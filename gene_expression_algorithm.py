import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.mlab import PCA
import numpy as np
import csv
from sklearn.cluster import KMeans
from sklearn import datasets, linear_model
from scipy.spatial.distance import euclidean


def get_patients_genes_and_gene_expression(cancer):
    line_num = 0
    gene_expression_all = []
    genes = []
    patients = []
    patient_ids_to_include = get_patients_to_death(cancer).keys()
    
    # open gene expression file and read in genes, and gene expression data
    # structured rows: genes, columns: patients
    max_val_genes = []
    min_val_genes = []
    print "READING THE DATA"
    with open(cancer + ".rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt") as tsv:
        for line in csv.reader(tsv, dialect="excel-tab"):
            if line_num == 0:
                patients = line[1:]
                line_num += 1
            elif line_num == 1:
                line_num += 1
            else:
                gene_data = [float(d) for d in line[1:]]
                if max(gene_data) != 0:
                    genes.append(line[0])
                    gene_expression_all.append(gene_data)
                    max_val_genes.append(max(gene_data))
                    min_val_genes.append(min(gene_data))
                line_num += 1

    # extract the unique patient label from the patient ID
    patients = [patient.split("-")[2].lower() for patient in patients]
    patients_filtered = []

    print "FILTER ONLY PATIENTS WITH DEATH "
    gene_expression_transposed = np.asarray(gene_expression_all).T.tolist()
    gene_expression_transposed_filtered = []
    for i in range(len(patients)):
        if patients[i] in patient_ids_to_include:
            patients_filtered.append(patients[i])
            gene_expression_transposed_filtered.append(gene_expression_transposed[i])

    gene_expression = np.asarray(gene_expression_transposed_filtered).T.tolist()
    print len(gene_expression), len(gene_expression[0])
    print "NORMALIZE THE DATA"
    # scale the genes so everything out of 1
    for a in range(len(gene_expression)):
        for c in range(len(gene_expression[a])):
            gene_expression[a][c] = (gene_expression[a][c] - min_val_genes[a]) / (max_val_genes[a] - min_val_genes[a])


    print "NUM PATIENTS ", len(patients)
    print "NUM GENES ", len(gene_expression)

    return patients_filtered, genes, gene_expression

def get_10_most_var_genes_inds(gene_expression, genes):
    print len(gene_expression[0])
    print len(gene_expression)
    
    std_dev = []
    print "CALCULATING ST DEV"
    for gene_data in gene_expression:
        std_dev.append(np.std(np.array(gene_data)))

    highest_frac_genes_ind = []
    highest_frac_genes = []

    print len(std_dev)

    #figure out the 10 most important fracs indices
    for i in range(20):
        max_ind = std_dev.index(max(std_dev))
        highest_frac_genes_ind.append(max_ind)
        highest_frac_genes.append(genes[max_ind])
        std_dev[max_ind] = 0

    print "IND ", highest_frac_genes_ind
    print highest_frac_genes
    return highest_frac_genes, highest_frac_genes_ind

def get_gene_exp_for_highest_var_genes(gene_expression, gene_names, genes):
    highest_frac_genes_ind = []
    for name in gene_names:
        highest_frac_genes_ind.append(genes.index(name))
    
    print "IND ", highest_frac_genes_ind
    # structured row: patients, columns: genes
    gene_expression_transposed = np.asarray(gene_expression).T.tolist()

    gene_expr_max_var = []
    for patient_info in gene_expression_transposed:
        gene_expr_max_var.append([patient_info[ind] for ind in highest_frac_genes_ind])
    return gene_expr_max_var
    #return PCA(np.asarray(gene_expr_max_var)).Y.tolist()

def get_patients_to_death(cancer):
    clinical_patients = []
    days_to_death = []
    line_num = 0
    with open(cancer + ".clin.merged.picked.txt") as tsv:
        for line in csv.reader(tsv, dialect="excel-tab"):
            if line_num == 0:
                clinical_patients = line[1:]
                line_num += 1
            elif line_num == 1:
                line_num += 1
            else:
                if line[0] == "days_to_death":
                    days_to_death = line[1:]

    clinical_patients = [c_patient.split("-")[2] for c_patient in clinical_patients]
    patients_to_death = {clinical_patients[i]: float(days_to_death[i]) for i in range(len(days_to_death)) if days_to_death[i] != "NA"}
    average = np.mean(np.array(patients_to_death.values()))
    std_dev = np.std(np.array(patients_to_death.values()))
    diff_val = (max(patients_to_death.values())- min(patients_to_death.values()) - average)/std_dev
    print "AVERAGE DEATH ", average
    print "STD DEATH", std_dev
    patients_to_death = {patient :  (patients_to_death[patient] - average)/std_dev for patient in patients_to_death.keys()}
    return patients_to_death

def get_cluster_stats(num_clusters, patients, patients_to_death, gene_expr_max_var):
    print "GENE EXPRS DIM ", len(gene_expr_max_var), " BY ", len(gene_expr_max_var[0])
    print "RUN KMEANS WITH N CLUSTERS: ", num_clusters
    kmean_test = KMeans(n_clusters=num_clusters, random_state=0).fit(gene_expr_max_var)
    patients_by_cluster = {}
    print "SPLIT PATIENTS ON CLUSTER"
    print "CHECK SIZE LABELS ", len(kmean_test.labels_)
    for i in range(num_clusters):
        patients_by_cluster[i] = []
    # add patients with death days to patients_by_cluster

    for i in range(len(kmean_test.labels_)):
        patient_cluster = kmean_test.labels_[i]
        patient_id = patients[i]
        if patient_id in patients_to_death.keys():
            patients_by_cluster[patient_cluster] = patients_by_cluster[patient_cluster] + [patients_to_death[patient_id]]

    print patients_by_cluster
    # do mean/median/st dev/anova (if num cluster > 2), else do t_test (if 2)
    means = {}
    stdvs = {}
    medians = {}
    vars = {}
    for i in range(num_clusters):
        means[i] = np.mean(np.array(patients_by_cluster[i]))
        stdvs[i] = np.std(np.array(patients_by_cluster[i]))
        medians[i] = np.median(np.array(patients_by_cluster[i]))

    print "MEANS : ", means
    print "STD DEV : ", stdvs
    print "MEDIANS : ", medians

def cluster_diameters_by_size(gene_expr_max_var):
    num_clusters_to_diameter = {}
    num_clusters_to_try = [2, 3, 4, 5, 6, 7, 8, 9, 10] #unsure?
    # try out a bunch of different cluster sizes?
    for num_clusters in num_clusters_to_try:
        kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(gene_expr_max_var)
        #creates list of lists, where each list represents a cluster and which point it belongs to.

        cluster_avgs = {}
        cluster_distances = {}
        num_patients_by_cluster = {}
        for i in range(num_clusters):
            cluster_distances[i] = 0
            num_patients_by_cluster[i] = 0
            cluster_avgs[i] = []
        
        #for i in range(len(kmeans.labels_)):
        #    patient_cluster_label = kmeans.labels_[i]
        #    patient_cluster_center = kmeans.cluster_centers_[patient_cluster_label]
        #    num_patients_by_cluster[patient_cluster_label] += 1
        #    if len(cluster_avgs[patient_cluster_label]) == 0:
        #        cluster_avgs[patient_cluster_label] = gene_expr_max_var[i]
        #    else:
        #        for j in gene_expr_max_var[i]:
        #            cluster_avgs[patient_cluster_label][j] += gene_expr_max_var[i][j]


        for i in range(len(kmeans.labels_)):
            patient_cluster_label = kmeans.labels_[i]
            num_patients_by_cluster[patient_cluster_label] += 1
            #patient_cluster_center = cluster_avgs[i]
            patient_cluster_center = kmeans.cluster_centers_[patient_cluster_label]
            # add to the cluster distances for that patient's cluster the euclidean distance to that center
            cluster_distances[patient_cluster_label] += euclidean(gene_expr_max_var[i], patient_cluster_center)

        for i in range(num_clusters):
            cluster_distances[i] /= float(num_patients_by_cluster[i])
        num_clusters_to_diameter[num_clusters] = sum(cluster_distances.values()) / float(num_clusters)
    print num_clusters_to_diameter
    return num_clusters_to_diameter

def plot_diameter_vs_num_clusters(num_clusters_to_diameter):
    x_is = num_clusters_to_diameter.keys()
    y_is = []

    for cluster_size in x_is:
        y_is.append(num_clusters_to_diameter[cluster_size])
    print x_is
    print y_is
    plt.plot(x_is, y_is)
    plt.xlabel("Number of Clusters")
    plt.ylabel("Diameter")
    plt.axis([0,8,0, max(num_clusters_to_diameter.values())])
    plt.title("Average Diameter Across Cluster Sizes")
    plt.grid(True)
    plt.show()


def regression_death_data(cancer):
    patients, genes, gene_expr = get_patients_genes_and_gene_expression(cancer)
    patients_to_death = get_patients_to_death(cancer)
    gene_transp = np.asarray(gene_expr).T.tolist()
    print len(gene_transp)
    print len(gene_transp[0])
    
    label_data = []
    gene_data = []
    for i in range(len(patients)):
        patient_id = patients[i]
        label_data = label_data + [patients_to_death[patient_id]]
        gene_data.append(gene_transp[i])

    #print gene_data[0]


    reg_model = linear_model.LinearRegression()
    reg_model.fit(gene_data, label_data)
    model_coef_abs_map = map(abs, reg_model.coef_.tolist())
    max_weights_gene_ind = []
    max_weights_gene = []
    for i in range(20):
        max_ind = model_coef_abs_map.index(max(model_coef_abs_map))
        max_weights_gene_ind.append(max_ind)
        max_weights_gene.append(genes[max_ind])
        model_coef_abs_map[max_ind] = 0
    print max_weights_gene

def main():
    cancer = "BRCA"
    #regression_death_data(cancer)
    patients, genes, gene_expression = get_patients_genes_and_gene_expression(cancer)
    #highest_frac_genes, highest_frac_genes_ind = get_10_most_var_genes_inds(gene_expression, genes)
    gene_names = ['TAS2R50|259296', 'C15orf5|81698', 'FLJ14107|80094', 'LOC201651|201651', 'USP17L2|377630', 'C11orf94|143678', 'LOC134466|134466', 'OR2AK2|391191', 'FAM71A|149647', 'SPTA1|6708', 'SNORA2A|677793', 'DNASE1L3|1776', 'DOC2B|8447', 'OR10V1|390201', 'OR52N1|79473', 'TSPO2|222642', 'OR52I2|143502', 'GPR4|2828', 'SPDYE7P|441251', 'IFITM1|8519']
    gene_expr_max_var = get_gene_exp_for_highest_var_genes(gene_expression, gene_names, genes)
    patients_to_death = get_patients_to_death(cancer)
    num_clusters_to_diameter = cluster_diameters_by_size(gene_expr_max_var)
    plot_diameter_vs_num_clusters(num_clusters_to_diameter)
    num_clusters = 4
    get_cluster_stats(num_clusters, patients, patients_to_death, gene_expr_max_var)

main()
