import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.mlab import PCA
import numpy as np
import csv
from sklearn.cluster import KMeans
from sklearn import datasets, linear_model
from scipy.spatial.distance import euclidean
from scipy import stats


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
    print "NORMALIZE THE DATA"
    # scale the genes so everything out of 1
    for a in range(len(gene_expression)):
        for c in range(len(gene_expression[a])):
            gene_expression[a][c] = (gene_expression[a][c] - min_val_genes[a]) / (max_val_genes[a] - min_val_genes[a])

    print "NUM PATIENTS ", len(patients)
    print "NUM PATIENTS FILTERED", len(patients_filtered)
    print "NUM GENES BEFORE " + str(line_num - 2)
    print "NUM GENES FILTERED", len(gene_expression)

    return patients_filtered, genes, gene_expression

def get_most_var_genes_inds(gene_expression, genes, num_genes):
    std_dev = []
    print "CALCULATING ST DEV"
    for gene_data in gene_expression:
        std_dev.append(np.std(np.array(gene_data)))

    highest_frac_genes_ind = []
    highest_frac_genes = []


    #figure out the num_genes most important fracs indices
    for i in range(num_genes):
        max_ind = std_dev.index(max(std_dev))
        highest_frac_genes_ind.append(max_ind)
        highest_frac_genes.append(genes[max_ind])
        std_dev[max_ind] = 0

    print "IND ", highest_frac_genes_ind
    print "HIGHEST FRAC GENES", highest_frac_genes
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

def get_cluster_stats(num_clusters, patients, patients_to_death, gene_expr_max_var, gene_names, genes):
    print "GENE EXPRS DIM ", len(gene_expr_max_var), " BY ", len(gene_expr_max_var[0])
    print "RUN KMEANS WITH N CLUSTERS: ", num_clusters
    kmean_test = KMeans(n_clusters=num_clusters, random_state=0).fit(gene_expr_max_var)
    patients_by_cluster = {}
    print "SPLIT PATIENTS ON CLUSTER"
    print "CHECK SIZE LABELS ", len(kmean_test.labels_)
    patient_gene_expression_data = {}
    for i in range(num_clusters):
        patients_by_cluster[i] = []
        patient_gene_expression_data[i] = []
    # add patients with death days to patients_by_cluster

    for i in range(len(kmean_test.labels_)):
        patient_cluster = kmean_test.labels_[i]
        patient_id = patients[i]
        if patient_id in patients_to_death.keys():
            patients_by_cluster[patient_cluster] = patients_by_cluster[patient_cluster] + [patients_to_death[patient_id]]
            patient_gene_expression_data[patient_cluster] = patient_gene_expression_data[patient_cluster] + [gene_expr_max_var[i]]

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

    anova_p_value = 0
    if num_clusters == 3:
        anova_p_value = stats.f_oneway(patients_by_cluster[0], patients_by_cluster[1], patients_by_cluster[2]).pvalue
    if num_clusters == 4:
        anova_p_value = stats.f_oneway(patients_by_cluster[0], patients_by_cluster[1], patients_by_cluster[2], patients_by_cluster[3]).pvalue
    if num_clusters == 5:
        anova_p_value =  stats.f_oneway(patients_by_cluster[0], patients_by_cluster[1], patients_by_cluster[2], patients_by_cluster[3], patients_by_cluster[4]).pvalue
    if num_clusters == 6:
        anova_p_value =  stats.f_oneway(patients_by_cluster[0], patients_by_cluster[1], patients_by_cluster[2], patients_by_cluster[3], patients_by_cluster[4], patients_by_cluster[5]).pvalue
    if num_clusters == 7:
        anova_p_value =  stats.f_oneway(patients_by_cluster[0], patients_by_cluster[1], patients_by_cluster[2], patients_by_cluster[3], patients_by_cluster[4], patients_by_cluster[5], patients_by_cluster[6]).pvalue
    if num_clusters == 8:
        anova_p_value =  stats.f_oneway(patients_by_cluster[0], patients_by_cluster[1], patients_by_cluster[2], patients_by_cluster[3],   patients_by_cluster[4], patients_by_cluster[5], patients_by_cluster[6], patients_by_cluster[7]).pvalue
    if num_clusters == 9:
        anova_p_value =  stats.f_oneway(patients_by_cluster[0], patients_by_cluster[1], patients_by_cluster[2], patients_by_cluster[3],   patients_by_cluster[4], patients_by_cluster[5], patients_by_cluster[6], patients_by_cluster[7], patients_by_cluster[8]).pvalue
    if num_clusters == 10:
        anova_p_value = stats.f_oneway(patients_by_cluster[0], patients_by_cluster[1], patients_by_cluster[2], patients_by_cluster[3],   patients_by_cluster[4], patients_by_cluster[5], patients_by_cluster[6], patients_by_cluster[7], patients_by_cluster[8], patients_by_cluster[9]).pvalue
    print "ANOVA P VALUE: ", anova_p_value

    if anova_p_value < 0.05:
        max_cluster = max(means)
        min_cluster = min(means)

        label_data = []
        gene_data = []
        cluster_use = max_cluster
        for i in range(num_clusters):
            for j in range(len(patient_gene_expression_data[i])):
                if i == cluster_use:
                    label_data = label_data + [1]
                else:
                    label_data = label_data + [0]
                gene_data.append(patient_gene_expression_data[i][j])

        reg_model = linear_model.LinearRegression()
        reg_model.fit(gene_data, label_data)
        model_coef_abs_map = map(abs, reg_model.coef_.tolist())
        max_weights_gene_weight = []
        max_weights_inds = []
        max_weights_gene = []
        for i in range(10):
            max_ind = model_coef_abs_map.index(max(model_coef_abs_map))
            max_weights_inds.append(max_ind)
            max_weights_gene_weight.append(max(model_coef_abs_map))
            max_weights_gene.append(gene_names[max_ind])
            model_coef_abs_map[max_ind] = 0

        max_gene_ind = max_weights_inds[0]
        cluster_val = 0
        all_others_val = 0
        for i in range(num_clusters):
            for j in range(len(patient_gene_expression_data[i])):
                if i == cluster_use:
                    cluster_val += patient_gene_expression_data[i][j][max_gene_ind]
                else:
                    all_others_val += patient_gene_expression_data[i][j][max_gene_ind]
        print "CLUSTER AVERAGE EXP OF GENE " +  str(max_weights_gene[0]) + " IS " + str(float(cluster_val)/len(patient_gene_expression_data[cluster_use]))
        print "ALL OTHER AVERAGE EXP OF GENE " +  str(max_weights_gene[0]) + " IS " + str(float(all_others_val)/(len(patients) -len(patient_gene_expression_data[cluster_use])))

        print max_weights_gene
        print max_weights_gene_weight




def cluster_diameters_by_size(gene_expr_max_var):
    num_clusters_to_diameter = {}
    num_clusters_to_try = [2, 3, 4, 5, 6, 7, 8, 9, 10]
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

        for i in range(len(kmeans.labels_)):
            patient_cluster_label = kmeans.labels_[i]
            num_patients_by_cluster[patient_cluster_label] += 1
            patient_cluster_center = kmeans.cluster_centers_[patient_cluster_label]
            # add to the cluster distances for that patient's cluster the euclidean distance to that center
            cluster_distances[patient_cluster_label] += euclidean(gene_expr_max_var[i], patient_cluster_center)

        for i in range(num_clusters):
            cluster_distances[i] /= float(num_patients_by_cluster[i])
        num_clusters_to_diameter[num_clusters] = sum(cluster_distances.values()) / float(num_clusters)
    return num_clusters_to_diameter

def plot_diameter_vs_num_clusters(num_clusters_to_diameter, cancer):
    x_is = num_clusters_to_diameter.keys()
    y_is = []

    for cluster_size in x_is:
        y_is.append(num_clusters_to_diameter[cluster_size])

    plt.plot(x_is, y_is)
    plt.xlabel("Number of Clusters")
    plt.ylabel("Diameter")
    plt.axis([0,10,0, max(num_clusters_to_diameter.values())])
    plt.title("Average Diameter Across Cluster Sizes " + cancer)
    plt.grid(True)
    plt.show()


def regression_death_data(cancer):
    patients, genes, gene_expr = get_patients_genes_and_gene_expression(cancer)
    patients_to_death = get_patients_to_death(cancer)
    gene_transp = np.asarray(gene_expr).T.tolist()
    
    label_data = []
    gene_data = []
    for i in range(len(patients)):
        patient_id = patients[i]
        label_data = label_data + [patients_to_death[patient_id]]
        gene_data.append(gene_transp[i])

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

def main():
    pass
    #cancer = "PAAD"
    #regression_death_data(cancer)
    #patients, genes, gene_expression = get_patients_genes_and_gene_expression(cancer)
    #highest_frac_genes, highest_frac_genes_ind = get_most_var_genes_inds(gene_expression, genes, 20)
    #gene_names_brca_var = ['C5orf41|153222', 'CACHD1|57685',  'FAM126A|84668', 'CBX7|23492', 'CXCL12|6387',  'BMX|660', 'RHOJ|57381', 'ZFP37|7539', 'PEAR1|375033', 'FER|2241']
    #gene_names_lusc_var = ['TNFSF12|8742', 'PTPRM|5797', 'CELF2|10659', 'RBPMS|11030', 'CACNA2D2|9254', 'RILPL2|196383', 'QSOX1|5768', 'GNA14|9630', 'SLC39A8|64116', 'MITF|4286']
    #gene_names_paad_var = ['UTY|7404', 'PARD3B|117583', 'PTPRG|5793', 'PRSS48|345062', 'RAD54L2|23132', 'ASXL2|55252', 'PTK6|5753', 'UPK1B|7348', 'CDK11A|728642', 'USP17L2|377630']
    #gene_names_brca_death = ['LOC201651|201651',  'USP17L2|377630', 'LOC134466|134466', 'FAM71A|149647', 'SPTA1|6708', 'DOC2B|8447', 'TSPO2|222642', 'GPR4|2828', 'SPDYE7P|441251', 'IFITM1|8519']
    #gene_names_lusc_death = ['LOC149620|149620', 'CD300LB|124599', 'NTM|50863', 'AMY2A|279', '?|442388', 'CELA3A|10136', 'AMY1A|276', 'SPDYE5|442590', 'TARSL2|123283', 'FAM21B|55747']
    #gene_names_paad_death = ['CACNG8|59283', 'PPP1R3A|5506', 'GP5|2814', 'T|6862', 'ISL2|64843', 'SPNS2|124976', 'FGF18|8817', 'S100A5|6276', 'KIAA0513|9764', 'FNDC8|54752']
    #gene_names = gene_names_paad_death
    #gene_expr_max_var = get_gene_exp_for_highest_var_genes(gene_expression, gene_names, genes)
    #patients_to_death = get_patients_to_death(cancer)
    #num_clusters_to_diameter = cluster_diameters_by_size(gene_expr_max_var)
    #plot_diameter_vs_num_clusters(num_clusters_to_diameter, cancer + " Death")
    #num_clusters = 6
    #get_cluster_stats(num_clusters, patients, patients_to_death, gene_expr_max_var, gene_names, genes)

main()
