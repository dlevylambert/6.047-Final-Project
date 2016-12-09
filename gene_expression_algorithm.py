import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.mlab import PCA
import numpy as np
import csv
from sklearn.cluster import KMeans
from sklearn import datasets, linear_model

line_num = 0
gene_expression = []
genes = []
patients = []

# open gene expression file and read in genes, and gene expression data
with open("BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt") as tsv:
    for line in csv.reader(tsv, dialect="excel-tab"):
        if line_num == 0:
            patients = line[1:]
            line_num += 1
        elif line_num == 1:
            line_num += 1
        else:
            genes.append(line[0])
            gene_expression.append(line[1:])

patients = [patient.split("-")[2].lower() for patient in patients]

# make everything into a float
for a in range(len(gene_expression)):
    for c in range(len(gene_expression[a])):
        gene_expression[a][c] = float(gene_expression[a][c])

dataMatrix = np.array(np.asarray(gene_expression).tolist())
result = PCA(dataMatrix)

highest_frac_genes_ind = []
highest_frac_genes = []

fracs = result.fracs.tolist()

#figure out the 10 most important fracs indices
for i in range(10):
    max_ind = fracs.index(max(fracs))
    highest_frac_genes_ind.append(max_ind)
    highest_frac_genes.append(genes[max_ind])
    fracs[max_ind] = 0

gene_expression_transposed = np.asarray(gene_expression).T.tolist()
max_columns = []
for patient_info in gene_expression_transposed:
    max_columns.append([patient_info[ind] for ind in highest_frac_genes_ind])

num_clusters_to_try = [2, 3, 4, 5, 6, 7] #unsure?
# try out a bunch of different cluster sizes?
for num_clusters in num_clusters_to_try:
    pass
    #kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(max_columns)
    # TODO KATY: collect stats on these clusters, look at kmeans.labels_, etc

clinical_patients = []
days_to_death = []
line_num = 0
with open("BRCA.clin.merged.picked.txt") as tsv:
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
num_clusters_test = 5
kmean_test = KMeans(n_clusters=num_clusters_test, random_state=0).fit(max_columns)
patients_by_cluster = {}
for i in range(num_clusters_test):
    patients_by_cluster[i] = []
# add patients with death days to patients_by_cluster
print patients_to_death
print patients_to_death.keys()
print "LEN PATIENTS ", len(patients)
print "LEN LABELS ", len(kmean_test.labels_)
for i in range(len(kmean_test.labels_)):
    patient_cluster = kmean_test.labels_[i]
    patient_id = patients[i]
    print patient_id
    if patient_id in patients_to_death.keys():
        patients_by_cluster[patient_cluster] = patients_by_cluster[patient_cluster] + [patients_to_death[patient_id]]

print patients_by_cluster
# do mean/median/st dev/anova (if num cluster > 2), else do t_test (if 2)
means = {}
stdvs = {}
medians = {}
vars = {}
for i in range(num_clusters_test):
    means[i] = np.mean(np.array(patients_by_cluster[i]))
    stdvs[i] = np.std(np.array(patients_by_cluster[i]))
    medians[i] = np.median(np.array(patients_by_cluster[i]))
    #vars[i] = np.var(np.array(patients_by_cluster[i]))

print means
print stdvs
print medians
print vars
