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

# open gene expression file and read in genes, and gene expression data
with open("BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt") as tsv:
    for line in csv.reader(tsv, dialect="excel-tab"):
        if line_num == 0:
            line_num += 1
        elif line_num == 1:
            line_num += 1
        else:
            genes.append(line[0])
            gene_expression.append(line[1:])

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


print highest_frac_genes
max_columns = []
for patient_info in gene_expression:
    max_columns.append([patient_info[ind] for ind in highest_frac_genes_ind])

num_clusters_to_try = [2, 3, 4, 5, 6, 7] #unsure?
# try out a bunch of different cluster sizes?
for num_clusters in num_clusters_to_try:
    pass
    #kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(max_columns)
    # TODO KATY: collect stats on these clusters, look at kmeans._labels, etc


