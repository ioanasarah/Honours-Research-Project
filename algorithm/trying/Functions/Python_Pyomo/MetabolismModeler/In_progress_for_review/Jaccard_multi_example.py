from sklearn.metrics import jaccard_score
from itertools import combinations
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import SpectralClustering, KMeans
from concurrent.futures import ThreadPoolExecutor
from scipy.sparse import lil_matrix
import numpy as np


test_set1 = {"rxn1", "rxn2","rxn33", "rxn3", "rxn4", "rxn5", "rxn11","rxn77","rxn88", "rxn99", "rxn111", "rxn12"} 
test_set2 = {                        "rxn3", "rxn4", "rxn5", "rxn11","rxn77","rxn88", "rxn99", "rxn111",          "rnx13"}

test_set1 = {"rxn1", "rexn", "rdd", "r", "2", "3"}
test_set2 = {"rxn1", "rddd"}

samples = [
		test_set1, 
		test_set2
		]

n_samples = len(samples)
jaccard_matrix = lil_matrix((n_samples, n_samples))
def jaccard_for_pair(i, j, average):
    # Create a union of the sets from samples[i] and samples[j]
    union = samples[i].union(samples[j])
    # Create binary vectors relative to the union
    vector_i = np.array([1 if item in samples[i] else 0 for item in union])
    vector_j = np.array([1 if item in samples[j] else 0 for item in union])
    # Compute the Jaccard index for the pair
    return jaccard_score(vector_i, vector_j, average=average)


print(jaccard_for_pair(0, 1, "micro"))
print(jaccard_for_pair(0, 1, "binary"))