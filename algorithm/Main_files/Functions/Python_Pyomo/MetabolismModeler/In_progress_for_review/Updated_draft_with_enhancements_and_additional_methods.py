import numpy as np
from sklearn.metrics import jaccard_score
from itertools import combinations
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import SpectralClustering, KMeans
from concurrent.futures import ThreadPoolExecutor
from scipy.sparse import lil_matrix

def create_jaccard_sample_x_sample_matrix(samples):
    n_samples = len(samples)
    jaccard_matrix = lil_matrix((n_samples, n_samples))
    
    def jaccard_for_pair(i, j):
        # Create a union of the sets from samples[i] and samples[j]
        union = samples[i].union(samples[j])
        # Create binary vectors relative to the union
        vector_i = np.array([1 if item in samples[i] else 0 for item in union])
        vector_j = np.array([1 if item in samples[j] else 0 for item in union])
        # Compute the Jaccard index for the pair
        return jaccard_score(vector_i, vector_j, average='binary')
    
    with ThreadPoolExecutor() as executor:
        futures = []
        for i in range(n_samples):
            for j in range(i, n_samples):
                futures.append(executor.submit(jaccard_for_pair, i, j))
        
        idx = 0
        for i in range(n_samples):
            for j in range(i, n_samples):
                jaccard_matrix[i, j] = futures[idx].result()
                jaccard_matrix[j, i] = jaccard_matrix[i, j]
                idx += 1

    return jaccard_matrix.toarray()

def optimize_order_hierarchical(jaccard_matrix):
    Z = linkage(jaccard_matrix, 'average')
    optimal_order = leaves_list(Z)
    optimized_jaccard_matrix = jaccard_matrix[optimal_order, :][:, optimal_order]
    return optimized_jaccard_matrix, optimal_order, Z

def optimize_order_greedy(jaccard_matrix):
    n_samples = len(jaccard_matrix)
    unvisited = set(range(n_samples))
    current_index = 0
    optimal_order = [current_index]
    unvisited.remove(current_index)

    while unvisited:
        # Select the next index with the maximum similarity (Maximum Jaccard index)
        next_index = max(unvisited, key=lambda x: jaccard_matrix[current_index, x])
        optimal_order.append(next_index)
        unvisited.remove(next_index)
        current_index = next_index

    optimized_jaccard_matrix = jaccard_matrix[np.ix_(optimal_order, optimal_order)]
    return optimized_jaccard_matrix, optimal_order

def optimize_order_sklearn(jaccard_matrix, n_clusters=2):
    clustering = SpectralClustering(n_clusters=n_clusters, affinity='precomputed', assign_labels='discretize', random_state=0)
    labels = clustering.fit_predict(1 - jaccard_matrix)
    sorted_indices = np.argsort(labels)
    optimized_jaccard_matrix = jaccard_matrix[sorted_indices, :][:, sorted_indices]
    return optimized_jaccard_matrix, sorted_indices

def optimize_order_kmeans(jaccard_matrix, n_clusters=2):
    clustering = KMeans(n_clusters=n_clusters, random_state=0)
    labels = clustering.fit_predict(jaccard_matrix)
    sorted_indices = np.argsort(labels)
    optimized_jaccard_matrix = jaccard_matrix[sorted_indices, :][:, sorted_indices]
    return optimized_jaccard_matrix, sorted_indices

def visualize_optimized_jaccard_matrix(optimized_jaccard_matrix, title, cmap='viridis'):
    plt.figure(figsize=(10, 10))
    sns.heatmap(optimized_jaccard_matrix, cmap=cmap, annot=False, cbar=True)
    plt.title(title)
    plt.xlabel('Sample Index (Reordered)')
    plt.ylabel('Sample Index (Reordered)')
    plt.show()

def visualize_dendrogram(Z, labels):
    plt.figure(figsize=(12, 8))
    dendrogram(Z, labels=[f'Sample {i}' for i in labels], leaf_rotation=90, leaf_font_size=8)
    plt.title('Dendrogram of Hierarchical Clustering')
    plt.xlabel('Sample Index')
    plt.ylabel('Distance')
    plt.show()

if __name__ == "__main__":
    # Example test sets (arrays of sets containing strings)
    # original_samples = [
    #     {"gene1", "gene2", "gene3"},
    #     {"gene2", "gene3"},
    #     # {"gene1", "gene4", "gene5"},
    #     # {"gene3", "gene5", "gene6"},
    #     # {"gene1", "gene2", "gene6"},
    #     # {"gene2", "gene3", "gene5"},
    # ]
    test_set1 = {"rxn1", "rxn2", "rxn3", "rxn4", "rxn11"} 
    test_set2 = {"rxn3", "rxn4", "rxn5", "rxn6"}
    test_set3 = {"rxn1", "rxn2", "rxn3", "rxn5"}
    test_set4 = {"rxn2", "rxn3", "rxn4", "rxn5", "rxn6"} 
    test_set5 = {"rxn3", "rxn2","rxn11", "rxn6"}
    test_set6 = {"rxn3", "rxn5","rxn11", "rxn6", "rxn88", "rxn5"}
    test_set7 = {"rxn11", "rxn1", "rxn2"}
    test_set8 = {"rxn11", "rxn5", "rxn2", "rxn99"}
    test_set9 = {"rxn11", "rxn5", "rxn88"}
    # original_samples = [test_set8, test_set9]
    original_samples = [test_set1, test_set2, test_set3, test_set4, test_set5, test_set6, test_set7, test_set8, test_set9]
    # Create 322 samples by duplicating the existing test sets
    num_samples = 322
    samples = original_samples
    # samples = original_samples * (num_samples // len(original_samples))
    # samples = samples[:num_samples]
    print(samples)
    
    print(f"Total number of samples created: {len(samples)}")

    print("\nCreating Jaccard matrix...")
    jaccard_matrix = create_jaccard_sample_x_sample_matrix(samples)
    print(jaccard_matrix)
    print("\nJaccard matrix created.")

    # Hierarchical method
    print("\nOptimizing with Hierarchical Clustering...")
    optimized_matrix_hierarchical, optimal_order_hierarchical, Z = optimize_order_hierarchical(jaccard_matrix)
    visualize_optimized_jaccard_matrix(optimized_matrix_hierarchical, 'Optimized Jaccard Matrix (Hierarchical)')
    visualize_dendrogram(Z, optimal_order_hierarchical)

    # Greedy method
    print("\nOptimizing with Greedy Approach...")
    optimized_matrix_greedy, optimal_order_greedy = optimize_order_greedy(jaccard_matrix)
    visualize_optimized_jaccard_matrix(optimized_matrix_greedy, 'Optimized Jaccard Matrix (Greedy)')
    
    # Sklearn method
    print("\nOptimizing with Sklearn Spectral Clustering...")
    optimized_matrix_sklearn, optimal_order_sklearn = optimize_order_sklearn(jaccard_matrix, n_clusters=2)
    visualize_optimized_jaccard_matrix(optimized_matrix_sklearn, 'Optimized Jaccard Matrix (Sklearn Spectral Clustering)')
    
    # KMeans method
    print("\nOptimizing with KMeans Clustering...")
    optimized_matrix_kmeans, optimal_order_kmeans = optimize_order_kmeans(jaccard_matrix, n_clusters=2)
    visualize_optimized_jaccard_matrix(optimized_matrix_kmeans, 'Optimized Jaccard Matrix (KMeans Clustering)')

    print("\nFinal optimal orders:")
    print(f"Hierarchical: {optimal_order_hierarchical}")
    print(f"Greedy: {optimal_order_greedy}")
    print(f"Sklearn: {optimal_order_sklearn}")
    print(f"KMeans: {optimal_order_kmeans}")
