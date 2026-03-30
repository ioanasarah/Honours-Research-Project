# # # ### Grant_figures.py


# # #################### heatmap
# import numpy as np
# import seaborn as sns
# import matplotlib.pyplot as plt
# import pandas as pd
# from scipy.cluster.hierarchy import linkage
# from matplotlib.colors import LinearSegmentedColormap




# # Set seed for reproducibility
# np.random.seed(42)

# # Parameters
# num_genes = 200
# num_samples_per_group = 1
# num_special_genes = 45  # The number of genes with lower expression in group 2

# # Group 1: Generate random gene expression data (from a normal distribution)
# group_1 = np.random.normal(loc=9.8, scale=1, size=(num_genes//2, num_samples_per_group))
# group_1 = np.vstack([group_1, np.random.normal(loc=11.8, scale=1.2, size=(num_genes//2, num_samples_per_group))])

# # Group 2: Generate random gene expression data with similar distribution, but 
# # lower mean for a subset of genes (lower by 2 units on average for 30 genes)
# group_2 = np.random.normal(loc=10.2, scale=1, size=(num_genes//2 - num_special_genes, num_samples_per_group))
# group_2 = np.vstack([group_2, np.random.normal(loc=13.8, scale=1.4, size=(num_special_genes, num_samples_per_group))])
# group_2 = np.vstack([group_2, np.random.normal(loc=7.8, scale=2, size=(num_genes//2, num_samples_per_group))])

# group_3 = np.random.normal(loc=11.2, scale=1.2, size=(num_genes//2 - num_special_genes, num_samples_per_group))
# group_3 = np.vstack([group_3, np.random.normal(loc=12.5, scale=1.4, size=(num_special_genes, num_samples_per_group))])
# group_3 = np.vstack([group_3, np.random.normal(loc=11.8, scale=2, size=(num_genes//2, num_samples_per_group))])



# # Combine data into one DataFrame for easier plotting
# combined_data = np.hstack([group_1, group_2, group_3])
# gene_names = [f"Gene_{i+1}" for i in range(num_genes)]
# sample_names = [f"G1_Sample_{i+1}" for i in range(num_samples_per_group)] + [f"G2_Sample_{i+1}" for i in range(num_samples_per_group)] +[f"G3_Sample_{i+1}" for i in range(num_samples_per_group)]

# # Convert to a DataFrame
# df = pd.DataFrame(combined_data, index=gene_names, columns=sample_names)
# df_zscored = df.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
# n_bins = 200  # Number of bins in the colormap.
# colors = ['maroon', 'lightslategrey', "chartreuse"]
# custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=n_bins)
# sns.clustermap(df_zscored, cmap=custom_cmap, center=0, method='average', metric='euclidean', 
#                row_cluster=False, col_cluster=False, figsize=(12, 14))

# plt.show()


# ###########################

# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt

# # Creating a DataFrame from your data
# data = {
#     'Flux': [19.5, 19.5, 0.5, 19, 19, 20, 204, 10, 10, 10, 34.5, 10, 10, 10, 1, 24, 26.75, 34, 1],
#     'Expression': [138.1710216, 41.70666686, -1, 22.73741773, 19.95689046, 17.38122708, 85.58669581, 231.0760223, 13.82558113, 
#              87.37680513, 204.8092443, 239.9041021, 87.37680513, 29.96926758, 15.8725959, 192.6596178, 0.497111644, 
#              278.7773185, 1.279456944],
#     'Reaction ID': ['MAR04388_f', 'MAR08514', 'MAR06606_f', 'MAR08609', 'MAR08611', 'MAR00710_f', 'MAR03957', 
#                       'MAR04145', 'MAR04152_f', 'MAR04209', 'MAR04410_f', 'MAR04589_f', 'MAR06413', 'MAR06414_f', 
#                       'MAR03977', 'MAR06911_f', 'MAR06914', 'MAR06918', 'MAR00357'],
#     'Reaction': ["(S)-Lactate:NAD+ oxidoreductase", "(S)-Lactate:ferricytochrome-c 2-oxidoreductase", "", 
#                  "L-Proline:NADP+ 5-oxidoreductase", "L-proline:(acceptor) oxidoreductase", 
#                  "Isocitrate:NADP+ oxidoreductase (decarboxylating)", "isocitrate:NAD+ oxidoreductase (decarboxylating)", 
#                  "acetyl-CoA:oxaloacetate C-acetyltransferase (thioester-hydrolysing)", 
#                  "Succinate:CoA ligase (ADP-forming)", "", "(S)-malate hydro-lyase (fumarate-forming)", 
#                  "isocitrate hydro-lyase (cis-aconitate-forming)", "", 
#                  "succinyl-CoA:enzyme N6-(dihydrolipoyl)lysine S-succinyltransferase", 
#                  "diphosphate phosphohydrolase", "", "ferrocytochrome-c:oxygen oxidoreductase", "", ""],
#     'RxnFormula': [
#         "H+ + NADH + pyruvate --> L-lactate + NAD+", 
#         "2.0 ferricytochrome C + L-lactate --> 2.0 ferrocytochrome C + 2.0 H+ + pyruvate", 
#         "ascorbate + urate radical --> monodehydroascorbate + urate", 
#         "1-pyrroline-5-carboxylate + 2.0 H+ + NADPH --> NADP+ + proline", 
#         "FAD + proline --> 1-pyrroline-5-carboxylate + FADH2 + H+", 
#         "isocitrate + NADP+ --> AKG + CO2 + NADPH", 
#         "isocitrate + NAD+ --> AKG + CO2 + NADH", 
#         "acetyl-CoA + H2O + OAA --> citrate + CoA + H+", 
#         "ADP + Pi + succinyl-CoA --> ATP + CoA + succinate", 
#         "AKG + H+ + thiamin-PP --> 3-carboxy-1-hydroxypropyl-ThPP + CO2", 
#         "fumarate + H2O --> malate", 
#         "cis-aconitate + H2O --> isocitrate", 
#         "3-carboxy-1-hydroxypropyl-ThPP + lipoamide --> S-succinyldihydrolipoamide + thiamin-PP", 
#         "CoA + S-succinyldihydrolipoamide --> dihydrolipoamide + succinyl-CoA", 
#         "H2O + PPi --> H+ + 2.0 Pi", 
#         "FADH2 + ubiquinone --> FAD + ubiquinol", 
#         "4.0 ferrocytochrome C + 8.0 H+ + O2 --> 4.0 ferricytochrome C + 4.0 H+ + 2.0 H2O", 
#         "2.0 ferricytochrome C + 2.0 H+ + ubiquinol --> 2.0 ferrocytochrome C + 4.0 H+ + ubiquinone", 
#         "ATP + CoA + omega-3-arachidonic acid --> (8Z,11Z,14Z,17Z)-eicosatetraenoyl-CoA + AMP + PPi"
#     ],
# }

# df = pd.DataFrame(data)

# styled_df = df.style.set_table_styles([
#     {'selector': 'thead th', 'props': [('background-color', '#4CAF50'), ('color', 'white'), ('font-weight', 'bold')]},
#     {'selector': 'tbody td', 'props': [('border', '1px solid black')]}
# ]).set_properties(**{'text-align': 'center'}).hide(axis = "index")
# styled_df
# styled_df.to_html("styled_table.html")


# ################################################################ UMAP AND PCA

# import numpy as np
# import seaborn as sns
# import matplotlib.pyplot as plt
# import pandas as pd
# import umap
# from scipy.cluster.hierarchy import linkage
# from sklearn.decomposition import PCA
# from matplotlib.colors import LinearSegmentedColormap

# # Set seed for reproducibility
# np.random.seed(42)

# # Parameters
# num_genes = 300
# num_samples_per_group = 8
# num_special_genes = 45  # The number of genes with lower expression in group 2

# # Group 1: Generate random gene expression data (from a normal distribution)
# group_1 = np.random.normal(loc=9.8, scale=1, size=(num_genes//2, num_samples_per_group))
# group_1 = np.vstack([group_1, np.random.normal(loc=11.8, scale=1.2, size=(num_genes//2, num_samples_per_group))])

# # Group 2: Generate random gene expression data with similar distribution, but 
# # lower mean for a subset of genes (lower by 2 units on average for 30 genes)
# group_2 = np.random.normal(loc=10.2, scale=1, size=(num_genes//2 - num_special_genes, num_samples_per_group))
# group_2 = np.vstack([group_2, np.random.normal(loc=13.8, scale=1.4, size=(num_special_genes, num_samples_per_group))])
# group_2 = np.vstack([group_2, np.random.normal(loc=7.8, scale=2, size=(num_genes//2, num_samples_per_group))])

# group_3 = np.random.normal(loc=11.2, scale=1.2, size=(num_genes//2 - num_special_genes, num_samples_per_group))
# group_3 = np.vstack([group_3, np.random.normal(loc=12.5, scale=1.4, size=(num_special_genes, num_samples_per_group))])
# group_3 = np.vstack([group_3, np.random.normal(loc=11.8, scale=2, size=(num_genes//2, num_samples_per_group))])



# # Combine data into one DataFrame for easier plotting
# combined_data = np.hstack([group_1, group_2, group_3])
# genes = [f"Gene_{i+1}" for i in range(num_genes)]
# samples = [f"G1_Sample_{i+1}" for i in range(num_samples_per_group)] + [f"G2_Sample_{i+1}" for i in range(num_samples_per_group)] +[f"G3_Sample_{i+1}" for i in range(num_samples_per_group)]

# # Convert to a DataFrame
# df = pd.DataFrame(combined_data, index=genes, columns=samples)
# df_zscored = df.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
# n_bins = 100  # Number of bins in the colormap.
# colors = ['green', 'black', "red"]

# colors = [ "LightBlue", "Magenta", "Orange"]
# # UMAP Clustering
# umap_model = umap.UMAP(n_neighbors=10, min_dist=0.3, random_state=42)
# umap_results = umap_model.fit_transform(df_zscored.T)  # Transpose to cluster samples

# # PCA
# pca = PCA(n_components=2)
# pca_results = pca.fit_transform(df_zscored.T)  # Transpose to cluster samples

# # Create DataFrame for plotting
# df_umap = pd.DataFrame(umap_results, columns=['UMAP1', 'UMAP2'])
# df_umap['Sample'] = samples
# df_umap["Group"] = ["Healthy"if sample.startswith("G1") else "Disease" if sample.startswith("G2") else "Disease + Compensated" for sample in samples ]

# df_pca = pd.DataFrame(pca_results, columns=['PC1', 'PC2'])
# df_pca['Sample'] = samples
# df_pca["Group"] = ["Healthy"if sample.startswith("G1") else "Disease" if sample.startswith("G2") else "Disease + Compensated" for sample in samples ]

# # Plot UMAP
# plt.figure(figsize=(10, 6))
# sns.scatterplot(x='UMAP1', y='UMAP2', hue='Group', data=df_umap, palette=colors, s=100)
# plt.title('UMAP Clustering of Samples')
# plt.show()

# # Plot PCA
# plt.figure(figsize=(10, 6))
# sns.scatterplot(x='PC1', y='PC2', hue='Group', data=df_pca, palette=colors, s=100)
# plt.title('PCA Clustering of Samples')
# plt.show()

# ################################################################ Volcano plot


# # Split into two groups
# import pandas as pd
# import numpy as np
# from scipy.stats import ttest_ind
# from statsmodels.stats.multitest import multipletests
# import matplotlib.pyplot as plt

# # Set seed for reproducibility
# np.random.seed(41)

# # Parameters
# num_genes = 501
# num_samples_per_group = 8
# num_special_genes = 45  # The number of genes with lower expression in group 2

# # Group 1: Generate random gene expression data (from a normal distribution)
# group_1 = np.random.normal(loc=7.8, scale=2.2, size=(num_genes//3, num_samples_per_group))
# group_1 = np.vstack([group_1, np.random.normal(loc=12.8, scale=1.3, size=(num_genes//3, num_samples_per_group))])
# group_1 = np.vstack([group_1, np.random.normal(loc=18, scale=3.4, size=(num_genes//3, num_samples_per_group))])


# np.random.seed(43)
# # Group 2: Generate random gene expression data with similar distribution, but 
# # lower mean for a subset of genes (lower by 2 units on average for 30 genes)
# group_2 = np.random.normal(loc=10.2, scale=1.8, size=(num_genes//3 - num_special_genes, num_samples_per_group))
# group_2 = np.vstack([group_2, np.random.normal(loc=13.8, scale=1.4, size=(num_special_genes, num_samples_per_group))])
# group_2 = np.vstack([group_2, np.random.normal(loc=6.8, scale=2, size=(num_genes//3, num_samples_per_group))])
# group_2 = np.vstack([group_2, np.random.normal(loc=17.8, scale=3.8, size=(num_genes//3, num_samples_per_group))])


# # Combine data into one DataFrame for easier plotting
# combined_data = np.hstack([group_1, group_2])
# genes = [f"Gene_{i+1}" for i in range(num_genes)]
# samples= [f"G1_Sample_{i+1}" for i in range(num_samples_per_group)] + [f"G2_Sample_{i+1}" for i in range(num_samples_per_group)]
# groups =  ["Control"if sample.startswith("G1") else "Disease" for sample in samples]

# # Convert to a DataFrame
# df = pd.DataFrame(combined_data, index=genes, columns=samples)
# # df_zscored = df.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
# df = df.T
# df['Group'] = groups
# df = df.T

# controls = df.columns[df.loc['Group'] == 'Control']
# diseases = df.columns[df.loc['Group'] == 'Disease']
# group1 = df[controls].drop('Group', axis=0)
# group2 = df[diseases].drop('Group', axis=0)

# fold_changes = []
# p_values = []

# # Perform t-test for each gene
# for gene in df.index:
#     if gene == "Group":
#         continue
#     sample1 = group1.loc[gene]
#     sample2 = group2.loc[gene]
#     # print(sample1)
#     # print(sample2)
    
#     # Compute fold change
#     mean1 = np.mean(sample1)
#     mean2 = np.mean(sample2)
#     # print(mean1, mean2)
#     fold_change = mean2 / mean1 if mean1 != 0 else np.nan
#     fold_changes.append(fold_change-1)
#     # print(fold_change)
#     sample1 = sample1.tolist()
#     sample2 = sample2.tolist()

#     # Perform t-test
#     _, p_value = ttest_ind(sample1, sample2, equal_var=False)
#     p_values.append(p_value)

# gene_index = df.index[:-1]
# fold_changes = pd.Series(fold_changes, index=gene_index)
# p_values = pd.Series(p_values, index=gene_index)

# _, p_values_corrected, _, _ = multipletests(p_values, method='fdr_bh')

# # Convert fold_changes and p_values_corrected to Series for better indexing
# results = pd.DataFrame({
#     'Fold Change': fold_changes,
#     'Log(-p)': -np.log10(p_values_corrected)
# })

# # Define criteria
# fold_change_threshold = 0.5
# p_value_threshold = 0.05

# # Determine which points to highlight
# highlight = (abs(results['Fold Change']) > fold_change_threshold) & (results['Log(-p)'] > -np.log10(p_value_threshold))

# # Create the volcano plot
# plt.figure(figsize=(10, 8))
# plt.scatter(
#     results['Fold Change'], 
#     results['Log(-p)'], 
#     alpha=0.5, 
#     edgecolors='w', 
#     s=60,
#     c=np.where(highlight, 'red', 'gray')  # Color points based on criteria
# )

# # Add significance threshold lines
# plt.axhline(y=-np.log10(p_value_threshold), color='r', linestyle='--',
#  # label='p-value threshold'
#  )
# plt.axvline(x=fold_change_threshold, color='b', linestyle='--', 
# 	#label='Fold Change Threshold (positive)'
# 	)
# plt.axvline(x=-fold_change_threshold, color='b', linestyle='--', 
# 	# label='Fold Ch ange Threshold (negative)'
# 	)
# plt.axvline(x=0, color='black', linestyle='-', 
# 	linewidth=2
# 	# label='Fold Ch ange Threshold (negative)'
# 	)
# # Add labels and title
# plt.xlabel('Fold Change')
# plt.ylabel('Log(-p-value)')
# plt.title('Volcano Plot')
# plt.grid(True)
# plt.legend()

# # Show the plot
# plt.show()

# ################################################################ Histogram data

# import numpy as np
# import matplotlib.pyplot as plt

# # Parameters for the log-normal distribution
# mean = 0     # Mean of the underlying normal distribution
# sigma = 1    # Standard deviation of the underlying normal distribution
# N = 3000
# # Generate log-normal distribution
# np.random.seed(42)  # For reproducibility
# data = np.random.gamma(4, 0.4, N)
# # Randomly select 20 indices to highlight
# highlight_indices = np.random.choice(len(data), 10, replace=False)
# highlight_indices = np.hstack([highlight_indices,int(N*0.9), int(N*0.95), int(N*0.99), int(N*0.99)+7])
# data.sort()
# highlight_values = data[highlight_indices]
# print(highlight_indices)
# print(highlight_values)

# # Create histogram
# plt.figure(figsize=(12, 6))

# # Plot histogram of the entire dataset
# counts, bin_edges, patches = plt.hist(data, bins=50, alpha=0.7, color='gray', edgecolor='black', 
# 	# label='Distribution'
# 	)

# # Add highlighted bars
# for i in range(len(bin_edges) - 1):
#     # Find indices of data points in the current bin
#     bin_indices = np.where((data >= bin_edges[i]) & (data < bin_edges[i + 1]))[0]
#     # Highlight values in the chosen bins
#     highlight_bin_indices = np.intersect1d(bin_indices, highlight_indices)
#     if len(highlight_bin_indices) > 0:
#         # Compute height of the bin to match the histogram height
#         bin_height = counts[i]
#         print(bin_height)
#         if bin_height < 50:
#         	bin_height = 50
#         for value in highlight_values:
#             if bin_edges[i] <= value < bin_edges[i + 1]:
#                 # pass
#                 print(value)
#                 plt.vlines(value, 0, bin_height, color='red', linewidth=2.5)
#             else:
#             	# pass
#                 print(value)
#                 # bin_height_to_use = 40
#                 # plt.vlines(value, 0, bin_height_to_use, color='red', linewidth=2.5)

# # Add labels and title
# plt.xlabel('Value')
# plt.ylabel('Frequency')
# # plt.title('Histogram of Log-Normal Distribution with Highlighted Values')
# plt.legend()
# plt.grid(True)

# # Show the plot
# plt.show()


################################################################ reaction list
# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt

# # Creating a DataFrame from your data
# data = {
#     'Reaction ID': ['MAR04388_f', 'MAR08514', 'MAR06606_f', 'MAR08609', 'MAR08611', 'MAR00710_f', 'MAR03957', 
#                       'MAR04145', 'MAR04152_f', 'MAR04209', 'MAR04410_f', 'MAR04589_f', 'MAR06413', 'MAR06414_f', 
#                       'MAR03977', 'MAR06911_f', 'MAR06914', 'MAR06918', 'MAR00357'],
# }

# df = pd.DataFrame(data)

# styled_df = df.style.set_table_styles([
#     {'selector': 'thead th', 'props': [('background-color', '#4CAF50'), ('color', 'white'), ('font-weight', 'bold')]},
#     {'selector': 'tbody td', 'props': [('border', '1px solid black')]}
# ]).set_properties(**{'text-align': 'center'}).hide(axis = "index")
# styled_df
# styled_df.to_html("styled_table_only_reactions.html")

################################################################ sample jaccards
# import numpy as np
# import matplotlib.pyplot as plt
# from sklearn.metrics import jaccard_score
# from itertools import combinations, permutations
# import time
# from scipy.cluster.hierarchy import linkage, leaves_list
# import os
# from scipy.sparse import lil_matrix
# import numpy as np
# from sklearn.metrics import jaccard_score
# from itertools import combinations
# from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
# import matplotlib.pyplot as plt
# import seaborn as sns
# from sklearn.cluster import SpectralClustering, KMeans
# from concurrent.futures import ThreadPoolExecutor
# from scipy.sparse import lil_matrix
# from matplotlib.colors import LinearSegmentedColormap


# def create_jaccard_sample_x_sample_matrix(samples):
#     n_samples = len(samples)
#     jaccard_matrix = lil_matrix((n_samples, n_samples))
    
#     def jaccard_for_pair(i, j):
#         # Create a union of the sets from samples[i] and samples[j]
#         union = samples[i].union(samples[j])
#         # Create binary vectors relative to the union
#         vector_i = np.array([1 if item in samples[i] else 0 for item in union])
#         vector_j = np.array([1 if item in samples[j] else 0 for item in union])
#         # Compute the Jaccard index for the pair
#         return jaccard_score(vector_i, vector_j, average='binary')
    
#     with ThreadPoolExecutor() as executor:
#         futures = []
#         for i in range(n_samples):
#             for j in range(i, n_samples):
#                 futures.append(executor.submit(jaccard_for_pair, i, j))
        
#         idx = 0
#         for i in range(n_samples):
#             for j in range(i, n_samples):
#                 jaccard_matrix[i, j] = futures[idx].result()
#                 jaccard_matrix[j, i] = jaccard_matrix[i, j]
#                 idx += 1

#     return jaccard_matrix.toarray()

# def optimize_order_hierarchical(jaccard_matrix):
#     Z = linkage(jaccard_matrix, 'average')
#     optimal_order = leaves_list(Z)
#     optimized_jaccard_matrix = jaccard_matrix[optimal_order, :][:, optimal_order]
#     return optimized_jaccard_matrix, optimal_order, Z

# def optimize_order_greedy(jaccard_matrix):
#     n_samples = len(jaccard_matrix)
#     unvisited = set(range(n_samples))
#     current_index = 0
#     optimal_order = [current_index]
#     unvisited.remove(current_index)

#     while unvisited:
#         # Select the next index with the maximum similarity (Maximum Jaccard index)
#         next_index = max(unvisited, key=lambda x: jaccard_matrix[current_index, x])
#         optimal_order.append(next_index)
#         unvisited.remove(next_index)
#         current_index = next_index

#     optimized_jaccard_matrix = jaccard_matrix[np.ix_(optimal_order, optimal_order)]
#     return optimized_jaccard_matrix, optimal_order

# def optimize_order_sklearn(jaccard_matrix, n_clusters=2):
#     clustering = SpectralClustering(n_clusters=n_clusters, affinity='precomputed', assign_labels='discretize', random_state=0)
#     labels = clustering.fit_predict(1 - jaccard_matrix)
#     sorted_indices = np.argsort(labels)
#     optimized_jaccard_matrix = jaccard_matrix[sorted_indices, :][:, sorted_indices]
#     return optimized_jaccard_matrix, sorted_indices

# def optimize_order_kmeans(jaccard_matrix, n_clusters=2):
#     clustering = KMeans(n_clusters=n_clusters, random_state=0)
#     labels = clustering.fit_predict(jaccard_matrix)
#     sorted_indices = np.argsort(labels)
#     optimized_jaccard_matrix = jaccard_matrix[sorted_indices, :][:, sorted_indices]
#     return optimized_jaccard_matrix, sorted_indices

# def visualize_optimized_jaccard_matrix(optimized_jaccard_matrix, title, cmap='viridis'):
#     # plt.figure(figsize=(10, 10))
#     n_bins = 500  # Number of bins in the colormap.
#     colors = ['orange','blue', "purple"]
#     custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=n_bins)
#     # sns.heatmap(optimized_jaccard_matrix, cmap=custom_cmap, annot=False, cbar=True,
#     # 	row_cluster = True, metric='euclidean')
#     sns.clustermap(optimized_jaccard_matrix, cmap=custom_cmap, center=0.8, method='average', metric='euclidean', 
#                row_cluster=True, col_cluster=False, figsize=(10, 12))
#     # plt.title(title)
#     # plt.xlabel('Sample Index (Reordered)')
#     # plt.ylabel('Sample Index (Reordered)')
#     # plt.show()

# def visualize_dendrogram(Z, labels):
#     plt.figure(figsize=(12, 8))
#     dendrogram(Z, labels=[f'Sample {i}' for i in labels], leaf_rotation=90, leaf_font_size=8)
#     plt.title('Dendrogram of Hierarchical Clustering')
#     plt.xlabel('Sample Index')
#     plt.ylabel('Distance')
#     plt.show()

# import numpy as np
# import pandas as pd
# import random
# # Seed for reproducibility
# np.random.seed(42)



# # Helper function to generate correlated integers

# def generate_correlated_samples(num_samples, min_unique_count, max_unique_count, correlation, common_integers):
#     """
#     Generate samples with a given correlation and varying unique counts.
    
#     Parameters:
#     - num_samples: Number of samples to generate.
#     - min_unique_count: Minimum number of unique integers per sample.
#     - max_unique_count: Maximum number of unique integers per sample.
#     - correlation: Correlation factor (between 0 and 1).
    
#     Returns:
#     - A list of samples, each containing correlated values.
#     """
#     currently_included = []  # To store the already included values
#     sample_values = [[] for _ in range(num_samples)]  # Initialize a list of lists for sample values
#     unique_count = [np.random.randint(min_unique_count, max_unique_count + 1) for _ in range(num_samples)]
    
#     not_all_done = True  # Flag to check if all samples are filled
    
#     while not_all_done:
#         for sample_iterator in range(num_samples):
#             # Check if the current sample still needs more unique integers
#             if len(sample_values[sample_iterator]) < unique_count[sample_iterator]:
#                 random_number = random.random()
                
#                 # Use correlation to either select from currently included or generate new
#                 if random_number < correlation and len(currently_included) > 0:
#                     # Select a value from already included values to maintain correlation
#                     value = np.random.choice(currently_included, 1, replace=False)[0]
#                 else:
#                     # Generate a new unique value outside the range of common values (e.g., from 31 to 500)
#                     value = np.random.choice(np.arange(31, 501), 1, replace=False)[0]
                
#                 # Ensure value is not duplicated in the current sample
#                 if value not in sample_values[sample_iterator]:
#                     sample_values[sample_iterator].append(value)
                    
#                     # Add the value to currently_included and keep it unique
#                     currently_included.append(value)
#                     currently_included = list(set(currently_included))

#         # Check if all samples have reached their required unique count
#         not_all_done = any(len(sample_values[i]) < unique_count[i] for i in range(num_samples))
    
#     for _ in sample_values:
#     	_.extend(common_integers)
#     return sample_values



# # Parameters
# n_samples = 322
# common_count = 50
# min_unique_count = 10
# max_unique_count = 20
# correlation_target1 = 0.995
# correlation_target2 = 0.7
# n_samples_group1 = 161
# n_samples_group2 = 161

# # Generate common integers
# common_integers = np.arange(1, common_count + 1)

# samples_group1 = generate_correlated_samples(
# 	n_samples,
# 	min_unique_count,
# 	max_unique_count,
# 	correlation_target1,
# 	common_integers
# 	)

# # print(samples_group1)

# #####
# # small_samples = [samples_group1[0], samples_group1[1]]
# # for list_ in small_samples:
# # 	for element in list_:
# # 		element = str(element)
# # small_samples = [set(list_)for list_ in small_samples]
# # n_samples = len(small_samples)
# # print(small_samples)
# # jaccard_matrix = lil_matrix((n_samples, n_samples))
# # def jaccard_for_pair_edit(i, j, average):
# #     # Create a union of the sets from samples[i] and samples[j]
# #     union = small_samples[i].union(small_samples[j])
# #     # Create binary vectors relative to the union
# #     vector_i = np.array([1 if item in small_samples[i] else 0 for item in union])
# #     vector_j = np.array([1 if item in small_samples[j] else 0 for item in union])
# #     # Compute the Jaccard index for the pair
# #     return jaccard_score(vector_i, vector_j, average=average)


# # print(jaccard_for_pair_edit(0, 1, "micro"))
# # print(jaccard_for_pair_edit(0, 1, "binary"))
# #####

# samples_group2 = samples_group1 = generate_correlated_samples(
# 	n_samples,
# 	min_unique_count,
# 	max_unique_count,
# 	correlation_target2,
# 	common_integers
# 	)

# samples_group1.extend(samples_group2)
# samples = samples_group1

# #samples = [list_.tolist() for list_ in samples]
# for list_ in samples:
# 	for element in list_:
# 		element = str(element)

# samples = [set(list_)for list_ in samples]
# print(samples)
# print("\nCreating Jaccard matrix...")
# jaccard_matrix = create_jaccard_sample_x_sample_matrix(samples)
# print("\nJaccard matrix created.")

# # Hierarchical method
# print("\nOptimizing with Hierarchical Clustering...")
# optimized_matrix_hierarchical, optimal_order_hierarchical, Z = optimize_order_hierarchical(jaccard_matrix)
# visualize_optimized_jaccard_matrix(optimized_matrix_hierarchical, 'Optimized Jaccard Matrix (Hierarchical)')
# visualize_dendrogram(Z, optimal_order_hierarchical)

# # optimized_matrix_kmeans, optimal_order_kmeans = optimize_order_kmeans(jaccard_matrix, n_clusters=2)
# # visualize_optimized_jaccard_matrix(optimized_matrix_kmeans, 'Optimized Jaccard Matrix (KMeans Clustering)')

##################################################### horizontal bar plots
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib import cm
# from matplotlib.patches import Patch
# import random
# random.seed(10)

# generate_specific_values =True
# amount_of_groups = 3
# amount_to_show = 3

# # 1. Generate 322 random values between 0 and 1
# np.random.seed(42)
# values = np.random.rand(122)


# # 2. Sort the values and select the top 30
# top_values = np.sort(values)[-amount_to_show:]

# task_numbers = np.random.choice(122, amount_to_show, replace=False)
# # group_names = [
# # "Aerobic Respiration",
# # "Anaerobic Respiration",
# # "Fatty Acid Biosynthesis",
# # "Amino Acid Degradation",
# # "Amino Acid Biosynthesis",
# # "Lactose Degradation",
# # "Ammonium Processes",
# # "Primary Alcohol Metabolism",
# # "Drug Catabolic Processes",
# # "Spermidine Biosynthesis",
# # "Cofactor metabolic Processes",
# # "Ion Transfer Processes",
# # "Organic Hydroxy Compound Processes"
# #  ]
# group_names = ["Healthy", "Diseased", "Diseased + Compensating"]
# groups = [f"{group_names[(i%len(group_names))]}" for i in task_numbers]
# if generate_specific_values:
# 	amount_of_tasks_to_show = amount_of_groups//amount_to_show
# 	np.random.seed(42)
# 	values = []
# 	values.extend(np.random.normal(12,1,1))
# 	values.extend(np.random.normal(10.5,1,1))
# 	values.extend(np.random.normal(4,1,1))
# 	# values = np.vstack([values1, values2, values3]
# 	top_values = np.sort(values)
# 	group_names = [ "Diseased", "Diseased + Compensating","Healthy",]
# 	groups = [f"{group_names[(i%len(group_names))]}" for i in range((amount_of_groups))]

# print(top_values, groups)
# # 1. Identify unique groups
# unique_groups = list(set(groups))

# # 2. Assign colors evenly using a colormap
# # Using the "tab20" colormap for up to 20 distinct colors
# colormap = cm.get_cmap('tab20', len(unique_groups))

# orange_purple_blue = [ "LightBlue","Orange", "Purple",]

# # Create a mapping of group to color
# group_to_color = {group: colormap(i) for i, group in enumerate(unique_groups)}

# group_to_color = {group: orange_purple_blue[i] for i, group in enumerate(unique_groups)}
# # 3. Map colors to groups in the same order as the 'groups' list
# colors = [group_to_color[group] for group in groups]


# # 4. Create a horizontal bar chart
# plt.figure(figsize=(10, 6))

# # Create the horizontal bar chart
# plt.barh(range(amount_to_show), top_values, color=colors)

# np.random.seed(42)
# task_numbers = np.random.choice(122, amount_to_show, replace=False)
# # Add labels and title
# plt.yticks(range(amount_to_show), [f'Sample {i+1}' for i in task_numbers])
# # plt.xlabel('Tasks')
# plt.title('Task score')
# print(unique_groups)
# legend_patches = [Patch(color=group_to_color[group], label=group) for group in unique_groups]
# plt.legend(handles=legend_patches, title="Metabolic function", bbox_to_anchor=(1.05, 1), loc='upper left')

# plt.tight_layout()  # Adjust layout to make space for the legend
# plt.show()

# #################################################### Network graph
# import networkx as nx
# import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches
# import matplotlib.patheffects as PathEffects

# # Create a directed graph
# G = nx.DiGraph()

# # Path 1 (red path)
# path1 = ["Input", "A", "B", "C", "D", "E","F", "Output"]

# # Path 2 (green path)
# path2 = ["Input", "A",  "G", "H", "I","J", "Output"]

# # Add nodes and edges for both paths
# for i in range(len(path1) - 1):
#     G.add_edge(path1[i], path1[i + 1])

# for i in range(len(path2) - 1):
#     G.add_edge(path2[i], path2[i + 1])

# # Set the positions of the nodes to reflect a top-to-bottom layout
# # Manually adjust positions to space the nodes clearly
# pos = {
#     "Input": (0, 6),
#     "A": (0, 5),
#     "B": (-1, 5),  # Overlapping node
#     "C": (-1, 4),
#     "D": (-1, 3),
#     "E": (-1, 2),
#     "F": (-1, 1),
#     "G": (1, 4),
#     "H": (1, 3),
#     "I": (1, 2),
#     "J": (1, 1),
#     "Output": (0, 0)
# }

# # Draw the nodes and edges for path1 (in red)
# nx.draw_networkx_edges(G, pos, edgelist=[(path1[i], path1[i+1]) for i in range(len(path1)-1)], edge_color="red", width=2)

# # Draw the nodes and edges for path2 (in green)
# nx.draw_networkx_edges(G, pos, edgelist=[(path2[i], path2[i+1]) for i in range(len(path2)-1)], edge_color="green", width=2)

# # Draw the nodes
# unique_nodes = set(path1 + path2)
# non_overlapping_nodes = unique_nodes - set(path1).intersection(set(path2))
# non_overlap_path1 = {node for node in path1 if node in non_overlapping_nodes}
# non_overlap_path2 = {node for node in path2 if node in non_overlapping_nodes}
# nx.draw_networkx_nodes(G, pos, nodelist=set(non_overlap_path1), node_color="red", edgecolors="black", node_size=1000)
# nx.draw_networkx_nodes(G, pos, nodelist=set(non_overlap_path2), node_color="green", edgecolors="black", node_size=1000)
# # Draw labels
# nx.draw_networkx_labels(G, pos, font_size=12)

# overlapping_nodes = set(path1).intersection(set(path2))

# for node in overlapping_nodes:
#     node_pos = pos[node]
    
#     # Draw the red half
#     red_half = plt.Circle(node_pos, radius=0.15, color='red', transform=plt.gca().transData, zorder=2, clip_on=False)
#     plt.gca().add_patch(red_half)
    
#     # Draw the green half, slightly overlapping the red half
#     green_half = plt.Circle((node_pos[0] + 0.075, node_pos[1]), radius=0.15, color='green', transform=plt.gca().transData, zorder=2, clip_on=False)
#     plt.gca().add_patch(green_half)

# # Set aspect ratio to equal to ensure circles are not distorted
# plt.gca().set_aspect('equal', adjustable='datalim')

# # Create legend
# red_patch = mpatches.Patch(color='red', label='Sample 1')
# green_patch = mpatches.Patch(color='green', label='Sample 2')
# plt.legend(handles=[red_patch, green_patch])

# # Display the plot
# plt.title("Two Paths with Overlapping Nodes")
# plt.axis('off')
# plt.show()
##################################################### Overlapping histograms

# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns

# # Generate two distributions of 300 values each
# np.random.seed(42)  # For reproducibility
# dist1 = np.random.normal(loc=0, scale=1, size=1000)  # Mean = 0, Std = 1
# dist2 = np.random.normal(loc=0.5, scale=1, size=1000)  # Mean = 0.5, Std = 1 (shifted)


# # Plot the KDE (smoothed distributions)
# plt.figure(figsize=(10, 6))
# sns.kdeplot(dist1, color='red', shade=True, label='Sample 1')
# sns.kdeplot(dist2, color='blue', shade=True, label='Sample 2')

# # Add labels and title
# plt.xlabel('Value')
# plt.ylabel('Density')

# # Add a legend
# plt.legend()

# # Display the plot
# plt.show()


######################### Bar plot of NOde centrality
# import matplotlib.pyplot as plt
# import networkx as nx
# import matplotlib.pyplot as plt

# # Create a graph
# G = nx.erdos_renyi_graph(20, 0.3)

# # Compute centrality
# centrality = nx.degree_centrality(G)

# # Draw the network with node sizes proportional to centrality
# pos = nx.spring_layout(G)
# plt.figure(figsize=(10, 8))
# nx.draw_networkx_nodes(G, pos, node_size=[v * 1000 for v in centrality.values()], node_color='lightblue', edgecolors='black')
# nx.draw_networkx_edges(G, pos)
# nx.draw_networkx_labels(G, pos)
# plt.title('Network with Node Sizes Based on Degree Centrality')
# plt.show()


# import seaborn as sns
# import pandas as pd

# # Compute centrality
# centrality = nx.degree_centrality(G)

# # Convert to DataFrame
# centrality_df = pd.DataFrame(list(centrality.items()), columns=['Node', 'Centrality'])

# # Create heatmap
# plt.figure(figsize=(10, 6))
# sns.heatmap(centrality_df.set_index('Node').T, cmap='YlGnBu', annot=True)
# plt.title('Heatmap of Node Centralities')
# plt.show()

# influential_nodes = sorted(centrality, key=centrality.get, reverse=True)[:5]  # Top 5 nodes

# plt.figure(figsize=(10, 8))
# nx.draw_networkx_nodes(G, pos, node_size=500, node_color='lightgrey', edgecolors='black')
# nx.draw_networkx_nodes(G, pos, nodelist=influential_nodes, node_size=1000, node_color='red')
# nx.draw_networkx_edges(G, pos)
# nx.draw_networkx_labels(G, pos)
# plt.title('Network with Influential Nodes Highlighted')
# plt.show()
# import matplotlib.pyplot as plt

# # Compute centrality
# centrality = nx.degree_centrality(G)

# # Plot bar chart
# plt.figure(figsize=(10, 6))
# nodes = list(centrality.keys())
# values = list(centrality.values())
# plt.bar(nodes, values, color='blue')
# plt.xlabel('Node')
# plt.ylabel('Centrality')
# plt.title('Bar Plot of Node Centrality')
# plt.xticks(rotation=90)
# plt.show()


# # ############################################################### Ridgeline plot
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# import random
# from matplotlib import cm

# amount_of_samples = 10
# amount_of_chosen_samples = 10
# # Step 1: Generate 250 random centers, most between 5 and 25
# np.random.seed(42)
# centers = np.random.normal(loc=5, scale=10, size=amount_of_samples)
# centers = np.clip(centers, 0.5, 250)  # Ensure centers stay within the range of 5 and 25

# # Step 2: Generate gamma distributions surrounding the centers
# distributions = []
# seedparam = 100


# random_choice_ = 0.0
# for center in centers:    
#     shape_param = center   # Shape parameter related to the center
#     scale_param = (max(np.random.normal(loc=3, scale=3, size=1)[0],0.2))
#     # Scale parameter can be fixed or adjusted
#     np.random.seed(seedparam)
#     # random_choice_ = random.random()
#     if random_choice_ < 0.2:
#     	  distribution = np.random.poisson(shape_param, 100) +9.5
#     elif random_choice_ < 0.4:        
#         X1 = np.random.normal(center + 4, 0.1, 500)
#         X2 = np.random.normal(center + 20, 0.1, 500)
#         distribution = np.concatenate([X1, X2])
#     elif random_choice_ < 0.6:
#     	distribution = np.random.gamma(shape_param+9.5, scale_param, 100) +10
#     elif random_choice_ < 0.8:
#     	distribution = np.random.gumbel(shape_param+9.5, scale_param, 100)
#     else:
#     	distribution = np.random.normal(5, 0.1, 100) +10 

#     distributions.append(distribution)
#     seedparam +=1
#     random_choice_ += 0.1
#     random_choice_ = random_choice_%1

# group_names = [
# "Aerobic Respiration",
# "Anaerobic Respiration",
# "Fatty Acid Biosynthesis",
# "Amino Acid Degradation",
# "Amino Acid Biosynthesis",
# "Lactose Degradation",
# "Ammonium Processes",
# "Primary Alcohol Metabolism",
# "Drug Catabolic Processes",
# "Spermidine Biosynthesis",
# "Cofactor metabolic Processes",
# "Ion Transfer Processes",
# "Organic Hydroxy Compound Processes"
#  ]

# # Step 3: Create groups and assign colors
# task_numbers = list(range(1,amount_of_samples+1))
# groups = [group_names[(i // (amount_of_samples//))] for i in task_numbers]  # Assign groups


# # Step 4: Prepare the DataFrame for plotting
# data = []
# for i, dist in enumerate(distributions):
#     data.append(pd.DataFrame({'Value': dist, 'Task': f'Task {i+1}', 'Group': groups[i]}))

# df = pd.concat(data).reset_index()
# # randomly take 25
# random_choices = np.random.choice(amount_of_samples, min(amount_of_samples, amount_of_chosen_samples), replace = False)
# chosen_tasks = [f'Task {i+1}'for i in random_choices]
# df_sampled = df[df["Task"].isin(chosen_tasks)]
# from ridgeplot import ridgeplot
# grouped_samples = []
# for Task in df_sampled["Task"].unique():    
#     values = df_sampled[df_sampled["Task"] == Task]["Value"].to_numpy()
#     grouped_samples.append(values)

# # print(len(grouped_samples))
# # Define unique groups
# unique_groups = df["Group"].unique()

# # Create a colormap and map groups to colors
# colormap = cm.get_cmap('Dark2' , len(unique_groups))

# group_to_color = {group: colormap(i) for i, group in enumerate(unique_groups)}

# # Map colors to tasks based on their group
# df_sampled['Color'] = df_sampled['Group'].map(group_to_color)

# # print(df_sampled['Color'])
# # Now `grouped_samples` contains arrays, each representing a distribution for each group
# # Create the ridge plot
# fig = ridgeplot(
#     samples=grouped_samples,            # Pass the list of arrays
#     bandwidth=4,
#     kde_points=np.linspace(-12.5, 112.5, 500),
#     coloralpha=0.65,
#     labels= df_sampled["Task"].unique(),  # Ensure the labels match the groups
#     linewidth=2,
#     spacing=5 / 9,
# )
# plotly_fig = fig
# # print(len(df_sampled["Task"].unique()))
# # print(len(plotly_fig.data))
# groups_in_sample = unique_groups = df_sampled["Group"].unique()
# legend_dict = {group: idx+1 for idx, group in enumerate(groups_in_sample)}
# show_legends = {group: False for idx, group in enumerate(groups_in_sample)}
# counter_i = 0

# for i, trace in enumerate(plotly_fig.data):
#     if i%2 ==0: 
#     	counter_i +=1
#     	continue
#     i = i - counter_i

#     task_name = df_sampled["Task"].unique()[i]
    
#     # Get the group number for this task
#     group_number = df_sampled[df_sampled["Task"] == task_name]["Group"].iloc[0]

#     # Get the RGB color value based on the group number
#     color_values = group_to_color[group_number]#df_sampled[df_sampled['Group'] == group_number]['Color'].iloc[0]
#     rgb_colors = f"rgb({int(color_values[0]*255)}, {int(color_values[1]*255)}, {int(color_values[2]*255)})"
#     rgba_colors = f"rgba({int(color_values[0]*255)}, {int(color_values[1]*255)}, {int(color_values[2]*255)}, 0.65)"
#     # Set the color for each trace
#     trace.line.color = rgb_colors  
    
#     trace.fillcolor = rgba_colors  # Color for the filled area
#     trace.legendgroup = group_number  # Assign traces to a legend group

#     # Show legend for only one trace per group
#     if trace.name == task_name:
#         trace.name = f"Group {group_number}"  # Set the name for the legend
#         if not show_legends[group_number]:
#             show_legends[group_number] = True
#             trace.showlegend = True  # Show legend for this trace
#         else:
#             trace.showlegend = False
#     else:
#         trace.showlegend = False  # Hide legend for other traces


# # fig.for_each_trace(lambda t: t.update(name = newnames[t.name],
# #                                       legendgroup = newnames[t.name],
# #                                       hovertemplate = t.hovertemplate.replace(t.name, newnames[t.name])
# #                                      )
# #                   )
# # plotly_fig.update_layout(legend=dict(
# #     orientation="v",
# #     # yanchor="bottom",
# #     # y=1.02,
# #     # xanchor="right",
# #     # x=1
# # ))
# plotly_fig.update_xaxes(range= [0,45])
# plotly_fig.update_layout(
#     plot_bgcolor='white'
# )
# plotly_fig.update_xaxes(
#     mirror=False,
#     ticks='outside',
#     showline=True,
#     linecolor='black',
#     gridcolor='lightgrey'
# )
# plotly_fig.update_yaxes(
#     mirror=False,
#     ticks='outside',
#     showline=True,
#     linecolor='black',
#     gridcolor='lightgrey'
# )
# plotly_fig.show()

###################################################### 3d parameter plot
import dash
from dash import dcc, html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
import os
import gurobi_logtools as glt_gurb
from io import StringIO
import json
import math
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
import numpy as np
import re
from sklearn.datasets import make_classification
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
import numpy as np
from sklearn import svm
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from cobra.io import save_json_model, load_json_model
from cobra import Model
from os import path, makedirs, listdir
from typing import List, Dict, Tuple
import matplotlib.pyplot as plt
import seaborn as sns



def check_SV_constraints_and_model_feasilibity(
        cobra_model: Model,
):
    solution = cobra_model.optimize()
    if solution.status == "optimal":
        print("Model is feasible")
        return True
    else:
        print("Model is infeasible")
        return False



def check_all_json_models_in_dir_for_feasibility(
        dir_with_models: str,
) -> Dict[str, bool]:
    is_feasible = {}
    for file in listdir(dir_with_models):
        if file.endswith(".json") and file.startswith("model_json_"):
            try:
                filename = file.split(".json")[0]
                filename = filename.split("model_json_")[1]
                model = load_json_model(path.join(dir_with_models, file))
                is_feasible[filename] = check_SV_constraints_and_model_feasilibity(model)
                print(f"Model {filename} checked")
            except:
                print(f"Model {file} could not be loaded")
                continue
    return is_feasible


from datetime import datetime
def find_closest_log_file(filename):
    files_dir = os.getcwd()
    log_file_name = f"{filename}_gurobi.log"
    filtered_time = os.path.getmtime(log_file_name)
    
    # Initialize variables to keep track of the closest log file and time difference
    closest_excel_file = None
    smallest_time_diff = float('inf')
    
    # Iterate over all .log files in the specified directory
    for excel_file in os.listdir(files_dir):
        if excel_file.endswith(".xlsx") and excel_file.startswith("filtered_"):
            excel_file_path = os.path.join(files_dir, excel_file)
            excel_file_time = os.path.getmtime(excel_file_path)
            
            # Calculate the time difference
            time_diff = abs(excel_file_time - filtered_time)
            
            # Update the closest log file if this one is closer in time
            if time_diff < smallest_time_diff:
                closest_excel_file = excel_file
                smallest_time_diff = time_diff
    
    print(closest_excel_file)
    return closest_excel_file

    
def read_xlsx_to_identify_integrality_violations(filename
    ) -> bool:
    try:
        file = f"filtered_{filename}.xlsx"
        df = pd.read_excel(file)
        print(df.iloc[:,7].unique())
        print(set(df.iloc[:,7].unique()).issubset({0,1}))
        if set(df.iloc[:,7].unique()).issubset({0,1}):
            return True
        else:
            return False
    except:
        print(f"could not find {file}")
        filename_2 = find_closest_log_file(filename)
        #.split("_gurobi.log")[0]
        #file_2 = f"filtered_{filename_2}.xlsx"
        try:
            df = pd.read_excel(filename_2)
            print(df.iloc[:,7].unique())
            print(set(df.iloc[:,7].unique()).issubset({0,1}))
            if set(df.iloc[:,7].unique()).issubset({0,1}):
                return True
            else:
                return False
        except:
            print(f'Also could not find {filename_2}')
            return None
        return None


def read_xlsx_to_identify_number_of_reactions(filename
    ) -> bool:
    try:
        file = f"filtered_{filename}.xlsx"
        df = pd.read_excel(file)
        return len(df.index)
    except:
        print(f"could not find {file}")
        filename_2 = find_closest_log_file(filename)
        try:
            df = pd.read_excel(filename_2)
            return len(df.index)
        except:
            print(f'Also could not find {filename_2}')
            return None
        return None
    


def run_plotly_dash(path, z_is_column = "Runtime", df_name = "*"):

    os.chdir(path)
    # print(os.getcwd())
    summary, timelines = glt_gurb.get_dataframe([f"{df_name}.log"], timelines=True)
    df = summary
    # print(df)
    is_feasible_dict = check_all_json_models_in_dir_for_feasibility(path)


    df = df[["Runtime",'FeasibilityTol (Parameter)','IntFeasTol (Parameter)', "Status", "LogFilePath", "ObjVal"]]
    df['isFeasible'] = df['LogFilePath'].apply(lambda x: is_feasible_dict.get(x.split("_gurobi.log")[0], False))
    df['notViolatedIntegrality'] = df['LogFilePath'].apply(lambda x: read_xlsx_to_identify_integrality_violations(x.split("_gurobi.log")[0]))
    df["numberReactions"] = df['LogFilePath'].apply(lambda x: read_xlsx_to_identify_number_of_reactions(x.split("_gurobi.log")[0]))
    print(df)
    # print(df['LogFilePath'])
    # print(is_feasible_dict.keys())
    # print(df.at[5, 'LogFilePath'].split("_gurobi.log")[0]) 
    # print(df["isFeasible"])
    gradient_column = "Epsilon"
    if z_is_column == "Epsilon":
        gradient_column = "Runtime"

    def extract_epsilon(text):
    # Find all the numbers that follow INIT_irrev
        matches = re.findall(r'(1e-?[0-9]?[0-9])', text)
        # print(matches)
        if len(matches) >= 3:
            # If there are 3 or more occurrences, take the first one
            return matches[0]
        elif len(matches) == 2:
            # If there are 2 occurrences, take the portion before the first 1e-
            match = re.search(r'INIT_irrev([0-9]+(?:\.[0-9]+)?)1e-', text)
            if match:
                return match.group(1)
        return None

    def is_correct(objv):
        if objv <= 600:
            return False
        return True
    def clip_runtime(time):
        if time>=3000:
            return 3000
        return time

    def apply_log(text, base= 10):
        if text is not None:
            text = float(text)
            return math.log(text,base)
        return None
    def get_value(value):
        return float(value)

    df["Epsilon"] = df["LogFilePath"].apply(extract_epsilon)
    df["Correct"] = df["ObjVal"].apply(is_correct)
    df["Runtime"] = df["Runtime"].apply(clip_runtime)
    df["FeasibilityTol (Parameter)"] = df["FeasibilityTol (Parameter)"].apply(apply_log)
    df["IntFeasTol (Parameter)"] = df["IntFeasTol (Parameter)"].apply(apply_log)
    df = df.dropna()
    epsilon = df["Epsilon"].apply(get_value)
    runtime = df["Runtime"].apply(get_value)

    # with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    #     print(df)

    if gradient_column == "Epsilon":
        gradient_color_col = df["Epsilon"].apply(apply_log) 
    else:
        gradient_color_col = df["Runtime"].apply(get_value)
        df["Epsilon"] = df["Epsilon"].apply(apply_log)
    # print(df["Epsilon"].apply(get_value) )
    # print(df["Runtime"].apply(get_value) )
    X = df[['FeasibilityTol (Parameter)', 'IntFeasTol (Parameter)', z_is_column]]
    Y = df['Correct']
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    svc = SVC(kernel='linear')
    svc.fit(X,Y)
    x_vals = np.linspace(df["FeasibilityTol (Parameter)"].min(), df["FeasibilityTol (Parameter)"].max(), 100)
    y_vals = np.linspace(df["IntFeasTol (Parameter)"].min(), df["IntFeasTol (Parameter)"].max(), 100)
    x, y = np.meshgrid(x_vals, y_vals)
    z = lambda x,y: (-svc.intercept_[0]-svc.coef_[0][0]*x-svc.coef_[0][1]*y) / svc.coef_[0][2]



   


    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(16, 7))

    # First subplot: Scatter plot of ObjVal by isFeasible
    for label, group in df.groupby('isFeasible'):
        color = 'red' if label == False else 'green'
        axes[0].scatter(group.index, group['ObjVal'], label=f'isFeasible={label}', color=color)
    axes[0].set_xlabel('Index')
    axes[0].set_ylabel('ObjVal')
    axes[0].set_title('ObjVal by isFeasible')
    axes[0].legend()

    # Second subplot: Scatter plot of ObjVal by notViolatedIntegrality
    for label, group in df.groupby('notViolatedIntegrality'):
        color = 'red' if label == False else 'green'
        axes[1].scatter(group.index, group['ObjVal'], label=f'notViolatedIntegrality={label}', color=color)
    axes[1].set_xlabel('Index')
    axes[1].set_ylabel('ObjVal')
    axes[1].set_title('ObjVal by notViolatedIntegrality')
    axes[1].legend()

    df['NoViolation_Feasible'] = df['notViolatedIntegrality'] & df['isFeasible']

      # Second subplot: Scatter plot of ObjVal by notViolatedIntegrality
    for label, group in df.groupby('NoViolation_Feasible'):
        color = 'red' if label == False else 'green'
        axes[2].scatter(group.index, group['ObjVal'], label=f'NoViolation_Feasible={label}', color=color)
    axes[2].set_xlabel('Index')
    axes[2].set_ylabel('ObjVal')
    axes[2].set_title('ObjVal by NoViolation_Feasible')
    axes[2].legend()
    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Show the combined plot
    plt.show()

    # Create subplots: 1 row and 2 columns
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16, 7))

    # Filter the DataFrame to only include rows where NoViolation_Feasible is True
    filtered_df = df[df['NoViolation_Feasible'] == True]
    filtered_df = filtered_df.reset_index(drop=True)
    print(filtered_df["Status"].unique())
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(filtered_df[["Runtime", "Epsilon", "FeasibilityTol (Parameter)", "IntFeasTol (Parameter)", "ObjVal"]].sort_values("Runtime"))
    filtered_only_optimal = filtered_df[filtered_df['Status'] == "OPTIMAL"]
    print(f'Unique objective values = {filtered_only_optimal["ObjVal"].unique()}')

    # Further filtering based on ObjVal
    filtered_only_optimal = filtered_only_optimal[
        (filtered_only_optimal["ObjVal"] >= 606.707) & 
        (filtered_only_optimal["ObjVal"] <= 606.80)
    ]
    print(filtered_only_optimal[["Runtime", "Epsilon", "FeasibilityTol (Parameter)", "IntFeasTol (Parameter)"]])
    filtered_df[["Runtime", "Epsilon", "FeasibilityTol (Parameter)", "IntFeasTol (Parameter)", "ObjVal"]].sort_values("Runtime").to_excel("filtered_data_frame.xlsx")
    # Ensure DataFrame has a reset index
    # Jitter for ObjVal and Runtime
    jitter_strength = 0.0  # Adjust this value to control the amount of jitter
    objval_jitter = filtered_df['ObjVal'] + np.random.uniform(-jitter_strength, jitter_strength, size=filtered_df.shape[0])
    runtime_jitter = filtered_df['Runtime'] + np.random.uniform(-jitter_strength, jitter_strength, size=filtered_df.shape[0])

   
    # Epsilon values and assign categorical colors
    filtered_df['Epsilon'] = pd.to_numeric(filtered_df['Epsilon'], errors='coerce')
    Epsilon_values = filtered_df['Epsilon'].values
    Epsilon_values_groups = filtered_df['Epsilon'].unique()

    # Generate a color palette for the unique Epsilon values
    palette = sns.color_palette("bright", len(Epsilon_values_groups))  # Use any color palette

    # Map each unique Epsilon to a color
    color_map = {eps: color for eps, color in zip(Epsilon_values_groups, palette)}
    colors = [color_map[eps] for eps in Epsilon_values]

    # Create the scatter plot for the first subplot (axes[0])
    sc = axes[0].scatter(objval_jitter, runtime_jitter, c=colors, label='NoViolation_Correct=True')

    # Create a custom legend for the Epsilon groups
    handles = [plt.Line2D([0], [0], marker='o', color=color, linestyle='', label=f'Epsilon={eps}')
               for eps, color in color_map.items()]
    axes[0].legend(handles=handles, title='Epsilon Groups')

    # Set labels and title for the first subplot (axes[0])
    axes[0].set_xlabel('ObjVal')
    axes[0].set_ylabel('Runtime')
    axes[0].set_title('Scatter Plot of Runtime vs ObjVal (NoViolation_Feasible=True)')

    # Annotate points with their index values and draw lines
    for i in range(len(filtered_df)):
        # Draw a line from the point to the annotation
        # i = 0
        axes[0].annotate(
            i,  # The index value
            (objval_jitter[i], runtime_jitter[i]),  # The coordinates of the point
            textcoords="offset points",  # Specify the position of the annotation text
            xytext=(5, 5),  # Offset from the point (so the text isn't on the point)
            arrowprops=dict(arrowstyle="->", lw=0.5),  # Arrow pointing to the point
            fontsize=8  # Font size for the annotation
        )

    # Jitter for ObjVal and numberReactions (for the second subplot)
    jitter_strength = 0.0  # Adjust this value to control the amount of jitter
    objval_jitter = filtered_df['ObjVal'] + np.random.uniform(-jitter_strength, jitter_strength, size=filtered_df.shape[0])
    number_reactions_jitter = filtered_df['numberReactions'] + np.random.uniform(-jitter_strength, jitter_strength, size=filtered_df.shape[0])

    # Create the scatter plot for the second subplot (axes[1])
    # axes[1].scatter(objval_jitter, number_reactions_jitter, color='green', label='NoViolation_Correct=True')

    # Create the scatter plot for the first subplot (axes[0])
    sc = axes[1].scatter(objval_jitter, number_reactions_jitter, c=colors, label='NoViolation_Correct=True')

    # Create a custom legend for the Epsilon groups
    handles = [plt.Line2D([1], [1], marker='o', color=color, linestyle='', label=f'Epsilon={eps}')
               for eps, color in color_map.items()]
    axes[1].legend(handles=handles, title='Epsilon Groups')

    # Annotate points with their index values and draw lines
    for i in range(len(filtered_df)):
        # Draw a line from the point to the annotation
        axes[1].annotate(
            i,  # The index value
            (objval_jitter[i], number_reactions_jitter[i]),  # The coordinates of the point
            textcoords="offset points",  # Specify the position of the annotation text
            xytext=(5, 5),  # Offset from the point (so the text isn't on the point)
            arrowprops=dict(arrowstyle="->", lw=0.5),  # Arrow pointing to the point
            fontsize=8  # Font size for the annotation
        )

    # Set labels and title for the second subplot (axes[1])
    axes[1].set_xlabel('ObjVal')
    axes[1].set_ylabel('numberReactions')
    axes[1].set_title('Scatter Plot of ObjVal vs numberReactions (NoViolation_Feasible=True)')
    axes[1].legend()

    # Show the plot
    plt.show()


    scatter = go.Scatter3d(
        x=df['FeasibilityTol (Parameter)'],
        y=df['IntFeasTol (Parameter)'],
        z=df[f'{z_is_column}'],
        mode='markers',
        marker=dict(
            size=10,  # Increase point size for better visibility
            color=gradient_color_col,
            colorscale=["green","red"],  
            opacity=0.8,
            colorbar=dict(
                title=f'{gradient_column}',
                # tickvals=[0, 0.5, 1],
                # ticktext=['Low', 'Medium', 'High']
            )
        ),
        text=df[f'{gradient_column}'],  # Tooltip shows the Runtime
    )
    highlight = go.Scatter3d(
        x=df.loc[df['Correct'], 'FeasibilityTol (Parameter)'],
        y=df.loc[df['Correct'], 'IntFeasTol (Parameter)'],
        z=df.loc[df['Correct'], f'{z_is_column}'],
        mode='markers',
        marker=dict(
            size=18,  # Increase the size of the red circles
            color='red',
            symbol='circle-open',
            line=dict(width=3, color='red'),  # Thicker red outline
        ),
        showlegend=False
    )
    highlight2 = go.Scatter3d(
        x=df.loc[df['isFeasible'], 'FeasibilityTol (Parameter)'],
        y=df.loc[df['isFeasible'], 'IntFeasTol (Parameter)'],
        z=df.loc[df['isFeasible'], f'{z_is_column}'],
        mode='markers',
        marker=dict(
            size=22,  # Increase the size of the red circles
            color='green',
            symbol='circle-open',
            line=dict(width=3, color='green'),  # Thicker red outline
        ),
        showlegend=False
    )
    highlight3 = go.Scatter3d(
        x=df.loc[df['notViolatedIntegrality'], 'FeasibilityTol (Parameter)'],
        y=df.loc[df['notViolatedIntegrality'], 'IntFeasTol (Parameter)'],
        z=df.loc[df['notViolatedIntegrality'], f'{z_is_column}'],
        mode='markers',
        marker=dict(
            size=26,  # Increase the size of the red circles
            color='yellow',
            symbol='circle-open',
            line=dict(width=3, color='yellow'),  # Thicker red outline
        ),
        showlegend=False
    )
    highlight4 = go.Scatter3d(
        x=df.loc[df['NoViolation_Feasible'], 'FeasibilityTol (Parameter)'],
        y=df.loc[df['NoViolation_Feasible'], 'IntFeasTol (Parameter)'],
        z=df.loc[df['NoViolation_Feasible'], f'{z_is_column}'],
        mode='markers',
        marker=dict(
            size=30,  # Increase the size of the red circles
            color='Purple',
            symbol='circle-open',
            line=dict(width=4, color='Purple'),  # Thicker red outline
        ),
        showlegend=False
    )


    fig = go.Figure(data=[scatter, 
        # highlight,
        # highlight2, 
        # highlight3,
        highlight4,
        ])
    fig.add_surface(x=x, y=y, z=z(x,y), colorscale='Greys', showscale=False)
    min_z = min(df[f'{z_is_column}'])
    max_z = max(df[f'{z_is_column}'])
    if gradient_column == "Runtime":
        fig.update_scenes(zaxis_range=[min_z, max_z])
    fig.update_layout(
        scene=dict(
            xaxis=dict(
                title='FeasibilityTol (Parameter)',
                # type='log',  # Log scale for x-axis
            ),
            yaxis=dict(
                title='IntFeasTol (Parameter)',
                # type='log',  # Log scale for y-axis
            ),
            zaxis=dict(
                title=f'{z_is_column}',
                # type='log',  # Log scale for z-axis
            ),
        ),
    )

    # Show the figure
    fig.show()
    fig.write_html(f"{z_is_column}_to_z.html")

# if __name__ == '__main__':
#     # main_folder = r"E:\Git\Metabolic_Task_Score\Data\Pyomo_python_files\Jelle\test_runs\consensus_NT_LT_INIT_LINEAR_MIN_Tarray_3samples_8tasks"
#     main_folder = r"E:\Git\Metabolic_Task_Score\Data\Pyomo_python_files\Jelle\test_runs\consensus_no_threshold_fractional_options_test_enum"
#     # if not os.getcwd().endswith("results"):
#     #     os.chdir(os.path.join(os.getcwd(), "results"))
#     run_plotly_dash(main_folder)

path = r"E:\Git\Metabolic_Task_Score\Data\Main_files\For_running\combinations\consensus_2024_DCM_test_metabolic_tasks_2024_v1_01"
expression = os.path.join(path, r"NT__NT__minMax__Y__rev__reaction_expression_dataset.feather")
contribution = os.path.join(path, r"NT__NT__minMax__Y__rev__gene_contribution_dataset.feather")
# df1 = pd.read_feather(expression)
df2 = pd.read_feather(contribution)

### add with file pd.read_feather
### add pyarrow.lib.ArrowInvalid file is too small error in the chekcing part of the preprocessing code

# df1.to_csv( os.path.join(path, r"NT__NT__minMax__Y__rev__reaction_expression_dataset.csv"))
df2.to_csv(os.path.join(path, r"NT__NT__minMax__Y__rev__gene_contribution_dataset.csv"))

################################################################### example for troubleshooting some code

##### Some stuff that isn't working

# A = 5
# # B = 6

# def some_function(x):
#     print("this is running now")
#     return str(x) + "hello"

# C = some_function('d')


# if __name__ == "__main__":
#     D = 5
#     

#### Should provide a list with integers
#    so we just test this part of the code by providing a list with integers

# list_with_integers = [2,43,5,12,7,4]

# def function_to_do_things_with_integers(integer_list):
#     for idx in range(len(integer_list)):
#         print(idx)

# function_to_do_things_with_integers(list_with_integers)

