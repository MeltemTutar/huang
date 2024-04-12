
import os
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score, silhouette_samples
import numpy as np

def evaluate_cluster_sil_score(df, row_clusters, col_clusters):
    # Calculate silhouette scores
    row_silhouette_score = silhouette_score(df, row_clusters)
    col_silhouette_score = silhouette_score(df.T, col_clusters)

    return row_silhouette_score, col_silhouette_score


def evaluate_cluster_sil_score_mutated_samples(df, col_clusters, mutated_samples):
    # Calculate silhouette scores
    col_silhouette_scores = silhouette_samples(df.T, col_clusters)

    # Filter the data to keep only mutated samples
    mutated_samples_set = set(mutated_samples)
    mutated_indices = [i for i, sample_id in enumerate(df.columns) if sample_id in mutated_samples_set]

    # Filter col_silhouette_scores for mutated samples
    mutated_silhouette_scores = [col_silhouette_scores[i] for i in mutated_indices]

    # Calculate the average silhouette score for mutated samples
    avg_mutated_silhouette_score = np.mean(mutated_silhouette_scores)

    return avg_mutated_silhouette_score


# TODO: modify this function so that it is for one feature at a time (row or col)
# and call it twice from the main file. 
def hierarchical_clustering(expression_df_heatmap, output_folder, row_threshold, col_threshold, mutated_samples, target_gene):
    # Cluster the rows and columns using hierarchical clustering
    row_linkage = linkage(expression_df_heatmap, method='ward')
    col_linkage = linkage(expression_df_heatmap.T, method='ward')

    # Assign cluster labels using cluster
    row_clusters = fcluster(row_linkage, t=row_threshold, criterion='maxclust')
    col_clusters = fcluster(col_linkage, t=col_threshold, criterion='maxclust')

    row_score, col_score = evaluate_cluster_sil_score(expression_df_heatmap, row_clusters, col_clusters)
    mutated_col_score = evaluate_cluster_sil_score_mutated_samples(expression_df_heatmap, col_clusters, mutated_samples)

    # Create dictionaries to store cluster information
    row_cluster_info = {f"Cluster {cluster}": expression_df_heatmap.index[row_clusters == cluster] for cluster in np.unique(row_clusters)}
    # col_cluster_info = {f"Cluster {cluster}": expression_df_heatmap.columns[col_clusters == cluster] for cluster in np.unique(col_clusters)}

    col_cluster_info = {}
    for cluster in np.unique(col_clusters):
        samples_in_cluster = expression_df_heatmap.columns[col_clusters == cluster]
        mutated_samples_in_cluster = [sample for sample in samples_in_cluster if sample in set(mutated_samples)]
        col_cluster_info[f"Cluster {cluster}"] = [
            f"samples: {samples_in_cluster}",
            f"mutated_samples_count': {len(mutated_samples_in_cluster)}"
        ]
        
    row_cluster_file = os.path.join(output_folder, target_gene, 'row_clusters.txt')
    with open(row_cluster_file, 'w') as file:
        for cluster, items in row_cluster_info.items():
            file.write(f"{cluster}:\n")
            file.write(f"{', '.join(items)}\n\n")

    # Save column cluster information to a text file
    col_cluster_file = os.path.join(output_folder, target_gene, 'col_clusters.txt')
    with open(col_cluster_file, 'w') as file:
        for cluster, items in col_cluster_info.items():
            file.write(f"{cluster}:\n")
            file.write(f"{', '.join(items)}\n\n")
    
    return row_linkage, col_linkage, row_cluster_info, col_cluster_info, row_score, col_score, mutated_col_score


def plot_and_save_dendrograms(row_linkage, col_linkage, expression_df_heatmap, output_folder, target_gene):
    # Plot the row dendrogram
    plt.figure(figsize=(12, 10))
    # play around with fig sizes 
    gene_plot = dendrogram(row_linkage, 
                           labels=expression_df_heatmap.index, 
                           truncate_mode='lastp', 
                           p=int(0.75 * len(row_linkage)), show_leaf_counts=True)
    plt.title('Row Dendrogram')
    
    # Save the row dendrogram plot to a file
    row_dendrogram_filename = f'{output_folder}/{target_gene}/row_dendrogram.png'
    plt.savefig(row_dendrogram_filename)
    plt.close()

    # Plot the column dendrogram
    # play around with fig sizes 
    plt.figure(figsize=(12, 10))
    sample_plot = dendrogram(col_linkage, labels=expression_df_heatmap.columns, orientation='top', 
                             truncate_mode='lastp', p=int(0.75 * len(row_linkage)), show_leaf_counts=True)
    plt.title('Column Dendrogram')
    
    # Save the column dendrogram plot to a file
    col_dendrogram_filename = f'{output_folder}/{target_gene}/col_dendrogram.png'
    plt.savefig(col_dendrogram_filename)
    plt.close()

    return gene_plot, sample_plot
