
import os
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
import numpy as np


def hierarchical_clustering(expression_df_heatmap, output_folder, row_threshold, col_threshold):
    # Cluster the rows and columns using hierarchical clustering
    row_linkage = linkage(expression_df_heatmap, method='ward')
    col_linkage = linkage(expression_df_heatmap.T, method='ward')

    # Assign cluster labels using cluster
    row_clusters = fcluster(row_linkage, t=row_threshold, criterion='maxclust')
    col_clusters = fcluster(col_linkage, t=col_threshold, criterion='maxclust')

    # Calculate silhouette scores
    row_silhouette_score = silhouette_score(expression_df_heatmap, row_clusters)
    col_silhouette_score = silhouette_score(expression_df_heatmap.T, col_clusters)

    print(f"Row Silhouette Score: {row_silhouette_score}")
    print(f"Column Silhouette Score: {col_silhouette_score}")

    # Create dictionaries to store cluster information
    row_cluster_info = {f"Cluster {cluster}": expression_df_heatmap.index[row_clusters == cluster] for cluster in np.unique(row_clusters)}
    col_cluster_info = {f"Cluster {cluster}": expression_df_heatmap.columns[col_clusters == cluster] for cluster in np.unique(col_clusters)}

    row_cluster_file = os.path.join(output_folder, 'row_clusters.txt')
    with open(row_cluster_file, 'w') as file:
        file.write(f"Row Silhouette Score: {row_silhouette_score}\n\n")
        for cluster, items in row_cluster_info.items():
            file.write(f"{cluster}:\n")
            file.write(f"{', '.join(items)}\n\n")

    # Save column cluster information to a text file
    col_cluster_file = os.path.join(output_folder, 'col_clusters.txt')
    with open(col_cluster_file, 'w') as file:
        file.write(f"Column Silhouette Score: {col_silhouette_score}\n\n")
        for cluster, items in col_cluster_info.items():
            file.write(f"{cluster}:\n")
            file.write(f"{', '.join(items)}\n\n")
    
    return row_linkage, col_linkage, row_cluster_info, col_cluster_info, row_silhouette_score, col_silhouette_score


def plot_and_save_dendrograms(row_linkage, col_linkage, expression_df_heatmap, output_folder):
    # Plot the row dendrogram
    plt.figure(figsize=(12, 10))
    # play around with fig sizes 
    gene_plot = dendrogram(row_linkage, 
                           labels=expression_df_heatmap.index, 
                           truncate_mode='lastp', 
                           p=int(0.75 * len(row_linkage)), show_leaf_counts=True)
    plt.title('Row Dendrogram')
    
    # Save the row dendrogram plot to a file
    row_dendrogram_filename = f'{output_folder}/row_dendrogram.png'
    plt.savefig(row_dendrogram_filename)
    plt.close()

    # Plot the column dendrogram
    # play around with fig sizes 
    plt.figure(figsize=(12, 10))
    sample_plot = dendrogram(col_linkage, labels=expression_df_heatmap.columns, orientation='top', 
                             truncate_mode='lastp', p=int(0.75 * len(row_linkage)), show_leaf_counts=True)
    plt.title('Column Dendrogram')
    
    # Save the column dendrogram plot to a file
    col_dendrogram_filename = f'{output_folder}/col_dendrogram.png'
    plt.savefig(col_dendrogram_filename)
    plt.close()

    return gene_plot, sample_plot
