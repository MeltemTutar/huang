import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def volcano_plot(input_file_path, 
                 yaxis, 
                 xaxis,
                 output_folder,
                 target_gene,
                 significance_threshold, 
                 logfold_positive_threshold, 
                 logfold_negative_threshold):
    df = pd.read_csv(input_file_path)

    # Apply -log10 transformation to the p-value
    df['-log10_pvalue'] = -np.log10(df[yaxis])

    # Highlight significant points with large log-fold changes
    significant_genes_positive = df[(df[yaxis] < significance_threshold) & (df['logFC'] > logfold_positive_threshold)]
    significant_genes_negative = df[(df[yaxis] < significance_threshold) & (df['logFC'] < logfold_negative_threshold)]

    # Clip outliers in 'logFC' column for the plot
    df['logFC'] = df['logFC'].clip(lower=-10, upper=10)

    plt.scatter(x=df[xaxis], y=df['-log10_pvalue'], s=1, label='All Genes', alpha=0.5)
    plt.scatter(x=significant_genes_positive[xaxis], y=significant_genes_positive['-log10_pvalue'], s=10, c='red', marker='^', label='Significant Genes (Positive LogFC)')
    plt.scatter(x=significant_genes_negative[xaxis], y=significant_genes_negative['-log10_pvalue'], s=10, c='blue', marker='v', label='Significant Genes (Negative LogFC)')

    plt.xlabel(xaxis)
    plt.ylabel(f'-log10({yaxis})')
    plt.title('Volcano Plot')
    plt.axhline(-np.log10(significance_threshold), color='gray', linestyle='--', label=f'Significance Threshold ({significance_threshold})')
    plt.axvline(logfold_positive_threshold, color='gray', linestyle='--', label=f'Positive Log-Fold Change Threshold ({logfold_positive_threshold})')
    plt.axvline(logfold_negative_threshold, color='gray', linestyle='--', label=f'Negative Log-Fold Change Threshold ({logfold_negative_threshold})')
    
     # Create a separate legend outside the plot
    fig = plt.gcf()
    handles, labels = plt.gca().get_legend_handles_labels()
    fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    
    
    output_filename = f'{output_folder}/{target_gene}/volcano_plot.png'
    plt.savefig(output_filename, bbox_inches='tight')
    plt.close()
    significant_genes_positive.to_csv(f'{output_folder}/{target_gene}/signif_genes_positive.csv')
    significant_genes_negative.to_csv(f'{output_folder}/{target_gene}/signif_genes_negative.csv')

    return significant_genes_positive, significant_genes_negative


def create_gene_expression_boxplot(expression_df, significant_genes_df, mutated_status_df, output_folder, target_gene, positive=1 ):

    if len(significant_genes_df)==0:
        if positive ==1:
            print(f"for gene {target_gene} no significant positive  genes found, so no box plot created")
        elif positive ==0:
            print(f"for gene {target_gene} no significant negative  genes found, so no box plot created")
        return
    # Filter for significant genes
    filtered_expression_df_genes_of_interest = expression_df.loc[significant_genes_df.gene]

    # Unstack the expression DataFrame
    unstack_expression = filtered_expression_df_genes_of_interest.unstack().reset_index()
    unstack_expression.columns = ['Sample', 'Gene', 'Expression']

    # Merge the expression matrix with the mutation information
    boxplot_df = pd.merge(unstack_expression, mutated_status_df, how='inner', on='Sample')

    # Rename columns for clarity
    boxplot_df.columns = ['Sample', 'Gene', 'Expression', 'Target Gene Mutation Status']

    # Set up the plot
    plt.figure(figsize=(12, 8))

    # Create a box plot with 'hue' for each gene and 'dodge' for Mutation Status 
    sns.boxplot(x='Gene', y='Expression', data=boxplot_df, hue='Target Gene Mutation Status', dodge=True)

    # Customize the plot
    plt.xlabel('Gene')
    plt.ylabel('Gene Expression')
    plt.title('Box Plots of Gene Expression for Mutated and Non-Mutated Individuals')

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')

    if positive == 1:
        output_filename = f'{output_folder}/{target_gene}/positive_genes_expression_boxplot.png'
    else:
        output_filename = f'{output_folder}/{target_gene}/negative_genes_expression_boxplot.png'
    plt.savefig(output_filename)
    plt.close()


def histogram_of_column_and_save(df, column, output_folder, target_gene):
    # Plot the distribution of p-values
    plt.figure(figsize=(10, 6))
    plt.hist(df[column], bins=30, color='blue', edgecolor='black')
    plt.title(f'Distribution of {column}')
    plt.xlabel('P-values')
    plt.ylabel('Frequency')

    output_filename = f'{output_folder}/{target_gene}/histogram_{column}.png'
    # Save the histogram plot to a file
    plt.savefig(output_filename)
    plt.close()

def create_heatmap_and_save(expression_df_heatmap,
                            output_folder,
                            output_file_name,
                            target_gene
                            ):

    sns.set(font_scale=1.0) 

    # Create a clustered heatmap
    ax = sns.heatmap(
        expression_df_heatmap,
        annot=False,
        xticklabels=False, 
        yticklabels=expression_df_heatmap.index,
        cmap='coolwarm'
    )

    #Set the size of the overall figure
    # sns_df.fig.set_size_inches(15, len(expression_df_heatmap) * 0.2)  

    output_filename = f'{output_folder}/{target_gene}/{output_file_name}'

    # sns_df.ax_heatmap.set_yticklabels(sns_df.ax_heatmap.get_yticklabels(), rotation=0)
    
    # Save the plot to a file
    plt.savefig(output_filename, bbox_inches='tight')

    # Close the plot to prevent displaying it in the notebook (optional)
    plt.close()

    return ax
    
    


def create_clustered_heatmap_and_save(expression_df_heatmap, 
                                      row_linkage, 
                                      col_linkage, 
                                      sample_mutation_df, 
                                      output_folder,
                                      output_file_name,
                                      target_gene):
    # Group labels and colors for color bar
    type_map = {1: 'red', 0: 'yellow'}

    sns.set(font_scale=1.0) 

    # Create a clustered heatmap
    clustered_df = sns.clustermap(
        expression_df_heatmap,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        cmap='viridis',
        annot=False,
        fmt=".1f",
        linewidths=.5,
        col_colors=sample_mutation_df.set_index('Sample')['Mutation Status'].map(type_map),
        z_score=0,
        cbar_kws={"shrink": 0.7, "aspect": 30}  # Adjust color bar size
    )

    #Set the size of the overall figure
    clustered_df.fig.set_size_inches(15, len(expression_df_heatmap) * 0.2)  

    output_filename = f'{output_folder}/{target_gene}/{output_file_name}'

    clustered_df.ax_heatmap.set_yticklabels(clustered_df.ax_heatmap.get_yticklabels(), rotation=0)
    
    # Save the plot to a file
    clustered_df.savefig(output_filename, bbox_inches='tight')

    # Close the plot to prevent displaying it in the notebook (optional)
    plt.close()

    return clustered_df
