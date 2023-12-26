

import plot, utils, data_load, clustering, constants 
import warnings
import random


# Suppress FutureWarnings from multiple modules
modules_to_ignore = ("seaborn", "sklearn")

for module in modules_to_ignore:
    warnings.filterwarnings("ignore", category=FutureWarning, module=module)

maf_df = data_load.load_maf_data(file_path=constants.maf_file_path)
expression_df = data_load.load_txt_file_into_dataframe(file_path=constants.expression_file_path)
expression_df_melted = data_load.reformat_expression_data(df=expression_df)
mutation_expression_df_melted = data_load.preprocess_and_combine_mutation_expression(maf_df= maf_df, expression_df = expression_df_melted)

sil_scores_genes = {}
for target_gene in constants.genes:
    # calculate logfc, pvalue
    volcano_plot_df, individuals_mutated_target_gene, volcano_input_filename = utils.generate_stats_per_gene(
        express_mut_genes_df=mutation_expression_df_melted, 
        target_gene=target_gene,
        output_folder=target_gene)

    # get top 100 genes with the lowest p-values
    heatmap_data = data_load.generate_expression_heatmap(expression_df=expression_df, volcano_plot_df=volcano_plot_df, n=100)

    mutated_status_df = utils.get_mutated_status(expression_df_heatmap=heatmap_data, 
                                             individuals_mutated_target_gene=individuals_mutated_target_gene, 
                                             output_folder=target_gene)

    # cluster top 100 genes and samples
    row_linkage, col_linkage, row_clusters, col_clusters, gene_sil_score, sample_sil_score = clustering.hierarchical_clustering(
                                                                                                expression_df_heatmap=heatmap_data, 
                                                                                                output_folder=target_gene,
                                                                                                row_threshold=7, 
                                                                                                col_threshold=2)

    sil_scores_genes[target_gene]=gene_sil_score


print(sil_scores_genes)

# # getting indviduals with mutation for heatmap
# # TODO: combine this func into other functions
# mutated_status_df = utils.get_mutated_status(expression_df_heatmap=heatmap_data, 
#                                             individuals_mutated_target_gene=individuals_mutated_target_gene, 
#                                             output_folder=target_gene)

# # plot clustered heatmap
# clustered_heatmap = plot.create_clustered_heatmap_and_save(heatmap_data, 
#                                                       row_linkage, 
#                                                       col_linkage,
#                                                       mutated_status_df,
#                                                       output_folder=target_gene)

# # plot truncated dendograms 
# row_dendrogram, col_dendrogram = clustering.plot_and_save_dendrograms(row_linkage, 
#                                                            col_linkage, 
#                                                            heatmap_data, 
#                                                            output_folder=target_gene)

# # plot histograms of adjusted and regular pvalues 
# plot.histogram_of_column_and_save(volcano_plot_df, 'pvalue', output_folder=target_gene)
# plot.histogram_of_column_and_save(volcano_plot_df, 'adjusted_pvalue', output_folder=target_gene)

# # calculate significant genes determined by cutoffs
# significant_genes_positive, significant_genes_negative  = plot.volcano_plot(input_file_path=volcano_input_filename, 
#                                                                         yaxis='pvalue', 
#                                                                         xaxis='logFC', 
#                                                                         output_folder=target_gene,
#                                                                         significance_threshold=0.05, 
#                                                                         logfold_positive_threshold=2, 
#                                                                         logfold_negative_threshold=-2)

# # create box plots of significant genes (uoregulated and downregulated genes)
# plot.create_gene_expression_boxplot(expression_df, significant_genes_positive, mutated_status_df, positive=1, output_folder=target_gene)
# plot.create_gene_expression_boxplot(expression_df, significant_genes_negative, mutated_status_df, positive=0, output_folder=target_gene)

# # Reset warnings to default behavior after the code block
# warnings.resetwarnings()