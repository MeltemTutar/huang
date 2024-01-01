
import os
import numpy as np
import pandas as pd
from statsmodels.stats import multitest
from scipy.stats import ranksums
from statistics import mean


def log2(value):
    return np.log2(value + 1)

def calculate_log_fold(list_mutated, list_non_mutated, smoothing_factor=1e-8): 
    mean_mutated = mean(list_mutated)  + smoothing_factor
    mean_non_mutated = mean(list_non_mutated) + smoothing_factor
        
    return np.log2(mean_mutated/mean_non_mutated)


def calculate_z_score(list_mutated, list_non_mutated):
    mean_mutated = mean(list_mutated) 
    mean_non_mutated = mean(list_non_mutated) 
    all_expression_values = list_non_mutated + list_mutated
    z_score = (mean_mutated - mean_non_mutated) / np.std(all_expression_values)
    return z_score


def calculate_pvalue_and_effect_size_wilcox_ranksum(list_mutated, list_non_mutated):
     test = ranksums(list_mutated, list_non_mutated)
     u_statistic = test.statistic 
     effect_size = u_statistic / (len(list_mutated) * len(list_non_mutated))
     pvalue=test.pvalue
     return pvalue, effect_size


def calculate_adjusted_pvalue(pvalues, method='fdr_bh'):
    _, corrected_pvalues, _, _ = multitest.multipletests(pvalues, method=method)
    return corrected_pvalues


def generate_stats_per_gene(express_mut_genes_df, target_gene, output_folder):
    """
    Given expression and mutation data calculates LogFC, pvalue, mean of expression per gene 
    for a given target_gene 
    
    """
    gene_df = express_mut_genes_df[express_mut_genes_df['gene']==target_gene]

    if len(gene_df) == 0:
        raise ValueError("This gene is not valid! No mutations with this gene exist") 

    mutated_samples = gene_df[gene_df['mutation'] == 1]['sample']
    non_mutated_samples = gene_df[gene_df['mutation'] == 0]['sample']
    
    mutated_individuals_expression = express_mut_genes_df[express_mut_genes_df['sample'].isin(mutated_samples)]
    
    non_mutated_individuals_expression = express_mut_genes_df[express_mut_genes_df['sample'].isin(non_mutated_samples)]

    # gather data of mutated and non-mutated genes into lists
    mutated_individuals_data = mutated_individuals_expression.groupby(['gene'])['gene_expression'].apply(
    lambda x: list(x)).to_frame().reset_index().rename(columns={'gene_expression': 'gene_expression_mutated'})
    non_mutated_individuals_data = non_mutated_individuals_expression.groupby(['gene'])['gene_expression'].apply(
    lambda x: list(x)).to_frame().reset_index().rename(columns={'gene_expression': 'gene_expression_non_mutated'})

    # combine mutated and unmutated data into one df
    combined_data= pd.merge(mutated_individuals_data, non_mutated_individuals_data, on='gene', how='inner')

    # calculate fold change and p-value and z-score
    combined_data['logFC']= combined_data.apply(
        lambda x: calculate_log_fold(x.gene_expression_mutated, x.gene_expression_non_mutated), axis=1)
    # calculate wilcox pvalue
    combined_data[['pvalue', 'effect_size']] = combined_data.apply(
    lambda x: pd.Series(calculate_pvalue_and_effect_size_wilcox_ranksum(x.gene_expression_mutated, x.gene_expression_non_mutated)),
    axis=1
    )
    combined_data['expression_mutated_mean'] = combined_data['gene_expression_mutated'].apply(mean)
    combined_data['expression_nonmutated_mean'] = combined_data['gene_expression_non_mutated'].apply(mean)


    # calculate adjusted p-value
    combined_data['adjusted_pvalue'] = calculate_adjusted_pvalue(combined_data['pvalue'].values)

    # Output data to csv 
    combined_data.drop(columns=['gene_expression_mutated', 'gene_expression_non_mutated'], inplace=True)
    
    os.makedirs(f'{output_folder}/{target_gene}', exist_ok=True)
    output_filename = f'{output_folder}/{target_gene}/{mutated_samples.count()}_{non_mutated_samples.count()}_logfc_pvalue.csv'
    
    print(f"outputting data to {output_filename}")
    combined_data.to_csv(output_filename, index=False)

    return combined_data, mutated_samples, output_filename


def get_mutated_status(expression_df_heatmap, individuals_mutated_target_gene, output_folder, target_gene):
    individuals_mutated_target_gene = list(individuals_mutated_target_gene)
    
    mutated_status = [1 if item in individuals_mutated_target_gene else 0 for item in expression_df_heatmap.columns]

    sample_categories_df = pd.DataFrame({
        'Sample': expression_df_heatmap.columns,
        'Mutation Status': mutated_status
    })

    output_filename = f'{output_folder}/{target_gene}/sample_mutation_status.csv'
    sample_categories_df.to_csv(output_filename,  index=False)
    return sample_categories_df

def output_genes_with_highest_score(gene_scores, n, sortby='mutated_sample_sil_score'):
    # Sort genes based on mutated_sample_sil_score in descending order
    sorted_genes = sorted(gene_scores.items(), key=lambda x: x[1][sortby], reverse=True)

    # return the top genes
    return [(gene, scores) for gene, scores in sorted_genes[:n]]
