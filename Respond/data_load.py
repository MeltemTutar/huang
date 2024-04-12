import utils
import pandas as pd
from pydeseq2 import preprocessing
from scipy.stats import zscore

def load_maf_data(file_path, columns = ["Hugo_Symbol", "Tumor_Sample_Barcode"]):
    df = pd.read_csv(file_path, sep='\t', comment="#", usecols=columns)
    df.rename(columns={'Hugo_Symbol': 'gene', 'Tumor_Sample_Barcode': 'sample'}, inplace=True)
    df['mutation'] = 1
    return df

def load_txt_file_into_dataframe(file_path):
    # Read the .txt file into a pandas DataFrame
    df = pd.read_csv(file_path, sep='\t')  # Adjust the separator if needed
    
    # normalize to account for RNA sequencing depth (samples having different totals of RNA expression)
    # this allows us to make comparisons for the same gene across samples
    df = preprocessing.deseq2_norm(df.T)[0].T
    
    # take log2 of expression data to scale expression data. Reduces the effect of outliers
    df = df.map(utils.log2)

    return df


def reformat_expression_data(df):
    # Combine column names and index names into rows for every element
    melted_df = pd.melt(df.reset_index(), id_vars=['index'], var_name='column', col_level=0)

    # Rename columns
    melted_df.rename(columns={'index': 'gene', 'column': 'sample', 'value': 'gene_expression'}, inplace=True)

    return melted_df


def preprocess_and_combine_mutation_expression(maf_df, expression_df):
    ''' filter the expression data to those that have whole genome sequencing 
     i.e. appear in the mutation data frame (maf)
     join mutation and expression data
     '''
    
    exon_seq_samples = maf_df['sample'].unique()
    filtered_expression_df = expression_df[expression_df['sample'].isin(exon_seq_samples)]
    
    all_rows = expression_df.shape[0]
    filt_rows = filtered_expression_df.shape[0]
    
    percentage_filtered = (all_rows - filt_rows) / all_rows
    
    print('fraction of rows filtered is', percentage_filtered)

    express_mut_genes_df = pd.merge(maf_df, filtered_expression_df, on=['gene', 'sample'], how='right')

    express_mut_genes_df['mutation'].fillna(0, inplace=True)

    return express_mut_genes_df

def generate_expression_heatmap(expression_df, volcano_plot_df, n, top):
    # Get the indices of the top n rows based on having smallest pvalues

    if top == True: 
        indices = volcano_plot_df['pvalue'].nsmallest(n).index
    else:
        indices = (-volcano_plot_df['pvalue']).nsmallest(n).index
    
    top_n_rows = volcano_plot_df.loc[indices].set_index('gene')

    # go back to the expression df and make a heatmap of the top n genes
    expression_df_heatmap = expression_df.loc[top_n_rows.index]
    
    return expression_df_heatmap


def generate_expression_heatmap_positive_negative_same(expression_df, volcano_plot_df, n, significance_threshold):
    # small p value, positive fold change
    positive_sig = volcano_plot_df[(volcano_plot_df['pvalue'] < significance_threshold)].sort_values('logFC', ascending=False).head(n)
    positive_sig_expression = expression_df.loc[positive_sig['gene']]

    # small p value, negative fold change
    negative_sig = volcano_plot_df[(volcano_plot_df['pvalue'] < significance_threshold)].sort_values('logFC', ascending=True).head(n)
    negative_sig_expression = expression_df.loc[negative_sig['gene']]

    # non significant. highest p valyes 
    non_sig= volcano_plot_df.sort_values(by=['pvalue'], ascending=[False]).head(n)
    non_sig_expression = expression_df.loc[non_sig['gene']]

    positive_negative_same_expression_df = pd.concat([positive_sig_expression, negative_sig_expression, non_sig_expression])

    # Calculate z-score by row (gene)
    positive_negative_same_expression_df = positive_negative_same_expression_df.apply(zscore, axis=1)

    return positive_negative_same_expression_df,  positive_sig, negative_sig, non_sig

