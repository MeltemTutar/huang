
import data_load, constants 
import pandas as pd
import numpy as np
import os

def create_gsea_expression_input(output_folder):
    expression_df = data_load.load_txt_file_into_dataframe(file_path=constants.expression_file_path)
    mutation_df = data_load.load_maf_data(file_path=constants.maf_file_path)

    # get expression data of individuals with at least one mutation (i.e. has sequencing data)
    expression_df = expression_df.T[expression_df.T.index.isin(mutation_df['sample'].unique())].T
    # Add a new column with 'NA' values at a specified position
    position = 0
    new_column_name = 'description'
    expression_df.insert(position, new_column_name, 'NA')

    # Reset the index and rename the index column
    expression_df = expression_df.reset_index().rename(columns={'index': 'NAME'})
    os.makedirs(output_folder, exist_ok=True)

    output_path = f'{output_folder}/GSEA_expression_input.gct'
    with open(output_path, 'w') as file:
        file.write("#1.2\n")
        file.write(f"{expression_df.shape[0]}\t{expression_df.shape[1]- 2}\n")

    # Save the DataFrame to a GCT file
    expression_df.to_csv(output_path, index=False, sep='\t',  mode='a')
    return output_path

def create_gsea_expression_input_preranked(output_folder, target_gene, stats_path, logFC=1, pvalue=0.05):
    stats_df=pd.read_csv(stats_path)

    # filter cols
    stats_df = stats_df[(abs(stats_df['logFC']) >= logFC) & (stats_df['pvalue'] <+ pvalue)]

    # Extract required columns
    stats_df = stats_df[['gene', 'logFC', 'pvalue']]
    
    # Calculate the transformed expression
    stats_df['metric'] = np.sign(stats_df['logFC']) * (-np.log10(stats_df['pvalue']))
    
    output_filename=f'{output_folder}/{target_gene}/preranked_genes.rnk'
    stats_df[['gene', 'metric']].to_csv(output_filename, sep='\t', index=False, header=False)



def create_mutation_label_gsea(gsea_expression_path, output_folder, target_gene):
    mutation_df = data_load.load_maf_data(file_path=constants.maf_file_path)

    expression_df =  pd.read_csv(gsea_expression_path, sep='\t', skiprows=range(1, 2), header=1).drop(columns=['NAME', 'description'])
    gsea_mutation_df= expression_df.T.reset_index(names='sample')[['sample']].merge(
        mutation_df[mutation_df['gene']==target_gene].drop_duplicates(), on='sample', how='left')[['sample', 'mutation']].fillna(0)

    os.makedirs(f'{output_folder}/{target_gene}', exist_ok=True)
    output_path=f'{output_folder}/{target_gene}/mutation_GSEA_input.cls'
    
    # check if first sample is 0.0 or 1.0 because this determines the order it needs to appear in cls file
    gsea_mutation_df
    if gsea_mutation_df['mutation'].iloc[0]==1.0:
        second_line = f"#\t1.0\t0.0\n"
    elif gsea_mutation_df['mutation'].iloc[0]==0.0:
        second_line = f"#\t0.0\t1.0\n"

    with open(output_path, 'w') as file:
        file.write(f"{gsea_mutation_df.shape[0]}\t{len(gsea_mutation_df['mutation'].unique())}\t1\n")
        file.write(second_line)
    
    gsea_mutation_df[['mutation']].T.to_csv(output_path,
                                        index=False, header=False, sep='\t', mode='a')
    
