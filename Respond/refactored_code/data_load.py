import utils
from pydeseq2 import preprocessing
import pandas as pd


def load_maf_data(file_path):
    # Define columns of interest
    columns = ["Hugo_Symbol", "Tumor_Sample_Barcode"]
    df = pd.read_csv(file_path, sep='\t', comment="#", usecols=columns)
    df.rename(columns={'Hugo_Symbol': 'gene', 'Tumor_Sample_Barcode': 'sample'}, inplace=True)
    df['mutation'] = 1
    return df


def load_txt_file(file_path):
    # Read the .txt file into a pandas DataFrame
    df = pd.read_csv(file_path, sep='\t')  # Adjust the separator if needed

    # taking DESeq2’s median of ratios to normalize for RNA counts per sample
    # makes comparisons comparable across samples
    # df = preprocessing.deseq2_norm(df.T)[0]

    # take log2 of expression data to scale expression data
    # df = df.map(utils.log2)

    return df

