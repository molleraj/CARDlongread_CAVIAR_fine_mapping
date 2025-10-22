# prepare CAVIAR fine mapping causal post probabilities with tensorQTL results
# based on Jupyter notebooks written for HBCC and NABEC QTLs by Kensuke Daida
# first import necessary libraries
import pandas as pd
import numpy as np
from pgenlib import PgenReader
import pgenlib as pgen
import itertools
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

# subroutine to parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Find most likely causal variant in CAVIAR outputs and append variant ID plus causal post probability to filtered tensorQTL cis-QTL map file.")
    
    
# main script subroutine
def main():
    # For SV only QTL
    # Initialize an empty DataFrame to store results
    results = pd.DataFrame(columns=['Phenotype', 'Chromosome', 'TOPSV_SNP_ID', 'TOPSV_Causal_Post._Prob.'])

    # Loop through each phenotype and chromosome in the cis_df DataFrame
    for index, row in cis_df.iterrows():
        # Load the data for the current phenotype
        pheno = index
        CHR = row['chr']
        try:
            temp = pd.read_csv(f'{CAVIAR_OUTPUT}/{set_name}/RESULTS/{pheno}_caviar_post', sep='\t')
        
            # Sort by 'Causal_Post._Prob.' in descending order
            temp = temp.sort_values('Causal_Post._Prob.', ascending=False)

            # Initialize placeholders for SV and SNV data
            TOPSV_SNP_ID = None
            TOPSV_Causal_Post_Prob = None

            # Extract top structural variant (SV) if available
            top_sv_candidates = temp[temp['SNP_ID'].str.contains('napu')]
            if not top_sv_candidates.empty:
                TOPSV = top_sv_candidates.iloc[0, :]
                TOPSV_SNP_ID = TOPSV['SNP_ID']
                TOPSV_Causal_Post_Prob = TOPSV['Causal_Post._Prob.']
            # Append the results in one row to the DataFrame
            results = results.append({
                'Phenotype': pheno,
                'Chromosome': CHR,
                'TOPSV_SNP_ID': TOPSV_SNP_ID,
                'TOPSV_Causal_Post._Prob.': TOPSV_Causal_Post_Prob
            }, ignore_index=True)

        except Exception as e:
            print(f"Error processing {pheno} for {CHR}: {e}")
            
    # For SNV only QTL
    import pandas as pd

    # Initialize an empty DataFrame to store results
    results = pd.DataFrame(columns=['Phenotype', 'Chromosome', 'TOPSV_SNP_ID', 'TOPSV_Causal_Post._Prob.', 
                                    'TOPSNV_SNP_ID', 'TOPSNV_Causal_Post._Prob.'])

    # Loop through each phenotype and chromosome in the cis_df DataFrame
    for index, row in cis_df.iterrows():
        # Load the data for the current phenotype
        pheno = index
        CHR = row['chr']
        try:
            temp = pd.read_csv(f'{CAVIAR_OUTPUT}/{set_name}/RESULTS/{pheno}_caviar_post', sep='\t')
        
            # Sort by 'Causal_Post._Prob.' in descending order
            temp = temp.sort_values('Causal_Post._Prob.', ascending=False)

            # Initialize placeholders for SV and SNV data
            TOPSV_SNP_ID = None
            TOPSV_Causal_Post_Prob = None

            # Extract top structural variant (SV) if available
            top_sv_candidates = temp[temp['SNP_ID'].str.contains('napu')]
            if not top_sv_candidates.empty:
                TOPSV = top_sv_candidates.iloc[0, :]
                TOPSV_SNP_ID = TOPSV['SNP_ID']
                TOPSV_Causal_Post_Prob = TOPSV['Causal_Post._Prob.']

            # Extract top single nucleotide variant (SNV)
            TOPSNV = temp[~temp['SNP_ID'].str.contains('napu')].iloc[0, :]

            # Append the results in one row to the DataFrame
            results = results.append({
                'Phenotype': pheno,
                'Chromosome': CHR,
                'TOPSV_SNP_ID': TOPSV_SNP_ID,
                'TOPSV_Causal_Post._Prob.': TOPSV_Causal_Post_Prob,
                'TOPSNV_SNP_ID': TOPSNV['SNP_ID'],
                'TOPSNV_Causal_Post._Prob.': TOPSNV['Causal_Post._Prob.']
            }, ignore_index=True)

        except Exception as e:
            print(f"Error processing {pheno} for {CHR}: {e}")

    # Display or save results as needed
        
    
# run main subroutine
if __name__ == "__main__":
    main()