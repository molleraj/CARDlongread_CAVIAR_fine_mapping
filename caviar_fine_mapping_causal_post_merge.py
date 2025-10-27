# prepare CAVIAR fine mapping causal post probabilities with tensorQTL results
# based on Jupyter notebooks written for HBCC and NABEC QTLs by Kensuke Daida
# first import necessary libraries
import pandas as pd
import numpy as np
# from pgenlib import PgenReader
# import pgenlib as pgen
import itertools
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
import argparse
# import rpy2.robjects as ro
# from rpy2.robjects.packages import importr

# subroutine to parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Find most likely causal variant in CAVIAR outputs and append variant ID plus causal post probability to filtered tensorQTL cis-QTL map file.")
    # import significance filtered tensorQTL cis-QTL map file
    parser.add_argument("--cis_map_file", required=True, help="Path to input *filtered* tensorQTL cis-QTL phenotype/variant association map file.")
    # set CAVIAR output directory
    parser.add_argument("--caviar_dir", required=True, help="Path to CAVIAR fine mapping analysis directory.")
    # set output prefix
    parser.add_argument("--output_prefix", required=True, help="Name prefix for output files (e.g., nabec_July_2024_rna_TPM_SV_harmonized_prun).")
    # set QTL type (SV or SV+SNV)
    parser.add_argument("--variant_type", choices=["SV", "SV+SNV"], required=True, help="tensorQTL run type (SV or SV+SNV).")
    # return parsed arguments
    return parser.parse_args()
    
    
# main script subroutine
def main():
    args=parse_args()
    # load cis_df
    cis_df = pd.read_csv(args.cis_map_file,index_col=0)
    
    if (args.variant_type == "SV"):
        # For SV only QTL
        # Initialize an empty DataFrame to store results
        results = pd.DataFrame(columns=['Phenotype', 'Chromosome', 'CAVIAR Fine Mapping Top SV', 'CAVIAR Fine Mapping Top SV Causal Post Probability'])

        # Loop through each phenotype and chromosome in the cis_df DataFrame
        for index, row in cis_df.iterrows():
            # Load the data for the current phenotype
            pheno = index
            CHR = row['chr']
            try:
                temp = pd.read_csv(f'{args.caviar_dir}/{args.output_prefix}/RESULTS/{pheno}_caviar_post', sep='\t')
                
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
                results = pd.concat([results,pd.DataFrame({
                    'Phenotype': pheno,
                    'Chromosome': CHR,
                    'CAVIAR Fine Mapping Top SV': TOPSV_SNP_ID,
                    'CAVIAR Fine Mapping Top SV Causal Post Probability': TOPSV_Causal_Post_Prob
                }, index=[0])], ignore_index=True)

            except Exception as e:
                print(f"Error processing {pheno} for {CHR}: {e}")
            
    if (args.variant_type == "SV+SNV"):
        # Initialize an empty DataFrame to store results
        results = pd.DataFrame(columns=['Phenotype', 'Chromosome', 'CAVIAR Fine Mapping Top SV', 'CAVIAR Fine Mapping Top SV Causal Post Probability', 
                                    'CAVIAR Fine Mapping Top SNV', 'CAVIAR Fine Mapping Top SNV Causal Post Probability'])

        # Loop through each phenotype and chromosome in the cis_df DataFrame
        for index, row in cis_df.iterrows():
            # Load the data for the current phenotype
            pheno = index
            CHR = row['chr']
            try:
                temp = pd.read_csv(f'{args.caviar_dir}/{args.output_prefix}/RESULTS/{pheno}_caviar_post', sep='\t')
        
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
            
                    # Preset TOPSNV to none
                    # Initialize empty TOPSNV dataframe
                    TOPSNV = pd.DataFrame(columns=['SNP_ID','Prob_in_pCausalSet','Causal_Post._Prob.'])
                    # Extract top single nucleotide variant (SNV) (not containing napu in the name)
                    TOPSNV = temp[~temp['SNP_ID'].str.contains('napu')].iloc[0, :]

                # Append the results in one row to the DataFrame
                results = pd.concat([results,pd.DataFrame({
                    'Phenotype': pheno,
                    'Chromosome': CHR,
                    'CAVIAR Fine Mapping Top SV': TOPSV_SNP_ID,
                    'CAVIAR Fine Mapping Top SV Causal Post Probability': TOPSV_Causal_Post_Prob,
                    'CAVIAR Fine Mapping Top SNV': TOPSNV['SNP_ID'],
                    'CAVIAR Fine Mapping Top SNV Causal Post Probability': TOPSNV['Causal_Post._Prob.']
                }, index=[0])], ignore_index=True)

            except Exception as e:
                print(f"Error processing {pheno} for {CHR}: {e}")
    
    # export initial CAVIAR results for each phenotype
    results.to_csv(f'{args.caviar_dir}/{args.output_prefix}/RESULTS/{args.output_prefix}_caviar_causal_post_probs_table.csv',index=False)        
            
    # join CAVIAR and tensorQTL cis-QTL results
    # rename results table columns
    results = results.rename(columns={'Phenotype':'phenotype_id'})
    # reset cis_df indices
    cis_df = cis_df.reset_index()
    # get cis_df columns of interest for final output table
    cis_df = cis_df.loc[:,['phenotype_id','variant_id','af','pval_nominal','slope','slope_se','pval_perm','bh_fdr','qval']]
    # join tables on phenotype id
    cis_caviar_join_df = cis_df.merge(results,on='phenotype_id',how='left')
    # below only if SV AND SNV QTL:
    if (args.variant_type == "SV+SNV"):
        # add additional columns to dataframe
        cis_caviar_join_df['cis-QTL Top Variant Type']=None
        cis_caviar_join_df['CAVIAR Top Variant Type']=None
        cis_caviar_join_df['cis-QTL/CAVIAR Top Variant Type Concordance']=None
        # get most significant variant type from cis map
        cis_caviar_join_df['cis-QTL Top Variant Type']=cis_caviar_join_df['variant_id'].str.contains('_').replace({True: 'SV', False: 'SNV'})
        # get most significant variant type from CAVIAR fine mapping
        cis_caviar_join_df['CAVIAR Top Variant Type']=(cis_caviar_join_df['CAVIAR Fine Mapping Top SV Causal Post Probability']>cis_caviar_join_df['CAVIAR Fine Mapping Top SNV Causal Post Probability']).replace({True: 'SV', False: 'SNV'})
        # determine whether the same variant type is most significant in both
        cis_caviar_join_df['cis-QTL/CAVIAR Top Variant Type Concordance']=(cis_caviar_join_df['CAVIAR Top Variant Type']==cis_caviar_join_df['cis-QTL Top Variant Type'])
    # Display or save results as needed
    cis_caviar_join_df.to_csv(f'{args.caviar_dir}/{args.output_prefix}/RESULTS/{args.output_prefix}_cisqtl_caviar_join_table.csv',index=False)
    
# run main subroutine
if __name__ == "__main__":
    main()