# prepare lists of variants and z-scores for CAVIAR fine mapping from tensorQTL results
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
    parser = argparse.ArgumentParser(description="Prepare variant lists and z-score tables for fine mapping 100 most significantly associated variants per phenotype for significantly associated phenotype-variant pairs.")
    # specify output prefix
    parser.add_argument("--output_prefix", required=True, help="Name prefix for output files (e.g., nabec_July_2024_rna_TPM_SV_harmonized_prun).")
    # argument for cis QTL tensorQTL map file
    parser.add_argument("--cis_map_file", required=True, help="Path to input tensorQTL cis-QTL phenotype/variant association map file.")
    # argument for parquet files
    parser.add_argument("--cis_parquet_dir", required=True, help="Path to parquet files containing association information for all cis-QTL phenotype/variant pairs.")
    # argument for CAVIAR analysis directory
    parser.add_argument("--caviar_dir", required=True, help="Path to CAVIAR fine mapping analysis directory.")
    # return parsed arguments
    return parser.parse_args()

# main script subroutine
def main():
    # load arguments
    args=parse_args()
    # load cis-QTL map
    cis_df = pd.read_csv(args.cis_map_file,index_col=0)
    # filter for significant variants based on FDR adjusted q-value (>0.05)
    cis_df = cis_df[cis_df['qval'] < 0.05]
    # calculate z score based on slope and slope se
    cis_df['z_score'] = cis_df['slope']/cis_df['slope_se']
    # remove napu prefix from variant ids where necessary
    cis_df['variant_id'] = cis_df['variant_id'].str.replace('napu_','')
    # output file for making CAVIAR input LD matrices
    cis_df[['chr']].to_csv(f'{args.caviar_dir}/{args.output_prefix}/CAVIAR_LD_{args.output_prefix}_calc.csv')
    # output filtered, renamed variant cis-map
    cis_df.to_csv(f'{args.caviar_dir}/{args.output_prefix}/{args.output_prefix}.qval_filtered.cis.map.csv')

    # make copy of cis-QTL map file
    temp = cis_df.copy()

    temp = temp.sort_values('chr')

    # Load the parquet files for all chromosomes upfront
    chromosome_data = {}

    for chr_num in temp['chr'].unique():
        chromosome_data[chr_num] = pd.read_parquet(f'{args.cis_parquet_dir}/{args.output_prefix}.cis_qtl_pairs.{chr_num}.parquet')
    
    # iterate through filtered cis map to get top 100 variants for each phenotype
    for index, row in temp.iterrows():
        pheno = index
        # define CHR (chromosome) variable
        CHR = row['chr']
        # Use the pre-loaded data for the current chromosome
        TENSOR = chromosome_data[CHR]
        # get chromosome data for which phenotype id matches that in filtered cis map
        NUM = TENSOR.loc[(TENSOR.phenotype_id == pheno)]


        # Use the pre-loaded data for the current chromosome
        TENSOR = chromosome_data[CHR]
        NUM = TENSOR.loc[(TENSOR.phenotype_id ==pheno)]
        NUM['zscore'] = NUM['slope']/NUM['slope_se']
        NUM = NUM.sort_values('pval_nominal')
        # Sort by pval_nominal to get the top 100
        NUM_sorted = NUM.sort_values('pval_nominal')

        # Get the top 100 variants
        top_100_variants = NUM_sorted.head(100)

        # Check if any variant_id starts with 'napu_'
        napu_variants = top_100_variants[top_100_variants['variant_id'].str.startswith('napu_')]

        # If no 'napu_' variant exists in the top 100
        if napu_variants.empty:
            # Find the variant that starts with 'napu_' and has the lowest pval_nominal (i.e., get at least ONE SV in fine mapping set)
            napu_lowest_variant = NUM_sorted[NUM_sorted['variant_id'].str.startswith('napu_')].nsmallest(1, 'pval_nominal')
    
            if not napu_lowest_variant.empty:
                # Replace the 100th variant with the 'napu_' variant
                top_100_variants = top_100_variants[:-1]  # Remove the last variant (100th)
                top_100_variants = top_100_variants.append(napu_lowest_variant)  # Append the napu_ variant

        # Output the top 100 variants
        #top_100_variants

        top_100_variants[['variant_id','zscore']].to_csv(f'{args.caviar_dir}/{args.output_prefix}/{pheno}_zscore.txt',index=False,header=None,sep='\t')
        top_100_variants['variant_id'].to_csv(f'{args.caviar_dir}/{args.output_prefix}/{pheno}_variant_id.txt',header=None,index=False)
        
# run main subroutine
if __name__ == "__main__":
    main()
    