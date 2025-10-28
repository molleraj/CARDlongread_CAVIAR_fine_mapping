# NIA CARD ANG CAVIAR QTL Fine Mapping Pipeline
This repository includes a pipeline for streamlined fine mapping of potentially causal QTL variants using [CAVIAR](http://genetics.cs.ucla.edu/caviar/index.html).

The pipeline is run in four steps:
1. Preparation of up to 100 most significantly associated variants per phenotype, making sure to include at least one structural variant when necessary, for table of phenotype/variant pairs passing FDR-corrected p-value (q-value) filter (typically 0.05).
2. Calculating linkage disequilibrium (LD) matrices for each variant set selected above.
3. Running CAVIAR on each variant set per phenotype.
4. Appending CAVIAR causal post probabilities to initial significant phenotype/variant pair table.

## Usage

```
usage: caviar_fine_mapping_variant_prep.py [-h] --output_prefix OUTPUT_PREFIX --cis_map_file CIS_MAP_FILE --cis_parquet_dir CIS_PARQUET_DIR --caviar_dir CAVIAR_DIR

Prepare variant lists and z-score tables for fine mapping 100 most significantly associated variants per phenotype for significantly associated phenotype-variant pairs.

optional arguments:
  -h, --help            show this help message and exit
  --output_prefix OUTPUT_PREFIX
                        Name prefix for output files (e.g., nabec_July_2024_rna_TPM_SV_harmonized_prun).
  --cis_map_file CIS_MAP_FILE
                        Path to input tensorQTL cis-QTL phenotype/variant association map file.
  --cis_parquet_dir CIS_PARQUET_DIR
                        Path to parquet files containing association information for all cis-QTL phenotype/variant pairs.
  --caviar_dir CAVIAR_DIR
                        Path to CAVIAR fine mapping analysis directory.
```
```
make_LD_matrices_for_caviar.sh CAVIAR_OUTPUT_DIR CAVIAR_OUTPUT_PREFIX BFILE_PREFIX_PATH
```
```
run_CAVIAR.sh CAVIAR_OUTPUT_DIR CAVIAR_OUTPUT_PREFIX
```
```
usage: caviar_fine_mapping_causal_post_merge.py [-h] --cis_map_file CIS_MAP_FILE --caviar_dir CAVIAR_DIR --output_prefix OUTPUT_PREFIX --variant_type {SV,SV+SNV}

Find most likely causal variant in CAVIAR outputs and append variant ID plus causal post probability to filtered tensorQTL cis-QTL map file.

optional arguments:
  -h, --help            show this help message and exit
  --cis_map_file CIS_MAP_FILE
                        Path to input *filtered* tensorQTL cis-QTL phenotype/variant association map file.
  --caviar_dir CAVIAR_DIR
                        Path to CAVIAR fine mapping analysis directory.
  --output_prefix OUTPUT_PREFIX
                        Name prefix for output files (e.g., nabec_July_2024_rna_TPM_SV_harmonized_prun).
  --variant_type {SV,SV+SNV}
                        tensorQTL run type (SV or SV+SNV).
```

## Example output
