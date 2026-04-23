# get LD matrices for CAVIAR in gnu parallel parallelized manner
# instead of sequentially obtaining LD matrices, get specified number (PARALLEL_JOB_COUNT, argument 4) at a time
# get input variables
export CAVIAR_OUTPUT_DIR=$1
export CAVIAR_OUTPUT_PREFIX=$2
export BFILE_PREFIX_PATH=$3
# set below to 0 to use all available threads
export PARALLEL_JOB_COUNT=$4
# prepare LD matrices with plink --extract --r2 square yes-really
# make results folder
mkdir -p ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/RESULTS
# replace while loop with gnu parallel --csv input option
parallel --header : \
	--env BFILE_PREFIX_PATH,CAVIAR_OUTPUT_DIR,BFILE_PREFIX_PATH,PARALLEL_JOB_COUNT \
	-j ${PARALLEL_JOB_COUNT} \
	--csv "plink --bfile ${BFILE_PREFIX_PATH} --extract ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/{phenotype_id}_variant_id.txt \
		--r2 square yes-really \
		--out ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/{phenotype_id}_LD" \
	:::: ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/CAVIAR_LD_${CAVIAR_OUTPUT_PREFIX}_calc.csv
# while IFS=',' read -r pheno chr
# do
#	plink --bfile ${BFILE_PREFIX_PATH} --extract ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/${pheno}_variant_id.txt \
#      --r2 square yes-really \
#      --out ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/${pheno}_LD
# done < <(tail -n +2 ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/CAVIAR_LD_${CAVIAR_OUTPUT_PREFIX}_calc.csv)