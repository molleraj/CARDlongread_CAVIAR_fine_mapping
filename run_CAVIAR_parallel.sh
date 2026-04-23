# run CAVIAR on z score tables and LD matrices in gnu parallel parallelized manner
# instead of sequentially obtaining LD matrices, get specified number (PARALLEL_JOB_COUNT, argument 4) at a time
# get input variables
export CAVIAR_OUTPUT_DIR=$1
export CAVIAR_OUTPUT_PREFIX=$2
# set below to 0 to use all available threads
export PARALLEL_JOB_COUNT=$3
# run CAVIAR
# replace while loop with gnu parallel --csv input option
parallel --header : \
	--env CAVIAR_OUTPUT_DIR,CAVIAR_OUTPUT_PREFIX,PARALLEL_JOB_COUNT \
	-j ${PARALLEL_JOB_COUNT} \
	--csv "CAVIAR -l ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/{phenotype_id}_LD.ld \
		-z ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/{phenotype_id}_zscore.txt \
		-o ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/RESULTS/{phenotype_id}_caviar -c 1" \
	:::: ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/CAVIAR_LD_${CAVIAR_OUTPUT_PREFIX}_calc.csv
# mkdir -p ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/RESULTS
# while IFS=',' read -r pheno chr
# do
# 	CAVIAR -l ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/${pheno}_LD.ld -z ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/${pheno}_zscore.txt -o ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/RESULTS/${pheno}_caviar -c 1
# done < <(tail -n +2 ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/CAVIAR_LD_${CAVIAR_OUTPUT_PREFIX}_calc.csv)