# run CAVIAR on z score tables and LD matrices
# get input variables
CAVIAR_OUTPUT_DIR=$1
CAVIAR_OUTPUT_PREFIX=$2
# run CAVIAR
mkdir -p ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/RESULTS
while IFS=',' read -r pheno chr
do
	CAVIAR -l ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/${pheno}_LD.ld -z ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/${pheno}_zscore.txt -o ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/RESULTS/${pheno}_caviar -c 1
done < <(tail -n +2 ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/CAVIAR_LD_${CAVIAR_OUTPUT_PREFIX}_calc.csv)