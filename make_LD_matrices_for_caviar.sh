# get LD matrices for CAVIAR
# get input variables
CAVIAR_OUTPUT_DIR=$1
CAVIAR_OUTPUT_PREFIX=$2
BFILE_PREFIX_PATH=$3
# prepare LD matrices with plink --extract --r2 square yes-really
mkdir -p ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/RESULTS
while IFS=',' read -r pheno chr
do
	plink --bfile ${BFILE_PREFIX_PATH} --extract ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/${pheno}_variant_id.txt \
      --r2 square yes-really \
      --out ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/${pheno}_LD
done < <(tail -n +2 ${CAVIAR_OUTPUT_DIR}/${CAVIAR_OUTPUT_PREFIX}/CAVIAR_LD_${CAVIAR_OUTPUT_PREFIX}_calc.csv)