# get LD matrices for CAVIAR
plink --bfile ${bfile_prefix_path} --extract ${CAVIAR_OUTPUT}/${set_name}/${pheno}_variant_id.txt \
      --r2 square yes-really \
      --out ${CAVIAR_OUTPUT}/${set_name}/${pheno}_LD