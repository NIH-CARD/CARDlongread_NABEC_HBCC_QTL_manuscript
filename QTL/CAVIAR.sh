#!/bin/bash
# SBATCH --time 3:00:00
# SBATCH --mem=80g
# SBATCH --cpus-per-task=40

ml caviar 
ml plink/1.9

CAVIAR_OUTPUT=/data/CARDPB/data/NABEC/projects/QTL_paper_2024/newSV-eQTL/analysis/CAVIAR

set_name=$1
bfile_prefix_path=$2

mkdir -p ${CAVIAR_OUTPUT}/${set_name}/RESULTS
while IFS=',' read -r pheno chr
do
    plink --bfile ${bfile_prefix_path} --extract ${CAVIAR_OUTPUT}/${set_name}/${pheno}_variant_id.txt \
      --r2 square yes-really \
      --out ${CAVIAR_OUTPUT}/${set_name}/${pheno}_LD

    #run CAVIAR
    CAVIAR -l ${CAVIAR_OUTPUT}/${set_name}/${pheno}_LD.ld -z ${CAVIAR_OUTPUT}/${set_name}/${pheno}_zscore.txt -o ${CAVIAR_OUTPUT}/${set_name}/RESULTS/${pheno}_caviar -c 1
done < <(tail -n +2 ${CAVIAR_OUTPUT}/${set_name}/CAVIAR_LD_${set_name}_calc.csv)
