#!/bin/bash
#SBATCH --time 6:00:00
#SBATCH --mem=60g

genotools \
  --bfile NABEC_HBCC_merged_filtered \
  --out NABEC_HBCC_merged_filtered_ancestry_v2 \
  --ref_panel /data/GP2/utils/ref_dir/ref_panel_gp2_prune \
  --ancestry \
  --ref_labels /data/GP2/utils/ref_dir/ref_panel_ancestry.txt