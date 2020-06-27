#!/bin/bash

for f in {1..22}; do
echo "LD_pruning '$f' "
R -q --vanilla --args CFF_sid_onlyGT.gds --out_prefix "chr_$f" --sample_id keep_samples.rds --variant_id "pre_LD_SNP_filter.rds --maf 0.05 --missing 0.05 --window_size 1 --autosome_only TRUE --chromosome "$f" < ld_pruning.R > 6_26ld_pruning.log &
sleep 30
done

#Run with: sh ld_pruning.sh
