#!/bin/bash

for f in {1..22}; do
echo "Assoc_testing '$f' "
R -q --vanilla --args CFF_sid_onlyGT.gds annot.rds 6_18pcr_mat.rds F508_count poisson --out_prefix "$f" --covars "PC1 PC2 PC3 site" --variant_id var_filter_SNVs_MAF0.05.rds --sample_id keep_samples_noTwins.rds --chromosome "$f" < assoc_test.R > assoc_test.log &
sleep 30
done

#Run with: sh assoc_test.sh

