#!/bin/bash

for f in {1..22}; do
echo "Assoc_testing '$f' "
R -q --vanilla --args CFF_sid_onlyGT.gds annot.rds 6_18pcr_mat.rds F508_count gaussian --out_prefix "$f" --covars "PC1 PC2 PC3" --variant_id var_filter_SNVs_MAF0.05.rds --sample_id keep_samples.rds --chromosome "$f" < assoc_test2.R &
sleep 30
done

#Run with: sh assoc_test.sh

