  
#!/bin/bash

for f in {1..22}; do
echo "Assoc_testing '$f' "
R -q --vanilla --args CFF_sid_onlyGT.gds annot.rds 6_18pcr_mat.rds sex_registry binomial --out_prefix "sex_$f" --covars "PC1 PC2 PC3 site age_cohort F508_count" --variant_id var_filter_SNVs_MAF0.05.rds --sample_id keep_samples_noTwins.rds --chromosome "$f" < assoc_test.R > sex_assoc_test.log &
sleep 30
done

#Run with: sh sex_assoc_test.sh
