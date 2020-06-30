
#!/bin/bash

for f in {1..22}; do
echo "Assoc_testing '$f' "
R -q --vanilla --args CFF_sid_onlyGT.gds annot.rds --out_prefix "sex_noPCs_chr_$f" --variant_id SNPS_bi_GATK_VQSR_MAF0.05_miss0.05.rds --sample_id keep_samples_noTwins.rds --chromosome "$f" --null_model sex_noPCs_nullmod.rds < assoc_test.R > sex_assoc_test_noPCs.log &

sleep 30
done

#Run with: sh sex_assoc_test_noPCs.sh
