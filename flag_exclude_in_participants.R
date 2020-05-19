d.delim("key_cffwgs.tsv", sep = "\t", header = TRUE)
nrow(key)
#[1] 5199
flag <- scan("samples_flag.txt", "character", sep = "\n")
#Read 79 items
exclude <- scan("samples_exclude.txt", "character", sep = "\n")
#Read 51 items
gds.id <- scan("sample_id_gds.txt", "character", sep = "\n")
#Read 5134 items

"%notin%" <- Negate("%in%")

sum(key$vcf_id %notin% gds.id)
#[1] 65
in_key_not_gds <- key[key$vcf_id %notin% gds.id,"vcf_id"]
#These should be removed!

sum(participants$pid %in% key$pid)
#[1] 5161 #All participants are in the key
sum(key$pid %in% participants$pid)
#[1] 5199 #All indiviudals in key are in participants, so 38 people must be repeated in the key
sum(duplicated(key$pid))
#[1] 38 #38 pairs of duplicates, so difference between key and participants file is just that duplicates were removed
#I only want to keep one of each duplicate... keep the person who is in the gds file
#Get all duplicates:
key[duplicated(key$pid) | duplicated(key$pid, fromLast = TRUE),]

#Because each duplicate has a different vcf_id and sid, I want to keep the one in the gds.id file
duplicates <- key[duplicated(key$pid) | duplicated(key$pid, fromLast = TRUE), "vcf_id"]
length(duplicates)
#[1] 72
sum(duplicates %in% gds.id)
#[1] 71
#all but one of the duplicates are in the gds file... since I can't determine which is the correct sample, I will remove them from the analysis entirely?
sum(in_key_not_gds %in% duplicates)
#[1] 1 #One of the IDs that are in the key but not the gds are in duplicates
sum(exclude %in% duplicates)
#[1] 0


#key[duplicated(key$pid) | duplicated(key$pid, fromLast = TRUE), c("sid", "vcf_id")] <- "NA"
nrow(key)
#[1] 5199
sum(key$vcf_id %notin% gds.id)
#[1] 65

#Add key and phenotype together so phenotype data contains the vcf_id
phenotype <- merge(participants, key, by = "pid", all = TRUE)

phenotype <- phenotype[order(phenotype$vcf_id),]
gds.id <- sort(gds.id)

sum(phenotype$vcf_id %notin% gds.id)
#[1] 65
not_shared_in_phenotype <- phenotype[phenotype$vcf_id %notin% gds.id,"vcf_id"]
sum(exclude %in% not_shared_in_phenotype)
#[1] 0
sum(flag %in% not_shared_in_phenotype)
#[1] 18


sum(gds.id %notin% phenotype$vcf_id)
#[1] 0 #There are no samples in the gds file without corresponding phenotype info

#anything in phenotype file without a corresponding entry in gds file should be removed, otherwise analysis won't run
phenotype_pruned <- phenotype[phenotype$vcf_id %in% gds.id,]
nrow(phenotype_pruned)
#[1] 5134
sum(phenotype_pruned$vcf_id %notin% gds.id)
#[1] 0



phenotype_pruned$site  <- sub("_.*", "", phenotype_pruned$vcf_id) #This may not be perfectly accurate because some individuals were included in multiple studies and some vcf_ids have changed
table(phenotype_pruned[!is.na(phenotype$pruned$include_in_analysis),"site"])
# JHU  UNC   UW
#1859 1774 1381


colnames(phenotype_pruned)

#To make more managable, I'm just selecting phenotypes I'm interested in
phenotype_pruned_temp  <- phenotype_pruned[,c("pid", "sex_wgs", "birthdate_year", "cftr_var_1_wgs", "cftr_var_2_wgs", "cftr_addl_vars_wgs", "cftr_gt_category_wgs", "age_dx", "year_dx", "age_death", "knorma", "vcf_id", "include_in_analysis")]

write.table(phenotype_pruned, "phenotype.txt", sep = "\t")
