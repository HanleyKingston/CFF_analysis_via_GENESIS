participants <- read.delim("participants_cffwgs.tsv", sep = "\t", header = TRUE)
key <- read.delim("key_cffwgs.tsv", sep = "\t", header = TRUE)
nrow(key)
#[1] 5199
flag <- scan("samples_flag.txt", "character", sep = "\n")
#Read 79 items
exclude <- scan("samples_exclude.txt", "character", sep = "\n")
#Read 51 items
gds.id <- scan("sample_id_gds.txt", "character", sep = "\n")
#Read 5134 items

"%notin%" <- Negate("%in%")

#Key has pid, sid, and vcf_id, participants only has pid... I will use the sid in analysis
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
not_shared_in_phenotype <- phenotype[phenotype$vcf_id %notin% gds.id, "vcf_id"]
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
identical(as.character(phenotype_pruned$vcf_id), gds.id)
#[1] TRUE

#find participants without necissary phenotype info
nrow(phenotype_pruned[is.na(phenotype_pruned$cftr_var_1_wgs) | is.na(phenotype_pruned$cftr_var_2_wgs),]) #Note: 609 people are missing a knorma value
#[1] 57
missing_F508 <- phenotype_pruned[is.na(phenotype_pruned$cftr_var_1_wgs) | is.na(phenotype_pruned$cftr_var_2_wgs), "vcf_id"] #Note: 609 people are missing a knorma value
sum(missing_F508 %in% exclude)
#[1] 14

phenotype_pruned$include_in_analysis <- ifelse(phenotype_pruned$vcf_id %in% duplicates | phenotype_pruned$vcf_id %in% exclude | phenotype_pruned$vcf_id %in% missing_F508,
                                               NA, ifelse(phenotype_pruned$vcf_id %in% flag, "flag", "keep"))

#Create a column for site (This may not be perfectly accurate because some individuals were included in multiple studies and some vcf_ids have changed)
phenotype_pruned$site  <- sub("_.*", "", phenotype_pruned$vcf_id)
table(phenotype_pruned[!is.na(phenotype_pruned$include_in_analysis),"site"])
# JHU  UNC   UW
#1841 1772 1358

colnames(phenotype_pruned)
nrow(phenotype_pruned)
nrow(phenotype_pruned[!is.na(phenotype_pruned$include_in_analysis),])
#[1] 5971

#Create column of deltaF508 count
phenotype_pruned$F508_count <- ifelse(phenotype_pruned$cftr_var_1_wgs == "F508del" & phenotype_pruned$cftr_var_2_wgs == "F508del", 2,
                                      ifelse(phenotype_pruned$cftr_var_1_wgs == "F508del" | phenotype_pruned$cftr_var_2_wgs == "F508del", 1, 0))

#Create a column of self-reported race:
phenotype_pruned$race_or_ethnicity <- NA

for(line in 1:nrow(phenotype_pruned)){
  if(rowSums(phenotype_pruned[line,c(6:12)], na.rm = TRUE) != 1){
      if(any(is.na(phenotype_pruned[line,c(6:12)]))){
        race_or_ethnicity <- NA
        }else if(phenotype_pruned[line,]$race_white == 1 & phenotype_pruned[line,]$race_natAm == 1){
          race_or_ethnicity <- "natAm" #Note: for purposes of plotting, people who are white and native american will be recorded as native american
        }else if(phenotype_pruned[line,]$race_white == 1 & phenotype_pruned[line,]$hispanic == 1){
          race_or_ethnicity <- "hispanic" #Note: for purposes of plotting, people who are white and native american will be recorded as native american
        }else{
          race_or_ethnicity <- "admixed_or_other"
          }                    
      }else if(phenotype_pruned[line,]$race_white == 1){
        race_or_ethnicity <- "white"
      }else if(phenotype_pruned[line,]$race_black == 1){
        race_or_ethnicity <- "black"
      }else if(phenotype_pruned[line,]$race_natAm == 1){
        race_or_ethnicity <- "natAm"
      }else if(phenotype_pruned[line,]$race_asian == 1){
        race_or_ethnicity <- "asian"
      }else if(phenotype_pruned[line,]$race_pac == 1){
        race_or_ethnicity <- "pac"
      }else if(!is.na(phenotype_pruned[line,]$hispanic) & phenotype_pruned[line,]$hispanic == 1){
        race_or_ethnicity <- "hispanic"
      }else{
        race_or_ethnicity <- "admixed_or_other"
        }
  phenotype_pruned$race_or_ethnicity[line] <- race_or_ethnicity
  }
              
table(phenotype_pruned$race_or_ethnicity)
#admixed_or_other            asian            black         hispanic
#              65               13               95              147
#           natAm            white
#              25             4750
        
#Plot count of deltaF508 per study site:
counts <- table(phenotype_pruned$F508_count, phenotype_pruned$site)

pdf("F508_count_by_study.pdf")
barplot(counts, main="Count of DeltaF508 per study Site",
  xlab="Study", legend = rownames(counts), beside = TRUE)
dev.off()

#To make more managable, I'm just selecting phenotypes I'm interested in
phenotype_pruned_selectCol  <- phenotype_pruned[,c("pid", "sid", "sex_wgs", "birthdate_year", "cftr_var_1_wgs", "cftr_var_2_wgs", "age_death", "knorma", "vcf_id", "include_in_analysis", "site", "F508_count", "race_or_ethnicity")]

saveRDS(phenotype_pruned_selectCol, "phenotype.rds")
write.table(phenotype_pruned_selectCol, "phenotype.txt", sep = "\t")

