participants <- read.delim("participants_cffwgs.tsv", sep = "\t", header = TRUE)
nrow(participants)
#[1] 5161

#Read in sample_key and completely remove any samples without matching VCF_IDs
gds.id <- readRDS("gds_id.rds")
nrow(participants[participants$VCF_ID %in% gds.id,])

#Filter based on sid_pid_to_keep list
sample_key <- read.table("sample_names_key.txt", header = TRUE)
sample_key <- sample_key[sample_key$vcf_id %in% gds.id,]
nrow(sample_key)
#[1] 5134

nrow(sample_key)
#[1] 5199

sum(participants$pid %in% sample_key$pid)
#[1] 5161

participants2 <- merge(participants, sample_key, by = "pid", all.x = TRUE, all.y = FALSE)
nrow(participants2)
#[1] 5199



#Save sample filter
keep_samples <- as.vector(sample_key$sid)


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

