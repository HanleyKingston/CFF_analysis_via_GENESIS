"%notin%" <- Negate("%in%")

participants <- read.delim("participants_cffwgs.tsv", sep = "\t", header = TRUE)
nrow(participants)
#[1] 5161

sum(duplicated(sample_key$pid))
#[1] 37
#Because sample key has duplicated pids, must remove the coroect duplicate (based on VCF) before mergeing with participants file based on pid
dups <- as.matrix(read.table("dups.txt", header = TRUE))
sample_key[sample_key$vcf_id %in% dups[,1:3], c("pid", "sid", "vcf_id")]
dups_df <- sample_key[sample_key$vcf_id %in% dups[,1:3] & sample_key$vcf_id %notin% dups[,5], c("vcf_id", "sid", "pid")]

sample_key <- sample_key[sample_key$vcf_id %notin% dups_df$vcf_id,]
sum(sample_key$sid %in% gds.id)

sample_key$pid <- ifelse(sample_key$vcf_id %in% dups_df$vcf_id, NA, as.character(sample_key$pid))
sum(is.na(sample_key$pid))
#[1] 37

sum(participants$pid %in% dups_df$pid)
#[1] 33

sum(participants$pid %in% sample_key$pid)
#[1] 5097

participants$pid <- ifelse(participants$pid %in% dups_df$pid, NA, as.character(participants$pid))
sum(is.na(participants$pid))
#[1] 33


#This only works because participants file has no duplicates
participants2 <- merge(participants, sample_key, by = "pid", all.x = FALSE, all.y = TRUE)
nrow(participants2)
#[1] 5134
sum(is.na(participants$pid))

sum(participants2$sid %in% gds.id)


#Read in sample_key and completely remove any samples without matching VCF_IDs
gds.id <- readRDS("gds_id.rds")
length(gds.id)
#[1] 5134

sample_key <- read.table("sample_names_key.txt", header = TRUE)

sample_key[duplicated(sample_key$pid),c("sid", "pid", "vcf_id")]
sample_key[duplicated(sample_key$pid,fromLast=TRUE),c("sid", "pid", "vcf_id")]


nrow(sample_key)
#[1] 5199
sum(sample_key$sid %in% gds.id)
#[1] 5134
sample_key <- sample_key[sample_key$sid %in% gds.id,]
#Filtering sample key fist on sid ensures that the right partiicpant samples of the duplicated pids are included  in the participant phenotype df
nrow(sample_key)
#[1] 5134



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
    
    
    
#Filter based on sid_pid_to_keep list
#drop_dups <- read.table("drop_dups.txt", header = TRUE)
#Check all drop identities are in sample_key and remove them
#sum(drop_dups$sid %in% sample_key$sid) == nrow(drop_dups)
#[1] TRUE
#sample_key_temp <- sample_key[sample_key$sid %notin% drop_dups$sid,]
#nrow(sample_key_temp)
#[1] 5162
#sum(duplicated(sample_key_temp$vcf_id)) #Must be 0!
#[1] 0


