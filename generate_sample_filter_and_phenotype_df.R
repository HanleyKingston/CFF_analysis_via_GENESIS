"%notin%" <- Negate("%in%")

participants <- read.delim("participants_cffwgs.tsv", sep = "\t", header = TRUE)
nrow(participants)
#[1] 5161

sample_key <- read.table("sample_names_key.txt", header = TRUE)
sum(participants$pid %in% sample_key$pid)
#[1] 5161


#This only works because participants file has no duplicates... each dulpicated pid from sample_key will appear twice, so can filter after the fact
participants2 <- merge(participants, sample_key, by = "pid", all = TRUE)
nrow(participants2)
#[1] 5199

sum(duplicated(participants2$pid))
#[1] 38


#Read in sample_key and completely remove any samples without matching VCF_IDs
gds.id <- readRDS("gds_id.rds")
length(gds.id)
#[1] 5134
sum(participants2$sid %in% gds.id)

#Filtering sample key fist on sid ensures that the right partiicpant samples of the duplicated pids are included  in the participant phenotype df & order to match gds file:
phenotype_pruned <- participants2[match(gds.id, participants2$sid),]
nrow(phenotype_pruned)
#[1] 5134
identical(as.character(phenotype_pruned$sid), gds.id)

#Create include in analysis filter:
#Read in keep_samples by sid fitler (base on QC and duplicates...)
keep_samples <- scan("keep_samples.txt", "character", sep = "\n")
#Read 4966 items
sum(keep_samples %in% gds.id)
#[1] 4966
saveRDS(keep_samples, "keep_sample.rds")

phenotype_pruned$include_in_analysis <- ifelse(phenotype_pruned$sid %in% keep_samples, "include", "exclude")


#Create a column for site (This may not be perfectly accurate because some individuals were included in multiple studies and some vcf_ids have changed)
phenotype_pruned$site  <- sub("_.*", "", phenotype_pruned$vcf_id)
table(phenotype_pruned[phenotype_pruned$include_in_analysis == "include","site"])
# JHU  UNC   UW
#1831 1785 1350

colnames(phenotype_pruned)

#Create column of deltaF508 count
phenotype_pruned$F508_count <- ifelse(phenotype_pruned$cftr_var_1_wgs == "F508del" & phenotype_pruned$cftr_var_2_wgs == "F508del", 2,
                                      ifelse(phenotype_pruned$cftr_var_1_wgs == "F508del" | phenotype_pruned$cftr_var_2_wgs == "F508del", 1, 0))

#Create a column of self-reported race:
phenotype_pruned$race_or_ethnicity <- NULL

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
              
table(phenotype_pruned[phenotype_pruned$include_in_analysis == "include",]$race_or_ethnicity)
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



