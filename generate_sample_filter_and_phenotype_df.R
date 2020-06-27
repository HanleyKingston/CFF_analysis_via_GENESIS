"%notin%" <- Negate("%in%")

participants <- read.delim("participants_cffwgs.tsv", sep = "\t", header = TRUE)
nrow(participants)
#[1] 5161

sample_key <- read.table("sample_names_key.txt", header = TRUE)
sum(participants$pid %in% sample_key$pid)
#[1] 5161

sum(duplicated(participants$pid))
#[1] 0

sum(duplicated(sample_key$sid))
#[1] 0
sum(duplicated(sample_key$vcf_id))
#[1] 0
sum(duplicated(sample_key$pid))
#[1] 38


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
#[1] TRUE

#Create include in analysis filter:
#Read in keep_samples by sid fitler (base on QC and duplicates...)
keep_samples <- scan("keep_samples.txt", "character", sep = "\n")
#Read 4966 items
sum(keep_samples %in% gds.id)
#[1] 4966
saveRDS(keep_samples, "keep_samples.rds")

phenotype_pruned$include_in_analysis <- ifelse(phenotype_pruned$sid %in% keep_samples, "include", "exclude")


#Create a column for site (This may not be perfectly accurate because some individuals were included in multiple studies and some vcf_ids have changed)
phenotype_pruned$site  <- sub("_.*", "", phenotype_pruned$vcf_id)
table(phenotype_pruned[phenotype_pruned$include_in_analysis == "include","site"])
# JHU  UNC   UW
#1831 1785 1350

colnames(phenotype_pruned)

#Create column of deltaF508 count
table(phenotype_pruned$cftr_var_1) #F508del is by far the most common and does not seem to appear under any other name
table(phenotype_pruned$cftr_var_2)

phenotype_pruned$F508_count <- ifelse(phenotype_pruned$cftr_var_1_wgs == "F508del" & phenotype_pruned$cftr_var_2_wgs == "F508del", 2,
                                      ifelse(phenotype_pruned$cftr_var_1_wgs == "F508del" | phenotype_pruned$cftr_var_2_wgs == "F508del", 1, 0))

#Create a column of self-reported race (need to verify if this is self-report):
phenotype_pruned$race_or_ethnicity <- NULL

for(line in 1:nrow(phenotype_pruned)){
  if(rowSums(phenotype_pruned[line,c(6:12)], na.rm = TRUE) != 1){
      if(any(is.na(phenotype_pruned[line,c(6:12)]))){
        race_or_ethnicity <- NA
        }else if(phenotype_pruned[line,]$race_white == 1 & phenotype_pruned[line,]$race_natAm == 1){
          race_or_ethnicity <- "natAm" #Note: for purposes of plotting, because the dataset is majority white, people who are white and native american will be recorded as native american
        }else if(phenotype_pruned[line,]$race_white == 1 & phenotype_pruned[line,]$hispanic == 1){
          race_or_ethnicity <- "white_hispanic"
        }else if(phenotype_pruned[line,]$race_black == 1 & phenotype_pruned[line,]$hispanic == 1){
          race_or_ethnicity <- "black_hispanic"
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
#admixed_or_other            asian            black            natAm
#              64               13               93               25
#           white   white_hispanic
#            4603              139

        
#Plot count of deltaF508 per study site:
counts <- table(phenotype_pruned$F508_count, phenotype_pruned$site)

png("F508_count_by_study.png")
barplot(counts, main="Count of DeltaF508 per study Site",
  xlab="Study", legend = rownames(counts), beside = TRUE)
dev.off()

  
colnames(phenotype_pruned)

    
#Create birthday 5-yr age cohorts:
range(phenotype_pruned[phenotype_pruned$include_in_analysis == "include", "birthdate_year"], na.rm = TRUE)
#[1] 1943 2011
phenotype_pruned$age_cohort <- cut(phenotype_pruned$birthdate_year, 14)

    
#Count NA's for covariates in study (should these be included in analyses?)    
table(phenotype_pruned[phenotype_pruned$include_in_analysis == "include", "sex_registry"], useNA = "ifany")
#   F    M <NA>
#2322 2634   10
    
table(phenotype_pruned$age_cohort, useNA = "ifany")
#          4          14          51          94         147         179
#(1972,1977] (1977,1982] (1982,1987] (1987,1992] (1992,1996] (1996,2001]
#        231         303         609         682         919        1084
#(2001,2006] (2006,2011]        <NA>
#        701          83          33

sum(is.na(phenotype_pruned[phenotype_pruned$include_in_analysis == "include", "birth_cohort"]))
#[1] 0
    
sum(is.na(phenotype_pruned[phenotype_pruned$include_in_analysis == "include","site"]))
#[1] 0

#Sex_registry needs to be 0-1 for association testing
unique(phenotype_pruned$sex_registry)
phenotype_pruned$sex_registry <- ifelse(phenotype_pruned$sex_registry == "F", 0, 1)

#To make more managable, I'm just selecting phenotypes I'm interested in
phenotype_pruned_selectCol  <- phenotype_pruned[,c("pid", "sid", "sex_wgs", "sex_registry", "birthdate_year", "cftr_var_1_wgs", "cftr_var_2_wgs", "age_death", "knorma", "vcf_id", "include_in_analysis", "site", "F508_count", "race_or_ethnicity", "race_white", "race_black", "race_natAm", "race_asian", "race_pac", "race_other", "hispanic", "age_cohort")]
    
saveRDS(phenotype_pruned_selectCol, "phenotype.rds")



