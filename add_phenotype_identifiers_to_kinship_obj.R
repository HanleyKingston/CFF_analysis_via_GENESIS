args <- commandArgs(trailingOnly = TRUE)

library(ggplot2)
library(SNPRelate)

#Read in phenotype and subset by keep_samples
phenotype <- readRDS(args[1])
keep_samples <- readRDS(file = "keep_samples.rds")
phenotype <- phenotype[phenotype$sid %in% keep_samples,]

pcrel <- readRDS(file = args[2])
kinship <- pcrel$kinBtwn

for(i in 1:nrow(kinship)){
  RE1 <- phenotype[phenotype$sid == kinship$ID1[i], "race_or_ethnicity"]
  RE2 <- phenotype[phenotype$sid == kinship$ID2[i], "race_or_ethnicity"]
  kinship$race_or_ethnicity[i] <- ifelse(RE1 == "white" & RE2 == "white", "white", ifelse(RE1 == "black" & RE2 == "black", "black", ifelse(RE1 == RE2, "match", "different_race_or_ethnicity")))
  if(i %% 100000 == 0)
     print("1,000,000 rows run")
     }
  }

saveRDS(kinship, "kinship_phenotype_added.rds")
