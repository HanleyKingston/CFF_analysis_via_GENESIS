#Load phenotype data:
phenotype <- readRDS("phenotype.rds")

##The null model requires a sample ID column:
phenotype$sample.id <- as.character(phenotype$sid)

##Read in PCA covariates
pca <- readRDS("6-18pcair.rds")
pcs.df <- as.data.frame(pca$vectors[,1:8])
colnames(pcs.df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")

keep_samples <- readRDS("keep_samples.rds")
pcs.df$sample.id <- row.names(pcs.df)
#check:
identical(keep_samples, pcs.df$sample.id)
#[1] TRUE

gds.id <- readRDS("gds_id.rds")
#Alternatively, can extract from gds file with:
#library(SeqVarTools)
#gds.id <- seqGetData(seqOpen("CFF_sid_onlyGT.gds"), "sample.id")
##Check that phenotype sample ids match pc IDs (after applying filter to samples)
identical(gds.id, phenotype$sample.id)
#[1] TRUE
identical(gds.id[gds.id %in% keep_samples], pcs.df$sample.id)
#[1] TRUE

#Add PCA covariates to phenotype data by subject nomber
merged_phen <- merge(phenotype, pcs.df, by = "sample.id", all.x=TRUE)
sum(is.na(merged_phen$PC1))
#[1] 163

#order annotated df to match order of gds ids
merged_phen <- merged_phen[match(gds.id, merged_phen$sample.id),]
identical(merged_phen$sample.id, gds.id)
#[1] TRUE

##annotate data:
library(Biobase)
metadata <- data.frame(labelDescription = c(
  "sample.id - use in analyses",
  "participant ID",
  "subject ID - matches sample.id",
  "chromosomal sex",
  "birthdate year",
  "fist CFTR variant - for deltaF508 carriers, F508del is listed first",
  "second CFTR variant",
  "age at death",
  "knorma score",
  "VCF ID - don't use",
  "whether participant should be included in analysis based on QC - note: use keep_samples filter to exclude, not this column",
  "recruitment site - may have errors",
  "count of F508 deletion",
  "self-reported race or ethnicity - note, hispanic whites are reported as hispanic and white Native Americans as Native American",
  "does participant self-report as white - not: need to check if this is actually self-report",
  "does participant self-report as black",
  "does participant self-report as Native American",
  "does participant self-report as Asian",
  "does participant self-report as Pacific Islander",
  "does participant self-report as a race not included above",
  "does participant self-report as Hispanic",
  "first principal component",
  "second principle component",
  "third principle component",
  "fourth principle component",
  "fifth principle component",
  "sixth principle component",
  "seventh principle component",
  "eighth principle component"
))


#Generate annotated dataframe
annot <- AnnotatedDataFrame(merged_phen, metadata)

#Save annot as an R object
saveRDS(annot, file = "annot.rds")

library(dplyr)
head(pData(annot))
