
#Load phenotype data:
phenotype <- readRDS("phenotype.rds")

##The null model requires a sample ID column:
phenotype$sample.id <- as.character(phenotype$sid)

##Read in PCA covariates
pca <- readRDS("CFF_LDsqrt0.1pcair.rds")
pcs.df <- as.data.frame(pca$vectors[,1:6])
colnames(pcs.df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")

keep_samples <- readRDS("keep_samples.rds")
pcs.df$sample.id <- row.names(pcs.df)
identical(as.character(phenotype[phenotype$sample.id %in% keep_samples, "sample.id"]), gds.id)

##Check that phenotype sample ids match pc IDs (after applying filter to samples)
#phenotype <- phenotype[match(pcs.df$sample.id, phenotype$sample.id),]
identical(as.character(phenotype[phenotype$sample.id %in% keep_samples, "sample.id"]), pcs.df$sample.id)


#Add PCA covariates to phenotype data by subject nomber
merged_phen <- merge(phenotype, pcs.df, by = "sample.id", all.x=TRUE)

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
  "first principal component",
  "second principle component",
  "third principle component",
  "fourth principle component",
  "fifth principle component",
  "sixth principle component"
))


#Generate annotated dataframe
annot <- AnnotatedDataFrame(merged_phen, metadata)

#Save annot as an R object
saveRDS(annot, file = "annot.rds")

library(dplyr)
head(pData(annot))
