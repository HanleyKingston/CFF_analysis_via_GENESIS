
#Load phenotype data:
phenotype <- read.table("phenotype.txt", sep = "\t", header = TRUE)

##The null model requires a sample ID column:
phenotype$sample.id <- as.character(phenotype$sid)

##Read in PCA covariates
pca <- readRDS("CFF_LDsqrt0.1pcair.rds")
pcs.df <- as.data.frame(pca$vectors[,1:6])
colnames(pcs.df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")

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
  "cftr_addl_vars_wgs",
  "cftr_gt_category_wgs",
  "age_dx",
  "year_dx",
  "age at death",
  "knorma score",
  "VCF ID - don't use",
  "whether participant should be included in analysis based on QC - note: use keep_samples filter to exclude, not this column",
  "recruitment site - may have errors",
  "evaluation  method",
  "age when sampled",
  "estimated European subpopulation",
  "AD age-at-onset or censoring age",
  "count of APOE E2 alleles",
  "count of APOE E4 alleles",
  "count of APOE E3 alleles",
  "E4_carrier_one_or_two_coppies",
  "heterozygous of homozygous for E4",
  "sample ID - character of SUBJ_NO",
  "cox proportional hazard regression: survival",
  "deviance residuals",
  "first principal component - no LD-prune",
  "second principle component - no LD-prune",
  "third principle component - no LD-prune",
  "first principal component - LD pruned for R^2 = 0.2",
  "second principal component - LD pruned for R^2 = 0.2",
  "third principal component - LD pruned for R^2 = 0.2",
  "first principal component - LD pruned for R^2 = 0.1",
  "second principal component - LD pruned for R^2 = 0.1",
  "third principal component - LD pruned for R^2 = 0.1"
))


#Generate annotated dataframe
annot <- AnnotatedDataFrame(merged_phen, metadata)

#Save annot as an R object
saveRDS(annot, file = "annot.rds")
