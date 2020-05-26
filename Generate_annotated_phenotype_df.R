
#Load phenotype data:
phenotype <- read.table("phenotype.txt", sep = "\t", header = TRUE)

##The null model requires a sample ID column:
phenotype$sample.id <- as.character(phenotype$sid)

##Check that columns match gds file
library(SeqVarTools)
gdsfmt::showfile.gds(closeall=TRUE)
gds <- seqOpen("CFF_sid_onlyGT.gds")
gds.id <- seqGetData(gds, "sample.id")
sum(gds.id %in% phenotype$sample.id) == sum(phenotype$sample.id %in% gds.id) #phenotype file isn't fitlered so will have more total samples
keep_samples <- readRDS("keep_samples.rds")

#phenotype <- phenotype[match(gds.id, phenotype$sample.id),]

identical(as.character(phenotype[phenotype$sample.id %in% keep_samples, "sample.id"]), gds.id)

pca <- readRDS("CFF_LDsqrt0.1pcair.rds")
pcs.df <- as.data.frame(pca$vectors[,1:6])
colnames(pcs.df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")

pcs.df$sample.id <- row.names(pcs.df)


#Add PCA covariates to phenotype data by subject nomber
merged_phen <- merge(merged_phen, pca_noLD, by = "SUBJ_NO", all.x=TRUE, all.y = FALSE)
merged_phen <- merge(merged_phen, pca_LDsqrt0.2, by = "SUBJ_NO", all.x=TRUE, all.y = FALSE)
merged_phen <- merge(merged_phen, pca_LDsqrt0.1, by = "SUBJ_NO", all.x=TRUE, all.y = FALSE)


##annotate data:
library(Biobase)
metadata <- data.frame(labelDescription = c(
  "family number",
  "subject number",
  "subject's father",
  "subject's mother",
  "subject's sex",
  "case-control status",
  "age at last evaluation",
  "age at death",
  "self-reported race",
  "self-reported Hispanic status",
  "estimated ethnicity",
  "age at onset of dementia",
  "age at diagnosis of dementia",
  "first APOE allele",
  "second APOE allele",
  "recruitment site",
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
