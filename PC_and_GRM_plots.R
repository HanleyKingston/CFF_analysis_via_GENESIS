#Generate PCs plots, percent variance explained (scree) plots, and relatedness plots

library(SeqArray)
library(SNPRelate)
library(SeqVarTools)
library(GENESIS)
library(MASS)
library(ggplot2)


#Read in phenotype and subset by keep_samples
phenotype <- read.table("phenotype.txt", header = TRUE)
keep_samples <- readRDS(file = "keep_samples.rds")
phenotype <- phenotype[phenotype$sid %in% keep_samples,]


pca <- readRDS("CFF_LDsqrt0.1pcair.rds")
pcs.df <- as.data.frame(pca$vectors[,1:7])
pcs.df$sid <- rownames(pcs.df)
pcs.df <- merge(pcs.df, phenotype[,c("sid","site", "race_or_ethnicity")], by="sid")

rels_V <- pca$rels
pcs.df$relate <- ifelse(pcs.df$sid %in% rels_V, "related", "unrelated")

pdf("/home/hkings/DATA/PC1and2.pdf")
ggplot() +
   geom_point(aes(pcs.df[,2], pcs.df[,3], col = factor(pcs.df$race_or_ethnicity), pch = pcs.df$relate)) +
   xlab("PC 1") + ylab("PC 2") +
   scale_shape(solid = FALSE) +
   labs(col= "Race or Ethnicity", shape="PC set")
dev.off()

pdf("/home/hkings/DATA/PC2and3.pdf")
ggplot() +
   geom_point(aes(pcs.df[,3], pcs.df[,4], col = factor(pcs.df$race_or_ethnicity), pch = pcs.df$relate)) +
   xlab("PC 2") + ylab("PC 3") +
   scale_shape(solid = FALSE) +
   labs(col= "Race or Ethnicity", shape="PC set")
dev.off()

pdf("/home/hkings/DATA/PC3and4.pdf")
ggplot() +
   geom_point(aes(pcs.df[,4], pcs.df[,5], col = factor(pcs.df$race_or_ethnicity), pch = pcs.df$relate)) +
   xlab("PC 3") + ylab("PC 4") +
   scale_shape(solid = FALSE) +
   labs(col= "Race or Ethnicity", shape="PC set")
dev.off()

#Check impact of site and PC relatedness assignment on results
pdf("/home/hkings/DATA/PC1and2.pdf")
ggplot() +
   geom_point(aes(pcs.df[,2], pcs.df[,3], col = factor(pcs.df$site), pch = pcs.df$relate)) +
   xlab("PC 1") + ylab("PC 2") +
   scale_shape(solid = FALSE) +
   labs(col= "Study", shape="PC set")
dev.off()

#Plot percent variance explained by each pc:
pca.df <- as.data.frame(pca$vectors[pca$unrels,])
var_prop_vect <- pca.df$varprop
PC_labs <- 1:length(var_prop_vect)

pdf("percent_var.pdf")
ggplot(dat, aes(x=PC_labs, y = 100*var_prop_vect)) +
    geom_point() + theme_bw() +
    xlab("PC") + ylab("percent variance accounted for")
dev.off()


#plotkinship estimates:
pcrel <- readRDS(file = "CFF_LDsqrt0.1pcr_obj.rds")
kinship <- pcrel$kinBtwn

#plot 2nd iteration kinship:
kinship$site  <- sub("_.*", "", kinship$ID1) #This may not be perfectly accurate because some individuals were included in multiple studies and some vcf_ids have changed
kinship$race_or_ethnicity  <- phenotype$race_or_ethnicity

png("kinship.png")
ggplot(kinship, aes(k0, kin, col = kinship$race_or_ethnicity, shape = factor(kinship$site))) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color = "grey") +
    geom_point(alpha=0.2) +
    ylab("kinship estimate") +
    ggtitle("kinship")
dev.off()


