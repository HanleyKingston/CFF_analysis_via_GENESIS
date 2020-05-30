#Generate PCs plots, percent variance explained (scree) plots, and relatedness plots

library(SeqArray)
library(SNPRelate)
library(SeqVarTools)
library(GENESIS)
library(MASS)
library(ggplot2)


#Read in phenotype and subset by keep_samples
phenotype <- readRDS("phenotype.rds")
keep_samples <- readRDS(file = "keep_samples.rds")
phenotype <- phenotype[phenotype$sid %in% keep_samples,]


pca <- readRDS("CFF_LDsqrt0.1pcair.rds")
pcs.df <- as.data.frame(pca$vectors[,1:7])
pcs.df$sid <- rownames(pcs.df)
pcs.df <- merge(pcs.df, phenotype[,c("sid","site", "race_or_ethnicity")], by="sid", all.x = TRUE, all.y = FALSE)

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
var_prop_vect <- 100*pca$varprop[1:12]
PC_labs <- 1:12

pdf("percent_var.pdf")
plot(x=PC_labs, y = var_prop_vect, xlab = "PC", ylab = "percent variance accounted for")
dev.off()


#plotkinship estimates:
pcrel <- readRDS(file = "CFF_LDsqrt0.1pcr_obj.rds")
kinship <- pcrel$kinBtwn

png("kinship.png")
ggplot(kinship, aes(k0, kin)) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color = "grey") +
    geom_point(alpha=0.2) +
    ylab("kinship estimate") +
    ggtitle("kinship")
dev.off()

#For KING plot:
KingRel <- readRDS(file = "king_obj.rds")

kinship <- snpgdsIBDSelection(KingRel)
head(kinship)

pdf("king_plot.pdf")
ggplot(kinship, aes(IBS0, kinship)) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
    geom_point(alpha=0.2) +
    ylab("kinship estimate") +
    theme_bw()
dev.off()


plot(KingRel$IBSO, KingRel$kinship)

