

##Plot second iteration PCs
library(SeqArray)
library(SNPRelate)
library(SeqVarTools)
library(GENESIS)
library(MASS)
library(ggplot2)


#Read in phenotype and subset by keep_samples
phenotype <- read.table("phenotype.txt", header = TRUE)
keep_samples <- readRDS(file = "keep_samples.rds")
phenotype <- phenotype[phenotype$vcf_id %in% keep_samples,]
pca <- readRDS("pcair_LDsqrt0.1.rds")
pcs.df <- as.data.frame(pca$vectors[,1:4])



pcs.df$vcf_id <- rownames(pcs.df)
pcs.df <- merge(pcs.df, phenotype[,c("vcf_id","site")], by="vcf_id")

rels_V <- pca$rels
pcs.df$relate <- ifelse(pcs.df$vcf_id %in% rels_V, "related", "unrelated")

pdf("/home/hkings/DATA/PC1and2.pdf")
plot(pcs.df[,2], pcs.df[,3], col = factor(pcs.df$site), pch = pcs.df$relate, xlab="PC1", ylab="PC2")
dev.off()

pdf("/home/hkings/DATA/PC2and3.pdf")
plot(pcs.df[,3], pcs.df[,4], col = factor(pcs.df$site), pch = pcs.df$relate, xlab="PC2", ylab="PC3")
dev.off()

pdf("/home/hkings/DATA/PC3and4.pdf")
plot(pcs.df[,4], pcs.df[,5], col = factor(pcs.df$site), pch = pcs.df$relate, xlab="PC3", ylab="PC4")
dev.off()

#plot second iteration kinship estimates:
pcrel2 <- readRDS(file = "pcr_grm_LDsqrt0.1.rds")
kinship <- pcrel2$kinBtwn

#plot 2nd iteration kinship:
png("kinship.png")
ggplot(kinship, aes(k0, kin)) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color = "grey") +
    geom_point(alpha=0.5) + #alpha doesn't work with Xming
    ylab("kinship estimate") +
    ggtitle("kinship") +
    theme_bw()
dev.off()
#How to add color to this plot?

#Plot percent variance explained by each pc:
pcs.df2 <- as.data.frame(pca$vectors[pca$unrels,])
names(pcs.df2) <- paste0("PC", 1:ncol(pcs.df2))
pcs.df2$sample.id <- rownames(pcs.df2)
dat <- data.frame(pc=1:(ncol(pcs.df2)-1), varprop = pca$varprop)

pdf("percent_var.pdf")
ggplot(dat, aes(x=factor(pc), y = 100*varprop)) +
    geom_point() + theme_bw() +
    xlab("PC") + ylab("percent variance accounted for")
dev.off()


