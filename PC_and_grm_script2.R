#Arguments:
#1. gds_file: the file path to the gds file (with .vcf.gds extension)
#2. LD-pruning R object (a list of variants to incldue)
#3. keep_variants: a file path to a list of variants to keep (must match corresponding rownames in phenotype and gds - saved as an R object
#4. keep_samples: a file path to a list of samples to keep (must match corresponding smaple IDs in phenotype and gds files - default is familyID_SUBJID) - saved as an R object
#5. text to uniquely identify plots and figures

#ex. run with command: Rscript PC_and_grm_script2.R CFF_5134_onlyGT.gds "pruned.rds" "keep_var_stringent.rds" "keep_samples.rds" "LDsqrt0.1"

#to test: args <- c("CFF_5134_onlyGT.gds", "0.316227766", "keep_var_stringent.rds", "keep_samples.rds", "LDsqrt0.1")

library(SeqArray)
library(SNPRelate)
library(SeqVarTools)
library(GENESIS)
library(MASS)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
print(args)

##Open GDS file:
gdsfmt::showfile.gds(closeall=TRUE) #makes sure gds file isn't already open - not really necissary in a script
gdsfile <- args[1]
gds <- seqOpen(gdsfile)
#Note: can access gds with (ex.): seqGetData(gds, "sample.id")

pruned <- readRDS(file = args[2])
keep_vars <- readRDS(file = args[3])
keep_samples <- readRDS(file = args[4])
text <- args[5]

cat(
"This generates 2 files and 4 plots:
Files
  PC: pcs_",text,".txt
  GRM: grm_",text,".txt (for no LD-pruning, also generates the PC-niave King-robust grm: King_robust_grm)
Plots:
  PCs:  PCs1and2_",text,".pdf, PCs2and3_",text,".pdf, PCs3and4_",text,".pdf
  Kinship:  kinship_",text,".png (K0 (IBS0) by kinship)

Where ",text, " files/plots used an LD threshold of ",text," for correlation (method=\"corr\")
All plots and figures are from running PC-air and PC-relate twice (ie. 2 iterations) with the following inputs:
  LD threshold for LD-pruning: ",text,"
  kinship threshold (kin.thresh in pcair) = 2^(-11/2) =~0.022 (used to subest into related and unrelated groups for PC calcualtions)
  div.thresh in pcair = (-2)^(-11/2)
  scaleKin in pcrelateMakeGRM = 2
  thresh in pcrelateToMatrix = 2^(-11/2) (This is threshold below which kinship values in the GRM are set to 0, this is for computational efficiency)
  500000 bp is the window used for LD pruning
  ", sep = "")


#Test
pdf("test_plot.pdf")
plot(c(1,2),c(1,3), main = "I am a test to make sure PC_and_grm_script.R is executing")
dev.off()

#Generate King reatedness matrix
print(paste("LD applied at threshold:", text))
#get initial relatedness estimates with KING:
king <- snpgdsIBDKING(gds, snp.id=intersect(pruned, keep_vars), sample.id = keep_samples)
#Working space: 3,547 samples, 109,292 SNVs

#Put kinship estimates in a matrix
kingMat <- king$kinship
colnames(kingMat) <- rownames(kingMat) <- king$sample.id


##Develop PC-matrix with PC-Air
print(paste("1st iteration PC-Air: LD applied at threshold:", text))
mypcair <- pcair(gds, kinobj = kingMat, kin.thresh=2^(-11/2), divobj = kingMat, snp.include = intersect(pruned, keep_vars)) #This includes the filtered set of samples without having to give it the keep_smaples object

#Save a copy of 1st iteration pc object:
saveRDS(mypcair, paste("pcair_", text, "_temp.rds", sep = ""))


##Generate relatedness estimation adjusted for PCs using PC_Relate
#filter gds file by keep_samples and keep_vars vectors and (for LD pruning) by pruned SNPs
seqResetFilter(gds)
seqSetFilter(seqData, variant.id=keep_vars, sample.id=keep_samples)
seqData <- SeqVarData(gds)

print("1st iteration PC-relate")
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
mypcrel <- pcrelate(iterator, pcs=mypcair$vectors[,1:3], training.set=mypcair$unrels)
pcrelate_matrix <- pcrelateToMatrix(mypcrel, scaleKin=2, thresh = 2^(-11/2)) #default = NULL, 2^(-11/2) suggested by bioconductor... should this be NULL for the first iteration???, scaleKin multiplies by 2 so average person's kinship with self is 1
#note:thresh appears to be applies after scalekin
#For some reason, there is one pair of individuals that are not coerced to 0 despite a kinship estimate of 3.689586e-07 in the LD pruning at sqrt(0.1) case... no idea why!

#Save a copy of 1st iteration PCrelate object
saveRDS(pcrel2, paste("pcr_grm_",text,"_temp.rds", sep = ""))


#Generate 2nd iteration PCs
print(paste("2nd iteration PC-Air: LD applied at threshold:", text))
pca <- pcair(seqData, kinobj=pcrelate_matrix, kin.thresh=2^(-11/2), divobj=kingMat, snp.include= intersect(pruned, keep_vars))

#Generate 2nd iteration kinship estimates
resetIterator(iterator, verbose = TRUE) #This doesn't seem necissary here. Included just in case???

print("2nd iteration PC-relate")
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
pcrel2 <- pcrelate(iterator, pcs=pca$vectors[,1:3], training.set=pca$unrels)

pcrelate_matrix <- pcrelateToMatrix(pcrel2, scaleKin=2, thresh = 2^(-11/2)) #Thresh is the threshold below which kinship values are coerced to 0

#Save entire pc object
saveRDS(mypcair, paste("pcair_", text, ".rds", sep = ""))

#save eigenvectors 1-3
pcs <- pca$vectors
write.table(pcs[,1:3], file = paste("pcs_",text,".txt", sep = ""), sep = "\t")

#Save entire relatedness object:
saveRDS(pcrel2, paste("pcr_grm_",text,".rds", sep = ""))

#Save GRM matrix
write.matrix(pcrelate_matrix, file = paste("grm_",text,".txt", sep = ""), sep=",")

#To read in, use:
#myGRM <- as.matrix(read.csv("grm_LDsqrt0.2.txt", header=T, na.strings="NA"))


##Plot second iteration PCs
#Read in merged_phen2 and subset by keep_samples
merged_phen2 <-read.table("merged_phen2", header=TRUE)
merged_phen2$SUBJ_NO <- merged_phen2$SUBJ_NO<-paste(merged_phen2$X.FAM_NO,merged_phen2$SUBJ_NO, sep = "_")
merged_phen2 <- merged_phen2[merged_phen2$SUBJ_NO %in% keep_samples,]
pcs.df <- as.data.frame(pcs)
pcs.df$SUBJ_NO <- rownames(pcs)
pcs.df <- merge(pcs.df, merged_phen2[,c("SUBJ_NO","CaucCode")], by="SUBJ_NO")

rels_V <- mypcair$rels
pcs.df$relate <- ifelse(pcs.df$SUBJ_NO %in% rels_V, "related", "unrelated")

pdf(paste("~/raw_data_analysis/plots_and_charts/PCs1and2_",text,".pdf", sep = ""))
plot(pcs.df[,2], pcs.df[,3], col=pcs.df$CaucCode, pch = pcs.df$relate, xlab="PC1", ylab="PC2", main = text)
dev.off()

pdf(paste("~/raw_data_analysis/plots_and_charts/PCs2and3_",text,".pdf", sep = ""))
plot(pcs.df[,3], pcs.df[,4], col=pcs.df$CaucCode, pch = pcs.df$relate, xlab="PC2", ylab="PC3", main = text)
dev.off()

pdf(paste("~/raw_data_analysis/plots_and_charts/PCs3and4_",text,".pdf", sep = ""))
plot(pcs.df[,4], pcs.df[,5], col=pcs.df$CaucCode, pch = pcs.df$relate, xlab="PC3", ylab="PC4", main = text)
dev.off()

#plot second iteration kinship estimates:
kinship <- pcrel2$kinBtwn

#plot 2nd iteration kinship:
png(paste("~/raw_data_analysis/plots_and_charts/kinship_",text,".png", sep = ""))
ggplot(kinship, aes(k0, kin)) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
    #geom_point(alpha=0.5) + #alpha doesn't work with Xming
    ylab("kinship estimate") +
    ggtitle(text) +
    theme_bw()
dev.off()

#Add: calculate % var explained by each PC
