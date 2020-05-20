#Identify identical twins & remove from analysis

#Read in PCrelate object and existing keep_samples vector
pcrel2 <- readRDS(file = "pcr_grm_LDsqrt0.1.rds")
kinship <- pcrel2$kinBtwn
keep_samples <- readRDS(file = "keep_samples.rds")

#Find likely identical twins
range(kinship$kin)
kinship[kinship$kin > 0.5, c("kin", "k2", "ID1", "ID2")]

#Remove identical twins from keep_samples vector
exclude_twins <- kinship[kinship$kin > 0.5, "ID1"] #ID1 is always the lowest of the 2 number IDs
"%notin%" <- Negate("%in%")
keep_samples <- keep_samples[keep_samples %notin% exclude_twins]

saveRDS(keep_samples, file = "keep_samples.rds")
