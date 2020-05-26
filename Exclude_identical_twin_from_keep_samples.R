#Identify identical twins & remove from analysis

#Read in PCrelate object and existing keep_samples vector
pcrel2 <- readRDS(file = "CFF_LDsqrt0.1pcr_obj.rds")
kinship <- pcrel2$kinBtwn
keep_samples <- readRDS(file = "keep_samples.rds")

#Find likely identical twins
range(kinship$kin)
kinship[kinship$kin > 0.40, c("kin", "k2", "ID1", "ID2")]

#Remove identical twins from keep_samples vector
exclude_twins <- kinship[kinship$kin > 0.5, "ID1"] #ID1 is always the lowest of the 2 number IDs
"%notin%" <- Negate("%in%")
keep_samples_noTwins <- keep_samples[keep_samples %notin% exclude_twins]

saveRDS(keep_samples_noTwins, file = "keep_samples_noTwins.rds")

