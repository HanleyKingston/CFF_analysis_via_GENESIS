library(SeqArray)

#Note: I made a copy of my GDS file before changing the sample_ids, which is why the name doesn't match the original gds file
gds <- openfn.gds("CFF_sid_onlyGT.gds", readonly=FALSE)

phenotype <- read.table("phenotype.txt", header = TRUE)

gds.id <- readRDS("sample_id_gds.rds")

sum(as.character(phenotype$vcf_id) %in% as.character(gds.id))
#[1] 5134

identical(as.character(phenotype$vcf_id), as.character(gds.id)) #This must be TRUE!
#[1] FALSE

#Can coerce to same order (if all the same IDs are shared) with this:
phenotype <- phenotype[match(gds.id, phenotype$vcf_id),]

sid <- as.character(phenotype$sid)

#Save re-ordered phenotype file:
write.table(phenotype, "phenotype.txt", sep = "\t")

add.gdsn(gds, "sample.id", sid, replace=TRUE, compress="LZMA_RA", closezip=TRUE)
closefn.gds(gds)

gds <- seqOpen("CFF_sid_onlyGT.gds")
head(seqGetData(gds, "sample.id"))

#check:
identical(seqGetData(gds, "sample.id"), as.character(phenotype$sid)) #must be TRUE!

seqClose(gds)
