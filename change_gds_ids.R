library(SeqArray)

#Note: I made a copy of my GDS file before changing the sample_ids, which is why the name doesn't match the original gds file
gds <- openfn.gds("CFF_sid_onlyGT.gds", readonly=FALSE)

phenotype <- read.table("phenotype.txt", header = TRUE)

gds.id <- scan("sample_id_gds.txt", "character", sep = "\n")
#Read 5134 items


sum(as.character(phenotype$vcf_id) %in% as.character(gds.id))
#[1] 5134
phenotype <- phenotype[match(gds.id, phenotype$vcf_id),]
identical(as.character(phenotype$vcf_id), as.character(gds.id)) #This must be TRUE!
#[1] TRUE

add.gdsn(gds, "sample.id", phenotype$pid, replace=TRUE, compress="LZMA_RA", closezip=TRUE)
closefn.gds(gds)

gds <- seqOpen("CFF_sid_onlyGT.gds")
head(seqGetData(gds, "sample.id"))

#check:
identical(as.character(seqGetData(gds, "sample.id")), as.character(phenotype$pid)) #must be TRUE!

seqClose(gds)
