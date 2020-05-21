library(SeqArray)

gds <- openfn.gds("CFF_pid.gds", readonly=FALSE)
#Note: I made a copy of my GDS file before changing the sample_ids, which is why the name doesn't match the original gds file

phenotype <- read.table("phenotype.txt", header = TRUE)

gds.id <- scan("sample_id_gds.txt", "character", sep = "\n")
#Read 5134 items

identical(as.character(phenotype$vcf_id), as.character(gds.id))
#[1] FALSE #No idea why this is false, but everything else I checked suggests they are identical

add.gdsn(gds, "sample.id", phenotype$pid, replace=TRUE, compress="LZMA_RA", closezip=TRUE)
closefn.gds(gds)

gds <- seqOpen("CFF_pid.gds")
head(seqGetData(gds, "sample.id"))

#check:
identical(as.character(seqGetData(gds, "sample.id")), as.character(phenotype$pid))

seqClose(gds)
