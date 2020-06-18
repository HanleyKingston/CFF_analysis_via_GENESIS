library(SeqArray)

gds.id <- seqGetData(seqOpen("CFF_5134_onlyGT.gds"), "sample.id")
sample_key <- read.table("sample_names_key.txt", header = TRUE)

sum(duplicated(sample_key$vcf_id)) #Must be 0, otherwise, can't compair to gds file
#[1] 0

sample_key <- sample_key[match(gds.id, sample_key$vcf_id),]
nrow(sample_key)
#[1] 5134

identical(as.character(sample_key$vcf_id), as.character(gds.id)) #This must be TRUE!

sid_list <- as.vector(as.character(sample_key[sample_key$vcf_id %in% gds.id, "sid"]))

class(sid_list)
#[1] "character"

#Close gds file and reopen in write mode. Note: make a copy of the gds file and rename it "CFF_sid_onlyGT.gds"
gdsfmt::showfile.gds(closeall=TRUE)
rm(gds.id)
gds <- openfn.gds("CFF_sid_onlyGT.gds", readonly=FALSE)

add.gdsn(gds, "sample.id", sid_list, replace=TRUE, compress="LZMA_RA", closezip=TRUE)
closefn.gds(gds)

gds <- seqOpen("CFF_sid_onlyGT.gds")

#Check:
gds.id <- seqGetData(gds, "sample.id")
class(gds.id)
#[1] "character"
is.vector(gds.id)
#[1] TRUE
#check:
identical(gds.id, as.character(sample_key$sid)) #must be TRUE!

#Save a vector of gds IDs:
saveRDS(gds.id, "gds_id.rds")

seqClose(gds)
