library(SeqArray)

gds.id <- seqGetData(seqOpen("CFF_5134_onlyGT.gds"), "sample.id")
sample_key <- read.table("sample_names_key.txt", header = TRUE)
nrow(sample_key)
#[1] 5199

sum(duplicated(sample_key$vcf_id)) #Must be 0, otherwise, can't compair to gds file
#[1] 0

sample_key <- sample_key[match(gds.id, sample_key$vcf_id),]
nrow(sample_key)
#[1] 5134

identical(as.character(sample_key$vcf_id), as.character(gds.id)) #This must be TRUE!

sid_list <- as.character(sample_key$sid)

class(sid_list)
#[1] "character"

#Spot check:
sid_list[1200]
#[1] "S38756"
gds.id[1200]
#[1] "JHU_TSS_592_2003_1"
sample_key[1200,]
#        pid    sid             vcf_id
#4203 P82544 S38756 JHU_TSS_592_2003_1


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
