library(SeqArray)

gds.id <- seqGetData(seqOpen("CFF_5134_onlyGT.gds"), "sample.id")
sample_key <- read.table("sample_names_key.txt", header = TRUE)

sum(duplicated(sample_key$vcf_id))
#[1] 0

sum(as.character(sample_key$vcf_id) %in% as.character(gds.id))
#[1] 5134

sample_key[match(sample_key$vcf_id,]
identical(as.character(sample_key[sample_key$vcf_id %in% gds.id, "vcf_id"]), as.character(gds.id))

sid_list <- as.vector(sample_key[sample_key$vcf_id %in% gds.id, "sid"])

class(sid_list)
#[1] "character"

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
