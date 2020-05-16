library(SeqArray)

gds_list = list.files(path = "/labdata12/CF_WGS2/Variants/CFF_5134_GDSs/seqArray_onlyGT/",
		      pattern = "onlyGT.gds")

gds_list = gds_list[-length(gds_list)]
gds_list = gds_list[-23]

seqMerge(gds.fn = gds_list, out.fn = "/home/hkings/DATA/CFF_5134_onlyGT.gds")
