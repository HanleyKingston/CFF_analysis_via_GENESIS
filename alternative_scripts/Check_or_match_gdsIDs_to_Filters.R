#WARNING:Soemthing withthe final check of same or different variant IDs is not working, needs revision

library(SeqVarTools)
library(dplyr)

flag.metric.df <- get(load("flag.metric.RData"))
head(flag.metric.df)
gds <- seqOpen("CFF_sid_onlyGT.gds")


        
# Flag metrics file *also* has variant.id - same number of variants but in different order. Same name of variable. Confusing.
# Want to map variant id in metrics file to variant id in gds.
res_list <- list()
for (chromosome in unique(flag.metric.df$chr)) {
  #Make a datframe of gds variant info
  seqSetFilterChrom(gds, chromosome)
  # Get variant info for that chromosome from the gds.
  variants_gds <- data.frame(
    position = seqGetData(gds, "position"),
    chromosome = seqGetData(gds, "chromosome"),
    variant_id_gds = seqGetData(gds, "variant.id"),
    stringsAsFactors = FALSE
  )
  # Subset flag metrics file to this chromosome.
  flag_metrics_chr <- flag.metric.df %>%
    filter(chr == chromosome)
  # Make sure that that no duplicate variants by chromosome and position
  stopifnot(all(!duplicated(flag_metrics_chr$pos)))
  stopifnot(all(!duplicated(variants_gds$position)))
  # Check that they have the same number of records
  stopifnot(nrow(variants_gds) == nrow(flag_metrics_chr))
  variant_map <- merge(flag_metrics_chr, variants_gds, by.x = "pos", by.y = "position")
  # Make sure no extra rows were added.
  stopifnot(nrow(variant_map) == nrow(flag_metrics_chr))
  res_list[[chromosome]] <- variant_map
  seqResetFilter(gds)
}
res <- bind_rows(res_list)


if(sum(res$variant.id != res$variant_id_gds) != 0)
   print("Original varaint ID's were not identical")


print(paste("number of variants that were not shared between files:", sum(res$variant.id != res$variant_id_gds)))
#[1] 105601274

#spot_check:
rand <- sample(1:nrow(res),1)
if(res[res$variant_id_gds == rand, "pos"] == seqGetData(gds, "position")[rand]){
           print("Matching WAS succesful")
}

#If original varaint ID's were not succesful, save the GDS ID's matched to the filter IDs as a new file AND run the "Generate_varaint_filter_when_gds_doesn't_match.R" script in alternative scripts
saveRDS(res, "flag.metric_with_gds_ids.rds")


