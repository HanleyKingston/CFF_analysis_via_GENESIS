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
  head(variants_gds)
  # Subset flag metrics file to this chromosome.
  flag_metrics_chr <- flag.metric.df %>%
    filter(chr == chromosome)
  # Make sure that that no variant chromosome and position
  stopifnot(all(!duplicated(variants_gds$position)))
  # Check that they have the same number of records
  stopifnot(nrow(variants_gds) == nrow(flag_metrics_chr))
  variant_map <- flag_metrics_chr %>%
    left_join(variants_gds, by = c("chr" = "chromosome", "pos" = "position"))
  # Make sure no extra rows were added.
  stopifnot(nrow(variant_map) == nrow(flag_metrics_chr))
  res_list[[chromosome]] <- variant_map
  seqResetFilter(gds)
}
res <- bind_rows(res_list)


if(sum(res$variant.id != res$variant_id_gds) != 0)
   print("Original varaint ID's were not identical")


print(paste("number of variants that were not shared between files:" sum(res$variant.id != res$variant_id_gds)))
#[1] 105601274

#spot_check:
rand <- sample(1:nrow(res),1)
if(res[res$variant_id_gds == rand, "pos"] == seqGetData(gds, "position")[rand]){
           print("Matching WAS succesful")
}

saveRDS(res, "flag.metric_with_gds_ids.rds")


