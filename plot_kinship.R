#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

#Read arguments
argp <- arg_parser("Generate kinship plot") %>%
  add_argument("pcrelate or king file", help = "PC-Relate file (.rds)") %>%
  add_argument("--out_prefix", help = "Prefix for output files", default = "") %>%
  add_argument("--group", help = "grouping variable - a column in king or pc-relate dataframe") %>%
  add_argument("--is_king", flag = TRUE, help = "Is input file format from King?") %>%
  add_argument("--x_axis", default = "k0", help = "x variable (note: with --is_king flag, default is 'IBSO'") %>%
  add_argument("--y_axis", default = "kin", help = "y variable (note: with --is_king flag, default is 'kinship'")
               
argv <- parse_args(argp)
print(argv)

library(ggplot2)
library(SNPRelate)
sessionInfo()
print(argv)


rel <- readRDS(argv$pcrelate_file)
out_prefix <- argv$out_prefix


if(argv$is_king == TRUE){
	kinship.df <- snpgdsIBDSelection(rel)
	} else {
	kinship.df <- rel$kinBtwn
	}

if(argv$is_king == TRUE & argv$x_axis == "k0" & argv$y_axis == "kin"){
  argv$x_axis <- "IBS0"
  argv$y_axis <- "kinship"
  }

if (!is.na(argv$group)) {
  group <- argv$group
  } else {
  group <- NULL
}

png(paste0(argv$out_prefix, "_kinship.png"))
ggplot(kinship.df, aes_string(argv$x_axis, argv$y_axis, color = group)) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color = "grey") +
    geom_point(alpha=0.2) + #Note: if you get a "partial transparancy is not supported..." error remove "alpha" argument
    ylab("kinship estimate") +
    ggtitle("kinship")
dev.off()
