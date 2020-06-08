#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

#Read arguments
argp <- arg_parser("Generate kinship plot") %>%
  add_argument("pcrelate_file", help = "PC-Relate file (.rds)") %>%
  add_argument("--out_prefix", help = "Prefix for output files", default = "") %>%
  add_argument("--group", help = "grouping variable") %>%
  add_argument("--is_king", flag = TRUE, help = "Is input file format from King?") %>%
  add_argument("--x_axis", default = "k0", help = "x variable") %>%
  add_argument("--y_axis", default = "kin", help = "y variable")
               
argv <- parse_args(argp)
print(argv)

library(ggplot2)
library(SNPRelate)
sessionInfo()
print(argv)


pcrel <- readRDS(argv$pcrelate_file)
out_prefix <- argv$out_prefix


if(argv$is_king == TRUE){
	kinship <- snpgdsIBDSelection(KingRel)
	} else {
	kinship <- pcrel$kinBtwn
	}

if(argv$is_king & argv$x_axis == "k0" & argv$y_axis == "kin"){
  argv$x_axis <- "IBS0"
  argv$y_axis <- "kinship"
  }

#if (!is.na(argv$group)) {
#  group <- argv$group
#  } else {
#  ## make up dummy group
#  group <- "group"
#  kinship$group <- "NA"
#}

png(paste0(argv$out_prefix, "_kinship.png"))
ggplot(kinship, aes_string(argv$x_axis, argv$y_axis, group = group)) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color = "grey") +
    geom_point(alpha=0.2) +
    ylab("kinship estimate") +
    ggtitle("kinship")
dev.off()
