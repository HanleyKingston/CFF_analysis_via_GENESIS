#! /usr/bin/env Rscript
library(argparser)
library(magrittr)
argp <- arg_parser("Generate kinship plot") %>%
  add_argument("pcrelate_file", help = "PC-AiR file (.rds)") %>%
  add_argument("--out_prefix", help = "Prefix for output files",
               default = "") %>%
  add_argument("--group", help = "grouping variable") %>%
  add_argument("--king", flag = TRUE, help = "Is input file format from King?") %>%
  add_argument("--x_axis", default = "k0") %>%
  add_arguemnt("--y_axis", default = "kin") %>%
               
argv <- parse_args(argp)

library(ggplot2)
library(SNPRelate)
sessionInfo()
print(argv)

out_prefix <- argv$out_prefix

pcrel <- readRDS(argv$pcrelate_file)

kinship <- ifelse(argv$keep_king, snpgdsIBDSelection(KingRel), pcrel$kinBtwn)

png(paste0(out_prefix, "_kinship.png"))
ggplot(kinship, aes(k0, kin)) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color = "grey") +
    geom_point(alpha=0.2) +
    ylab("kinship estimate") +
    ggtitle("kinship")
dev.off()

#For KING plot:
KingRel <- readRDS(file = "king_obj.rds")


head(kinship)

png("king_plot.png")
ggplot(kinship, aes(IBS0, kinship)) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
    geom_point(alpha=0.2) +
    ylab("kinship estimate") +

dev.off()
