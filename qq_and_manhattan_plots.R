#R -q --vanilla --args assoc.rds --out_prefix CFF_F508 --ymin 1 --ymax 10 < assoc_plots.R &

library(qqman)
library(utils)
library(ggplot2)

argp <- arg_parser("Create qq-plot and Manhattan plot from association results") %>%
  add_argument("assoc_file", help = "file of association test results") %>%
  add_argument("--out_prefix", help = "Prefix for output files",
               default = "") %>%
  add_argument("--ymin", help = "minimum log(y value) to show on Manhattan plot", default = NULL) %>%
  add_argument("--ymax", help = "maximum log(y value) to show on Manhattan plot", default = NULL)


argv <- parse_args(argp)

assoc <- readRDS(argv$assoc_file)


#Generate lambda
lambda <- round(median(-2*log(assoc$Score.pval))/1.39, 3)
print(lambda)

#Generate qq-plot
png(paste0(argv$out_prefix, "_qq_plot.png"))
	qq(Assoc$Score.pval) #Need to add a title, but qq doesn't acept "main ="
  legend("topleft", legend = paste("lambda =", lambda))
  dev.off()
  
  
#generate Manhattan plot
Assoc <- type.convert(Assoc) #Chromosomes must be numeric for manhattan function

png(paste0(argv$out_prefix, "_manhattan_plot.png"))
	manhattan(Assoc, snp = "variant.id", chr = "chr", bp = "pos", p = "Score.pval", main = paste(title), ylim = ylim, annotatePval = 10^(-4))
  legend("topleft", legend = paste("lambda =", lambda))
  dev.off()


#Return a list of SNPs out of range based on ymax
#if(!is.null(ymax) & any(assoc$Score.pval <= 10^(-ymax)){
#  top_SNP <- Assoc[assoc$Score.pval <= 10^(-ymax), c("Score.pval", "variant.id")]
#  print("SNPs out of plot bounds:")
#  print(top_SNP)
