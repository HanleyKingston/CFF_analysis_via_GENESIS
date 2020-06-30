#! /usr/bin/env Rscript
library(argparser)

argp <- arg_parser("Association plots")
argp <- add_argument(argp, "assoc_file",
                     help = "Association test results file (.rds)")
argp <- add_argument(argp, "--out_prefix", help = "Prefix for output files",
                     default = "")
argp <- add_argument(argp, "--qq_ylim", type="integer",  help = "y limit for qq plot (10^-ylim)")
argp <- add_argument(argp, "--Man_ylim", type="integer", help = "y limit for manhattan plot (10^-ylim)")
argv <- parse_args(argp)

sessionInfo()
print(argv)
library(dplyr)
library(ggplot2)
library(RColorBrewer)


assoc <- readRDS(argv$assoc_file)

lambda <- median(-2*log(assoc$Score.pval))/1.39

print(lambda)

## qq plot
n <- nrow(assoc)
x <- seq(n)
dat <- data.frame(obs=sort(assoc$Score.pval),
                  exp=x/n,
                  upper=qbeta(0.025, x, rev(x)),
                  lower=qbeta(0.975, x, rev(x)))

p <- ggplot(dat, aes(-log10(exp), -log10(obs))) +
    geom_line(aes(-log10(exp), -log10(upper)), color="gray") +
    geom_line(aes(-log10(exp), -log10(lower)), color="gray") +
    geom_point() +
    geom_abline(intercept=0, slope=1, color="red") +
    xlab(expression(paste(-log[10], "(expected P)"))) +
    ylab(expression(paste(-log[10], "(observed P)"))) +
    ylim(0, argv$qq_ylim) +
    ggtitle(paste("lambda =", format(lambda, digits=4, nsmall=3))) +
    theme_bw() +
    theme(plot.title = element_text(size = 22))
ggsave(paste0(argv$out_prefix, "qq.png"), plot=p, width=6, height=6)

rm(dat)


## manhattan plot
chr <- as.factor(assoc$chr) %>% levels()
cmap <- setNames(rep_len(brewer.pal(8, "Dark2"), length(chr)), chr)

# significance level
## genome-wide significance
signif <- c(5e-8, 5e-9, 1e-9)

p <- ggplot(assoc, aes(reorder(chr, sort(as.numeric(chr))),
                       -log10(Score.pval), group=interaction(chr, pos), color=chr)) +
    geom_point(position=position_dodge(0.8)) +
    scale_color_manual(values=cmap, breaks=names(cmap)) +
    geom_hline(yintercept=-log10(signif), linetype='dashed') +
    ylim(0, argv$Man_ylim) +
    theme_bw() +
    theme(legend.position="none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("chromosome") + ylab(expression(-log[10](p)))
ggsave(paste0(argv$out_prefix, "manh.png"), plot=p, width=10, height=5)
