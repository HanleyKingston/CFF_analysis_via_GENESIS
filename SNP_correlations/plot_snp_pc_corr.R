.libPaths("/labdata12/CF_WGS2/amstilp/analysis_pipeline_cff_wgs/R_library/")

library(argparser)
library(dplyr)
library(tidyr)
library(ggplot2)

argp <- arg_parser("Correlation of variants with PCs")
argp <- add_argument(argp, "--infile-prefix", help = "prefix for input PC-snp correlation files, eg <infile-prefix>_chr22.rds")
argp <- add_argument(argp, "--outfile", help = "output plot file")
argp <- add_argument(argp, "--npcs", default = 3, help = "number of pcs to plot")
argv <- parse_args(argp)
print(argv)

infile_prefix <- argv$infile_prefix
outfile <- argv$outfile
n_pcs <- argv$npcs

infile_pattern = sprintf("%s_chr\\d+?.rds", basename(infile_prefix))
files <- list.files(dirname(infile_prefix), pattern = infile_prefix)

corr_list <- lapply(files, function(x) readRDS(x))

# Put in chr/pos order.
corr <- bind_rows(corr_list) %>%
  mutate(chromosome = ordered(chromosome, levels = c(1:22, "X"))) %>%
  arrange(chromosome, position)

# Add variant filtering.

# Reshape data to long format for plotting.
dat <- corr %>%
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "corr") %>%
  mutate(
    abs_corr = abs(corr)
  ) %>%
  arrange(chromosome, position) %>%
  mutate(xpos = 1:n()) %>%
  filter(!is.na(corr))

# Thin points - up to 10,000 points from evenly spaced bins.
dat_thinned <- dat %>%
  group_by(PC) %>%
  mutate(bin = ntile(corr, 10)) %>%
  ungroup() %>%
  group_by(PC, bin) %>%
  sample_n(min(n(), 10000))

chrs <- sort(unique(dat_thinned$chromosome))
cmap <- setNames(rep(c("grey20", "grey50"), length.out = length(chrs)), chrs)
chr_labels <- dat_thinned %>% group_by(chromosome) %>% summarise(xpos = mean(xpos))

p <-
  ggplot(dat_thinned, aes( x = xpos, y = abs_corr)) +
  geom_point(aes(color = chromosome), size = 0.2) +
  scale_color_manual(values = cmap, breaks = names(cmap)) +
  facet_grid(rows = vars(PC)) +
  theme_bw() +
  ylim(c(0, 1)) +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
    #panel.grid.minor = element_blank()
  ) +
  # Chromosome labels at midpoint of the chromosome.
  scale_x_continuous(label = chr_labels$chromosome, breaks = chr_labels$xpos) +
  xlab("Chromosome") +
  ylab("abs(r)")


ggsave(outfile, plot = p)
