## ---------------------------
## Script name: genome_vis.R
##
## Purpose of script: Calculate relative depth and SNP density for a given sample,
## and then plot a genome-scale view.
##
## Author: Nancy Scott
##
## Date Created: 2023-03-29
##
## Email: scot0854@umn.edu
## ---------------------------
## Notes: Adapted from https://github.com/berman-lab/ymap
## and https://github.com/stajichlab/C_lusitaniae_popseq.
## Inputs are tab-delim files from samtools depth and samtools mpileup, with headers.
## This script uses those headers, so change carefully.
## For samtools depth, recommend using bam files that have been corrected for gc bias (optional, reduces copy number noise).
## Script order is candida_gc_correct.sh -> candida_ymap.sh (uses berman_count_snps_v5.py) -> genome_vis.R
## Can run this R script from candida_ymap.sh or interactively, see below.
## ---------------------------
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# Input file variables
genome_df_file <-  args[1]
sample_id <- args[2] # or "YourID"
feature_file <- args[3] # or "path/to/features.txt"
label_file <- args[4] # or "path/to/chr_labels.txt"

## ---------------------------
# Load packages
library(readxl) 
library(tidyverse)
library(ggplot2)

## ---------------------------
# Set variables
window <- 5000 # size of window used for rolling mean and snp density
ploidy <- 2

y_axis_labels <- c(1,2,3,4)  # manual y-axis labels, adjust as needed
inter_chr_spacing <- 150000 # size of space between chrs

save_dir <- paste0("images/Calbicans/",Sys.Date(),"_genome_plots/") # path with trailing slash, or  "" to save locally 
ref <- "sc5314" # short label for file name or "" to leave out

# Plotting variables
# X-axis labels overwrite input scaffold names in final plot
chr_ids <- scan(label_file, what = character())

snp_low <- "white"  # snp LOH colors, plot function uses 2-color gradient scale
snp_high <- "black"  # snp LOH colors, plot function uses 2-color gradient scale
cnv_color <- "dodgerblue4"  # copy number color
    
ploidy_multiplier <- 2  # this multiplied by ploidy sets the max-y scale

chrom_outline_color <- "gray15"  # color of chromosome outlines
    
chrom_line_width <- 0.2  # line width of chromosome outlines

## ---------------------------
# Final dataframe of joined copy number, snps, and plotting positions per window
genome_depth <- read_xlsx(genome_df_file, sheet=1)

## ---------------------------
# Small dataframes for chrom. outlines and features
chroms <- genome_depth %>%
  group_by(index) %>%
  summarise(xmin=min(plot_pos), xmax=max(plot_pos), ymin=0, ymax=Inf)

features <- read_tsv(feature_file, show_col_types = FALSE)

features <- features %>% 
     group_by(chr, index=consecutive_id(chr)) %>% 
     left_join(chroms, by=join_by(index)) 
 
features <- features %>% 
     mutate(plot_start = start + xmin, plot_end = end + xmin)

# Tick marks to center chromosome ID label
ticks <- tapply(genome_depth$plot_pos, genome_depth$index, quantile, probs =
                0.5, na.remove = TRUE)

## ---------------------------
# Plot linear genome
p <- ggplot(genome_depth) +
  scale_color_gradient(low=snp_low,high=snp_high, na.value = "white", guide = "none") +
  geom_segment(aes(x = plot_pos, y = 0, color = snp_count, xend = plot_pos, yend = Inf)) +
  geom_segment(aes(x = plot_pos, 
                   y = ifelse(copy_number <= ploidy*ploidy_multiplier, copy_number, Inf),
                   xend = plot_pos, yend = ploidy), alpha = 0.9, color = cnv_color) +
  geom_rect(data=chroms, aes(group=index, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            linewidth = chrom_line_width, fill = NA, 
            colour = chrom_outline_color, linejoin = "round", inherit.aes = FALSE) +
  geom_point(data = features, size = 2,
               aes(group=index, x=plot_start, y=ymin, shape = Feature, fill = Feature),
               position = position_nudge(y=0.07)) +
    scale_fill_manual(values = c("white", "grey26", "deepskyblue")) +
    scale_shape_manual(values = c(24,21,22)) +
  ylab(sample_id) +
  scale_x_continuous(name = NULL, expand = c(0, 0), breaks = ticks, labels=chr_ids) +
  scale_y_continuous(limits = c(0, ploidy*ploidy_multiplier), breaks = y_axis_labels) +
  theme_classic() +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.ticks = element_line(color = NA),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size=12))

## ---------------------------
# Save plot
ggsave(sprintf("%s%s_%s_%s_%sbp.png", save_dir, Sys.Date(), sample_id, ref, window),
       p, width = 18, height = 1.7, units = "in", device = png, dpi = 300, bg = "white")
