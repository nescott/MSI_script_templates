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
## Notes: Adapted significantly from https://github.com/berman-lab/ymap
## and https://github.com/stajichlab/C_lusitaniae_popseq.
## Inputs are tab-delim files from samtools depth and samtools mpileup, with headers.
## This script uses those headers, so change carefully.
## For samtools depth, recommend using bam files that have been corrected for gc bias (optional, reduces copy number noise).
## Script order is gc_correct.sh -> ymap_prep.sh (uses berman_count_snps_v5.py) -> genome_vis.R
## ---------------------------
## load packages
library(RcppRoll)
library(tidyverse)
library(ggplot2)
library(writexl)
## ---------------------------
## input file variables
read_depth_file <- "depth/AMS5643_sc5314_gccorrected_depth.txt"
snp_file <- "alleles/AMS5643_sc5314_putative_SNPs.txt"
sample_id <- "AMS5643"
mito <- "Ca19-mtDNA" #scaffold ID - will use for subsetting

## data-wrangling variables
window <- 5000 # size of window used for rolling mean and snp density
ploidy <- 2 # for proper plotting of y-axis

## plotting variables
chr_ids <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "ChrR") # manual x-axis labels overwrite input scaffold names in final plot
y_axis_labels <- c(1,2,3,4)  # manual y-axis labels for now, to remove "0" from axis if wanted
inter_chr_spacing <- 150000 # size of space between chrs
snp_low <- "white"  # snp LOH colors, plot function uses 2-color gradient scale
snp_high <- "black"  # snp LOH colors, plot function uses 2-color gradient scale
copy_number <- "steelblue4"  # copy number color
ploidy_multiplier <- 2  # this number multiplied by ploidy sets the max-y scale
chrom_outline_color <- "gray15"  # color of chromosome outlines
chrom_line_width <- 0.2  # line width of chromosome outlines

## output variables
save_dir <- "plots/" # path with trailing slash, or just "" to save in same folder
ref <- "SC5314_a21" # short label for generating file name 

## ---------------------------
## Base R doesn't have a mode calculation
# nice to compare this to median but not essential
# (this function is used to make chr_mode column below, delete if necessary)
Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

## ---------------------------
## relative copy number calcs from samtools depth input
genome_raw <- read.table(read_depth_file, header = TRUE)
genome_raw <- genome_raw %>%
  filter(chr != mito)
genome_raw$rolling_mean <- roll_mean(genome_raw$depth, window)[seq_len(length(genome_raw$chr))]

raw_genome_median <- median(genome_raw$depth) #includes all chromosomes, may need correcting

chr_median <- genome_raw %>%  # checking each chromosome for outliers relative to genome
  group_by(chr) %>%
  summarise(chr_mode = Modes(depth), chr_med = median(depth))  # can manually compare mode and median if questioning median

subset_chr_median <- chr_median %>%  # moderate filtering to avoid aneuploidy skew of "normal" genome depth
  filter(chr_med <= raw_genome_median *1.15 & chr_med >= raw_genome_median * 0.85)

genome_median <- median(subset_chr_median$chr_med)  # filtered median used to calculate relative depth

genome_window <- genome_raw %>%
  group_by(chr, index=consecutive_id(chr)) %>%
  reframe(position=unique((pos %/% window)*window+1))

genome_window <- genome_window %>%
  group_by(index) %>%
  mutate(chr_length = max(position))

chrs <- as.vector(unique(genome_window$chr_length))  # for plotting
chr_plot <- c()
for(i in 1:length(chrs)){chr_plot[i] <- sum(chrs[1:i-1])}  # plotting

genome_depth <- genome_raw %>%
  group_by(chr, index=consecutive_id(chr)) %>%
  filter(pos %in% genome_window$position) %>%
  mutate(relative_depth = rolling_mean/genome_median) %>%
  mutate(copy_number= relative_depth * ploidy) %>%
  mutate(chr_sums=chr_plot[index]) %>%  # for proper x-axis plotting
  mutate(plot_pos=ifelse(index==1, pos, (pos+chr_sums+(inter_chr_spacing*(index-1))))) # for proper x-axis plotting
## ---------------------------
## SNP freq calcs (pulls in position data from genome_depth dataframe)
genome_snp <- read.table(snp_file, header = TRUE)
genome_snp <- genome_snp %>%
  filter(chr != mito) %>%
  mutate(reads=rowSums(pick(A,T,G,C))) %>%
  mutate(across(c(A,T,C,G), ~ .x /reads, .names = "{.col}_freq"))

genome_snp <- genome_snp %>%
  mutate(snp_bin=(pos %/% window) * window +1) %>%
  left_join(genome_depth, by=c("chr","snp_bin"="pos"))

gaf <- genome_snp %>%  # sets a limit for allele frequency for heterozygosity and sums those within limit, per bin
  group_by(chr, snp_bin) %>%
  summarize(snp_count = sum((A_freq >= (1/copy_number)*0.5 & A_freq <=(1-(1/copy_number)*0.5)) |
                              (T_freq >= (1/copy_number)*0.5 & T_freq <=(1-(1/copy_number)*0.5)) |
                              (G_freq >= (1/copy_number)*0.5 & G_freq <=(1-(1/copy_number)*0.5)) |
                              (C_freq >= (1/copy_number)*0.5 & C_freq <=(1-(1/copy_number)*0.5))
)
)
## ---------------------------
## final dataframe of joined copy number, snps, and plotting positions per window
genome_depth <- genome_depth %>%
  left_join(gaf, by=c("chr", "pos"="snp_bin"))
## ---------------------------
## chromosome outlines and tick locations to add to final plot (tick marks are for Chr ID)
chroms <- genome_depth %>%
group_by(index) %>%
summarise(xmin=min(plot_pos), xmax=max(plot_pos), ymin=0, ymax=Inf)

ticks <- tapply(genome_depth$plot_pos, genome_depth$index, quantile, probs = 0.5)
## ---------------------------
## plot linear genome
p <- ggplot(genome_depth) +
  scale_color_gradient(low=snp_low,high=snp_high, na.value = "white", guide = "none") +
  geom_segment(aes(x = plot_pos, y = 0, color = snp_count, xend = plot_pos, yend = Inf)) +
  geom_segment(aes(x = plot_pos, y = ifelse(copy_number <= ploidy*ploidy_multiplier, copy_number, Inf),
                   xend = plot_pos, yend = ploidy), alpha = 0.9, color = copy_number) +
  geom_rect(data=chroms, aes(group=index, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            linewidth = chrom_line_width, fill = NA, colour = chrom_outline_color, linejoin = "round", inherit.aes = FALSE) +
  ylab(sample_id) +
  scale_x_continuous(expand = c(0, 0), breaks = ticks, labels = chr_ids) +
  scale_y_continuous(name = NULL, limits = c(0, ploidy*ploidy_multiplier), breaks = y_axis_labels) +
  theme_classic() +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.ticks = element_line(color = NA),
        axis.line = element_blank(),
        axis.text = element_text(size = 12))
## ---------------------------
## save plot as jpg. Height to width ratio is eyeballed for now
ggsave(sprintf("%s%s_%s_%s_%sbp.jpg", save_dir, Sys.Date(), sample_id, ref, window),
       p, width = 18, height = 1.7, units = "in")
## ---------------------------
## Save dataframes as excel
outfiles <- list(plotting_data=genome_depth,
                 read_depth_summary=chr_median,
                 raw_genome_median=as.data.frame(raw_genome_median),
                 corrected_genome_median=as.data.frame(genome_median))

write_xlsx(outfiles, path = paste0(save_dir,Sys.Date(),"_",sample_id,".xlsx"))
