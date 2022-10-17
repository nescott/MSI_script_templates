## ---------------------------
## Script name: PAFR_genome_comparison_plot.R
##
## Purpose of script: Plot genome-genome alignments
##
## Author: Nancy Scott
##
## Date Created: 2022-04-19
##
## Email: scot0854@umn.edu
## ---------------------------
## Notes: minimap2 for aligning genome fasta files (default ouput is paf)
## Consider trying minimap options for max sequence divergence (-x asm5/asm10/asm20)
## Original genome scripts and plots are gone; original gaps in ref genome plots 
## not reproducible so far
## and using pafr to filter by primary alignment or mapping quality
## ---------------------------
## load packages
library(pafr)
library(ggplot2)

# load alignment files
#NAME_ME <- read_paf("FILE_PATH")

# put all plotting details in one place in a function for dot plots
# can probably mess around with ggplot for prettier naming on output
genome_dotplot <- function(paf){
    dotplot(paf, label_seqs = T, order_by = "qstart") + theme_bw()
}

# same for coverage
genome_coverage_plot <- function(paf){
    plot_coverage(paf, fill = 'qname') + scale_fill_brewer(palette = "Set1")
}
