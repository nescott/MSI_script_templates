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
##  minimap2 -x asm10 genome1.fa genome2.fa > align.paf
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
