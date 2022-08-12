library(ggplot2)
library(RColorBrewer)

# Input file: gc-corrected bam (from deeptools) <- mosdepth calculates mean per-window depth <- 
# generate tab-delimited file with updated chromosomes numbers (using sed) for ggplot-friendly file

strain <- "AMS5200"
tile <- "1500bp"
ref <- "UCRiverside"

genome_windows <- read.table('data/AMS5200_deeptools_ucr_1500_chrnum_mosdepth.gg.tab')
colnames(genome_windows) <-  c("Chr","Start","End","Depth")
chrlist <- c(1:8)
genome_depth <- genome_windows[genome_windows$Chr %in% chrlist, ]
genome_depth <- genome_depth[order(genome_depth$Chr,genome_depth$Start), ]
genome_depth$index <- rep.int(seq_along(unique(genome_depth$Chr)), times=tapply(genome_depth$Start,genome_depth$Chr,length))
genome_depth$Chr <- as.character(genome_depth$Chr)

genome_depth$pos = NA

ncrh <- length(unique(genome_depth$Chr))
lastbase=0
ticks = NULL
minor = vector(,8)
for (i in 1:8 ) {
  if (i==1) {
    genome_depth[genome_depth$index==i, ]$pos=genome_depth[genome_depth$index==i, ]$Start
  } else {
    ## chromosome position maybe not start at 1, eg. 9999. So gaps may be produced. 
    lastbase = lastbase + max(genome_depth[genome_depth$index==(i-1),"Start"])
    minor[i] = lastbase
    genome_depth[genome_depth$index == i,"Start"] =
      genome_depth[genome_depth$index == i,"Start"]-min(genome_depth[genome_depth$index==i,"Start"]) +1
    genome_depth[genome_depth$index == i,"End"] = lastbase
    genome_depth[genome_depth$index == i, "pos"] = genome_depth[genome_depth$index == i,"Start"] + lastbase
  }
}
ticks <-tapply(genome_depth$pos,genome_depth$index,quantile,probs=0.5)
ticks
minorB <- tapply(genome_depth$End,genome_depth$index,max,probs=0.5)
minorB
minor
xmax = ceiling(max(genome_depth$pos) * 1.03)
xmin = floor(max(genome_depth$pos) * -0.03)
#Title="Depth of Sequencing Coverage"

p <- ggplot(genome_depth, aes(x=pos,y=Depth, xend=pos, yend=1.0, color=Chr))  + 
  scale_colour_brewer(palette = "Set2") +
  geom_segment(alpha=0.4, size=0.5) +
  geom_point(alpha=0.9,size=0.5,shape=20) +
  labs(title=sprintf("Depth of sequencing coverage, %s, %s, %s",strain, ref, tile),xlab="Position",y="Normalized Read Depth") +
  scale_x_continuous(name="Chromosome", expand = c(0, 0),
                     breaks=ticks,
                     labels=(unique(genome_depth$CHR))) +
  scale_y_continuous(name="Normalized Read Depth", expand = c(0, 0),
                     limits = c(0,3)) + theme_classic() +
  guides(fill = guide_legend(keywidth = 3, keyheight = 1))
ggsave(sprintf("~/plots/%s%s%s.pdf", strain, tile, ref),p,width=7,height=2.5)

manualColors = c("dodgerblue2","red1","grey20")
for (n in chrlist ) {
  Title=sprintf("Chr%s, %s, %s reference, %s", n, strain, ref, tile)
  print(Title)
  l <- subset(genome_depth,genome_depth$Chr==n)
  l$bp <- l$Start
  p<-ggplot(l,
            aes(x=bp,y=Depth, xend=bp, yend=1.0)) +
    geom_segment(alpha=0.6,size=0.3) +
    geom_point(alpha=0.6,size=0.3,shape=20) +
    scale_color_manual(values = manualColors) +
    labs(title=Title,xlab="Position",ylab="Normalized Read Depth") +
    scale_x_continuous(expand = c(0, 0), name="Position") +
    scale_y_continuous(name="Normalized Read Depth", expand = c(0, 0),
                       limits = c(0,3)) + theme_classic() +
    guides(fill = guide_legend(keywidth = 3, keyheight = 1))
  ggsave(sprintf("~/plots/%s_%s_%s_Chr%s.pdf",strain, ref, tile, n),p,width=7,height=2.5)
  p
}

