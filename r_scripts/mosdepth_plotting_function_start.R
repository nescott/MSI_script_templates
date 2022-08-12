library(ggplot2)
library(RColorBrewer)

fda_windows <- read.table("")
colnames(fda_windows) <-  c("Chr","Start","End","Depth")
chrlist <- c(1:8)
fda <- fda_windows[fda_windows$Chr %in% chrlist, ]
fda <- d[order(fda$Chr,fda$Start),]
fda$index <- rep.int(seq_along(unique(fda$Chr)), times=tapply(fda$Start,fda$Chr,length))
fda$pos = NA

read_depth <-  function(mosdepth_formatted) {
  windows <- read.table(mosdepth_formatted)
  colnames(windows) <- c("Chr","Start","End","Depth")
  chrlist <- c(1:8)
  chr_windows <- windows[windows$Chr %in% chrlist]
  chr_windows <- chr_windows[order(chr_windows$Chr, chr_windows$Start)]
  chr_windows$index <- rep.int(seq_along(unique(chr_windows$Chr)), times = tapply(chr_windows$Start, chr_windows$Chr, length))
  chr_windows$pos = NA
}