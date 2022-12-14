#!/bin/bash

set -ue
set -o pipefail

# how to get effective genome size:
# ~/bin/faCount $ref.fna.gz -summary
# how to get 2bit file format for deeptools commands:
#  ~/bin/faToTwoBit GCA_003675555.2_ASM367555v2_genomic.fna.gz GCA_00367555.2bit

window=2500  #bp window for tiling across genome
ref=fda  #short ID
in_dir=../../align/fda/bam/
bam=$(find $in_dir -mindepth 1 -maxdepth 1 -type f -name "*.bam")
genome_size=12152819
ref2bit=~/umn/C_lus_refs/ASM1463611v1/GCA014636115.2bit

for b in $bam; do
  strain=$(basename $b | cut -d "_" -f 1)  

  #~/bin/mosdepth -T 1,10,50,100,200 -n --by "${window}" -t 2 \
  #"${strain}_${ref}_deeptools_${window}" \
  #gc_corrected_bams/"${strain}_${ref}_deeptools.bam"
  
  #mean=$(zcat "${strain}_${ref}_deeptools_${window}".regions.bed.gz \
  #| awk '{total += $4} END { print total/NR}')

  #zcat "${strain}_${ref}_deeptools_${window}".regions.bed.gz \
  #| awk -v val="$mean" 'BEGIN{OFS="\t"} {print $1,$2,$3,$4/val}' \
  #> "${strain}_${ref}_deeptools_${window}"_mosdepth.gg.tab
  
  #sed -f "sed_${ref}_chr" < "${strain}_${ref}_deeptools_${window}"_mosdepth.gg.tab \
  #> "${strain}_${ref}_deeptools_${window}"_chrnum_mosdepth.gg.tab

  sed i "1 s/^/chr\tstart\tstop\t${strain}\n/" \
  "${strain}_"${ref}_deeptools_${window}"_chrnum_mosdepth.gg.tab
done
