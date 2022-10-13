#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=15gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=3:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %x_%u_%j.out
#SBATCH -e %x_%u_%j.err

set -ue
set -o pipefail

# local modules

module use "$HOME"/modulefiles.local
module load deeptools/20221013 # computes and corrects GC bias
# Requires effective genome size and reference genome in 2bit format
module load mosdepth/20221013 # Calculates depth per tiled window

window=2500  #bigger bp window smooths visualization
ref=fda  #short ID
in_dir=../bam
bam=$(find $in_dir -mindepth 1 -maxdepth 1 -type f -name "AMS5200*.bam")
genome_size=12152819
ref2bit="$HOME"/faToTwoBit/GCA01436115.1.2bit

# output dirs
arr=("gc_corrected_bams" "mosdepth_txt" "mosdepth_bed" "tab")
for d in "${arr[@]}"; do
  if [ ! -d "$d" ]; then
    mkdir "$d"
  fi
done

# for each bam, compute and correct GC bias, calculate depth,
# calculate mean and perform mean-correction,
# generate tab-delimited table for graphing,
# switch to chr number instead of accession scaffold id,
# generate headers:
# ready for use with candida karyoploter R scripts

for b in $bam; do
  strain=$(basename "$b" | cut -d "_" -f 1)

  computeGCBias -b "$b" --effectiveGenomeSize "${genome_size}" -g "${ref2bit}" \
  -o gc_corrected_bams/"${strain}_${ref}_freq.txt"

  correctGCBias -b "$b" --effectiveGenomeSize "${genome_size}" -g "${ref2bit}" \
  -freq gc_corrected_bams/"${strain}_${ref}_freq.txt" \
  -o gc_corrected_bams/"${strain}_${ref}_deeptools.bam"

  mosdepth -T 1,10,50,100,200 -n --by "${window}" -t 4 \
  "${strain}_${ref}_deeptools_${window}" gc_corrected_bams/"${strain}_${ref}_deeptools.bam"

  mean=$(zcat "${strain}_${ref}_deeptools_${window}".regions.bed.gz \
  | awk '{total += $4} END { print total/NR}')

  zcat "${strain}_${ref}_deeptools_${window}".regions.bed.gz \
  | awk -v val="$mean" 'BEGIN{OFS="\t"} {print $1,$2,$3,$4/val}' \
  > "${strain}_${ref}_deeptools_${window}"_mosdepth.gg.tab

  sed -f "sed_${ref}_chr" < "${strain}_${ref}_deeptools_${window}"_mosdepth.gg.tab \
  > "${strain}_${ref}_deeptools_${window}"_chrnum_mosdepth.gg.tab

  sed -i "1s/^/chr\tstart\tstop\t${strain}\n/" \
  "${strain}_${ref}_deeptools_${window}"_chrnum_mosdepth.gg.tab

  mv AMS*.txt mosdepth_text/
  mv AMS*.tab tab/
  mv AMS*.bed* mosdepth_bed/
done
