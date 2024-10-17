#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=1:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --array=1
set -ue
set -o pipefail

# local modules

module use "$HOME"/modulefiles.local
module load mosdepth/20221013 # Calculates depth per tiled window

line=${SLURM_ARRAY_TASK_ID}
window=500  # bigger bp window smooths visualization
ref=sc5314  # short ID
bam_file=calbicans.bam # text file listing bams with paths
mtdna= # contig/scaffold ID for removal so mean correction is only of nuclear genome

# output dirs
mkdir -p "mosdepth_txt" "mosdepth_bed" "tab"

# calculate uncorrected nuclear mean and perform mean-correction,
# generate tab-delimited table for graphing,
# ready for use with candida karyoploter R scripts or custom ggplot2 scripts

in_bam=$(awk -v val="$line" 'NR == val {print $0}' $bam_file)
strain=$(basename "$in_bam" | cut -d "_" -f 1,2)

mosdepth -T 1,10,50,100,200 -n --by "${window}" -t 4 \
"${strain}_${ref}_deeptools_${window}" gc_corrected_bams/"${strain}_${ref}_deeptools.bam"

zcat "${strain}_${ref}_deeptools_${window}".regions.bed.gz \
| grep -v "${mtdna}" > "${strain}_${ref}_deeptools_${window}".nuclear.bed

gzip "${strain}_${ref}_deeptools_${window}".nuclear.bed

mean=$(zcat "${strain}_${ref}_deeptools_${window}".nuclear.bed.gz \
| awk '{total += $4} END { print total/NR}')

zcat "${strain}_${ref}_deeptools_${window}".nuclear.bed.gz \
| awk -v val="$mean" 'BEGIN{OFS="\t"} {print $1,$2,$3,$4/val}' \
> "${strain}_${ref}_deeptools_${window}"_mosdepth.gg.tab

sed -i "1s/^/chr\tstart\tstop\t${strain}\n/" \
"${strain}_${ref}_deeptools_${window}"_mosdepth.gg.tab

mv ./"${strain}"_*.txt mosdepth_txt
mv ./"${strain}"_*.tab tab
mv ./"${strain}"_*.bed* mosdepth_bed
