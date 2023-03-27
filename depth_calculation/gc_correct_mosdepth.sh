#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=2:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %x_%u_%j.out
#SBATCH -e %x_%u_%j.err
#SBATCH --array=1-
set -ue
set -o pipefail

# local modules

module use "$HOME"/modulefiles.local
module load deeptools/20221013 # computes and corrects GC bias
# Requires effective genome size and reference genome in 2bit format,
# generated from faToTwoBit and faCount in ~/bin, see README in same dir
module load mosdepth/20221013 # Calculates depth per tiled window

line=${SLURM_ARRAY_TASK_ID}
window=500  # bigger bp window smooths visualization
ref=sc5314  # short ID
bam_file= # to iterate over as array
genome_size=14320608  # see above to get this number
ref2bit="$HOME"/bin/faToTwoBit/C_albicans_SC5314_version_A21-s02-m09-r08_chromosomes.2bit

# output dirs
arr=("gc_corrected_bams" "mosdepth_txt" "mosdepth_bed" "tab")
for d in "${arr[@]}"; do
  if [ ! -d "$d" ]; then
    mkdir "$d"
  fi
done

# for each bam, compute and correct GC bias, calculate depth,
# generate tab-delimited table for plotting and downstream analysis

in_bam=$(awk -v val="$line" 'NR == val {print $0}' $bam_file)
strain=$(basename "$in_bam" | cut -d "_" -f 1,2)

computeGCBias -b "$in_bam" --effectiveGenomeSize "${genome_size}" -g "${ref2bit}" \
-o gc_corrected_bams/"${strain}_${ref}_freq.txt"

correctGCBias -b "$in_bam" --effectiveGenomeSize "${genome_size}" -g "${ref2bit}" \
-freq gc_corrected_bams/"${strain}_${ref}_freq.txt" \
-o gc_corrected_bams/"${strain}_${ref}_deeptools.bam"

mosdepth -T 1,10,50,100,200 -n --by "${window}" -t 4 \
"${strain}_${ref}_deeptools_${window}" gc_corrected_bams/"${strain}_${ref}_deeptools.bam"

zcat "${strain}_${ref}_deeptools_${window}".regions.bed.gz \
| awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' \
> "${strain}_${ref}_deeptools_${window}"_mosdepth.gg.tab

sed -i "1s/^/chr\tstart\tstop\t${strain}\n/" \
"${strain}_${ref}_deeptools_${window}"_mosdepth.gg.tab

mv ./*.txt mosdepth_txt
mv ./*.tab tab
mv ./*.bed* mosdepth_bed
