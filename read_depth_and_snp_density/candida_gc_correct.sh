#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --time=30
#SBATCH -p msismall,msilarge
#SBATCH -o %x_%u_%j.out
#SBATCH -e %x_%u_%j.err
#SBATCH --array=1
set -ue
set -o pipefail

# local modules

module use /home/selmecki/shared/software/modulefiles.local
module load deeptools/20221013 # computes and corrects GC bias
# Requires effective genome size (14320608 for C. albicans SC5314 A21) and reference genome in 2bit format
# which is in the Reference_Genomes/SC5314_A21 directory

line=${SLURM_ARRAY_TASK_ID}
ref=sc5314  # short ID
bam_file=bam_files.txt  # a text file of bam files with path if necessary, to iterate over as array
gc_dir=/scratch.global/YOUR_X500/ # with trailing slash, for output of GC-corrected BAMs, consider using /scratch.global/UMN_ID/
genome_size=14320608 # C. albicans effective genome size
ref2bit=/home/selmecki/shared/disaster_recovery/Reference_Genomes/SC5314_A21/C_albicans_SC5314_version_A21-s02-m09-r08_chromosomes.2bit

# for each bam, compute and correct GC bias, calculate depth,
# generate tab-delimited table for plotting and downstream analysis

in_bam=$(awk -v val="$line" 'NR == val {print $0}' $bam_file)
strain=$(basename "$in_bam" | cut -d "_" -f 1)
if [ "$strain" == "AMS" ]; then
    strain=$(basename "$in_bam" | cut -d "_" -f 1,2)
fi

computeGCBias -b "$in_bam" --effectiveGenomeSize "${genome_size}" -g "${ref2bit}" \
-o "${gc_dir}${strain}_${ref}_freq.txt"

correctGCBias -b "$in_bam" --effectiveGenomeSize "${genome_size}" -g "${ref2bit}" \
-freq "${gc_dir}${strain}_${ref}_freq.txt" \
-o "${gc_dir}${strain}_${ref}_deeptools.bam"