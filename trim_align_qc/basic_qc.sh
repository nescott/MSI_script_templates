#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=8:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %x_%u_%j.out
#SBATCH -e %x_%u_%j.err

set -ue
set -o pipefail
unset DISPLAY

fastq_dir=
bam_dir=
logs_dir=
temp_dir=/scratch.global/scot0854  #scratch directory
out_dir=
file_name="$(date +%F)_"

# Use local modules
module use /home/selmecki/scot0854/modulefiles.local

# Load modules
module load fastqc/0.11.9
module load qualimap/20221111
module load multiqc/20221111

# raw fastq file qc
find "$fastq_dir" -name "*trimmed_*P.fq" -exec fastqc -o -t 8 "$temp_dir" {} \;

# bam qc
find "$bam_dir" -name "*.bam" \
-exec qualimap bamqc -bam {} -outdir $temp_dir/{} --java-mem-size=4G  \;

# multiqc

multiqc "$temp_dir" "$logs_dir" -o "$out_dir" -n "$file_name"
