#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=4:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %x_%u_%j.out
#SBATCH -e %x_%u_%j.err

set -ue
set -o pipefail

fastq_dir=
bam_dir=
out_dir=
file_name=

# Use local modules
module use /home/selmecki/scot0854/modulefiles.local

# Load modules
module load fastqc/0.11.9
module load qualimap/20221111
module load multiqc/20221111

# raw fastq file qc
find "$fastq_dir" -name "*trimmed_*P.fq" -exec fastqc -o -t 8 "$out_dir" {} \;

# bam qc
find "$bam_dir" -name "*.bam" \
-exec qualimap bamqc -bam -outdir "$out_dir" --java-mem-size=3G {} \;

# multiqc

multiqc logs/ "$outdir" -o "$out_dir" -n "$file_name"

