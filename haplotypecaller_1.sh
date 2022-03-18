#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --mem=10gb
#SBATCH -o %x_%u_%A_%a.out
#SBATCH -e %x_%u_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=6:00:00
#SBATCH -p amdsmall,amdlarge,amd512,amd2tb,small,large,max
#SBATCH --array=1-

set -ue
set -o pipefail

# Sample and reference variables
sample_file=
line=${SLURM_ARRAY_TASK_ID}
ploidy_number=1
reference_fasta=/home/selmecki/shared/Reference_Genomes/C_lusitaniae/ATCC42720/GCA_000003835.1_ASM383v1_genomic.fna
bam_dir=

# make directory for gvcf files
# Get strain ID froms sample file line equal to array task ID
strain=$(awk -v val="${line}" 'NR == val { print $1}' "${sample_file}")

# Check for/create output dir
if [ ! -d "gvcf" ]; then
  mkdir gvcf
fi

# Get bam file from strain and path
bam="${bam_dir}""${strain}"_trimmed_bwa_sorted_markdup.bam

# Load modules
module load gatk

# HaplotypeCaller and Genotype
gatk --java-options "-Xmx6g -XX:ParallelGCThreads=2" HaplotypeCaller \
 -R "${reference_fasta}" \
 -I "${bam}" \
 -O gvcf/"${strain}"_haplotypecaller.g.vcf \
 -ploidy $((ploidy_number)) \
 -ERC GVCF 
