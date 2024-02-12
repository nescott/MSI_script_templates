#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --time=48:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %x_%u_%A_%a.out
#SBATCH -e %x_%u_%A_%a.err
#SBATCH --array=

# call variants for all samples in a population using freebayes, subsetting by region (chromosome)
# array size == number of chromosomes/contigs in reference file
# time needed depends on number of samples (C. albicans 100 samples required >36 hours)
# input files: a list of bam files (with path if needed) and a list of regions
# output: multi-sample vcf file per each region, file name is the region name (may want to change this) 
# NOTE: freebayes options are -C (number of supporting reads), -F (min alt allele frequency for calling)
#  and -p (ploidy) - change as needed, should probably change these to vars at top of script
set -ue
set -o pipefail

line=${SLURM_ARRAY_TASK_ID}
region_list=
ref_fasta=   # FREEBAYES REQUIRES UNZIPPED REF so annoying
bam_list=bam.files
species=   # no spaces
ref=  # abbreviation for reference genome

mkdir -p chr_vcf

chr=$(awk -v val="$line" 'NR == val { print $0}' $region_list)

#Load modules
module load samtools/1.10
module load freebayes/20180409

freebayes -f "${ref_fasta}" \
  -C 10 \
  -F 0.4 \
  -p 2 \
  -r "${chr}" \
  -L "${bam_list}" \
  -v "chr_vcf/${species}_${ref}_${chr}.vcf"
