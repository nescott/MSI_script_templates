#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=3gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=3:00:00
#SBATCH -p amdsmall,amdlarge,amd512,amd2tb
#SBATCH -o %x_%u_%j.out
#SBATCH -e %x_%u_%j.err

#Perform variant calling on a sorted, markduped bam using Mutect2

set -ue
set -o pipefail

ref_fasta=
sample_bam=
normal_bam=
normal_name= #use sample name from bam header
strain=

#Load GATK module 
module load gatk/4.1.2

# Run and filter calls
if [ ! -d "vcf" ]; then
  mkdir vcf
fi

gatk Mutect2 -R "${ref_fasta}" -I "${sample_bam}" -I "${normal_bam}" -normal "${normal_name}" -O "${strain}"_unfiltered.vcf
gatk FilterMutectCalls -R "${ref_fasta}" -V "${strain}"_unfiltered.vcf -O "${strain}"_filtered.vcf

