#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=12:00:00
#SBATCH -p amdsmall,amdlarge,amd512,amd2tb,small,large,max
#SBATCH -o %x_%u_%A_%a.out
#SBATCH -e %x_%u_%A_%a.err
#SBATCH --array=1-8

# call variants for all samples in a population using freebayes, subsetting by region
set -ue
set -o pipefail

species=Clusitaniae 
ref=fda
ref_fasta=/home/selmecki/shared/Reference_Genomes/C_lusitaniae/ASM1463611v1/GCA_014636115.1_ASM1463611v1_genomic.fna
line=${SLURM_ARRAY_TASK_ID} 
bam_list=bam.files
region_list=regions.txt

chr=$(awk -v val="$line" 'NR == val { print $0}' $region_list)

#Load modules 
module load samtools/1.10
module load freebayes/20180409

freebayes -f "${ref_fasta}" \
  -C 10 \
  -F 0.4 \
  -p 1 \
  -r "${chr}" \
  -L "${bam_list}" \
  -v "${species}_${ref}_${chr:0:15}.vcf"
