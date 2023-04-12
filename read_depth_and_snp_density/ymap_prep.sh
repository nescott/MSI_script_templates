#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=3gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH -t 30
#SBATCH -p msismall,msilarge
#SBATCH -o %x_%u_%j.out
#SBATCH -e %x_%u_%j.err
#SBATCH --array=3-8
set -ue
set -o pipefail

module load samtools/1.10

line=${SLURM_ARRAY_TASK_ID}
ref=sc5314  # short ID
snp_bam_file=bam.files
depth_bam_file=gc_corrected.files
fasta=/home/selmecki/shared/disaster_recovery/Reference_Genomes/SC5314_A21/C_albicans_SC5314_version_A21-s02-m09-r08_chromosomes.fasta

# depth and allele counts per bam file from samtools

snp_bam=$(awk -v val="$line" 'NR == val {print $0}' $snp_bam_file)
depth_bam=$(awk -v val="$line" 'NR == val {print $0}' $depth_bam_file)
snp_strain=$(basename "$snp_bam" | cut -d "_" -f 1)
depth_strain=$(basename "$depth_bam" | cut -d "_" -f 1)

samtools depth -aa -o "${depth_strain}"_"${ref}"_gc_corrected_depth.txt "${depth_bam}"

sed -i "1s/^/chr\tpos\tdepth\n/" "${depth_strain}"_"${ref}"_gc_corrected_depth.txt

samtools mpileup -f "${fasta}" "${snp_bam}" | awk '{print $1, $2, $3, $4, $5}' > "${snp_strain}".pileup

python3 berman_count_snps_v5.py "${snp_strain}".pileup > "${snp_strain}"_"${ref}"_putative_SNPs.txt

sed -i "1s/^/chr\tpos\tref\tA\tT\tC\tG\n/" "${snp_strain}"_"${ref}"_putative_SNPs.txt

rm "${snp_strain}".pileup
