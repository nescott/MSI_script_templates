#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --time=1:30:00
#SBATCH -p msismall,msilarge
#SBATCH -o %x_%u_%j.out
#SBATCH -e %x_%u_%j.err
#SBATCH --array=1
set -ue
set -o pipefail

module load samtools/1.10
module load R/4.1.0

line=${SLURM_ARRAY_TASK_ID}
ref=sc5314  # short ID
snp_bam_file=bam_files.txt # text list of bam files with path
depth_bam_file=gc_corrected_files.txt # text list of GC-corrected bams with path
fasta=/home/selmecki/shared/disaster_recovery/Reference_Genomes/SC5314_A21/C_albicans_SC5314_version_A21-s02-m09-r08_chromosomes.fasta

# get rid of big intermediate file including if slurm script fails before finishing
function finish {
  if [ -e alleles/"${snp_strain}".pileup ]; then
    rm alleles/"${snp_strain}".pileup
  fi
}
trap finish EXIT

mkdir -p alleles depth plots

# depth and allele counts per bam file from samtools

# process input files, using array ID for line number:
# get the SNP bam file, the depth bam file, and the strain IDs
snp_bam=$(awk -v val="$line" 'NR == val {print $0}' $snp_bam_file)
depth_bam=$(awk -v val="$line" 'NR == val {print $0}' $depth_bam_file)
snp_strain=$(basename "$snp_bam" | cut -d "_" -f 1)
if [ "$snp_strain" == "AMS" ]; then
    snp_strain=$(basename "$snp_bam" | cut -d "_" -f 1,2)
fi
depth_strain=$(basename "$depth_bam" | cut -d "_" -f 1)
if [ "$depth_strain" == "AMS" ]; then
    depth_strain=$(basename "$depth_bam" | cut -d "_" -f 1,2)
fi

# if the snp strain and depth strain don't match (meaning the input files were
# sorted differently), then exit with error
[[ "$depth_strain" != "$snp_strain" ]] && { >&2 echo "Depth and snp IDs don't match"; exit 1; }

# read depth per position
samtools depth -aa -o depth/"${depth_strain}"_"${ref}"_gc_corrected_depth.txt "${depth_bam}"

# add a standard header to make the R script happy
sed -i "1s/^/chr\tpos\tdepth\n/" depth/"${depth_strain}"_"${ref}"_gc_corrected_depth.txt

# alleles per position (filters on mapping quality of 60, assuming BWA MEM
# alignment), grab the columns of interest
samtools mpileup -q 60 -f "${fasta}" "${snp_bam}" | awk '{print $1, $2, $3, $4, $5}' > alleles/"${snp_strain}".pileup

# use a borrowed python script to parse mpileup into neat columns of A, T, G, C
python3 /home/selmecki/shared/software/berman_count_snps_v5.py alleles/"${snp_strain}".pileup > alleles/"${snp_strain}"_"${ref}"_putative_SNPs.txt

# add standard header
sed -i "1s/^/chr\tpos\tref\tA\tT\tG\tC\n/" alleles/"${snp_strain}"_"${ref}"_putative_SNPs.txt

# OPTIONAL - run R script and output plots automatically
# Make sure you've set all the variables in the  R script first!
Rscript --vanilla ../genome_vis.R depth/"${depth_strain}"_"${ref}"_gc_corrected_depth.txt alleles/"${snp_strain}"_"${ref}"_putative_SNPs.txt "${depth_strain}"
