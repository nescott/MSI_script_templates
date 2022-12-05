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

# Sample and reference variables --------------------------
sample_file=
line=${SLURM_ARRAY_TASK_ID}
ploidy_number=1
reference_fasta=/home/selmecki/shared/Reference_Genomes/C_lusitaniae/ATCC42720/GCA_000003835.1_ASM383v1_genomic.fna
bam_dir=

# make directory for gvcf files
# Get strain ID froms sample file line equal to array task ID
strain=$(awk -v val="${line}" 'NR == val { print $1}' "${sample_file}")

# Get bam file from strain and path
bam="${bam_dir}""${strain}"_trimmed_bwa_sorted_markdup.bam

# Load modules
module load gatk

# HaplotypeCaller and Genotype

srun gatk --java-options "-Xmx6g -XX:ParallelGCThreads=2" HaplotypeCaller \
 -R "${reference_fasta}" \
 -I "${bam}" \
 -O "${strain}"_haplotypecaller.g.vcf \
 -ploidy $((ploidy_number)) \
 -ERC GVCF 

gatk --java-options "-Xmx6g" GenotypeGVCFs \
 -R "${reference_fasta}" \
 -V "${strain}"_haplotypecaller.g.vcf \
 -ploidy $((ploidy_number)) \
 -O "${strain}"_haplotypecaller_genotype.g.vcf

# Start of Hard-filtering to SNP-only and indels-only sets
# Subset to SNPs-only callset with SelectVariants
gatk --java-options "-Xmx6g"  SelectVariants \
    -V "${strain}"_haplotypecaller_genotype.g.vcf \
    -select-type SNP \
    -O "${strain}"_snps.vcf

# Subset to indels-only callset with SelectVariants
gatk --java-options "-Xmx6g" SelectVariants \
    -V "${strain}"_haplotypecaller_genotype.g.vcf \
    -select-type INDEL \
    -O "${strain}"_indels.vcf

# Filter SNPS (recommend experimentation with filter parameters)
gatk --java-options "-Xmx6g" VariantFiltration \
-V "${strain}"_snps.vcf \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O "${strain}"_snps_filtered.vcf

# Filter indels (recommend experimentation with filter parameters) - note different parameters for SNPS vs indels
gatk VariantFiltration \
-V "${strain}"_indels.vcf \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "FS > 200.0" --filter-name "FS200" \
-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
-O "${strain}"_indels_filtered.vcf

# Merge SNPs and indels into one vcf
java -jar /panfs/roc/msisoft/picard/2.18.16/picard.jar MergeVcfs \
I="${strain}"_snps_filtered.vcf \
I="${strain}"_indels_filtered.vcf \
O="${strain}"_snps_indels_filtered.vcf
