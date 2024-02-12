#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --mem=10gb
#SBATCH -o %x_%u.out
#SBATCH -e %x_%u.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=6:00:00
#SBATCH -p amdsmall,amdlarge,amd512,amd2tb,small,large,max

set -ue
set -o pipefail

# Sample and reference variables --------------------------
ploidy_number=
reference_fasta=
species=
ref=

mkdir -p "$PWD/genotyped_vcf" "$PWD/filtered_vcf"

# Load modules
module load gatk
module load htslib
module load bcftools

# HaplotypeCaller and Genotype
gatk --java-options "-Xmx6g" GenotypeGVCFs \
 -R "${reference_fasta}" \
 -V gendb://db/"${species}_${ref}" \
 -ploidy $((ploidy_number)) \
 -O genotyped_vcf/"${species}"_"${ref}"_haplotypecaller_genotype.g.vcf

# Start of Hard-filtering to SNP-only and indels-only sets
# Subset to SNPs-only callset with SelectVariants
gatk --java-options "-Xmx6g"  SelectVariants \
    -V genotyped_vcf/"${species}"_"${ref}"_haplotypecaller_genotype.g.vcf \
    -select-type SNP \
    -O genotyped_vcf/"${species}"_"${ref}"_snps.vcf

# Subset to indels-only callset with SelectVariants
gatk --java-options "-Xmx6g" SelectVariants \
    -V genotyped_vcf/"${species}"_"${ref}"_haplotypecaller_genotype.g.vcf \
    -select-type INDEL \
    -O genotyped_vcf/"${species}"_"${ref}"_indels.vcf

# Filter SNPS (recommend experimentation with filter parameters)
gatk --java-options "-Xmx6g" VariantFiltration \
-V genotyped_vcf/"${species}"_"${ref}"_snps.vcf \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O filtered_vcf/"${species}"_"${ref}"_snps_filtered.vcf


# Filter indels (recommend experimentation with filter parameters) - note different parameters for SNPS vs indels
gatk VariantFiltration \
-V genotyped_vcf/"${species}"_"${ref}"_indels.vcf \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "FS > 200.0" --filter-name "FS200" \
-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
-O filtered_vcf/"${species}"_"${ref}"_indels_filtered.vcf

# Merge SNPs and indels into one vcf, compress and index
bgzip filtered_vcf/"${species}"_"${ref}"_snps_filtered.vcf
bgzip filtered_vcf/"${species}"_"${ref}"_indels_filtered.vcf
bcftools concat -a -o "${species}"_"${ref}"_merged_filtered.vcf
bgzip "${species}"_"${ref}"_merged_filtered.vcf
bcftools index -t "${species}"_"${ref}"_merged_filtered.vcf.gz
