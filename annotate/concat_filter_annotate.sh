#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --time=1:00:00
#SBATCH -p msilarge,msismall
#SBATCH -o %x_%u_%j.out
#SBATCH -e %x_%u_%j.err

# Concatenate chromosome vcfs into genome-wide vcf from Freebayes variant calling
# Bcftools Filter: remove complex variants, remove fixed variants, mapping quality >40,
#  strand balance probability of alt > 0, number of reads on reverse strand > 0,
#  number of reads right/left of alt > 1
# Snpeff annotate using manually built DB
# SnpSift for high and moderate impact variants, zip and index

set -ue
set -o pipefail

chr_dir=chr_vcf # location of vcf files to combine
raw_vcf=  # name for unfiltered vcf, include .vcf
sorted_vcf= # name for sorted, contatenated vcf, include .vcf
allele_num=20 # number of samples in vcf multiplied ploidy of organism
bcftools_out=  # intermediate vcf file, include .vcf
snpeff=/home/selmecki/shared/software/snpEff/snpEff.jar
snpeff_config=/home/selmecki/shared/software/snpEff/snpEff.config
snpeff_db=  # must be in snpeff.config file and must be name of directory in snpeff data subdir
annotate_vcf= # another intermediate vcf file, include .vcf
snpsift=/home/selmecki/shared/software/snpEff/SnpSift.jar
final_vcf=  # fully annotated and filtered output, include vcf
genotype_table=  # tab delimited table for use with R scripts (MCA using FactoMineR)

#Load modules
module load bcftools/1.10.2
module load htslib/1.9

# concatenate regions, leave uncompressed for faster piping below
find "$chr_dir" -name "*.vcf" > chr_files.txt
bcftools concat -f chr_files.txt -o "${raw_vcf}"

# sort samples alphanumerically (original freebayes sorting is random)
bcftools query -l "${raw_vcf}" | sort > samples.txt
bcftools view -S samples.txt "${raw_vcf}" > "${sorted_vcf}"

# view -e option excludes (here is complex variants and fixed SNPs) 
# view -i option includes (mapping quality at least 40, alt reads forward and reverse, reads on both sides of alt allele)
bcftools view -e "INFO/TYPE='complex'" "${sorted_vcf}" \
| bcftools view -e "INFO/AC=${allele_num}" \
| bcftools view -i \
"INFO/MQM>=40 && INFO/SAR>=1 && INFO/SAP>0 && INFO/RPL>1 && INFO/RPR>1" \
 -o "${bcftools_out}"

# annotate using snpeff with manually built database
java -Xmx4g -jar "${snpeff}" -c "${snpeff_config}" "${snpeff_db}" \
"${bcftools_out}" >  "${annotate_vcf}"

# OPTIONAL: filter with snpsift by predicted impact - this removes intergenic and modifiers
java -Xmx4g -jar "${snpsift}" filter "ANN[*].IMPACT has 'HIGH' \
| ANN[*].IMPACT has 'MODERATE' | ANN[*].IMPACT has 'LOW' " "${annotate_vcf}" > "${final_vcf}"

# subset annotated (but not snpsift filtered) to just SNPs
# and output tab-delimited file for use in R MCA script for preliminary clustering
bcftools view -e 'GT="mis"' "${annotate_vcf}" \
| bcftools view -e 'TYPE!="indel"' \
| bcftools query -H -f '%CHROM\t%POS[\t%GT]\n' > "${genotype_table}"

# fix header of genotype table
sed -i -e '1s;\[[0-9]\{1,3\}\];;g' \
    -e '1s/^# //' \
    -e '1s/:GT//g' "${genotype_table}"

# zip and index output (facilitate IGV loading)
# save annotated vcf for additional clustering scripts
bgzip "${annotate_vcf}"
bcftools index -t "${annotate_vcf}".gz
bgzip "${final_vcf}"
bcftools index -t "${final_vcf}".gz
bcftools stats "${final_vcf}.gz" > ../logs/"${final_vcf}.stats"

# remove intermediate files
rm "${raw_vcf}"
rm "${sorted_vcf}"
rm "${bcftools_out}"