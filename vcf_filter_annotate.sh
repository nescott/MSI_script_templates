#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=2:00:00
#SBATCH -p msilarge,msismall
#SBATCH -o %x_%u_%j.out
#SBATCH -e %x_%u_%j.err

# Concatenate chromosome vcfs into genome-wide vcf from Freebayes variant calling
# Bcftools Filter: remove complex variants, remove fixed variants, mapping quality >40,
#  strand balance probability of alt > 0, number of reads on reverse strand > 0,
#  number of reads right/left of alt > 1, filter indels separated by =<10 bp
# Vcftools thin by SNPs =< 5 bp
# Snpeff annotate using manually built DB, zip and index

set -ue
set -o pipefail

chr_dir=chr_vcf
raw_vcf=  # include .vcf
sample_num=
bcftools_out=  # include .vcf
vcftools=/home/selmecki/shared/software/VCFtools/bin/vcftools
vcftools_out=  # do not include .vcf
snpeff=/home/selmecki/shared/software/snpEff/snpEff.jar
snpeff_config=/home/selmecki/shared/software/snpEff/snpEff.config
snpeff_db=
annotate_vcf=

#Load modules
module load bcftools/1.10.2
module load htslib/1.9

find "$chr_dir" -name "*.vcf" > chr_files.txt
bcftools concat -f chr_files.txt -o "${raw_vcf}"

bcftools view -e "INFO/TYPE='complex'" "${raw_vcf}" \
| bcftools view -e "INFO/AC=${sample_num}" \
| bcftools view -i \
"INFO/MQM>=40 && INFO/SAR>=1 && INFO/SAP>0 && INFO/RPL>1 && INFO/RPR>1" \
| bcftools filter -G 10 -o "${bcftools_out}"

"${vcftools}" --vcf "${bcftools_out}" --thin 5 --recode --recode-INFO-all \
--out "${vcftools_out}"

# annotate
#Annotate using snpeff with built database (remember alternate codons!)
java -Xmx4g -jar "${snpeff}" -c "${snpeff_config}" "${snpeff_db}" \
"${vcftools_out}.recode.vcf" >  "${annotate_vcf}"

#zip and index output (facilitate IGV loading)
bgzip "${annotate_vcf}"
bcftools index -t "${annotate_vcf}".gz
