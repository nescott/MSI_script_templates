#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=1:00:00
#SBATCH -p amdsmall,amdlarge,amd512,amd2tb
#SBATCH -o %x_%u_%j.out
#SBATCH -e %x_%u_%j.err

set -ue
set -o pipefail

species=
ref=  #recognizable abbrevation for reference genome variants were called against
vcf=  #input file name with relative path as needed (e.g., filtered_vcf/file.vcf)
annotate_vcf=  #output annotated file name including .vcf
sift_vcf=   #output final filtered name including .vcf
snpeff_config=/home/selmecki/shared/software/snpEff/snpEff.config #path to snpeff config file
snpeff_db=  #name of built database (no path needed, included in config file)
snpeff=/home/selmecki/shared/software/snpEff/snpEff.jar #path to snpeff jar file
snpSift=/home/selmecki/shared/software/snpEff/SnpSift.jar

# Check for/create output directories
if [ ! -d "annotate" ]; then
  mkdir "annotate"
fi

#Load modules
module load htslib/1.9
module load bcftools/1.10.2

#Annotate using snpeff with built database (remember alternate codons!)
java -Xmx4g -jar "${snpeff}" -v -c "${snpeff_config}" -csvStats \
annotate/"${species}"_"${ref}".csv "${snpeff_db}" "${vcf}" \
> annotate/"${species}"_"${ref}"_"${annotate_vcf}"

# Filter using SnpSift (in addition to GATK hard filters). Change as needed
# Think about homozygous and heterozygous variants
# Joint vs single samples probably need different approaches
java -Xmx4g -jar "${snpSift}" filter \
 "(FILTER = 'PASS') && ((ANN[*].IMPACT has 'HIGH') \
| (ANN[*].IMPACT has 'MODERATE') | (ANN[*].IMPACT has 'LOW'))" \
annotate/"${annotate_vcf}" > \
annotate/"${sift_vcf}"

#zip and index output (facilitate IGV loading)
bgzip annotate_vcf/"${sift_vcf}"
bcftools index -t annotate_vcf/"${sift_vcf}".gz



