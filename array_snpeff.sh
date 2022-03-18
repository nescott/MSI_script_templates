#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=1:00:00
#SBATCH -p amdsmall,amdlarge,amd512,amd2tb
#SBATCH -o job_out/%x_%u_%A_%a.out
#SBATCH -e job_out/%x_%u_%A_%a.err
#SBATCH --array=1-20

set -ue
set -o pipefail

line=${SLURM_ARRAY_TASK_ID} 
sample_file=  #tab delimited file: sample read1 read2
vcf=  #input file name not including strain, e.g. filtered.vcf
annotate_vcf=  #output file name not including strain, e.g. filtered_ann.vcf
snpeff=/home/selmecki/shared/software/snpEff/snpEff.jar #path to snpeff jar file
snpeff_config=/home/selmecki/shared/software/snpEff/snpEff.config #path to snpeff config file
snpeff_db=  #name of built database (no path needed, included in config file)

# Check for/create output directories
if [ ! -d "annotate" ]; then
  mkdir "annotate"
fi

#Load modules
module load htslib/1.9
module load bcftools/1.10.2

#Get strain ID froms sample file line equal to array task ID
strain=$(awk -v val="${line}" 'NR == val { print $1}' "${sample_file}")

#Annotate using snpeff with built database (remember alternate codons!)
java -Xmx4g -jar "${snpeff}" -v -c "${snpeff_config}" -csvStats annotate/"${strain}".csv "${snpeff_db}" vcf/"${strain}"_"${vcf}" > annotate/"${strain}"_"${annotate_vcf}"

#zip and index output (facilitate IGV loading)
bgzip vcf/"${strain}"_"${annotate_vcf}"
bcftools index -t vcf/"${strain}"_"${annotate_vcf}".gz
