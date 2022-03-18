#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=1:00:00
#SBATCH -p amdsmall,amdlarge,amd512,amd2tb
#SBATCH -o job_out/%x_%u_%A_%a.out
#SBATCH -e job_out/%x_%u_%A_%a.err
#SBATCH --array=1-20

#Specific to C. lusitaniae GCA_000003835:
#Annotate using snpeff and downloaded GCA_000003835 database
#Downloaded database required correcting chromosome name (Snpeff database  uses scaffold_x, vcf has genbank accession) and didn't use correct codon tables
#Built database using genbank accession fna and gtf, but specifying codon tables

set -ue
set -o pipefail

line=${SLURM_ARRAY_TASK_ID} 
cur=$PWD
sample_file=MEC_Clusitaniae_paths.txt
snpeff=/home/selmecki/shared/software/snpEff/snpEff.jar
config=/home/selmecki/shared/software/snpEff/snpEff.config
DB=ATCC42720_GCA000003835
snpSift=/home/selmecki/shared/software/snpEff/SnpSift.jar

#Load modules
module load htslib/1.9
module load bcftools/1.10.2

#Get strain ID froms sample file line equal to array task ID
strain=$(awk -v val="$line" 'NR == val { print $1}' $sample_file)

#check for sample-specific output directory, move into or fail if not found
if [ -d "${strain}_out" ]; then
  cd "${strain}"_out
else
  echo >&2 "output folder and files don't exist for variant calling"
  exit 1
fi

#annotate (vcf from mutect2 with ATCC normal) and filter (passed Mutect2 and has impact other than just modifier)

srun java -jar "${snpeff}" -v -c "${config}" -csvStats "${strain}"_ATCCnormal.csv -s "${strain}"_ATCCnormal.html \
-ud 250 "${DB}" "${strain}"_ATCCnormal_filtered.vcf > "${strain}"_ATCCnormal_filtered_ann.vcf    

srun java -jar "${snpSift}"  filter \
"(FILTER = 'PASS') & ((ANN[*].IMPACT has 'HIGH') | (ANN[*].IMPACT has 'MODERATE') | (ANN[*].IMPACT has 'LOW'))" \
"${strain}"_ATCCnormal_filtered_ann.vcf > "${strain}"_ATCCnormal_filtered_ann_pass_impacts.vcf

srun bgzip "${strain}"_ATCCnormal_filtered_ann_pass_impacts.vcf
srun bcftools index -t "${strain}"_ATCCnormal_filtered_ann_pass_impacts.vcf.gz

#annotate (vcf from mutect2 with AMS5200 normal)
if [[ ${line} != 1 ]]; then
  srun java -jar "${snpeff}" -v -c "${config}" -csvStats "${strain}"_AMSnormal.csv -s "${strain}"_AMSnormal.html -ud 250 \
  "${DB}"  "${strain}"_AMSnormal_filtered.vcf > "${strain}"_AMSnormal_filtered_ann.vcf
  
  srun java -jar "${snpSift}" filter \
 "(FILTER = 'PASS') & ((ANN[*].IMPACT has 'HIGH') | (ANN[*].IMPACT has 'MODERATE') | (ANN[*].IMPACT has 'LOW'))" \
  "${strain}"_AMSnormal_filtered_ann.vcf > "${strain}"_AMSnormal_filtered_ann_pass_impacts.vcf
fi

#back to parent directory
cd "$cur"
