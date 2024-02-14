#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=3gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=3:00:00
#SBATCH -p amdlarge,large
#SBATCH -o job_out/%x_%u_%A_%a.out
#SBATCH -e job_out/%x_%u_%A_%a.err
#SBATCH --array=1-

#Perform variant calling on sorted, markduped bams using Mutect2
#Depends on finding existing output dir and bam per strain listed in sample file

set -ue
set -o pipefail

ref_fasta=
sample_file=
line=${SLURM_ARRAY_TASK_ID}

mkdir -p "$PWD/vcf" "$PWD/logs" "$PWD/bam"

#Load GATK module
module load gatk/4.1.2

#Get strain ID froms sample file line equal to array task ID
strain=$(awk -v val="${line}" 'NR == val { print $1}' "${sample_file}")

#variant calling and filtering, with normal sample
gatk Mutect2 -R "${ref_fasta}" -I bam/"${strain}"_trimmed_bwa_sorted_markdup.bam \
-O vcf/"${strain}"_unfiltered.vcf

gatk FilterMutectCalls -R "$ref_fasta" -V vcf/"${strain}"_unfiltered.vcf \
-O vcf/"${strain}"_filtered.vcf
