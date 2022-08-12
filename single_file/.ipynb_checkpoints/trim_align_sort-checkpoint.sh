#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=15gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=2:00:00
#SBATCH -p amdsmall,amdlarge,amd512,amd2tb
#SBATCH -o %x_%u_%j.out
#SBATCH -e %x_%u_%j.err

set -ue
set -o pipefail

species= #no spaces in name
ref_fasta=
read1= #with path if needed

#Parse the fastq file for strain, verify it matches the strain entered manually 
strain=$(basename "${read1}" | cut -d '_' -f 1)

# Load modules for trimming and aligning
module load trimmomatic/0.39
module load bwa/0.7.17
module load samtools/1.10

# Check for/create output directories
arr=("$PWD/trimmed_fastq" "$PWD/logs" "$PWD/bam")
for d in "${arr[@]}"; do
  if [ ! -d "$d" ]; then
    mkdir "$d"
  fi
done

# Quality trimming 
java -jar /panfs/roc/msisoft/trimmomatic/0.39/trimmomatic.jar PE -threads 8 \
-phred33 -trimlog logs/"${strain}".trimlog -basein "${read1}" \
-baseout trimmed_fastq/"${strain}"_trimmed.fastq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33

# Alignment, fix mate-pair errors from alignment, sort, mark duplicates
bwa mem -t 128 -R "@RG\tID:${species}_${strain}\tPL:ILLUMINA\tPM:NextSeq\tSM:${strain}" \
"${ref_fasta}" trimmed_fastq/"${strain}"_trimmed_1P.fastq.gz \
trimmed_fastq/"${strain}"_trimmed_2P.fastq.gz \
| samtools fixmate -m - - \
| samtools sort -l 0 -T "${species}" -@8 - \
| samtools markdup -@8 - bam/"${strain}"_trimmed_bwa_sorted_markdup.bam

#reindex
samtools index bam/"${strain}"_trimmed_bwa_sorted_markdup.bam

#basic stats
samtools flagstat bam/"${strain}"_trimmed_bwa_sorted_markdup.bam \
> logs/"${strain}"_trimmed_bwa_sorted_markdup.stdout

