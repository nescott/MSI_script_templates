#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=15gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=2:00:00
#SBATCH -p msismall,msilarge 
#SBATCH -o %x_%u_%j.out
#SBATCH -e %x_%u_%j.err

set -ue
set -o pipefail

species= #no spaces in name
ref_fasta=
raw_fq1=
raw_fq2=
strain=
adapters=/home/selmecki/shared/software/bbmap/resources/adapters.fa
phix=/home/selmecki/shared/software/bbmap/resources/adapters.fa
bbduk=/home/selmecki/shared/software/bbmap/bbduk.sh
read1= #with path if needed

#Parse the fastq file for strain, verify it matches the strain entered manually 
strain=$(basename "${read1}" | cut -d '_' -f 1)

# Load modules for trimming and aligning
module load bwa/0.7.17
module load samtools/1.10

# Check for/create output directories
arr=("$PWD/trimmed_fastq" "$PWD/logs" "$PWD/bam")
for d in "${arr[@]}"; do
  if [ ! -d "$d" ]; then
    mkdir "$d"
  fi
done

# Adapter and quality trimming using JGI BBTools data preprocessing guidelines 
## trim adapters 
"${bbduk}" in1="${raw_fq1}" in2="${raw_fq2}" out1="${strain}"_trim_adapt1.fq out2="${strain}"_trim_adapt2.fq ref="${adapters}" ktrim=r k=23 mink=11 hdist=1 ftm=5 tpe tbo

## contaminant (phix) filtering
"${bbduk}" in1="${strain}"_trim_adapt1.fq in2="${strain}"_trim_adapt2.fq out1="${strain}"_unmatched1.fq out2="${strain}"_unmatched2.fq outm1="${strain}"_matched1.fq outm2="${strain}"_matched2.fq ref="${phix}" k=31 hdist=1 stats=phistats.txt

## quality trimming (bbduk user guide recommends this as separate step from adapter trimming)
"${bbduk}" in1="${strain}"_unmatched1.fq in2="${strain}"_unmatched2.fq out1="${strain}"_trimmed_1P.fq out2="${strain}"_trimmed_2P.fq qtrim=rl trimq=10

# Reference alignment, fix mate-pair errors from alignment, sort, mark duplicates
# Including sample information to ensure unique read groups if freebayes is used
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

