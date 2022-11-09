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
#SBATCH --array=1-

set -ue
set -o pipefail

species= #no spaces in name
ref_fasta=
sample_file=  #tab delimited sampleID read1 read2
adapters=/home/selmecki/shared/software/bbmap/resources/adapters.fa
phix=/home/selmecki/shared/software/bbmap/resources/adapters.fa
bbduk=/home/selmecki/shared/software/bbmap/bbduk.sh
line=${SLURM_ARRAY_TASK_ID}

#Read sample file line corresponding to array task ID and get variables
strain=$(awk -v val="$line" 'NR == val { print $1}' $sample_file)
read1=$(awk -v val="$line" 'NR == val { print $2}' $sample_file)
read2=$(awk -v val="$line" 'NR == val { print $3}' $sample_file)

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
"${bbduk}" in1="${read1}" in2="${read2}" \
out1=trimmed_fastq/"${strain}"_trim_adapt1.fq out2=trimmed_fastq/"${strain}"_trim_adapt2.fq \
ref="${adapters}" ktrim=r k=23 mink=11 hdist=1 ftm=5 tpe tbo

## contaminant (phix) filtering
"${bbduk}" in1=trimmed_fastq/"${strain}"_trim_adapt1.fq in2=trimmed_fastq/"${strain}"_trim_adapt2.fq \
out1=trimmed_fastq/"${strain}"_unmatched1.fq out2=trimmed_fastq/"${strain}"_unmatched2.fq \
outm1=trimmed_fastq/"${strain}"_matched1.fq outm2=trimmed_fastq/"${strain}"_matched2.fq \
ref="${phix}" k=31 hdist=1 stats=phistats.txt

## quality trimming (bbduk user guide recommends this as separate step from adapter trimming)
"${bbduk}" in1=trimmed_fastq/"${strain}"_unmatched1.fq in2=trimmed_fastq/"${strain}"_unmatched2.fq \
out1=trimmed_fastq/"${strain}"_trimmed_1P.fq out2=trimmed_fastq/"${strain}"_trimmed_2P.fq qtrim=rl trimq=10

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

