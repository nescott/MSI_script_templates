#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=25gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=
#SBATCH --time=4:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --array=

# For a set of samples, perform adapter and quality trimming, align to selected reference genome,
# sort and mark duplicates in BAM files, and generate basic alignment stats
# input: tab delimited file of sample ID, read 1 and read 2
# array size == number of samples

set -ue
set -o pipefail

line=${SLURM_ARRAY_TASK_ID}
sample_file=
# make a global scratch directory with your x500 id, then a relevant subdirectory
tempdir=  # with trailing slash
species=  # no spaces in name
ref_fasta=  # can be zipped or unzipped, must be indexed

# Read sample file line corresponding to array task ID and get variables
strain=$(awk -v val="$line" 'NR == val { print $1}' $sample_file)
read1=$(awk -v val="$line" 'NR == val { print $2}' $sample_file)
read2=$(awk -v val="$line" 'NR == val { print $3}' $sample_file)

# Load modules for trimming and aligning
module use /home/selmecki/shared/software/modulefiles.local

module load bbmap
module load bwa/0.7.17
module load samtools/1.10

# Check for/create output directories
arr=("${tempdir}trimmed_fastq" "logs" "bam")
for d in "${arr[@]}"; do
  if [ ! -d "$d" ]; then
    mkdir -p "$d"
  fi
done

# JGI BBTools data preprocessing guidelines:
## trim adapters
bbduk.sh in1="${read1}" in2="${read2}" \
out1="${tempdir}"trimmed_fastq/"${strain}"_trim_adapt1.fq \
out2="${tempdir}"trimmed_fastq/"${strain}"_trim_adapt2.fq \
ref=adapters ktrim=r k=23 mink=11 hdist=1 ftm=5 tpe tbo

## contaminant filtering per bbduk user guide
bbduk.sh in1="${tempdir}"trimmed_fastq/"${strain}"_trim_adapt1.fq \
in2="${tempdir}"trimmed_fastq/"${strain}"_trim_adapt2.fq \
out1="${tempdir}"trimmed_fastq/"${strain}"_unmatched1.fq \
out2="${tempdir}"trimmed_fastq/"${strain}"_unmatched2.fq \
outm1="${tempdir}"trimmed_fastq/"${strain}"_matched1.fq \
outm2="${tempdir}"trimmed_fastq/"${strain}"_matched2.fq \
ref=phix,artifacts k=31 hdist=1 stats=logs/"${strain}"_phistats.txt

## quality trimming (bbduk user guide recommends this as separate step from adapter trimming)
bbduk.sh in1="${tempdir}"trimmed_fastq/"${strain}"_unmatched1.fq \
in2="${tempdir}"trimmed_fastq/"${strain}"_unmatched2.fq \
out1="${tempdir}"trimmed_fastq/"${strain}"_trimmed_1P.fq \
out2="${tempdir}"trimmed_fastq/"${strain}"_trimmed_2P.fq \
qtrim=rl trimq=10

# Reference alignment, fix mate-pair errors from alignment, sort, mark duplicates
# Including sample information to ensure unique read groups if freebayes is used
bwa mem -t 8 -R "@RG\tID:${species}_${strain}\tPL:ILLUMINA\tPM:NextSeq\tSM:${strain}" \
"${ref_fasta}" "${tempdir}"trimmed_fastq/"${strain}"_trimmed_1P.fq \
"${tempdir}"trimmed_fastq/"${strain}"_trimmed_2P.fq \
| samtools fixmate -m - - \
| samtools sort -l 0 -T "${species}" -@8 - \
| samtools markdup -@8 - bam/"${strain}"_trimmed_bwa_sorted_markdup.bam

# reindex
samtools index bam/"${strain}"_trimmed_bwa_sorted_markdup.bam

# basic stats
samtools flagstat bam/"${strain}"_trimmed_bwa_sorted_markdup.bam \
> logs/"${strain}"_trimmed_bwa_sorted_markdup.stdout
