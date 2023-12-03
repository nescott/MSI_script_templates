#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=
#SBATCH --time=2:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --array=1-

# Performs quality trimming with trimmomatic, alignment with bwa mem, and
# post-alignment clean-up and basic stats with samtools.

set -ue
set -o pipefail

sample_file=  #tab delimited sampleID read1 read2
species=  #no spaces in name
instrument=
ref_fasta=  #include path; indices should be same directory
line=${SLURM_ARRAY_TASK_ID}

# Check for/create output directories
arr=("$PWD/trimmed_fastq" "$PWD/logs" "$PWD/bam")
for d in "${arr[@]}"; do
  if [ ! -d "$d" ]; then
    mkdir "$d"
  fi
done

# Load modules for trimming and aligning
module load trimmomatic/0.39
module load bwa/0.7.17
module load samtools/1.10

# Read sample file line corresponding to array task ID and get variables
strain=$(awk -v val="$line" 'NR == val { print $1}' $sample_file)
read1=$(awk -v val="$line" 'NR == val { print $2}' $sample_file)

# Quality trimming
java -jar /panfs/roc/msisoft/trimmomatic/0.39/trimmomatic.jar PE -threads 8 \
-phred33 -trimlog logs/"${strain}".trimlog -basein \
"${read1}" -baseout trimmed_fastq/"${strain}"_trimmed.fastq.gz LEADING:3 TRAILING:3 \
 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33

# Alignment, fix mate-pair errors from alignment, sort, mark duplicates
bwa mem -t 8 -R "@RG\tID:${species}_${strain}\tPL:ILLUMINA\tPM:${instrument}\tSM:${strain}" \
${ref_fasta} trimmed_fastq/"${strain}"_trimmed_1P.fastq.gz trimmed_fastq/"${strain}"_trimmed_2P.fastq.gz \
| samtools fixmate -m - - \
| samtools sort -l 0 -T ${species} -@8 - \
| samtools markdup -@8 - bam/"${strain}"_trimmed_bwa_sorted_markdup.bam

# Reindex
samtools index bam/"${strain}"_trimmed_bwa_sorted_markdup.bam

# Basic stats
samtools flagstat bam/"${strain}"_trimmed_bwa_sorted_markdup.bam \
> logs/"${strain}"_trimmed_bwa_sorted_markdup.stdout
