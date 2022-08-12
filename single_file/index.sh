#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=index_%j.out
#SBATCH --error=index_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH -t 20
#SBATCH -p amdsmall,small

reference_fasta=

#Load modules
module load bwa/0.7.17
module load picard/2.25.6
module load samtools/1.10
module load htslib

#bgzip if needed; samtools doesn't like gzip
if [[ $reference_fasta =~ \.gz$ ]]; then
  gunzip "${reference_fasta}"
  bgzip  "$(basename "${reference_fasta}" .gz)"
fi

#Generate indices
bwa index -a bwtsw "${reference_fasta}"

java -jar /panfs/roc/msisoft/picard/2.25.6/picard.jar \
 CreateSequenceDictionary -R "${reference_fasta}"

samtools faidx "${reference_fasta}"




