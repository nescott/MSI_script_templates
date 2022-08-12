#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=2:00:00
#SBATCH -p amdsmall,amdlarge,amd512,amd2tb
#SBATCH -o job_out/%x_%u_%A_%a.out
#SBATCH -e job_out/%x_%u_%A_%a.err
#SBATCH --array=2-20

#run fastqc per dir
 
set -ue
set -o pipefail

sample_file=MEC_Calbicans_samples_paths.txt  #tab delimited sampleID read1 read2
line=${SLURM_ARRAY_TASK_ID} 


#Load modules 
module load fastqc   

#Read sample file line corresponding to array task ID and get variables
sample=$(awk -v val=$line 'NR == val { print $1}' $sample_file)
read1=$(awk -v val=$line 'NR == val { print $2}' $sample_file)
read2=$(awk -v val=$line 'NR == val { print $3}' $sample_file)

cd qc

fastqc -o . $read1 
fastqc -o . $read2

cd ..

