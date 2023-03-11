#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=10:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %x_%u_%j.out
#SBATCH -e %x_%u_%j.err

set -ue
set -o pipefail

# local modules

module use "$HOME"/modulefiles.local
module load flye/20230228 

in_reads=/home/selmecki/shared/disaster_recovery/Sequencing_Runs/AnnaSelmecki221005/MEC225/MEC225/MEC225.hifi_reads.fastq.gz
out_dir=flye

flye --pacbio-hifi "${in_reads}" --out-dir "${out_dir}" --threads 4 -g 13m
