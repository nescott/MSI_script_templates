#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=20 
#SBATCH -p amdsmall,amdlarge,amd512,amd2tb
#SBATCH -o %x_%u_%j.out
#SBATCH -e %x_%u_%j.err

job1=$(sbatch --parsable haplotypecaller_1.sh)
job2=$(sbatch --parsable --dependency=afterok:${job1} samplemap.sh)
job3=$(sbatch --parsable --dependency=afterok:${job1}:${job2} genomicsdb.sh)
job4=$(sbatch --parsable --dependency=afterok:${job1}:${job2}:${job3} genotype_filter.sh)
