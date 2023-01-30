#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=260gb
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH -p msismall,msilarge
#SBATCH --array=1-

sample_file=  #tab delimited sampleID read1 read2
line=${SLURM_ARRAY_TASK_ID} 
centrifuge_home=/home/selmecki/shared/software/centrifuge
db=/home/selmecki/shared/centrifuge_fungi_refseq/ncbi_nt/nt

#Read sample file line corresponding to array task ID and get variables
strain=$(awk -v val="$line" 'NR == val { print $1}' $sample_file)
read1=$(awk -v val="$line" 'NR == val { print $2}' $sample_file)
read2=$(awk -v val="$line" 'NR == val { print $3}' $sample_file)

srun "${centrifuge_home}"/centrifuge -p 10 -x "${db}" -1 "${read1}" -2 "${read2}" \
-S "${strain}"_centrifuge.txt

srun "${centrifuge_home}"/centrifuge-kreport -x "${db}" "${strain}"_centrifuge.txt > "${strain}"_kreport.txt

