#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=10:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o _%j.out
#SBATCH -e %j.err

set -ue
set -o pipefail

# local modules

module use /home/selmecki/shared/software/modulefiles.local
module load flye/20230927

in_reads=
out_dir=
genome_size=
coverage=

flye --pacbio-hifi "${in_reads}" --out-dir "${out_dir}" --threads 4 -g "${genome_size}" --asm-coverage "${coverage}"
