#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --mem=10gb
#SBATCH -o %x_%u.out
#SBATCH -e %x_%u.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=6:00:00
#SBATCH -p amdsmall,amdlarge,amd512,amd2tb,small,large,max
#SBATCH --array=1-

set -ue
set -o pipefail

species=
ref=

mkdir -p db

module load gatk

gatk --java-options "-Xmx6g" GenomicsDBImport \
  --genomicsdb-workspace-path db/"${species}_${ref}" \
  -L intervals.list \
  --sample-name-map "${species}".sample_map
  --validate-sample-name-map true
