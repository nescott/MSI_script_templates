#!/bin/bash
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH -o %x_%u.out
#SBATCH -e %x_%u.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=1:00:00
#SBATCH -p amdsmall,amdlarge,amd512,amd2tb,small,large,max

set -ue
set -o pipefail
species=

# Create sample map
find gvcf -name '*.g.vcf' -type f -printf '%p\n' \
 | awk -F'[/_]' 'BEGIN{OFS="\t"} {print $2, $0}'> "${species}".sample_map
