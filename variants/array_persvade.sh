#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=66gb
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --output=%A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH -p msismall,msibigmem
#SBATCH --array=25,39,42,45,50,55,56,64,71,73,86,87,93,97

set -ue
set -o pipefail

line=${SLURM_ARRAY_TASK_ID}
sample_file=bam.files
ploidy=2
cnv_window=500
ref_fasta=C_albicans_SC5314_version_A21-s02-m09-r08_chromosomes.fasta
gff_file=C_albicans_SC5314_A21_current_features.gff
mito=Ca19-mtDNA
genetic_code=12
m_code=3
sv_sims=breakpoints_aroundHomRegions_wsize=500bp_maxEval=1e-05_minQcovS=50.bedpe
sim_p=diploid_hetero
repeat_file=combined_repeats.tab
threads=32
mem_frac=0.128

bam=$(awk -v val="$line" 'NR == val { print $1}' $sample_file)
strain=$(basename "$bam" | cut -d "_" -f 1)

mkdir -p "${strain}"/{params,sv,cnv,genecov,integrate,annotate}
mkdir -p /scratch.global/scot0854/persvade/"${strain}"

export PERSVADE_TMPDIR=/scratch.global/scot0854/persvade/"${strain}"

module load singularity

singularity exec -e mikischikora_persvade_v1.02.6.sif \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate \
    perSVade_env && python /perSVade/scripts/perSVade \
    optimize_parameters  \
    --ref ${ref_fasta} \
    -sbam ${bam} \
    -mchr ${mito}\
    --repeats_file ${repeat_file} \
    --regions_SVsimulations ${sv_sims} \
    --simulation_ploidies ${sim_p} \
    -o ${strain}/params \
    --thr ${threads} \
    --fraction_available_mem ${mem_frac}"

singularity exec -e mikischikora_persvade_v1.02.6.sif \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate \
    perSVade_env && python /perSVade/scripts/perSVade \
    call_SVs  \
    --SVcalling_parameters ${strain}/params/optimized_parameters.json \
    --ref ${ref_fasta} \
    -sbam ${bam} \
    -mchr ${mito} \
    --repeats_file ${repeat_file} \
    -o ${strain}/sv \
    --thr ${threads} \
    --fraction_available_mem ${mem_frac}"

singularity exec -e mikischikora_persvade_v1.02.6.sif \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate \
    perSVade_env && python /perSVade/scripts/perSVade \
    call_CNVs  \
    --cnv_calling_algs HMMcopy,AneuFinder \
    -p ${ploidy} \
    --window_size_CNVcalling ${cnv_window} \
    --ref ${ref_fasta} \
    -sbam ${bam} \
    -mchr ${mito} \
    -o ${strain}/cnv \
    --thr ${threads} \
    --fraction_available_mem ${mem_frac}"

singularity exec -e mikischikora_persvade_v1.02.6.sif \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate \
    perSVade_env && python /perSVade/scripts/perSVade \
    get_cov_genes \
    --ref ${ref_fasta} \
    -gff  ${gff_file} \
    -sbam ${bam} \
    -o ${strain}/genecov \
    --thr ${threads} \
    --fraction_available_mem ${mem_frac}"

singularity exec -e mikischikora_persvade_v1.02.6.sif \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate \
    perSVade_env && python /perSVade/scripts/perSVade \
    integrate_SV_CNV_calls  \
    -p ${ploidy} \
    --ref ${ref_fasta} \
    -sbam ${bam} \
    -mchr ${mito} \
    --repeats_file ${repeat_file} \
    --outdir_callSVs ${strain}/sv \
    --outdir_callCNVs ${strain}/cnv \
    -o ${strain}/integrate \
    --thr ${threads} \
    --fraction_available_mem ${mem_frac}"

singularity exec -e mikischikora_persvade_v1.02.6.sif \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate \
    perSVade_env && python /perSVade/scripts/perSVade \
    annotate_SVs \
    --ref ${ref_fasta} \
    -gff  ${gff_file} \
    -mchr ${mito} \
    -mcode ${m_code} \
    -gcode ${genetic_code} \
    --SV_CNV_vcf  ${strain}/integrate/SV_and_CNV_variant_calling.vcf \
    -o ${strain}/annotate \
    --thr ${threads} \
    --fraction_available_mem ${mem_frac}"
