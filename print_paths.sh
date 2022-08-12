#!/bin/bash

#Creates sample input file for trim_align_sort script
#Output is "sample_id\tpath/to/read1.fastq.gz\tpath/to/read2.fastq.gz"

dir=("/home/selmecki/shared/CUSOM/SharedFiles/Sequencing_Data/2021_02_16_MiGS/2021_02_26_MiGS/" "/home/selmecki/shared/CUSOM/SharedFiles/Sequencing_Data/2021_03_01_MiGS/parapsilosis_fastq" "/home/selmecki/shared/CUSOM/SharedFiles/Sequencing_Data/2021_03_01_MiGS/albicans_fastq/2021_MEC_clinical_isolates" "/home/selmecki/shared/Sequencing_runs/AnnaSelmecki21104/2021_MEC_clinical_isolates/" "/home/selmecki/shared/Sequencing_runs/AnnaSelmecki220323/2022_MEC_clinical_isolates/" "/home/selmecki/shared/Sequencing_runs/AnnaSelmecki220406/2022_MEC_clinical_isolates/" "/home/selmecki/shared/Sequencing_runs/AnnaSelmecki220618/2022_MEC_clinical_isolates/" "/home/selmecki/shared/Sequencing_runs/AnnaSelmecki220625/2022_MEC_clinical_isolates/") 
samplefile=all_paths.txt
todo_list=samples.txt

for d in "${dir[@]}"; do
  files=$(find "${subdir}" -mindepth 1 -maxdepth 1 -type f | sort) 
done

for f in $files; do
  sample=$(basename "${files}" | cut -d '_' -f 1)

  if grep -Fq  "$sample" "${todo_list}"; then
    printf "%s\t%s\n" "${sample}" "${files}" >> "${samplefile}"
  fi
done
