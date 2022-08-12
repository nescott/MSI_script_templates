#!/bin/bash

#Creates sample input file for trim_align_sort script
#Output is "sample_id\tpath/to/read1.fastq.gz\tpath/to/read2.fastq.gz"

dir=$1 #requires trailing slash on directory
samplefile=$2 #name for output
todo_list=$3 #list of sample IDs to find

subdircount=$(find $dir -maxdepth 1 -type d | wc -l)

if [[ "$subdircount" -eq 1 ]]
then
  files=$(find . -mindepth 1 -maxdepth 1 -type f | sort | tr '\n' '\t')
  sample=$(basename "${files}" | cut -d '_' -f 1)
        
  if [[ $sample =~ [0-9] ]]; then
    sample_id=$sample
  else
    sample_id=$(basename "${files}" | cut -d '_' -f 2)
  fi
  
  if grep -Fq  "$sample_id" "${todo_list}"; then
    printf "%s\t%s\n" "${sample_id}" "${files}" >> "${samplefile}"
  fi
else
  for subdir in "${dir}"*
    do
      files=$(find "${subdir}" -mindepth 1 -maxdepth 1 -type f | sort | tr '\n' '\t')
      sample=$(basename "${files}" | cut -d '_' -f 1)
        
    if [[ $sample =~ [0-9] ]]; then
      sample_id=$sample
    else
      sample_id=$(basename "${files}" | cut -d '_' -f 2)
    fi 
    
    if grep -Fq  "$sample_id" "${todo_list}"; then
      printf "%s\t%s\n" "${sample_id}" "${files}" >> "${samplefile}"
    fi
  done
fi