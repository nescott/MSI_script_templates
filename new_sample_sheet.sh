#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH -t 50
#SBATCH -p msismall,msilarge

#Creates sample ID and sequencing path file
#Output is "sample_id path/to/read1.fastq.gz path/to/read2.fastq.gz path/to/any/reads.fastq.gz"

dir=/home/selmecki/shared/disaster_recovery/Sequencing_Runs/
out_file=AMS_sequencing_paths.txt  #name for output

find $dir  -type f -name "*fastq*" | sort | grep -Eiv "RNA|SRA|MinION" \
|tee fastq.txt | xargs basename -s ".fastq.gz" |tee basenames.txt \
| grep -Eo "(AMS|MEC)+_?[0-9]{3,5}" | sort | uniq > samples.txt

awk '
{
    if (substr($NF,1,4)=="_")
        str1=substr($NF,1,8)
    else if (substr($NF,1,3)=="MEC")
        str1=substr($NF,1,6)
    else
        str1=substr($NF,1,7)
}
FNR==NR {
    a[$0]
    next
}
(str1 in a) {
a[str1]=(a[str1]?a[str1] OFS:" ")$0
}
END {
  for (i in a)
      if (a[i]!="")
          print i,a[i]
}
' samples.txt FS="/" fastq.txt | tr -s ' ' | sort -k1,1 > "${out_file}"
