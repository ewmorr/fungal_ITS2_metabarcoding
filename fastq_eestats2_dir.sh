#!/usr/bin/bash

#A directory of directories is provided to this script as input
# usearch fastq_eestats will be run on each file in those directories
dir=$1
for d in $dir/*/*pair.fastq

do(
file=${d##*/} #filename
base=${file%.*} #base name
dir2=${d%/*} #path to file

echo $base
echo $dir2

#usearch10 -fastq_eestats2 $d -output $dir2/$base._eestats.txt -length_cutoffs 100,300,25
#usearch10 -fastx_truncate $d  -trunclen 250 -fastqout $base_250.fastq
usearch10 -fastq_filter $d -fastq_truncqual 2 -fastq_maxee 2 -fastqout $dir2/$base.maxee2.fastq
)
done

