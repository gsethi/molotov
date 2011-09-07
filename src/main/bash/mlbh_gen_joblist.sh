#!/bin/bash
# NOTE: This version computes a split factor without overlapping segments. 
# As a result, any motifs crossing split boundaries WILL NOT BE COMPUTED.
DATASET=$1
CHROMOSOME=$2

number_of_lines=`wc -l $DATASET/genomes/masked/$CHROMOSOME.fa | awk '{print $1}'`
line_split=$3
line_start=2
line_end=$line_split

while [ $line_end -le $number_of_lines ]
do  
    echo "1 mlbh_scanMotif.sh" $DATASET $CHROMOSOME $line_start $line_end
    line_start=$line_end
    line_end=$[line_end+line_split]
done

echo "1 mlbh_scanMotif.sh" $DATASET $CHROMOSOME $line_start $number_of_lines

