#!/bin/bash
## Author: Hector Rovira & Tom Robinson 
## Golem Wrapper

function write_line {
   while read CHR_FILE CHR_START CHR_END STRAND PROB NORM SCORE MOTIF SEQ
   do
      if [ ${SEQ} ]; then
         ACTUAL_START="`expr ${CHR_START} + ${BP_OFFSET}`"
         ACTUAL_END="`expr ${CHR_END} + ${BP_OFFSET}`"
         echo -e "$CHROMOSOME\t$ACTUAL_START\t$ACTUAL_END\t$STRAND\t$PROB\t$NORM\t$SCORE\t$MOTIF\t$SEQ"
      fi
   done
}

DATASET_HOME=$1
BACKGROUND_FILE=$DATASET_HOME/background/mouse_promoters_1mb_2010.fa
THRESHOLD_FILE=$DATASET_HOME/thresholds/mouse_promoters_1m_least_stringent.tsv

CHROMOSOME=$2
OFFSET_LINE_START=$3
OFFSET_LINE_END=$4
BP_OFFSET=$[OFFSET_LINE_START*50]

SEQ_FILE=$CHROMOSOME"_"$OFFSET_LINE_START"_"$OFFSET_LINE_END".fa"
SOURCE_FILE=$DATASET_HOME"/genomes/masked/"$CHROMOSOME".fa"

echo ">"$CHROMOSOME > $SEQ_FILE
awk -v start=$OFFSET_LINE_START -v end=$OFFSET_LINE_END ' NR>=start && NR <= end' $SOURCE_FILE >> $SEQ_FILE

mlbh_linux $SEQ_FILE $BACKGROUND_FILE $THRESHOLD_FILE $DATASET_HOME/matrices -quiet -filter | write_line
mlbh_linux $SEQ_FILE $BACKGROUND_FILE $THRESHOLD_FILE $DATASET_HOME/matrices -quiet -filter -reverse | write_line

rm $SEQ_FILE

