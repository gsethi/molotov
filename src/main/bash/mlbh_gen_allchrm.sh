#!/bin/bash
DATASET_HOME=$1
SPLIT_SIZE=$2
FILE_NAME=joblist_$SPLIT_SIZE.txt

mlbh_gen_joblist.sh $DATASET_HOME chr1 $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chr2 $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chr3 $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chr4 $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chr5 $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chr6 $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chr7 $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chr8 $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chr9 $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chr10 $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chr11 $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chr12 $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chr13 $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chr14 $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chr15 $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chr16 $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chr17 $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chr18 $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chr19 $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chrM $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chrX $SPLIT_SIZE >> $FILE_NAME
mlbh_gen_joblist.sh $DATASET_HOME chrY $SPLIT_SIZE >> $FILE_NAME
