#!/bin/bash

INPUT=$1
CORES=$2

BASE_IN=`echo $INPUT | sed -e "s/\.fastq.gz//g"`

unpigz -p $CORES -c ${INPUT} | awk '{print (NR%4 == 1) ? "@" ++i : $0}' | pigz -p $CORES -c > ${BASE_IN}_renamed.fastq.gz
