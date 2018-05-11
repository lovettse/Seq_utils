#!/bin/bash

TARGET=$1

paste - - - - - - - - < $TARGET | tee >(cut -f 1-4 | tr "\t" "\n" > `basename $TARGET .${TARGET##*.}`_R1.fastq) | cut -f 5-8 | tr "\t" "\n" > `basename $TARGET .${TARGET##*.}`_R2.fastq;
