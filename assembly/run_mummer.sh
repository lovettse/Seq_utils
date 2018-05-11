#!/bin/bash

REF=$1
QUERY=$2
OUT_BASE=$3

nucmer --maxmatch -c 100 -p $OUT_BASE $REF $QUERY 
show-coords -c $OUT_BASE.delta > $OUT_BASE.coords
mummerplot -p $OUT_BASE --filter --png $OUT_BASE.delta -R $REF -Q $QUERY --layout
