#!/bin/bash

IN_FASTA=$1

cat $IN_FASTA \
| paste - - \
| perl -ne '@x=split m/\t/; unshift @x, length($x[1]); print join "\t",@x;' \
| sort -n \
| cut -f2- \
| tr "\t" "\n"