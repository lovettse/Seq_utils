#!/bin/bash

awk '/^>/ {OUT=substr($0,2) ".fasta"}; {print >> OUT; close(OUT)}' $1
