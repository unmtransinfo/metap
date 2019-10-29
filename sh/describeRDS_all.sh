#!/bin/bash
###

cwd=$(pwd)

IDIR="/home/oleg/workspace/metap/data/input"
DATADIR="${cwd}/data/info"

i=0
N=$(ls $IDIR/*.rds |wc -l)
for fpath in $(ls $IDIR/*.rds) ; do
	i=$(($i + 1))
	fname=$(basename $fpath)
	ofile=${DATADIR}/${fname}.info
	printf "%d/%d. %s => %s\n" $i $N "$fname" "$ofile"
	${cwd}/R/describeRDS.R ${IDIR}/${fname} >&${ofile}
done
