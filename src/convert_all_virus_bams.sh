#!/bin/bash

DIR=$1
BAMTOOLS=$2

cd $DIR
BAMS=`ls -d *bam`

for (( i=0; i<${#BAMS[@]}; i++ )); do
	CMD="$BAMTOOLS convert -f bed -in ${BAMS[$i]}"
	echo $CMD
	$CMD
done

