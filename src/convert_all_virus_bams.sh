#!/bin/bash

DIR=$1
BAMTOOLS=$2

cd $DIR
BAMS=`ls -d *bam`

for (( i=0; i<${#BAMS[@]}; i++ )); do
	BAM=${BAMS[$i]}
	BED=`echo $BAM | sed 's/bam$/bed/g'`
	CMD="$BAMTOOLS convert -format bed -in $BAM -out $BED"
	echo $CMD
	$CMD
done

