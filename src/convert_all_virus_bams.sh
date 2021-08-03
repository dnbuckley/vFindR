#!/bin/bash

DIR=$1
BAMTOOLS=$2

cd $DIR
BAMS=`ls -d *bam`

for BAM in $BAMS; do
	BED=`echo $BAM | sed 's/bam$/bed/g'`
	CMD="$BAMTOOLS convert -format bed -in $BAM -out $BED"
	echo $CMD
	$CMD
done

