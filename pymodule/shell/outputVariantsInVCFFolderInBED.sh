#!/bin/bash

defaultFormatType=1
if test $# -lt 2 ; then
	echo "$0 INPUT_VCF_Folder OutputFile [FormatType]"
	echo
	echo "This script takes all the variants in all VCF files in the VCF folder and then outputs them in BED format, tab-delimited, chr, start, stop."
	echo "	NOTE: BED file is 0-based and stop-th base is not included. i.e. start=0, stop=100 selects bases from 0 to 99."
	echo
	echo "FormatType, what type of format. default is $defaultFormatType."
	echo "	1: BED"
	echo "	2: tab-delimited, 3-column: chr start stop. coordinates are 1-based. start=stop for SNPs. could be used for vcftools to keep sites."
	exit 1
fi
shellDir=`dirname $0`
source $shellDir/common.sh

inputFolder=$1
outputFilename=$2
formatType=$3
if test -z "$formatType"; then
	formatType=$defaultFormatType
fi
echo FormatType is $formatType.
n=0
exitIfFileExists $outputFilename

for i in `ls $inputFolder/*vcf.gz`; do
	if test $formatType -eq 1; then
		zcat $i |grep -v ^#|awk -F ' ' '{print $1 "\t" $2-1 "\t" $2}' >> $outputFilename
	elif test $formatType -eq 2; then
		zcat $i |grep -v ^#|awk -F ' ' '{print $1 "\t" $2 "\t" $2}' >> $outputFilename
	fi
done
