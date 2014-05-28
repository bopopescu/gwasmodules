#!/bin/sh

if test $# -le 3; then
	echo ""
	echo "Usage: $0 inputFname1 inputFname2 outputFname"
	echo ""
	echo "Note:"
	echo "	bedtools intersect creates a .bed file with lots of redundant entries."
	echo "	This script clears that redundancy."
	echo "	Two temporary files (.intersect.bed, .sort.bed) will be created and deleted in the end."
	echo ""
	echo "Examples:"
	echo "	$0 ~/RefGenomes/dustPlus10_M1-22XY.bed.gz genomicSuperDups_hg19.bed ~/RefGenomes/dustPlus10_M1-22XY.overlap.genomicSuperDups_hg19.bed "
	echo ""
	exit
fi
inputFname1=$1
inputFname2=$2
outputFname=$3
intersectTmpFname=$outputFname.intersect.bed
sortTmpFname=$outputFname.sort.bed

bedtoolsPath=bedtools
date
echo -n "Intersecting ..."
$bedtoolsPath intersect -a $inputFname1 -b $inputFname2 > $intersectTmpFname
echo "Done"

echo -n "Sorting ..."
sort -k 1,1 -k 2,2n $intersectTmpFname > $sortTmpFname
echo "Done"
echo -n "Merging sorted .bed files ..."
$bedtoolsPath merge -i $sortTmpFname > $outputFname 
echo "Done"

echo -n "Removing temp files ..."
rm $intersectTmpFname
rm $sortTmpFname

echo "Done"
date
