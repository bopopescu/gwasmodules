#!/bin/sh
if test $# -lt 2 ; then
	echo ""
	echo "Usage: $0 inputFname outputFname"
	echo ""
	echo "Note:"
	echo "	This script clears that redundancy."
	echo "	Two temporary files (.intersect.bed, .sort.bed) will be created and deleted in the end."
	echo ""
	echo "Examples:"
	echo "	$0 ~/RefGenomes/dustPlus10_M1-22XY.bed.gz ~/RefGenomes/dustPlus10_M1-22XY.uniqueMerged.bed "
	echo ""
	exit
fi
inputFname=$1
outputFname=$2
sortTmpFname=$outputFname.sort.bed

bedtoolsPath=bedtools
date

echo -n "Sorting ..."
#2014.05.28 problem with sort is that it can't tell if input has header or not
#sort -k 1,1 -k 2-3n $intersectTmpFname > $sortTmpFname
#2024.05.28 problem with bedtools sort is that it doesnot use the 3rd column (stop) in sorting, only column 1,2
$bedtoolsPath sort -i $inputFname > $sortTmpFname
echo "Done"

echo -n "Merging sorted .bed files ..."
$bedtoolsPath merge -i $sortTmpFname > $outputFname 
echo "Done"

echo -n "Removing temp files ..."
rm $sortTmpFname
echo "Done"
date
