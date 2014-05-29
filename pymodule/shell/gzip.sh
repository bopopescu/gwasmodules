#!/bin/bash
# 2012-7.19 wrap around gzip so that input of gzip won't be deleted
if test $# -lt 2
then
	echo "Usage: $0 inputFile gzipOutputFile"
	echo
	echo "Note:"
	echo "	the gzipOutputFile should have .gz as suffix."
	echo
	echo "Example:"
	echo "	$0 input.txt input.txt.gz"
exit
fi

shellDir=`dirname $0`
source $shellDir/common.sh

inputFname=$1
outputFname=$2

exitIfFileExists $outputFname

#2013.11.22 use bgzip in case tabix indexing is needed.
~/bin/bgzip -c $inputFname > $outputFname
