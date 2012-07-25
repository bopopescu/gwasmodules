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

inputFname=$1
outputFname=$2
gzip -c $inputFname > $outputFname
