#!/bin/bash
# 2013.2.4

trunkSizeDefault=50000

if test $# -lt 3
then
	echo "Usage: $0 psmcFolderPath psmcInputFname psmcSplitInputFname [trunkSize]"
	echo
	echo "Notes:"
	echo "	#. psmcInputFname is input of psmc (Li Durbin, 2011)."
	echo "	#. psmcSplitInputFname is split input of psmc."
	echo "	#. trunkSize is $trunkSizeDefault by default."
	echo "	#.   It is the number of characters per chunk. Since each character is a '0','1' or '.' per 100 bp, each trunk contains 100bp*trunkSize bases."
	echo
	echo "Example:"
	echo "	$0 ~/script/psmc/ diploid.psmcfa diploidSplit.psmcfa 50000"
exit
fi

source $HOME/.bash_profile

shellDir=~/script/shell/
source $shellDir/common.sh
psmcFolderPath=$1
psmcInputFname=$2
psmcSplitInputFname=$3
trunkSize=$4
if [ -z $trunkSize ]
then
	trunkSize=$trunkSizeDefault
fi

splitfaPath=$psmcFolderPath/utils/splitfa
###
#Usage: splitfa <in.fa> [trunk_size=500000]

$splitfaPath $psmcInputFname $trunkSize > $psmcSplitInputFname
