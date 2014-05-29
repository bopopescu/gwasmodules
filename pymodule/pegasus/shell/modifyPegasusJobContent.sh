#!/bin/bash

source ~/.bash_profile
attributeNameDefault=TimeToLive
runTypeDefault=0
jobFilenameSuffixDefault=.sub
if test $# -lt 3; then
	echo "Usage:"
	echo "  $0 pegasusWorkFolder jobFilenamePrefix attributeNewValue [runType] [attributeName] [jobFilenameSuffix]"
	echo ""
	echo "Note:"
	echo "	#. This program modifies pegasus job scripts (jobFilenamePrefix*.sub, .sub is controlled by argument jobFilenameSuffix) in the working folder." 
	echo "	#. Right now it is only suited for modifying numerical attributes (re pattern: [0-9]*), like changing 82800 in 'TimeToLive>=82800' to 178000. '>=' could any combination of <>= (at least one) and could have optional space in the front or end, i.e. 'TimeToLive = 82800' or 'TimeToLive=   82800'."
	echo " 	#. Default attributeName is $attributeNameDefault"
	echo "  #. If runType is 0, this program leaves the new content in jobFilename.tmp but do not overwrite the original one."
	echo "  #. If runType is 1, this program overwrites (by mv command) the original FILEPATH with FILEPATH.tmp (and FILEPATH.tmp will be gone)."
	echo "  #. If runType is 2, this program deletes FILEPATH.tmp left by runType 0."
	echo "  #. If runType is 3, same effect as runType 1 but no user-asking, just execute."
	echo "  #. Default runType is $runTypeDefault."
	echo " 	#. Default jobFilenameSuffix is $jobFilenameSuffixDefault."
	echo
	echo "Examples:"
	echo "	# replace old value TimeToLive with 178000 and overwrite original file"
	echo "	$0 work/BaseQualityRecalibration/LocalRealignmentBQSR_AlnID2828_2847_vsMethod87.2013.Apr.17T230821/ merge_pegasus-addOrReplaceReadGroupsJava 178000 1"
	echo
	echo "  # dry-run: replace old MaxPermSize in java jobs (original files not affected), .in files."
	echo "	$0 work/BaseQualityRecalibration/LocalRealignmentBQSR_AlnID2828_2847_vsMethod87.2013.Apr.17T230821/ merge_pegasus-addOrReplaceReadGroupsJava 4000 0 MaxPermSize .in"
	echo
	echo "	# dry-run, change request_memory"
	echo "  $0 work/InspectAlignment/InspectPopulationMonkeyAlignment_RefSeq3488_AlnMethod6_GATKDOC.2013.Jun.21T105136/  DOCWalkerJava_ID00 9040 0 request_memory"
	echo
	echo "  # dry-run, change memory"
	echo ' 	$0 work/InspectAlignment/InspectPopulationMonkeyAlignment_RefSeq3488_AlnMethod6_GATKDOC.2013.Jun.21T105136/  DOCWalkerJava_ID00 9040 0 "(memory"'
	echo
	echo
	exit 1
fi

#debug mode
#set -vx
pegasusWorkFolder=$1
jobFilenamePrefix=$2
attributeNewValue=$3
runType=$4
attributeName=$5
jobFilenameSuffix=$6

if [ -z $runType ]; then
	runType=$runTypeDefault
fi

if [ -z $jobFilenameSuffix ]; then
	jobFilenameSuffix=$jobFilenameSuffixDefault
fi
if test -z "$attributeName"; then
	attributeName=$attributeNameDefault
fi

echo attributeName is $attributeName
echo attributeNewValue is $attributeNewValue
echo runType is $runType


affectedFiles=`ls $pegasusWorkFolder/$jobFilenamePrefix*$jobFilenameSuffix`

echo "These files have matches: "
echo  $affectedFiles

counter=0

if test "$runType" = "1"; then
	echo "File overwriting option is on."
	echo -n "Continue? (y/n): "
	read answer
	if [ -z $answer ]; then
		exit 1
	fi
	if [ $answer != 'y' ]; then
		exit 1
	fi
fi

for f in `ls $pegasusWorkFolder/$jobFilenamePrefix*$jobFilenameSuffix`; do
	sed 's/'$attributeName'\([ ]*[<>=][<>=]*[ ]*\)[0-9]*/'$attributeName'\1'$attributeNewValue'/' $f >$f.tmp;
	counter=`echo $counter+1|bc`
	if test "$runType" = "1" || test "$runType" = "3"; then
		mv $f.tmp $f
	elif test "$runType" = "2"; then
		rm $f.tmp
	fi
done

echo "Modified $counter files."
