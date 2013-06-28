#!/bin/bash

source ~/.bash_profile
runTypeDefault=0
if test $# -lt 6; then
	echo "Usage:"
	echo "  $0 pegasusWorkFolder jobFilenamePrefix XmxOldValue XmxNewValue MaxPermSizeNewValue totalMemoryNewValue [runType]"
	echo ""
	echo "Note:"
	echo "	#. This program modifies pegasus job scripts (jobFilenamePrefix*.sub) in a workflow working folder."
	echo " 	#. XmxOldValue has to match the original value in order for XmxNewValue to be set."
	echo "  #. if XmxNewValue is <=0, its modification will not be carried out. Same for MaxPermSize and totalMemoryNewValue."
	echo "  #. If runType is 0, this program leaves the new content in jobFilename.tmp but do not overwrite the original one."
	echo "  #. If runType is 1, this program overwrites (by mv command) the original FILEPATH with FILEPATH.tmp (and FILEPATH.tmp will be gone)."
	echo "  #. If runType is 2, this program deletes FILEPATH.tmp left by runType 0."
	echo "  #. Default runType is $runTypeDefault."
	echo
	echo "Examples:"
	echo "	# dry-run: Leave MaxPermSize alone, modify Xmx and total memory request"
	echo "	$0 work/BaseQualityRecalibration/LocalRealignmentBQSR_AlnID2828_2847_vsMethod87.2013.Apr.17T230821/ merge_pegasus-addOrReplaceReadGroupsJava 3500 7000 -1 12000 0"
	echo
	echo "  # dry-run: set MaxPermSize to 8000 (Mb) and total memory request to 12000, leave Xmx unchanged."
	echo "	$0 work/BaseQualityRecalibration/LocalRealignmentBQSR_AlnID2828_2847_vsMethod87.2013.Apr.17T230821/ merge_pegasus-addOrReplaceReadGroupsJava 4000 0 8000 12000 0"
	echo
	echo
	echo
	exit 1
fi

#debug mode
#set -vx
pegasusWorkFolder=$1
jobFilenamePrefix=$2
XmxOldValue=$3
XmxNewValue=$4
MaxPermSizeNewValue=$5
totalMemoryNewValue=$6
runType=$7

if [ -z $runType ]; then
	runType=$runTypeDefault
fi

echo runType is $runType

if test $XmxNewValue -gt 0; then
	~/script/shell/processFile/substitutePatternInAllMatchedFiles.sh $pegasusWorkFolder "$jobFilenamePrefix*.sub"  Xmx$XmxOldValue\m Xmx$XmxNewValue\m 1 $runType 
	~/script/shell/processFile/substitutePatternInAllMatchedFiles.sh $pegasusWorkFolder "$jobFilenamePrefix*.in"  Xmx$XmxOldValue\m Xmx$XmxNewValue\m 1 $runType 
fi

if test $MaxPermSizeNewValue -gt 0; then
	./modifyPegasusJobContent.sh $pegasusWorkFolder $jobFilenamePrefix $MaxPermSizeNewValue $runType MaxPermSize
fi

if test $totalMemoryNewValue -gt 0; then
	./modifyPegasusJobContent.sh $pegasusWorkFolder $jobFilenamePrefix $totalMemoryNewValue $runType request_memory
	
	~/script/pymodule/pegasus/shell/modifyPegasusJobContent.sh $pegasusWorkFolder $jobFilenamePrefix $totalMemoryNewValue $runType "(memory"
fi

date
