#!/bin/sh

noOfWaitingSecondsDefault=5
if test $# -lt 1
	then
	echo "Usage:"
	echo " $0 WORKFLOW1_FOLDER [WORKFLOW2_FOLDER] ...."
	echo
	echo "Note:"
	echo "	#. WORKFLOW1_FOLDER could also be a file that contains workflow folder names (one on each line)"
	echo
	echo "This program shows last 10 lines of pegasus-status for each workflow every $noOfWaitingSecondsDefault seconds."
	echo
	exit
fi

noOfWaitingSeconds=$noOfWaitingSecondsDefault

workflowFolders=$*

while [ "1" = "1" ]; do
	echo
	echo "=================================================================="
	echo -n "==============  current date is  "
	date
	for i in $workflowFolders; do
		if test -d $i; then
			# it is a folder
			workdir=$i
			echo $workdir
			pegasus-status $workdir|tail -n 10
			sleep $noOfWaitingSeconds
		elif test -r $i; then
			#it is a file that contains workflow folder names
			fileWithWorkflowFolderNames=$i
			for workdir in `cat $fileWithWorkflowFolderNames`; do
				echo $workdir
				pegasus-status $workdir|tail -n 10
				sleep $noOfWaitingSeconds
			done
		fi
	done
done
