#!/bin/sh

noOfWaitingSecondsDefault=5
if test $# -lt 1
	then
	echo "Usage:"
	echo " $0 WORKFLOW1_FOLDER [WORKFLOW2_FOLDER] ...."
	echo
	echo "This program shows last 10 lines of pegasus-status for each workflow every $noOfWaitingSecondsDefault seconds."
	echo
	exit
fi

noOfWaitingSeconds=$noOfWaitingSecondsDefault

workflowFolders=$*

while [ "1" = "1" ]; do
	echo -n "=========current date is  "
	date
	for i in $workflowFolders; do
		echo $i
		pegasus-status $i|tail -n 10
		sleep $noOfWaitingSeconds
	done
done
