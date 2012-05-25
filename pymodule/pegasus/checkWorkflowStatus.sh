#!/bin/sh

noOfLinesDefault=30
workflowDir=$1
noOfLines=$2
if [ -z $noOfLines ]
then
	noOfLines=$noOfLinesDefault
fi

watch -d -n 8 "pegasus-status --noutf8  $workflowDir |tail -n $noOfLines"
