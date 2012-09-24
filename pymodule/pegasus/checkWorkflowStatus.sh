#!/bin/sh

noOfLinesDefault=30
workflowDir=$1
noOfLines=$2
if [ -z $noOfLines ]
then
	noOfLines=$noOfLinesDefault
fi

#--noutf8 
watch -d -n 13 "pegasus-status  $workflowDir |tail -n $noOfLines"
