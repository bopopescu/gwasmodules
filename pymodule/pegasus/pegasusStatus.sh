#!/bin/sh

noOfLinesDefault=30
noOfLines=$1
if [ -z $noOfLines ]
then
	noOfLines=$noOfLinesDefault
fi

watch -d -n 6 "pegasus-status |tail -n $noOfLines"
