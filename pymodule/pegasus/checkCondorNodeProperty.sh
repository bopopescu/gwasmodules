#!/bin/sh

propertyDefault=sshDBTunnel
property=$1
if [ -z $property ]
then
	property=$propertyDefault
fi

#watch -d -n 12 "condor_status |tail -n $noOfLines"
#while [ 0==0 ];do
	fname1=/tmp/condorMachine
	fname2=/tmp/condorSSHDBTunnel
	condor_status -long|grep "^Machine " > $fname1
	condor_status -long|grep "^$property" > $fname2
	paste $fname1 $fname2|less
	#sleep 5
#done
