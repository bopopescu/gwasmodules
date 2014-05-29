#!/bin/bash

runningWorkflowLogFnameDefault=runningWorkflows.txt
failedWorkflowLogFnameDefault=failedWorkflows.txt
TOPDIR=`pwd`

if test $# -lt 1 ; then
	echo "Usage:"
	echo "  $0 workflowFolder [runningWorkflowLogFname] [failedWorkflowLogFname]"
	echo ""
	echo "Note:"
	echo "	#. runningWorkflowLogFname will be opened in append mode and will contain the folder names of successfully submitted workflows. One on each line. Default is $runningWorkflowLogFnameDefault."
	echo "	#. failedWorkflowLogFname will be opened in append mode and will contain the folder names of UN-successfully submitted workflows. One on each line. Default is $failedWorkflowLogFnameDefault."
	echo ""
	echo "Examples:"
	echo "	#"
	echo "	$0 work/BaseQualityRecalibration/LocalRealignmentBQSR_AlnID2828_2847_vsMethod87.2013.Apr.17T230821/"
	echo
	echo "	# use non-default files to record running and failed workflows"
	echo "	$0 work/BaseQualityRecalibration/LocalRealignmentBQSR_AlnID2828_2847_vsMethod87.2013.Apr.17T230821/ runningBQSRWorkflows.txt failedBQSRWorkflows.txt"
	echo
	exit 1
fi
workdir=$1
runningWorkflowLogFname=$2
failedWorkflowLogFname=$3

if test -z "$runningWorkflowLogFname"
then
	runningWorkflowLogFname=$runningWorkflowLogFnameDefault
fi

if test -z "$failedWorkflowLogFname"
then
	failedWorkflowLogFname=$failedWorkflowLogFnameDefault
fi

#pegasus-run --conf $workdir/pegasus*.properties $workdir
pegasus-run $workdir
exitCode=$?
if test $exitCode = "0"; then
	echo $workdir >> $runningWorkflowLogFname
elif test $exitCode = "2"; then
	#error due to existence of .dag.condor.sub; .dag.lib.out; .dag.lib.out; .dag.dagman.log files (the workflow was terminated by external force, and thus no proper folder cleanup)
	pushd $workdir
	for suffix in .dag.condor.sub .dag.lib.err .dag.lib.out .dag.dagman.log; do
		fname=`ls *$suffix`;
		if test -w $fname; then
			mv $fname $fname.old
		fi
	done
	popd
	#resubmit it again
	$0 $workdir $runningWorkflowLogFname $failedWorkflowLogFname
else
	echo $workdir >>$failedWorkflowLogFname
fi
