#!/bin/bash
# 2013.2.4

if test $# -lt 3
then
	echo "Usage: $0 psmcFolderPath psmcOutputFname msCommandFname"
	echo
	echo "Notes:"
	echo "	#. psmcOutputFname is output of psmc (Li Durbin, 2011)."
	echo "	#. msCommandFname is for ms-simulation"
	echo
	echo "Example:"
	echo "	$0 ~/script/psmc/ diploid.psmc ms-cmd.sh"
exit
fi

source $HOME/.bash_profile

shellDir=~/script/shell/
source $shellDir/common.sh
psmcFolderPath=$1
psmcOutputFname=$2
msCommandFname=$3

psmc2historyPath=$psmcFolderPath/utils/psmc2history.pl
history2msPath=$psmcFolderPath/utils/history2ms.pl

#
#Usage: psmc2history.pl [-n 20] [-u "-1"] <in.psmc.par>
#default: (n=>20, u=>-1). n is probably the which time interval of psmc output
# u is probably the mutation rate or inverse mutation rate. if it's negative (default), estimate from psmc output.
#
#Usage:   history2ms.pl [options] <in.psmc.par>
#
#Options: -n INT    number of chromosome to simulate [2]
#         -L INT    length of each chromosome [30000000]
#         -s INT    skip used in psmc run [100]
#         -u FLOAT  neutral mutation rate [2.5e-08]
#         -R FLOAT  recomb. rate in hotspots are FLOAT times larger [10]
#         -g INT    years per generation [25]
#         -d INT    divergence time [0]
#         -r INT    # replicates [1]
#         -M        output macs command line
#
#

$psmc2historyPath $psmcOutputFname | $history2msPath > $msCommandFname


exitCodeAll="${PIPESTATUS[0]} ${PIPESTATUS[1]}"
exitCode1=`echo $exitCodeAll|awk -F ' ' '{print $1}'`
exitCode2=`echo $exitCodeAll|awk -F ' ' '{print $2}'`

echo exit codes: $exitCode1, $exitCode2

if test "$exitCode1" = "0" && test "$exitCode2" = "0"
then
	exit 0
else
	exit 3
fi
