#!/bin/bash
# 2013.2.4

noOfChromosomesDefault=2
chromosomeLengthDefault=30000000
mutationRatePerNucleotidePerGenerationDefault=5e-9
noOfYearsPerGenerationDefault=25
recombinationHotSpotRateMultiplerDefault=10
divergenceTimeDefault=0


if test $# -lt 3
then
	echo "Usage: $0 psmcFolderPath psmcOutputFname msCommandFname [history2msArguments]"
	echo
	echo "Notes:"
	echo "	#. This program generates ms commandline that could simulate sequences subject to the PSMC demographic output."
	echo "	#. psmcOutputFname is output of psmc (Li Durbin, 2011)."
	echo "	#. msCommandFname is for ms-simulation"
	echo "	#. history2msArguments is optional used to change the defaults in history2ms.pl. Here are potential arguments."
	echo "		Usage:   history2ms.pl [options] <in.psmc.par>"
	echo "		"
	echo "		Options: -n INT    number of chromosome to simulate [2]"
	echo "		         -L INT    length of each chromosome [30000000]"
	echo "		         -s INT    skip used in psmc run [100]"
	echo "		         -u FLOAT  neutral mutation rate [2.5e-08]"
	echo "		         -R FLOAT  recomb. rate in hotspots are FLOAT times larger [10]"
	echo "		         -g INT    years per generation [25]"
	echo "		         -d INT    divergence time [0]"
	echo "		         -r INT    # replicates [1]"
	echo "		         -M        output macs command line"
	echo
	echo "Example:"
	echo "	$0 ~/script/psmc/ diploid.psmc ms-cmd.sh"
	echo
	echo "	# simulate 2 chromosomes (-n 2), 20Mb (-L 20000000)"
	echo "	$0 ~/script/psmc/ diploid.psmc ms-cmd.sh -n 2 -L 20000000"
	echo
	
exit
fi

source $HOME/.bash_profile

shellDir=~/script/shell/
source $shellDir/common.sh
psmcFolderPath=$1
psmcOutputFname=$2
msCommandFname=$3

shift
shift
shift
#after shifting out all the previous arguments times, all arguments left are passed to history2ms.pl
history2msArguments=$*

echo "Arguments passed to history2ms.pl is $history2msArguments."

psmc2historyPath=$psmcFolderPath/utils/psmc2history.pl
history2msPath=$psmcFolderPath/utils/history2ms.pl

#
#Usage: psmc2history.pl [-n 20] [-u "-1"] <in.psmc.par>
#default: (n=>20, u=>-1). n is probably the which time interval of psmc output
# u is probably the mutation rate or inverse mutation rate. if it's negative (default), estimate from psmc output.
#
#
#


$psmc2historyPath $psmcOutputFname | $history2msPath $history2msArguments > $msCommandFname


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
