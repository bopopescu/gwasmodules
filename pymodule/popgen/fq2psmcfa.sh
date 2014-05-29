#!/bin/bash
# 2013.2.4
minDPDefault=10
maxDPDefault=500
minMapQDefault=30
minBaseQDefault=20

if test $# -lt 3
then
	echo "Usage: $0 psmcFolderPath consensusSequenceFname consensusSequencePSMCFname [minBaseQ]"
	echo
	echo "Notes:"
	echo "	#. consensusSequenceFname is output of pymodule/pegasus/mapper/ExtractConsensusSequenceFromAlignment.sh ."
	echo "	#. consensusSequencePSMCFname will be input of psmc (Li Durbin, 2011)."
	echo "	#. default minBaseQ is $minBaseQDefault."
	echo
	echo "Example:"
	echo "	$0 ~/script/psmc/ diploid.fq.gz diploid.psmcfa"
exit
fi

source $HOME/.bash_profile

shellDir=~/script/shell/
source $shellDir/common.sh
psmcFolderPath=$1
consensusSequenceFname=$2
consensusSequencePSMCFname=$3
minBaseQ=$4
if [ -z $minBaseQ ]
then
	minBaseQ=$minBaseQDefault
fi

#
# Usage: fq2psmcfa [-nvx] [-q 10] [-g 10000] [-s 100] <in.fq>
#
# case 'q': min_qual = atoi(optarg); break;
# case 'x': mask_pseudo = 1; break;
# case 'v': tv_only = 1; break;
# case 'n': ts_only = 1; break;
# case 'g': n_min_good = atoi(optarg); break;
# case 's': BLOCK_LEN = atoi(optarg); break;
# 

fq2psmcfaPath=$psmcFolderPath/utils/fq2psmcfa
echo "minBaseQ is $minBaseQ."

$fq2psmcfaPath -q $minBaseQ $consensusSequenceFname > $consensusSequencePSMCFname
