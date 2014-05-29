#!/bin/bash
# 2013.3.6 bwa mem, piped to gzipped sam for paired end short-read
if test $# -lt 4
then
	echo "Usage: $0 bwaPath refFastaFname outputSamGzipFname fastqF1 [fastqF2] [otherMemArguments]"
	echo
	echo "Note:"
	echo "	1. refFastaFname must have been indexed by bwa."
	echo
	echo "Example:"
	echo "	$0 ~/bin/bwa 120_480K_supercontigs.fasta 435_6136_2004030_GA_5/gerald_81GPBABXX_8_CGATGT.sam.gz gerald_81GPBABXX_8_CGATGT_1.fastq.gz gerald_81GPBABXX_8_CGATGT_2.fastq.gz"
	echo
	echo "	$0 ~/bin/bwa 120_480K_supercontigs.fasta 435_6136_2004030_GA_5/gerald_81GPBABXX_8_CGATGT.sam.gz gerald_81GPBABXX_8_CGATGT_1.fastq.gz gerald_81GPBABXX_8_CGATGT_2.fastq.gz -M -a -t 4"
	echo
	echo "# set the fastqF2 as empty space, single-end read alignment with additional arguments"
	echo "	$0 ~/bin/bwa 120_480K_supercontigs.fasta 435_6136_2004030_GA_5/gerald_81GPBABXX_8_CGATGT.sam.gz gerald_81GPBABXX_8_CGATGT_1.fastq.gz ' ' -M -a -t 4"
	echo
exit
fi

bwaPath=$1
refFastaFname=$2
outputSamGzipFname=$3
fastqF1=$4
fastqF2=$5
shift
shift
shift
shift
shift
otherMemArguments=$*
$bwaPath mem $otherMemArguments $refFastaFname $fastqF1 $fastqF2 | gzip > $outputSamGzipFname
exitCodeAll="${PIPESTATUS[0]} ${PIPESTATUS[1]}"
exitCode=`echo $exitCodeAll|awk -F ' ' '{print $1}'`
exitCode2=`echo $exitCodeAll|awk -F ' ' '{print $2}'`
echo exit code: $exitCode, $exitCode2
if test "$exitCode" = "0" && test "$exitCode2" = "0"
then
	#normal exit
	exit
else
	if test "$exitCode" = "0"
	then
		finalExitCode=$exitCode2
	else
		finalExitCode=$exitCode
	fi
	exit $finalExitCode
fi
