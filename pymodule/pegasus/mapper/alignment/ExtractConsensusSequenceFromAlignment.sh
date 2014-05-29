#!/bin/bash
# 2013.2.4
minDPDefault=10
maxDPDefault=500
minMapQDefault=30
minBaseQDefault=20
minRMSMapQDefault=10
minDistanceToIndelDefault=5

if test $# -lt 3
then
	echo "Usage: $0 refFastaFname inputAlignmentFname outputFname [minDP] [maxDP] [minMapQ] [minBaseQ] [minRMSMapQ] [minDistanceToIndel] [mpileupInterval]"
	echo
	echo "Note:"
	echo "	#. This is a program that generates a consensus sequence out of inputAlignmentFname for psmc (Li Durbin 2011) analysis."
	echo "	#. refFastaFname must be faidx-indexed."
	echo "	#. shell has to be bash. PIPESTATUS is used."
	echo "	#. default minDP is $minDPDefault."
	echo "	#. default maxDP is $maxDPDefault."
	echo "	#. default minMapQ is $minMapQDefault."
	echo "	#. default minBaseQ is $minBaseQDefault."
	echo "	#. default minRMSMapQ is $minRMSMapQDefault."
	echo "	#. default minDistanceToIndel is $minDistanceToIndelDefault."
	echo "	#. default mpileupInterval is nothing, no restriction."
	echo
	echo "Example:"
	echo "	$0 individual_sequence/524_superContigsMinSize2000.fasta individual_alignment/555_15_vs_524_by_2.bam /tmp/555_15_vs_524_by_2.fq.gz 24 48"
exit
fi

source $HOME/.bash_profile

shellDir=~/script/shell/
source $shellDir/common.sh

bwaPath=$HOME/bin/bwa
samtoolsPath=$HOME/bin/samtools
bcftoolsPath=$HOME/bin/bcftools
vcfutilsPath=$HOME/bin/vcfutils.pl

picardPath=$HOME/script/picard/dist/

refFastaFname=$1
inputAlignmentFname=$2
outputFname=$3

minDP=$4
if [ -z $minDP ]
then
	minDP=$minDPDefault
fi

maxDP=$5
if [ -z $maxDP ]
then
	maxDP=$maxDPDefault
fi

minMapQ=$6
if [ -z $minMapQ ]
then
	minMapQ=$minMapQDefault
fi

minBaseQ=$7
if [ -z $minBaseQ ]
then
	minBaseQ=$minBaseQDefault
fi

minRMSMapQ=$8
if [ -z $minRMSMapQ ]
then
	minRMSMapQ=$minRMSMapQDefault
fi

minDistanceToIndel=$9
if [ -z $minDistanceToIndel ]
then
	minDistanceToIndel=$minDistanceToIndel
fi

shift	#no double-digit (or more digits) argument number
mpileupInterval=$9
if [ -n $mpileupInterval ]; then
	if test -r "$mpileupInterval"
	then
		#it's a file, BED format.
		echo "$mpileupInterval is treated as a BED file."
		mpileupIntervalArgument="-l $mpileupInterval"
	else
		echo "$mpileupInterval" is treated as a region-defining string.
		mpileupIntervalArgument="-r $mpileupInterval"
	fi
else
	mpileupIntervalArgument=""
fi

echo "mpileupIntervalArgument is $mpileupIntervalArgument."

###################
#
#Usage: samtools mpileup [options] in1.bam [in2.bam [...]]
#
#Input options:
#
#       -6           assume the quality is in the Illumina-1.3+ encoding
#       -A           count anomalous read pairs
#       -B           disable BAQ computation
#       -b FILE      list of input BAM files [null]
#       -C INT       parameter for adjusting mapQ; 0 to disable [0]
#       -d INT       max per-BAM depth to avoid excessive memory usage [250]
#       -f FILE      faidx indexed reference sequence file [null]
#       -G FILE      exclude read groups listed in FILE [null]
#       -l FILE      list of positions (chr pos) or regions (BED) [null]
#       -M INT       cap mapping quality at INT [60]
#       -r STR       region in which pileup is generated [null]
#       -R           ignore RG tags
#       -q INT       skip alignments with mapQ smaller than INT [0]
#       -Q INT       skip bases with baseQ/BAQ smaller than INT [13]
#
#Output options:
## 2012.5.4 "-S -D" from mpileup increases running time significantly
#
#       -D           output per-sample DP in BCF (require -g/-u)
#       -g           generate BCF output (genotype likelihoods)
#       -O           output base positions on reads (disabled by -g/-u)
#       -s           output mapping quality (disabled by -g/-u)
#       -S           output per-sample strand bias P-value in BCF (require -g/-u)
#       -u           generate uncompress BCF output
#
#SNP/INDEL genotype likelihoods options (effective with `-g' or `-u'):
#
#       -e INT       Phred-scaled gap extension seq error probability [20]
#       -F FLOAT     minimum fraction of gapped reads for candidates [0.002]
#       -h INT       coefficient for homopolymer errors [100]
#       -I           do not perform indel calling
#       -L INT       max per-sample depth for INDEL calling [250]
#       -m INT       minimum gapped reads for indel candidates [1]
#       -o INT       Phred-scaled gap open sequencing error probability [40]
#       -P STR       comma separated list of platforms for indels [all]
#
#Notes: Assuming diploid individuals.

##############
#Usage:   vcfutils.pl vcf2fq [options] <all-site.vcf>
#
#Options: -d INT    minimum depth          [3]
#         -D INT    maximum depth          [100000]
#         -Q INT    min RMS mapQ           [10]
#         -l INT    INDEL filtering window [5]
#
#
####################
#
#Usage: bcftools view [options] <in.bcf> [reg]
#
#Input/output options:
#
#       -A        keep all possible alternate alleles at variant sites
#       -b        output BCF instead of VCF
#       -D FILE   sequence dictionary for VCF->BCF conversion [null]
#       -F        PL generated by r921 or before (which generate old ordering)
#       -G        suppress all individual genotype information
#       -l FILE   list of sites (chr pos) or regions (BED) to output [all sites]
#       -L        calculate LD for adjacent sites
#       -N        skip sites where REF is not A/C/G/T
#       -Q        output the QCALL likelihood format
#       -s FILE   list of samples to use [all samples]
#       -S        input is VCF
#       -u        uncompressed BCF output (force -b)
#
#Consensus/variant calling options:
#
#       -c        SNP calling (force -e)
#       -d FLOAT  skip loci where less than FLOAT fraction of samples covered [0]
#       -e        likelihood based analyses
#       -g        call genotypes at variant sites (force -c)
#       -i FLOAT  indel-to-substitution ratio [-1]
#       -I        skip indels
#       -p FLOAT  variant if P(ref|D)<FLOAT [0.5]
#       -P STR    type of prior: full, cond2, flat [full]
#       -t FLOAT  scaled substitution mutation rate [0.001]
#       -T STR    constrained calling; STR can be: pair, trioauto, trioxd and trioxs (see manual) [null]
#       -v        output potential variant sites only (force -c)
#
#Contrast calling and association test options:
#
#       -1 INT    number of group-1 samples [0]
#       -C FLOAT  posterior constrast for LRT<FLOAT and P(ref|D)<0.5 [1]
#       -U INT    number of permutations for association testing (effective with -1) [0]
#       -X FLOAT  only perform permutations for P(chi^2)<FLOAT [0.01]
#

###### aruments for $vcfutilsPath varFilter
#Usage:   vcfutils.pl varFilter [options] <in.vcf>
#
#Options: -Q INT    minimum RMS mapping quality for SNPs [10]
#         -d INT    minimum read depth [2]
#         -D INT    maximum read depth [10000000]
#         -a INT    minimum number of alternate bases [2]
#         -w INT    SNP within INT bp around a gap to be filtered [3]
#         -W INT    window size for filtering adjacent gaps [10]
#         -1 FLOAT  min P-value for strand bias (given PV4) [0.0001]
#         -2 FLOAT  min P-value for baseQ bias [1e-100]
#         -3 FLOAT  min P-value for mapQ bias [0]
#         -4 FLOAT  min P-value for end distance bias [0.0001]
#		 -e FLOAT  min P-value for HWE (plus F<0) [0.0001]
#         -p        print filtered variants
#
#Note: Some of the filters rely on annotations generated by SAMtools/BCFtools.
#
#-D100: maximum depth=100
#-w 10:	SNP within 10 bp around a gap to be filtered.
#-d 3:	minimum read depth is 3.
#-e 0:	min P-value for HWE (plus F<0). setting it to zero removes this filter.


bcftoolsArguments=""

echo "bcftoolsArguments is $bcftoolsArguments."

#2013.05.08 trace what commands are being executed....
set -vx


echo "outputFname is $outputFname"
#2013.05.08 need "bash -c" if  "-" or "|" or ">" is put into a command "" string, but you will lose PIPESTATUS
#bash -c "$samtoolsCommandline | $bcftoolsCommandline" 

$samtoolsPath mpileup -C50 -q $minMapQ -Q $minBaseQ $mpileupIntervalArgument -u -f $refFastaFname $inputAlignmentFname | $bcftoolsPath view $bcftoolsArguments -c - | $vcfutilsPath vcf2fq -d $minDP -D $maxDP -Q $minRMSMapQ -l $minDistanceToIndel | gzip > $outputFname

exitCodeAll="${PIPESTATUS[0]} ${PIPESTATUS[1]} ${PIPESTATUS[2]} ${PIPESTATUS[3]}"
exitCode1=`echo $exitCodeAll|awk -F ' ' '{print $1}'`
exitCode2=`echo $exitCodeAll|awk -F ' ' '{print $2}'`
exitCode3=`echo $exitCodeAll|awk -F ' ' '{print $3}'`
exitCode4=`echo $exitCodeAll|awk -F ' ' '{print $4}'`

echo exit codes: $exitCode1, $exitCode2, $exitCode3, $exitCode4

if test "$exitCode1" = "0" && test "$exitCode2" = "0" && test "$exitCode3" = "0" && test "$exitCode4" = "0"
then
	exit 0
else
	exit 3
fi
