#!/bin/sh

realSAMtoolsPath=~/script/samtools.git/samtools

if test $# -lt 1
	then
	echo "2013.3.18 This temporary program is used to switch-&-bait the samtools used in workflows to convet sam to bam. The symlink ~/bin/samtools will be pointed to this program."
	echo " The samtools sam-2-bam syntax looks like: /u/home/eeskin/polyacti/bin/samtools  view -bSh -o  1675_5908_2000001_GA_0_1/93754_79836_1669_gerald_C0TKHACXX_3_ACTTGA_1_.bam   1675_5908_2000001_GA_0_1/93754_79836_1669_gerald_C0TKHACXX_3_ACTTGA_1_.sam.gz"
	echo " if 2nd argument is not view, then call the real samtools, $realSAMtoolsPath with all arguments."
	exit
fi
samtoolsCommand=$1
if test "$samtoolsCommand" = "view"; then
	input=$5
	output=$4
	java -jar ~/script/picard/dist/SamFormatConverter.jar VALIDATION_STRINGENCY=LENIENT I=$input O=$output
else
	$realSAMtoolsPath $*
fi
