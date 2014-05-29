#!/bin/sh

shellDir=`dirname $0`
if test $# -lt 3 ; then
	echo "  $0 InfoFieldKey outputFname inputFile1 inputFile2 ..."
	echo ""
	echo "Note:"
	echo " 	#. This program calls pegasus/mapper/extractor/ExtractInfoFromVCF.py to extract VCF info fields, bgzip, tabix index and creates a database for searching."
	echo "	#. InfoFieldKey is one of keys in VCF info field whose value is to be extracted, like AF, HaplotypeScore, DP."
	echo " 	#. outputFname is the un-gzipped version of the final database."
	echo ""
	echo "Examples:"
	echo " 	#extract AF from VCF files. The final databse is tmp/method_232_AF/VCF_AF.tsv.gz with a .tbi index file."
	echo " 	$0 AF tmp/method_232_AF/VCF_AF.tsv ~/NetworkData/vervet/db/genotype_file/method_232/*vcf.gz"
	echo ""
	exit 1
fi


InfoFieldKey=$1
outputFname=$2
shift
shift
inputFiles=$*

~/script/pymodule/pegasus/mapper/extractor/ExtractInfoFromVCF.py -o $outputFname -k $InfoFieldKey $inputFiles
echo -n "bgzipping $outputFname ..."
bgzip $outputFname
echo ""
echo -n "tabix indexing $outputFname.gz ..."
tabix -s 1  2 -e 2 -S 1 -f $outputFname.gz
echo ""

