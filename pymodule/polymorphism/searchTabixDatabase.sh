#!/bin/sh

shellDir=`dirname $0`
chromosomeColumnInQueryDefault=1
startColumnInQueryDefault=2
if test $# -lt 2 ; then
	echo "  $0 inputQueryFname databaseFname [chromosomeColumnInQuery] [startColumnInQuery]"
	echo ""
	echo "Note:"
	echo " 	#. This program searches a tabix database given SNPs in inputQueryFname."
	echo "	#. inputQueryFname is a tab delimited file with chromosome and start positions."
	echo " 	#. databaseFname is a gzipped tabix-indexed database."
	echo " 	#. default chromosomeColumnInQuery is $chromosomeColumnInQueryDefault (column order starts from 1, not 0)."
	echo " 	#. default startColumnInQuery is $startColumnInQueryDefault"
	echo ""
	echo "Examples:"
	echo " 	$0 /tmp/damaging.short.tsv method_232_AF.tsv.gz > output.tsv"
	echo " 	$0 tmp/damaging.short.tsv 94506_VCF_94286_VCF_CAE8.vcf.gz >> tmp/damaging.short.vcf"
	echo ""
	exit 1
fi


inputQueryFname=$1
databaseFname=$2
chromosomeColumnInQuery=$chromosomeColumnInQueryDefault
startColumnInQuery=$startColumnInQueryDefault

echo "inputQueryFname=$inputQueryFname" 1>&2
echo "databaseFname=$databaseFname" 1>&2

for i in `awk -F '\t' '{print $'$chromosomeColumnInQuery'":"$'$startColumnInQuery'"-"$'$startColumnInQuery' }' $inputQueryFname`; do
	tabix $databaseFname $i
done
