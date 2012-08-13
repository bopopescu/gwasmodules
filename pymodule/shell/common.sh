#!/bin/bash

findValueGivenAnOptionName () {
	if [ -z "$1" ]
	then
		echo "Option Name is not provided."
		echo ;
	else
		optionNamePosition=`echo $arguments|awk -F ' ' '{i=1; while (i<=NF){if ($i=="'$1'") {print i}; ++i}}'`
		#echo $1 position: $optionNamePosition
		if [ -z "$optionNamePosition" ]
		then
			echo;
		else
			optionValuePosition=`echo $optionNamePosition+1|bc`
			#echo optionValuePosition $optionValuePosition
			optionValue=`echo $arguments|awk -F ' ' '{ if ('$optionValuePosition'<=NF) {print $'$optionValuePosition'} else print }'`
			echo $optionValue
		fi
	fi
}

checkVCFFileIfEmpty () {
	fname=$1
	if test -r $fname
	then
		suffix=`echo $fname|awk -F '.' '{print $NF}'`
		if test "$suffix" = "gz"
		then
			numberOfLoci=`gunzip -c $fname|grep -v "^#"|wc -l|awk -F ' ' '{print $1}'`
		else
			numberOfLoci=`grep -v "^#" $fname|wc -l|awk -F ' ' '{print $1}'`
		fi
		if test $numberOfLoci -gt 0
		then
			echo 0;
		else
			echo 1;
		fi
	else
		echo 1;
	fi
}



outputEmptyVCFWithInputHeader () {
	if [ -n "$vcfInputFname" ]
	then
		egrep "^#" $vcfInputFname 1>$outputVCFFname
	else
		if [ -n "$gzvcfInputFname" ]
		then
			gunzip -c $gzvcfInputFname | egrep "^#" 1>$outputVCFFname
		else
			touch $outputVCFFname
		fi
	fi
	#this echo is to avoid non-zero exit by egrep (nothing matches)
	echo "empty vcf with header from $vcfInputFname is created."
}

checkIfFileExists () {
	fname=$1
	if test -r $fname
	then
		echo 0;
	else
		echo 1;
	fi
}
exitIfFileExists () {
	fname=$1
	if test -r $fname
	then
		echo "$fname exists. Do not want to append/overwrite it. Check and rename it."
		echo 0;
		exit 2;
	else
		echo 1;
	fi
}
