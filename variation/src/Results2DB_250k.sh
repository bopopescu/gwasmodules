#!/bin/bash

if test $# -lt 6
then
	echo "Usage:"
	echo "    Results2DB_250k.sh FILE_PREFIX CALL_METHOD_ID ANALYSIS_METHOD_ID DB_PASSWORD START_PHENOTYPE_ID END_PHENOTYPE_ID [DB_USER] [DB_HOST] [RESULTS_METHOD_TYPE_ID]"
	echo
	echo "This script calls Results2DB_250k.py to submit results into database."
	echo	"DB_USER is optional. If not supplied ,it's assumed to be yh."
	echo	"DB_HOST is optional. If not supplied ,it's assumed to be papaya.usc.edu."
	echo	"RESULTS_METHOD_TYPE_ID is optional. If not supplied ,it's assumed to be 1."
	echo
	echo "Example:"
	echo " Results2DB_250k.sh /Network/Data/250k/tmp-bvilhjal/marg_results/Marg_newDataset_ 17 5 secret 1 180 yh banyan.usc.edu 1"
	exit
fi

file_prefix=$1
call_method_id=$2
analysis_method_id=$3
db_password=$4
start_phenotype_id=$5
end_phenotype_id=$6
if test -n "$7" ; then
	db_user=$7
else
	db_user=yh
fi
if test -n "$8" ; then
	db_host=$8
else
	db_host=papaya.usc.edu
fi
if test -n "$9" ; then
	results_method_type_id=$9
else
	results_method_type_id=1
fi

for((i=$start_phenotype_id;i<=$end_phenotype_id;i++))
	do echo $i;
	#~/script/variation/src/Results2DB_250k.py -a $call_method_id -e $i -i $file_prefix$i\_*[^s][^r].tsv -l $analysis_method_id -u $db_user -p $db_password -c -z $db_host -s $results_method_type_id
	~/script/variation/src/Results2DB_250k.py -a $call_method_id -e $i -i $file_prefix$i\_*.tsv -l $analysis_method_id -u $db_user -p $db_password -c -z $db_host -s $results_method_type_id
done

