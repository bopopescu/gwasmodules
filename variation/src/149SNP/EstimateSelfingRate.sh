#!/bin/sh
input_dir='/tmp/'
cd $input_dir
pop_size=$1
filename_prefix=stock_149SNP_popid2ecotypeid_$pop_size
which_method=$2
output_dir='/tmp/est_selfing_out'
mkdir $output_dir
for input_file in `ls $filename_prefix.*`
	do echo $input_file, which_method=$which_method
	~/script/variation/src/EstimateSelfingRate.py -i $input_file -o $output_dir/$input_file\.y$which_method -y $which_method -s popid2s_$pop_size -c -u yh -p secret
done

