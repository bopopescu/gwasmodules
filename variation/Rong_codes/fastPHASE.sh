#!/bin/sh
input_dir=$1
cd $input_dir
date
for i in `ls *chr?`; do echo $i;  echo ~/script/variation/genotyping/fastphase.1.2.3.linux/fastPHASE -o$i.out $i;
~/script/variation/genotyping/fastphase.1.2.3.linux/fastPHASE -o$i.out $i; done
date
