#!/bin/csh
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=8 
#PBS -l mem=5000mb
#PBS -q cmb
#PBS -N p1_emma_test_emma
(python /home/cmbpanfs-01/bvilhjal/src/gwa.py --id=/home/cmbpanfs-01/bvilhjal/results/test_pid1_LD_emma  --region_plots=0   --specific_methods=emma --debug_filter=0.005  ../Data/250k/250K_t54.csv ../Data/phenotypes/phen_all_051710.tsv 1 > /home/cmbpanfs-01/bvilhjal/results/test_pid1_LD_emma_job.out) >& /home/cmbpanfs-01/bvilhjal/results/test_pid1_LD_emma_job.err
