#!/bin/csh
#PBS -q cmb
#PBS -l walltime=24:00:00
#PBS -l mem=3000m 
#PBS -N TLS_400_test
(python /home/GMI/bjarni.vilhjalmsson/Projects/py_src/ThreeLocusTest.py -o /home/GMI/bjarni.vilhjalmsson/Projects/gwa_results/TLS_test_400_100_pm4  --numberPerRun=100  --latent_corr=1  --phenotypeFile=/tmp/test.phen  --phenotype_model=4  --pvalueThreshold=1.010000 --mapping_method=kw  --phenotype_error=0.500000  --kinship_error=0.000000 -t 54 400 > /home/GMI/bjarni.vilhjalmsson/Projects/gwa_results/TLS_test_400_100_pm4_job.out) >& /home/GMI/bjarni.vilhjalmsson/Projects/gwa_results/TLS_test_400_100_pm4_job.err
