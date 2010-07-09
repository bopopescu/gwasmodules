"""
This is a configuration file used by various scripts, to set up environment specific paths, etc.

In particular used for various GWAS..
"""

import os
#These paths need to be specified for each user!

user = os.getenv("USER")
home_dir = os.getenv("HOME")+"/"
results_dir = home_dir+"tmp/"
rf_dir=None

if user=="bjarni":   # bamboo.usc.edu
	java_dir = "/System/Library/Frameworks/JavaVM.framework/Versions/1.6/Commands/"
	margarita_dir =  home_dir+"Projects/"
	#programDir="/home/cmb-01/bvilhjal/Projects/Python-snps/"

	
elif user =="bjarnivilhjalmsson": #betula (laptop)
	java_dir = "/System/Library/Frameworks/JavaVM.framework/Versions/1.6.0/Home/bin/"
	margarita_dir =  home_dir+"Projects/gwas_programs/"
	script_dir=home_dir+"Projects/src_2009/py_src/" 	#The dir with all the python scripts.
	

# Paths needed to run the scripts on the cluster.  
elif user=="bvilhjal":   #hpc-cmb.usc.edu  
	java_dir = "/usr/bin/"
	margarita_dir =  home_dir+"Projects/"
	rf_dir="/home/cmb-01/bvilhjal/Projects/randomForest/"
	home_dir = "/home/cmbpanfs-01/bvilhjal/"
	results_dir=home_dir+"results/"  	#The dir where all the results and output files are placed.
	script_dir=home_dir+"src/" 	#The dir with all the python scripts.
	data_dir=home_dir+"data/"		#The dir where both seq and phen data is located.


elif user=="bjarni.vilhjalmsson":   #hpc-cmb.usc.edu  
	results_dir=home_dir+"Projects/gwa_results/"  	#The dir where all the results and output files are placed.
	script_dir=home_dir+"Projects/py_src/" 	#The dir with all the python scripts.
	data_dir=home_dir+"Projects/data/"		#The dir where both seq and phen data is located.


elif user=="arthur.korte":   #hpc-cmb.usc.edu  
	script_dir="/home/GMI/bjarni.vilhjalmsson/Projects/py_src/" 	#The dir with all the python scripts.
	results_dir=home_dir+"gwas_results/"  	#The dir where all the results and output files are placed.
	data_dir=home_dir+"gwas_data/"		#The dir where both seq and phen data is located.

