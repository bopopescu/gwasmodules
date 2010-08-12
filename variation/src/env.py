"""
This is a configuration file used by various scripts, to set up environment specific paths, etc.

In particular used for various GWAS..
"""

import os
#These paths need to be specified for each user!

user = os.getenv("USER")
home_dir = os.getenv("HOME")+"/"
config_file = home_dir+'.gwa_config'
results_dir ='/tmp/'
data_dir = '/Network/Data/250k/'
rf_dir=None
env = dict()
env['home_dir']=home_dir
env['results_dir']=results_dir
env['data_dir']=data_dir
env['rf_dir']=rf_dir
env['default_lookup_db']='gmi-ara-devel-be.gmi.oeaw.ac.at'
env['default_insert_db']='gmi-ara-devel-be.gmi.oeaw.ac.at'
env['db_results_dir']='/Network/Data/250k/db/results/type_1/'


try:
	import configobj
	try:
		config = configobj.ConfigObj(config_file)
		print config
		for c in config:
			env[c] = config[c]
	except Exception, err_str:
		print 'Configuration file:',config_file,'is missing?:',err_str
		print 'Warning! Using default configurations!'
#		print 'Creating',config_file,'with default configurations.'	
#		print 'Please update this file appropriately!  '
#		config = configobj.ConfigObj()
#		for c in env:
#			config[c]=env[c]
#		config.filename=config_file
#		config.write()
except Exception, err_str:
	print 'Failed importing the configobj module:',err_str	
	print 'Warning! Using default configurations!'


#
#
#if user=="bjarni":   # bamboo.usc.edu
#	java_dir = "/System/Library/Frameworks/JavaVM.framework/Versions/1.6/Commands/"
#	margarita_dir =  home_dir+"Projects/"
#	#programDir="/home/cmb-01/bvilhjal/Projects/Python-snps/"
#
#	
#elif user =="bjarnivilhjalmsson": #betula (laptop)
#	java_dir = "/System/Library/Frameworks/JavaVM.framework/Versions/1.6.0/Home/bin/"
#	margarita_dir =  home_dir+"Projects/gwas_programs/"
#	script_dir=home_dir+"Projects/src_2009/py_src/" 	#The dir with all the python scripts.
#	
#
## Paths needed to run the scripts on the cluster.  
#elif user=="bvilhjal":   #hpc-cmb.usc.edu  
#	java_dir = "/usr/bin/"
#	margarita_dir =  home_dir+"Projects/"
#	rf_dir="/home/cmb-01/bvilhjal/Projects/randomForest/"
#	home_dir = "/home/cmbpanfs-01/bvilhjal/"
#	results_dir=home_dir+"results/"  	#The dir where all the results and output files are placed.
#	script_dir=home_dir+"src/" 	#The dir with all the python scripts.
#	data_dir=home_dir+"data/"		#The dir where both seq and phen data is located.
#
#
#elif user=="bjarni.vilhjalmsson":   #hpc-cmb.usc.edu  
#	results_dir=home_dir+"Projects/gwa_results/"  	#The dir where all the results and output files are placed.
#	script_dir=home_dir+"Projects/py_src/" 	#The dir with all the python scripts.
#	data_dir=home_dir+"Projects/data/"	#The dir where both seq and phen data is located.
#
#
#elif user=="arthur.korte":   #hpc-cmb.usc.edu  
#	script_dir="/home/GMI/bjarni.vilhjalmsson/Projects/py_src/" 	#The dir with all the python scripts.
#	results_dir=home_dir+"gwas_results/"  	#The dir where all the results and output files are placed.
#	data_dir=home_dir+"gwas_data/"		#The dir where both seq and phen data is located.
