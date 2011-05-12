"""
This is a configuration file used by various scripts, to set up environment specific paths, etc.

In particular used for various GWAS..
"""

import os
#These paths need to be specified for each user!

user = os.getenv("USER")
home_dir = os.getenv("HOME") + "/"
config_file = home_dir + '.gwa_config'
results_dir = '/home/GMI/bjarni.vilhjalmsson/gwa_results/'
data_dir = '/home/GMI/bjarni.vilhjalmsson/data/'
data_1001_dir = '/home/GMI/bjarni.vilhjalmsson/Projects/data/1001genomes/'
data_quan_dir = '/home/GMI/bjarni.vilhjalmsson/Projects/data/quan_seq_data/'
rf_dir = None
env = dict()
env['home_dir'] = home_dir
env['results_dir'] = results_dir
env['data_dir'] = data_dir
env['data_1001_dir'] = data_1001_dir
env['data_quan_dir'] = data_quan_dir
env['rf_dir'] = rf_dir
env['default_lookup_db'] = 'gmi-ara-devel-be.gmi.oeaw.ac.at'
env['default_insert_db'] = 'gmi-ara-devel-be.gmi.oeaw.ac.at'
env['db_results_dir'] = '/Network/Data/250k/db/results/type_1/'
env['tmp_dir'] = '/tmp/'
env['phen_dir'] = '/Users/bjarni.vilhjalmsson/Projects/Data/phenotypes/'


#This should be changed in the .gwa_config file.
env['db_user'] = "bvilhjal"
env['db_passwd'] = "*rri_bjarni@usc"

try:
	import configobj
	try:
		config = configobj.ConfigObj(config_file)
		print 'GWAs configuration file loaded:', config_file
		for c in config:
			env[c] = config[c]
	except Exception, err_str:
		print 'Configuration file:', config_file, 'is missing?:', err_str
		print 'Warning! Using default configurations!'
#		print 'Creating',config_file,'with default configurations.'	
#		print 'Please update this file appropriately!  '
#		config = configobj.ConfigObj()
#		for c in env:
#			config[c]=env[c]
#		config.filename=config_file
#		config.write()
except Exception, err_str:
	print 'Failed importing the configobj module:', err_str
	print 'Warning! Using default configurations!'

