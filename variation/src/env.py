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
env['tmp_dir']='/tmp/'

try:
	import configobj
	try:
		config = configobj.ConfigObj(config_file)
		print 'GWAs configuration file loaded:',config_file
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

