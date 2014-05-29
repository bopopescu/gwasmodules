#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s 

Description:
	2012.6.5
		abstract class for db-interacting pegasus jobs

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule import AbstractMapper
from pymodule import AbstractDBInteractingClass

class AbstractDBInteractingJob(AbstractDBInteractingClass):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	#option_default_dict.pop(('inputFname', 1, ))
	option_default_dict.update(AbstractMapper.db_option_dict.copy())
	option_default_dict.update({
							('logFilename', 0, ): [None, '', 1, 'file to contain logs. use it only when this job is at the end of pegasus workflow \
					and has no other output file. Because otherwise, pegasus optimization will get rid of this job (no output file, why need it?)'],\
							('commit', 0, int):[0, '', 0, 'commit db transaction'],\
							})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractDBInteractingClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
	
	def connectDB(self):
		"""
			split out of __init__() so that derived classes could overwrite this function
		"""
		pass
	
	