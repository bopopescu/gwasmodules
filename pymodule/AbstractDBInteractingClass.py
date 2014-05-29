#!/usr/bin/env python
"""
Examples:

Description:
	2012.2.10
		an abstract class that includes db-related arguments and etc
"""

import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from ProcessOptions import ProcessOptions
import utils
from utils import PassingData
from pegasus.mapper.AbstractMapper import AbstractMapper

class AbstractDBInteractingClass(AbstractMapper):
	__doc__ = __doc__
	
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict = {
						('sshTunnelCredential', 0, ): ['', '', 1, 'a ssh credential to allow machine to access db server. \
										polyacti@login3, yuhuang@hpc-login2. if empty or port is empty, no tunnel', ],\
						('logFilename', 0, ): [None, '', 1, 'file to contain logs. use it only if this program is at the end of pegasus workflow \
		and has no output file'],\
						}
	option_default_dict.update(AbstractMapper.db_option_dict)
	
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		2011-7-11
		"""
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		
		#2012.7.4 keep track of all the source&destination files, used by moveNewISQFileIntoDBStorage()
		self.srcFilenameLs = []
		self.dstFilenameLs = []
		#2012.11.18
		if getattr(self, "logFilename", None):
			self.logF = open(self.logFilename, 'w')
		else:
			self.logF = None
	
	def parseInputFile(self, inputFname=None, **keywords):
		"""
		2013.08.23
			if a program is adding a file to db-affiliated storage, this is used for parsing.
		"""
		pass
		return PassingData()
	
	def connectDB(self):
		"""
		2012.5.11
			place holder. AbstractVervetMapper.py will use it 
		"""
		pass
	
	def rmGivenFiles(self, filenameLs=[], rmCommand='rm -rf'):
		"""
		2012.7.4
			delete all files in filenameLs
		"""
		sys.stderr.write("Deleting %s files ...\n"%(len(filenameLs)))
		for filename in filenameLs:
			commandline = '%s %s'%(rmCommand, filename)
			return_data = utils.runLocalCommand(commandline, report_stderr=True, report_stdout=True)
			if return_data.stderr_content:
				#something wrong.
				sys.stderr.write("commandline %s failed: %s\n"%(commandline, return_data.stderr_content))
		sys.stderr.write(".\n")
	
	def cleanUpAndExitOnFailure(self, exitCode=1):
		"""
		2012.7.13 an exit function when the program failed somewhere
		"""
		#delete all target files.
		self.rmGivenFiles(filenameLs=self.dstFilenameLs)
		sys.exit(exitCode)
	
	def cleanUpAndExitOnSuccess(self, exitCode=0):
		"""
		2012.7.13  an exit function when the program succeeded in the end
		"""
		sys.exit(exitCode)
	
	def outputLogMessage(self, logMessage=None, logFilename=None, logF=None):
		"""
		2012.11.18
			1. do not close _logF
			2. use self.logF if None available.
			3. append to logFilename if the latter is given.
		2012.7.17
		"""
		_logF = None
		if logF is None:
			if logFilename:
				_logF = open(logFilename, 'a')
			elif self.logF:
				_logF = self.logF
		if _logF:
			_logF.write(logMessage)
	
	def closeLogF(self):
		"""
		2012.11.18
		"""
		if getattr(self, 'logF', None):
			self.logF.close()
	
	def sessionRollback(self, session=None):
		"""
		2013.04.05
			wrap try...except around it because sometimes db connection is gone halfway through
			and it causes program to exit without proper cleanup
		"""
		try:
			session.rollback()
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()

	def run(self):
		pass
	
if __name__ == '__main__':
	main_class = AbstractDBInteractingClass
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()