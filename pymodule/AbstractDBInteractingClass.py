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

class AbstractDBInteractingClass(object):
	__doc__ = __doc__
	option_default_dict = {
						('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('port', 0, ):[None, '', 1, 'database port number. must be non-empty if need ssh tunnel'],\
						('sshTunnelCredential', 0, ): ['', '', 1, 'a ssh credential to allow machine to access db server. \
										polyacti@login3, yuhuang@hpc-login2. if empty or port is empty, no tunnel', ],\
						('logFilename', 0, ): [None, '', 1, 'file to contain logs. use it only if this program is at the end of pegasus workflow \
		and has no output file'],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']
						}

	def __init__(self, inputFnameLs=None, **keywords):
		"""
		2011-7-11
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		self.inputFnameLs = inputFnameLs
		self.connectDB()
		
		#2012.7.4 keep track of all the source&destination files, used by moveNewISQFileIntoDBStorage()
		self.srcFilenameLs = []
		self.dstFilenameLs = []
	
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
	
	def outputLogMessage(self, logMessage=None, logFilename=None):
		"""
		2012.7.17
		"""
		if logFilename is None:
			logFilename = self.logFilename
		if logFilename:
			outf = open(logFilename, 'w')
			outf.write(logMessage)
			outf.close()
	
	def run(self):
		pass
	
if __name__ == '__main__':
	main_class = AbstractDBInteractingClass
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()