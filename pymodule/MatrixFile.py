#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s

Description:
	2012.3.1

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from ProcessOptions import  ProcessOptions
import csv
from pymodule import utils, figureOutDelimiter

class MatrixFile(object):
	__doc__ = __doc__
	option_default_dict = {
						('inputFname', 0, ): [None, 'i', 1, 'filename.'],\
						('inputFile', 0, ): [None, '', 1, 'opened file handler'],\
						('openMode', 1, ): ['r', '', 1, 'r'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']
						}
	def __init__(self, inputFname=None, **keywords):
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		if not self.inputFname:
			self.inputFname = inputFname
		self.isCSVReader = True
		if self.inputFname:
			self.inputFile = utils.openGzipFile(self.inputFname, openMode=self.openMode)
		
		delimiter = figureOutDelimiter(self.inputFile)
		if delimiter=='\t' or delimiter==',':
			reader = csv.reader(self.inputFile, delimiter=delimiter)
		else:
			reader = self.inputFile
			self.isCSVReader = False
		self.reader = reader
		self.col_name2index = None
		self.header = None
	
	def constructColName2IndexFromHeader(self):
		"""
		2012.8.23
		"""
		self.header = self.next()
		self.col_name2index = utils.getColName2IndexFromHeader(self.header)
	
	def getColIndexGivenColHeader(self, colHeader=None):
		"""
		2012.8.23
		
		"""
		if self.col_name2index is None:
			return None	#no header
		else:
			return self.col_name2index.get(colHeader)
	
	def __iter__(self):
		return self
	
	def next(self):
		try:
			row = self.reader.next()
		except:
			raise StopIteration
		if not self.isCSVReader:
			row = row.strip().split()
		return row
	
	def close(self):
		del self.reader
		self.inputFile.close()
	

	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		

if __name__ == '__main__':
	main_class = MatrixFile
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()