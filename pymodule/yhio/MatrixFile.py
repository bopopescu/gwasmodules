#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s

Description:
	2012.3.1 a matrix file API like csv.reader, it can deal with any delimiter: coma, tab, space or multi-space.
	Example:
	
		reader = MatrixFile(inputFname='/tmp/input.txt', openMode='r')
		reader = MatrixFile('/tmp/input.txt', openMode='r')
		reader.constructColName2IndexFromHeader()
		for row in reader:
			row[reader.getColName2IndexFromHeader('KID')]
		
		inf = utils.openGzipFile(inputFname, openMode='r')
		reader = MatrixFile(inputFile=inf)
		
		#2013.2.1 writing
		writer = MatrixFile('/tmp/output.txt', openMode='w', delimiter='\t')
		writer.writeHeader(...)
		writer.writerow(row)
		writer.close()
	
	

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule.ProcessOptions import  ProcessOptions
import csv
from pymodule import utils, figureOutDelimiter

class MatrixFile(object):
	__doc__ = __doc__
	option_default_dict = {
						('inputFname', 0, ): [None, 'i', 1, 'filename.'],\
						('inputFile', 0, ): [None, '', 1, 'opened file handler'],\
						('openMode', 1, ): ['r', '', 1, 'mode to open the inputFname if inputFile is not presented.'],\
						('delimiter', 0, ): [None, '', 1, ''],\
						('header', 0, ): [None, '', 1, 'the header to be in output file, openMode=w'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']
						}
	def __init__(self, inputFname=None, **keywords):
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		if not self.inputFname:
			self.inputFname = inputFname
		if self.inputFname and self.inputFile is None:
			self.inputFile = utils.openGzipFile(self.inputFname, openMode=self.openMode)
		
		self.csvFile = None
		self.isRealCSV = False
		if self.openMode=='r':	#reading mode
			if self.delimiter is None:
				self.delimiter = figureOutDelimiter(self.inputFile)
			
			if self.delimiter=='\t' or self.delimiter==',':
				self.csvFile = csv.reader(self.inputFile, delimiter=self.delimiter)
				self.isRealCSV = True
			else:
				self.csvFile = self.inputFile
				self.isRealCSV = False
		else:	#writing mode
			if self.delimiter:
				self.csvFile = csv.writer(self.inputFile, delimiter=self.delimiter)
				self.isRealCSV = True
			else:
				self.csvFile = self.inputFile
				self.isRealCSV = False
		self.col_name2index = None
	
	def writeHeader(self):
		"""
		2012.11.16
		"""
		if self.header:
			if self.isRealCSV:
				self.csvFile.writerow(self.header)
			else:
				self.csvFile.write("%s\n"%(self.delimiter.join(self.header)))
	
	def constructColName2IndexFromHeader(self):
		"""
		2012.8.23
		"""
		self.header = self.next()
		self.col_name2index = utils.getColName2IndexFromHeader(self.header)
		return self.col_name2index
	
	def getHeader(self):
		"""
		2012.11.22
		"""
		return self.header
	
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
			row = self.csvFile.next()
		except:
			raise StopIteration
		if not self.isRealCSV:
			row = row.strip().split()
		return row
	
	def close(self):
		del self.csvFile
		self.inputFile.close()
	
	def writerow(self, row=None):
		"""
		mimic csv's writerow()
		"""
		if row:
			if self.isRealCSV:
				self.csvFile.writerow(row)
			else:
				self.csvFile.write("%s\n"%(self.delimiter.join(row)))
	
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