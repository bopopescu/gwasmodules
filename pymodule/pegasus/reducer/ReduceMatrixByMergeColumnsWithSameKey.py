#!/usr/bin/env python
"""
Examples:
	#testing merge three identical genotype files
	%s -o /tmp/ccc.tsv /tmp/call_1.tsv /tmp/call_1.tsv /tmp/call_1.tsv
	
	%s 
	
Description:
	2011-12-21
		This program merge lines keyed by keyColumnLs. The non-key columns are appended next to each other.
		If one input misses lines for some keys, those lines will have empty data over there.
		
		All input files must have the keys at the same column(s).
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])


sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from pymodule import ProcessOptions, figureOutDelimiter, utils, PassingData
from pymodule.yhio.MatrixFile import MatrixFile
from AbstractReducer import AbstractReducer

class ReduceMatrixByMergeColumnsWithSameKey(AbstractReducer):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractReducer.option_default_dict)
	option_default_dict.update({
						("keyColumnLs", 1, ): [0, 'k', 1, 'index(es) of the key in each input file. must be same. comma/dash-separated. i.e. 0-2,4 '],\
						("keyHeaderLs", 0, ): [0, '', 1, 'header(s) of the key. comma-separated'],\
						('valueColumnLs', 1, ):["1", 'v', 1, 'comma/tab-separated list, specifying columns from which to aggregate total value by key'],\
						})

	def __init__(self, inputFnameLs=None, **keywords):
		"""
		2011-7-12
		"""
		AbstractReducer.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		if self.keyColumnLs:
			self.keyColumnLs = utils.getListOutOfStr(self.keyColumnLs, data_type=int)
		else:
			self.keyColumnLs = []
		if self.keyHeaderLs:
			self.keyHeaderLs = self.keyHeaderLs.split(',')
		else:
			self.keyHeaderLs = []
		
		self.keyColumnSet = set(self.keyColumnLs)
	
	def appendSelectedCellIntoGivenList(self, givenLs=[], inputLs=[], indexLs=[]):
		"""
		2012.1.9
		"""
		for columnIndex in indexLs:
			if columnIndex<len(inputLs):
				givenLs.append(inputLs[columnIndex])
		return givenLs
	
	def generateKey(self, row, keyColumnLs):
		"""
		2012.1.17
			make sure columnIndex is >=0
		2012.1.9
		"""
		keyLs = []
		for columnIndex in keyColumnLs:
			if columnIndex<len(row) and columnIndex>=0:
				keyLs.append(row[columnIndex])
		key = tuple(keyLs)
		return key
	
	def outputFinalData(self, outputFname, key2dataLs=None, delimiter=None, header=None):
		"""
		2013.3.3 bugfix , added openMode='w' for MatrixFile()
		2013.2.12 replace csv.writer with MatrixFile
		2012.7.30
			open the outputFname regardless whether there is data or not.
		2012.1.9
		"""
		writer = MatrixFile(inputFname=outputFname, delimiter=delimiter, openMode='w')
		if key2dataLs and delimiter and header:
			if header:
				writer.writerow(header)
			keyLs = key2dataLs.keys()
			keyLs.sort()
			for key in keyLs:
				dataLs = key2dataLs.get(key)
				writer.writerow(list(key) + dataLs)
		writer.close()
	
	def handleNewHeader(self, oldHeader=None, newHeader=None, keyColumnLs=None, valueColumnLs=None, keyColumnSet=None):
		"""
		2012.1.9
		"""
		originalHeaderLength = len(oldHeader)
		if len(newHeader)==0:	#add the key columns into the new header
			self.appendSelectedCellIntoGivenList(newHeader, oldHeader, keyColumnLs)
		for i in xrange(originalHeaderLength):
			if i not in keyColumnSet:
				valueColumnLs.append(i)
		self.appendSelectedCellIntoGivenList(newHeader, oldHeader, valueColumnLs)
		return newHeader
	
	def handleValueColumns(self, row, key2dataLs=None, keyColumnLs=[], valueColumnLs=[], noOfDataColumnsFromPriorFiles=None, \
						visitedKeySet=None):
		"""
		2012.1.9
		"""
		key = self.generateKey(row, keyColumnLs)
		if key not in key2dataLs:
			key2dataLs[key] = ['']*noOfDataColumnsFromPriorFiles
		visitedKeySet.add(key)
		
		for columnIndex in valueColumnLs:
			key2dataLs[key].append(row[columnIndex])

	def traverse(self):
		"""
		self.noHeader:	#2012.8.10
		2012.1.9
		"""
		newHeader = []
		key2dataLs = {}	#key is the keyColumn, dataLs corresponds to the sum of each column from valueColumnLs 
		delimiter = None
		noOfDataColumnsFromPriorFiles = 0
		for inputFname in self.inputFnameLs:
			if not os.path.isfile(inputFname):
				if self.exitNonZeroIfAnyInputFileInexistent:
					sys.exit(3)
				else:
					continue
			reader = None
			try:
				inputFile = utils.openGzipFile(inputFname)
				delimiter = figureOutDelimiter(inputFile)
				reader = MatrixFile(inputFile=inputFile, delimiter=delimiter)
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
			
			valueColumnLs = []
			try:
				header = reader.next()
				self.handleNewHeader(header, newHeader, self.keyColumnLs, valueColumnLs, keyColumnSet=self.keyColumnSet)
				if self.noHeader:	#2012.8.10
					inputFile.seek(0)
					reader = MatrixFile(inputFile=inputFile, delimiter=delimiter)
			except:	#in case something wrong (i.e. file is empty)
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
			
			if reader is not None and valueColumnLs:
				visitedKeySet = set()
				for row in reader:
					try:
						self.handleValueColumns(row, key2dataLs=key2dataLs, keyColumnLs=self.keyColumnLs, \
								valueColumnLs=valueColumnLs, noOfDataColumnsFromPriorFiles=noOfDataColumnsFromPriorFiles, \
								visitedKeySet=visitedKeySet)
					except:	#in case something wrong (i.e. file is empty)
						sys.stderr.write('Ignore this row: %s.\n'%repr(row))
						sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
						import traceback
						traceback.print_exc()
				del reader
				#append empty data to keys who are not present in this current "reader" file
				totalKeySet = set(key2dataLs.keys())
				unvisitedKeySet = totalKeySet - visitedKeySet
				for key in unvisitedKeySet:
					for i in valueColumnLs:
						key2dataLs[key].append('')
			noOfDataColumnsFromPriorFiles += len(valueColumnLs)
		if self.noHeader:	#2012.8.10
			newHeader = None
		returnData = PassingData(key2dataLs=key2dataLs, delimiter=delimiter, header=newHeader)
		return returnData
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		returnData = self.traverse()
		self.outputFinalData(self.outputFname, key2dataLs=returnData.key2dataLs, delimiter=returnData.delimiter, header=returnData.header)

if __name__ == '__main__':
	main_class = ReduceMatrixByMergeColumnsWithSameKey
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	import copy
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
