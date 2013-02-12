#!/usr/bin/env python
"""
Examples:
	#testing merge three identical genotype files
	%s -k 0 -v 3 -o test.tsv trio_inconsistency.tsv trio_inconsistency1.tsv trio_inconsistency2.tsv
	
	%s 
	
Description:
	2012.1.9
		This program takes mean/median/stdev of values of chosen columns (all input files) with same keys from the keyColumnLs.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])


sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, figureOutDelimiter, utils, PassingData
import csv, numpy
from ReduceMatrixByChosenColumn import ReduceMatrixByChosenColumn

class ReduceMatrixByAverageColumnsWithSameKey(ReduceMatrixByChosenColumn):
	__doc__ = __doc__
	option_default_dict = ReduceMatrixByChosenColumn.option_default_dict.copy()

	def __init__(self, inputFnameLs, **keywords):
		"""
		2011-7-12
		"""
		ReduceMatrixByChosenColumn.__init__(self, inputFnameLs, **keywords)
	
	
	def handleValueColumns(self, row, key2dataLs=None, keyColumnLs=[], valueColumnLs=[], **keywords):
		"""
		2012.1.9
		"""
		key = self.generateKey(row, keyColumnLs)
		if key not in key2dataLs:
			key2dataLs[key] = []	#
			for i in xrange(len(valueColumnLs)):	#[[]]*len(valueColumnLs) will fail. every [] is a reference to one [].
				key2dataLs[key].append([])
		for i in xrange(len(valueColumnLs)):
			columnIndex = valueColumnLs[i]
			if columnIndex<len(row):
				value = float(row[columnIndex])
				key2dataLs[key][i].append(value)
	
	def avgKey2DataLs(self, key2dataLs, no_of_key_columns=1, header=[]):
		"""
		2012.1.9
			1. take mean/median/stdev of every cell in dataLs,
			2. modify newHeader to reflect that
		"""
		sys.stderr.write("Averaging key2dataLs (%s entries ) ..."%(len(key2dataLs)))
		newKey2DataLs = {}
		newHeader = []
		keyColHeader = header[:no_of_key_columns]
		valueColHeader = header[no_of_key_columns:]
		newValueColHeader = []
		no_of_value_columns = len(valueColHeader)
		for i in xrange(no_of_value_columns):
			valueColName = valueColHeader[i]
			newValueColHeader += ['mean_%s'%(valueColName), 'median_%s'%(valueColName), 'stdev_%s'%(valueColName)]
		
		for key, dataLs in key2dataLs.iteritems():
			if key not in newKey2DataLs:
				newKey2DataLs[key] = []
			no_of_value_columns = len(dataLs)
			for i in xrange(no_of_value_columns):
				meanValue = numpy.mean(dataLs[i])
				medianValue = numpy.median(dataLs[i])
				stdev = numpy.std(dataLs[i])
				newKey2DataLs[key] += [meanValue, medianValue, stdev]
		sys.stderr.write("Done.\n")
		return PassingData(key2dataLs= newKey2DataLs, header=keyColHeader + newValueColHeader)
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		returnData = self.traverse()
		newReturnData = self.avgKey2DataLs(returnData.key2dataLs, no_of_key_columns=len(self.keyColumnLs), header=returnData.header)
		self.outputFinalData(self.outputFname, newReturnData.key2dataLs, returnData.delimiter, header=newReturnData.header)

if __name__ == '__main__':
	main_class = ReduceMatrixByAverageColumnsWithSameKey
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	import copy
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
