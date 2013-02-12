#!/usr/bin/env python
"""
Examples:
	#testing merge three identical genotype files
	%s -k 0 -v 4,5 -o /tmp/test.tsv trio_inconsistency_summary_hist_homo_het.tsv
	
	%s
	
Description:
	2012.1.9
		This program first sums values of chosen columns (all input files) with same keys from the keyColumnLs.
		In the end, it divides values from first two chosen columns and appends it in the output as an extra column.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])


sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, figureOutDelimiter, utils, PassingData
import csv, numpy
from ReduceMatrixByChosenColumn import ReduceMatrixByChosenColumn

class ReduceMatrixBySumSameKeyColsAndThenDivide(ReduceMatrixByChosenColumn):
	__doc__ = __doc__
	option_default_dict = ReduceMatrixByChosenColumn.option_default_dict.copy()
	option_default_dict.update({
						("operatorType", 1, int): [1, 'p', 1, 'For the last column, 1: firstValue/2ndValue. 2: 1stValue-2ndValue.'],\
						})
	def __init__(self, inputFnameLs, **keywords):
		"""
		2011-7-12
		"""
		ReduceMatrixByChosenColumn.__init__(self, inputFnameLs, **keywords)
	
	
	def divideTwoColumnsInKey2DataLs(self, key2dataLs, no_of_key_columns=1, header=[]):
		"""
		2012.1.9
			1. generate a new column by dividing the values from 1st two cells in dataLs,
			2. modify newHeader to reflect that
		"""
		sys.stderr.write("Averaging key2dataLs (%s entries ) ..."%(len(key2dataLs)))
		keyColHeader = header[:no_of_key_columns]
		valueColHeader = header[no_of_key_columns:]
		if len(valueColHeader)>1:
			header.append('%s_by_%s'%(valueColHeader[0], valueColHeader[1]))
		for key, dataLs in key2dataLs.iteritems():
			no_of_value_columns = len(dataLs)
			if no_of_value_columns>1:
				if self.operatorType==2:
					ratio = float(dataLs[0]) - float(dataLs[1])
				else:
					if dataLs[1]!=0:
						ratio = dataLs[0]/float(dataLs[1])
					else:
						ratio = -1
				key2dataLs[key].append(ratio)
		sys.stderr.write("Done.\n")
		return PassingData(key2dataLs= key2dataLs, header=header)
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		returnData = self.traverse()
		newReturnData = self.divideTwoColumnsInKey2DataLs(returnData.key2dataLs, no_of_key_columns=len(self.keyColumnLs), header=returnData.header)
		self.outputFinalData(self.outputFname, newReturnData.key2dataLs, returnData.delimiter, header=newReturnData.header)

if __name__ == '__main__':
	main_class = ReduceMatrixBySumSameKeyColsAndThenDivide
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	import copy
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
