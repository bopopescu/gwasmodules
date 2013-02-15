#!/usr/bin/env python
"""
Examples:
	%s -s 1 -o /tmp/Contig685_largeSiteGap.tsv -V 5000  -W distanceToNextSite  /tmp/Contig685_siteGap.tsv
	
	# 2012.8.6 select rows based on the fraction of heterozygotes per individual [0.2, 0.8].
	%s -s 1 -o /tmp/hetPerMonkey_hist.png
		-W NoOfHet_by_NoOfTotal -V 0.2 -x 0.8 /tmp/homoHetCountPerSample.tsv
	

Description:
	2012.8.13 select certain rows from a matrix based on the whichColumn (chosen column)
		it accepts both input from "-i oneFile.txt" or the trailing standalone arguments.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])


#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:	   #64bit
#	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
#	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
#else:   #32bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))


from pymodule import ProcessOptions
from pymodule.yhio.AbstractMatrixFileWalker import AbstractMatrixFileWalker


class SelectRowsFromMatrix(AbstractMatrixFileWalker):
	__doc__ = __doc__
#						
	option_default_dict = AbstractMatrixFileWalker.option_default_dict
	option_default_dict.update({
			('minWhichColumnValue', 0, float): [None, 'V', 1, 'minimum value of the which column, after log transformation if needed.\
			Default is no such filter.'],\
			('maxWhichColumnValue', 0, float): [None, 'x', 1, 'maximum value of the which column, after log transformation if needed.\
			Default is no such filter.'],\
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractMatrixFileWalker.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def processRow(self, row=None, pdata=None):
		"""
		2012.8.13
			handles each row in each file
		"""
		col_name2index = getattr(pdata, 'col_name2index', None)
		y_ls = getattr(pdata, 'y_ls', None)	#don't add anything to it, then self.plot won't be called in self.fileWalker
		if col_name2index and y_ls is not None:
			if self.whichColumnHeader:
				whichColumn = col_name2index.get(self.whichColumnHeader, None)
			else:
				whichColumn = self.whichColumn
			
			yValue = self.handleYValue(row[whichColumn])
			if self.minWhichColumnValue is not None and yValue<self.minWhichColumnValue:
				return
			if self.maxWhichColumnValue is not None and yValue>self.maxWhichColumnValue:
				return
			self.invariantPData.writer.writerow(row)
	
	
if __name__ == '__main__':
	main_class = SelectRowsFromMatrix
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
