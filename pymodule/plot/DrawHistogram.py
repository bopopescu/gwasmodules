#!/usr/bin/env python
"""
Examples:
	%s 
	
	# 2012.8.6 draw histogram of fraction of heterozygotes per individual.
	%s -i 10 -s 1 -o /tmp/hetPerMonkey_hist.png
		-W NoOfHet_by_NoOfTotal /tmp/homoHetCountPerSample.tsv
	

Description:
	2012.8.6
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

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter
from pymodule import yh_matplotlib, GenomeDB
import numpy, random, pylab
from AbstractPlot import AbstractPlot

class DrawHistogram(AbstractPlot):
	__doc__ = __doc__
#						
	option_default_dict = AbstractPlot.option_default_dict
	option_default_dict.pop(('xColumnHeader', 1, ))
	option_default_dict.pop(('xColumnPlotLabel', 0, ))
	
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractPlot.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def plot(self, x_ls=None, y_ls=None, pdata=None, min_no_of_data_points=10, needLog=False, max_no_of_bins=20, **kwargs):
		"""
		2011-8.6
			get called by the end of fileWalker() for each inputFname.
			needLog doesn't need to be True. as self.logWhichColumn and self.positiveLog will decide on how to log the y-value.
		"""
		no_of_data_points = len(y_ls)
		if no_of_data_points>=min_no_of_data_points:
			no_of_bins = max(10, min(max_no_of_bins, no_of_data_points/10))
			n, bins, patches = pylab.hist(y_ls, no_of_bins, log=needLog)
	
	def processRow(self, row=None, pdata=None):
		"""
		2012.8.2
			handles each row in each file
		"""
		col_name2index = getattr(pdata, 'col_name2index', None)
		y_ls = getattr(pdata, 'y_ls', None)
		if col_name2index and y_ls is not None:
			if self.whichColumnHeader:
				whichColumn = col_name2index.get(self.whichColumnHeader, None)
			else:
				whichColumn = self.whichColumn
			
			yValue = self.handleYValue(row[whichColumn])
			y_ls.append(yValue)
	
	def handleXLabel(self,):
		"""
		2012.8.6
			the whichColumnPlotLabel is usually reserved for y-axis, but in histogram, it's on x-axis.
		"""
		if getattr(self, 'whichColumnPlotLabel', None):
			pylab.xlabel(self.whichColumnPlotLabel)
	
	def handleYLabel(self,):
		"""
		2012.8.6
		"""
		pass
	
if __name__ == '__main__':
	main_class = DrawHistogram
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
