#!/usr/bin/env python
"""
Examples:
	%s  -s 1 -o 5566_association_locus_support_vs_overlap.png --whichColumnHeader no_of_peaks
		--minNoOfTotal=0  /tmp/5566_association_locus.h5  /tmp/5566_association_locus_overlap_3.h5
	
	%s  -s 1 -o 5566_association_locus_connectivity_vs_overlap.png
		--minNoOfTotal=0  /tmp/5566_association_locus.h5  /tmp/5566_association_locus_overlap_3.h5
	

Description:
	2012/11/22
		this program draws box plots (multiple) of association-locus attributes (connectivity, no_of_peaks, no_of_results, etc.) ,
			each boxplot with different min_overlap_ratio.
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
import numpy, random, pylab
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter
from pymodule import yh_matplotlib
from pymodule.plot.PlotBoxPlot import PlotBoxPlot

class BoxPlotAssociationLocusAttributeVsOverlapThreshold(PlotBoxPlot):
	__doc__ = __doc__
#						
	option_default_dict = PlotBoxPlot.option_default_dict.copy()
	option_default_dict.pop(('xColumnHeader', 1, ))
	# change default of whichColumnHeader
	option_default_dict[('whichColumnHeader', 0, )][0] = 'connectivity'
	option_default_dict[('xColumnPlotLabel', 0, )][0] = 'minimum overlap ratio'
	# change default of file-format
	option_default_dict[('inputFileFormat', 0, int)][0] = 2
	option_default_dict[('minNoOfTotal', 1, int)][0] = 1
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		PlotBoxPlot.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def processRow(self, row=None, pdata=None):
		"""
		2012.8.2
			handles each row in each file
		"""
		col_name2index = getattr(pdata, 'col_name2index', None)
		xValue2yValueLs = getattr(pdata, 'xValue2yValueLs', None)
		y_ls = getattr(pdata, 'y_ls', None)
		
		reader = getattr(pdata, 'reader', None)
		tableObject =reader.getTableObject()
		xValue = tableObject.getAttribute('min_overlap_ratio')
		if xValue2yValueLs is not None:
			#if self.whichColumnHeader:
			#	whichColumn = col_name2index.get(self.whichColumnHeader, None)
			#else:
			#	whichColumn = self.whichColumn
			xValue = self.processValue(xValue, processType=self.logX)
			yValue = row[self.whichColumnHeader]
			yValue = self.processValue(yValue, processType=self.logY)
			if xValue not in xValue2yValueLs:
				xValue2yValueLs[xValue] = []
			xValue2yValueLs[xValue].append(yValue)
			y_ls.append(yValue)
	

if __name__ == '__main__':
	main_class = BoxPlotAssociationLocusAttributeVsOverlapThreshold
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()