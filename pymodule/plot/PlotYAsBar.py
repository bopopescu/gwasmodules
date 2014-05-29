#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s  -s 1 -o ../plot/AAC_tally_barChart.png -W NumberOfLoci -l AAC -x AAC -i 20 -D NumberOfLoci ./AAC_tally.tsv.gz
	

Description:
	2011-11-28
		this program draws a manhattan plot (gwas plot) and a histogram for some vcftools outputted windowed statistics.
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
from AbstractPlot import AbstractPlot

class PlotYAsBar(AbstractPlot):
	__doc__ = __doc__
#						
	option_default_dict = AbstractPlot.option_default_dict.copy()
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractPlot.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
		
	def initiatePassingData(self, ):
		"""
		2012.8.2
			this function gets called in the beginning of each fileWalker() (for each inputFname)
		"""
		pdata = PassingData(xValue2yValue={}, x_ls=[], y_ls=[], invariantPData=self.invariantPData)
		#2012.8.16 pass to global data
		self.invariantPData.y_ls = pdata.y_ls
		self.invariantPData.x_ls = pdata.x_ls
		return pdata
	
	def getNumberOfData(self, pdata):
		"""
		2012.8.6
		"""
		return len(pdata.xValue2yValue)
	
	def plot(self, x_ls=None, y_ls=None, pdata=None, color_ls=['r', 'g', 'b', 'y','c', 'k'],width = 0.5, **kwargs):
		"""
		2011-9-30
			get called by the end of fileWalker() for each inputFname.
		"""
		xValue2yValue = getattr(pdata, 'xValue2yValue', None)
		xValue_list = xValue2yValue.keys()
		xValue_list.sort()
		yValue_list = [xValue2yValue[item] for item in xValue_list]
		
		#width = 0.75/len(yValue_list)	#these are for multi-series of bar charts.
		#c = color_ls[i%len(color_ls)]
		c = color_ls[0]
		ind = numpy.array(xValue_list)
		plotObject = pylab.bar(ind, yValue_list, width=width, bottom=0, color=c, log=False, **kwargs)
		self.plotObjectLs.append(plotObject)
		self.plotObjectLegendLs.append(os.path.basename(pdata.filename))
		
		#if len(rects_ls)==2:
		#	pylab.legend( (rects_ls[0][0], rects_ls[1][0]), ('version 0', 'version 1') )
		
		#pylab.plot(x_ls, y_ls, self.formatString)
	
	def processRow(self, row=None, pdata=None):
		"""
		2012.8.2
			handles each row in each file
		"""
		col_name2index = getattr(pdata, 'col_name2index', None)
		xValue2yValue = getattr(pdata, 'xValue2yValue', None)
		y_ls = getattr(pdata, 'y_ls', None)
		if col_name2index and xValue2yValue is not None:
			if self.whichColumnHeader:
				whichColumn = col_name2index.get(self.whichColumnHeader, None)
			else:
				whichColumn = self.whichColumn
			x_index = col_name2index.get(self.xColumnHeader, None)
			
			xValue = float(row[x_index])
			xValue = self.processValue(xValue, processType=self.logX)
			yValue = float(row[whichColumn])
			yValue = self.processValue(yValue, processType=self.logY)
			xValue2yValue[xValue] = yValue
			y_ls.append(yValue)
	

if __name__ == '__main__':
	main_class = PlotYAsBar
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()