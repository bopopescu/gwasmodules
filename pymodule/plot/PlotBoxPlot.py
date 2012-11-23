#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s  -s 1 -o ../plot/AAC_tally_barChart.png -W NumberOfLoci -l AAC -x AAC -i 20 -D NumberOfLoci ./AAC_tally.tsv.gz
	

Description:
	2012.11.22
		this program draws multiple box plots of yValue. each boxplot's data has the same xValue (on X-axis).
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

class PlotBoxPlot(AbstractPlot):
	__doc__ = __doc__
#						
	option_default_dict = AbstractPlot.option_default_dict.copy()
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractPlot.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		#2012.11.22
		self.boxPlotXLabelList = []
		
	def initiatePassingData(self, ):
		"""
		2012.11.22
			this function gets called in the beginning of each fileWalker() (for each inputFname)
		"""
		pdata = PassingData(xValue2yValueLs={}, x_ls=[], y_ls=[], invariantPData=self.invariantPData)
		#2012.8.16 pass to global data
		self.invariantPData.y_ls = pdata.y_ls
		self.invariantPData.x_ls = pdata.x_ls
		return pdata
	
	def getNumberOfData(self, pdata):
		"""
		2012.11.22
		"""
		return len(pdata.xValue2yValueLs)
	
	def handleXLim(self,):
		"""
		2012.11.22
		"""
		pylab.xlim(xmin=self.xMin-1, xmax=self.xMax+1)
		
	def handleYLim(self,):
		"""
		2012.11.22
		"""
		self._handleYLim()
	
	def plot(self, x_ls=None, y_ls=None, pdata=None, color_ls=['r', 'g', 'b', 'y','c', 'k'],width = 0.5, **kwargs):
		"""
		2012.11.22
			get called by the end of fileWalker() for each inputFname.
		"""
		xValue2yValueLs = getattr(pdata, 'xValue2yValueLs', None)
		xValueList = []
		y2DList = []
		noOfBoxPlotsBefore = len(self.boxPlotXLabelList)
		for xValue, yValueList in xValue2yValueLs.iteritems():
			xValueList.append(xValue)
			self.boxPlotXLabelList.append(xValue)
			y2DList.append(yValueList)
			yMin = min(yValueList)
			self.setGlobalMinVariable(extremeVariableName='yMin', givenExtremeValue=yMin)
			yMax = max(yValueList)
			self.setGlobalMaxVariable(extremeVariableName='yMax', givenExtremeValue=yMax)
		rects_ls = []
		
		#width = 0.75/len(yValue_list)	#these are for multi-series of bar charts.
		#c = color_ls[i%len(color_ls)]
		#pylab.boxplot(y2DList, positions=xValueList)	#this changes the x-location of the box plot
		noOfBoxPlotsNow = len(self.boxPlotXLabelList)
		pylab.boxplot(y2DList, positions=range(noOfBoxPlotsBefore, noOfBoxPlotsNow))
		
		self.setGlobalMinVariable(extremeVariableName='xMin', givenExtremeValue=0)
		self.setGlobalMaxVariable(extremeVariableName='xMax', givenExtremeValue=noOfBoxPlotsNow-1)
	
	def processRow(self, row=None, pdata=None):
		"""
		2012.11.22
			handles each row in each file
		"""
		col_name2index = getattr(pdata, 'col_name2index', None)
		xValue2yValueLs = getattr(pdata, 'xValue2yValueLs', None)
		y_ls = getattr(pdata, 'y_ls', None)
		if col_name2index and xValue2yValueLs is not None:
			if self.whichColumnHeader:
				whichColumn = col_name2index.get(self.whichColumnHeader, None)
			else:
				whichColumn = self.whichColumn
			x_index = col_name2index.get(self.xColumnHeader, None)
			
			xValue = float(row[x_index])
			yValue = float(row[whichColumn])
			yValue = self.handleYValue(yValue)
			if xValue not in xValue2yValueLs:
				xValue2yValueLs[xValue] = []
			xValue2yValueLs[xValue].append(yValue)
			y_ls.append(yValue)
	
	def saveFigure(self, invariantPData=None, **keywords):
		"""
		2012.11.22
			before saving the figure, mark the xticks
		"""
		pylab.xticks(range(len(self.boxPlotXLabelList)), self.boxPlotXLabelList)
		AbstractPlot.saveFigure(self, invariantPData=invariantPData, **keywords)

if __name__ == '__main__':
	main_class = PlotBoxPlot
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()