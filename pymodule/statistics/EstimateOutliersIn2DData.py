#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i /tmp/Contig315_StKitts_vs_Nevis.tsv --xColumnHeader=StKitts --whichColumnHeader=Nevis
		-s 1.0 -o /tmp/Contig315_StKitts_vs_Nevis.2D.png
	

Description:
	2012.10.12
		program to estimate how many outliers off the y=x axis.
		1. hard cutoff. abs(y-x)<=minDelta
		2. model y-x as a normal distribution, estimate its mean/variance
			then add them up as chi-squared statistic.
	If "-i ..." is given, it is regarded as one of the input files (plus the ones in trailing arguments). 
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
from pymodule import yh_matplotlib
import numpy, random
from pymodule.AbstractMatrixFileWalker import AbstractMatrixFileWalker
from pymodule.plot.AbstractPlot import AbstractPlot
from pymodule import statistics

class EstimateOutliersIn2DData(AbstractPlot):
	__doc__ = __doc__
	option_default_dict = AbstractPlot.option_default_dict.copy()
	option_default_dict.update({
						('minAbsDelta', 1, float): [None, '', 1, 'minimum of abs(y-x) for a sample to be declared an outlier'],\
						})
	#('columnForX', 1, ): ["", 'x', 1, 'index of the column to be x-axis'],\
	#('columnForY', 1, ): ["", 'y', 1, 'index of the column to be y-axis'],\
	#('whichColumnPlotLabel', 0, ), ('xColumnHeader', 1, ), ('xColumnPlotLabel', 0, )
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractPlot.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		
		
		self.x_value_ls = []
		self.y_value_ls = []
		self.z_value_ls = []
	
	def setup(self, **keywords):
		"""
		2012.10.15
			run before anything is run
		"""
		AbstractPlot.setup(self, **keywords)
		
		self.writer = self.invariantPData.writer
		self.noOfOutliers = 0
		self.chiSqStat = 0
		self.noOfNonMissing = 0 
		self.chiSqMinusLogPvalue = None
		self.dataLs = []
		
	def processHeader(self, header=None, pdata=None):
		"""
		2012.10.15
			called everytime the header of an input file is derived in fileWalker()
			Have to change this. The output header is different from input.
		"""
		if self.invariantPData.writer and not self.invariantPData.headerOutputted:
			newHeader = ["outputID", 'noOfOutliers', 'noOfNonMissing', 'outlierFraction', 'chiSqStat', 'chiSqMinusLogPvalue',\
						'xMedianValue', 'yMedianValue', 'corr']
			self.invariantPData.writer.writerow(newHeader)
			self.invariantPData.headerOutputted = True
		
	def processRow(self, row=None, pdata=None):
		"""
		2012.10.7
		"""
		returnValue = 0
		col_name2index = getattr(pdata, 'col_name2index', None)
		x_ls = getattr(pdata, 'x_ls', None)
		y_ls = getattr(pdata, 'y_ls', None)
		
		if col_name2index:
			x_index = col_name2index.get(self.xColumnHeader, None)
			if self.whichColumnHeader:
				y_index = col_name2index.get(self.whichColumnHeader, None)
			else:
				y_index = self.whichColumn
			
			xValue = row[x_index]
			yValue = row[y_index]
			
			if yValue!=self.missingDataNotation and xValue!=self.missingDataNotation:
				xValue = self.processValue(value=xValue, takeLogarithm=self.logX, positiveLog=self.positiveLog, \
										valueForNonPositiveValue=self.valueForNonPositiveYValue)
				yValue = self.processValue(yValue, takeLogarithm=self.logY, positiveLog=self.positiveLog, \
										valueForNonPositiveValue=self.valueForNonPositiveYValue)
				absDelta = abs(xValue - yValue)
				self.noOfNonMissing += 1
				self.dataLs.append(absDelta)
				self.x_value_ls.append(xValue)
				self.y_value_ls.append(yValue)
				
				if absDelta>=self.minAbsDelta:
					self.noOfOutliers += 1
				returnValue = 1
		return returnValue
	
	def afterFileFunction(self, **keywords):
		"""
		2012.10.7
		"""
		pass
		
	
	def reduce(self, **keywords):
		"""
		2012.10.15
			run after all files have been walked through
		"""
		meanStdData = statistics.estimateMeanStdFromData(dataVector=self.dataLs, excludeTopFraction=0.2)
		chiSqData = statistics.calculateChiSqStatOfDeltaVector(dataVector=self.dataLs, mean=meanStdData.mean, \
												std=meanStdData.std)
		xMedianValue =  numpy.median(self.x_value_ls)
		yMedianValue =  numpy.median(self.y_value_ls)
		corr = numpy.corrcoef(self.x_value_ls, self.y_value_ls)[0,1]
		
		thisOutputID = os.path.basename(self.outputFname)
		outlierFraction = float(self.noOfOutliers)/self.noOfNonMissing
		dataRow = [thisOutputID, self.noOfOutliers, self.noOfNonMissing, outlierFraction, chiSqData.chiSqStat,\
				chiSqData.chiSqMinusLogPvalue, xMedianValue, yMedianValue, corr]
		self.writer.writerow(dataRow)
		
		sys.stderr.write("%s/%s (%.3f) outliers, chiSqStat=%.3f, chiSqMinusLogPvalue=%.3f.\n"%\
						(self.noOfOutliers, self.noOfNonMissing, outlierFraction,\
						chiSqData.chiSqStat, chiSqData.chiSqMinusLogPvalue))
		
		#close the self.invariantPData.writer
		AbstractMatrixFileWalker.reduce(self, **keywords)
		
	
if __name__ == '__main__':
	main_class = EstimateOutliersIn2DData
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
