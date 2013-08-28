#!/usr/bin/env python
"""
Examples:
	%s 
	
	# 2012.8.2 draw data column NumberOfLoci (y-axis, -W) vs. AAC (x-axis -l ...).
	# sample all data (-s 1), generate svg figure (-n)
	# take positive (-p) log (-g) of the whichColumn value (y-axis)
	%s -p -g -l AAC -W NumberOfLoci -D NoOfLoci -O ~/NoOfLoci_vs_AAC -n -s 1 -x AAC
			VCFStat_Method8_L800000P4000000m1000000.2012.8.1T0331/11Contigs_AAC_tally.tsv

Description:
	2012.8.2
		abstract class for plot classes, can plot XY scatter/line (pending self.formatString) plot.
	If you specify --outputFname, make sure its suffix is .png.
	If "-i ..." is given, it is regarded as one of the input files (plus the ones in trailing arguments). 
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import pylab
import csv, random, numpy
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, getColName2IndexFromHeader, figureOutDelimiter,\
	yh_matplotlib
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
from pymodule.yhio.AbstractMatrixFileWalker import AbstractMatrixFileWalker

class AbstractPlot(AbstractMatrixFileWalker):
	__doc__ = __doc__
	option_default_dict = AbstractMatrixFileWalker.option_default_dict.copy()
	#option_default_dict.update(AbstractMapper.db_option_dict.copy())
	option_default_dict.update({
						('title', 0, ): [None, 't', 1, 'title for the figure.'],\
						('figureDPI', 1, int): [200, 'f', 1, 'dpi, dots per inch, for the output figures (png)'],\
						('formatString', 1, ): ['.', '', 1, 'formatString passed to matplotlib plot'],\
						('markerSize', 1, int): [10, '', 1, 'size of plotting marker (dot size, usually), matplotlib default is 5'],\
						('ylim_type', 1, int): [1, 'y', 1, 'y-axis limit type, 1: whatever matplotlib decides. 2: min to max'],\
						('whichColumnPlotLabel', 0, ): ['', 'D', 1, 'plot label for data of the whichColumn', ],\
						('xColumnHeader', 1, ): ['', 'l', 1, 'header of the x-axis data column, ' ],\
						('xColumnPlotLabel', 0, ): ['', 'x', 1, 'x-axis label (posColumn) in manhattan plot', ],\
						
						('logX', 0, int): [0, '', 1, 'value 0: nothing; 1: log(X), 2: -log(X)'],\
						('xScaleLog', 0, int): [0, '', 1, 'regarding the x-axis scale. 0: nothing; 1: scale of log(10), 2: scale of log(2),\
	if this is non-zero, then logX should be toggled off. otherwise mess up the data'],\
						('yScaleLog', 0, int): [0, '', 1, 'regarding the y-axis scale. 0: nothing; 1: scale of log(10), 2: scale of log(2),\
	if this is non-zero, then logY should be toggled off. otherwise data will be messed up.'],\
						('need_svg', 0, ): [0, 'n', 0, 'whether need svg output', ],\
						('defaultFontLabelSize', 1, int): [12, '', 1, 'default font & label size on the plot'], \
						('defaultFigureWidth', 1, float): [10, '', 1, 'default figure width (in inch) in the pedigree plot'], \
						('defaultFigureHeight', 1, float): [10, '', 1, 'default figure height (in inch) in the pedigree plot'], \
						('xmargin', 0, float): [0, '', 1, 'this is margin within the plotting area. how far data points should be away from axis. 0-1.0'], \
						('ymargin', 0, float): [0, '', 1, 'this is margin within the plotting area. how far data points should be away from axis. 0-1.0'], \
						('plotLeft', 0, float): [0.1, '', 1, 'the left side of the subplots of the figure'], \
						('plotRight', 0, float): [0.9, '', 1, 'the right side of the subplots of the figure'], \
						('plotBottom', 0, float): [0.1, '', 1, ' the bottom of the subplots of the figure'], \
						('plotTop', 0, float): [0.9, '', 1, ' the top of the subplots of the figure'], \
						('legendType', 0, int): [0, '', 1, '0: no legend; 1: legend'], \
						
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractMatrixFileWalker.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
		#if user wants to preserve data in a data structure that is visisble throughout reading different files.
		# then use this self.invariantPData.
		#AbstractMatrixFileWalker has initialized a structure like below.
		#self.invariantPData = PassingData()
		#2012.11.22
		self.xMin = None
		self.xMax = None
		self.yMin = None
		self.yMax = None
		#2012.11.25 for legend
		self.plotObjectLs = []
		self.plotObjectLegendLs = []
		
		#2013.07.01 set these attributes if the class or its descendants have valid values
		if getattr(self, 'defaultFontLabelSize', None):
			yh_matplotlib.setFontAndLabelSize(self.defaultFontLabelSize)
		if getattr(self, "defaultFigureWidth", None) and  getattr(self, "defaultFigureHeight", None):
			yh_matplotlib.setDefaultFigureSize((self.defaultFigureWidth, self.defaultFigureHeight))
		if getattr(self, "plotLeft", None) and  getattr(self, "plotBottom", None):
			yh_matplotlib.setPlotDimension(left=self.plotLeft, right=self.plotRight, bottom=self.plotBottom, top=self.plotTop)

	def addPlotLegend(self, plotObject=None, legend=None, pdata=None, **keywords):
		"""
		2012.11.25
		"""
		self.plotObjectLs.append(plotObject)
		self.plotObjectLegendLs.append(legend)
	
	def plot(self, x_ls=None, y_ls=None, pdata=None, **keywords):
		"""
		2011-9-30
			get called by the end of fileWalker() for each inputFname.
		"""
		plotObject = pylab.plot(x_ls, y_ls, self.formatString, markersize=self.markerSize, )[0]
		self.addPlotLegend(plotObject=plotObject, legend=os.path.basename(pdata.filename), pdata=pdata)
		
		self.setGlobalMinVariable(extremeVariableName='xMin', givenExtremeValue=min(x_ls))
		self.setGlobalMaxVariable(extremeVariableName='xMax', givenExtremeValue=max(x_ls))
		self.setGlobalMinVariable(extremeVariableName='yMin', givenExtremeValue=min(y_ls))
		self.setGlobalMaxVariable(extremeVariableName='yMax', givenExtremeValue=max(y_ls))
	
	
	def processRow(self, row=None, pdata=None):
		"""
		2012.8.31
			deal with missing data
		2012.8.2
			handles each row in each file
		"""
		col_name2index = getattr(pdata, 'col_name2index', None)
		x_ls = getattr(pdata, 'x_ls', None)
		y_ls = getattr(pdata, 'y_ls', None)
		if col_name2index and x_ls is not None and y_ls is not None:
			if self.whichColumnHeader:
				whichColumn = col_name2index.get(self.whichColumnHeader, None)
			else:
				whichColumn = self.whichColumn
			x_index = col_name2index.get(self.xColumnHeader, None)
			
			xValue = row[x_index]
			yValue = row[whichColumn]
			if xValue and yValue and yValue not in self.missingDataNotation and xValue not in self.missingDataNotation:
				xValue = self.processValue(xValue, processType=self.logX, valueForNonPositiveValue=self.valueForNonPositiveYValue)
				yValue = self.processValue(yValue, processType=self.logY, \
							valueForNonPositiveValue=self.valueForNonPositiveYValue)
				x_ls.append(xValue)
				y_ls.append(yValue)
	
	def setGlobalExtremeVariable(self, extremeVariableName=None, givenExtremeValue=None, extremeType=1):
		"""
		2012.11.22
			extremeType:
				1: min
				2: max
			
		"""
		if givenExtremeValue is not None:
			extremeVariable = getattr(self, extremeVariableName)
			if extremeVariable is None:
				setattr(self, extremeVariableName, givenExtremeValue)
			else:
				if extremeType==1:	#replace a min variable
					if givenExtremeValue<extremeVariable:
						setattr(self, extremeVariableName, givenExtremeValue)
				else:	#set a max variable
					if givenExtremeValue>extremeVariable:
						setattr(self, extremeVariableName, givenExtremeValue)
	
	def setGlobalMinVariable(self, extremeVariableName=None, givenExtremeValue=None):
		"""
		2012.11.22
		"""
		self.setGlobalExtremeVariable(extremeVariableName=extremeVariableName, givenExtremeValue=givenExtremeValue, \
									extremeType=1)
	
	def setGlobalMaxVariable(self, extremeVariableName=None, givenExtremeValue=None):
		"""
		2012.11.22
		"""
		self.setGlobalExtremeVariable(extremeVariableName=extremeVariableName, givenExtremeValue=givenExtremeValue, \
									extremeType=2)
		
	def processHeader(self, header=None, pdata=None):
		"""
		2012.8.13
			called everytime the header of an input file is derived in fileWalker()
		"""
		pass
	
	def _handleXLim(self, **keywords):
		"""
		2012.11.22
			could be used by descendants, but not called by default
		"""
		if self.xMin is not None and self.xMax is not None:
			delta = abs(self.xMax-self.xMin)/10.0
			if delta<=0:
				delta = 0.5
			pylab.xlim(xmin=self.xMin-delta, xmax=self.xMax+delta)
		
	def _handleYLim(self, **keywords):
		"""
		2012.11.22
			could be used by descendants, but not called by default
		"""
		if self.yMin is not None and self.yMax is not None:
			originalYMin, originalYMax = pylab.ylim()
			delta = abs(self.yMax-self.yMin)/10.0
			if delta<=0:
				delta = 0.5
			pylab.ylim(ymin=self.yMin-delta, ymax=self.yMax+delta)
	
	def handleXLim(self, **keywords):
		"""
		2012.11.22
		"""
		pass
		
	def handleYLim(self, **keywords):
		"""
		2012.11.22
		"""
		pass
	
	def handleXLabel(self, **keywords):
		"""
		2012.8.6
		"""
		if getattr(self, 'xColumnPlotLabel', None):
			xlabel = self.xColumnPlotLabel
		else:
			xlabel = getattr(self, "xColumnHeader", "")
		pylab.xlabel(xlabel)
		return xlabel
		
	def handleYLabel(self, **keywords):
		"""
		2012.8.6
		"""
		if getattr(self, 'whichColumnPlotLabel', None):
			ylabel = self.whichColumnPlotLabel
		else:
			ylabel = getattr(self, "whichColumnHeader", "")
		pylab.ylabel(ylabel)
		return ylabel
	
	def handleTitle(self, **keywords):
		"""
		2012.8.16
		"""
		if self.title:
			title = self.title
		else:
			title = yh_matplotlib.constructTitleFromTwoDataSummaryStat(self.invariantPData.x_ls, self.invariantPData.y_ls)
		pylab.title(title)
		return title
	
	def changeFigureScaleToLog(self, xScaleLog=None, yScaleLog=None):
		"""
		2012.11.24
			not used in default mode
		"""
		
		if xScaleLog is None:
			xScaleLog = self.xScaleLog
		if yScaleLog is None:
			yScaleLog = self.yScaleLog
		if xScaleLog==1:	#only change the axis into log-scale if it's positive log
			pylab.gca().set_xscale('log', basex=10)
		elif xScaleLog==2:
			pylab.gca().set_xscale('log', basex=2)
		
		if yScaleLog==1:	#only change the axis into log-scale if it's positive log
			pylab.gca().set_yscale('log', basey=10)
		elif yScaleLog==2:
			pylab.gca().set_yscale('log', basey=2)
		
	def saveFigure(self, invariantPData=None, **keywords):
		"""
		2013.07.16 add legend or not , depending on self.legendType
		2012.11.24
		2012.10.7
		"""
		sys.stderr.write("Saving figure ")
		#2012.12.4 change scale
		self.changeFigureScaleToLog()
		#2012.12.4 change margin
		self.setMargins()
		
		if self.legendType!=0 and self.plotObjectLegendLs and self.plotObjectLs:
			#add the legend
			pylab.legend(self.plotObjectLs, self.plotObjectLegendLs, shadow = True)
		
		if self.outputFnamePrefix:
			pngOutputFname = '%s.png'%self.outputFnamePrefix
			svgOutputFname = '%s.svg'%self.outputFnamePrefix
		elif self.outputFname:
			pngOutputFname = self.outputFname
			svgOutputFname = '%s.svg'%(self.outputFname[:-4])
		else:
			sys.stderr.write("could not get outputFnamePrefix from self.outputFnamePrefix %s or self.outputFname %s.\n"%\
							(self.outputFnamePrefix, self.outputFname))
			sys.exit(1)
		sys.stderr.write("to %s ..."%(pngOutputFname))
		pylab.savefig(pngOutputFname, dpi=self.figureDPI)
		if self.need_svg:
			pylab.savefig(svgOutputFname, dpi=self.figureDPI)
		sys.stderr.write("  .\n")

	def setup(self, **keywords):
		"""
		2012.10.15
			run before anything is run
		"""
		AbstractMatrixFileWalker.setup(self, **keywords)
		pylab.clf()
	
	def reduce(self, **keywords):
		"""
		2012.10.15
			run after all files have been walked through
		"""
		self.saveFigure(invariantPData=self.invariantPData)
		#delete self.invariantPData.writer if it exists
		AbstractMatrixFileWalker.reduce(self, **keywords)
	
	def setMargins(self):
		"""
		2012.12.4
		"""
		if self.xmargin is not None:
			pylab.gca().set_xmargin(self.xmargin)
		if self.ymargin is not None:
			pylab.gca().set_ymargin(self.ymargin)
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		self.setup()
		
		for inputFname in self.inputFnameLs:
			if os.path.isfile(inputFname):
				self.fileWalker(inputFname, afterFileFunction=None, run_type=1, processRowFunction=self.processRow)
				#afterFileFunction = None means self.plot
		
		self.handleTitle()
		self.handleXLabel()
		self.handleYLabel()
		self.handleXLim()
		self.handleYLim()
		self.reduce()
		
		self.closeFiles()
		

if __name__ == '__main__':
	main_class = AbstractPlot
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()