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
from pymodule.AbstractMatrixFileWalker import AbstractMatrixFileWalker

class AbstractPlot(AbstractMatrixFileWalker):
	__doc__ = __doc__
	option_default_dict = AbstractMatrixFileWalker.option_default_dict.copy()
	#option_default_dict.update(AbstractMapper.db_option_dict.copy())
	option_default_dict.update({
						('title', 0, ): [None, 't', 1, 'title for the figure.'],\
						('figureDPI', 1, int): [200, 'f', 1, 'dpi for the output figures (png)'],\
						('formatString', 1, ): ['.', 'm', 1, 'formatString passed to matplotlib plot'],\
						('ylim_type', 1, int): [1, 'y', 1, 'y-axis limit type, 1: whatever matplotlib decides. 2: min to max'],\
						('whichColumnPlotLabel', 0, ): ['', 'D', 1, 'plot label for data of the whichColumn', ],\
						('xColumnHeader', 1, ): ['', 'l', 1, 'header of the x-axis data column, ' ],\
						('xColumnPlotLabel', 0, ): ['', 'x', 1, 'x-axis label (posColumn) in manhattan plot', ],\
						
						('logX', 0, int): [0, '', 0, 'whether to take -log of X'],\
						('need_svg', 0, ): [0, 'n', 0, 'whether need svg output', ],\
						
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractMatrixFileWalker.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
		#if user wants to preserve data in a data structure that is visisble throughout reading different files.
		# then use this self.invariantPData.
		#AbstractMatrixFileWalker has initialized a structure like below.
		#self.invariantPData = PassingData()
	
	
	def plot(self, x_ls=None, y_ls=None, pdata=None):
		"""
		2011-9-30
			get called by the end of fileWalker() for each inputFname.
		"""
		pylab.plot(x_ls, y_ls, self.formatString)
	
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
			if yValue not in self.missingDataNotation and xValue not in self.missingDataNotation:
				xValue = float(xValue)
				yValue = self.handleYValue(yValue)
				x_ls.append(xValue)
				y_ls.append(yValue)
	
	def processHeader(self, header=None, pdata=None):
		"""
		2012.8.13
			called everytime the header of an input file is derived in fileWalker()
		"""
		pass
	
	def handleXLabel(self,):
		"""
		2012.8.6
		"""
		if getattr(self, 'xColumnPlotLabel', None):
			xlabel = self.xColumnPlotLabel
		else:
			xlabel = getattr(self, "xColumnHeader", "")
		pylab.xlabel(xlabel)
		return xlabel
		
	def handleYLabel(self,):
		"""
		2012.8.6
		"""
		if getattr(self, 'whichColumnPlotLabel', None):
			ylabel = self.whichColumnPlotLabel
		else:
			ylabel = getattr(self, "whichColumnHeader", "")
		pylab.ylabel(ylabel)
		return ylabel
	
	def handleTitle(self,):
		"""
		2012.8.16
		"""
		if self.title:
			title = self.title
		else:
			title = yh_matplotlib.constructTitleFromTwoDataSummaryStat(self.invariantPData.x_ls, self.invariantPData.y_ls)
		pylab.title(title)
		return title

	def saveFigure(self, invariantPData=None, **keywords):
		"""
		2012.10.7
		
		"""
		sys.stderr.write("Saving figure ...")
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
		
		self.reduce()
		

if __name__ == '__main__':
	main_class = AbstractPlot
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()