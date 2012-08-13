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
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import csv, random
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, getColName2IndexFromHeader, figureOutDelimiter
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
import pylab

class AbstractPlot(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 1, ))
	#option_default_dict.update(AbstractMapper.db_option_dict.copy())
	option_default_dict.update({
						('minNoOfTotal', 1, int): [100, 'i', 1, 'minimum no of total variants (denominator of inconsistent rate)'],\
						('title', 0, ): [None, 't', 1, 'title for the figure.'],\
						('figureDPI', 1, int): [200, 'f', 1, 'dpi for the output figures (png)'],\
						('formatString', 1, ): ['.', 'm', 1, 'formatString passed to matplotlib plot'],\
						('ylim_type', 1, int): [1, 'y', 1, 'y-axis limit type, 1: 0 to max. 2: min to max'],\
						('samplingRate', 1, float): [1, 's', 1, 'how often you include the data, a probability between 0 and 1.'],\
						('whichColumn', 0, int): [3, 'w', 1, 'data from this column (index starting from 0) is plotted as y-axis value'],\
						('whichColumnHeader', 0, ): ["", 'W', 1, 'column header for the data to be plotted as y-axis value, substitute whichColumn'],\
						('whichColumnPlotLabel', 0, ): ['', 'D', 1, 'plot label for data of the whichColumn', ],\
						('logWhichColumn', 0, int): [0, 'g', 0, 'whether to take -log of the whichColumn'],\
						('positiveLog', 0, int): [0, 'p', 0, 'toggle to take log, rather than -log(), \
				only effective when logWhichColumn is toggled. '],\
						('valueForNonPositiveYValue', 1, float): [50, 'v', 1, 'if the whichColumn value is not postive and logWhichColumn is on,\
			what yValue should be.'],\
						('xColumnHeader', 1, ): ['', 'l', 1, 'header of the x-axis data column, ' ],\
						('xColumnPlotLabel', 0, ): ['', 'x', 1, 'x-axis label (posColumn) in manhattan plot', ],\
						('need_svg', 0, ): [0, 'n', 0, 'whether need svg output', ],\
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
		#if user wants to preserve data in a data structure that is visisble throughout reading different files.
		# then use this self.invariantPData.
		self.invariantPData = PassingData()
		
	def connectDB(self):
		"""
			split out of __init__() so that derived classes could overwrite this function
		"""
		pass
	
	def initiatePassingData(self, ):
		"""
		2012.8.2
			this function gets called in the beginning of each fileWalker() (for each inputFname)
		"""
		pdata = PassingData(x_ls = [], y_ls = [], invariantPData=self.invariantPData)
		return pdata
	
	def plot(self, x_ls=None, y_ls=None, pdata=None):
		"""
		2011-9-30
			get called by the end of fileWalker() for each inputFname.
		"""
		pylab.plot(x_ls, y_ls, self.formatString)
	
	def handleYValue(self, yValue=None):
		"""
		2012.8.2
		"""
		yValue = float(yValue)
		if self.logWhichColumn:
			if yValue>0:
				if self.positiveLog:
					yValue = math.log10(yValue)
				else:
					yValue = -math.log10(yValue)
				
			else:
				yValue = self.valueForNonPositiveYValue
		return yValue
	
	def processRow(self, row=None, pdata=None):
		"""
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
			
			xValue = float(row[x_index])
			yValue = self.handleYValue(row[whichColumn])
			
			y_ls.append(yValue)
			x_ls.append(xValue)
	
	def getNumberOfData(self, pdata):
		"""
		2012.8.6
		"""
		return len(pdata.y_ls)
	
	def fileWalker(self, inputFname=None, plotFunction=None, processRowFunction=None , run_type=1):
		"""
		2012.8.1
		"""
		sys.stderr.write("walking through %s ..."%(inputFname))
		counter =0
		pdata = self.initiatePassingData()
		if processRowFunction is None:
			processRowFunction = self.processRow
		if plotFunction is None:
			plotFunction = self.plot
		try:
			inf = utils.openGzipFile(inputFname)
			delimiter = figureOutDelimiter(inputFname)
			isCSVReader = True
			if delimiter=='\t' or delimiter==',':
				reader = csv.reader(inf, delimiter=figureOutDelimiter(inputFname))
				header = reader.next()
			else:
				reader = inf
				header = inf.readline().strip().split()	#whatever splits them
				isCSVReader = False
			col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
			pdata.col_name2index = col_name2index
			
			for row in reader:
				if not isCSVReader:
					row = row.strip().split()
				if self.samplingRate<1 and self.samplingRate>=0:
					r = random.random()
					if r>self.samplingRate:
						continue
				processRowFunction(row=row, pdata=pdata)
				counter += 1
			if self.getNumberOfData(pdata)>self.minNoOfTotal:
				plotFunction(pdata.x_ls, pdata.y_ls, pdata=pdata)
			del reader
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
		sys.stderr.write("%s data.\n"%(counter))
	
	def handleXLabel(self,):
		"""
		2012.8.6
		"""
		if getattr(self, 'xColumnPlotLabel', None):
			pylab.xlabel(self.xColumnPlotLabel)
	
	def handleYLabel(self,):
		"""
		2012.8.6
		"""
		if getattr(self, 'whichColumnPlotLabel', None):
			pylab.ylabel(self.whichColumnPlotLabel)
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		pylab.clf()
		
		for inputFname in self.inputFnameLs:
			if os.path.isfile(inputFname):
				self.fileWalker(inputFname, plotFunction=self.plot, run_type=1, processRowFunction=self.processRow)
		
		if self.title:
			pylab.title(self.title)
		self.handleXLabel()
		self.handleYLabel()
		
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
		sys.stderr.write("Done.\n")

if __name__ == '__main__':
	main_class = AbstractPlot
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()