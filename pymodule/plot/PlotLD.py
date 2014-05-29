#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s  -s 0.01 -o ./LDPlot.png -W "R^2" --fitCurve /tmp/5988_VCF_Contig966.geno.ld
	

Description:
	2012.8.18
		this programs draw LD plot out of a matrix like output, i.e. vcftools's LD output. 
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
from pymodule import yh_matplotlib, GenomeDB, statistics
from pymodule.yhio.AbstractMatrixFileWalker import AbstractMatrixFileWalker
from AbstractPlot import AbstractPlot

class PlotLD(AbstractPlot):
	__doc__ = __doc__
	option_default_dict = AbstractPlot.option_default_dict.copy()
	option_default_dict.update({
			('chrLengthColumnHeader', 1, ): ['chrLength', 'c', 1, 'label of the chromosome length column', ],\
			('chrColumnHeader', 1, ): ['CHR', 'C', 1, 'label of the chromosome column', ],\
			('minChrLength', 1, int): [1000000, 'L', 1, 'minimum chromosome length for one chromosome to be included', ],\
			('pos2ColumnHeader', 1, ): ['POS2', 'u', 1, 'label of the 2nd position column, xColumnHeader is the 1st position column', ],\
			('maxDist', 0, int): [None, '', 1, 'if given, pairs beyond this distance are tossed.', ],\
			('minDist', 0, int): [None, '', 1, 'if given, pairs below this distance are tossed.', ],\
			('movingAverageType', 0, int): [2, '', 1, '1: median r2 within each step, 2: mean, 3: fraction that is >0.8.', ],\
#			('fitCurve', 0, ): [0, '', 0, 'toggle to fit an exponential decay function to the data', ],\
			})
	option_default_dict[('missingDataNotation', 0, )][0] = '-nan'
	
	#option_default_dict.pop(('outputFname', 1, ))
	
	"""
	option_for_DB_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgresql', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ]}
	option_default_dict.update(option_for_DB_dict)
	"""
	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		AbstractPlot.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def plot(self, x_ls=None, y_ls=None, pdata=None):
		"""
		2012.9.6 overwrite the ancestral function
		2011-9-30
			get called by the end of fileWalker() for each inputFname.
		"""
		#pylab.plot(x_ls, y_ls, self.formatString, alpha=0.3)
		#if self.fitCurve:
		movingAverageData = statistics.movingAverage(listOfList=[x_ls, y_ls], no_of_steps=20, needReorderData=True, \
										reduceType=self.movingAverageType, minValueForFraction=0.8)
		n_x_ls = movingAverageData.listOfList[0]
		n_y_ls = movingAverageData.listOfList[1]
		pylab.plot(n_x_ls, n_y_ls, self.formatString)
		pylab.ylim(ymin=0)
		
		
		#output the reduced x and y
		if self.outputFnamePrefix:
			prefix = self.outputFnamePrefix
		else:
			prefix = os.path.splitext(self.outputFname)[0]
		outputFname = '%s.tsv'%(prefix)
		outf = open(outputFname, 'a')
		writer = csv.writer(outf, delimiter='\t')
		header = ['distance', 'reducedR2']
		writer.writerow(header)
		for i in xrange(len(n_x_ls)):
			data_row = [n_x_ls[i], n_y_ls[i]]
			writer.writerow(data_row)
		del writer
		
		
		#pylab.ylim(ymax=1.0)
		#splineFitData = statistics.splineFit(x_ls=x_ls, y_ls=y_ls, no_of_steps=100, needReorderData=True)
		#pylab.plot(splineFitData.x_ls, splineFitData.y_ls)
		"""
		x_ar = numpy.array(x_ls, numpy.float)
		y_ar = numpy.array(y_ls, numpy.float)
		#sort x_ar and y_ar must be in the order of x_ar
		indexOfOrderList = numpy.argsort(x_ar)
		x_ar = x_ar[indexOfOrderList]
		y_ar = y_ar[indexOfOrderList]
		#from Bio.Statistics import lowess
		#lowessY = lowess.lowess(x_ar, y_ar)
		"""
		
		"""
		### error from this code:
		#	Warning: overflow encountered in multiply
		#			return a*numpy.exp(-b*x+d) + c
		
		from scipy.optimize import curve_fit
		def func(x, a, b, c, d):
			return a*numpy.exp(-b*x+d) + c
		
		popt, pcov = curve_fit(func, numpy.array(x_ls), numpy.array(y_ls))
		x_ar = numpy.linspace(0, 1, 25)
		functor = lambda x: func(x, popt[0], popt[1], popt[2], popt[3])
		y_ar = map(functor, x_ar)
		pylab.plot(x_ar, y_ar, 'r-')
		"""
	
	
	def processRow(self, row=None, pdata=None):
		"""
		2012.10.25
			handle self.minDist and self.maxDist, make sure pairs within this distance.
		2012.8.31 skip missing data via self.missingDataNotation
		2012.8.18
		"""
		col_name2index = getattr(pdata, 'col_name2index', None)
		chr_id_index = col_name2index.get(self.chrColumnHeader, None)
		chrLength_index = col_name2index.get(self.chrLengthColumnHeader, None)
		pos1_index = col_name2index.get(self.xColumnHeader, None)
		pos2_index = col_name2index.get(self.pos2ColumnHeader, None)
		x_ls = getattr(pdata, 'x_ls', None)
		y_ls = getattr(pdata, 'y_ls', None)
		if col_name2index and y_ls is not None:
			if self.whichColumnHeader:
				whichColumn = col_name2index.get(self.whichColumnHeader, None)
			elif self.whichColumn:
				whichColumn = self.whichColumn
			else:
				whichColumn = None
			if chrLength_index:
					chrLength = int(row[chrLength_index])
					if chrLength<self.minChrLength:
						return
			if whichColumn is not None and pos1_index is not None and pos2_index is not None:
				yValue = row[whichColumn]
				if yValue not in self.missingDataNotation:
					pos1 = int(float(row[pos1_index]))
					pos2 = int(float(row[pos2_index]))
					xValue = abs(pos2-pos1)
					if self.minDist and xValue<self.minDist:
						return
					if self.maxDist and xValue>self.maxDist:
						return
					xValue = self.processValue(xValue, processType=self.logX)
					x_ls.append(xValue)
					yValue = self.processValue(yValue, processType=self.logY)
					y_ls.append(yValue)
		return 1
	
if __name__ == '__main__':
	main_class = PlotLD
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()