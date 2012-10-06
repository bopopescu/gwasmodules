#!/usr/bin/env python
"""
Examples:
	#-1.0 is NA, rather than -1
	%s -i -N "-1.0" -o tmp/Contig97_variantSiteStatAndTrioInconsistencyRate -x mismatchRate -y missingRate -z inconsistencyRateInTrio 
		-s 0.5 tmp/Contig97_variantSiteStatAndTrioInconsistencyRate.tsv
	
	%s -i  -o tmp/Contig97_variantSiteStatAndTrioInconsistencyRate
		-x mismatchRate -y missingRate -z inconsistencyRateInTrio
		-s 0.5 tmp/Contig97_variantSiteStatAndTrioInconsistencyRate.tsv
	

Description:
	2011-11-2
		The input is tab/coma-delimited, with a header and has at least 3 columns.
		The three designated columns must be of float value.
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

class Draw2DHistogramOfMatrix(object):
	__doc__ = __doc__
	option_default_dict = {('outputFnamePrefix', 1, ): [None, 'o', 1, 'output filename prefix for output figures.'],\
						('columnForX', 1, ): ["", 'x', 1, 'index of the column to be x-axis'],\
						('columnForY', 1, ): ["", 'y', 1, 'index of the column to be y-axis'],\
						('columnForZ', 1, ): ["", 'z', 1, 'index of the column to be z-axis'],\
						('logX', 0, int): [0, 'l', 0, 'whether to take -log of X'],\
						('logY', 0, int): [0, 'm', 0, 'whether to take -log of Y'],\
						('logZ', 0, int): [0, 'n', 0, 'whether to take -log of Z'],\
						('title', 0, ): [None, 't', 1, 'title for the figure.'],\
						('samplingRate', 1, float): [0.01, '', 1, 'how often you include the data'],\
						('figureDPI', 1, int): [200, '', 1, 'dpi for the output figures (png)'],\
						('NANotation', 1, ):["-1", '', 1, 'which represents NA , apart from empty string'],\
						('ignoreNA', 0, int):[0, 'i', 0, 'toggle this to ignore pairs with either missing'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	#option_default_dict.update({('maxDepth', 1, int): [50, 'm', 1, 'genotype depth of every member of the trio should be below this number']})
	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		self.inputFnameLs = inputFnameLs
		
		self.x_value_ls = []
		self.y_value_ls = []
		self.z_value_ls = []
		
	def trioInconsistentRateFileWalker(self, inputFname, processFunc=None, columnForX="", columnForY="", columnForZ="", \
									logX=0, logY=0, logZ=0, run_type=1):
		"""
		2011-11-2
			remove the maxDepth filter. apply afterwards through filterDataByDepth().
		2011-9-30
		
		"""
		reader = csv.reader(open(inputFname), delimiter=figureOutDelimiter(inputFname))
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		x_index = col_name2index.get(columnForX)
		y_index = col_name2index.get(columnForY)
		z_index = col_name2index.get(columnForZ)
		counter = 0
		real_counter = 0
		for row in reader:
			counter += 1
			if counter%50000==0:
				sys.stderr.write("%s\t%s\t%s"%('\x08'*120, counter, real_counter))
			r = random.random()
			# exit if one row has less data to support all 3 columns
			if len(row)<max(x_index, y_index, z_index)+1:
				sys.stderr.write("%s\t%s\t%s\n"%('\x08'*120, counter, real_counter))
				sys.stderr.write("Warning: the length of this row %s is beyond any of the indices (%s, %s, %s) here.\n"%\
								(len(row), x_index, y_index, z_index))
				sys.stderr.write(repr(row))
				sys.exit(3)
				continue
			if row[x_index] and row[y_index] and row[z_index] and r<=self.samplingRate:
				if self.ignoreNA:
					if row[x_index]==self.NANotation or row[y_index]==self.NANotation or row[z_index] ==self.NANotation:
						continue
				real_counter += 1
				x_value = float(row[x_index])
				y_value = float(row[y_index])
				z_value = float(row[z_index])
				if logX:	#handled by the xscale
					if x_value>0:
						x_value = -math.log10(x_value)
					else:
						x_value = 100
				if logY:
					if y_value>0:
						y_value = -math.log10(y_value)
					else:
						y_value = 100
				if logZ:
					if z_value>0:
						z_value = -math.log10(z_value)
					else:
						z_value = 100
				self.x_value_ls.append(x_value)
				self.y_value_ls.append(y_value)
				self.z_value_ls.append(z_value)
			
		sys.stderr.write("%s\t%s\t%s\n"%('\x08'*120, counter, real_counter))
		del reader
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		import pylab
		pylab.clf()
		
		if len(self.inputFnameLs)>1:	#multiple input file (multiple contigs). use lines.
			self.formatString = '-'
		else:	#only one input file (one contig). use dot.
			self.formatString = '.'
		
		sys.stderr.write("Reading in data ...\n")
		for inputFname in self.inputFnameLs:
			if os.path.isfile(inputFname):
				self.trioInconsistentRateFileWalker(inputFname, processFunc=None, columnForX=self.columnForX, columnForY=self.columnForY, \
												columnForZ=self.columnForZ, logX=self.logX, logY=self.logY, logZ=self.logZ)
		sys.stderr.write("Done.\n")
		
		
		if self.title is None:
			title = " %s input files"%(len(self.inputFnameLs))
		else:
			title = self.title
		
		outputFnamePrefix = self.outputFnamePrefix
		
		loci_count_ls = [1]*len(self.x_value_ls)
		colorBarLabelForZ = self.columnForZ
		colorBarLabelForLociCount = 'log10(no. of loci)'
		
		#if uniformly distributed over 2D, each hexagon has ~50 points.
		# the maximum gridsize is 30.
		gridsize = min(30, int(math.sqrt(len(loci_count_ls)/50.0)))
		
		#if self.logX:
		#	xscale = "log"	#it raises exception when negative value or zero is encountered
		#else:
		xscale = "linear"
		outputFname = '%s_%s_vs_%s_%s_as_Z.png'%(outputFnamePrefix, self.columnForX, self.columnForY, self.columnForZ)
		yh_matplotlib.drawHexbin(self.x_value_ls, self.y_value_ls, self.z_value_ls, \
								fig_fname=outputFname, gridsize=gridsize, title=title, xlabel=self.columnForX, ylabel=self.columnForY,\
								colorBarLabel=colorBarLabelForZ, reduce_C_function=numpy.median, dpi=self.figureDPI,\
								mincnt=50, marginals=False, xscale=xscale)	#at least 50 sites in one hexagon
		outputFname = '%s_%s_vs_%s_loci_count.png'%(outputFnamePrefix, self.columnForX, self.columnForY,)
		yh_matplotlib.drawHexbin(self.x_value_ls, self.y_value_ls, loci_count_ls, \
								fig_fname=outputFname, gridsize=gridsize, title=title, xlabel=self.columnForX, ylabel=self.columnForY,\
								colorBarLabel=colorBarLabelForLociCount, reduce_C_function=yh_matplotlib.logSum, dpi=self.figureDPI,\
								xscale=xscale)
		
if __name__ == '__main__':
	main_class = Draw2DHistogramOfMatrix
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
