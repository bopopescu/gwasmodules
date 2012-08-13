#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s  -s 0.01 -o ./LDPlot.png -W "R^2" /tmp/5988_VCF_Contig966.geno.ld
	

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
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter
from pymodule import yh_matplotlib, GenomeDB
import numpy, random, pylab
from AbstractPlot import AbstractPlot

class PlotLD(AbstractPlot):
	__doc__ = __doc__
#						
	option_default_dict = {
			('outputFname', 1, ): [None, 'o', 1, 'output file for the figure.'],\
			('minNoOfTotal', 1, int): [100, 'i', 1, 'minimum no of total variants (denominator of inconsistent rate)'],\
			('title', 0, ): [None, 't', 1, 'title for the figure.'],\
			('figureDPI', 1, int): [200, 'f', 1, 'dpi for the output figures (png)'],\
			('formatString', 1, ): ['.', 'a', 1, 'formatString passed to matplotlib plot'],\
			('ylim_type', 1, int): [1, 'y', 1, 'y-axis limit type, 1: 0 to max. 2: min to max'],\
			('samplingRate', 1, float): [0.001, 's', 1, 'how often you include the data'],\
			('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
			('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.'],\
			('whichColumn', 0, int): [3, 'w', 1, 'data from this column (index starting from 0) is plotted as y-axis value'],\
			('whichColumnLabel', 0, ): ["R^2", 'W', 1, 'column label (in the header) for the data to be plotted as y-axis value, substitute whichColumn'],\
			('logWhichColumn', 0, int): [0, 'g', 0, 'whether to take -log of the whichColumn'],\
			('whichColumnPlotLabel', 1, ): ['r2', 'D', 1, 'plot label for data of the whichColumn', ],\
			('chrLengthColumnLabel', 1, ): ['chrLength', 'c', 1, 'label of the chromosome length column', ],\
			('chrColumnLabel', 1, ): ['CHR', 'C', 1, 'label of the chromosome column', ],\
			('minChrLength', 1, int): [1000000, 'm', 1, 'minimum chromosome length for one chromosome to be included', ],\
			('pos1ColumnLabel', 1, ): ['POS1', 'l', 1, 'label of the 1st position column', ],\
			('pos2ColumnLabel', 1, ): ['POS2', 'p', 1, 'label of the 2nd position column', ],\
			('posColumnPlotLabel', 1, ): ['distance', 'x', 1, 'x-axis label in  plot', ],\
			}
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
	
	
	def fileWalker(self, inputFname, plotFunction=None, run_type=1, \
					chrColumnLabel='CHR', minChrLength=1000000, chrLengthColumnLabel='chrLength',\
					pos1ColumnLabel="POS1", pos2ColumnLabel='POS2'):
		"""
		2012.8.1
		"""
		sys.stderr.write("walking through %s ..."%(inputFname))
		counter =0
		try:
			reader = csv.reader(open(inputFname), delimiter=figureOutDelimiter(inputFname))
			header = reader.next()
			col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
			chr_id_index = col_name2index.get(chrColumnLabel, None)
			pos1_index = col_name2index.get(pos1ColumnLabel, None)
			pos2_index = col_name2index.get(pos2ColumnLabel, None)
			chrLength_index = col_name2index.get(chrLengthColumnLabel, None)
			if self.whichColumnLabel:
				whichColumn = col_name2index.get(self.whichColumnLabel, None)
			else:
				whichColumn = self.whichColumn
			
			x_ls = []
			y_ls = []
			for row in reader:
				if self.samplingRate<1 and self.samplingRate>=0:
					r = random.random()
					if r>self.samplingRate:
						continue
				if chrLength_index:
					chrLength = int(row[chrLength_index])
					if chrLength<minChrLength:
						continue
				chr_id = row[chr_id_index]
				pos1 = int(float(row[pos1_index]))
				pos2 = int(float(row[pos2_index]))
				yValue = float(row[whichColumn])
				if self.logWhichColumn:
					if yValue>0:
						yValue = -math.log10(yValue)
					else:
						yValue = 50
				
				y_ls.append(yValue)
				x_data = abs(pos2-pos1)
				x_ls.append(x_data)
				counter += 1
			if len(x_ls)>self.minNoOfTotal:
				plotFunction(x_ls, y_ls)
			del reader
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
		sys.stderr.write("%s data.\n"%(counter))
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		pylab.clf()
		
		for inputFname in self.inputFnameLs:
			if os.path.isfile(inputFname):
				self.fileWalker(inputFname, plotFunction=self.plot, run_type=1, \
					chrColumnLabel=self.chrColumnLabel, minChrLength=self.minChrLength, \
					chrLengthColumnLabel=self.chrLengthColumnLabel,\
					pos1ColumnLabel=self.pos1ColumnLabel, pos2ColumnLabel=self.pos2ColumnLabel)
		
		if self.title:
			pylab.title(self.title)
		pylab.xlabel(self.posColumnPlotLabel)
		pylab.ylabel(self.whichColumnPlotLabel)
		pylab.savefig(self.outputFname, dpi=self.figureDPI)
		sys.stderr.write("Done.\n")

if __name__ == '__main__':
	main_class = PlotLD
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()