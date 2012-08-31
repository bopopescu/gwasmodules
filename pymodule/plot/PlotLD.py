#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s  -s 0.01 -o ./LDPlot.png -W "R^2" /tmp/5988_VCF_Contig966.geno.ld
	

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
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter
from pymodule import yh_matplotlib, GenomeDB
import numpy, random, pylab
from AbstractPlot import AbstractPlot
from pymodule.AbstractMatrixFileWalker import AbstractMatrixFileWalker

class PlotLD(AbstractPlot):
	__doc__ = __doc__
	option_default_dict = AbstractPlot.option_default_dict.copy()
	option_default_dict.update({
			('chrLengthColumnHeader', 1, ): ['chrLength', 'c', 1, 'label of the chromosome length column', ],\
			('chrColumnHeader', 1, ): ['CHR', 'C', 1, 'label of the chromosome column', ],\
			('minChrLength', 1, int): [1000000, 'L', 1, 'minimum chromosome length for one chromosome to be included', ],\
			('pos2ColumnHeader', 1, ): ['POS2', 'u', 1, 'label of the 2nd position column, xColumnHeader is the 1st position column', ],\
			})
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
	
	
	
	def processRow(self, row=None, pdata=None):
		"""
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
				if yValue!=self.missingDataNotation:
					pos1 = int(float(row[pos1_index]))
					pos2 = int(float(row[pos2_index]))
					xValue = abs(pos2-pos1)
					x_ls.append(xValue)
					yValue = self.handleYValue(yValue)
					y_ls.append(yValue)
	
	
if __name__ == '__main__':
	main_class = PlotLD
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()