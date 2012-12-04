#!/usr/bin/env python
"""
Examples:
	%s 
	
	# 2012.8.6 draw histogram of fraction of heterozygotes per site.
	%s  -i 10 -s 0.0001 -o /tmp/hetPerSite_hist.png
		-W "OBS(HOM1/HET/HOM2)"
		/Network/Data/vervet/vervetPipeline/VCFStat_Method9_L0P2Mil_m500K.2012.8.3T0136/VCFStatGzip/hweMerge.tsv.gz
	

Description:
	2012.8.6
		derivative of DrawHistogram.py, parse the OBS(HOM1/HET/HOM2) (i.e. 13/3/0 ) the column.
			This won't work in pegasus workflow as the column header has parenthesis ( ) in it and 
			i do not know how to escape them in pegasus. 
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
from DrawHistogram import DrawHistogram

class DrawVCFtoolsHWEHetFractionHistogram(DrawHistogram):
	__doc__ = __doc__
#						
	option_default_dict = DrawHistogram.option_default_dict
	
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		DrawHistogram.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def handleYValue(self, yValue=None):
		"""
		2012.8.6
			parse the OBS(HOM1/HET/HOM2) (i.e. 13/3/0 ) to get het fraction
		"""
		vector = yValue.split('/')
		vector = map(int, vector)
		noOfHomo1, noOfHet, noOfHomo2 = vector
		noOfTotal = sum(vector)
		if noOfTotal>0:
			yValue = float(noOfHet/float(noOfTotal))
		else:
			yValue = -1
		yValue = self.processValue(yValue, processType=self.logY)
		return yValue
	
	
if __name__ == '__main__':
	main_class = DrawVCFtoolsHWEHetFractionHistogram
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
