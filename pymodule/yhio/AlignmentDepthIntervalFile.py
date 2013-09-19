#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s

Description:
	2013.08.30 a child class of MatrixFile. used to describe AlignmentDepthIntervalFile in VervetDB
	 
	Example:
	
		
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy, numpy
from pymodule.ProcessOptions import  ProcessOptions
from pymodule import utils, PassingData
from MatrixFile import MatrixFile

class AlignmentDepthIntervalFile(MatrixFile):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(MatrixFile.option_default_dict)
	option_default_dict.update({
						('minDepth', 0, float): [None, '', 1, ''],\
							
						})
	def __init__(self, inputFname=None, **keywords):
		MatrixFile.__init__(self, inputFname=inputFname, **keywords)
		
		#summary data
		self.no_of_intervals = 0
		self.interval_value_ls = []
		self.interval_length_ls = []
		self.chromosome_size = 0
		
		self.min_interval_value = None
		self.max_interval_value = None
		self.median_interval_value = None
		
		self.min_interval_length = None
		self.max_interval_length = None
		self.median_interval_length = None
		
		#
		self.interval_ls = []
		
	def parseRow(self, row):
		"""
		2013.08.30
		"""
		start, stop, length, depth = row[:4]
		start = int(start)
		stop = int(stop)
		length = int(length)
		depth = float(depth)
		return PassingData(start=start, stop=stop, length=length, depth=depth)
	
	def _postProcessParsedRowDataForSummary(self, pdata=None):
		"""
		2013.08.30
		"""
		self.no_of_intervals += 1
		self.chromosome_size += pdata.length
		self.interval_value_ls.append(pdata.depth)
		self.interval_length_ls.append(pdata.length)
		
		if self.min_interval_value is None or pdata.depth <self.min_interval_value:
			self.min_interval_value = pdata.depth
		if self.max_interval_value is None or pdata.depth >self.max_interval_value:
			self.max_interval_value = pdata.depth
	
	def readThroughAndProvideSummary(self):
		"""
		2013.08.30
			called by vervet/src/db/input/AddAlignmentDepthIntervalFile2DB.py
		"""
		
	
		
		col_name2index= self.smartReadHeader()
		if col_name2index is None:
			pdata = self.parseRow(self._row)
			self._postProcessParsedRowDataForSummary(pdata)
		
		for row in self:
			pdata = self.parseRow(row)
			self._postProcessParsedRowDataForSummary(pdata)
		
		
		self.min_interval_length = numpy.min(self.interval_length_ls)
		self.max_interval_length = numpy.max(self.interval_length_ls)
		self.median_interval_length = numpy.median(self.interval_length_ls)
		
		self.mean_interval_value=numpy.mean(self.interval_value_ls)
		self.median_interval_value=numpy.median(self.interval_value_ls)
		return PassingData(no_of_intervals=self.no_of_intervals, chromosome_size=self.chromosome_size, \
						mean_interval_value=self.mean_interval_value,\
						median_interval_value=self.median_interval_value,\
						min_interval_value=self.min_interval_value,\
						max_interval_value=self.max_interval_value,\
						
						min_interval_length=self.min_interval_length,\
						max_interval_length=self.max_interval_length ,\
						median_interval_length=self.median_interval_length)
	
	def _postProcessParsedRowDataForInterval(self, pdata=None, minDepth=None, maxDepth=None):
		"""
		2013.08.30
		"""
		toInclude = True
		if minDepth is not None and pdata.depth<minDepth:
			toInclude = False
		if maxDepth is not  None and pdata.depth>maxDepth:
			toInclude = False
		if toInclude:
			if len(self.interval_ls)>0:
				previousInterval = self.interval_ls[-1]
			else:
				previousInterval = None
			if previousInterval is not None and pdata.start==previousInterval.stop+1:	#adjacent, merge the two
				previousIntervalArea = previousInterval.length*previousInterval.depth
				currentIntervalArea = pdata.length*pdata.depth
				newLength = previousInterval.length + pdata.length
				newDepth = (previousIntervalArea+currentIntervalArea)/newLength
				previousInterval.stop = pdata.stop
				previousInterval.length = newLength
				previousInterval.depth = newDepth
				if self.max_interval_length is None or newLength<self.max_interval_length:
					self.max_interval_length = newLength
			else:
				self.interval_ls.append(pdata)
				if self.min_interval_length is None or pdata.length<self.min_interval_length:
					self.min_interval_length = pdata.length
		
	def getAllIntervalWithinDepthRange(self, minDepth=None, maxDepth=None):
		"""
		2013.08.30
		"""
		col_name2index= self.smartReadHeader()
		if col_name2index is None:
			pdata = self.parseRow(self._row)
			self._postProcessParsedRowDataForInterval(pdata, minDepth=minDepth, maxDepth=maxDepth)
		
		for row in self:
			pdata = self.parseRow(row)
			self._postProcessParsedRowDataForInterval(pdata, minDepth=minDepth, maxDepth=maxDepth)
		return self.interval_ls
		
if __name__ == '__main__':
	main_class = AlignmentDepthIntervalFile
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()