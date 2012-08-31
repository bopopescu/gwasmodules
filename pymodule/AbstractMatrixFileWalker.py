#!/usr/bin/env python
"""
Examples:
	%s
	
Description:
	2012.8.15
		abstract class for programs that walk through a list of matrix-like files.
		Running it will reformat input into tab-delimited tsv matrix file.
		If whichColumn/whichColumnHeader is given, it'll convert it into float or log or -log.
		It combines both input from "-i oneFile.txt" or the trailing standalone arguments.
		It could act as both a mapper (one input) or reducer (multiple input).

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv, random
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, getColName2IndexFromHeader, figureOutDelimiter
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper

class AbstractMatrixFileWalker(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	#option_default_dict.update(AbstractMapper.db_option_dict.copy())
	option_default_dict.update({
						('inputFname', 0, ):option_default_dict[('inputFname', 1, )],\
						('minNoOfTotal', 1, int): [10, 'M', 1, 'minimum no of data from one file for afterFileFunction() to run'],\
						('samplingRate', 1, float): [1, 's', 1, 'how often you include the data, a probability between 0 and 1.'],\
						('whichColumn', 0, int): [None, 'w', 1, 'data from this column (index starting from 0) is plotted as y-axis value'],\
						('whichColumnHeader', 0, ): [None, 'W', 1, 'column header for the data to be plotted as y-axis value, substitute whichColumn'],\
						('logWhichColumn', 0, int): [0, 'g', 0, 'whether to take -log of the whichColumn'],\
						('positiveLog', 0, int): [0, 'p', 0, 'toggle to take log, rather than -log(), \
				only effective when logWhichColumn is toggled. '],\
						('valueForNonPositiveYValue', 1, float): [50, 'v', 1, 'if the whichColumn value is not postive and logWhichColumn is on,\
			what yValue should be.'],\
						('missingDataNotation', 0, ): ['NA', '', 1, 'notation for missing data. will be skipped.'],\
								})
	#pop in the end after its value is used above
	option_default_dict.pop(('inputFname', 1, ))
	
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
		#if user wants to preserve data in a data structure that is visible throughout reading different files.
		# then use this self.invariantPData.
		self.invariantPData = PassingData(writer=None, headerOutputted=False, x_ls = [], y_ls = [])
		
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
		#2012.8.16 pass to global data
		self.invariantPData.y_ls = pdata.y_ls
		self.invariantPData.x_ls = pdata.x_ls
		return pdata
	
	
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
			handles each row in each file, here it replaces the yValue
		"""
		col_name2index = getattr(pdata, 'col_name2index', None)
		y_ls = getattr(pdata, 'y_ls', None)
		if col_name2index and y_ls is not None:
			if self.whichColumnHeader:
				whichColumn = col_name2index.get(self.whichColumnHeader, None)
			elif self.whichColumn:
				whichColumn = self.whichColumn
			else:
				whichColumn = None
			if whichColumn is not None:
				yValue = row[whichColumn]
				if yValue!=self.missingDataNotation:
					yValue = self.handleYValue(yValue)
				row[whichColumn] = yValue
			if self.invariantPData.writer:
				self.invariantPData.writer.writerow(row)
	
	def afterFileFunction(self, pdata=None, **keywords):
		"""
		2012.8.15
			executed after the file is read through
		"""
		pass
	
	def processHeader(self, header=None, pdata=None):
		"""
		2012.8.13
			called everytime the header of an input file is derived in fileWalker()
		"""
		if self.invariantPData.writer and not self.invariantPData.headerOutputted:
			self.invariantPData.writer.writerow(header)
			self.invariantPData.headerOutputted = True
	
	
	def getNumberOfData(self, pdata):
		"""
		2012.8.6
		"""
		return len(pdata.y_ls)
	
	def fileWalker(self, inputFname=None, afterFileFunction=None, processRowFunction=None , run_type=1):
		"""
		2012.8.1
		"""
		sys.stderr.write("walking through %s ..."%(inputFname))
		counter =0
		pdata = self.initiatePassingData()
		if processRowFunction is None:
			processRowFunction = self.processRow
		if afterFileFunction is None:
			afterFileFunction = self.plot
		try:
			inf = utils.openGzipFile(inputFname)
			delimiter = figureOutDelimiter(inf)
			isCSVReader = True
			if delimiter=='\t' or delimiter==',':
				reader = csv.reader(inf, delimiter=delimiter)
				header = reader.next()
			else:
				reader = inf
				header = inf.readline().strip().split()	#whatever splits them
				isCSVReader = False
			self.processHeader(header=header, pdata=pdata) #2012.8.13
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
				afterFileFunction(x_ls=pdata.x_ls, y_ls=pdata.y_ls, pdata=pdata)
			del reader
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
			sys.exit(3)
		sys.stderr.write("%s data.\n"%(counter))
	
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		if self.outputFname:
			writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		else:
			writer = None
		#pass it to the invariantPData
		self.invariantPData.writer = writer
		
		for inputFname in self.inputFnameLs:
			if os.path.isfile(inputFname):
				self.fileWalker(inputFname, afterFileFunction=self.afterFileFunction, run_type=1, processRowFunction=self.processRow)
		if writer:
			del writer
		
			
if __name__ == '__main__':
	main_class = AbstractMatrixFileWalker
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()