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
from pymodule.yhio.MatrixFile import MatrixFile
from pymodule.yhio.HDF5MatrixFile import HDF5MatrixFile

class AbstractMatrixFileWalker(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	#option_default_dict.update(AbstractMapper.db_option_dict.copy())
	option_default_dict.update({
						('minNoOfTotal', 1, int): [10, 'M', 1, 'minimum no of data from one file for afterFileFunction() to run'],\
						('maxNoOfTotal', 0, int): [None, '', 1, 'maximum no of data to sample from one file. if not set, no limit.'],\
						('samplingRate', 1, float): [1, 's', 1, 'how often you include the data, a probability between 0 and 1.'],\
						('whichColumn', 0, int): [None, 'w', 1, 'data from this column (index starting from 0) is plotted as y-axis value'],\
						('whichColumnHeader', 0, ): [None, 'W', 1, 'column header for the data to be plotted as y-axis value, substitute whichColumn'],\
						('logY', 0, int): [0, '', 1, 'value 0: nothing; 1: log(), 2: -log(). replacing self.logWhichColumn.'],\
						('valueForNonPositiveYValue', 1, float): [-1, 'v', 1, 'default value when log-transformation fails (when value is negative)'],\
						('missingDataNotation', 0, ): ['NA', '', 1, 'coma-separated list of notations for missing data.\n\
	missing data will be skipped.'],\
						('inputFileFormat', 0, int): [1, '', 1, '1: csv-like plain text file; 2: HDF5MatrixFile.'],\
						('outputFileFormat', 0, int): [1, '', 1, '1: csv-like plain text file; 2: HDF5MatrixFile.'],\
						})
	
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		2012.10.25 self.missingDataNotation will be processed into a set.
		"""
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
		#if user wants to preserve data in a data structure that is visible throughout reading different files.
		# then use this self.invariantPData.
		self.invariantPData = PassingData(writer=None, headerOutputted=False, x_ls = [], y_ls = [], z_ls=[])
		if getattr(self, 'missingDataNotation', None):
			self.missingDataNotation = set(utils.getListOutOfStr(self.missingDataNotation, data_type=str, separator2=None))
		else:
			self.missingDataNotation = set()
	
	def connectDB(self):
		"""
			split out of __init__() so that derived classes could overwrite this function
		"""
		pass
	
	def preFileFunction(self, **keywords):
		"""
		2012.11.23
		"""
		if hasattr(self, 'initiatePassingData'):
			return self.initiatePassingData(**keywords)
	
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
	
	def processValue(self, value=None, processType=None, valueForNonPositiveValue=None, **keywords):
		"""
		2012.11.282 change positiveLog to processType
			processType 0: nothing; 1: log(), 2: -log().
		2012.10.15
			the default value of takeLogarithm depends on (self.logWhichColumn or self.logY)
		2012.10.7
		"""
		if valueForNonPositiveValue is None:
			valueForNonPositiveValue = self.valueForNonPositiveYValue
		value = float(value)
		if processType is not None and processType>=1:
			if value>0:
				if processType==1:
					value = math.log10(value)
				elif processType==2:
					value = -math.log10(value)
			else:
				value = valueForNonPositiveValue
		return value
	
	def handleYValue(self, yValue=None, processType=None, valueForNonPositiveValue=None, **keywords):
		"""
		2012.11.28 change positiveLog to processType
		2012.10.7
			add argument takeLogarithm, positiveLog, valueForNonPositiveValue
		2012.8.2
		"""
		return self.processValue(value=yValue, processType=processType, \
						valueForNonPositiveValue=valueForNonPositiveValue, **keywords)
	
	def processHeader(self, header=None, pdata=None, dtypeList=None):
		"""
		2012.11.22
		2012.8.13
			called right after the header of an input file is derived in fileWalker()
		"""
		if not self.invariantPData.headerOutputted:
			if self.outputFileFormat==1:
				if self.invariantPData.writer and header:
					self.invariantPData.writer.writerow(header)
			elif self.outputFileFormat==2:	#for HDF5MatrixFile 
				if not dtypeList and header:	#generate a dtypeList based on header
					dtypeList = []
					for colID in header:
						dtypeList.append((colID, HDF5MatrixFile.varLenStrType))
				#dtypeList = [('locus_id','i8'),('chromosome', HDF5MatrixFile.varLenStrType), ('start','i8'), ('stop', 'i8'), \
				#	('score', 'f8'), ('MAC', 'i8'), ('MAF', 'f8')]
				writer = HDF5MatrixFile(self.outputFname, openMode='w', dtypeList=dtypeList)
				self.invariantPData.writer = writer
		self.invariantPData.headerOutputted = True
	
	
	def processRow(self, row=None, pdata=None):
		"""
		2012.10.15
			return returnValue: 0 if not included or nothing is done on it.
				1 if included or something is carried out on it.
		2012.8.2
			handles each row in each file, here it replaces the yValue
		"""
		returnValue = 0
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
				if yValue not in self.missingDataNotation:
					yValue = self.processValue(yValue, processType=self.logY, valueForNonPositiveValue=self.valueForNonPositiveYValue)
				row[whichColumn] = yValue
			if self.invariantPData.writer:
				self.invariantPData.writer.writerow(row)
				returnValue = 1
		return returnValue
	
	def getNumberOfData(self, pdata):
		"""
		2012.8.6
		"""
		return len(pdata.y_ls)
		
	def afterFileFunction(self, **keywords):
		"""
		2012.10.7
		"""
		if hasattr(self, 'plot'):
			return self.plot(**keywords)
	
	
	def fileWalker(self, inputFname=None, preFileFunction=None, afterFileFunction=None, processRowFunction=None , run_type=1):
		"""
		2012.11.23 added argument preFileFunction
		2012.11.22 inputFileFormat==2, support HDF5MatrixFile
		2012.8.1
		"""
		sys.stderr.write("walking through %s ..."%(inputFname))
		counter = 0
		real_counter = 0
		noOfSampled = 0
		if preFileFunction is None:
			preFileFunction = self.preFileFunction
		pdata = preFileFunction()
		if processRowFunction is None:
			processRowFunction = self.processRow
		if afterFileFunction is None:
			afterFileFunction = self.afterFileFunction
		try:
			if self.inputFileFormat==1:	#2012.11.22 support HDF5MatrixFile
				reader = MatrixFile(inputFname)
			else:
				reader = HDF5MatrixFile(inputFname, openMode='r')
			
			col_name2index = reader.constructColName2IndexFromHeader()
			header = reader.getHeader()
			
			#output the header to the output file if necessary 
			self.processHeader(header=header, pdata=pdata) #2012.8.13
			pdata.reader = reader
			pdata.col_name2index = col_name2index
			pdata.filename = inputFname
			
			for row in reader:
				counter += 1
				if self.samplingRate<1 and self.samplingRate>=0:
					r = random.random()
					if r>self.samplingRate:
						continue
				rowReturnValue = processRowFunction(row=row, pdata=pdata)
				if rowReturnValue is None:	#2012.10.15
					rowReturnValue = 0
				real_counter += rowReturnValue
				noOfSampled += 1
				if self.maxNoOfTotal and real_counter>self.maxNoOfTotal:
					break
			if self.getNumberOfData(pdata)>=self.minNoOfTotal:
				afterFileFunction(x_ls=pdata.x_ls, y_ls=pdata.y_ls, pdata=pdata)
			del reader
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
			sys.exit(3)
		if counter>0:
			fraction = float(noOfSampled)/float(counter)
		else:
			fraction = 0
		sys.stderr.write("%s/%s (%.3f) data sampled. real_counter=%s.\n"%(noOfSampled, counter, fraction, real_counter))
	
	def setup(self, **keywords):
		"""
		2012.11.22
		2012.10.25
			do not open the file if it's a png file
		2012.10.15
			run before anything is run
		"""
		writer = None
		if self.outputFileFormat==1:
			suffix = os.path.splitext(self.outputFname)[1]
			if self.outputFname and suffix!='.png':
				writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		else:	#HDF5MatrixFile
			#can't generate HDF5MatrixFile, because it needs dtypeList
			pass
		#pass it to the invariantPData
		self.invariantPData.writer = writer
	
	def reduce(self, **keywords):
		"""
		2012.10.15
			run after all files have been walked through
		"""
		if getattr(self.invariantPData, 'writer', None):
			self.invariantPData.writer.close()
			del self.invariantPData.writer
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		self.setup()
		
		for inputFname in self.inputFnameLs:
			if os.path.isfile(inputFname):
				self.fileWalker(inputFname, afterFileFunction=self.afterFileFunction, run_type=1, processRowFunction=self.processRow)
		
		self.reduce()
			
if __name__ == '__main__':
	main_class = AbstractMatrixFileWalker
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()