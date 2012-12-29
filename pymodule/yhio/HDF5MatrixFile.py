#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s

Description:
	2012.11.15 a csv file format, stored in HDF5.
		This HDF5 file is composed of groups. At least one group, group0.
	i.e.
		reader = HDF5MatrixFile(inputFname=filename, openMode='r')
		reader = HDF5MatrixFile(filename, openMode='r')
		for row in reader:
			...
		tableObject = reader.getTableObject(tableName=tableName)
		for row in tableObject:
			...
		
		rowDefinition = [('locus_id','i8'),('chr', HDF5MatrixFile.varLenStrType), ('start','i8'), ('stop', 'i8'), \
					('score', 'f8'), ('MAC', 'i8'), ('MAF', 'f8')]
		headerList = [row[0] for row in rowDefinition]
		dtype = numpy.dtype(rowDefinition)
		
		writer = HDF5MatrixFile(inputFname=filename, openMode='w', dtype=dtype)
		writer = HDF5MatrixFile(filename, openMode='w', dtype=dtype)
		
		if writer:
			tableObject = writer.createNewTable(tableName=tableName, dtype=dtype)
			tableObject.setColIDList(headerList)
		elif outputFname:
			writer = HDF5MatrixFile(outputFname, openMode='w', dtype=dtype, tableName=tableName)
			writer.writeHeader(headerList)
			tableObject = writer.getTableObject(tableName=tableName)
		cellList = []
		for data_obj in self.data_obj_ls:
			dataTuple = self._extractOutputRowFromDataObject(data_obj=data_obj)
			cellList.append(dataTuple)
		tableObject.writeCellList(cellList)
		if closeFile:
			writer.close()
			del writer
		
		
		#each number below is counting bytes, not bits
		rowDefinition = [('locus_id','i8'),('chromosome', HDF5MatrixFile.varLenStrType), ('start','i8'), ('stop', 'i8'), \
					('score', 'f8'), ('MAC', 'i8'), ('MAF', 'f8')]
		if writer is None and filename:
			writer = HDF5MatrixFile(filename, openMode='w', rowDefinition=rowDefinition, tableName=tableName)
			tableObject = writer.getTableObject(tableName=tableName)
		elif writer:
			tableObject = writer.createNewTable(tableName=tableName, rowDefinition=rowDefinition)
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import csv
import h5py
import numpy
from pymodule.utils import PassingData, PassingDataList
from pymodule.ProcessOptions import  ProcessOptions
from pymodule.yhio.MatrixFile import MatrixFile
varLenStrType = h5py.new_vlen(str)

class YHTableInHDF5Group(object):
	option_default_dict = {
							('h5Group', 0, ): [None, '', 1, "the h5py group ojbect."],\
							('newGroup', 1, int): [0, '', 1, "whether this is a new group or an existing group in a file"],\
							('dataMatrixDtype', 0, ): ['f', '', 1, 'data type in the dataMatrix. candidates are i, f8, compound type, etc.'],\
							('compression', 0, ): [None, '', 1, 'the compression engine for all underlying datasets'],\
							('compression_opts', 0, ): [None, '', 1, 'option for the compression engine, gzip level, or tuple for szip'],\
							
							}
	def __init__(self, **keywords):
		"""
		dataMatrixDtype could be a compound type:
			http://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html
			http://docs.scipy.org/doc/numpy/reference/generated/numpy.dtype.html
				
				#A record data type containing a 16-character string (in field name)
					#and a sub-array of two 64-bit floating-point number (in field grades):
				dt = numpy.dtype([('name', numpy.str_, 16), ('grades', numpy.float64, (2,))])
				
				my_dtype = numpy.dtype([('field1', 'i'), ('field2', 'f'), ('field3', varLenStrType)])
				
				#Using array-protocol type strings:
				#each number below is counting bytes, not bits
				>>> numpy.dtype([('a','f8'),('b','S10')])
				dtype([('a', '<f8'), ('b', '|S10')])
				
				#Using tuples. int is a fixed type, 3 the field's shape. void is a flexible type, here of size 10:
				numpy.dtype([('hello',(numpy.int,3)),('world',numpy.void,10)])
				
				#Using dictionaries. Two fields named 'gender' and 'age':
				numpy.dtype({'names':['gender','age'], 'formats':['S1',numpy.uint8]})
				
				#Offsets in bytes, here 0 and 25:
				numpy.dtype({'surname':('S25',0),'age':(numpy.uint8,25)})
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		
		self.dataMatrixDSName = "dataMatrix"
		self.rowIDListDSName = "rowIDList"
		self.colIDListDSName = "colIDList"
		if not self.newGroup:
			self._readInData()
		else:
			self._createDatasetSkeletonForOneGroup(h5Group=self.h5Group, dtype=self.dataMatrixDtype)
		
		self.newWrite = True	#a flag used to control whether it's first time to write stuff (first time=set whole matrix)
		
		self.rowIndexCursor = 0
		
	def _createDatasetSkeletonForOneGroup(self, h5Group=None, dtype=None):
		"""
		2012.11.18
		"""
		if h5Group is None:
			h5Group = self.h5Group
		if dtype is None:
			dtype = self.dataMatrixDtype
		self.dataMatrix = h5Group.create_dataset(self.dataMatrixDSName, shape=(1,), dtype=dtype, \
										maxshape=(None, ), compression=self.compression, compression_opts=self.compression_opts)
		self.dataMatrix.resize((0,))	#by default it contains one "null" data point.
		self.rowIDList = h5Group.create_dataset(self.rowIDListDSName, shape=(1,), dtype=varLenStrType, maxshape=(None,),\
											compression=self.compression, compression_opts=self.compression_opts)
		self.rowIDList.resize((0,))
		self.colIDList = h5Group.create_dataset(self.colIDListDSName, shape=(1,), dtype=varLenStrType, maxshape=(None,),\
											compression=self.compression, compression_opts=self.compression_opts)
		self.colIDList.resize((0,))
		self.rowID2rowIndex = {}
		self.colID2colIndex = {}
	
	def _readInData(self):
		"""
		2012.11.16
		"""
		self.dataMatrix = self.h5Group[self.dataMatrixDSName]
		self.rowIDList = self.h5Group[self.rowIDListDSName]
		self.colIDList = self.h5Group[self.colIDListDSName]
		self._processRowIDColID()
		#pass the HDF5Group attributes to this object itself, it ran into "can't set attribute error".
		# conflict with existing property
		#for attributeName, attributeValue in self.h5Group.attrs.iteritems():
		#	object.__setattr__(self, attributeName, attributeValue)
	
	def _processRowIDColID(self):
		"""
		2012.11.16 similar to SNPData.processRowIDColID()
		"""
			
		rowIDList = self.rowIDList
		colIDList = self.colIDList
		rowID2rowIndex = {}
		if rowIDList:	#2008-12-03
			for i in range(len(rowIDList)):
				rowID = rowIDList[i]
				rowID2rowIndex[rowID] = i
		
		
		colID2colIndex = {}
		if colIDList:	#2008-12-03
			for i in range(len(colIDList)):
				colID = colIDList[i]
				colID2colIndex[colID] = i
		
		self.rowID2rowIndex = rowID2rowIndex
		self.colID2colIndex = colID2colIndex
	
	def extendDataMatrix(self, dataMatrix=None, **keywords):
		"""
		2012.11.16
			dataMatrix has to be 2D
		"""
		if dataMatrix:
			m = len(dataMatrix)
			s = self.dataMatrix.shape[0]
			if self.newWrite:	#defaultData in self.dataMatrix is of shape (1,)
				self.dataMatrix.resize((m,))
				self.dataMatrix[0:m] = dataMatrix
			else:
				self.dataMatrix.resize((s+m,))
				self.dataMatrix[s:s+m] = dataMatrix
			self.newWrite = False
	
	def writeCellList(self, cellList=None, **keywords):
		"""
		2012.11.19
			call self.extendDataMatrix()
		"""
		self.extendDataMatrix(dataMatrix=cellList, **keywords)
	
	def _appendID(self, idDataset=None, idValue=None):
		"""
		2012.11.8
		"""
		self._extendIDList(idDataset=idDataset, idList=[idValue])

	def _extendIDList(self, idDataset=None, idList=None):
		"""
		2012.11.16
		"""
		if idDataset is None:
			idDataset = self.colIDList
		
		if idList:
			m = len(idList)
			s = idDataset.shape[0]
			idDataset.resize((s+m,))
			idDataset[s:s+m] = idList

	def _setIDList(self, idDataset=None, idList=None):
		"""
		2012.11.16
		"""
		if idDataset is None:
			idDataset = self.colIDList
		
		if idList:
			m = len(idList)
			s = idDataset.shape[0]
			idDataset.resize((m,))
			idDataset[:] = idList
		
	def setColIDList(self, colIDList=None):
		"""
		"""
		self._setIDList(idDataset=self.colIDList, idList=colIDList)
		
	def setRowIDList(self, rowIDList=None):
		"""
		"""
		self._setIDList(idDataset=self.rowIDList, idList=rowIDList)
	
	def getColIndex(self, colID=None):
		"""
		
		"""
		colID2colIndex = self.colID2colIndex
		colIndex = None
		if colID2colIndex:
			colIndex = colID2colIndex.get(colID, None)
		return colIndex
	
	def getRowIndex(self, rowID=None):
		"""
		"""
		rowID2rowIndex = self.rowID2rowIndex
		rowIndex = None
		if rowID2rowIndex:
			rowIndex = rowID2rowIndex.get(rowID, None)
		return rowIndex
	
	def getCellDataGivenRowColID(self, rowID=None, colID=None):
		"""
		"""
		
		rowIndex = self.getRowIndex(rowID)
		colIndex = self.getColIndex(colID)
		
		cellData = None
		if rowIndex is not None and colIndex is not None:
			cellData = self.dataMatrix[rowIndex][colIndex]
		return cellData
	
	@property
	def name(self):
		return self.h5Group.name
	
	def __iter__(self):
		return self
	
	def next(self):
		"""
		2012.11.19 part of faking as a file object
		"""
		pdata = PassingDataList()
		if self.rowIndexCursor<self.dataMatrix.shape[0]:
			row = self.dataMatrix[self.rowIndexCursor]
			for colID in self.colIDList:	#iteration over colIDList is in the same order as the ascending order of colIndex
				#but iteration over self.colID2colIndex is not in the same order as the ascending order of colIndex
				colIndex = self.colID2colIndex.get(colID)
				setattr(pdata, colID, row[colIndex])
			self.rowIndexCursor += 1
		else:
			raise StopIteration
		return pdata
	
	def reset(self):
		"""
		2012.11.19 part of faking as a file object
		"""
		self.rowIndexCursor = 0
	
	def constructColName2IndexFromHeader(self):
		"""
		2012.11.22 added so that YHSingleTableFile could be used as HDF5MatrixFile
		"""
		return self.colID2colIndex
	
	
	def getColIndexGivenColHeader(self, colHeader=None):
		"""
		2012.11.15
			this is from the combined column header list.
		"""
		return self.colID2colIndex.get(colHeader)
	
	def addAttribute(self, name=None, value=None, overwrite=True):
		"""
		
		"""
		if name in self.h5Group.attrs:
			sys.stderr.write("Warning: h5Group %s already has attribute %s=%s.\n"%(self.name, name, value))
			if overwrite:
				self.h5Group.attrs[name] = value
			else:
				return False
		else:
			self.h5Group.attrs[name] = value
		#pass the HDF5Group attributes to this object itself , it ran into "can't set attribute error". conflict with existing property
		#object.__setattr__(self, name, value)
		#setattr(self, name, value)
		return True
	
	def addAttributeDict(self, attributeDict=None):
		"""
		2012.12.18
		"""
		addAttributeDictToYHTableInHDF5Group(tableObject=self, attributeDict=attributeDict)
	
	def getAttribute(self, name=None, defaultValue=None):
		return self.h5Group.attrs.get(name, defaultValue)
	
	def getAttributes(self):
		return self.h5Group.attrs
	
	def getListAttributeInStr(self, name=None):
		"""
		2012.11.22
			this attribute must be a list or array
		"""
		attr_in_str = ''
		attributeValue = self.getAttribute(name=name)
		if attributeValue is not None and type(attributeValue)==numpy.ndarray:
			if hasattr(attributeValue, '__len__') and attributeValue.size>0:
				ls = map(str, attributeValue)
				attr_in_str = ','.join(ls)
		return attr_in_str
	
	def addObjectAttributeToSet(self, attributeName=None, setVariable=None):
		"""
		2012.12.3
		2012.11.22
			do not add an attribute to the set if it's not available or if it's none
		"""
		attributeValue = self.getAttribute(attributeName, None)
		if attributeValue is not None and setVariable is not None:
			setVariable.add(attributeValue)
		return setVariable
	
	def addObjectListAttributeToSet(self, attributeName=None, setVariable=None):
		"""
		2012.12.3
		2012.12.2 bugfix
		2012.11.23
		"""
		attributeValue = self.getAttribute(attributeName, None)
		flag = False
		if type(attributeValue)==numpy.ndarray:	#"if attributeValue" fails for numpy array
			if hasattr(attributeValue, '__len__') and attributeValue.size>0:
				flag = True
		elif attributeValue or attributeValue == 0:
			flag = True
		if flag and setVariable is not None:
			if type(attributeValue)==str:
				attributeValueList = getListOutOfStr(attributeValue, data_type=data_type, separator1=',', separator2='-')
			else:
				attributeValueList = attributeValue
			setVariable |= set(list(attributeValueList))
		return setVariable
	
	def close(self):
		"""
		2012.12.4
			deprecated. dimension set to 0 upon creation.
			resize the data matrix to zero if it's in write mode and self.newWrite is still True
		"""
		if self.newGroup and self.newWrite:	#in write mode and no writing has happened.
			self.dataMatrix.resize((0,))
			#or you can delete it
			#del self.h5Group[self.dataMatrixDSName]
			#self.h5Group.delete_dataset( ...
	

class HDF5MatrixFile(MatrixFile):
	varLenStrType = varLenStrType	#convenient for outside program to access this variable length string type
	__doc__ = __doc__
	option_default_dict = MatrixFile.option_default_dict.copy()
	option_default_dict.update({
							('tableName', 0, ): [None, '', 1, "name for the first table, default is $tableNamePrefix\0. "],\
							('tableNamePrefix', 0, ): ['table', '', 1, "prefix for all table's names"],\
							('rowDefinition', 0, ): [None, '', 1, "data type list for a compound dtype. It overwrites dtype. i.e. a list of i.e. ('start','i8')"],\
							('dtype', 0, ): [None, '', 1, 'data type in the first table to be created. candidates are i, f8, etc.'],\
							('compression', 0, ): ['gzip', '', 1, 'the compression engine for all underlying datasets, gzip, szip, lzf'],\
							('compression_opts', 0, ): [None, '', 1, 'option for the compression engine, gzip level, or tuple for szip'],\
							})
	def __init__(self, inputFname=None, **keywords):
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		if not self.inputFname:
			self.inputFname = inputFname
		
		self.header = None
		self.combinedColIDList = None	#same as header
		self.combinedColID2ColIndex = None
		
		self.hdf5File = h5py.File(self.inputFname, self.openMode)
		self.tableObjectList = []
		self.tablePath2Index = {}
		
		if self.openMode=='r':
			self._readInData()
		elif self.openMode=='w':
			self.createNewTable(tableName=self.tableName, dtype=self.dtype, rowDefinition=self.rowDefinition)
		
		self.rowIndexCursor = 0	#2012.11.16 for iteration
	
	def createNewTable(self, tableName=None, dtype=None, rowDefinition=None):
		"""
		2012.11.20 add argument rowDefinition
		2012.11.19
		"""
		colIDList = None
		if rowDefinition:
			colIDList = [row[0] for row in rowDefinition]
			dtype = numpy.dtype(rowDefinition)
		if dtype is None:
			dtype = self.dtype
		if tableName is None:
			tableName = self._getNewTableName()
		tableObject = YHTableInHDF5Group(h5Group=self.hdf5File.create_group(tableName), \
								newGroup=True, dataMatrixDtype=dtype, \
								compression=self.compression, compression_opts=self.compression_opts)
		if colIDList:
			tableObject.setColIDList(colIDList)
		self._appendNewTable(tableObject)
		return tableObject
	
	def _appendNewTable(self, tableObject=None):
		if tableObject:
			tableName = tableObject.name	##the tableName preceds with "/"
			if tableName in self.tablePath2Index:	 
				sys.stderr.write("ERROR, table %s already in self.tablePath2Index, index=%s.\n"%(tableName,\
																			self.tablePath2Index.get(tableName)))
				sys.exit(3)
			self.tablePath2Index[tableName] = len(self.tablePath2Index)
			self.tableObjectList.append(tableObject)
	
	def _getNewTableName(self, tableNamePrefix=None):
		"""
		2012.11.19
		"""
		if not tableNamePrefix:
			tableNamePrefix = self.tableNamePrefix
		i = len(self.tableObjectList)
		tableName = "%s%s"%(tableNamePrefix, i)
		while "/"+tableName in self.tablePath2Index:	#stop until a unique name shows up
			i += 1
			tableName = "%s%s"%(tableNamePrefix, i)
		return tableName
	
	def _readInData(self):
		"""
		2012.11.16
		"""
		for tableName, h5Group in self.hdf5File.iteritems():
			tableObject = YHTableInHDF5Group(h5Group=h5Group, newGroup=False,\
										dataMatrixDtype=self.dtype)
			self._appendNewTable(tableObject)
		self._setupCombinedColumnIDMapper()
		
	def _setupCombinedColumnIDMapper(self,):
		"""
		2012.11.16
		"""
		self.combinedColID2ColIndex = {}
		self.header = []
		self.combinedColIDList = self.header
		for tableObject in self.tableObjectList:
			colIDList = tableObject.colIDList
			for colID in colIDList:
				if colID in self.combinedColID2ColIndex:
					sys.stderr.write("Error: column ID %s already used in column %s.\n"%(colID, self.combinedColID2ColIndex.get(colID)))
					sys.exit(3)
				else:
					self.combinedColID2ColIndex[colID] = len(self.combinedColID2ColIndex)
					self.header.append(colID)
	
	def constructColName2IndexFromHeader(self):
		"""
		2012.11.22 overwrite parent function
		"""
		return self.combinedColID2ColIndex
	
	def getColIndexGivenColHeader(self, colHeader=None):
		"""
		2012.11.15
			this is from the combined column header list.
		"""
		return self.combinedColID2ColIndex.get(colHeader)
	
	def getTableObject(self, tableIndex=None, tableName=None):
		"""
		"""
		tableObject = None
		if tableIndex is not None:
			tableObject = self.tableObjectList[tableIndex]
		elif tableName:
			if tableName[0]!='/':	#bugfix, add root in front if not there.
				tableName = '/'+tableName
			tableIndex = self.tablePath2Index.get(tableName)
			if tableIndex is not None:
				tableObject = self.tableObjectList[tableIndex]
		else:	#return first table
			tableObject = self.tableObjectList[0]
		return tableObject
	
	def getColIndex(self, colID=None, tableIndex=None, tableName=None):
		"""
		
		"""
		tableObject = self.getTableObject(tableIndex=tableIndex, tableName=tableName)
		return tableObject.getColIndex(colID=colID)
	
	def getRowIndex(self, rowID=None, tableIndex=None, tableName=None):
		"""
		"""
		tableObject = self.getTableObject(tableIndex=tableIndex, tableName=tableName)
		return tableObject.getRowIndex(rowID=rowID)
	
	def getCellDataGivenRowColID(self, rowID=None, colID=None, tableIndex=None, tableName=None):
		"""
		"""
		tableObject = self.getTableObject(tableIndex=tableIndex, tableName=tableName)
		return tableObject.getCellDataGivenRowColID(rowID=rowID, colID=colID)
	
	def __iter__(self):
		return self
	
	def next(self):
		"""
		2012.11.16
			combine the numeric and str row (if both exist). and return the combined row
		"""
		
		row = None
		pdata = PassingDataList()
		for tableObject in self.tableObjectList:
			if self.rowIndexCursor<tableObject.dataMatrix.shape[0]:
				if row is None:
					row = list(tableObject.dataMatrix[self.rowIndexCursor])
				else:
					rowAppend = tableObject.dataMatrix[self.rowIndexCursor]
					row += list(rowAppend)
			else:
				raise StopIteration
			for colID  in self.header:	#iteration over header is in the same order as the ascending order of colIndex
				#but iteration over self.colID2colIndex is not in the same order as the ascending order of colIndex
				colIndex = self.combinedColID2ColIndex.get(colID)
				if colIndex < len(row):
					setattr(pdata, colID, row[colIndex])
		self.rowIndexCursor += 1
		return pdata
	
	def close(self):
		"""
		#2012.12.4 close the table objects first
		"""
		self.hdf5File.close()
		del self.hdf5File
	
	def writeHeader(self, headerList=None, tableIndex=None, tableName=None):
		"""
		2012.11.16
			only the first table
		"""
		self.setColIDList(colIDList=headerList, tableIndex=tableIndex, tableName=tableName)
	
	def writeOneCell(self, oneCell=None, tableIndex=None, tableName=None):
		"""
		2012.11.19
			mimic csv's writerow()
			each cell must be a tuple (readable only). list is not acceptable. 
			it's not very efficient as it's resizing the dataMatrix all the time.
		"""
		tableObject = self.getTableObject(tableIndex=tableIndex, tableName=tableName)
		tableObject.extendDataMatrix(dataMatrix=[oneCell])

	def writeCellList(self, cellList=None, tableIndex=None, tableName=None):
		"""
		2012.11.18
			for bulk writing. more efficient.
			each cell must be a tuple (readable only). list is not acceptable. 
		"""
		tableObject = self.getTableObject(tableIndex=tableIndex, tableName=tableName)
		tableObject.extendDataMatrix(dataMatrix=cellList)
	
	
	def setColIDList(self, colIDList=None, tableIndex=None, tableName=None):
		"""
		"""
		tableObject = self.getTableObject(tableIndex=tableIndex, tableName=tableName)
		tableObject.setColIDList(colIDList=colIDList)
		
	def setRowIDList(self, rowIDList=None, tableIndex=None, tableName=None):
		"""
		"""
		tableObject = self.getTableObject(tableIndex=tableIndex, tableName=tableName)
		tableObject.setRowIDList(rowIDList=rowIDList)
	
	def addAttribute(self, name=None, value=None, overwrite=True, tableIndex=None, tableName=None):
		"""
		2012.11.28 find the tableObject and let it do the job
		"""
		tableObject = self.getTableObject(tableIndex=tableIndex, tableName=tableName)
		return tableObject.addAttribute(name=name, value=value, overwrite=overwrite)
	
	def addAttributeDict(self, attributeDict=None, tableObject=None):
		if tableObject is None:
			tableObject = self.getTableObject()
		
		addAttributeDictToYHTableInHDF5Group(tableObject=tableObject, attributeDict=attributeDict)
	
	def getAttribute(self, name=None, defaultValue=None, tableIndex=None, tableName=None):
		"""
		2012.11.28
		"""
		tableObject = self.getTableObject(tableIndex=tableIndex, tableName=tableName)
		return tableObject.getAttribute(name=name, defaultValue=defaultValue)
	
	def getAttributes(self, tableIndex=None, tableName=None):
		"""
		2012.11.28
		"""
		tableObject = self.getTableObject(tableIndex=tableIndex, tableName=tableName)
		return tableObject.getAttributes()
	
	def getListAttributeInStr(self, name=None, tableIndex=None, tableName=None):
		"""
		2012.11.22
			this attribute must be a list or array
		"""
		tableObject = self.getTableObject(tableIndex=tableIndex, tableName=tableName)
		return tableObject.getListAttributeInStr(name=name)
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()

def addAttributeDictToYHTableInHDF5Group(tableObject=None, attributeDict=None):
	"""
	2012.11.22 convenient function
		attributeValue could not be high-level python objects, such as list, set.
		numpy.array could replace list.
	"""
	if tableObject is not None and attributeDict is not None:
		for attributeName, attributeValue in attributeDict.iteritems():
			doItOrNot = False
			if type(attributeValue)==numpy.ndarray:
				if hasattr(attributeValue, '__len__') and attributeValue.size>0:
					doItOrNot = True
			elif attributeValue or attributeValue == 0:	#empty array will be ignored but not 0
				doItOrNot = True
			if doItOrNot:
				tableObject.addAttribute(name=attributeName, value=attributeValue)
	

if __name__ == '__main__':
	main_class = HDF5MatrixFile
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
