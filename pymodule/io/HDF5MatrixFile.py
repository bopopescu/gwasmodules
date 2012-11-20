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
		groupObject = reader.getGroupObject(groupName=groupName)
		for row in groupObject:
			...
		
		dtypeList = [('locus_id','i8'),('chr', HDF5MatrixFile.varLenStrType), ('start','i8'), ('stop', 'i8'), \
					('score', 'f8'), ('MAC', 'i8'), ('MAF', 'f8')]
		headerList = [row[0] for row in dtypeList]
		dtype = numpy.dtype(dtypeList)
		
		writer = HDF5MatrixFile(inputFname=filename, openMode='w', dtype=dtype)
		writer = HDF5MatrixFile(filename, openMode='w', dtype=dtype)
		
		if writer:
			groupObject = writer.createNewGroup(groupName=groupName, dtype=dtype)
			groupObject.setColIDList(headerList)
		elif outputFname:
			writer = HDF5MatrixFile(outputFname, openMode='w', dtype=dtype, firstGroupName=groupName)
			writer.writeHeader(headerList)
			groupObject = writer.getGroupObject(groupName=groupName)
		cellList = []
		for data_obj in self.data_obj_ls:
			dataTuple = self._extractOutputRowFromDataObject(data_obj=data_obj)
			cellList.append(dataTuple)
		groupObject.writeCellList(cellList)
		if closeFile:
			writer.close()
			del writer
		
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import utils, figureOutDelimiter
from pymodule.ProcessOptions import  ProcessOptions
from pymodule.io.MatrixFile import MatrixFile
import csv
import h5py
import numpy
varLenStrType = h5py.new_vlen(str)

class HDF5GroupWrapper(object):
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
		self.rowIDList = h5Group.create_dataset(self.rowIDListDSName, shape=(1,), dtype=varLenStrType, maxshape=(None,),\
											compression=self.compression, compression_opts=self.compression_opts)
		self.colIDList = h5Group.create_dataset(self.colIDListDSName, shape=(1,), dtype=varLenStrType, maxshape=(None,),\
											compression=self.compression, compression_opts=self.compression_opts)
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
	
	def extendDataMatrix(self, dataMatrix=None):
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
	
	def writeCellList(self, cellList):
		"""
		2012.11.19
			call self.extendDataMatrix()
		"""
		self.extendDataMatrix(dataMatrix=cellList)
	
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
		
		if self.rowIndexCursor<self.dataMatrix.shape[0]:
			row = self.dataMatrix[self.rowIndexCursor]
			self.rowIndexCursor += 1
		else:
			raise StopIteration
		return row
	
	def reset(self):
		"""
		2012.11.19 part of faking as a file object
		"""
		self.rowIndexCursor = 0
	
	def getColIndexGivenColHeader(self, colHeader=None):
		"""
		2012.11.15
			this is from the combined column header list.
		"""
		return self.colID2colIndex.get(colHeader)
	
	def addAttribute(self, name=None, value=None, overwrite=False):
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
		return True
	
	def getAttribute(self, name=None):
		return self.h5Group.attrs.get(name)
	
	def getAttributes(self):
		return self.h5Group.attrs

class HDF5MatrixFile(MatrixFile):
	varLenStrType = varLenStrType	#convenient for outside program to access this variable length string type
	__doc__ = __doc__
	option_default_dict = MatrixFile.option_default_dict.copy()
	option_default_dict.update({
							('firstGroupName', 0, ): [None, '', 1, "name for the first group, default is $groupNamePrefix\0."],\
							('groupNamePrefix', 0, ): ['group', '', 1, "prefix for all group's names"],\
							('dtype', 0, ): ['f', '', 1, 'data type in the first group to be created. candidates are i, f8, etc.'],\
							('compression', 0, ): ['gzip', '', 1, 'the compression engine for all underlying datasets, gzip, szip, lzf'],\
							('compression_opts', 0, ): [None, '', 1, 'option for the compression engine, gzip level, or tuple for szip'],\
							})
	def __init__(self, inputFname=None, **keywords):
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		if not self.inputFname:
			self.inputFname = inputFname
		
		self.combinedColID2ColIndex = None
		
		self.hdf5File = h5py.File(self.inputFname, self.openMode)
		self.h5GroupList = []
		self.groupName2Index = {}
		
		if self.openMode=='r':
			self._readInData()
		elif self.openMode=='w':
			self.createNewGroup(groupName=self.firstGroupName, dtype=self.dtype)
		
		self.rowIndexCursor = 0	#2012.11.16 for iteration
	
	def createNewGroup(self, groupName=None, dtype=None):
		"""
		2012.11.19
		"""
		if dtype is None:
			dtype = self.dtype
		if groupName is None:
			groupName = self.getNewGroupName()
		elif groupName in self.groupName2Index:	#use itself as prefix
			groupName = self.getNewGroupName(groupNamePrefix=groupName)
		groupObject = HDF5GroupWrapper(h5Group=self.hdf5File.create_group(groupName), \
								newGroup=True, dataMatrixDtype=dtype, \
								compression=self.compression, compression_opts=self.compression_opts)
		self._appendNewGroup(groupObject)
		return groupObject
	
	def _appendNewGroup(self, groupObject=None):
		if groupObject:
			groupName = groupObject.name
			if groupName in self.groupName2Index:
				sys.stderr.write("ERROR, group %s already in self.groupName2Index, index=%s.\n"%(groupName,\
																			self.groupName2Index.get(groupName)))
				sys.exit(3)
			self.groupName2Index[groupName] = len(self.groupName2Index)
			self.h5GroupList.append(groupObject)
	
	def getNewGroupName(self, groupNamePrefix=None):
		"""
		2012.11.19
		"""
		if not groupNamePrefix:
			groupNamePrefix = self.groupNamePrefix
		i = len(self.h5GroupList)
		groupName = "%s%s"%(groupNamePrefix, i)
		while groupName in self.groupName2Index:	#stop until a unique name shows up
			i += 1
			groupName = "%s%s"%(groupNamePrefix, i)
		return groupName
	
	def _readInData(self):
		"""
		2012.11.16
		"""
		for groupName, h5Group in self.hdf5File.iteritems():
			groupObject = HDF5GroupWrapper(h5Group=h5Group, newGroup=False,\
										dataMatrixDtype=self.dtype)
			self._appendNewGroup(groupObject)
		self._setupCombinedColumnIDMapper()
		
	def _setupCombinedColumnIDMapper(self,):
		"""
		2012.11.16
		"""
		self.combinedColID2ColIndex = {}
		for groupObject in self.h5GroupList:
			colIDList = groupObject.colIDList
			for colID in colIDList:
				if colID in self.combinedColID2ColIndex:
					sys.stderr.write("Error: column ID %s already used in column %s.\n"%(colID, self.combinedColID2ColIndex.get(colID)))
					sys.exit(3)
				else:
					self.combinedColID2ColIndex[colID] = len(self.combinedColID2ColIndex)
	
	def getColIndexGivenColHeader(self, colHeader=None):
		"""
		2012.11.15
			this is from the combined column header list.
		"""
		return self.combinedColID2ColIndex.get(colHeader)
	
	def getGroupObject(self, groupIndex=None, groupName=None):
		"""
		"""
		groupObject = None
		if groupIndex is not None:
			groupObject = self.h5GroupList[groupIndex]
		elif groupName:
			if groupName[0]!='/':	#bugfix, add root in front if not there.
				groupName = '/'+groupName
			groupIndex = self.groupName2Index.get(groupName)
			if groupIndex is not None:
				groupObject = self.h5GroupList[groupIndex]
		else:	#return first group
			groupObject = self.h5GroupList[0]
		return groupObject
	
	def getColIndex(self, colID=None, groupIndex=None, groupName=None):
		"""
		
		"""
		groupObject = self.getGroupObject(groupIndex=groupIndex, groupName=groupName)
		return groupObject.getColIndex(colID=colID)
	
	def getRowIndex(self, rowID=None, groupIndex=None, groupName=None):
		"""
		"""
		groupObject = self.getGroupObject(groupIndex=groupIndex, groupName=groupName)
		return groupObject.getRowIndex(rowID=rowID)
	
	def getCellDataGivenRowColID(self, rowID=None, colID=None, groupIndex=None, groupName=None):
		"""
		"""
		groupObject = self.getGroupObject(groupIndex=groupIndex, groupName=groupName)
		return groupObject.getCellDataGivenRowColID(rowID=rowID, colID=colID)
	
	def __iter__(self):
		return self
	
	def next(self):
		"""
		2012.11.16
			combine the numeric and str row (if both exist). and return the combined row
		"""
		
		row = None
		for groupObject in self.h5GroupList:
			if self.rowIndexCursor<groupObject.dataMatrix.shape[0]:
				if row is None:
					row = list(groupObject.dataMatrix[self.rowIndexCursor])
				else:
					rowAppend = groupObject.dataMatrix[self.rowIndexCursor]
					row += list(rowAppend)
			self.rowIndexCursor += 1
		else:
			raise StopIteration
		return row
	
	def close(self):
		self.hdf5File.close()
		del self.hdf5File
	
	def writeHeader(self, headerList=None, groupIndex=None, groupName=None):
		"""
		2012.11.16
			only the first group
		"""
		self.setColIDList(colIDList=headerList, groupIndex=groupIndex, groupName=groupName)
	
	def writeOneCell(self, oneCell=None, groupIndex=None, groupName=None):
		"""
		2012.11.19
			mimic csv's writerow()
			each cell must be a tuple (readable only). list is not acceptable. 
			it's not very efficient as it's resizing the dataMatrix all the time.
		"""
		groupObject = self.getGroupObject(groupIndex=groupIndex, groupName=groupName)
		groupObject.extendDataMatrix(dataMatrix=[oneCell])

	def writeCellList(self, cellList=None, groupIndex=None, groupName=None):
		"""
		2012.11.18
			for bulk writing. more efficient.
			each cell must be a tuple (readable only). list is not acceptable. 
		"""
		groupObject = self.getGroupObject(groupIndex=groupIndex, groupName=groupName)
		groupObject.extendDataMatrix(dataMatrix=cellList)
	
	
	def setColIDList(self, colIDList=None, groupIndex=None, groupName=None):
		"""
		"""
		groupObject = self.getGroupObject(groupIndex=groupIndex, groupName=groupName)
		groupObject.setColIDList(colIDList=colIDList)
		
	def setRowIDList(self, rowIDList=None, groupIndex=None, groupName=None):
		"""
		"""
		groupObject = self.getGroupObject(groupIndex=groupIndex, groupName=groupName)
		groupObject.setRowIDList(rowIDList=rowIDList)
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		

if __name__ == '__main__':
	main_class = HDF5MatrixFile
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()