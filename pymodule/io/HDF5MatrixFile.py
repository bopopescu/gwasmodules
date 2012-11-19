#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s

Description:
	2012.11.15 a csv file format, stored in HDF5.
		This HDF5 file has two default datasets. One is a numerical 2D numericGroup.
			The other is 2D "str_matrix" of type string.
		The two datasets match the rows to each other. The number of columns may differ.
			This setup is due to the fact that within one dataset, only one data type could exist.
	i.e.
		reader = HDF5MatrixFile(inputFname=filename, openMode='r')
		reader = HDF5MatrixFile(filename, openMode='r')
		
		writer = HDF5MatrixFile(inputFname=filename, openMode='w')
		writer = HDF5MatrixFile(filename, openMode='w')
		writer.writeHeader(numericGroupHeader=['locus_id', 'start', 'stop', 'score', 'mac', 'maf'], strGroupHeader=['chr'])
		writer.writerow(numericGroupRow=[1,2,3, 0.5, 15, 0.15], strGroupRow=['X'])
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule.ProcessOptions import  ProcessOptions
from pymodule import utils, figureOutDelimiter
from pymodule import MatrixFile
import csv
import h5py
import numpy
str_type = h5py.new_vlen(str)

class HDF5GroupWrapper(object):
	def __init__(self, h5Group = None, openMode='r', dataMatrixDtype='f'):
		self.h5Group = h5Group
		self.openMode = openMode
		self.dataMatrixDtype = dataMatrixDtype
		
		self.dataMatrixDSName = "dataMatrix"
		self.rowIDListDSName = "rowIDList"
		self.colIDListDSName = "colIDList"
		if self.openMode=='r':
			self._readInData()
		elif self.openMode =='w':
			self._createDatasetSkeletonForOneGroup(h5Group=self.h5Group, dtype=self.dataMatrixDtype)

	def _createDatasetSkeletonForOneGroup(self, h5Group=None, dtype=None):
		"""
		2012.11.18
		"""
		if h5Group is None:
			h5Group = self.h5Group
		if dtype is None:
			dtype = self.dataMatrixDtype
		if dtype==str_type:
			defaultData = [['0']]
		else:
			defaultData = [[]]
		self.dataMatrix = h5Group.create_dataset(self.dataMatrixDSName, data=defaultData, dtype=dtype, \
										maxshape=(None, None), compression='gzip', compression_opts=4)
		#for str_type, there must be an element inside
		self.rowIDList = h5Group.create_dataset(self.rowIDListDSName, data=['0'], dtype=str_type, maxshape=(None,),\
											compression='gzip', compression_opts=4)
		self.colIDList = h5Group.create_dataset(self.colIDListDSName, data=['0'], dtype=str_type, maxshape=(None,),\
											compression='gzip', compression_opts=4)
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
	
	def _extendDataMatrix(self, dataMatrix=None):
		"""
		2012.11.16
			dataMatrix has to be 2D
		"""
		if dataMatrix:
			m = len(dataMatrix)
			n = len(dataMatrix[0])
			s, t = self.dataMatrix.shape
			if t!=0 and n!=t:
				sys.stderr.write("Error: dataset is of shape (%s,%s). The extension is of shape (%s,%s). Different number of columns.\n"%\
								(s, t, m, n))
				sys.exit(3)
			else:
				self.dataMatrix.resize((s+m, n))
				self.dataMatrix[s:s+m] = dataMatrix

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
	
class HDF5MatrixFile(MatrixFile):
	__doc__ = __doc__
	option_default_dict = MatrixFile.option_default_dict.copy()
	option_default_dict.update({
							('dtype', 0, ): ['f', '', 1, 'data type in the main numeric numericGroup. candidates are i, f8, etc.'],\
							})
	def __init__(self, inputFname=None, **keywords):
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		if not self.inputFname:
			self.inputFname = inputFname
		
		self.combinedColID2ColIndex = None
		
		self.hdf5File = h5py.File(self.inputFname, self.openMode)
		self.numericGroupName = 'numeric'
		self.strGroupName = 'str'
		
		if self.openMode=='r':
			self._readInData()
		elif self.openMode=='w':
			self.numericGroup = HDF5GroupWrapper(self.hdf5File.create_group(self.numericGroupName), \
										openMode=self.openMode, dataMatrixDtype=self.dtype)
			self.strGroup = HDF5GroupWrapper(self.hdf5File.create_group(self.strGroupName), openMode=self.openMode,\
											dataMatrixDtype=str_type)
			
		self.rowIndexCursor = 0	#2012.11.16 for iteration
	
	def _readInData(self):
		"""
		2012.11.16
		"""
		self.numericGroup = HDF5GroupWrapper(self.hdf5File[self.numericGroupName], openMode=self.openMode,\
											dataMatrixDtype=self.dtype)
		self.strGroup = HDF5GroupWrapper(self.hdf5File[self.strGroupName], openMode=self.openMode, \
										dataMatrixDtype=self.dtype)
		self._setupCombinedColumnIDMapper()
		
	def _setupCombinedColumnIDMapper(self,):
		"""
		2012.11.16
			the strGroup is arranged after numericGroup.
		"""
		self.combinedColID2ColIndex = {}
		colIDList = self.numericGroup.colIDList + self.strGroup.colIDList
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
	
	def getColIndex(self, colID=None, groupObject=None):
		"""
		
		"""
		if groupObject is None:
			groupObject = self.numericGroup
		
		return groupObject.getColIndex(colID=colID)
	
	def getRowIndex(self, rowID=None, groupObject=None):
		"""
		"""
		if groupObject is None:
			groupObject = self.numericGroup
		return groupObject.getRowIndex(rowID=rowID)
	
	def getNumericColIndex(self, colID=None):
		"""
		"""
		return self.getColIndex(colID=colID, groupObject=self.numericGroup)
	
	def getNumericRowIndex(self, rowID=None):
		"""
		"""
		return self.getRowIndex(rowID=rowID, groupObject=self.numericGroup)
	
	def getStrColIndex(self, colID=None):
		"""
		"""
		return self.getColIndex(colID=colID, groupObject=self.strGroup)
	
	def getStrRowIndex(self, rowID=None):
		"""
		"""
		return self.getRowIndex(rowID=rowID, groupObject=self.strGroup)
	
	def getCellDataGivenRowColID(self, rowID=None, colID=None, groupObject=None):
		"""
		"""
		if groupObject is None:
			groupObject = self.numericGroup
		return groupObject.getCellDataGivenRowColID(rowID=rowID, colID=colID)
	
	def getNumericCellDataGivenRowColID(self, rowID=None, colID=None):
		"""
		"""
		return self.getCellDataGivenRowColID(rowID, colID, groupObject=self.numericGroup)
	
	def getStrCellDataGivenRowColID(self, rowID=None, colID=None):
		"""
		"""
		return self.getCellDataGivenRowColID(rowID, colID, groupObject=self.strGroup)
	
	def __iter__(self):
		return self
	
	def next(self):
		"""
		2012.11.16
			combine the numeric and str row (if both exist). and return the combined row
		"""
		if self.rowIndexCursor<self.numericGroup.dataMatrix.shape[0]:
			row = self.numericGroup.dataMatrix[self.rowIndexCursor]
			if self.rowIndexCursor<self.strGroup.dataMatrix.shape[0]:
				strRow = self.strGroup.dataMatrix[self.rowIndexCursor]
				row = row + strRow
			self.rowIndexCursor += 1
		else:
			raise StopIteration
		return row
	
	def close(self):
		self.hdf5File.close()
		del self.hdf5File
	
	def writeHeader(self, numericGroupHeader=None, strGroupHeader=None):
		"""
		2012.11.16
		"""
		if numericGroupHeader:
			self.setNumericGroupColIDList(colIDList=numericGroupHeader)
		if strGroupHeader:
			self.setStrGroupColIDList(colIDList=strGroupHeader)
	
	def writerow(self, numericGroupRow=None, strGroupRow=None):
		"""
		mimic csv's writerow()
			it's not very efficient as it's resizing the dataMatrix all the time.
		"""
		if numericGroupRow:
			self.numericGroup._extendDataMatrix(dataMatrix=[numericGroupRow])
		if strGroupRow:
			self.strGroup._extendDataMatrix(dataMatrix=[strGroupRow])

	def writeMatrix(self, numericGroupMatrix=None, strGroupMatrix=None):
		"""
		2012.11.18
			for bulk writing, both arguments are 2D list or array
		"""
		if numericGroupMatrix:
			self.numericGroup._extendDataMatrix(dataMatrix=numericGroupMatrix)
		if strGroupMatrix:
			self.strGroup._extendDataMatrix(dataMatrix=strGroupMatrix)
				
	
	def _extendDataMatrix(self, dataMatrix=None, groupObject=None):
		"""
		2012.11.16
			dataMatrix has to be 2D
		"""
		if groupObject is None:
			groupObject = self.numericGroup
		if dataMatrix:
			groupObject._extendDataMatrix(dataMatrix=dataMatrix)
	
	def _setColIDList(self, colIDList=None, groupObject=None):
		"""
		"""
		if groupObject is None:
			groupObject = self.numericGroup
		groupObject.setColIDList(colIDList=colIDList)
		
	def _setRowIDList(self, rowIDList=None, groupObject=None):
		"""
		"""
		if groupObject is None:
			groupObject = self.numericGroup
		groupObject.setRowIDList(rowIDList=rowIDList)
	
	def setNumericGroupColIDList(self, colIDList=None):
		"""
		"""
		self._setColIDList(colIDList=colIDList, groupObject=self.numericGroup)
	
	def setNumericGroupRowIDList(self, rowIDList=None):
		"""
		"""
		self._setRowIDList(rowIDList=rowIDList, groupObject=self.numericGroup)
		
	def setStrGroupColIDList(self, colIDList=None):
		"""
		"""
		self._setColIDList(colIDList=colIDList, groupObject=self.strGroup)
	
	def setStrGroupRowIDList(self, rowIDList=None):
		"""
		"""
		self._setRowIDList(rowIDList=rowIDList, groupObject=self.strGroup)
	
	
	
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