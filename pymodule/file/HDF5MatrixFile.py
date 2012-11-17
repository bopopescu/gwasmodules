#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s

Description:
	2012.11.15 a csv file format, stored in HDF5.
		This HDF5 file has two default datasets. One is a numerical 2D numericDataset.
			The other is 2D "str_matrix" of type string.
		The two datasets match the rows to each other. The number of columns may differ.
			This setup is due to the fact that within one dataset, only one data type could exist.

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

class HDF5MatrixFile(MatrixFile):
	__doc__ = __doc__
	option_default_dict = MatrixFile.option_default_dict.copy()
	option_default_dict.update({
							('dtype', 0, ): ['f', '', 1, 'data type in the main numeric numericDataset. candidates are i, f8, etc.'],\
							})
	def __init__(self, inputFname=None, **keywords):
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		if not self.inputFname:
			self.inputFname = inputFname
		
		self.combinedColID2ColIndex = None
		
		self.hdf5File = h5py.File(self.inputFname, self.openMode)
		if self.openMode=='r':
			self.readInData()
		elif self.openMode=='w':
			self.numericDataset = self.hdf5File.create_dataset("numericDataset", data=[[]], dtype=self.dtype, maxshape=(None, None))
			self.numericDataset.attrs['rowIDList'] = []
			self.numericDataset.attrs['colIDList'] = []
			self.strDataset = self.hdf5File.create_dataset("strDataset", data=[[]], dtype=str_type, maxshape=(None, None))
			self.strDataset.attrs['rowIDList'] = []
			self.strDataset.attrs['colIDList'] = []
		self.rowIndexCursor = 0	#2012.11.16 for iteration
	
	def readInData(self):
		"""
		2012.11.16
		"""
		self.numericDataset = self.hdf5File['numericDataset']
		self.strDataset = self.hdf5File['strDataset']
		self.processRowIDColID(self.numericDataset)
		self.processRowIDColID(self.strDataset)
		self._setupCombinedColumnIDMapper()
	
	def processRowIDColID(self, h5Dataset=None):
		"""
		2012.11.16 similar to SNPData.processRowIDColID()
		"""
		rowIDList = h5Dataset.attrs['rowIDList']
		colIDList = h5Dataset.attrs['colIDList']
		
		
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
		
		h5Dataset.attrs['rowID2rowIndex'] = rowID2rowIndex
		h5Dataset.attrs['colID2colIndex'] = colID2colIndex
	
	def _setupCombinedColumnIDMapper(self,):
		"""
		2012.11.16
			the strDataset is arranged after numericDataset.
		"""
		self.combinedColID2ColIndex = {}
		colIDList = self.numericDataset.attrs['colIDList'] + self.strDataset.attrs['colIDList']
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
	
	def getColIndex(self, colID=None, h5Dataset=None):
		"""
		
		"""
		if h5Dataset is None:
			h5Dataset = self.numericDataset
		colID2colIndex = h5Dataset.attrs.get('colID2colIndex', None)
		colIndex = None
		if colID2colIndex:
			colIndex = colID2colIndex.get(colID, None)
		return colIndex
	
	def getRowIndex(self, rowID=None, h5Dataset=None):
		"""
		"""
		if h5Dataset is None:
			h5Dataset = self.numericDataset
		rowID2rowIndex = h5Dataset.attrs.get('rowID2rowIndex', None)
		rowIndex = None
		if rowID2rowIndex:
			rowIndex = rowID2rowIndex.get(rowID, None)
		return rowIndex
	
	def getNumericColIndex(self, colID=None):
		"""
		"""
		return self.getColIndex(colID=colID, h5Dataset=self.numericDataset)
	
	def getNumericRowIndex(self, rowID=None):
		"""
		"""
		return self.getRowIndex(rowID=rowID, h5Dataset=self.numericDataset)
	
	def getStrColIndex(self, colID=None):
		"""
		"""
		return self.getColIndex(colID=colID, h5Dataset=self.strDataset)
	
	def getStrRowIndex(self, rowID=None):
		"""
		"""
		return self.getRowIndex(rowID=rowID, h5Dataset=self.strDataset)
	
	def getCellDataGivenRowColID(self, rowID=None, colID=None, h5Dataset=None):
		"""
		"""
		if h5Dataset is None:
			h5Dataset = self.numericDataset
		
		rowIndex = self.getRowIndex(rowID, h5Dataset=h5Dataset)
		colIndex = self.getColIndex(colID, h5Dataset=h5Dataset)
		
		cellData = None
		if rowIndex is not None and colIndex is not None:
			cellData = h5Dataset[rowIndex][colIndex]
		return cellData
	
	def getNumericCellDataGivenRowColID(self, rowID=None, colID=None):
		"""
		"""
		return self.getCellDataGivenRowColID(rowID, colID, h5Dataset=self.numericDataset)
	
	def getStrCellDataGivenRowColID(self, rowID=None, colID=None):
		"""
		"""
		return self.getCellDataGivenRowColID(rowID, colID, h5Dataset=self.strDataset)
	
	def __iter__(self):
		return self
	
	def next(self):
		"""
		2012.11.16
			combine the numeric and str row (if both exist). and return the combined row
		"""
		if self.rowIndexCursor<self.numericDataset.shape[0]:
			row = self.numericDataset[self.rowIndexCursor]
			if self.rowIndexCursor<self.strDataset.shape[0]:
				strRow = self.strDataset[self.rowIndexCursor]
				row = row + strRow
			self.rowIndexCursor += 1
		else:
			raise StopIteration
		return row
	
	def close(self):
		self.hdf5File.close()
		del self.hdf5File
	
	def writeHeader(self, numericDatasetHeader=None, strDatasetHeader=None):
		"""
		2012.11.16
		"""
		if numericDatasetHeader:
			self.setNumericDatasetColIDList(colIDList=numericDatasetHeader)
		if strDatasetHeader:
			self.setStrDatasetColIDList(colIDList=strDatasetHeader)
	
	def writerow(self, numericDatasetRow=None, strDatasetRow=None):
		"""
		mimic csv's writerow()
		"""
		if numericDatasetRow:
			self.extendDataset(data_matrix=[numericDatasetRow], h5Dataset=self.numericDataset)
		if strDatasetRow:
			self.extendDataset(data_matrix=[strDatasetRow], h5Dataset=self.strDataset)
		
	
	def extendDataset(self, data_matrix=None, h5Dataset=None):
		"""
		2012.11.16
			data_matrix has to be 2D
		"""
		if h5Dataset is None:
			h5Dataset = self.numericDataset
		if data_matrix:
			m = len(data_matrix)
			n = len(data_matrix[0])
			s, t = h5Dataset.shape
			if n!=t:
				sys.stderr.write("Error: dataset is of shape (%s,%s). The extension is of shape (%s,%s). Different number of columns.\n"%\
								(s, t, m, n))
				sys.exit(3)
			else:
				h5Dataset.resize((s+m, t))
				h5Dataset[s:s+m] = data_matrix
		
	
	def setDatasetColIDList(self, colIDList=None, h5Dataset=None):
		"""
		"""
		if h5Dataset is None:
			h5Dataset = self.numericDataset
		h5Dataset.attrs['colIDList'] = colIDList
		
	def setDatasetRowIDList(self, rowIDList=None, h5Dataset=None):
		"""
		"""
		if h5Dataset is None:
			h5Dataset = self.numericDataset
		h5Dataset.attrs['rowIDList'] = rowIDList
	
	def setNumericDatasetColIDList(self, colIDList=None):
		"""
		"""
		self.setDatasetColIDList(colIDList=colIDList, h5Dataset=self.numericDataset)
	
	def setNumericDatasetRowIDList(self, rowIDList=None):
		"""
		"""
		self.setDatasetRowIDList(rowIDList=rowIDList, h5Dataset=self.numericDataset)
		
	def setStrDatasetColIDList(self, colIDList=None):
		"""
		"""
		self.setDatasetColIDList(colIDList=colIDList, h5Dataset=self.strDataset)
	
	def setStrDatasetRowIDList(self, rowIDList=None):
		"""
		"""
		self.setDatasetRowIDList(rowIDList=rowIDList, h5Dataset=self.strDataset)
	
	
	
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