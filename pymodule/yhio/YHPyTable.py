#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s

Description:
	2012.12.15 table-data stored in pytables.
	i.e.
		reader = PyTablesMatrixFile(inputFname=filename, openMode='r')
		reader = PyTablesMatrixFile(filename, openMode='r')
		for row in reader:
			...
		tableObject = reader.getTableObject(tableName=tableName)
		for row in tableObject:
			...
		
		dtypeList = [('locus_id','i8'),('chr', HDF5MatrixFile.varLenStrType), ('start','i8'), ('stop', 'i8'), \
					('score', 'f8'), ('MAC', 'i8'), ('MAF', 'f8')]
		headerList = [row[0] for row in dtypeList]
		dtype = numpy.dtype(dtypeList)
		
		writer = PyTablesMatrixFile(inputFname=filename, openMode='w', dtype=dtype)
		writer = PyTablesMatrixFile(filename, openMode='w', dtype=dtype)
		
		if writer:
			tableObject = writer.createNewTable(tableName=tableName, dtype=dtype)
			tableObject.setColIDList(headerList)
		elif outputFname:
			writer = PyTablesMatrixFile(outputFname, openMode='w', dtype=dtype, tableName=tableName)
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
		dtypeList = [('locus_id','i8'),('chromosome', HDF5MatrixFile.varLenStrType), ('start','i8'), ('stop', 'i8'), \
					('score', 'f8'), ('MAC', 'i8'), ('MAF', 'f8')]
		if writer is None and filename:
			writer = PyTablesMatrixFile(filename, openMode='w', dtypeList=dtypeList, tableName=tableName)
			tableObject = writer.getTableObject(tableName=tableName)
		elif writer:
			tableObject = writer.createNewTable(tableName=tableName, dtypeList=dtypeList)

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import csv
import tables
import numpy
from pymodule.utils import PassingData, PassingDataList
from pymodule.ProcessOptions import  ProcessOptions
from pymodule.yhio.MatrixFile import MatrixFile
from pymodule.yhio.HDF5MatrixFile import HDF5MatrixFile, YHTableInHDF5Group, addAttributeDictToYHTableInHDF5Group 



class YHPyTable(tables.Table):
	__doc__ = __doc__
	"""
	2012.12.16 API is very similar to HDF5MatrixFile
	2012.12.16 adapted from http://pytables.github.com/cookbook/simple_table.html
	"""
	#mimics the sqlalchemy
	query = tables.Table.readWhere
	
	def __init__(self, inputFname=None, openMode='r', \
				groupName=None, tableName=None,\
				description=None,
				title='', filters=None, rowDefinition=None,\
				expectedrows=512000, **keywords):
		"""
		argument rowDefinition is skipped, just to make it compatible with HDF5MatrixFile
		"""
		self.inputFname = inputFname
		self.openMode = openMode
		self.groupName = groupName
		self.tableName = tableName
		
 		self.hdf5File = tables.openFile(inputFname, openMode)
 		self.uservars = None
 		
		if groupName is None:
			groupName = 'default'
		
		self.parentNode = self.hdf5File._getOrCreatePath('/' + groupName, True)
        
		if tableName in self.parentNode: # existing table
			description = None
		elif description is None: # pull the description from the attrs
			description = dict(self._get_description())
		
		if filters is None:
			filters = tables.Filters(complib="blosc", complevel=5, shuffle=True)

		#if self.openMode=='w':
		tables.Table.__init__(self, self.parentNode, tableName,
						description=description, title=title,
						filters=filters,
						expectedrows=expectedrows,
						_log=False, **keywords)
		
		self._c_classId = self.__class__.__name__
		
		self._processRowIDColID()
		self.no_of_rows = self.nrows	#a counter used in self.writeOneCell, self.nrows does not get updated until flush.

	def _processRowIDColID(self):
		"""
		2012.12.16 similar to SNPData.processRowIDColID()
		"""
		self.rowIDList = []
		self.rowID2rowIndex = {}
		self.colIDList = []
		self.colID2colIndex = {}
		for i in xrange(len(self.colnames)):
			colID = self.colnames[i]
			self.colIDList.append(colID)
			self.colID2colIndex[colID] = i
		
		
	
	def _get_description(self):
		# pull the description from the attrs
		for attr_name in dir(self):
			if attr_name[0] == '_': continue
			try:
				attr = getattr(self, attr_name)
			except:
				continue
			if isinstance(attr, tables.Atom):
				yield attr_name, attr

	def insert_many(self, data_generator, attr=False):
		row = self.row
		cols = self.colnames
		if not attr:
			for d in data_generator:
				for c in cols:
					row[c] = d[c]
				row.append()
		else:
			for d in data_generator:
				for c in cols:
					row[c] = getattr(d, c)
				row.append()
		self.flush()
	
	def getTableObject(self, tableIndex=None, tableName=None, groupName=None):
		"""
		2012.12.18 compatible with HDF5MatrixFile
		"""
		return self
		
	def setColIDList(self, colIDList=None, **keywords):
		"""
		"""
		pass
	
	def setRowIDList(self, rowIDList=None, **keywords):
		pass
	
	def writeHeader(self, headerList=None, tableIndex=None, tableName=None):
		"""
		2012.11.16
			only the first group
		"""
		pass
	
	def writeOneCell(self, oneCell=None, cellType=1):
		"""
		2012.12.16
			each cell is either a tuple/list or object with column-name attributes
			cellType
				1: list or tuple, in the order of table columns, not including first "id" column
				2: object with attributes whose names are same as those of the table columns
				3: a dictionary with key=column-name, value=column-value
		"""
		row = self.row
		for i in xrange(len(self.colnames)):	#assuming data in oneCell is in the same order as tableObject.colnames
			colname = self.colnames[i]
			if colname=='id':	#auto-increment the ID column if it exists
				row[colname] = self.no_of_rows + 1
				self.no_of_rows += 1
			else:
				if cellType==1:	#list/tuple
					if "id" in self.colinstances:	#colinstances maps the name of a column to its Column (see The Column class) 
							#or Cols (see The Cols class) instance.
						cellColIndex = i-1					#oneCell does not include id
					else:	#no "id" column, so same index
						cellColIndex = i
					row[colname] = oneCell[cellColIndex]
				elif cellType==2:	#object with attributes
					row[colname] = getattr(oneCell, colname, None)
				elif cellType==3:	#dictionary
					row[colname] = oneCell.get(colname, None)
		row.append()

	def writeCellList(self, cellList=None, cellType=1):
		"""
		"""
		for oneCell in cellList:
			self.writeOneCell(oneCell, cellType=cellType)
		self.flush()
	
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
		pass
		"""
		rowIndex = self.getRowIndex(rowID)
		colIndex = self.getColIndex(colID)
		
		cellData = None
		if rowIndex is not None and colIndex is not None:
			cellData = self.[rowIndex][colIndex]
		return cellData
		"""
	
	def constructColName2IndexFromHeader(self):
		"""
		2012.11.22 overwrite parent function
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
		if hasattr(self._v_attrs, name):
			sys.stderr.write("Warning: h5Group %s already has attribute %s=%s.\n"%(self.name, name, value))
			if overwrite:
				setattr(self._v_attrs, name, value)
			else:
				return False
		else:
			setattr(self._v_attrs, name, value)
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
		return getattr(self._v_attrs, name, defaultValue)
	
	def getAttributes(self):
		dc = {}
		for attributeName in self._v_attrs._f_list():
			dc[attributeName] = getattr(self._v_attrs, attributeName)
		return dc
	
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



class YHMultiPyTable(tables.File):
	"""
	2012.12.18 not ready. a version of YHPyTable which could contain multiple tables.
	2012.12.16 API is very similar to HDF5MatrixFile
	"""
	def __init__(self, inputFname=None, openMode='r', \
				firstGroupName=None, tableName=None, groupNamePrefix='group', tableNamePrefix='table',\
				rowDefinition=None, debug=0, report=0, **keywords):
		self.inputFname = inputFname
		self.header = None
		self.openMode = openMode
		self.firstGroupName = firstGroupName
		self.tableName = tableName
		self.groupNamePrefix = groupNamePrefix
		self.tableNamePrefix = tableNamePrefix
		self.rowDefinition = rowDefinition
		self.debug = debug
		self.report = report
		
		self.combinedColIDList = None	#same as header
		self.combinedColID2ColIndex = None
		
		#self.hdf5File = tables.openFile(self.inputFname, self.openMode)
		#self.root = self.hdf5File.root
		
		tables.File.__init__(self, self.inputFname, mode=self.openMode, title='', rootUEP='/', filters=_filter,\
							**keywords)
		
		if self.openMode=='r':
			self._readInData()
		elif self.openMode=='w':
			self.createNewTable(tableName=self.tableName, rowDefinition=self.rowDefinition)
		
		self._c_classId = self.__class__.__name__
	
	def _createNewGroup(self, groupName=None, rowDefinition=None, tableName=None):
		"""
		2012.12.16 this would also create a new table object within the newly-created group
			return tableObject
			** deprecated
		"""
		if groupName is None or self.root.__contains__(groupName):
			groupName = self._getNewGroupName()
		groupObject = self.createGroup(self.root, groupName)
		if tableName is None:
			tableName = self._getNewTableName()
		tableObject = self.createTable(groupObject, tableName, rowDefinition, title=None)
		return groupObject
	
	def createNewTable(self, tableName=None, rowDefinition=None):
		"""
		2012.12.16 this would also create a new group object if no group exists.
		"""
		groupObject = self._getGroupObject(groupIndex=None, groupName=None)
		if not groupObject:	#create a new one if it's not available.
			groupName = self._getNewGroupName()
			groupObject = self.createGroup(self.root, groupName)
		
		if tableName is None:
			tableName = self._getNewTableName(groupObject=groupObject)
		
		tableObject = self.createTable(groupObject, tableName, rowDefinition, title=None)
		return tableObject
	
	def _getNewGroupName(self, groupNamePrefix=None):
		"""
		2012.11.19
		"""
		if not groupNamePrefix:
			groupNamePrefix = self.groupNamePrefix
		i = len(self.root._v_groups)
		groupName = "%s%s"%(groupNamePrefix, i)
		while self.root.__contains__(groupName):	#stop until a unique name shows up
			i += 1
			groupName = "%s%s"%(groupNamePrefix, i)
		return groupName

	def _getNewTableName(self, tableNamePrefix=None, groupObject=None):
		"""
		2012.12.15
		"""
		if groupObject is None:
			groupObject = self._getGroupObject()
		if not tableNamePrefix:
			tableNamePrefix = self.tableNamePrefix
		i = len(groupObject._v_leaves)
		tableName = "%s%s"%(tableNamePrefix, i)
		while tableName in groupObject.__contains__(tableName):	#stop until a unique name shows up
			i += 1
			tableName = "%s%s"%(tableNamePrefix, i)
		return tableName
	
	def _readInData(self):
		"""
		2012.12.16
		"""
		self._setupCombinedColumnIDMapper()
	
	def _setupCombinedColumnIDMapper(self,):
		"""
		2012.12.16
			use the first table only
		"""
		self.combinedColID2ColIndex = {}
		tableObject = self.getTableObject()
		self.header = []
		self.combinedColIDList = self.header
		for colID in tableObject.colnames:
			if colID in self.combinedColID2ColIndex:
				sys.stderr.write("Error: column ID %s already used in column %s.\n"%(colID, self.combinedColID2ColIndex.get(colID)))
				sys.exit(3)
			else:
				self.combinedColID2ColIndex[colID] = len(self.combinedColID2ColIndex)
				self.header.append(colID)
		return self.combinedColID2ColIndex
	
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
	
	def _getGroupObject(self, groupIndex=None, groupName=None):
		"""
		"""
		groupObject = None
		if groupName:
			groupObject = self.getNode("/", name=groupName, classname='Group')
		else:	#return first group
			if len(self.root._v_groups.values())>0:
				groupObject = self.root._v_groups.values()[0]	#first group
		return groupObject
	
	def getTableObject(self, tableIndex=None, tableName=None, groupName=None):
		"""
		2012.12.16
		"""
		groupObject = self._getGroupObject(groupIndex=None, groupName=groupName)
		tableObject = None
		if tableName:
			tableObject = groupObject._f_getChild(tableName)
		else:
			tableObject = groupObject._v_leaves.values()[0]	#first leave
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
		2012.12.16
			go through each leaf (table & array), iteratively.
			usually the file has only one leaf (=table)
		"""
		
		row = None
		pdata = PassingDataList()
		for groupObject in self.walkGroups("/"):
			for leafNode in  groupObject._f_walkNodes('Leaf'):
				for row in leafNode:
					yield row
		"""
			for colID  in self.header:	#iteration over header is in the same order as the ascending order of colIndex
				colIndex = self.combinedColID2ColIndex.get(colID)
				if colIndex < len(row):
					setattr(pdata, colID, row[colIndex])
		self.rowIndexCursor += 1
		return pdata
		"""
	
	def writeHeader(self, headerList=None, tableIndex=None, tableName=None):
		"""
		2012.11.16
			only the first group
		"""
		pass
	
	def writeOneCell(self, oneCell=None, tableIndex=None, tableName=None, cellType=1):
		"""
		2012.12.16
			mimic csv's writerow()
			cellType
				1: list or tuple, in the order of table columns, not including first "id" column
				2: object with attributes whose names are same as those of the table columns
				3: a dictionary with key=column-name, value=column-value
		"""
		tableObject = self.getTableObject(tableIndex=tableIndex, tableName=tableName)
		tableObject.writeOneCell(oneCell, cellType=cellType)
		"""
		row = tableObject.row
		for i in xrange(len(tableObject.colnames)):	#assuming data in oneCell is in the same order as tableObject.colnames
			colname = tableObject.colnames[i]
			row[colname] = oneCell[i]
		row.append()
		"""
	def writeCellList(self, cellList=None, tableIndex=None, tableName=None, cellType=1):
		"""
		"""
		tableObject = self.getTableObject(tableIndex=tableIndex, tableName=tableName)
		tableObject.writeCellList(cellList, cellType=cellType)
		"""
		for oneCell in cellList:
			self.writeOneCell(oneCell, tableIndex=tableIndex, tableName=tableName)
		tableObject.flush()
		"""
	
	def setColIDList(self, colIDList=None, tableIndex=None, tableName=None):
		"""
		"""
		pass
		
	def setRowIDList(self, rowIDList=None, tableIndex=None, tableName=None):
		"""
		"""
		pass
	
	def addAttributeDict(self, attributeDict=None, tableObject=None):
		if tableObject is None:
			tableObject = self.getTableObject()
		
		addAttributeDictToYHTableInHDF5Group(tableObject=tableObject, attributeDict=attributeDict)
	
	def addAttribute(self, name=None, value=None, overwrite=True, tableIndex=None, tableName=None):
		"""
		2012.11.28 find the tableObject and let it do the job
		"""
		tableObject = self.getTableObject(tableIndex=tableIndex, tableName=tableName)
		return tableObject.addAttribute(name=name, value=value, overwrite=overwrite)
	
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

