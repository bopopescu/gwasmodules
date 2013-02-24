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
from tables.file import _checkfilters

class YHTable(tables.Table, YHTableInHDF5Group):
	"""
	2012.12.16 API is very similar to HDF5MatrixFile
	2012.12.16 adapted from http://pytables.github.com/cookbook/simple_table.html
	"""
	#mimics the sqlalchemy
	query = tables.Table.readWhere
	
	def __init__(self, parentNode=None, tableName=None,\
				description=None, rowDefinition=None,\
				title='', filters=None, \
				expectedrows=52000, chunkshape=None, byteorder=None, **keywords):
		"""
	description
		  An IsDescription subclass or a dictionary where the keys are the field
		  names, and the values the type definitions. In addition, a pure NumPy
		  dtype is accepted.  And it can be also a recarray NumPy object,
		  RecArray numarray object or NestedRecArray. If None, the table metadata
		  is read from disk, else, it's taken from previous parameters.
	rowDefinition is a backup of description, to make the class compatible with HDF5MatrixFile
		"""
		self.tableName = tableName
		self.parentNode = parentNode
		
		if tableName in self.parentNode: # existing table
			description = None
		elif description is None: # pull the description from the attrs
			if rowDefinition is not None:	#a backup of description
				description = rowDefinition
			else:
				description = dict(self._get_description())
		
		if filters is None:
			filters = tables.Filters(complib="blosc", complevel=5, shuffle=True)
		
		#YHTableInHDF5Group.__init__(self, newGroup=0)
		#if self.openMode=='w':
		tables.Table.__init__(self, self.parentNode, tableName,
						description=description, title=title,
						filters=filters,
						expectedrows=expectedrows,
						_log=False, **keywords)
		
		self._c_classId = self.__class__.__name__
		
		self._processRowIDColID()

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
		
		self.no_of_rows = self.nrows	#a counter used in self.writeOneCell, self.nrows does not get updated until flush.
		
		
	
	def _get_description(self):
		# pull the description from the attrs
		for attr_name in dir(self):
			if attr_name[0] == '_':
				continue
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
		#2013.1.9 make sure enough columns in oneCell
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
				self.no_of_rows += 1
				row[colname] = self.no_of_rows
			else:
				if cellType==1:	#list/tuple
					if "id" in self.colinstances:	#colinstances maps the name of a column to its Column (see The Column class) 
							#or Cols (see The Cols class) instance.
						cellColIndex = i-1					#oneCell does not include id
					else:	#no "id" column, so same index
						cellColIndex = i
					if cellColIndex<len(oneCell):	#2013.1.9 make sure enough columns in oneCell
						value = oneCell[cellColIndex]
					else:
						continue
				elif cellType==2:	#object with attributes
					value = getattr(oneCell, colname, None)
				elif cellType==3:	#dictionary
					value = oneCell.get(colname, None)
				if value is not None:
					row[colname] = value
		row.append()
	
	def appendOneRow(self, oneCell=None, cellType=1):
		"""
		2013.1.9
		"""
		self.writeOneCell(oneCell=oneCell, cellType=cellType)
	
	def writeCellList(self, cellList=None, cellType=1):
		"""
		"""
		for oneCell in cellList:
			self.writeOneCell(oneCell, cellType=cellType)
		self.flush()
	
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
	
	def getAttribute(self, name=None, defaultValue=None):
		return getattr(self._v_attrs, name, defaultValue)
	
	def getAttributes(self):
		dc = {}
		for attributeName in self._v_attrs._f_list():
			dc[attributeName] = getattr(self._v_attrs, attributeName)
		return dc
	
	#2102.12.20 overwrite YHTableInHDF5Group.next()
	#next = tables.Table.next	tables.Table doesn't have next
	
	
	def close(self):
		"""
		2012.12.20 do nothing
		"""
		pass


class YHSingleTableFile(YHTable):
	"""
	2012.12.16 API is very similar to HDF5MatrixFile
	2012.12.16 adapted from http://pytables.github.com/cookbook/simple_table.html
	"""
	#mimics the sqlalchemy	
	def __init__(self, inputFname=None, openMode='r', \
				groupName=None, tableName=None,\
				description=None,
				title='', filters=None, rowDefinition=None,\
				expectedrows=512000, **keywords):
		"""
		rowDefinition is backup of description, to make it compatible with HDF5MatrixFile
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
		
		YHTable.__init__(self, self.parentNode, tableName,
						description=description, rowDefinition=rowDefinition, \
						title=title,
						filters=filters,
						expectedrows=expectedrows,
						_log=False, **keywords)
		
		self._c_classId = self.__class__.__name__
	
	def getTableObject(self, tableIndex=None, tableName=None, groupName=None):
		"""
		2012.12.18 compatible with HDF5MatrixFile
		"""
		return self

class transplant:
	"""
	2012.12.20 from http://code.activestate.com/recipes/81732/
	"""
	def __init__(self, method, host, method_name=None):
		self.host = host
		self.method = method
		setattr(host, method_name or method.__name__, self)	#calling method_name is effectively calling this transplant class

	def __call__(self, *args, **kwargs):
		nargs = [self.host]
		nargs.extend(args)
		return apply(self.method, nargs, kwargs)

def funcToMethod(func, newClass, method_name=None, instance=None):
	"""
	2012.12.20 from http://code.activestate.com/recipes/81732/,
		to attach func to newClass with new method_name and instance
	"""
	import new
	method = new.instancemethod(func, instance, newClass)
	if not method_name:
		method_name=func.__name__
	instance.__dict__[method_name]=method

#def turnTableIntoYHTable():
#	for functionName in ['getAttributes', 'getAttribute', 'addAttribute', 'addAttributeDict']:
#		function = getattr(YHTable, functionName)
#		funcToMethod(function, tables.Table, functionName)

#turnTableIntoYHTable()


def armTableObjectWithYHTableMethods(tableObject=None):
	"""
	2012.12.20 function to arm a tableObject (read from a pytables file) with functions from YHTable and its parental class YHTableInHDF5Group
		almost like type-cast (python doesn't support type cast) 
	"""
	import types
	for functionName in ['_processRowIDColID', 'addAttribute', 'addAttributeDict', 'addObjectAttributeToSet', \
						'addObjectListAttributeToSet', 'close', 'constructColName2IndexFromHeader', \
						'getAttribute', 'getAttributes',  'getCellDataGivenRowColID', \
						'getColIndex', 'getColIndexGivenColHeader', \
						'getListAttributeInStr', \
						'getRowIndex', \
						'writeOneCell', 'writeCellList', 'setColIDList', 'setRowIDList']:
		if functionName in tableObject.__dict__:	#make sure not overwriting tables.Table's own function
			continue 
		if functionName in YHTable.__dict__:	#YHTable has precedence
			ClassToBorrow = YHTable
		elif functionName in YHTableInHDF5Group.__dict__:
			ClassToBorrow = YHTableInHDF5Group
		else:
			continue
		function = ClassToBorrow.__dict__[functionName]
		
		#function = getattr(YHTable, functionName)
		## use getattr() from above to get the function will result in this error:
		###"TypeError: unbound method getAttribute() must be called with YHTable instance as first argument (got Table instance instead)"
		
		tableObject.__dict__[functionName] = types.MethodType(function, tableObject, tableObject.__class__)
		
		#transplant(function, tableObject)	#another way of doing this, won't work because these functions are still bound to their old class
		
		#function.im_class = tables.Table	#read-only attributes for bound methods
		#func.im_func=fun
		#function.im_self = tableObject
		#funcToMethod(function, tables.Table, functionName, instance=tableObject)
	tableObject._processRowIDColID()
	return tableObject

class YHFile(tables.File, HDF5MatrixFile):
	"""
	2012.12.18 not ready. a pytables file that could contain multiple tables.
	2012.12.16 API is very similar to HDF5MatrixFile
	i.e.
	
	#open a file with a defined table structure
	writer = YHFile(self.outputFname, openMode='w', rowDefinition=CountAssociationLocusTable)
	...
	
	#read a file
	reader = YHFile(inputFname, openMode='r')
	...
	
	#open a file without passing the table structure.
	writer = YHFile(self.outputFname, openMode='w')
	# create a table with proper structure
	writer.createNewTable(rowDefinition=CountAssociationLocusTable)
	...
	"""
	def __init__(self, inputFname=None, openMode='r', \
				tableName=None, groupNamePrefix='group', tableNamePrefix='table',\
				rowDefinition=None, filters=None, expectedrows=500000, \
				autoRead=True, autoWrite=True, \
				debug=0, report=0, **keywords):
		self.inputFname = inputFname
		self.header = None
		self.openMode = openMode
		self.tableName = tableName
		self.groupNamePrefix = groupNamePrefix
		self.tableNamePrefix = tableNamePrefix
		self.rowDefinition = rowDefinition
		self.expectedrows = expectedrows
		self.autoRead = autoRead	#whether to automatically read in data if openMode is in reading mode
		self.autoWrite = autoWrite	#whether to automatically create table objects if openMode is in writing mode
		
		self.debug = debug
		self.report = report
		
		self.combinedColIDList = None	#same as header
		self.combinedColID2ColIndex = None
		
		#to get table objects fast, without retrieving them from their parental groups.
		self.tableObjectList = []
		self.tablePath2Index = {}
		
		#self.hdf5File = tables.openFile(self.inputFname, self.openMode)
		#self.root = self.hdf5File.root
		if filters is None:
			filters = tables.Filters(complib="blosc", complevel=5, shuffle=True)
		
		tables.File.__init__(self, self.inputFname, mode=self.openMode, title='', rootUEP='/', filters=filters,\
							**keywords)
		
		if self.openMode=='r' or self.openMode=='a':
			self._setupCombinedColumnIDMapper()
			if self.autoRead:
				self._readInData(tableName=self.tableName)
		
		if self.autoWrite and self.openMode=='w' and self.rowDefinition is not None:
			self.createNewTable(tableName=self.tableName, rowDefinition=self.rowDefinition, expectedrows=self.expectedrows)
		
		self._c_classId = self.__class__.__name__
	
	def _createNewGroup(self, groupName=None, **keywords):
		"""
		2012.12.16 this would also create a new table object within the newly-created group
			return tableObject
			** deprecated
		"""
		if groupName is None or self.root.__contains__(groupName):
			groupName = self._getNewGroupName()
		groupObject = self.createGroup(self.root, groupName)
		return groupObject
	
	def _appendNewTable(self, tableObject=None):
		"""
		2012.12.20
		"""
		if tableObject is not None:
			tablePath = tableObject._v_pathname	##the tableName preceds with "/"
			if tablePath in self.tablePath2Index:	 
				sys.stderr.write("ERROR, table %s already in self.tablePath2Index, index=%s.\n"%(tablePath,\
																			self.tablePath2Index.get(tablePath)))
				sys.exit(3)
			self.tablePath2Index[tablePath] = len(self.tablePath2Index)
			self.tableObjectList.append(tableObject)
	
	def _fetchTableObject(self, tablePath):
		"""
		2012.12.20
		"""
		tableObject = None
		if tablePath in self.tablePath2Index:
			tableObject = self.tableObjectList[self.tablePath2Index.get(tablePath)]
		return tableObject
	
	def createNewTable(self, tableName=None, rowDefinition=None, title="",
					filters=None, expectedrows=50000,
					chunkshape=None, byteorder=None):
		"""
		2012.12.20 this would also create a new group object if no group exists.
			symlink to createTable()
		"""
		groupObject = self._getGroupObject()
		if not groupObject:	#create a new one if it's not available.
			groupName = self._getNewGroupName()
			groupObject = self.createGroup(self.root, groupName)
		
		if tableName is None:
			tableName = self._getNewTableName(groupObject=groupObject)
		if rowDefinition is None:
			rowDefinition = self.rowDefinition
		return self.createTable(where=groupObject, tableName=tableName, description=rowDefinition, \
							rowDefinition=rowDefinition, title=title, filters=filters, expectedrows=expectedrows, \
							chunkshape=chunkshape, byteorder=byteorder, createparents=False)
	
	def createTable(self, where=None, tableName=None, description=None, rowDefinition=None, title="",
					filters=None, expectedrows=50000,
					chunkshape=None, byteorder=None,
					createparents=False):
		"""
		2012.12.20 overwrite the parent (stock) function
		
		Create a new table with the given name in where location.

		Parameters
		----------
		where : str or Group
			The parent group from which the new table will hang. It can be a
			path string (for example '/level1/leaf5'), or a Group instance
			(see :ref:`GroupClassDescr`).
		name : str
			The name of the new table.
		description : Description
			This is an object that describes the table, i.e. how
			many columns it has, their names, types, shapes, etc.  It
			can be any of the following:

				* *A user-defined class*: This should inherit from the
				  IsDescription class (see :ref:`IsDescriptionClassDescr`) where
				  table fields are specified.
				* *A dictionary*: For example, when you do not know beforehand
				  which structure your table will have).
				* *A Description instance*: You can use the description
				  attribute of another table to create a new one with the same
				  structure.
				* *A NumPy dtype*: A completely general structured NumPy dtype.
				* *A NumPy (structured) array instance*: The dtype of this
				  structured array will be used as the description.  Also, in
				  case the array has actual data, it will be injected into the
				  newly created table.
				* *A RecArray instance (deprecated)*: Object from the numarray
				  package.  This does not give you the possibility to create a
				  nested table.  Array data is injected into the new table.
				* *A NestedRecArray instance (deprecated)*: If you want to have
				  nested columns in your table and you are using numarray, you
				  can use this object. Array data is injected into the new
				  table.

		title : str
			A description for this node (it sets the TITLE HDF5 attribute
			on disk).
		filters : Filters
			An instance of the Filters class (see :ref:`FiltersClassDescr`) that
			provides information about the desired I/O filters to be applied
			during the life of this object.
		expectedrows : int
			A user estimate of the number of records that will be in the table.
			If not provided, the default value is EXPECTED_ROWS_TABLE (see
			:file:`tables/parameters.py`). If you plan to create a bigger table
			try providing a guess; this will optimize the HDF5 B-Tree creation
			and management process time and memory used.
		chunkshape
			The shape of the data chunk to be read or written in a single HDF5
			I/O operation. Filters are applied to those chunks of data. The rank
			of the chunkshape for tables must be 1. If None, a sensible value
			is calculated based on the expectedrows parameter (which is
			recommended).
		byteorder : str
			The byteorder of data *on disk*, specified as 'little' or 'big'.
			If this is not specified, the byteorder is that of the platform,
			unless you passed an array as the description, in which case
			its byteorder will be used.
		createparents : bool
			Whether to create the needed groups for the parent path to exist
			(not done by default).

		See Also
		--------
		Table : for more information on tables

		"""
		if rowDefinition is None:
			rowDefinition = self.rowDefinition
		if description is None:	#2012.12.20 rowDefinition is backup of description
			description = rowDefinition
		parentNode = self._getOrCreatePath(where, createparents)
		_checkfilters(filters)
		tableObject =  YHTable(parentNode=parentNode, tableName=tableName,
					 description=description, title=title,
					 filters=filters, expectedrows=expectedrows,
					 chunkshape=chunkshape, byteorder=byteorder)
		if tableObject._v_pathname not in self.tablePath2Index:
			self._appendNewTable(tableObject)
		else:
			sys.stderr.write("Error: table path name %s already in this PyTables file.\n"%(tableObject._v_pathname))
			sys.exit(3)
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
		if self.tableName:
			tableName = self.tableName
		else:
			if not tableNamePrefix:
				tableNamePrefix = self.tableNamePrefix
			i = len(groupObject._v_leaves)
			tableName = "%s%s"%(tableNamePrefix, i)
		while groupObject.__contains__(tableName):	#stop until a unique name shows up
			i += 1
			tableName = "%s%s"%(tableNamePrefix, i)
		return tableName
	
	def _readInData(self, tableName=None, tableObject=None):
		"""
		2012.12.16
		"""
		pass
	
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
	"""
	def constructColName2IndexFromHeader(self):
		#2012.11.22 overwrite parent function
		return self.combinedColID2ColIndex
	
	def getColIndexGivenColHeader(self, colHeader=None):
		#2012.11.15
		#	this is from the combined column header list.
		return self.combinedColID2ColIndex.get(colHeader)
	"""
	
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
		if tableName is None and self.tableName:
			tableName = self.tableName
		if tableIndex is None and tableName is None and self.tableObjectList:	#return the 1st one.
			return self.tableObjectList[0]
		groupObject = self._getGroupObject(groupIndex=None, groupName=groupName)
		tableObject = None	
		if tableName:
			tablePath = groupObject._v_pathname + '/' + tableName
			if tablePath in self.tablePath2Index:	#already cached
				return self._fetchTableObject(tablePath)
			else:
				tableObject = groupObject._f_getChild(tableName)
		else:
			tableObject = groupObject._v_leaves.values()[0]	#first leaf
		if not hasattr(tableObject, 'getAttribute'):
			tableObject = armTableObjectWithYHTableMethods(tableObject)
			self._appendNewTable(tableObject)
		return tableObject
	
	@property
	def nrows(self):
		"""
		2012.12.21
		"""
		tableObject = self.getTableObject()
		return tableObject.nrows
	
	def __iter__(self):
		"""
		2012.12.23 don't just "return self" and remove the next() because using "yield" is making a generator itself
			i.e.
			for row in yhFile:
				print row
			
			Python first calls obj=yhFile.__iter__(), then every "for" step, it calls obj.next().
			 
			In this __iter__(), "yield something" returns a newObj = iter(something) (AKA a generator), the in each "for" step, python calls newObj.next(). 
		2012.12.23 have to overwrite tables.File's __iter__() otherwise it's returning hdf5 group objects (/, /group0, etc.) iteratively
		return self
	
	def next(self):
		2012.12.20 iterate through just the first table
		
		"""
		
		row = None
		pdata = PassingDataList()
		tableObject = self.getTableObject(tableName=self.tableName)
		for row in tableObject:
			try:
				yield row	#yield row returns a generator
			except:
				raise StopIteration
		"""
		#2012.12.16
		#	go through each leaf (table & array), iteratively.
		#	usually the file has only one leaf (=table)
		for groupObject in self.walkGroups("/"):
			for leafNode in  groupObject._f_walkNodes('Leaf'):
				for row in leafNode:
					yield row
		"""
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
	
	def appendOneRow(self, oneCell=None, tableIndex=None, tableName=None, cellType=1):
		"""
		2013.1.9
		"""
		self.writeOneCell(oneCell=oneCell, tableIndex=tableIndex, tableName=tableName, cellType=cellType)
	
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

def castPyTablesRowIntoPassingData(rowPointer=None):
	"""
	2012.12.21 rowPointer from PyTables iteration is like a C pointer to the current row (no real content).
		need to convert it to a real object if you try to store its content in memory and use it later. 
	
	rowPointer has these methods: 
		['__class__', '__contains__', '__delattr__', '__delitem__', '__doc__', '__format__', '__getattribute__', 
		'__getitem__', '__hash__', '__init__', '__iter__', '__new__', '__next__', '__pyx_vtable__', '__reduce__', 
		'__reduce_ex__', '__repr__', '__setattr__', '__setitem__', '__sizeof__', '__str__', '__subclasshook__', 
		'_fillCol', '_flushBufferedRows', '_flushModRows', '_getUnsavedNrows', '_iter', 'append', 'fetch_all_fields', 
		'next', 'nrow', 'table', 'update']
		
		However, errors like this crop up during copying of these rowPointer. 

  File "/usr/lib/python2.7/copy_reg.py", line 93, in __newobj__
    return cls.__new__(cls, *args)
  File "tableExtension.pyx", line 706, in tables.tableExtension.Row.__cinit__ (tables/tableExtension.c:6910)
TypeError: __cinit__() takes exactly 1 positional argument (0 given)

	2012.12.21 could not use PassingDataList. because of these errors.

  File "/usr/lib/python2.7/copy.py", line 257, in _deepcopy_dict
    y[deepcopy(key, memo)] = deepcopy(value, memo)
  File "/usr/lib/python2.7/copy.py", line 182, in deepcopy
    rv = reductor(2)
TypeError: 'NoneType' object is not callable

	"""
	pdata = PassingData()	
	for colname in rowPointer.table.colnames:
		setattr(pdata, colname, rowPointer[colname])
	return pdata