#!/usr/bin/env python
"""
2012.12.15 data structure for AssociationLandscape in PyTablesMatrixFile format
 
"""
import sys, os, math
__doc__ = __doc__%()

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import csv
import tables
import numpy
import networkx as nx
from tables import UInt64Col, Float64Col
from pymodule.utils import PassingData, PassingDataList
from pymodule.ProcessOptions import  ProcessOptions
from pymodule.yhio.YHPyTables import YHTable, YHFile
from Association import AssociationTable, AssociationTableFile

class AssociationLandscapeTable(tables.IsDescription):
	id = UInt64Col(pos=0)
	start_locus_id = UInt64Col(pos=1)
	stop_locus_id = UInt64Col(pos=2)
	no_of_loci = UInt64Col(pos=3)
	deltaX = UInt64Col(pos=4)
	deltaY = Float64Col(pos=5)
	#no beta0, beta1, beta2

	
class AssociationLandscapeTableFile(AssociationTableFile):
	"""
	2012.12.18 usage examples:
		
		#for writing
		landscapeFile = AssociationLandscapePyTable(self.outputFname, openMode='w')
		landscapeFile.addAttributeDict(attributeDict)
		landscapeFile.appendAssociationLandscapeBridgeList(bridge_ls=landscapeData.bridge_ls)
		
		#for read-only,
		landscapeTable = AssociationLandscapePyTable(self.inputFname, openMode='r')
	"""
	def __init__(self, inputFname=None, openMode='r', \
				tableName='association_landscape', groupNamePrefix='group', tableNamePrefix='table',\
				filters=None,\
				min_MAF=0.1, associationTableName='association', **keywords):
		
		rowDefinition = AssociationLandscapeTable
		YHFile.__init__(self, inputFname=inputFname, openMode=openMode, \
				tableName=tableName, groupNamePrefix=groupNamePrefix, tableNamePrefix=tableNamePrefix,\
				rowDefinition=rowDefinition, filters=filters, debug=0, report=0)

		self.associationTableName = associationTableName
		self.min_MAF = min_MAF
		
		self.bridge_ls = None
		self.locusLandscapeNeighborGraph = None
		if openMode=='r':
			self.associationLandscapeTable = self.getTableObject(tableName=self.tableName)
			self.associationTable = self.getTableObject(tableName=self.associationTableName)
			self.genome_wide_result = self._readInGWR(min_MAF=self.min_MAF, tableObject=self.associationTable)
			self._constructLandscapeGraph(genome_wide_result=self.genome_wide_result)
		elif openMode =='w':
			self.associationLandscapeTable = self.createNewTable(tableName=self.tableName, rowDefinition=AssociationLandscapeTable,\
														expectedrows=50000)
			self.associationTable = self.createNewTable(tableName=self.associationTableName, rowDefinition=AssociationTable,\
													expectedrows=300000)
	
	def _constructLandscapeGraph(self, genome_wide_result=None, tableObject=None):
		"""
		2012.12.17
		"""
		if genome_wide_result is None:
			genome_wide_result = self.genome_wide_result
		if tableObject is None:
			tableObject = self.associationLandscapeTable
		
		sys.stderr.write("Reading landscape from table %s ..."%(self.tableName))
		
		current_obj = None
		self.bridge_ls = []
		self.locusLandscapeNeighborGraph = nx.Graph()
		self.HDF5AttributeNameLs = []
		#landscapeTableObject = self.getTableObject(tableName=tableName)
		
		"""
		for attributeName, value in landscapeTableObject.getAttributes().iteritems():
			self.HDF5AttributeNameLs.append(attributeName)
			setattr(returnData, attributeName, value)
		"""
		for row in tableObject:
			start_locus_id = row['start_locus_id']	#fastest accessing
			#start_locus_id = row[self.associationLandscapeTable.getColIndex('start_locus_id')]	#2nd fastest,
			#start_locus_id = row[self.getColIndex('start_locus_id')]	#slowest
			if start_locus_id==0:	#empty data. happens when inputFname contains no valid landscape, but one default null data point.
				continue
			#fastest accessing
			stop_locus_id = row['stop_locus_id']
			no_of_loci = row['no_of_loci']
			deltaX = row['deltaX']
			deltaY = row['deltaY']
			
			"""
			#2nd fastest, one layer in between
			stop_locus_id = row[self.associationLandscapeTable.getColIndex('stop_locus_id')]
			no_of_loci = row[self.associationLandscapeTable.getColIndex('no_of_loci')]
			deltaX = row[self.associationLandscapeTable.getColIndex('deltaX')]
			deltaY = row[self.associationLandscapeTable.getColIndex('deltaY')]
			#should be slowest
			stop_locus_id = row[self.getColIndex('stop_locus_id')]
			no_of_loci = row[self.getColIndex('no_of_loci')]
			deltaX = row[self.getColIndex('deltaX')]
			deltaY = row[self.getColIndex('deltaY')]
			"""
			
			start_obj = genome_wide_result.get_data_obj_by_locus_db_id(start_locus_id)
			stop_obj = genome_wide_result.get_data_obj_by_locus_db_id(stop_locus_id)
			
			self.bridge_ls.append([start_obj, stop_obj, no_of_loci, deltaX, deltaY])
			
			source_index = start_obj.index
			#genome_wide_result.get_data_obj_index_by_locus_db_id(start_locus_id)
			target_index = stop_obj.index
			
			self.locusLandscapeNeighborGraph.add_edge(source_index, target_index, \
										weight=None)
			self.locusLandscapeNeighborGraph[source_index][target_index]['no_of_loci'] = no_of_loci
			self.locusLandscapeNeighborGraph[source_index][target_index]['deltaX'] = deltaX
			self.locusLandscapeNeighborGraph[source_index][target_index]['deltaY'] = deltaY
		
		sys.stderr.write("%s bridges.\n"%(len(self.bridge_ls)))
	
	def appendAssociationLandscapeBridgeList(self, bridge_ls=None):
		"""
		2012.11.18
		"""
		sys.stderr.write("Outputting the %s bridges from the landscape to %s ..."%(len(bridge_ls), self.inputFname))
		#output the data_object.id in bridge_ls to outputFname
		#each number below is counting bytes, not bits
		previous_locus_id = None
		cellList = []
		for bridge in bridge_ls:
			current_obj = bridge[0]
			obj_with_fastest_score_increase = bridge[1]
			no_of_loci, deltaX, deltaY = bridge[2:5]
			dataTuple = (current_obj.db_id, obj_with_fastest_score_increase.db_id, no_of_loci, deltaX, deltaY)
			self.associationLandscapeTable.writeOneCell(dataTuple)
			#cellList.append(dataTuple)
		#tableObject.writeCellList(cellList)
		self.flush()
		sys.stderr.write(" \n")
