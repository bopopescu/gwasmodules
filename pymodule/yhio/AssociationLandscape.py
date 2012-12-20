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
from tables import *
import numpy
from pymodule.utils import PassingData, PassingDataList
from pymodule.ProcessOptions import  ProcessOptions
from pymodule.yhio.YHPyTable import YHPyTable
import networkx as nx

class AssociationLandscapePyTable(YHPyTable):
	id = UInt64Col(pos=0)
	start_locus_id = UInt64Col(pos=1)
	stop_locus_id = UInt64Col(pos=2)
	no_of_loci = UInt64Col(pos=3)
	deltaX = UInt64Col(pos=4)
	deltaY = Float64Col(pos=5)
	
	"""
	2012.12.18 usage examples:
		
		#for writing
		landscapeFile = AssociationLandscapePyTable(self.outputFname, openMode='w')
		landscapeFile.addAttributeDict(attributeDict)
		landscapeFile.appendAssociationLandscapeBridgeList(bridge_ls=landscapeData.bridge_ls)
		
		#for read-only, genome_wide_result is required if locusLandscapeNeighborGraph is needed.
		landscapeTable = AssociationLandscapePyTable(self.inputFname, openMode='r', \
							tableName='landscape', genome_wide_result=associationTable.genome_wide_result)
	"""
	
	def __init__(self, inputFname=None, openMode='r', \
				groupName=None, tableName='landscape',\
				description=None,
				title='', filters=None, rowDefinition=None,\
				expectedrows=512000, genome_wide_result=None, **keywords):
		YHPyTable.__init__(self, inputFname=inputFname, openMode=openMode, \
				groupName=groupName, tableName=tableName,\
				description=description,
				title=title, filters=filters,
				expectedrows=expectedrows, **keywords)

		self.genome_wide_result = genome_wide_result
		self.bridge_ls = None
		self.locusLandscapeNeighborGraph = None
		if openMode=='r' and self.genome_wide_result is not None:
			self._constructLandscapeGraph(tableName=self.tableName, genome_wide_result=self.genome_wide_result)
	
	def _constructLandscapeGraph(self, tableName=None, genome_wide_result=None):
		"""
		2012.12.17
		"""
		if tableName is None:
			tableName = self.tableName
		if genome_wide_result is None:
			genome_wide_result = self.genome_wide_result
		
		sys.stderr.write("Reading landscape from table %s ..."%(tableName))
		
		current_obj = None
		self.bridge_ls = []
		self.locusLandscapeNeighborGraph = nx.Graph()
		self.HDF5AttributeNameLs = []
		landscapeTableObject = self.getTableObject(tableName=tableName)
		
		"""
		for attributeName, value in landscapeTableObject.getAttributes().iteritems():
			self.HDF5AttributeNameLs.append(attributeName)
			setattr(returnData, attributeName, value)
		"""
		for row in landscapeTableObject:
			start_locus_id = row[self.getColIndex('start_locus_id')]
			if start_locus_id==0:	#empty data. happens when inputFname contains no valid landscape, but one default null data point.
				continue
			stop_locus_id = row[self.getColIndex('stop_locus_id')]
			no_of_loci = row[self.getColIndex('no_of_loci')]
			deltaX = row[self.getColIndex('deltaX')]
			deltaY = row[self.getColIndex('deltaY')]
			
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
		sys.stderr.write("Outputting the %s bridges from the landscape ..."%(len(bridge_ls)))
		#output the data_object.id in bridge_ls to outputFname
		#each number below is counting bytes, not bits
		previous_locus_id = None
		cellList = []
		for bridge in bridge_ls:
			current_obj = bridge[0]
			obj_with_fastest_score_increase = bridge[1]
			no_of_loci, deltaX, deltaY = bridge[2:5]
			dataTuple = (current_obj.db_id, obj_with_fastest_score_increase.db_id, no_of_loci, deltaX, deltaY)
			self.writeOneCell(dataTuple)
			#cellList.append(dataTuple)
		#tableObject.writeCellList(cellList)
		self.flush()
		sys.stderr.write(" \n")
