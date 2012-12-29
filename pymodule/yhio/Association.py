#!/usr/bin/env python
"""
2012.11.20 data structures for IO
"""
import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import numpy
import networkx as nx
from pymodule.utils import PassingData
from pymodule.yhio.HDF5MatrixFile import HDF5MatrixFile, addAttributeDictToYHTableInHDF5Group
from pymodule.yhio.YHPyTables import YHTable, YHFile
from pymodule.yhio.SNP import getGenomeWideResultFromHDF5MatrixFile

import tables
from tables import UInt64Col, Float64Col, StringCol

class AssociationTable(tables.IsDescription):
	"""
	2012.12.18 pytable class to store the genome-wide association result
	"""
	id = UInt64Col(pos=0)
	locus_id = UInt64Col(pos=1)
	chromosome = StringCol(64, pos=2)	#64 byte-long
	start = UInt64Col(pos=3)
	stop = UInt64Col(pos=4)
	score = Float64Col(pos=5)
	MAC = UInt64Col(pos=6)
	MAF = Float64Col(pos=7)
	genotype_var_perc = Float64Col(pos=8)

class AssociationTableFile(YHFile):
	#no beta0, beta1, beta2
	def __init__(self, inputFname=None, openMode='r', \
				tableName='association', groupNamePrefix='group', tableNamePrefix='table',\
				filters=None,\
				min_MAF=0.1, **keywords):
		YHFile.__init__(self, inputFname=inputFname, openMode=openMode, \
				tableName=tableName, groupNamePrefix=groupNamePrefix, tableNamePrefix=tableNamePrefix,\
				rowDefinition=None, filters=filters, debug=0, report=0)
		
		self.min_MAF = min_MAF
		self.genome_wide_result = None
		self.associationTable = None
		if openMode=='r':
			self.associationTable = self.getTableObject(tableName=self.tableName)
			self._readInGWR(min_MAF=self.min_MAF, tableObject=self.associationTable)
		elif openMode=='w':
			self.associationTable = self.createNewTable(tableName=self.tableName, rowDefinition=AssociationTable, \
											expectedrows=300000)
	
	def _readInGWR(self, inputFname=None, tableName=None, min_MAF=None, tableObject=None):
		"""
		"""
		if inputFname is None:
			inputFname = self.inputFname
		if tableName is None:
			tableName = self.tableName
		if min_MAF is None:
			min_MAF = self.min_MAF
		pdata = PassingData(min_MAF=min_MAF)
		self.genome_wide_result = getGenomeWideResultFromHDF5MatrixFile(reader=self, tableName=tableName, tableObject=tableObject,\
							min_value_cutoff=None, do_log10_transformation=False, pdata=pdata,\
							construct_chr_pos2index=False, construct_data_obj_id2index=False, \
							construct_locus_db_id2index=True,\
							report=True)
		return self.genome_wide_result

def getAssociationLandscapeDataFromHDF5File(inputFname=None, associationTableName='association', \
										landscapeTableName='landscape', min_MAF=0.1):
	"""
	2012.11.20
		input is in HDF5MatrixFile format (which is output of variation/src/association_peak/DefineAssociationLandscape.py)
		contains two hdf5 groups. one is by associationTableName. the other is by landscapeTableName.
	"""
	pdata = PassingData(min_MAF=min_MAF)
	genome_wide_result = getGenomeWideResultFromHDF5MatrixFile(inputFname=inputFname, \
						min_value_cutoff=None, do_log10_transformation=False, pdata=pdata,\
						construct_chr_pos2index=False, construct_data_obj_id2index=False, \
						construct_locus_db_id2index=True,\
						report=True, tableName=associationTableName)
	
	returnData = PassingData(genome_wide_result=genome_wide_result)
	
	sys.stderr.write("Reading landscape from %s ..."%(inputFname))
	current_obj = None
	bridge_ls = []
	locusLandscapeNeighborGraph = nx.Graph()
	reader = HDF5MatrixFile(inputFname, openMode='r')
	landscapeTableObject = reader.getTableObject(tableName=landscapeTableName)
	returnData.HDF5AttributeNameLs = []
	for attributeName, value in landscapeTableObject.getAttributes().iteritems():
		returnData.HDF5AttributeNameLs.append(attributeName)
		setattr(returnData, attributeName, value)
	
	for row in landscapeTableObject:
		if row.start_locus_id==0:	#empty data. happens when inputFname contains no valid landscape, but one default null data point.
			continue
		start_locus_id = row.start_locus_id
		stop_locus_id = row.stop_locus_id
		no_of_loci = row.no_of_loci
		deltaX = row.deltaX
		
		start_obj = genome_wide_result.get_data_obj_by_locus_db_id(start_locus_id)
		stop_obj = genome_wide_result.get_data_obj_by_locus_db_id(stop_locus_id)
		
		bridge_ls.append([start_obj, stop_obj, no_of_loci, deltaX])
		
		source_index = start_obj.index
		#genome_wide_result.get_data_obj_index_by_locus_db_id(start_locus_id)
		target_index = stop_obj.index
		
		locusLandscapeNeighborGraph.add_edge(source_index, target_index, \
									weight=None)
		locusLandscapeNeighborGraph[source_index][target_index]['no_of_loci'] = no_of_loci
		locusLandscapeNeighborGraph[source_index][target_index]['deltaX'] = deltaX
		
	del reader
	sys.stderr.write("%s bridges.\n"%(len(bridge_ls)))
	returnData.bridge_ls = bridge_ls
	returnData.locusLandscapeNeighborGraph = locusLandscapeNeighborGraph
	return returnData

def outputAssociationLandscapeInHDF5(bridge_ls=None, outputFname=None, writer=None, closeFile=False, tableName='landscape',\
							attributeDict=None,):
	"""
	2012.11.18
	"""
	sys.stderr.write("Outputting the %s bridges from the landscape ..."%(len(bridge_ls)))
	#output the data_object.id in bridge_ls to outputFname
	#each number below is counting bytes, not bits
	rowDefinition = [('start_locus_id','i8'),('stop_locus_id', 'i8'), ('no_of_loci','i8'), ('deltaX', 'i8')]
	if writer:
		tableObject = writer.createNewTable(tableName=tableName, rowDefinition=rowDefinition)
	elif outputFname:
		writer = HDF5MatrixFile(outputFname, openMode='w', rowDefinition=rowDefinition, tableName=tableName)
		tableObject = writer.getTableObject(tableName=tableName)
	else:
		sys.stderr.write("Error: no writer(%s) or filename(%s) to dump.\n"%(writer, filename))
		sys.exit(3)
	addAttributeDictToYHTableInHDF5Group(tableObject=tableObject, attributeDict=attributeDict)
	
	previous_locus_id = None
	cellList = []
	for bridge in bridge_ls:
		current_obj = bridge[0]
		obj_with_fastest_score_increase = bridge[1]
		no_of_loci, deltaX = bridge[2:4]
		dataTuple = (current_obj.db_id, obj_with_fastest_score_increase.db_id, no_of_loci, deltaX)
		cellList.append(dataTuple)
	tableObject.writeCellList(cellList)
	if closeFile:
		writer.close()
	sys.stderr.write("%s objects.\n"%(len(cellList)))
	return writer

def constructAssociationPeakRBDictFromHDF5File(inputFname=None, peakPadding=10000, tableName='association_peak'):
	"""
	2012.11.12
		similar to Stock_250kDB.constructRBDictFromResultPeak(), but from HDF5MatrixFile-like file
	"""
	from pymodule.algorithm.RBTree import RBDict
	from pymodule.yhio.CNV import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
	
	sys.stderr.write("Constructing association-peak RBDict from HDF5 file %s, (peakPadding=%s) ..."%(inputFname, peakPadding))
	reader = HDF5MatrixFile(inputFname, openMode='r')
	associationPeakRBDict = RBDict()
	associationPeakRBDict.result_id = None	#2012.6.22
	associationPeakRBDict.peakPadding = peakPadding
	associationPeakRBDict.HDF5AttributeNameLs = []
	
	tableObject = reader.getTableObject(tableName=tableName)
	for attributeName, value in tableObject.getAttributes().iteritems():
		associationPeakRBDict.HDF5AttributeNameLs.append(attributeName)
		setattr(associationPeakRBDict, attributeName, value)
	
	counter = 0
	real_counter = 0
	for row in tableObject:
		if not row.chromosome:	#empty chromosome, which happens when inputFname contains no valid peaks, but the default null peak (only one).
			continue
		counter += 1
		segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=row.chromosome, \
						span_ls=[max(1, row.start - peakPadding), row.stop + peakPadding], \
						min_reciprocal_overlap=1, result_peak_id=None)
						#2010-8-17 overlapping keys are regarded as separate instances as long as they are not identical.
		if segmentKey not in associationPeakRBDict:
			associationPeakRBDict[segmentKey] = []
		else:
			sys.stderr.write("Warning: segmentKey of %s already in associationPeakRBDict with this row: %s.\n"%\
							(row, associationPeakRBDict[segmentKey][0]))
		associationPeakRBDict[segmentKey].append(row)
	sys.stderr.write("%s peaks in %s spans.\n"%(counter, len(associationPeakRBDict)))
	return associationPeakRBDict


def outputAssociationPeakInHDF5(association_peak_ls=None, filename=None, writer=None, tableName='association_peak', closeFile=True,\
							attributeDict=None,):
	"""
	2012.11.20
	"""
	sys.stderr.write("Dumping association peaks into %s (HDF5 format) ..."%(filename))
	#each number below is counting bytes, not bits
	rowDefinition = [('chromosome', HDF5MatrixFile.varLenStrType), ('start','i8'), ('stop', 'i8'), \
				('start_locus_id','i8'), ('stop_locus_id','i8'), \
				('no_of_loci', 'i8'), ('peak_locus_id', 'i8'), ('peak_score', 'f8')]
	if writer is None and filename:
		writer = HDF5MatrixFile(filename, openMode='w', rowDefinition=rowDefinition, tableName=tableName)
		tableObject = writer.getTableObject(tableName=tableName)
	elif writer:
		tableObject = writer.createNewTable(tableName=tableName, rowDefinition=rowDefinition)
	else:
		sys.stderr.write("Error: no writer(%s) or filename(%s) to dump.\n"%(writer, filename))
		sys.exit(3)
	#add neighbor_distance, max_neighbor_distance, min_MAF, min_score, ground_score as attributes
	addAttributeDictToYHTableInHDF5Group(tableObject=tableObject, attributeDict=attributeDict)
	cellList = []
	#2012.11.28 sort it
	association_peak_ls.sort()
	for association_peak in association_peak_ls:
		dataTuple = (association_peak.chromosome, association_peak.start, association_peak.stop, \
					association_peak.start_locus_id, association_peak.stop_locus_id, \
					association_peak.no_of_loci,\
					association_peak.peak_locus_id, association_peak.peak_score)
		cellList.append(dataTuple)
	
	if tableObject is None:
		sys.stderr.write("Error: tableObject (name=%s) is None. could not write.\n"%(tableName))
		sys.exit(3)
	tableObject.writeCellList(cellList)
	if closeFile:
		writer.close()
	sys.stderr.write("%s objects.\n"%(len(cellList)))
	return writer


def constructAssociationLocusRBDictFromHDF5File(inputFname=None, locusPadding=0, tableName='association_locus'):
	"""
	2012.11.25
		similar to constructAssociationPeakRBDictFromHDF5File
	"""
	from pymodule.algorithm.RBTree import RBDict
	from pymodule.yhio.CNV import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
	
	sys.stderr.write("Constructing association-locus RBDict from HDF5 file %s, (locusPadding=%s) ..."%(inputFname, locusPadding))
	reader = HDF5MatrixFile(inputFname, openMode='r')
	associationLocusRBDict = RBDict()
	associationLocusRBDict.locusPadding = locusPadding
	associationLocusRBDict.HDF5AttributeNameLs = []
	tableObject = reader.getTableObject(tableName=tableName)
	for attributeName, value in tableObject.getAttributes().iteritems():
		associationLocusRBDict.HDF5AttributeNameLs.append(attributeName)
		setattr(associationLocusRBDict, attributeName, value)
	
	counter = 0
	real_counter = 0
	for row in tableObject:
		if not row.chromosome:	#empty chromosome, which happens when inputFname contains no valid locus, but the default null locus (only one).
			continue
		counter += 1
		segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=row.chromosome, \
						span_ls=[max(1, row.start - locusPadding), row.stop + locusPadding], \
						min_reciprocal_overlap=1, no_of_peaks=row.no_of_peaks, \
						no_of_results=row.no_of_results, connectivity=row.connectivity)
						#2010-8-17 overlapping keys are regarded as separate instances as long as they are not identical.
		if segmentKey not in associationLocusRBDict:
			associationLocusRBDict[segmentKey] = []
		associationLocusRBDict[segmentKey].append(row)
	sys.stderr.write("%s peaks in %s spans.\n"%(counter, len(associationLocusRBDict)))
	return associationLocusRBDict

def outputAssociationLociInHDF5(associationLocusList=None, filename=None, writer=None, tableName='association_locus', \
					closeFile=True,\
					attributeDict=None):
	"""
	2012.12.10
		for each locus, output the association peaks that fall into the locus.
			for each association peak, include 
				* result-id 
				* phenotype id
				* chromosome
				* start
				* stop
				* start_locus
				* stop_locus
				* no_of_loci
				* peak_locus
				* peak-score
	2012.11.20
	"""
	sys.stderr.write("Dumping association loci into %s (HDF5 format) ..."%(filename))
	#each number below is counting bytes, not bits
	rowDefinition = [('chromosome', HDF5MatrixFile.varLenStrType), \
				('start','i8'), ('stop', 'i8'), \
				('no_of_peaks', 'i8'), ('connectivity', 'f8'), ('no_of_results', 'i8')]
	if writer is None and filename:
		writer = HDF5MatrixFile(filename, openMode='w', rowDefinition=rowDefinition, tableName=tableName)
		tableObject = writer.getTableObject(tableName=tableName)
	elif writer:
		tableObject = writer.createNewTable(tableName=tableName, rowDefinition=rowDefinition)
	else:
		sys.stderr.write("Error: no writer(%s) or filename(%s) to dump.\n"%(writer, filename))
		sys.exit(3)
	#add neighbor_distance, max_neighbor_distance, min_MAF, min_score, ground_score as attributes
	addAttributeDictToYHTableInHDF5Group(tableObject=tableObject, attributeDict=attributeDict)
	cellList = []
	#2012.11.28 sort it
	associationLocusList.sort()
	for associationLocus in associationLocusList:
		dataTuple = (associationLocus.chromosome, associationLocus.start, associationLocus.stop, associationLocus.no_of_peaks,\
					associationLocus.connectivity, associationLocus.no_of_results)
		cellList.append(dataTuple)
	
	if tableObject is None:
		sys.stderr.write("Error: tableObject (name=%s) is None. could not write.\n"%(tableName))
		sys.exit(3)
	tableObject.writeCellList(cellList)
	if closeFile:
		writer.close()
	sys.stderr.write("%s objects.\n"%(len(cellList)))
	return writer
