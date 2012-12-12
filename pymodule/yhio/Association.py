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
from pymodule.yhio.HDF5MatrixFile import HDF5MatrixFile, addAttributeDictToHDF5GroupObject
from pymodule.yhio.SNP import getGenomeWideResultFromHDF5MatrixFile


def getAssociationLandscapeDataFromHDF5File(inputFname=None, associationGroupName='association', \
										landscapeGroupName='landscape', min_MAF=0.1):
	"""
	2012.11.20
		input is in HDF5MatrixFile format (which is output of variation/src/association_peak/DefineAssociationLandscape.py)
	"""
	pdata = PassingData(min_MAF=min_MAF)
	genome_wide_result = getGenomeWideResultFromHDF5MatrixFile(inputFname=inputFname, \
						min_value_cutoff=None, do_log10_transformation=False, pdata=pdata,\
						construct_chr_pos2index=False, construct_data_obj_id2index=False, \
						construct_locus_db_id2index=True,\
						report=True, groupName=associationGroupName)
	
	returnData = PassingData(genome_wide_result=genome_wide_result)
	
	sys.stderr.write("Reading landscape from %s ..."%(inputFname))
	current_obj = None
	bridge_ls = []
	locusLandscapeNeighborGraph = nx.Graph()
	reader = HDF5MatrixFile(inputFname, openMode='r')
	landscapeGroupObject = reader.getGroupObject(groupName=landscapeGroupName)
	returnData.HDF5AttributeNameLs = []
	for attributeName, value in landscapeGroupObject.getAttributes().iteritems():
		returnData.HDF5AttributeNameLs.append(attributeName)
		setattr(returnData, attributeName, value)
	
	for row in landscapeGroupObject:
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

def outputAssociationLandscapeInHDF5(bridge_ls=None, outputFname=None, writer=None, closeFile=False, groupName='landscape',\
							attributeDict=None,):
	"""
	2012.11.18
	"""
	sys.stderr.write("Outputting the %s bridges from the landscape ..."%(len(bridge_ls)))
	#output the data_object.id in bridge_ls to outputFname
	#each number below is counting bytes, not bits
	dtypeList = [('start_locus_id','i8'),('stop_locus_id', 'i8'), ('no_of_loci','i8'), ('deltaX', 'i8')]
	if writer:
		groupObject = writer.createNewGroup(groupName=groupName, dtypeList=dtypeList)
	elif outputFname:
		writer = HDF5MatrixFile(outputFname, openMode='w', dtypeList=dtypeList, firstGroupName=groupName)
		groupObject = writer.getGroupObject(groupName=groupName)
	else:
		sys.stderr.write("Error: no writer(%s) or filename(%s) to dump.\n"%(writer, filename))
		sys.exit(3)
	addAttributeDictToHDF5GroupObject(groupObject=groupObject, attributeDict=attributeDict)
	
	previous_locus_id = None
	cellList = []
	for bridge in bridge_ls:
		current_obj = bridge[0]
		obj_with_fastest_score_increase = bridge[1]
		no_of_loci, deltaX = bridge[2:4]
		dataTuple = (current_obj.db_id, obj_with_fastest_score_increase.db_id, no_of_loci, deltaX)
		cellList.append(dataTuple)
	groupObject.writeCellList(cellList)
	if closeFile:
		writer.close()
	sys.stderr.write("%s objects.\n"%(len(cellList)))
	return writer

def constructAssociationPeakRBDictFromHDF5File(inputFname=None, peakPadding=10000, groupName='association_peak'):
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
	
	groupObject = reader.getGroupObject(groupName=groupName)
	for attributeName, value in groupObject.getAttributes().iteritems():
		associationPeakRBDict.HDF5AttributeNameLs.append(attributeName)
		setattr(associationPeakRBDict, attributeName, value)
	
	counter = 0
	real_counter = 0
	for row in groupObject:
		if not row.chromosome:	#empty chromosome, which happens when inputFname contains no valid peaks, but the default null peak (only one).
			continue
		counter += 1
		segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=row.chromosome, \
						span_ls=[max(1, row.start - peakPadding), row.stop + peakPadding], \
						min_reciprocal_overlap=1, result_peak_id=None)
						#2010-8-17 overlapping keys are regarded as separate instances as long as they are not identical.
		if segmentKey not in associationPeakRBDict:
			associationPeakRBDict[segmentKey] = []
		associationPeakRBDict[segmentKey].append(row)
	sys.stderr.write("%s peaks in %s spans.\n"%(counter, len(associationPeakRBDict)))
	return associationPeakRBDict


def outputAssociationPeakInHDF5(association_peak_ls=None, filename=None, writer=None, groupName='association_peak', closeFile=True,\
							attributeDict=None,):
	"""
	2012.11.20
	"""
	sys.stderr.write("Dumping association peaks into %s (HDF5 format) ..."%(filename))
	#each number below is counting bytes, not bits
	dtypeList = [('chromosome', HDF5MatrixFile.varLenStrType), ('start','i8'), ('stop', 'i8'), \
				('start_locus_id','i8'), ('stop_locus_id','i8'), \
				('no_of_loci', 'i8'), ('peak_locus_id', 'i8'), ('peak_score', 'f8')]
	if writer is None and filename:
		writer = HDF5MatrixFile(filename, openMode='w', dtypeList=dtypeList, firstGroupName=groupName)
		groupObject = writer.getGroupObject(groupName=groupName)
	elif writer:
		groupObject = writer.createNewGroup(groupName=groupName, dtypeList=dtypeList)
	else:
		sys.stderr.write("Error: no writer(%s) or filename(%s) to dump.\n"%(writer, filename))
		sys.exit(3)
	#add neighbor_distance, max_neighbor_distance, min_MAF, min_score, ground_score as attributes
	addAttributeDictToHDF5GroupObject(groupObject=groupObject, attributeDict=attributeDict)
	cellList = []
	#2012.11.28 sort it
	association_peak_ls.sort()
	for association_peak in association_peak_ls:
		dataTuple = (association_peak.chromosome, association_peak.start, association_peak.stop, \
					association_peak.start_locus_id, association_peak.stop_locus_id, \
					association_peak.no_of_loci,\
					association_peak.peak_locus_id, association_peak.peak_score)
		cellList.append(dataTuple)
	
	if groupObject is None:
		sys.stderr.write("Error: groupObject (name=%s) is None. could not write.\n"%(groupName))
		sys.exit(3)
	groupObject.writeCellList(cellList)
	if closeFile:
		writer.close()
	sys.stderr.write("%s objects.\n"%(len(cellList)))
	return writer


def constructAssociationLocusRBDictFromHDF5File(inputFname=None, locusPadding=0, groupName='association_locus'):
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
	groupObject = reader.getGroupObject(groupName=groupName)
	for attributeName, value in groupObject.getAttributes().iteritems():
		associationLocusRBDict.HDF5AttributeNameLs.append(attributeName)
		setattr(associationLocusRBDict, attributeName, value)
	
	counter = 0
	real_counter = 0
	for row in groupObject:
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

def outputAssociationLociInHDF5(associationLocusList=None, filename=None, writer=None, groupName='association_locus', \
					closeFile=True,\
					attributeDict=None):
	"""
	2012.12.10
		for each locus, output the association peaks that fall into the locus.
			for each association peak, include result-id, phenotype id, peak-score, start-locus,
				chromosome
				start
				stop
				start_locus
				stop_locus
				no_of_loci
				peak_locus
				peak_score
	2012.11.20
	"""
	sys.stderr.write("Dumping association loci into %s (HDF5 format) ..."%(filename))
	#each number below is counting bytes, not bits
	dtypeList = [('chromosome', HDF5MatrixFile.varLenStrType), \
				('start','i8'), ('stop', 'i8'), \
				('no_of_peaks', 'i8'), ('connectivity', 'f8'), ('no_of_results', 'i8')]
	if writer is None and filename:
		writer = HDF5MatrixFile(filename, openMode='w', dtypeList=dtypeList, firstGroupName=groupName)
		groupObject = writer.getGroupObject(groupName=groupName)
	elif writer:
		groupObject = writer.createNewGroup(groupName=groupName, dtypeList=dtypeList)
	else:
		sys.stderr.write("Error: no writer(%s) or filename(%s) to dump.\n"%(writer, filename))
		sys.exit(3)
	#add neighbor_distance, max_neighbor_distance, min_MAF, min_score, ground_score as attributes
	addAttributeDictToHDF5GroupObject(groupObject=groupObject, attributeDict=attributeDict)
	cellList = []
	#2012.11.28 sort it
	associationLocusList.sort()
	for associationLocus in associationLocusList:
		dataTuple = (associationLocus.chromosome, associationLocus.start, associationLocus.stop, associationLocus.no_of_peaks,\
					associationLocus.connectivity, associationLocus.no_of_results)
		cellList.append(dataTuple)
	
	if groupObject is None:
		sys.stderr.write("Error: groupObject (name=%s) is None. could not write.\n"%(groupName))
		sys.exit(3)
	groupObject.writeCellList(cellList)
	if closeFile:
		writer.close()
	sys.stderr.write("%s objects.\n"%(len(cellList)))
	return writer