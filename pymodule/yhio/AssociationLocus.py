#!/usr/bin/env python
"""
2012.12.15 data structure for AssociationLocus in PyTablesMatrixFile format
 
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
from AssociationPeak import AssociationPeakPyTable

class AssociationLocusPyTable(YHPyTable):
	id = UInt64Col(pos=0)
	chromosome = StringCol(64, pos=1)	#64 byte-long
	start = UInt64Col(pos=2)
	stop = UInt64Col(pos=3)
	
	no_of_peaks = UInt64Col(pos=3)
	connectivity = Float64Col(pos=4)
	no_of_results = UInt64Col(pos=5)
	phenotype_id_list = UInt64Col(shape=(1000,), pos=6)
	
	association_peak = AssociationPeakPyTable()
	
	
	def __init__(self, inputFname=None, openMode='r', \
				groupName=None, tableName='association_locus',\
				description=None,
				title='', filters=None, rowDefinition=None,\
				expectedrows=512000, genome_wide_result=None, **keywords):
		YHPyTable.__init__(self, inputFname=inputFname, openMode=openMode, \
				groupName=groupName, tableName=tableName,\
				description=description,
				title=title, filters=filters,
				expectedrows=expectedrows, **keywords)

		if openMode=='r':
			self._constructLandscapeGraph(tableName=self.tableName, genome_wide_result=self.genome_wide_result)
	
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
	
	def appendAssociationLoci(associationLocusList=None):
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
		sys.stderr.write("Dumping %s association loci into %s (HDF5 format) ..."%(len(associationLocusList), self.inputFname))
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
