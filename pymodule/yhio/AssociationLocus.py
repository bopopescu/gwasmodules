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
from tables import UInt64Col, StringCol, Float64Col
import numpy
from pymodule.utils import PassingData, PassingDataList
from pymodule.ProcessOptions import  ProcessOptions
from pymodule.yhio.YHPyTables import YHTable, YHFile, castPyTablesRowIntoPassingData
from pymodule.algorithm.RBTree import RBDict
from pymodule.yhio.CNV import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
from AssociationPeak import AssociationPeakTable

class AssociationLocusTable(tables.IsDescription):
	"""
	2013.1.28 bugfix, was pos=3 for no_of_peaks (same as stop), now change it to pos=4, and increment others accordingly
	"""
	id = UInt64Col(pos=0)
	chromosome = StringCol(64, pos=1)	#64 byte-long
	start = UInt64Col(pos=2)
	stop = UInt64Col(pos=3)
	no_of_peaks = UInt64Col(pos=4)
	connectivity = Float64Col(pos=5)
	no_of_results = UInt64Col(pos=6)
	phenotype_id_ls_in_str = StringCol(1000, pos=7)

class AssociationLocus2PeakTable(tables.IsDescription):
	id = UInt64Col(pos=0)
	association_locus_id = UInt64Col(pos=1)
	association_peak = AssociationPeakTable()


class AssociationLocusTableFile(YHFile):

	"""
	usage examples:
	
		associationLocusTableFile = AssociationLocusTableFile(self.outputFname, openMode='w')
		associationLocusTableFile.addAttributeDict(attributeDict)
		associationLocusTableFile.appendAssociationPeak(association_peak_ls=association_peak_ls)
		
		#for read-only
		associationLocusTableFile = AssociationLocusTableFile(inputFname, openMode='r')
		rbDict = associationLocusTableFile.associationLocusRBDict
		
		associationLocusTableFile = AssociationLocusTableFile(inputFname, openMode='r', constructLocusRBDict=False)	#don't need it
		call_method_id_ls = associationLocusTableFile.getAttribute('call_method_id_ls')
	"""
	def __init__(self, inputFname=None, openMode='r', \
				tableName='association_locus', groupNamePrefix='group', tableNamePrefix='table',\
				filters=None, autoRead=True, autoWrite=True, \
				locus2PeakTableName='association_locus2peak', locusPadding=0, constructLocusRBDict=True,\
				**keywords):
		
		self.constructLocusRBDict = constructLocusRBDict
		self.locus2PeakTableName = locus2PeakTableName
		self.locusPadding = locusPadding
		self.associationLocusRBDict = None
		
		YHFile.__init__(self, inputFname=inputFname, openMode=openMode, \
				tableName=tableName, groupNamePrefix=groupNamePrefix, tableNamePrefix=tableNamePrefix,\
				rowDefinition=None, filters=filters, debug=0, report=0,\
				autoRead=False, autoWrite=False)
		
		#to overwrite self.autoRead that is set by YHFile.__init__
		self.autoRead = autoRead
		self.autoWrite = autoWrite
		
		if self.autoRead and (self.openMode=='r' or self.openMode=='a'):
			self.associationLocusTable = self.getTableObject(tableName=self.tableName)
			self.associationLocus2PeakTable = self.getTableObject(tableName=self.locus2PeakTableName)
			if self.constructLocusRBDict:
				self.associationLocusRBDict = self._readInData(tableName=self.tableName, tableObject=self.associationLocusTable)
		elif openMode == 'w':
			self.associationLocusTable = self.createNewTable(tableName=self.tableName, rowDefinition=AssociationLocusTable,\
													expectedrows=50000)
			self.associationLocus2PeakTable = self.createNewTable(tableName=self.locus2PeakTableName, \
													rowDefinition=AssociationLocus2PeakTable, expectedrows=500000)
	
	def _readInData(self, tableName=None, tableObject=None, bugfixType=None):
		"""
		2013.1.28 added argument bugfixType (default is None)
			1: swap stop & no_of_peaks, an earlier bug exchanged the positions of the two.
		2013.1.26 added phenotype_id_set in the node
		2012.11.25
			similar to constructAssociationPeakRBDictFromHDF5File
		"""
		if tableName is None:
			tableName = self.tableName
		YHFile._readInData(self, tableName=tableName, tableObject=tableObject)
		if not self.constructLocusRBDict:
			return
		
		locusPadding = self.locusPadding
		sys.stderr.write("Constructing association-locus RBDict (locusPadding=%s) ..."%(locusPadding))
		if tableObject is None:
			tableObject = self.getTableObject(tableName=tableName)
		associationLocusRBDict = RBDict()
		associationLocusRBDict.locusPadding = locusPadding
		associationLocusRBDict.HDF5AttributeNameLs = []
		
		for attributeName, value in tableObject.getAttributes().iteritems():
			associationLocusRBDict.HDF5AttributeNameLs.append(attributeName)
			setattr(associationLocusRBDict, attributeName, value)
		
		counter = 0
		real_counter = 0
		for rowPointer in tableObject:
			row = castPyTablesRowIntoPassingData(rowPointer)
			if not row.chromosome:	#empty chromosome, which happens when inputFname contains no valid locus, but the default null locus (only one).
				continue
			counter += 1
			phenotype_id_ls = row.phenotype_id_ls_in_str.split(',')
			phenotype_id_set = set(map(int, phenotype_id_ls))
			if bugfixType==1:
				#2013.1.28 old association-loci file have two columns swapped. run this to correct it.
				# a function in variation/src/misc.py is written:
				#	DB250k.correctAssociationLocusFileFormat(db_250k=db_250k, data_dir=None)
				rowPointer['stop'] = row.no_of_peaks
				rowPointer['no_of_peaks'] = row.stop
				rowPointer.update()
				row.no_of_peaks = rowPointer['no_of_peaks']
				row.stop = rowPointer['stop']
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=row.chromosome, \
							span_ls=[max(1, row.start - locusPadding), row.stop + locusPadding], \
							min_reciprocal_overlap=1, no_of_peaks=row.no_of_peaks, \
							no_of_results=row.no_of_results, connectivity=row.connectivity,\
							phenotype_id_set=phenotype_id_set, locus_id=row.id)
							#2010-8-17 overlapping keys are regarded as separate instances as long as they are not identical.
			if segmentKey not in associationLocusRBDict:
				associationLocusRBDict[segmentKey] = []
			associationLocusRBDict[segmentKey].append(row)
		sys.stderr.write("%s peaks in %s spans.\n"%(counter, len(associationLocusRBDict)))
		self.associationLocusRBDict = associationLocusRBDict
		return associationLocusRBDict
	
	def appendAssociationLoci(self, associationLocusList=None):
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
		sys.stderr.write("Saving %s association loci into file %s ..."%(len(associationLocusList), self.inputFname))
		#add neighbor_distance, max_neighbor_distance, min_MAF, min_score, ground_score as attributes
		#2012.11.28 sort it
		associationLocusList.sort()
		for associationLocus in associationLocusList:
			#dataTuple = (associationLocus.chromosome, associationLocus.start, associationLocus.stop, associationLocus.no_of_peaks,\
			#			associationLocus.connectivity, associationLocus.no_of_results, associationLocus.phenotype_id_ls_in_str)
			self.associationLocusTable.writeOneCell(associationLocus, cellType=2)
			associationLocus2PeakTableRow = self.associationLocus2PeakTable.row
			for association_peak in associationLocus.association_peak_ls:
				self.associationLocus2PeakTable.no_of_rows += 1
				associationLocus2PeakTableRow['id'] = self.associationLocus2PeakTable.no_of_rows
				associationLocus2PeakTableRow['association_locus_id'] = self.associationLocusTable.no_of_rows
				associationLocus2PeakTableRow['association_peak/id'] = association_peak['id']
				associationLocus2PeakTableRow['association_peak/chromosome'] = association_peak['chromosome']
				associationLocus2PeakTableRow['association_peak/start'] = association_peak['start']
				associationLocus2PeakTableRow['association_peak/stop'] = association_peak['stop']
				associationLocus2PeakTableRow['association_peak/start_locus_id'] = association_peak['start_locus_id']
				associationLocus2PeakTableRow['association_peak/stop_locus_id'] = association_peak['stop_locus_id']
				associationLocus2PeakTableRow['association_peak/no_of_loci'] = association_peak['no_of_loci']
				associationLocus2PeakTableRow['association_peak/peak_locus_id'] = association_peak['peak_locus_id']
				associationLocus2PeakTableRow['association_peak/peak_score'] = association_peak['peak_score']
				associationLocus2PeakTableRow.append()
				#self.associationLocus2PeakTable.writeOneCell(oneCell, cellType=3)	#doesn't work. nested type shows up as one column in colnames.
		sys.stderr.write("\n")
	
	def fetchOneLocusGivenID(self, locus_id=None):
		"""
		2013.1.27 query the file by the AssociationLocus.id
		"""
		rows = self.readWhere('(id=%s)'%(locus_id))
		if len(rows)>1:
			sys.stderr.write("Warning: %s (>1) rows found with id=%s. Return 1st one.\n"%(len(rows), locus_id))
		elif len(rows)==0:
			#sys.stderr.write("Warning: nothing found with id=%s.\n"%(locus_id))
			return None
		return rows[0]
	
	def fetchOneLocusGivenChrStartStop(self, chromosome=None, start=None, stop=None):
		"""
		2013.1.27
		"""
		rows = self.readWhere('(chromosome=%s) & (start=%s) & (stop=%s)'%(chromosome, start, stop))
		if len(rows)>1:
			sys.stderr.write("Warning: %s (>1) rows found with chr=%s, start=%s, stop=%s. Return 1st one.\n"%\
							(len(rows), chromosome, start, stop))
		elif len(rows)==0:
			#sys.stderr.write("Warning: nothing found with id=%s.\n"%(locus_id))
			return None
		return rows[0]
