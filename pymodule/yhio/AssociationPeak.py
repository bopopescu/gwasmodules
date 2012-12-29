#!/usr/bin/env python
"""
2012.12.15 data structure for AssociationPeak in PyTablesMatrixFile format
 
"""
import sys, os, math
__doc__ = __doc__%()

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import csv
import tables
from tables import UInt64Col, Float64Col, StringCol
import numpy
from pymodule.utils import PassingData, PassingDataList
from pymodule.ProcessOptions import ProcessOptions
from pymodule.yhio.YHPyTables import YHTable, YHFile, castPyTablesRowIntoPassingData

class AssociationPeakTable(tables.IsDescription):
	id = UInt64Col(pos=0)
	chromosome = StringCol(64, pos=1)	#64 byte-long
	start = UInt64Col(pos=2)
	stop = UInt64Col(pos=3)
	start_locus_id = UInt64Col(pos=4)
	stop_locus_id = UInt64Col(pos=5)
	no_of_loci = UInt64Col(pos=6)
	peak_locus_id = UInt64Col(pos=7)
	peak_score = Float64Col(pos=8)

class AssociationPeakTableFile(YHFile):

	"""
	usage examples:
	
		peakFile = AssociationPeakTableFile(self.outputFname, openMode='w')
		peakFile.addAttributeDict(attributeDict)
		peakFile.appendAssociationPeak(association_peak_ls=association_peak_ls)
		
		#for read-only
		peakFile = AssociationPeakTableFile(inputFname, openMode='r', peakPadding=0)
		rbDict = peakFile.associationPeakRBDict
	"""
	def __init__(self, inputFname=None, openMode='r', \
				tableName='association_peak', groupNamePrefix='group', tableNamePrefix='table',\
				filters=None, peakPadding=0, \
				**keywords):
		
		YHFile.__init__(self, inputFname=inputFname, openMode=openMode, \
				tableName=tableName, groupNamePrefix=groupNamePrefix, tableNamePrefix=tableNamePrefix,\
				rowDefinition=None, filters=filters, debug=0, report=0)
		
		self.peakPadding = peakPadding
		self.associationPeakRBDict = None
		if openMode=='r':
			self.associationPeakTable = self.getTableObject()
			self._constructAssociationPeakRBDict(tableObject=self.associationPeakTable)
		elif openMode == 'w':
			self.associationPeakTable = self.createNewTable(tableName=self.tableName, rowDefinition=AssociationPeakTable,\
												expectedrows=50000)

	def _constructAssociationPeakRBDict(self, tableObject=None):
		"""
		2012.11.12
			similar to Stock_250kDB.constructRBDictFromResultPeak(), but from HDF5MatrixFile-like file
		"""
		from pymodule.algorithm.RBTree import RBDict
		from pymodule.yhio.CNV import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
		if tableObject is None:
			tableObject = self.associationPeakTable
		sys.stderr.write("Constructing association-peak RBDict from HDF5 file %s, (peakPadding=%s) ..."%(self.inputFname, self.peakPadding))
		associationPeakRBDict = RBDict()
		associationPeakRBDict.result_id = None	#2012.6.22
		associationPeakRBDict.peakPadding = self.peakPadding
		associationPeakRBDict.HDF5AttributeNameLs = []
		
		for attributeName, value in self.getAttributes().iteritems():
			associationPeakRBDict.HDF5AttributeNameLs.append(attributeName)
			setattr(associationPeakRBDict, attributeName, value)
		
		counter = 0
		real_counter = 0
		for row in tableObject:
			if not row['chromosome']:	#empty chromosome, which happens when inputFname contains no valid peaks, but the default null peak (only one).
				continue
			counter += 1
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=row['chromosome'], \
							span_ls=[max(1, row['start'] - self.peakPadding), row['stop'] + self.peakPadding], \
							min_reciprocal_overlap=1, result_peak_id=None)
							#2010-8-17 overlapping keys are regarded as separate instances as long as they are not identical.
			if segmentKey not in associationPeakRBDict:
				associationPeakRBDict[segmentKey] = []
			else:
				sys.stderr.write("Warning: segmentKey of %s already in associationPeakRBDict with this row: %s.\n"%\
								(row, associationPeakRBDict[segmentKey][0]))
			associationPeakRBDict[segmentKey].append(castPyTablesRowIntoPassingData(row))	#row is a pointer to the current row.
		sys.stderr.write("%s peaks in %s spans.\n"%(counter, len(associationPeakRBDict)))
		
		self.associationPeakRBDict = associationPeakRBDict
		return self.associationPeakRBDict
	
	def appendAssociationPeak(self, association_peak_ls=None):
		"""
		2012.11.20
		"""
		sys.stderr.write("Dumping %s association peaks into %s ..."%(len(association_peak_ls), self.inputFname))
		#each number below is counting bytes, not bits
		cellList = []
		#2012.11.28 sort it
		association_peak_ls.sort()
		for association_peak in association_peak_ls:
			dataTuple = (association_peak.chromosome, association_peak.start, association_peak.stop, \
						association_peak.start_locus_id, association_peak.stop_locus_id, \
						association_peak.no_of_loci,\
						association_peak.peak_locus_id, association_peak.peak_score)
			self.associationPeakTable.writeOneCell(dataTuple)
		self.flush()
		sys.stderr.write(" \n")
	
