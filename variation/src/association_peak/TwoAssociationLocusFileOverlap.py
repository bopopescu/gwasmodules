#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s  -o /tmp/two_association_locus_overlap.h5 /tmp/association_locus_overlap_9.h5 /tmp/association_locus_overlap_5.h5

Description:
	2012.11.22 This program calculates the distance of association loci from 1st input to 2nd input or
		overlapping fraction if loci overlap.
	If "-i ..." is given, it is regarded as one of the input files (plus the ones in trailing arguments).
	The total number of input should be 2. If more than 2, only the first 2 are used.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import pylab
import csv, random, numpy
import tables
from tables import UInt64Col, Float64Col, StringCol
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, getColName2IndexFromHeader, figureOutDelimiter,\
	yh_matplotlib, HDF5MatrixFile, YHFile
from pymodule import RBDict, CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
from pymodule import yhio, AssociationLocusTableFile
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
from CountAssociationLocus import CountAssociationLocus

class TwoAssociationLocusOverlapTable(tables.IsDescription):
	"""
	2012.12.24 new PyTables-based table definition
	"""
	id = UInt64Col(pos=0)
	chromosome = StringCol(64, pos=1)
	start = UInt64Col(pos=2)
	stop = UInt64Col(pos=3)
	fractionCoveredByAssociation2 = Float64Col(pos=4)

class TwoAssociationLocusFileOverlap(CountAssociationLocus):
	__doc__ = __doc__
	option_default_dict = CountAssociationLocus.option_default_dict.copy()
	# change default of file-format
	option_default_dict[('inputFileFormat', 0, int)][0] = 2
	option_default_dict[('outputFileFormat', 0, int)][0] = 2
	#option_default_dict.update(AbstractMapper.db_option_dict.copy())
	option_default_dict.update({
					('locusPadding', 0, int): [0, '', 1, 'the padding around each locus (used only to extend the loci)' ],\
				})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		CountAssociationLocus.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		#self.connectDB() called within its __init__()
		self.associationLocusRBDictList = []
	
	def setup(self, **keywords):
		"""
		2012.11.25 setup the output
		"""
		writer = None
		if self.outputFileFormat==1:
			suffix = os.path.splitext(self.outputFname)[1]
			if self.outputFname and suffix!='.png':
				writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		else:	#HDF5MatrixFile
			writer = YHFile(self.outputFname, openMode='w', rowDefinition=TwoAssociationLocusOverlapTable)
		self.writer = writer
		
	def fileWalker(self, inputFname=None, afterFileFunction=None, processRowFunction=None , run_type=1, **keywords):
		"""
		2012.11.22
		"""
		sys.stderr.write("walking through %s ..."%(inputFname))
		counter = 0
		real_counter = 0
		noOfSampled = 0
		pdata = self.initiatePassingData()
		if afterFileFunction is None:
			afterFileFunction = self.afterFileFunction
		try:
			associationLocusTableFile = AssociationLocusTableFile(inputFname, openMode='r')
			rbDict = associationLocusTableFile.associationLocusRBDict
			#rbDict = yhio.Association.constructAssociationLocusRBDictFromHDF5File(inputFname=inputFname, \
			#							locusPadding=self.locusPadding, tableName='association_locus')
			self.associationLocusRBDictList.append(rbDict)
			associationLocusTableFile.close()
			
			counter += 1
			real_counter += 1
			noOfSampled += 1
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
			sys.exit(3)
		if counter>0:
			fraction = float(noOfSampled)/float(counter)
		else:
			fraction = 0
		sys.stderr.write("%s/%s (%.3f) data sampled. real_counter=%s.\n"%(noOfSampled, counter, fraction, real_counter))
	
	def reduce(self, **keywords):
		"""
		2012.11.25 derive overlap & distance-to-closest-locus statistics between two input and output the statistics as well.
		"""
		compareIns = CNVCompare(min_reciprocal_overlap=0.0000001)	#to detect any overlap
		#for i in xrange(len(self.associationLocusRBDictList)):
		#	for j in xrange(i+1, len(self.associationLocusRBDictList)):
		associationLocusRBDict1 = self.associationLocusRBDictList[0]
		associationLocusRBDict2 = self.associationLocusRBDictList[1]
		#output attributes from its HDF5 file
		for attributeName in associationLocusRBDict1.HDF5AttributeNameLs:
			attributeValue = getattr(associationLocusRBDict1, attributeName)
			self.writer.addAttribute(name='association1_%s'%(attributeName), value=attributeValue)
		for attributeName in associationLocusRBDict2.HDF5AttributeNameLs:
			attributeValue = getattr(associationLocusRBDict2, attributeName)
			self.writer.addAttribute(name='association2_%s'%(attributeName), value=attributeValue)
		
		cellList = []
		for locus1Node in associationLocusRBDict1:
			targetNodeLs = []
			associationLocusRBDict2.findNodes(locus1Node.key, node_ls=targetNodeLs, compareIns=compareIns)
			fractionCoveredByAssociation2 = 0.0
			if targetNodeLs:
				for locus2Node in targetNodeLs:
					overlap1 = get_overlap_ratio(locus1Node.key.span_ls, locus2Node.key.span_ls)[0]
					fractionCoveredByAssociation2 += overlap1
			else:
				closestNodeData = associationLocusRBDict2.findClosestNode(key=locus1Node.key)
				distanceToSmallerNode = None
				distanceToBiggerNode = None
				if closestNodeData is not None:
					smallerNode = closestNodeData.smallerNode
					biggerNode = closestNodeData.biggerNode
					if smallerNode == biggerNode:
						sys.stderr.write("Warning: should not happen. Found no overlapping nodes for %s but now found one, %s.\n"%\
										(repr(locus1Node), repr(smallerNode)))
						sys.exit(4)
					else:
						if smallerNode and smallerNode.key.chromosome==locus1Node.key.chromosome:	#make sure it's same chromosome
							distanceToSmallerNode = locus1Node.key.start - smallerNode.key.stop
						if biggerNode and biggerNode.key.chromosome==locus1Node.key.chromosome:
							distanceToBiggerNode = biggerNode.key.start - locus1Node.key.stop
				else:
					sys.stderr.write("Error: impossible to happen. Found no close nodes for %s (%s:%s-%s).\n"%\
									(repr(locus1Node.key), locus1Node.key.chromosome, locus1Node.key.start, \
									locus1Node.key.stop))
					#continue	#2012.12.3 temporary bugfix to skip loci with empty chromosomes.
					# bug from upstream AssociationPeak2AssociationLocus.py
					sys.exit(6)
				if distanceToBiggerNode is not None and distanceToSmallerNode is not None:
					fractionCoveredByAssociation2 = min(distanceToBiggerNode, distanceToSmallerNode)
				elif distanceToBiggerNode is not None:
					fractionCoveredByAssociation2 = distanceToBiggerNode
				elif distanceToSmallerNode is not None:
					fractionCoveredByAssociation2 = distanceToSmallerNode
				else:	#on a new chromosome
					fractionCoveredByAssociation2 = -1	#
				
				#if fractionCoveredByAssociation2<-10:
				#	import pdb
				#	pdb.set_trace()
			oneCell = (locus1Node.key.chromosome, locus1Node.key.start, locus1Node.key.stop, float(fractionCoveredByAssociation2))
			cellList.append(oneCell)
		self.writer.writeCellList(cellList)
		
		if getattr(self, 'writer', None):
			del self.writer

if __name__ == '__main__':
	main_class = TwoAssociationLocusFileOverlap
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()