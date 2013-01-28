#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s  -o /tmp/compare_association_locus_overlap.h5
		/Network/Data/250k/db/genome_wide_association_locus/type_4/24_251results_call57_ana32_type4.h5
		/Network/Data/250k/db/genome_wide_association_locus/type_4/4_251results_call75_ana32_type4.h5

Description:
	2013.1.26 This program finds correspondence between the two association loci from the input 
		and checks how different they are in terms of phenotypes in which they are significant.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import pylab
import csv, random, numpy
import tables
from tables import UInt64Col, Float64Col, StringCol, UInt32Col, Float32Col
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, getColName2IndexFromHeader, figureOutDelimiter,\
	yh_matplotlib, HDF5MatrixFile, YHFile
from pymodule import RBDict, CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
from pymodule import yhio, AssociationLocusTableFile, CNV
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
from TwoAssociationLocusFileOverlap import TwoAssociationLocusFileOverlap

class TwoGenomeWideAssociationLocusMapTable(tables.IsDescription):
	"""
	2013.1.26 further attributes associated with this table:
		input1_fname, input2_fname, 
		gw_association_locus1 (id, call-method, analysis-method-ls), gw_association_locus2 (id, cm, am)
		
		if "pos=.." is not added, they are sorted alphabetically by their names.
	"""
	id = UInt64Col(pos=0)
	
	input1_locus_id = UInt64Col()
	input1_chromosome = StringCol(64, )
	input1_start = UInt64Col()
	input1_stop = UInt64Col()
	
	input2_locus_id = UInt64Col()
	input2_chromosome = StringCol(64,)
	input2_start = UInt64Col()
	input2_stop = UInt64Col()
	
	locusOverlapFraction = Float64Col()
	
	no_of_total_phenotypes = UInt32Col()	#all significant phenotypes
	total_phenotype_ls_in_str = StringCol(1000)
	fraction_of_total_phenotypes = Float32Col()	#divided by all phenotypes with association
	
	no_of_overlap_phenotypes = UInt32Col()
	overlap_phenotype_ls_in_str = StringCol(1000)
	fraction_of_overlap_phenotypes = Float32Col()	#divided by the no_of_total_phenotypes (3 cols above, with significant hits)
	
	no_of_input1_only_phenotypes = UInt32Col()
	input1_only_phenotype_ls_in_str = StringCol(1000)
	fraction_of_input1_only_phenotypes = Float32Col()
	
	no_of_input2_only_phenotypes = UInt32Col()
	input2_only_phenotype_ls_in_str = StringCol(1000)
	fraction_of_input2_only_phenotypes = Float32Col()


class CompareTwoGWAssociationLocusByPhenotypeVector(TwoAssociationLocusFileOverlap):
	__doc__ = __doc__
	option_default_dict = TwoAssociationLocusFileOverlap.option_default_dict.copy()
	# change default of file-format
	option_default_dict[('inputFileFormat', 0, int)][0] = 2
	option_default_dict[('outputFileFormat', 0, int)][0] = 2
	option_default_dict.update({
					('locusPadding', 0, int): [0, '', 1, 'the padding around each locus (used only to extend the loci)' ],\
				})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		TwoAssociationLocusFileOverlap.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
		
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
			writer = YHFile(self.outputFname, openMode='w', rowDefinition=TwoGenomeWideAssociationLocusMapTable)
		self.writer = writer
	
	def constructLocusCorrespondanceBetweenTwoFiles(self):
		"""
		2013.1.27
		"""
	
	def reduce(self, **keywords):
		"""
		2013.1.27
		"""
		associationLocusRBDict1 = self.associationLocusFileList[0].associationLocusRBDict
		associationLocusRBDict2 = self.associationLocusFileList[1].associationLocusRBDict
		nodePairList = CNV.findCorrespondenceBetweenTwoCNVRBDict(associationLocusRBDict1, associationLocusRBDict2)
		
		
		#output attributes from its HDF5 file
		for attributeName in associationLocusRBDict1.HDF5AttributeNameLs:
			attributeValue = getattr(associationLocusRBDict1, attributeName)
			self.writer.addAttribute(name='association1_%s'%(attributeName), value=attributeValue)
		for attributeName in associationLocusRBDict2.HDF5AttributeNameLs:
			attributeValue = getattr(associationLocusRBDict2, attributeName)
			self.writer.addAttribute(name='association2_%s'%(attributeName), value=attributeValue)
		
		all_tested_phenotype_list = list(set(associationLocusRBDict1.phenotype_method_id_ls)|set(associationLocusRBDict2.phenotype_method_id_ls))
		
		for nodePair in nodePairList:
			input1Node, input2Node, overlapFraction = nodePair[:3]
			if input1Node is not None and input2Node is not None:
				input1_locus_id = input1Node.key.locus_id
				input1_chromosome = input1Node.key.chromosome
				input1_start = input1Node.key.start
				input1_stop = input1Node.key.stop
				
				input2_locus_id = input2Node.key.locus_id
				input2_chromosome = input2Node.key.chromosome
				input2_start = input2Node.key.start
				input2_stop = input2Node.key.stop
				
				total_phenotype_id_set = input1Node.key.phenotype_id_set | input2Node.key.phenotype_id_set
				no_of_total_phenotypes = len(total_phenotype_id_set)
				
				overlapping_phenotype_id_set = input1Node.key.phenotype_id_set & input2Node.key.phenotype_id_set
				no_of_overlap_phenotypes = len(overlapping_phenotype_id_set)
				
				input1_only_phenotype_id_set = input1Node.key.phenotype_id_set - input2Node.key.phenotype_id_set
				no_of_input1_only_phenotypes = len(input1_only_phenotype_id_set)
				input2_only_phenotype_id_set = input2Node.key.phenotype_id_set - input1Node.key.phenotype_id_set
				no_of_input2_only_phenotypes = len(input2_only_phenotype_id_set)
			elif input1Node is None and input2Node is not None:
				input1_locus_id = None
				input1_chromosome = None 
				input1_start = None
				input1_stop = None
				
				input2_locus_id = input2Node.key.locus_id
				input2_chromosome = input2Node.key.chromosome
				input2_start = input2Node.key.start
				input2_stop = input2Node.key.stop
				
				total_phenotype_id_set = input2Node.key.phenotype_id_set
				no_of_total_phenotypes = len(total_phenotype_id_set)
				
				overlapping_phenotype_id_set = set()
				no_of_overlap_phenotypes = len(overlapping_phenotype_id_set)
				
				input1_only_phenotype_id_set = set()
				no_of_input1_only_phenotypes = len(input1_only_phenotype_id_set)
				input2_only_phenotype_id_set = input2Node.key.phenotype_id_set
				no_of_input2_only_phenotypes = len(input2_only_phenotype_id_set)
				
			elif input1Node is not None and input2Node is None:
				input1_locus_id = input1Node.key.locus_id
				input1_chromosome = input1Node.key.chromosome
				input1_start = input1Node.key.start
				input1_stop =  input1Node.key.stop
				
				input2_locus_id = None
				input2_chromosome = None
				input2_start = None
				input2_stop = None
				
				total_phenotype_id_set = input1Node.key.phenotype_id_set
				no_of_total_phenotypes = len(total_phenotype_id_set)
				
				overlapping_phenotype_id_set = set()
				no_of_overlap_phenotypes = len(overlapping_phenotype_id_set)
				
				input1_only_phenotype_id_set = input1Node.key.phenotype_id_set
				no_of_input1_only_phenotypes = len(input1_only_phenotype_id_set)
				input2_only_phenotype_id_set = set()
				no_of_input2_only_phenotypes = len(input2_only_phenotype_id_set)
			else:
				sys.stderr.write("Error: Both input1Node and input2Node are None.\n")
				sys.exit(3)
			oneCell = PassingData(
					input1_locus_id = input1_locus_id,
					input1_chromosome = input1_chromosome,
					input1_start = input1_start,
					input1_stop = input1_stop,
					
					input2_locus_id = input2_locus_id,
					input2_chromosome = input2_chromosome,
					input2_start = input2_start,
					input2_stop = input2_stop,
					
					locusOverlapFraction = overlapFraction,
					no_of_total_phenotypes = no_of_total_phenotypes,
					total_phenotype_ls_in_str = utils.getSuccinctStrOutOfList(list(total_phenotype_id_set)),
					fraction_of_total_phenotypes = float(len(total_phenotype_id_set))/len(all_tested_phenotype_list),
					
					no_of_overlap_phenotypes = no_of_overlap_phenotypes,
					overlap_phenotype_ls_in_str = utils.getSuccinctStrOutOfList(list(overlapping_phenotype_id_set)),
					fraction_of_overlap_phenotypes = float(no_of_overlap_phenotypes)/no_of_total_phenotypes,
					
					no_of_input1_only_phenotypes = no_of_input1_only_phenotypes,
					input1_only_phenotype_ls_in_str = utils.getSuccinctStrOutOfList(list(input1_only_phenotype_id_set)),
					fraction_of_input1_only_phenotypes = float(no_of_input1_only_phenotypes)/no_of_total_phenotypes,
					
					no_of_input2_only_phenotypes = no_of_input2_only_phenotypes,
					input2_only_phenotype_ls_in_str = utils.getSuccinctStrOutOfList(list(input2_only_phenotype_id_set)),
					fraction_of_input2_only_phenotypes = float(no_of_input2_only_phenotypes)/no_of_total_phenotypes)
			
			self.writer.writeOneCell(oneCell, cellType=2)
		

if __name__ == '__main__':
	main_class = CompareTwoGWAssociationLocusByPhenotypeVector
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
