#!/usr/bin/env python
"""
Examples:
	%s 
	
	
	%s  -w GenomeDBPassword -F /Network/Data/250k/tmp-yh/SNPInfo_LocusType1.pickle -l 1 -u yh

Description:
	2012.3.16
		find the overlapping peaks (form a connected component) among the sets. => mark their plot filenames that they form one cc.
		for each cc, find out its combined genomic span, whether a candidate gene is near/within this span.
		near-candidate and non-near-candidate cc's have their plots separate.
		  cc.size()>1 => two gwas have overlapping peaks
		output: regionOfInterestFile for each span (including multiple peaks, give 10kb on each side for the span.), 
			make chr, start, stop, phenotype_id, list of peak_id, candidate gene id/name,
				candidate or not in the beginning of the filename, 
				whether overlap between two different call methods:
				overlap_candidate_..., overlap_noncandidate_..., onlyCall32_candidate..., onlyCall80_noncandidate..., 
			part of output plot filename 
		how to identify the overlapping peaks? check whether one edge within cc is connecting two different call methods.
		the output file will be fed to two Drawing jobs with each taking different call method id.

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])
import matplotlib; matplotlib.use("Agg")

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from AbstractVariationMapper import AbstractVariationMapper
import h5py, numpy
from pymodule import GenomeDB
from pymodule import SNP 
from variation.src import Stock_250kDB
from variation.src.TwoGWASPeakOverlap import TwoGWASPeakOverlap 
import networkx as nx

class PickleSNPInfo(AbstractVariationMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVariationMapper.option_default_dict.copy()
	option_default_dict.update({
							('snpInfoPickleFname', 1, ): ['', 'F', 1, 'The file to contain pickled SNPInfo.'],\
							('call_method_id', 0, int): [None, 'm', 1, "provide if locus_type_id is not."],\
							('locus_type_id', 0, int): [None, 'l', 1, ''],\
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		AbstractVariationMapper.__init__(self, inputFnameLs, **keywords)

	def run(self):
		"""
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_250k = self.db_250k
		
		if not self.locus_type_id and self.call_method_id:
			rm = Stock_250kDB.ResultsMethod.query.filter_by(call_method_id=self.call_method_id).first()
			self.locus_type_id = rm.locus_type_id
		if not self.locus_type_id:
			sys.stderr.write("self.locus_type_id %s is still not valid. Quit.\n"%(self.locus_type_id))
			sys.exit(3)
		
		db_250k.dealWithSNPInfo(self.snpInfoPickleFname, locus_type_id=self.locus_type_id)
	

if __name__ == '__main__':
	main_class = PickleSNPInfo
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()