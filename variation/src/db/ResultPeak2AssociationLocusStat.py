#!/usr/bin/env python
"""
Examples:
	#2012.9.30
	%s  -i 1-35,39-48,57-82,158-159,161-179,182-186,272-283,314-351,362-380,418-589 -l 80,75 -u yh -z banyan -s 2,3,1,3 -c
	-o ~/script/variation/data/AssociationPeakStat.tsv

	#2012.9.25
	%s  -i 1-35,39-48,57-82,158-159,161-179,182-186,272-274,277-283,314-351,362-380 -l 32,80,57,75 -u yh -z banyan -s 2,3,1,3 -c

Description:
	2012.6.1
		program to combine association result peaks into association locus and its stats

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule import SNPData, SNP
#from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
from variation.src.mapper.AbstractVariationMapper import AbstractVariationMapper
from variation.src import Stock_250kDB
import networkx as nx
from pymodule.CNV import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
from pymodule.RBTree import RBDict
import numpy

class ResultPeak2AssociationLocusStat(AbstractVariationMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVariationMapper.option_default_dict.copy()
	option_default_dict.update({
							('phenotype_id_ls', 1, ): ['', 'i', 1, 'loci from these peaks will be outputted', ],\
							('call_method_id_ls', 1, ): ['', 'l', 1, 'loci from these peaks will be outputted', ],\
							('analysis_method_id_ls', 1, ): ['1,32', 'a', 1, 'loci from these peaks will be outputted', ],\
							('min_overlap_ratio', 1, float): ['0.05', 'm', 1, 'minimum overlap ratio, overlap length/total' ],\
							('result_peak_type_id_ls', 1, ): ['', 's', 1, 'list of result peak type IDs, in this order:\
								snp-KW-gwas,snp-EMMA-gwas,deletion-KW-gwas,deletion-EMMA-gwas', ],\
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		AbstractVariationMapper.__init__(self, inputFnameLs, **keywords)
		self.phenotype_id_ls = getListOutOfStr(self.phenotype_id_ls, data_type=int)
		self.phenotype_id_ls.sort()
		self.result_peak_type_id_ls = getListOutOfStr(self.result_peak_type_id_ls, data_type=int)
		self.call_method_id_ls = getListOutOfStr(self.call_method_id_ls, data_type=int)
		self.analysis_method_id_ls = getListOutOfStr(self.analysis_method_id_ls, data_type=int)
		
	
	def constructGraph(self, resultPeakRBDictList=[], min_overlap_ratio=0.1):
		"""
		2012.6.24
		"""
		sys.stderr.write("Constructing graph between all result peaks ...")
		rpg = resultPeakGraph = nx.Graph()
		n = len(resultPeakRBDictList)
		counter = 0
		for i in xrange(n):
			for j in xrange(i+1, n):
				peakRBDict1 = resultPeakRBDictList[i]
				peakRBDict2 = resultPeakRBDictList[j]
				compareIns = CNVCompare(min_reciprocal_overlap=min_overlap_ratio)
				for node in peakRBDict1:
					counter += 1
					segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=node.key.chromosome, \
									span_ls=[node.key.start, node.key.stop], \
									min_reciprocal_overlap=0.0000001, )	#min_reciprocal_overlap doesn't matter here.
						# it's decided by compareIns.
					node_ls = []
					peakRBDict2.findNodes(segmentKey, node_ls=node_ls, compareIns=compareIns)
					total_perc_overlapped_by_result2 = 0.
					rpg.add_node(node.key.result_peak_id, span=segmentKey.span_ls)
					for result2_node in node_ls:
						if result2_node.key.result_peak_id not in rpg:
							rpg.add_node(result2_node.key.result_peak_id, span = [result2_node.key.start, result2_node.key.stop])
						
						rpg.add_edge(node.key.result_peak_id, result2_node.key.result_peak_id)
		sys.stderr.write("%s nodes. %s edges. %s connected components.\n"%(rpg.number_of_nodes(), rpg.number_of_edges(), \
															nx.number_connected_components(rpg) ) )
		return rpg
	
	def discoverAssociationLocus(self, db_250k=None, resultPeakGraph=None, min_overlap_ratio=0.1):
		"""
		2012.6.24
		"""
		sys.stderr.write("Discovering association loci ...")
		cc_graph_list = nx.connected_component_subgraphs(resultPeakGraph)
		counter = 0
		for cc_graph in cc_graph_list:
			#calculate connectivity of this component
			ne = cc_graph.number_of_edges()
			nn = cc_graph.number_of_nodes()
			if nn>1:
				connectivity = ne/float(nn*(nn-1)/2)
			else:
				connectivity = 1
			start_ls = []
			stop_ls = []
			result_peak_ls = []
			#get span of each node, then take median of all its start/stop
			for n in cc_graph:
				span = resultPeakGraph.node[n]['span']
				start_ls.append(span[0])
				stop_ls.append(span[1])
				result_peak = Stock_250kDB.ResultPeak.get(n)
				result_peak_ls.append(result_peak)
			median_start = numpy.median(start_ls)
			median_stop = numpy.median(stop_ls)
			no_of_peaks = nn
			
			association_locus = db_250k.getAssociationLocus(chromosome=result_peak.chromosome, start=median_start, stop=median_stop, \
										no_of_peaks=nn, connectivity=connectivity,\
						threshold=min_overlap_ratio, result_peak_ls=result_peak_ls)
			counter += 1
		db_250k.session.flush()
		sys.stderr.write("%s association loci added into db.\n"%(counter))
		
	
	def run(self):
		"""
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_250k = self.db_250k
		
		db_250k.session.begin()
		
		result_ls = db_250k.getResultLs(call_method_id=None, analysis_method_id_ls=self.analysis_method_id_ls, \
								phenotype_method_id_ls=self.phenotype_id_ls, \
								call_method_id_ls=self.call_method_id_ls, cnv_method_id=None)
		
		
		
		resultPeakRBDictList =[]
		
		noCorrectionAnalysisMethodIDSet = set([1,2,3,5,6,16])
		EMMAAnalysisMethodIDSet = set([4,7,8,32,33])
		
		for result in result_ls:
			if result.call_method.locus_type_id ==1 and result.analysis_method_id in noCorrectionAnalysisMethodIDSet:	#KW on snp
				result_peak_type_id = self.result_peak_type_id_ls[0]
			elif result.call_method.locus_type_id==1 and result.analysis_method_id in EMMAAnalysisMethodIDSet:	#EMMA/EMMAX on snp
				result_peak_type_id = self.result_peak_type_id_ls[1]
			elif result.call_method.locus_type_id==2 and result.analysis_method_id in noCorrectionAnalysisMethodIDSet:	#KW on deletion
				result_peak_type_id = self.result_peak_type_id_ls[2]
			elif result.call_method.locus_type_id==2 and result.analysis_method_id in EMMAAnalysisMethodIDSet:	#EMMA on deletion 
				result_peak_type_id = self.result_peak_type_id_ls[3]
			else:
				sys.stderr.write("Error: ")
				sys.exit(2)
				result_peak_type_id = None
			rbDict = db_250k.constructRBDictFromResultPeak(result.id, result_peak_type_id, peakPadding=0)
			resultPeakRBDictList.append(rbDict)
		resultPeakGraph = self.constructGraph(resultPeakRBDictList=resultPeakRBDictList, \
									min_overlap_ratio=self.min_overlap_ratio)
		
		self.discoverAssociationLocus(db_250k=db_250k, resultPeakGraph=resultPeakGraph, min_overlap_ratio=self.min_overlap_ratio)
		
		if self.commit:
			db_250k.session.flush()
			db_250k.session.commit()
		else:
			db_250k.session.rollback()


if __name__ == '__main__':
	main_class = ResultPeak2AssociationLocusStat
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()