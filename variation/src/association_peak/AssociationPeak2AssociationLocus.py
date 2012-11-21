#!/usr/bin/env python
"""
Examples:
	#2012.9.30
	%s  -i 1-35,39-48,57-82,158-159,161-179,182-186,272-283,314-351,362-380,418-589 -l 80,75 -u yh -z banyan -s 2,3,1,3 -c
	-o ~/script/variation/data/AssociationPeakStat.tsv

	#2012.11.20
	%s -o /tmp/association_locus.h5 -m 0.05 /tmp/5566_association_peak.h5 /tmp/5566_association_peak.h5

Description:
	2012.6.1
		program to combine association result peaks into association locus and its stats

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import numpy
import networkx as nx
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule import AbstractMapper, io
from pymodule import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
from pymodule import RBDict
from AssociationLandscape2Peak import AssociationLandscape2Peak

class AssociationPeak2AssociationLocus(AssociationLandscape2Peak):
	__doc__ = __doc__
	option_default_dict = AssociationLandscape2Peak.option_default_dict.copy()
	#get rid of db-related options
	ProcessOptions.removeCertainOptions(option_default_dict=option_default_dict, option_dict_to_remove=AbstractMapper.db_option_dict)
	ProcessOptions.removeCertainOptions(option_default_dict, AssociationLandscape2Peak.my_option_dict)
	
	option_default_dict.update({
							('min_overlap_ratio', 1, float): ['0.05', '', 1, 'minimum overlap ratio, overlap length/total' ],\
							('peakPadding', 0, int): [0, '', 1, 'the padding around each peak (use only to extend the overlap between peaks)' ],\
							})
	
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AssociationLandscape2Peak.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def constructGraph(self, associationPeakRBDictList=[], min_overlap_ratio=0.1):
		"""
		2012.6.24
		"""
		sys.stderr.write("Constructing graph between all result peaks ...")
		rpg = resultPeakGraph = nx.Graph()
		n = len(associationPeakRBDictList)
		counter = 0
		for i in xrange(n):
			for j in xrange(i+1, n):
				peakRBDict1 = associationPeakRBDictList[i]
				peakRBDict2 = associationPeakRBDictList[j]
				compareIns = CNVCompare(min_reciprocal_overlap=min_overlap_ratio)
				for node in peakRBDict1:
					counter += 1
					""""
					segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=node.key.chromosome, \
									span_ls=[node.key.start, node.key.stop], \
									min_reciprocal_overlap=0.0000001, )	#min_reciprocal_overlap doesn't matter here.
						# it's decided by compareIns.
					"""
					targetNodeLs = []
					peakRBDict2.findNodes(node.key, node_ls=targetNodeLs, compareIns=compareIns)
					total_perc_overlapped_by_result2 = 0.
					rpg.add_node(node.key.getKey(), span=node.key.span_ls)
					for result2_node in targetNodeLs:
						result2_peak_id = result2_node.key.getKey()
						if result2_peak_id not in rpg:
							rpg.add_node(result2_peak_id, span = [result2_node.key.chromosome, result2_node.key.start, \
																result2_node.key.stop])
						
						rpg.add_edge(node.key.getKey(), result2_peak_id)
		sys.stderr.write("%s nodes. %s edges. %s connected components.\n"%(rpg.number_of_nodes(), rpg.number_of_edges(), \
															nx.number_connected_components(rpg) ) )
		return rpg
	
	def discoverAssociationLocus(self, resultPeakGraph=None, min_overlap_ratio=0.1):
		"""
		2012.11.20
		2012.6.24
		"""
		rpg = resultPeakGraph
		sys.stderr.write("Discovering association loci from graph of %s nodes. %s edges. %s connected components..."%\
						(rpg.number_of_nodes(), rpg.number_of_edges(), \
						nx.number_connected_components(rpg) ))
		cc_graph_list = nx.connected_component_subgraphs(resultPeakGraph)
		counter = 0
		associationLocusList = []
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
				chromosome = n[0]
				span = resultPeakGraph.node[n]['span']
				start_ls.append(span[0])
				stop_ls.append(span[1])
				result_peak_ls.append(n)
			median_start = numpy.median(start_ls)
			median_stop = numpy.median(stop_ls)
			no_of_peaks = nn
			
			associationLocus = PassingData(chromosome=chromosome, start=median_start, stop=median_stop, \
								no_of_peaks=nn, connectivity=connectivity,\
								result_peak_ls=result_peak_ls)
			associationLocusList.append(associationLocus)
			counter += 1
		sys.stderr.write("%s association loci.\n"%(counter))
		return associationLocusList

	
	def run(self):
		"""
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		
		associationPeakRBDictList =[]
		for inputFname in self.inputFnameLs:
			rbDict = io.constructAssociationPeakRBDictFromHDF5File(inputFname=inputFname, \
											peakPadding=self.peakPadding, groupName='association_peak')
			associationPeakRBDictList.append(rbDict)
		"""
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
			associationPeakRBDictList.append(rbDict)
		"""
		resultPeakGraph = self.constructGraph(associationPeakRBDictList=associationPeakRBDictList, \
									min_overlap_ratio=self.min_overlap_ratio)
		
		associationLocusList = self.discoverAssociationLocus(resultPeakGraph=resultPeakGraph, min_overlap_ratio=self.min_overlap_ratio)
		rbDict = associationPeakRBDictList[-1]	#use the last one to fetch metadata
		io.Association.outputAssociationLociInHDF5(associationLocusList=associationLocusList, filename=self.outputFname, \
										writer=None, groupName='association_locus', closeFile=True,\
					min_MAF=rbDict.min_MAF, neighbor_distance=rbDict.neighbor_distance, \
					max_neighbor_distance=rbDict.max_neighbor_distance, \
					min_score=rbDict.min_score,\
					ground_score=rbDict.ground_score, min_overlap_ratio=self.min_overlap_ratio, peakPadding=self.peakPadding)


if __name__ == '__main__':
	main_class = AssociationPeak2AssociationLocus
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()