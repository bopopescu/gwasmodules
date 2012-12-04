#!/usr/bin/env python
"""
Examples:
	#2012.9.30
	%s  -i 1-35,39-48,57-82,158-159,161-179,182-186,272-283,314-351,362-380,418-589 -l 80,75 -u yh -z banyan -s 2,3,1,3 -c
	-o ~/script/variation/data/AssociationPeakStat.tsv

	#2012.11.20
	%s -o /tmp/association_locus.h5 -m 0.05 /tmp/5566_association_peak.h5 /tmp/5567_association_peak.h5

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
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, PassingDataList
from pymodule import AbstractMapper, yhio
from pymodule import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
from pymodule import RBDict
from AssociationLandscape2Peak import AssociationLandscape2Peak

class AssociationPeak2AssociationLocus(AssociationLandscape2Peak):
	__doc__ = __doc__
	option_default_dict = AssociationLandscape2Peak.option_default_dict.copy()
	#get rid of db-related options
	ProcessOptions.removeCertainOptions(option_default_dict=option_default_dict, option_dict_to_remove=AbstractMapper.db_option_dict)
	ProcessOptions.removeCertainOptions(option_default_dict, AssociationLandscape2Peak.my_option_dict)
	option_default_dict.pop(('min_MAF', 0, float))
	
	option_default_dict.update({
					('min_overlap_ratio', 1, float): ['0.05', '', 1, 'minimum overlap ratio, overlap length/total' ],\
					('peakPadding', 0, int): [0, '', 1, 'the padding around each peak (use only to extend the overlap between peaks)' ],\
					})
	
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AssociationLandscape2Peak.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		#2012.11.22 a sell-auto-increment value to keep track of all peaks this program is handling. 
		self.peakID = 0
	
	def getAssociationPeakKeyFromRBTreeKey(self, RBDict=None, RBNodeKey=None):
		"""
		2012.11.22
			each RBNodeKey is constructed via yhio.association.constructAssociationPeakRBDictFromHDF5File
			
			result_id has to be baked-in, otherwise, same location from different results would all be collapsed as one peak.
		"""
		return (RBDict.result_id, RBNodeKey.chromosome, RBNodeKey.start, RBNodeKey.stop)
	
	def constructGraph(self, associationPeakRBDictList=[], min_overlap_ratio=0.1):
		"""
		2012.6.24
		"""
		sys.stderr.write("Constructing graph between all result peaks ...")
		rpg = resultPeakGraph = nx.Graph()
		n = len(associationPeakRBDictList)
		counter = 0
		compareIns = CNVCompare(min_reciprocal_overlap=min_overlap_ratio)
		for i in xrange(n):
			for j in xrange(i+1, n):
				peakRBDict1 = associationPeakRBDictList[i]
				peakRBDict2 = associationPeakRBDictList[j]
				for result1Node in peakRBDict1:
					counter += 1
					""""
					segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=result1Node.key.chromosome, \
									span_ls=[result1Node.key.start, result1Node.key.stop], \
									min_reciprocal_overlap=0.0000001, )	#min_reciprocal_overlap doesn't matter here.
						# it's decided by compareIns.
					"""
					targetNodeLs = []
					peakRBDict2.findNodes(result1Node.key, node_ls=targetNodeLs, compareIns=compareIns)
					total_perc_overlapped_by_result2 = 0.
					result1NodeHashKey = self.getAssociationPeakKeyFromRBTreeKey(peakRBDict1, result1Node.key)
					rpg.add_node(result1NodeHashKey, chromosome=result1Node.key.chromosome, \
								span=result1Node.key.span_ls, \
								result_id=peakRBDict1.result_id)
					for result2Node in targetNodeLs:
						result2NodeHashKey = self.getAssociationPeakKeyFromRBTreeKey(peakRBDict2, result2Node.key)
						if result2NodeHashKey not in rpg:
							rpg.add_node(result2NodeHashKey, chromosome=result2Node.key.chromosome, \
										span = [result2Node.key.start, result2Node.key.stop],\
										result_id=peakRBDict2.result_id)
						
						rpg.add_edge(result1NodeHashKey, result2NodeHashKey)
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
			result_id_set = set()
			chromosome_set = set()	#should be only one chromosome
			for n in cc_graph:
				nodeObject = resultPeakGraph.node[n]
				chromosome_set.add(nodeObject['chromosome'])
				span = nodeObject['span']
				start_ls.append(span[0])
				stop_ls.append(span[1])
				result_peak_ls.append(n)
				result_id_set.add(nodeObject['result_id'])
			if len(chromosome_set)>1:
				sys.stderr.write("Error: %s chromosomes (%s) in one connected component.\n"%(len(chromosome_set), repr(chromosome_set)))
				sys.exit(7)
			median_start = numpy.median(start_ls)
			median_stop = numpy.median(stop_ls)
			no_of_results = len(result_id_set)
			associationLocus = PassingDataList()
			#assign each value separately to impose the order of variables in associationLocus's internal list
			associationLocus.chromosome = chromosome_set.pop()
			associationLocus.start=median_start
			associationLocus.stop=median_stop
			associationLocus.no_of_peaks=nn
			associationLocus.connectivity=connectivity
			associationLocus.no_of_results=no_of_results
			associationLocus.result_peak_ls=result_peak_ls
			#PassingDataList is sortable via (chromosome, start, stop ...)
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
		call_method_id_set = set()
		cnv_method_id_set = set()
		phenotype_method_id_set = set()
		analysis_method_id_set = set()
		for inputFname in self.inputFnameLs:
			rbDict = yhio.Association.constructAssociationPeakRBDictFromHDF5File(inputFname=inputFname, \
											peakPadding=self.peakPadding, groupName='association_peak')
			associationPeakRBDictList.append(rbDict)
			utils.addObjectAttributeToSet(objectVariable=rbDict, attributeName='call_method_id', \
										setVariable=call_method_id_set)
			utils.addObjectAttributeToSet(objectVariable=rbDict, attributeName='cnv_method_id', \
										setVariable=cnv_method_id_set)
			utils.addObjectAttributeToSet(objectVariable=rbDict, attributeName='phenotype_method_id', \
										setVariable=phenotype_method_id_set)
			utils.addObjectAttributeToSet(objectVariable=rbDict, attributeName='analysis_method_id', \
										setVariable=analysis_method_id_set)
		
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
		attributeDict = {'min_MAF':rbDict.min_MAF, 'neighbor_distance':rbDict.neighbor_distance, \
				'max_neighbor_distance':rbDict.max_neighbor_distance, 'min_score':rbDict.min_score,\
				'ground_score':getattr(rbDict, 'ground_score', None), \
				'call_method_id_ls': numpy.array(list(call_method_id_set)),\
				'cnv_method_id_ls': numpy.array(list(cnv_method_id_set)),\
				'phenotype_method_id_ls': numpy.array(list(phenotype_method_id_set)),\
				'analysis_method_id_ls': numpy.array(list(analysis_method_id_set)),\
				'min_overlap_ratio':self.min_overlap_ratio, \
				'peakPadding':self.peakPadding, 'total_no_of_results':len(associationPeakRBDictList)}
		yhio.Association.outputAssociationLociInHDF5(associationLocusList=associationLocusList, filename=self.outputFname, \
							writer=None, groupName='association_locus', closeFile=True,\
							attributeDict=attributeDict)


if __name__ == '__main__':
	main_class = AssociationPeak2AssociationLocus
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()