#!/usr/bin/env python
"""
Examples:
	%s 
	
	#2012.3.20 phenotype LD's call 32 (3395, peak type 2) vs call 80 (4885, peak type 3).
	%s -l 129 -m /Network/Data/250k/tmp-yh/genomeRBDict_tax3702_padding20kb.pickle -p PASSWORD -s 3395:2,3813:3,4885:3,4992:3 -u yh 
		-w genome_db_passwd -o /tmp/result_3395_2_3813_3_4885_3_4992_3.peak.tsv


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
from pymodule import GenomeDB
from variation.src import Stock_250kDB
import networkx as nx
from pymodule.CNV import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
from pymodule.RBTree import RBDict

class OutputMultiGWASOverlapPeakSpan(AbstractVariationMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVariationMapper.option_default_dict.copy()
	option_default_dict.update({
							('genome_drivername', 1,):['mysql', '', 1, 'which type of database is the genome database? mysql or postgresql', ],\
							('genome_hostname', 1, ): ['banyan', '', 1, 'hostname of the genome db server', ],\
							('genome_dbname', 1, ): ['genome', '', 1, 'genome database name', ],\
							('genome_schema', 0, ): ['', '', 1, 'genome database schema name', ],\
							('genome_db_user', 1, ): ['yh', 'g', 1, 'genome database username', ],\
							('genome_db_passwd', 1, ): [None, 'w', 1, 'genome database password', ],\
							
							('peakPadding', 1, int): [10000, '', 1, 'the extension for each peak on both sides. Rationale is if two peaks are ...'],\
							('result_id_peak_type_id_ls', 1, ): ['', 's', 1, 'result_id:peak_type_id in pairs.Each pair is separated by ",". i.e. 3158:3,3159:1,', ],\
							("list_type_id_list", 1, ): [None, 'l', 1, 'comma/dash separated list of Gene list type. must be in table gene_list_type.'],\
							
							('genomeRBDictPickleFname', 1, ): ['', 'm', 1, 'The file to contain pickled genomeRBDict.'],\
							('genePadding', 0, int): [20000, 'x', 1, "the extension around a gene on both sides to allow association between a locus and a gene. Proper distance is LD-decay."],\
							('tax_id', 0, int): [3702, '', 1, 'Taxonomy ID to get gene position and coordinates.'],\
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		AbstractVariationMapper.__init__(self, inputFnameLs, **keywords)
		if self.list_type_id_list:
			self.list_type_id_list = getListOutOfStr(self.list_type_id_list, data_type=int)
		else:
			self.list_type_id_list = []
		
		self.result_id_peak_type_id_ls = getListOutOfStr(self.result_id_peak_type_id_ls, data_type=str)
		result_id_peak_type_id_ls = []
		for result_id_peak_type_id in self.result_id_peak_type_id_ls:
			result_id_peak_type_id = result_id_peak_type_id.split(':')
			result_id_peak_type_id = map(int, result_id_peak_type_id)
			result_id_peak_type_id_ls.append(result_id_peak_type_id)
		self.result_id_peak_type_id_ls = result_id_peak_type_id_ls
	
	def constructPeakOverlapGraph(self, resultPeakRBDictList=[], genomeRBDict=None, candidate_gene_set=None, outputFname=None):
		"""
		2012.3.16
			make sure each edge is marked with a flag whether it's across two different call methods.
			
			for each component
				1. get the final span (chr, start, stop)
				2. check if any candidate gene is within or touches upon.
				3. check if any edge is across two different call methods,
					which means they are overlapping peaks via two different call methods.
			
		"""
		sys.stderr.write("Constructing result peak overlap graph ...")
		g = nx.Graph()
		
		compareIns = CNVCompare(min_reciprocal_overlap=0.0000001)	#any overlap is an overlap
		no_of_peaks_not_in_result2 = 0
		overlap_ls = []
		counter = 0
		no_of_results = len(resultPeakRBDictList)
		for i in xrange(no_of_results):
			for j in xrange(i+1, no_of_results):
				result1_peakRBDict = resultPeakRBDictList[i]
				result2_peakRBDict = resultPeakRBDictList[j]
				for queryNode in result1_peakRBDict:
					g.add_node(queryNode.value[0].id, chromosome=queryNode.key.chromosome, \
								span_ls=[queryNode.key.start, queryNode.key.stop], \
								call_method_id_ls=[queryNode.value[0].result.call_method_id],\
								phenotype_method_id_ls = [queryNode.value[0].result.phenotype_method_id])	#add this node first, could be singleton
					counter += 1
					segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=queryNode.key.chromosome, \
									span_ls=[queryNode.key.start, queryNode.key.stop], \
									min_reciprocal_overlap=0.0000001, )	#min_reciprocal_overlap doesn't matter here.
						# it's decided by compareIns.
					node_ls = []
					result2_peakRBDict.findNodes(segmentKey, node_ls=node_ls, compareIns=compareIns)
					total_perc_overlapped_by_result2 = 0.
					for node in node_ls:
						overlap1, overlap2, overlap_length, overlap_start_pos, overlap_stop_pos = get_overlap_ratio(segmentKey.span_ls, \
												[node.key.start, node.key.stop])[:5]
						total_perc_overlapped_by_result2 += overlap1
						g.add_edge(queryNode.value[0].id, node.value[0].id, chromosome=queryNode.key.chromosome, \
								span_ls=[min(queryNode.key.start, node.key.start), max(queryNode.key.stop, node.key.stop)], \
								call_method_id_ls=[queryNode.value[0].result.call_method_id, node.value[0].result.call_method_id],\
								phenotype_method_id_ls = [queryNode.value[0].result.phenotype_method_id, \
														node.value[0].result.phenotype_method_id])
					if total_perc_overlapped_by_result2==0:
						no_of_peaks_not_in_result2 += 1
						overlap_ls.append(-0.5)
					else:
						overlap_ls.append(total_perc_overlapped_by_result2)
		sys.stderr.write("%s nodes. %s edges. %s connected components.\n"%(g.number_of_nodes(), g.number_of_edges(), \
															nx.number_connected_components(g)))
		
		sys.stderr.write("Outputting overlap regions ...")
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		header = ['chromosome', 'start', 'stop', 'phenotype_id', 'fileNamePrefix']
		writer.writerow(header)
		no_of_output = 0
		for cc in nx.connected_components(g):
			chromosome= None
			min_start = None
			max_stop = None
			call_method_id_set = set()
			phenotype_method_id_set = set()
			sg = nx.subgraph(g, cc)
			if len(cc)==1:	#only one node, no edges
				node_id = cc[0]
				nodeData = sg.node[node_id]
				min_start, max_stop = nodeData['span_ls']
				chromosome = nodeData['chromosome']
				call_method_id_set = set(nodeData['call_method_id_ls'])
				phenotype_method_id_set = set(nodeData['phenotype_method_id_ls'])
			else:
				for e in sg.edges_iter(data=True):	#data=True, return edge attribute dict in 3-tuple (u,v,data).
					edge_data = e[2]
					chromosome = edge_data['chromosome']
					call_method_id_set = call_method_id_set.union(set(edge_data['call_method_id_ls']))
					phenotype_method_id_set = phenotype_method_id_set.union(set(edge_data['phenotype_method_id_ls']))
					span_ls = edge_data['span_ls']
					if min_start is None:
						min_start = span_ls[0]
					else:
						min_start = min(min_start, span_ls[0])
					if max_stop is None:
						max_stop = span_ls[1]
					else:
						max_stop = max(max_stop, span_ls[1])
			#2012.3.27 don't extend the box before checking for overlap with candidate genes.
			#min_start = max(1, min_start-genomeRBDict.genePadding)	#to extend so that candidate gene could be seen
			#max_stop = max_stop + genomeRBDict.genePadding	#to extend so that candidate gene could be seen
			
			#check whether a candidate gene is within this
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=str(chromosome), \
							span_ls=[min_start, max_stop], \
							min_reciprocal_overlap=0.0000001, )	#min_reciprocal_overlap doesn't matter here.
				# it's decided by compareIns.
			node_ls = []
			genomeRBDict.findNodes(segmentKey, node_ls=node_ls, compareIns=compareIns)
			nearCandidateGene = False
			near_peak_candidate_gene_id_list = []
			for node in node_ls:
				geneSegKey = node.key
				for oneGeneData in node.value:
					if oneGeneData.gene_id in candidate_gene_set:
						nearCandidateGene = True
						near_peak_candidate_gene_id_list.append(oneGeneData.ncbi_gene_id)	#use ncbi gene id instead
						min_start = min(min_start, geneSegKey.span_ls[0])	#2012.3.27 adjust to include the full length of the gene
						max_stop = max(max_stop, geneSegKey.span_ls[1])
			near_peak_candidate_gene_id_list.sort()
			near_peak_candidate_gene_id_list = map(str, near_peak_candidate_gene_id_list)
			
			fileNamePrefixLs = []
			if len(call_method_id_set)>1:
				fileNamePrefixLs.append('olp')
			else:
				call_method_id = call_method_id_set.pop()
				fileNamePrefixLs.append('onlyCall%s'%(call_method_id))
			if nearCandidateGene:
				fileNamePrefixLs.append('cand_%s'%('_'.join(near_peak_candidate_gene_id_list)))
			else:
				fileNamePrefixLs.append("nonCand")
			peak_id_ls_str = map(str, cc)
			fileNamePrefixLs.append("peak_id_%s"%('_'.join(peak_id_ls_str)))
			
			fileNamePrefix = '_'.join(fileNamePrefixLs)
			for phenotype_id in phenotype_method_id_set:
				data_row = [chromosome, min_start, max_stop, phenotype_id, 'pheno_%s_%s'%(phenotype_id, fileNamePrefix)]
				writer.writerow(data_row)
				no_of_output += 1
		del writer
		sys.stderr.write("%s lines outputted.\n"%(no_of_output))
		
	def run(self):
		"""
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_250k = self.db_250k
		
		#need to setup a different db setting
		db_genome = GenomeDB.GenomeDatabase(drivername=self.genome_drivername, username=self.genome_db_user,
						password=self.genome_db_passwd, hostname=self.genome_hostname, database=self.genome_dbname, \
						schema=self.genome_schema)
		db_genome.setup(create_tables=False)
		
		genomeRBDict = db_genome.dealWithGenomeRBDict(self.genomeRBDictPickleFname, tax_id=self.tax_id, \
													max_distance=self.genePadding, debug=self.debug)
		
		candidate_gene_set = set()
		for list_type_id in self.list_type_id_list:
			candidate_gene_list = db_250k.getGeneList(list_type_id)
			candidate_gene_set |= set(candidate_gene_list)
		
		resultPeakRBDictList = []
		for result_id, result_peak_type_id in self.result_id_peak_type_id_ls:
			result_peakRBDict = db_250k.constructRBDictFromResultPeak(result_id, result_peak_type_id, peakPadding=self.peakPadding)
			resultPeakRBDictList.append(result_peakRBDict)
		
		self.constructPeakOverlapGraph(resultPeakRBDictList=resultPeakRBDictList, genomeRBDict=genomeRBDict, \
									candidate_gene_set=candidate_gene_set, outputFname=self.outputFname)
	
if __name__ == '__main__':
	main_class = OutputMultiGWASOverlapPeakSpan
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()