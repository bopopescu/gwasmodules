#!/usr/bin/env python
"""

Examples:
	# 2010-7-29 merge deletions from CNVMethod 16 to 18 with min_overlap_ratio set at 0.9.
	CNVMergeAcrossArrays.py -u yh -a 16 -m 18 -z banyan.usc.edu -n 0.9 -r -c
	
	# 2010-8-13 use betweenness centrality to cut graph
	~/script/variation/src/CNVMergeAcrossArrays.py -u  yh -a 16 -m 25 -n 0.8 -z banyan.usc.edu -i 0.8 -c -g3
	
	# 2010-9-13 merging potential-bi-allelic deletions only
	CNVMergeAcrossArrays.py -u  yh -a 16 -m 28 -z banyan.usc.edu -i 0.95 -n 0.03 -c -g4
	
Description:
	If two segments from two arrays share much of the overlap (>=min_overlap_ratio), they are regarded as one segment
		and connected in a graph.
	All connected components in that graph are collapsed into one segment (median start, stop, probability/score)
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import getopt, numpy
from variation.src.common import nt2number, number2nt
from pymodule import ProcessOptions, PassingData
import Stock_250kDB
from pymodule.CNV import get_overlap_ratio
import networkx as nx


class CNVMergeAcrossArrays(object):
	__doc__ = __doc__
	"""
	2009-2-12
	"""
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
						('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('raw_cnv_method_id', 1, int): [8, 'a', 1, 'CNV method id of the training data in CNVCall'],\
						('cnv_method_id', 1, int): [10, 'm', 1, 'CNV method id for the predicted deletions to be stored in CNVCall.'],\
						('cnv_type_id', 1, int): [1, 'y', 1, 'CNV type id. table CNVType', ],\
						('min_overlap_ratio', 1, float): [0.8, 'n', 1, 'overlap ratios (overlap/segment-length) for two segments have to be above this.'],\
						('minConnectivity', 1, float): [0.8, 'i', 1, 'minimum connectivity for the deletion-graph, only for betweenness-centrality method'],\
						('graphClusterMethod', 1, int): [1, 'g', 1, 'use which method to cluster graph. 1: connected component 2: max-clique, \
							3: betweenness-centrality, 4: bi-allelic connected components (argument min_overlap_ratio is set low, 0.05, to create edge. implicit 0.85 is used.)'],\
						('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2009-2-12
		"""
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
												class_to_have_attr=self)
		self.graphClusterMethod2func = {
									1: self.clearOutGraphAndSaveByConnectedComponents,\
									2: self.clearOutGraphAndSaveByClique,\
									3: self.clearOutGraphAndSaveByBetweennessCentrality,\
									4: self.clearOutGraphByRetainingBiAllelicConnectedComponents}
		self.clearOutGraphAndSave = self.graphClusterMethod2func[self.graphClusterMethod]
		
	cnv_unique_key2cnv = {}	#2010-8-8 to avoid duplicated CNVs stored into db
	cnv_array_call_unique_key2cnv = {}	#2010-8-10
	def saveSegmentsIntoCNVAndCNVArrayCall(self, segments, param_obj=None):
		"""
		
		"""
		if not hasattr(param_obj, 'no_of_total'):
			setattr(param_obj, 'no_of_total', 0)
		if not hasattr(param_obj, 'no_of_into_db'):
			setattr(param_obj, 'no_of_into_db', 0)
		if not hasattr(param_obj, 'duplicated_cnv_unique_key2count'):
			setattr(param_obj, 'duplicated_cnv_unique_key2count', {})
		if not hasattr(param_obj, 'already_in_db_cnv_unique_key2count'):
			setattr(param_obj, 'already_in_db_cnv_unique_key2count', {})
		if not hasattr(param_obj, 'duplicated_cnv_array_call_unique_key2count'):
			setattr(param_obj, 'duplicated_cnv_array_call_unique_key2count', {})
		if not hasattr(param_obj, 'already_in_db_cnv_array_call_unique_key2count'):
			param_obj.already_in_db_cnv_array_call_unique_key2count = {}
		
		cnv_type_id = getattr(param_obj, 'cnv_type_id', 1)
		cnv_method_id = getattr(param_obj, 'cnv_method_id', 1)
		
		session = getattr(param_obj, 'session', None)
		if session is None:
			sys.stderr.write("Error: db session is not available.\n")
			return
		param_obj.no_of_total += 1
		
		import numpy
		segment_matrix = numpy.array(segments, numpy.float32)
		start = int(numpy.median(segment_matrix[:,1]))
		stop = int(numpy.median(segment_matrix[:,2]))
		chromosome = int(segment_matrix[0][0])
		
		cnv_unique_key = (chromosome, start, stop, cnv_method_id, cnv_type_id)
		newCNVObject = False
		if cnv_unique_key in self.cnv_unique_key2cnv:
			cnv = self.cnv_unique_key2cnv.get(cnv_unique_key)
			if cnv_unique_key not in param_obj.duplicated_cnv_unique_key2count:
				param_obj.duplicated_cnv_unique_key2count[cnv_unique_key] = 1	#started from 1, rather than 0
			param_obj.duplicated_cnv_unique_key2count[cnv_unique_key] += 1
		else:
			rows = Stock_250kDB.CNV.query.filter_by(chromosome=chromosome).filter_by(start=start).\
				filter_by(stop = stop).filter_by(cnv_type_id = cnv_type_id).\
				filter_by(cnv_method_id=cnv_method_id)
			if rows.count()==0:	#make sure it's not in db yet.
				newCNVObject = True
			else:
				cnv = rows.first()
				if cnv_unique_key not in param_obj.already_in_db_cnv_unique_key2count:
					param_obj.already_in_db_cnv_unique_key2count[cnv_unique_key] = 1	#started from 1, rather than 0
				param_obj.already_in_db_cnv_unique_key2count[cnv_unique_key] += 1
				#return	#2010-8-11 temporary. if CNV is in db, all CNVArrayCall associated with it shall be in db as well unless bug.
		
		frequency = float(len(param_obj.array_id2data))/param_obj.no_of_total_arrays
		if newCNVObject:
			size_affected = stop - start + 1
			score = numpy.median(segment_matrix[:,3])
			score_std = numpy.std(segment_matrix[:,3])
		
			cnv = Stock_250kDB.CNV(chromosome = chromosome, start=start, stop=stop, size_affected=size_affected,\
								score=score, score_std=score_std,\
								frequency=frequency, cnv_method_id=cnv_method_id, cnv_type_id=cnv_type_id,)
			self.cnv_unique_key2cnv[cnv_unique_key] = cnv	#2010-8-9
		else:
			cnv.frequency += frequency	#2010-8-6 add the frequency. a new CNV is found to be identical to an old one.
		session.add(cnv)
		
		for cnv_call_id in segment_matrix[:,-1]:	#record the sources for each CNV
			cnv_call = Stock_250kDB.CNVCall.get(int(cnv_call_id))	#2010-8-3 bugfix. int() is a must even though difference is .0.
			cnv.cnv_call_ls.append(cnv_call)
		
		# 2010-8-9 flush after the for loop is over to avoid the situation that cnv object gets saved half way
		needFlush = False
		for array_id, array_data in param_obj.array_id2data.iteritems():
			probability = numpy.median(array_data.probability_ls)
			saveNewObj = False
			cnv_array_call_unique_key = (array_id, cnv_unique_key, cnv_method_id)
			if cnv_array_call_unique_key in self.cnv_array_call_unique_key2cnv:
				cnv_array_call = self.cnv_array_call_unique_key2cnv.get(cnv_array_call_unique_key)
				if cnv_array_call_unique_key not in param_obj.duplicated_cnv_array_call_unique_key2count:
					param_obj.duplicated_cnv_array_call_unique_key2count[cnv_array_call_unique_key] = 1	#started from 1, rather than 0
				param_obj.duplicated_cnv_array_call_unique_key2count[cnv_array_call_unique_key] += 1
			elif cnv.id:
				rows = Stock_250kDB.CNVArrayCall.query.filter_by(array_id=array_id).filter_by(cnv_id=cnv.id).\
					filter_by(cnv_method_id=cnv_method_id)
				if rows.count()==0:	
					saveNewObj = True
				else:
					if cnv_array_call_unique_key not in param_obj.already_in_db_cnv_array_call_unique_key2count:
						param_obj.already_in_db_cnv_array_call_unique_key2count[cnv_array_call_unique_key] = 1	#started from 1, rather than 0
					param_obj.already_in_db_cnv_array_call_unique_key2count[cnv_array_call_unique_key] += 1
				
			else:
				saveNewObj = True
			if saveNewObj:
				cnv_array_call = Stock_250kDB.CNVArrayCall(array_id=array_id, cnv_method_id=cnv_method_id,\
														score = probability)
				cnv_array_call.cnv = cnv
				self.cnv_array_call_unique_key2cnv[cnv_array_call_unique_key] = cnv_array_call	#2010-8-10
				session.add(cnv_array_call)
				param_obj.no_of_into_db += 1
				if param_obj.no_of_into_db%5000==0:
					needFlush = True
		
		if needFlush:
			try:
				session.flush()
				#session.expunge_all()
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
				import pdb
				pdb.set_trace()
		
		report = getattr(param_obj, 'report', False)
		if report and param_obj.no_of_total%1000==0:
			no_of_duplicated_cnv_unique_keys = len(param_obj.duplicated_cnv_unique_key2count)
			no_of_already_in_db_cnv_unique_keys = len(param_obj.already_in_db_cnv_unique_key2count)
			k = len(param_obj.duplicated_cnv_array_call_unique_key2count)
			l = len(param_obj.already_in_db_cnv_array_call_unique_key2count)
			sys.stderr.write('%s%s into db out of %s. %s duplicated cnv_unique_key(s). %s already_in_db cnv_unique_key(s).\n\
						%s duplicated cnv_array_call_unique_key(s). %s already_in_db_cnv_array_call_unique_key(s).\n'%\
						('\x08'*100, param_obj.no_of_into_db, param_obj.no_of_total, no_of_duplicated_cnv_unique_keys,\
						no_of_already_in_db_cnv_unique_keys, k, l))
	
	def cutGraphByBetweenCentrality(self, G, cluster_ls=[], minConnectivity=0.8, centralityMeasure=1):
		"""
		2010-8-10
			cut the edge with biggest centrality until connectivity of all connected components fall are above minConnectivity
		"""
		H = nx.connected_component_subgraphs(G)
		for subG in H:
			n = subG.number_of_nodes()
			if n<=2:
				cluster_ls.append(subG.nodes())
				continue
			else:
				e = subG.number_of_edges()
				connectivity = 2*e/float(n*(n-1))
				if connectivity>minConnectivity:
					cluster_ls.append(subG.nodes())
					continue
				else:
					from annot.bin.graph.cc_from_edge_list import ClusterByEBC
					#2010-8-15 for debug purpose, to see which graph caused segmentation fault
					sys.stderr.write("Working on graph ...\n")
					for e in subG.edges():
						sys.stderr.write("\t edge (%s, %s)\n"%(e[0], e[1]))
					cf_instance =ClusterByEBC(subG.edges(), 1, minConnectivity)
					cf_instance.run()
					cluster_ls.extend(cf_instance.cc_vertex_list)
					"""
					if centralityMeasure==1:
						edge2centrality = nx.edge_betweenness_centrality(subG, normalized=True, weighted_edges=True)
					else:
						edge2centrality = nx.edge_current_flow_betweenness_centrality(subG, normalized=True)
						# weight is taken into account if it's set. 
					centrality_edge_ls = []
					for edge, centrality in edge2centrality.iteritems():
						centrality_edge_ls.append((centrality, edge))
					centrality_edge_ls.sort()
					edgeWithTopCentrality = centrality_edge_ls[-1][1]
					subG.remove_edge(edgeWithTopCentrality[0], edgeWithTopCentrality[1])
					self.cutGraphByBetweenCentrality(subG, cluster_ls, minConnectivity=minConnectivity, \
													centralityMeasure=centralityMeasure)
					"""
		return
	
	def clearOutGraphAndSaveByConnectedComponents(self, G, segment_ls, param_obj=None):
		"""
		2010-7-31
		"""
		H = nx.connected_component_subgraphs(G)
		for subG in H:
			segments = []
			array_id2data = {}	#different CNVCall s from the same array could be included in one component.
			for i in subG.nodes():
				segment = segment_ls[i]
				probability = segment[-3]
				array_id = segment[-2]
				if array_id not in array_id2data:
					array_id2data[array_id] = PassingData(probability_ls = [])
				array_id2data[array_id].probability_ls.append(probability)
				segments.append(segment_ls[i])
			param_obj.array_id2data = array_id2data
			self.saveSegmentsIntoCNVAndCNVArrayCall(segments, param_obj=param_obj)
		#2010-8-11 to clear it up
		self.cnv_unique_key2cnv = {}
		self.cnv_array_call_unique_key2cnv = {}
	
	def clearOutGraphAndSaveByBetweennessCentrality(self, G, segment_ls, param_obj=None,):
		"""
		2010-8-10
		
		"""
		cluster_ls = []
		self.cutGraphByBetweenCentrality(G, cluster_ls=cluster_ls, minConnectivity=getattr(param_obj, 'minConnectivity', 0.8),\
										centralityMeasure=1)
		saved_node_set = set()
		for node_list in cluster_ls:
			segments = []
			array_id2data = {}	#different CNVCall s from the same array could be included in one component.
			for i in node_list:
				if i in saved_node_set:	#already saved in a bigger clique.
					continue
				saved_node_set.add(i)
				segment = segment_ls[i]
				probability = segment[-3]
				array_id = segment[-2]
				if array_id not in array_id2data:
					array_id2data[array_id] = PassingData(probability_ls = [])
				array_id2data[array_id].probability_ls.append(probability)
				segments.append(segment_ls[i])
			if len(segments)>0:	#could be empty because nodes all went to the bigger cliques.
				param_obj.array_id2data = array_id2data
				self.saveSegmentsIntoCNVAndCNVArrayCall(segments, param_obj=param_obj)
		#2010-8-11 to clear it up
		self.cnv_unique_key2cnv = {}
		self.cnv_array_call_unique_key2cnv = {}
		
	def clearOutGraphAndSaveByClique(self, G, segment_ls, param_obj=None):
		"""
		2010-8-6
			similar to clearOutGraphAndSaveByConnectedComponents, but use clique-finding algorithm
		"""
		saved_node_set = set()
		for node_list in nx.find_cliques(G):
			segments = []
			array_id2data = {}	#different CNVCall s from the same array could be included in one component.
			for i in node_list:
				if i in saved_node_set:	#already saved in a bigger clique.
					continue
				saved_node_set.add(i)
				segment = segment_ls[i]
				probability = segment[-3]
				array_id = segment[-2]
				if array_id not in array_id2data:
					array_id2data[array_id] = PassingData(probability_ls = [])
				array_id2data[array_id].probability_ls.append(probability)
				segments.append(segment_ls[i])
			if len(segments)>0:	#could be empty because nodes all went to the bigger cliques.
				param_obj.array_id2data = array_id2data
				self.saveSegmentsIntoCNVAndCNVArrayCall(segments, param_obj=param_obj)
		#2010-8-11 to clear it up
		self.cnv_unique_key2cnv = {}
		self.cnv_array_call_unique_key2cnv = {}
	
	def clearOutGraphByRetainingBiAllelicConnectedComponents(self, G, segment_ls, param_obj=None):
		"""
		2010-9-13
			Two conditions for a CC to be regarded as representing just one allele.
				1. connectivity is >param_obj.minConnectivity
				2. edge weight >0.85 for each edge 
		"""
		H = nx.connected_component_subgraphs(G)
		for subG in H:
			minEdgeWeight = None
			for e in subG.edges_iter():
				weight = subG.get_edge_data(e[0], e[1])['weight']
				if minEdgeWeight is None or weight<minEdgeWeight:
					minEdgeWeight = weight
				if minEdgeWeight<0.85:
					break
			if minEdgeWeight<0.85:
				continue
			
			n = subG.number_of_nodes()
			if n<=2:
				connectivity = 1
			else:
				e = subG.number_of_edges()
				connectivity = 2*e/float(n*(n-1))
			
			if connectivity < param_obj.minConnectivity:
				continue
			
			segments = []
			array_id2data = {}	#different CNVCall s from the same array could be included in one component.
			
			
			for i in subG.nodes():
				segment = segment_ls[i]
				probability = segment[-3]
				array_id = segment[-2]
				if array_id not in array_id2data:
					array_id2data[array_id] = PassingData(probability_ls = [])
				array_id2data[array_id].probability_ls.append(probability)
				segments.append(segment_ls[i])
			param_obj.array_id2data = array_id2data
			self.saveSegmentsIntoCNVAndCNVArrayCall(segments, param_obj=param_obj)
		#2010-8-11 to clear it up
		self.cnv_unique_key2cnv = {}
		self.cnv_array_call_unique_key2cnv = {}
	
	@classmethod
	def getNonDuplicateArraysWithHighestMedianIntensity(cls, db_250k, raw_cnv_method_id=None, table_name=None):
		"""
		2010-9-30
			add argument table_name
		2010-7-29
		"""
		sys.stderr.write("Getting non-duplicate arrays with highest median intensity ...")
		if table_name is None:
			table_name = Stock_250kDB.CNVCall.table.name
		rows = db_250k.metadata.bind.execute("select a.* from view_array a, (select distinct array_id from %s where cnv_method_id=%s) \
				t where t.array_id=a.array_id"%(table_name, raw_cnv_method_id))
		ecotype_id2median_intensity_array_id_ls = {}
		no_of_arrays = 0
		for row in rows:
			no_of_arrays += 1
			ecotype_id = row.maternal_ecotype_id
			if ecotype_id not in ecotype_id2median_intensity_array_id_ls:
				ecotype_id2median_intensity_array_id_ls[ecotype_id] = []
			ecotype_id2median_intensity_array_id_ls[ecotype_id].append((row.median_intensity, row.array_id))
		
		non_duplicate_array_id_ls = []
		for ecotype_id, median_intensity_array_id_ls in ecotype_id2median_intensity_array_id_ls.iteritems():
			median_intensity_array_id_ls.sort()
			non_duplicate_array_id_ls.append(median_intensity_array_id_ls[-1][1])
		sys.stderr.write(" reduced from %s arrays to %s.\n"%(no_of_arrays, len(non_duplicate_array_id_ls)))
		return non_duplicate_array_id_ls
	
	def mergeSegmentsAcrossArrays(self, db_250k, raw_cnv_method_id = None, min_overlap_ratio=0.3, \
							maxDeletionLength=50000, param_obj=None, clearOutGraphAndSave=None):
		"""
		2010-7-29
		"""
		sys.stderr.write("Merging segments across arrays ... \n")
		#query = Stock_250kDB.CNVCall.query.filter_by(cnv_method_id=raw_cnv_method_id).\
		#	filter_by(cnv_type_id=param_obj.cnv_type_id).\
		#	order_by(Stock_250kDB.CNVCall.chromosome).order_by(Stock_250kDB.CNVCall.start).order_by(Stock_250kDB.CNVCall.stop)
		
		base_sql_string = 'select * from %s where cnv_method_id=%s and cnv_type_id=%s order by chromosome, start, stop'%\
				(Stock_250kDB.CNVCall.table.name, raw_cnv_method_id, param_obj.cnv_type_id)
		
		
		segment_index = 0
		segment_ls = []
		leftOverlapSegmentIndexLs = []	#this is the index of the segments on the left that are still overlapping with the current one,
		# which would be the segments potentially overlapping with the next segment. 
		
		G=nx.Graph()
		
		counter = 0
		real_counter = 0
		block_size = 10000
		rows = db_250k.metadata.bind.execute(base_sql_string)
		
		#rows = query.offset(counter).limit(block_size)
		if self.debug:
			minNumberOfNodesToClear = 1000
		else:	#set it higher if not in debug
			minNumberOfNodesToClear = 1000
		#while rows.count()!=0:
		for row in rows:
			counter += 1
			segment = [row.chromosome, row.start, row.stop, row.probability, row.array_id, row.id]
			current_segment_index = len(segment_ls)
			segment_ls.append(segment)
			G.add_node(current_segment_index)
			
			#calcualte the overlap ratio between the current segment and all previous potentially overlapping segments
			newLeftOverlapSegmentIndexLs = []
			for aheadSegmentIndex in leftOverlapSegmentIndexLs:
				aheadSegment = segment_ls[aheadSegmentIndex]
				if aheadSegment[0]!=segment[0]:		#on different chromosomes, ignore
					#leftMostOverlapSegmentIndex = aheadSegmentIndex + 1
					continue
				overlap1, overlap2, overlap_length = get_overlap_ratio([aheadSegment[1], aheadSegment[2]], \
															[segment[1], segment[2]])[:3]
				#if overlap1==0 or overlap2==0:	#2010-8-12 this is wrong. this could cut out earlier segments
				#	that could still be overlapping with future segments.
				#	leftMostOverlapSegmentIndex = aheadSegmentIndex + 1
				
				if overlap1>=min_overlap_ratio or overlap2>=min_overlap_ratio:
					# 2010-8-12 doesn't require both to be >=min_overlap_ratio because
					# scenario 1. overlap1>=min_overlap_ratio but overlap2 not
					# 	the future segments might be shorter than current segment while with similar overlap length,
					# 	thus its overlap2 could be
					# 	>=min_overlap_ratio even though the current one's is <min_overlap_ratio
					# scenario 2. opposite of scenario 1.
					# 	the current segment is embedded in the previous segment,
					#	future segments might extend over this aheadSegment and overlap1 reaches above min_overlap_ratio
					newLeftOverlapSegmentIndexLs.append(aheadSegmentIndex)
				
				if overlap1<0 or overlap2<0:
					sys.stderr.write("aheadSegment %s.\n"%repr(aheadSegment))
					sys.stderr.write("segment %s.\n"%repr(segment))
					sys.stderr.write("Error: overlap1 %s or overlap2 %s is zero.\n"%(overlap1, overlap2)) 
				elif overlap1>=min_overlap_ratio and overlap2>=min_overlap_ratio:
					G.add_edge(aheadSegmentIndex, current_segment_index, weight=min(overlap1, overlap2))
					real_counter += 1
			
			# always include the current_segment_index for the next segment to compare
			newLeftOverlapSegmentIndexLs.append(current_segment_index)
			leftOverlapSegmentIndexLs = newLeftOverlapSegmentIndexLs
			
			if counter%1000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
			if len(leftOverlapSegmentIndexLs)==1 and G.number_of_nodes()>=minNumberOfNodesToClear:
				#2010-8-12 current_segment_index is not connected to any prior node.
				G.remove_node(current_segment_index)	#2010-8-12 bugfix. remove current_segment_index before graph-clustering
				#2010-8-12  current_segment_index belongs to the next graph.
				clearOutGraphAndSave(G, segment_ls, param_obj = param_obj)
				
				# release memory
				del G
				#for k in xrange(leftMostOverlapSegmentIndex):
				#	segment_ls[k] = None
				
				# construct a new graph
				G=nx.Graph()
				G.add_node(current_segment_index)	#don't forget to add current node
			#2010-7-31 for debug purpose, break out after one bactch
			if self.debug and counter>=50000:
				break
			#rows = db_250k.metadata.bind.execute('%s offset %s limit %s'%(base_sql_string, counter, block_size))
			#rows = query.offset(counter).limit(block_size)
		sys.stderr.write("%s%s\t%s\n"%('\x08'*80, counter, real_counter))
		# clear out remaining nodes
		clearOutGraphAndSave(G, segment_ls, param_obj = param_obj)
		
		no_of_duplicated_cnv_unique_keys = len(param_obj.duplicated_cnv_unique_key2count)
		no_of_already_in_db_cnv_unique_keys = len(param_obj.already_in_db_cnv_unique_key2count)
		k = len(param_obj.duplicated_cnv_array_call_unique_key2count)
		l = len(param_obj.already_in_db_cnv_array_call_unique_key2count)
		sys.stderr.write('%s%s into db out of %s. %s duplicated cnv_unique_key(s). %s already_in_db cnv_unique_key(s).\n\
					%s duplicated cnv_array_call_unique_key(s). %s already_in_db_cnv_array_call_unique_key(s).\n'%\
					('\x08'*100, param_obj.no_of_into_db, param_obj.no_of_total, no_of_duplicated_cnv_unique_keys,\
					no_of_already_in_db_cnv_unique_keys, k, l))
			
		sys.stderr.write("Items in duplicated_cnv_unique_key2count:\n")
		for cnv_unique_key, count in param_obj.duplicated_cnv_unique_key2count.iteritems():
			sys.stderr.write("\t CNV unique key %s, count=%s.\n"%(repr(cnv_unique_key), count))
		sys.stderr.write("Items in already_in_db_cnv_unique_key2count:\n")
		for cnv_unique_key, count in param_obj.already_in_db_cnv_unique_key2count.iteritems():
			sys.stderr.write("\t CNV unique key %s, count=%s.\n"%(repr(cnv_unique_key), count))
		
		sys.stderr.write("Items in duplicated_cnv_array_call_unique_key2count:\n")
		for key, count in param_obj.duplicated_cnv_array_call_unique_key2count.iteritems():
			sys.stderr.write("\t CNVArrayCall unique key %s, count=%s.\n"%(repr(key), count))
		sys.stderr.write("Items in already_in_db_cnv_array_call_unique_key2count:\n")
		for key, count in param_obj.already_in_db_cnv_array_call_unique_key2count.iteritems():
			sys.stderr.write("\t CNVArrayCall unique key %s, count=%s.\n"%(repr(key), count))
		
		#self.getAndDrawConnectedComponentsStat(G, raw_cnv_method_id, param_obj)
	
	def getAndDrawConnectedComponentsStat(self, G, raw_cnv_method_id, param_obj=None):
		"""
		2010-7-31
			to check what summary stats look like
		"""
		sys.stderr.write("Getting the connected components info ...")
		H = nx.connected_component_subgraphs(G)
		density_ls = []
		no_of_nodes_ls = []
		no_of_edges_ls = []
		num_edges_num_nodes_ratio_ls = []
		for subG in H:
			n = subG.number_of_nodes()
			no_of_edges = subG.number_of_edges()
			no_of_nodes_ls.append(math.log10(n))
			if no_of_edges==0:
				no_of_edges_ls.append(-1)
			else:
				no_of_edges_ls.append(math.log10(no_of_edges))
			if n>1:
				density = float(no_of_edges)/(n*(n-1)/2.0)
				num_edges_num_nodes_ratio = math.log10(float(no_of_edges)/n)
			else:
				density = 0
				num_edges_num_nodes_ratio = -1
			density_ls.append(density)
			num_edges_num_nodes_ratio_ls.append(num_edges_num_nodes_ratio)
		sys.stderr.write(".\n")
		
		sys.stderr.write("drawing ...")
		x_ls = density_ls
		dataType = 3
		fileNamePrefix = 'OverlapDeletionGraphDensityHist'
		xlim_in_1D = None
		output_dir = os.path.expanduser('~/script/variation/data/CNV/AcrossArrayOverlapDeletionGraph/')
		xlabel = 'density'
		
		if dataType==1:
			y_ls = no_of_nodes_ls
			ylabel_in_2D = 'log10(no_of_nodes)'
			hist_x_ls = density_ls
			hist_xlabel = 'density'
		elif dataType==2:
			y_ls = no_of_edges_ls
			ylabel_in_2D = 'log10(no_of_edges)'
			hist_x_ls = no_of_edges_ls
			hist_xlabel = ylabel_in_2D
		elif dataType==3:
			y_ls = num_edges_num_nodes_ratio_ls
			ylabel_in_2D = 'log10(no_of_edges/no_of_nodes)'
			hist_x_ls = num_edges_num_nodes_ratio_ls
			hist_xlabel = ylabel_in_2D
		else:
			y_ls = no_of_nodes_ls
			ylabel_in_2D = 'log10(no_of_nodes)'
			hist_x_ls = no_of_nodes_ls
			hist_xlabel = ylabel_in_2D
		import pylab
		from pymodule.utils import addExtraLsToFilenamePrefix
		pylab.clf()
		output_fname = os.path.join(output_dir, '%s_CNVMethod%s_CNVType%s_dataType%s.png'%(fileNamePrefix, \
									raw_cnv_method_id,	param_obj.cnv_type_id, dataType, ))
		title = '%s components, %s nodes, %s edges.'%(len(x_ls), G.number_of_nodes(), G.number_of_edges())
		
		pylab.title(title)
		"""
		import statistics	# 2010-5-30 package from Michiel De Hoon
		y, x = statistics.pdf(x_ls)
		pylab.loglog(x, y, alpha=0.7)
		pylab.grid(True, alpha=0.6)
		"""
		pylab.hist(hist_x_ls, 20, log=True)
		if xlim_in_1D:
			pylab.xlim(xlim_in_1D)
		
		pylab.xlabel(hist_xlabel)
		pylab.ylabel('Count')
		pylab.savefig(output_fname, dpi=300)
		
		
		C_ls = [1]*len(y_ls)
		fig_fname = addExtraLsToFilenamePrefix(output_fname, ['2D'])
		from misc import CNV
		CNV.drawHexbin(x_ls, y_ls, C_ls, fig_fname=fig_fname, gridsize=20, \
				title=title, \
				xlabel=xlabel, \
				ylabel=ylabel_in_2D, \
				colorBarLabel='log10(count)', reduce_C_function=CNV.logSum)
		
		sys.stderr.write("Done.\n")
		#merged_segment_ls = self.mergeOverlappingORCloseSegmentsByGraph(segment_ls, max_reciprocal_gap_ratio=max_gap_ratio, \
		#															max_gap_len=max_gap_len, mergeFunc=self.mergeTwoSegments,\
		#															maxDeletionLength=maxDeletionLength)
		
		"""
		from CNVPredictDeletionBySVM import CNVPredictDeletionBySVM
		for merged_segment in merged_segment_ls:
			chromosome, start, stop, start_probe_id, stop_probe_id, no_of_probes_covered, size_affected,\
					amplitude, probability = merged_segment[:9]
			cnv_segment_obj = PassingData(array_id=array_id, start_probe_id=start_probe_id, stop_probe_id=stop_probe_id, \
								no_of_probes=no_of_probes_covered, amplitude=amplitude, segment_length=size_affected, \
								segment_chromosome=chromosome, \
								segment_start_pos=start, segment_stop_pos=stop, \
								median_intensity=None, probability=probability)
			#CNVPredictDeletionBySVM.saveSegmentObj(param_obj, cnv_segment_obj)
		"""
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		2010-7-29
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_250k = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				  			password=self.db_passwd, hostname=self.hostname, database=self.dbname, 
				   			schema=self.schema)
		db_250k.setup(create_tables=False)
		session = db_250k.session
		session.expire_on_commit = False	#otherwise a flush() would cause:	
		"""
		sqlalchemy.orm.exc.DetachedInstanceError: Instance <...x770eb6d0> is not bound to a Session;
		attribute refresh operation cannot proceed
		
		# 2010-7-31 	session.expire_on_commit = False doesn't help
		"""
		non_duplicate_array_id_ls = self.getNonDuplicateArraysWithHighestMedianIntensity(db_250k, self.raw_cnv_method_id)
		param_obj = PassingData(no_of_valid_deletions=0, session=db_250k.session, cnv_type_id=self.cnv_type_id, \
					cnv_method_id=self.cnv_method_id, no_of_total=0, no_of_into_db=0, report=self.report, \
					no_of_total_arrays = len(non_duplicate_array_id_ls), minConnectivity=self.minConnectivity)
		self.mergeSegmentsAcrossArrays(db_250k, self.raw_cnv_method_id, min_overlap_ratio=self.min_overlap_ratio,\
										param_obj=param_obj, clearOutGraphAndSave=self.clearOutGraphAndSave)
		session.flush()
		session.expunge_all()
		if self.commit:
			session.commit()
	
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = CNVMergeAcrossArrays
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
