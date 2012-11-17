#!/usr/bin/env python
"""

Examples:
	# 2010-7-29 merge deletions from CNVMethod 12 to 13 (CNVCall)
	~/script/variation/src/CNVMerge.py -u yh -a 12 -m 13 -z banyan.usc.edu -c
	
	# 2010-11-24 merge deletions from CNVQCCall into CNVQCCall
	~/script/variation/src/CNVMerge.py -u yh -z banyan.usc.edu -c -a 9 -m 31 -n 2
	
Description:
	Merge adjacent deletions within one array. After two adjacent deletions are merged,
		neighbors within 5kb will be checked if they could be merged with the newly-merged deletion.
	
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



class CNVMerge(object):
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
						('max_gap_ratio', 1, float): [0.2, 'x', 1, 'maximum min(gap/segment1_len, gap/segment2_len) for two segments to be merged.'],\
						('max_gap_len', 1, float): [200, 'g', 1, 'maximum gap length for two segments to be merged. OR with max_gap_ratio'],\
						('maxNeighborDist', 1, float): [5000, '', 1, 'maximum distance for one deletion to be neighbor of the other'],\
						('maxDeletionLength', 1, float): [5000000, 'D', 1, 'maximum deletion length allowed. two deletions will not be merged if their combined length exceeds this'],\
						('run_type', 0, int):[1, 'n', 1, '1: work on CNVCall, 2: work on CNVQCCall'],\
						('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2009-2-12
		"""
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
												class_to_have_attr=self)
	
	@classmethod
	def mergeTwoSegments(cls, segment, next_segment):
		"""
		2010-7-29
			modified from CNV.LerContig.mergeTwoSegments of variation/src/misc.py
			
			segment = [row.chromosome, row.start, row.stop, row.start_probe_id, row.stop_probe_id, row.no_of_probes_covered, \
					row.size_affected,\
					row.amplitude, row.probability]
		"""
		import numpy
		if segment[2]<=next_segment[2]:		# next_segment ends behind segment
			stop_segment = next_segment
			if segment[5]!=None and next_segment[5]!=None:
				no_of_probes_covered = segment[5] + next_segment[5]	#ignore the case that they are overlapping and the middle gap
			else:
				no_of_probes_covered = None
		else:	#next_segment ends before segment. which means segment wholly encompass next_segment.
			stop_segment = segment
			no_of_probes_covered = segment[5]
		size_affected = stop_segment[2] - segment[1] + 1
		
		# deal with amplitude
		if segment[7]!=None and next_segment[7]!=None:
			amplitude = numpy.mean([segment[7], next_segment[7]])	#mean of two. not very accurate.
		else:
			amplitude = None
		#deal with score/probability
		if segment[-1] is not None and next_segment[-1] is not None:
			score = numpy.mean([segment[-1], next_segment[-1]])
		else:
			score = None
		
		new_segment = [segment[0], segment[1], stop_segment[2], segment[3], stop_segment[4],\
					no_of_probes_covered, size_affected, amplitude, score, ]
		return new_segment
	
	@classmethod
	def calGapRatioOfTwoSegments(cls, segment, next_segment):
		"""
		2010-7-29
			copied from CNV.LerContig of variation/src/misc.py
			
			segment = [row.chromosome, row.start, row.stop, row.start_probe_id, row.stop_probe_id, row.no_of_probes_covered, \
					row.size_affected,\
					row.amplitude, row.probability]
		"""
		current_chr = segment[0]
		current_start = segment[1]
		current_end = segment[2]
			
		next_chr = next_segment[0]
		next_start = next_segment[1]
		next_end = next_segment[2]
		
		gap_len = None
		mergedLength = None
		if next_chr==current_chr:
			gap_len = next_start - current_end - 1
			len1 = current_end - current_start + 1
			len2 = next_end - next_start + 1
			mergedLength = next_end - current_start + 1
			perc1 = gap_len/float(len1)
			perc2 = gap_len/float(len2)
			gap_ratio = min(perc1, perc2)
		else:
			gap_ratio = None
		return (gap_len, gap_ratio, mergedLength,)
	
	@classmethod
	def isTwoSegmentsMergable(cls, segment, next_segment, max_reciprocal_gap_ratio=0.1, max_gap_len=10000,
							maxDeletionLength=50000):
		"""
		2010-6-11
		"""
		current_chr = segment[0]
		next_chr = next_segment[0]
		
		mergable = False
		gap_ratio = None
		if next_chr==current_chr:
			gap_len, gap_ratio, mergedLength = cls.calGapRatioOfTwoSegments(segment, next_segment)[:3]
			if mergedLength<=maxDeletionLength:
				if gap_len<=0:
					#+1 because if they are next to each other, merge them as well
					mergable = True
					gap_ratio = 0
				elif gap_ratio<=max_reciprocal_gap_ratio and gap_len<=max_gap_len:	#2010-8-19 flip the "OR" to "AND"
					mergable = True
		return (gap_ratio, mergable)
	
	@classmethod
	def findChromosomeNeighbor(cls, i, segment_ls, direction=1, maxNeighborDist=10000, oneNeighborEnough=True):
		"""
		2010-7-29
			add argument oneNeighborEnough. jump out if one adjacent neighbor is found.
		2010-6-11
			finding neighbors of segment_ls[i] on the chromosome in certain direction within maxNeighborDist.
			max direction
				-1: go left
				+1: go right
		"""
		segment = segment_ls[i]
		current_chr = segment[0]
		current_start = segment[1]
		current_end = segment[2]
		
		count = 0
		next_segment = None
		no_of_segments = len(segment_ls)
		neighbor_index_ls = []
		while 1:
			count += 1
			new_index = i + count*direction
			if new_index<0 or new_index>=no_of_segments:
				next_segment = None
				break
			next_segment = segment_ls[new_index]
			if next_segment:
				next_chr = next_segment[0]
				next_start = next_segment[1]
				next_end = next_segment[2]
				if next_chr!=current_chr:	#beyond this chromosome, no neighbor. break
					next_segment = None
					break
				
				if direction>0:	#to the right 
					gap_len = (next_start - current_end) - 1
				else:	#to the left
					gap_len = current_start - next_end - 1
				if gap_len>maxNeighborDist:	# too far. stop here
					next_segment = None
					break
				neighbor_index_ls.append(new_index)
				
				if oneNeighborEnough and len(neighbor_index_ls)>=1:	#2010-7-29
					break
		return neighbor_index_ls
		
		
	@classmethod
	def mergeOverlappingORCloseSegmentsByGraph(cls, segment_ls, max_reciprocal_gap_ratio=0.1, max_gap_len=100, \
									mergeFunc=None, maxDeletionLength=50000, maxNeighborDist=5000):
		"""
		2010-7-29
			copied from CNV.LerContig. modified to merge only adjacent segments
		2010-6-11
			a graph-based counterpart to mergeOverlappingORCloseMummerCoords()
			
			
			segment = [row.chromosome, row.start, row.stop, row.start_probe_id, row.stop_probe_id, row.size_affected,\
					row.amplitude, row.probability]
		"""
		sys.stderr.write("Merging overlapping or close-by %s segments by graph ... \n"%(len(segment_ls)))
		if mergeFunc is None:
			mergeFunc=cls.mergeTwoSegments
		segment_ls.sort()
		no_of_segments = len(segment_ls)
		import networkx as nx
		G=nx.Graph()
		
		edgeHeap = []
		sys.stderr.write("Constructing graph ... \n")
		counter  = 0
		real_counter = 0
		for i in xrange(no_of_segments-1):
			segment = segment_ls[i]
			current_chr = segment[0]
			current_start = segment[1]
			current_end = segment[2]
			G.add_node(i)
			
			next_segment = segment_ls[i+1]
			next_chr = next_segment[0]
			next_start = next_segment[1]
			next_end = next_segment[2]
			gap_ratio, mergable = cls.isTwoSegmentsMergable(segment, next_segment, \
								max_reciprocal_gap_ratio=max_reciprocal_gap_ratio, max_gap_len=max_gap_len,\
								maxDeletionLength=maxDeletionLength)[:2]
			if mergable:
				G.add_edge(i, i+1, weight=gap_ratio)
				#heappush(edgeHeap, (gap_ratio, i,j))
				real_counter += 1
			counter += 1
			if counter%1000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
		sys.stderr.write(" . %s nodes, %s edges. Done.\n"%(G.number_of_nodes(), G.number_of_edges()))
		
		sys.stderr.write("Random graph algorithm to merge segments ...\n")
		counter = 0
		real_counter = 0
		import random
		while G.number_of_edges()>0:
			no_of_edges = G.number_of_edges()
			edge_index = random.randint(1, no_of_edges)
			#edge_index = 1
			i, j = G.edges()[edge_index-1]
			try:
				neighbor_set = set(G.neighbors(i))
				neighbor_set |= set(G.neighbors(j))
				
				neighbor_set |= set(cls.findChromosomeNeighbor(i, segment_ls, direction=-1, maxNeighborDist=maxNeighborDist))	#5kb is hard-coded number
				
				# 2010-7-30 only to merge adjacent segments. no need to look at the right side of i.
				#neighbor_set |= set(cls.findChromosomeNeighbor(i, segment_ls, direction=+1, maxNeighborDist=maxNeighborDist))
				# 2010-7-30 only to merge adjacent segments. no need to look at the left side of i.
				#neighbor_set |= set(cls.findChromosomeNeighbor(j, segment_ls, direction=-1, maxNeighborDist=maxNeighborDist))
				
				neighbor_set |= set(cls.findChromosomeNeighbor(j, segment_ls, direction=+1, maxNeighborDist=maxNeighborDist))	#5kb is hard-coded number
				
				segment = segment_ls[i]
				next_segment = segment_ls[j]
				
				segment = mergeFunc(segment, next_segment)
				segment_ls[i] = segment
				segment_ls[j] = None
				
				# remove the two old nodes and affiliated edges
				G.remove_node(i)
				G.remove_node(j)
				
				# add the merged node back
				G.add_node(i)
				
				#i or j could be added when they are next to each other on chromosome
				neighbor_set.remove(i)
				neighbor_set.remove(j)
				
				for node in neighbor_set:
					if i<node:
						first_segment_index = i
						snd_segment_index = node
					elif i>node:
						first_segment_index = node
						snd_segment_index = i
					else:
						continue
					segment = segment_ls[first_segment_index]
					next_segment = segment_ls[snd_segment_index]
					gap_ratio, mergable = cls.isTwoSegmentsMergable(segment, next_segment, \
									max_reciprocal_gap_ratio=max_reciprocal_gap_ratio, max_gap_len=max_gap_len,\
									maxDeletionLength=maxDeletionLength)[:2]
					if mergable:
						G.add_edge(first_segment_index, snd_segment_index, weight=gap_ratio)
						#heappush(edgeHeap, (gap_ratio, first_segment_index, snd_segment_index))
			except:
				sys.stderr.write('Except type at counter=%s, real_counter=%s: %s\n'%(counter, real_counter, repr(sys.exc_info())))
				import traceback
				traceback.print_exc()
				for eH in edgeHeap:
					print eH
				print "%s nodes:"%(G.number_of_nodes())
				for node in G.nodes():
					print node
				print "%s edges:"%(G.number_of_edges())
				for edge in G.edges():
					print edge
				raise
			counter  += 1
			if counter%1000==0:
				sys.stderr.write("%s %s nodes %s edges"%('\x08'*80, G.number_of_nodes(), G.number_of_edges()))
		sys.stderr.write("%s %s nodes %s edges. Done.\n"%('\x08'*80, G.number_of_nodes(), G.number_of_edges()))
		
		return_ls = []
		for segment in segment_ls:
			if segment is not None:
				return_ls.append(segment)
		sys.stderr.write("shrinked to %s segments. Done.\n"%(len(return_ls)))
		return return_ls
	

	def getNonDuplicateArraysWithHighestMedianIntensity(self, db_250k, raw_cnv_method_id=None,):
		"""
		2010-7-29
		"""
		sys.stderr.write("Getting non-duplicate arrays with highest median intensity ...")
		rows = db_250k.metadata.bind.execute("select a.* from view_array a, (select distinct array_id from %s where cnv_method_id=%s) \
				t where t.array_id=a.array_id"%(Stock_250kDB.CNVCall.table.name, raw_cnv_method_id))
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
	
	def mergeSegmentsForOneArray(self, db_250k, array_id, raw_cnv_method_id = None, max_gap_ratio=0.3, \
							max_gap_len=None, maxDeletionLength=50000, param_obj=None):
		"""
		2010-7-29
		"""
		sys.stderr.write("Merging segments for array %s ... \n"%array_id)
		query = Stock_250kDB.CNVCall.query.filter_by(cnv_method_id=raw_cnv_method_id).\
			filter_by(array_id=array_id).filter_by(cnv_type_id=param_obj.cnv_type_id).\
			order_by(Stock_250kDB.CNVCall.chromosome).order_by(Stock_250kDB.CNVCall.start).order_by(Stock_250kDB.CNVCall.stop)
		
		segment_ls = []
		for row in query:
			segment = [row.chromosome, row.start, row.stop, row.start_probe_id, row.stop_probe_id, row.no_of_probes_covered, \
					row.size_affected,\
					row.amplitude, row.probability]
			segment_ls.append(segment)
		
		merged_segment_ls = self.mergeOverlappingORCloseSegmentsByGraph(segment_ls, max_reciprocal_gap_ratio=max_gap_ratio, \
																max_gap_len=max_gap_len, mergeFunc=self.mergeTwoSegments,\
																maxDeletionLength=maxDeletionLength,\
																maxNeighborDist=getattr(param_obj, 'maxNeighborDist', 5000))
		
		from CNVPredictDeletionBySVM import CNVPredictDeletionBySVM
		for merged_segment in merged_segment_ls:
			chromosome, start, stop, start_probe_id, stop_probe_id, no_of_probes_covered, size_affected,\
					amplitude, probability = merged_segment[:9]
			cnv_segment_obj = PassingData(array_id=array_id, start_probe_id=start_probe_id, stop_probe_id=stop_probe_id, \
								no_of_probes=no_of_probes_covered, amplitude=amplitude, segment_length=size_affected, \
								segment_chromosome=chromosome, \
								segment_start_pos=start, segment_stop_pos=stop, \
								median_intensity=None, probability=probability)
			CNVPredictDeletionBySVM.saveSegmentObj(param_obj, cnv_segment_obj)
		sys.stderr.write("Done.\n")
	
	def mergeSegmentsForOneQCAccession(self, db_250k, accession_id, raw_cnv_method_id = None, max_gap_ratio=0.3, \
							max_gap_len=None, maxDeletionLength=50000, param_obj=None):
		"""
		2010-11-24
			merging segments for CNVQCCall
		"""
		sys.stderr.write("Merging segments for accession %s ... \n"%accession_id)
		query = Stock_250kDB.CNVQCCall.query.filter_by(cnv_method_id=raw_cnv_method_id).\
			filter_by(accession_id=accession_id).filter_by(cnv_type_id=param_obj.cnv_type_id).\
			order_by(Stock_250kDB.CNVQCCall.chromosome).order_by(Stock_250kDB.CNVQCCall.start).order_by(Stock_250kDB.CNVQCCall.stop)
		
		segment_ls = []
		start_probe_id = None
		stop_probe_id = None
		for row in query:
			amplitude = row.copy_number
			segment = [row.chromosome, row.start, row.stop, start_probe_id, stop_probe_id, row.no_of_probes_covered, \
					row.size_affected, amplitude, row.score]
			segment_ls.append(segment)
		
		merged_segment_ls = self.mergeOverlappingORCloseSegmentsByGraph(segment_ls, max_reciprocal_gap_ratio=max_gap_ratio, \
																	max_gap_len=max_gap_len, mergeFunc=self.mergeTwoSegments,\
																	maxDeletionLength=maxDeletionLength,\
																	maxNeighborDist=getattr(param_obj, 'maxNeighborDist', 5000))
		
		for merged_segment in merged_segment_ls:
			chromosome, start, stop, start_probe_id, stop_probe_id, no_of_probes_covered, size_affected,\
					amplitude, probability = merged_segment[:9]
			copy_number = amplitude
			cnv_segment_obj = PassingData(accession_id=accession_id, start_probe_id=start_probe_id, stop_probe_id=stop_probe_id, \
								no_of_probes_covered=no_of_probes_covered, copy_number=copy_number, segment_length=size_affected, \
								segment_chromosome=chromosome, segment_start_pos=start, segment_stop_pos=stop, \
								median_intensity=None, probability=probability)
			self.saveSegmentObjInCNVQCCall(param_obj, cnv_segment_obj)
		sys.stderr.write("Done.\n")
	
	def saveSegmentObjInCNVQCCall(self, param_obj, cnv_segment_obj):
		"""
		2010-11-24
			save the segment object
		"""
		if not hasattr(param_obj, 'no_of_total'):
			setattr(param_obj, 'no_of_total', 0)
		if not hasattr(param_obj, 'no_of_into_db'):
			setattr(param_obj, 'no_of_into_db', 0)
		session = getattr(param_obj, 'session', None)
		if session is None:
			sys.stderr.write("Error: db session is not available.\n")
			return
		param_obj.no_of_total += 1
		
		cnv_method_id = getattr(param_obj, "cnv_method_id", None)
		cnv_type_id = getattr(param_obj, "cnv_type_id", None)
		accession_id = cnv_segment_obj.accession_id
		rows = Stock_250kDB.CNVQCCall.query.filter_by(accession_id=accession_id).filter_by(chromosome=cnv_segment_obj.segment_chromosome).\
				filter_by(start = cnv_segment_obj.segment_start_pos).filter_by(stop = cnv_segment_obj.segment_stop_pos).\
				filter_by(cnv_type_id = cnv_type_id).\
				filter_by(cnv_method_id = cnv_method_id)
		if rows.count()==0:	#make sure it's not in db yet.
			probability = getattr(cnv_segment_obj, 'probability', None)
			cnv_qc_call = Stock_250kDB.CNVQCCall(accession_id=accession_id, chromosome = cnv_segment_obj.segment_chromosome,\
										start = cnv_segment_obj.segment_start_pos, stop = cnv_segment_obj.segment_stop_pos, \
										no_of_probes_covered = cnv_segment_obj.no_of_probes_covered, \
										size_affected = cnv_segment_obj.segment_length,\
										copy_number = cnv_segment_obj.copy_number,\
										score = probability)
			cnv_qc_call.cnv_method_id = cnv_method_id
			cnv_qc_call.cnv_type_id = cnv_type_id
			cnv_qc_call.comment = getattr(cnv_segment_obj, 'comment', None)	#2010-7-25
			session.add(cnv_qc_call)
			param_obj.no_of_into_db += 1
		if param_obj.no_of_into_db%10000==0:
			session.flush()
			session.expunge_all()
		
		report = getattr(param_obj, 'report', False)
		if report and param_obj.no_of_total%10000==0:
			sys.stderr.write('%s%s into db out of %s\n'%('\x08'*100, param_obj.no_of_into_db, param_obj.no_of_total))
	
	
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
		
		param_obj = PassingData(no_of_valid_deletions=0, session=db_250k.session, cnv_type_id=self.cnv_type_id, \
						cnv_method_id=self.cnv_method_id, no_of_total=0, no_of_into_db=0, report=self.report,\
						maxNeighborDist=self.maxNeighborDist)
		if self.run_type==1:
			non_duplicate_array_id_ls = self.getNonDuplicateArraysWithHighestMedianIntensity(db_250k, self.raw_cnv_method_id)
			for array_id in non_duplicate_array_id_ls:
				self.mergeSegmentsForOneArray(db_250k, array_id, self.raw_cnv_method_id, max_gap_ratio=self.max_gap_ratio,\
											max_gap_len=self.max_gap_len, maxDeletionLength=self.maxDeletionLength,\
											param_obj=param_obj)
		elif self.run_type==2:
			rows = db_250k.metadata.bind.execute("select distinct c.accession_id from %s c where c.cnv_method_id=%s and c.cnv_type_id=%s"\
				%(Stock_250kDB.CNVQCCall.table.name, self.raw_cnv_method_id, self.cnv_type_id))
			for row in rows:
				self.mergeSegmentsForOneQCAccession(db_250k, row.accession_id, raw_cnv_method_id =self.raw_cnv_method_id, \
							max_gap_ratio=self.max_gap_ratio, max_gap_len=self.max_gap_len, \
							maxDeletionLength=self.maxDeletionLength, param_obj=param_obj)
		session.flush()
		session.expunge_all()
		if self.commit:
			session.commit()
		"""
		2009-2-12
		Program to 
			1. discretize CNV amplitude file outputted by GADA into 1 (normal & duplication) or -1 (deletion).
			2. make segmental plots in which each plot contains consecutive 1000 probes. y-axis is probe value summarizing across all accessions.
				input could be output of RunGADA.py or from step 1.
			
		import pylab
		cnvIntensityData = SNPData(input_fname=self.input_fname, turn_into_array=1, ignore_2nd_column=1, matrix_data_type=float)
		probe_pos_ls = []
		avg_intensity_ls = []
		
		if self.run_type == 1:
			newDataMatrix = numpy.ones(cnvIntensityData.data_matrix.shape, numpy.int)
		
		for j in range(cnvIntensityData.data_matrix.shape[1]):
			probe_id = cnvIntensityData.col_id_ls[j]
			probe_id = probe_id.split('_')
			probe_id = map(int, probe_id)
			probe_pos_ls.append(probe_id[1])
			avg_intensity_ls.append(numpy.sum(cnvIntensityData.data_matrix[:,j]))
			if self.run_type==1:
				for i in range(cnvIntensityData.data_matrix.shape[0]):
					if cnvIntensityData.data_matrix[i][j]<=self.max_del_intensity:
						newDataMatrix[i][j] = -1
		
		if self.run_type==1:
			newData = SNPData(row_id_ls=cnvIntensityData.row_id_ls, col_id_ls=cnvIntensityData.col_id_ls, data_matrix=newDataMatrix)
			newData.tofile(self.output_fname)
		elif self.run_type==2:
			block_size = 1000
			no_of_probes = len(probe_pos_ls)
			no_of_blocks = no_of_probes/block_size
			for i in range(no_of_blocks):
				if i*block_size>no_of_probes:
					break
				start_index = i*block_size
				end_index = min((i+1)*block_size, no_of_probes)
				fname = '%s_%s_%s.png'%(self.output_fname, probe_pos_ls[start_index], probe_pos_ls[end_index])
				pylab.clf()
				pylab.plot(probe_pos_ls[start_index:end_index], avg_intensity_ls[start_index:end_index], '.', markersize=4, alpha=0.4)
				pylab.xlabel('chromosome position')
				pylab.ylabel('sum intensity')
				pylab.savefig(fname, dpi=300)
		"""
	
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = CNVMerge
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
