#!/usr/bin/env python
"""

Examples:
	# define landscape for result 4634.
	DefineAssociationLandscape.py -z banyan -u yh -e 4634
	
Description:
	Program to find the landscape/peaks of a genome wide association result.
	It fills data into db table ResultPeak, but not table ResultLandscape.

"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, getopt
import warnings, traceback
from pymodule import ProcessOptions, PassingData, GenomeDB
import Stock_250kDB
from pymodule import SNP, getListOutOfStr, yh_matplotlib, PassingData
from pymodule.db import formReadmeObj
from pymodule.CNV import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
from pymodule.RBTree import RBDict
import networkx as nx

import CGAL	#computational geometry algorithm python binding

class DefineAssociationLandscape(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['stock_250k', 'd', 1, 'stock_250k database name', ],\
						('genome_dbname', 1, ): ['genome', 'g', 1, 'genome database name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('results_directory', 0, ):[None, 't', 1, 'The results directory. Default is None. use the one given by db.'],\
						('result_id_ls', 1, ): [None, 'e', 1, 'comma or dash-separated list of result ids, i.e. 3431-3438,4041'],\
						('neighbor_distance', 0, int): [5000, '', 1, "within this distance, a locus that increases the association score \
									the fastest is chosen as bridge end. outside this distance, whatever the next point is will be picked."],\
						('max_neighbor_distance', 0, int): [20000, '', 1, "beyond this distance, no bridge would be formed."],\
						('min_MAF', 0, float): [0.1, 'n', 1, 'minimum Minor Allele Frequency.'],\
						('min_score', 0, float): [4, 'f', 1, 'minimum score to call a peak'],\
						('ground_score', 0, float): [0, '', 1, 'minimum score possible in this test'],\
						('tax_id', 1, int): [3702, 'x', 1, 'to get the number of total genes from database, which species.'],\
						('output_fname', 0, ): [None, 'o', 1, 'if given, output the landscape result.'],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-08-19
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
	
		if getattr(self, 'result_id_ls', None):
			self.result_id_ls = getListOutOfStr(self.result_id_ls, data_type=int)
			self.result_id_ls.sort()
		else:
			self.result_id_ls = []
	
	
	def getXDistanceBetweenTwoDataObjs(self, current_obj=None, next_obj=None):
		"""
		2011-4-20
			the two arguments are of type DataObject in pymodule/SNP.py.
			It calculates the distance from current_obj to next_obj.
		"""
		
		deltaX = next_obj.position - current_obj.stopPosition
		return deltaX
	
	def findLandscape(self, data_obj_ls=None, neighbor_distance=5000, max_neighbor_distance=20000,\
					ground_score=0):
		"""
		2011-4-19
			add max_neighbor_distance, ground_score
		2011-4-18
			
		"""
		sys.stderr.write("Finding landscape for %s data objects ..."%(len(data_obj_ls)))
		current_obj = None
		bridge_ls = []
		start_index = 0	#starting index in each step of bridge building
		locusLandscapeNeighborGraph = nx.Graph()
		while start_index<len(data_obj_ls):
			current_obj = data_obj_ls[start_index]
			current_obj.index = start_index
			obj_with_fastest_score_increase = None
			for j in xrange(start_index + 1, len(data_obj_ls)):
				next_obj = data_obj_ls[j]
				next_obj.index = j
				if current_obj.chromosome==next_obj.chromosome:
					deltaY = next_obj.value - current_obj.value
					deltaX = self.getXDistanceBetweenTwoDataObjs(current_obj, next_obj)
					#watch the distance calculation
					rate = deltaY/deltaX
					next_obj.rate = rate
					if deltaX<=max_neighbor_distance:
						if obj_with_fastest_score_increase is None:	#regardless of the distance and rate, pick this one.
							obj_with_fastest_score_increase = next_obj
						else:
							if deltaX>neighbor_distance:	#this object is over the minimum distance
								break
							if next_obj.rate>=obj_with_fastest_score_increase.rate:
								#if it's equal, still pick this next one because it spans further.
								obj_with_fastest_score_increase = next_obj
					else:	#over the max_neighbor_distance and data objects are ordered. so safely break
						break
				else:
					break
			if obj_with_fastest_score_increase is not None:
				no_of_loci = obj_with_fastest_score_increase.index - current_obj.index - 1
				bridge_ls.append([current_obj, obj_with_fastest_score_increase, no_of_loci, deltaX])
				source_index = current_obj.index
				target_index = obj_with_fastest_score_increase.index
				locusLandscapeNeighborGraph.add_edge(source_index, target_index, \
										weight=obj_with_fastest_score_increase.rate)
				locusLandscapeNeighborGraph[source_index][target_index]['no_of_loci'] = no_of_loci
				deltaX = self.getXDistanceBetweenTwoDataObjs(current_obj, obj_with_fastest_score_increase)
				locusLandscapeNeighborGraph[source_index][target_index]['deltaX'] = deltaX
				start_index = obj_with_fastest_score_increase.index		#new starting point
			else:
				
				start_index = start_index+1	#no bridge is built. at chromosome boundary.
		
		sys.stderr.write("%s bridges. Done.\n"%(len(bridge_ls)))
		
		returnData = PassingData(bridge_ls=bridge_ls, locusLandscapeNeighborGraph=locusLandscapeNeighborGraph)
		return returnData
	
	def drawBridgeChromosomalLengthHist(self, bridge_ls):
		"""
		2011-4-18
		"""
		no_of_bridges = len(bridge_ls)
		sys.stderr.write("Drawing histogram of chromosomal length for %s bridges ... \n"%(no_of_bridges))
		bridge_chr_length_ls = []
		no_of_loci_per_bridge_ls = []
		for i in xrange(no_of_bridges):
			bridge = bridge_ls[i]
			bridge_chr_length_ls.append(bridge[3])
			no_of_loci_per_bridge_ls.append(bridge[2])
			
		yh_matplotlib.drawHist(bridge_chr_length_ls, title='Histogram of bridge chromosomal length', \
							xlabel_1D='chromosomal length',\
							outputFname='/tmp/chromosomal_length_hist.png', min_no_of_data_points=50, needLog=True)
		yh_matplotlib.drawHist(no_of_loci_per_bridge_ls, title='Histogram of no-of-loci per bridge', \
							xlabel_1D='no-of-loci',\
							outputFname='/tmp/no_of_loci_hist.png', min_no_of_data_points=50, needLog=True)
		sys.stderr.write("Done.\n")
	
	def addXYFromObject(self, data_obj, x_ls=None, y_ls=None, cumu_start=0):
		"""
		2011-4-21
			used in drawBridgeLs()
		"""
		x_ls.append(cumu_start + data_obj.position)
		y_ls.append(data_obj.value)
		if data_obj.stop_position:	#not None
			x_ls.append(cumu_start + data_obj.stop_position)
			y_ls.append(data_obj.value)
		
	
	def drawBridgeLs(self, bridge_ls=[], outputFname=None, oneGenomeData=None):
		"""
		2011-4-20
			draw a plot showing the landscape in bridge_ls, returned by findLandscape()
		"""
		sys.stderr.write("Drawing a whole genome plot for %s bridges ...\n"%(len(bridge_ls)))
		import pylab
		pylab.clf()
		pylab.figure(figsize=(25, 4))
		color_ls = ['r','b']
		chr_id_ls = oneGenomeData.chr_id2cumu_start.keys()
		chr_id_ls.sort()
		for bridge in bridge_ls:
			current_obj, next_obj = bridge[:2]
			oneGenomeData.chr_id2cumu_start
			try:
				chr_id_index = chr_id_ls.index(current_obj.chromosome)
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
				continue
			color = color_ls[chr_id_index%2]
			cumu_start = oneGenomeData.chr_id2cumu_start.get(current_obj.chromosome)
			x_ls = []
			y_ls = []
			self.addXYFromObject(current_obj, x_ls=x_ls, y_ls=y_ls, cumu_start=cumu_start)
			self.addXYFromObject(next_obj, x_ls=x_ls, y_ls=y_ls, cumu_start=cumu_start)
			pylab.plot(x_ls, y_ls, color=color)
			
		pylab.savefig(outputFname, dpi=300)
		sys.stderr.write("Done.\n")
	
	def findPeakDataObjWithinPeak(self, data_obj_ls=[], connected_component_node_list=[], peak_start_data_obj=None, \
							peak_stop_data_obj=None):
		"""
		2011-4-20
			given a peak starting from peak_start_data_obj and ending at peak_stop_data_obj, which locus in the middle
			has the maximum score. return the first one if there are multiple ones.
		"""
		data_obj_with_max_value = None
		for i in xrange(1, len(connected_component_node_list)-1):
			data_obj_index = connected_component_node_list[i]
			if data_obj_index>=peak_start_data_obj.index and data_obj_index<=peak_stop_data_obj.index:
				data_obj = data_obj_ls[data_obj_index]
				if data_obj_with_max_value is None:
					data_obj_with_max_value = data_obj
				elif data_obj.value>data_obj_with_max_value.value:
					data_obj_with_max_value=data_obj
		return data_obj_with_max_value
	
	
	def findPeaks(self, db_250k, locusLandscapeNeighborGraph, bridge_ls=None, data_obj_ls=None, ground_score=0, min_score=4,\
				rm=None, result_peak_type=None):
		"""
		2011-4-19
		"""
		no_of_bridges = len(bridge_ls)
		sys.stderr.write("Finding peaks among the %s bridges ..."%(no_of_bridges))
		no_of_edges = locusLandscapeNeighborGraph.number_of_edges()
		if no_of_edges!=no_of_bridges:
			sys.stderr.write("Warning: number of edges % in the graph is different from the number of bridges %s.\n"%\
							(no_of_edges, no_of_bridges))
		cc_list = nx.connected_components(locusLandscapeNeighborGraph)
		no_of_peaks = 0
		for cc in cc_list:
			peak_start_data_obj = None
			cc.sort()	#each node is identified by index in data_obj_ls
			cc.insert(0, -1)	#-1 is fake index. standing for the ground point
			cc.append(-1)	#-1 is fake index. standing for the ground point
			no_of_nodes = len(cc)
			for i in xrange(no_of_nodes-1):
				start_data_obj_index = cc[i]
				stop_data_obj_index = cc[i+1]
				if start_data_obj_index!=-1:
					start_data_obj = data_obj_ls[start_data_obj_index]
				else:
					start_data_obj = None
				bridge_x_start = getattr(start_data_obj, 'stopPosition',None)
				bridge_y_start = getattr(start_data_obj, 'value', None)
				if stop_data_obj_index!=-1:
					stop_data_obj = data_obj_ls[stop_data_obj_index]
				else:
					stop_data_obj = None
				bridge_x_stop = getattr(stop_data_obj, 'position', None)
				bridge_y_stop = getattr(stop_data_obj, 'value', None)
				
				if bridge_x_start is None:	#start_data_obj_index is -1
					bridge_x_start = bridge_x_stop
					bridge_y_start = ground_score
				if bridge_x_stop is None:	#stop_data_obj_index is -1
					bridge_x_stop = bridge_x_start
					bridge_y_stop = ground_score
				
				p1 = CGAL.Kernel.Point_2(bridge_x_start, bridge_y_start)
				p2 = CGAL.Kernel.Point_2(bridge_x_stop, bridge_y_stop)
				segment = CGAL.Kernel.Segment_2(p1, p2)
				
				cutoff_segment = CGAL.Kernel.Segment_2(CGAL.Kernel.Point_2(bridge_x_start, min_score),\
										CGAL.Kernel.Point_2(bridge_x_stop, min_score))
				intersection_point = CGAL.intersectionPointOfTwoSegments(segment, cutoff_segment)
				if hasattr(intersection_point, 'is_degenerate'):
					pass	#it's an empty segment. no intersection.
				else:
					if peak_start_data_obj is None:
						locus_index = stop_data_obj.index	#backwards starting from the stop_data_obj
						# start_data_obj could be None.
						while locus_index>=0:
							data_obj = data_obj_ls[locus_index]
							stop_position = getattr(data_obj, 'stopPosition', None)
							if stop_position<intersection_point.x():
								break
							locus_index -= 1
						peak_start_data_obj = data_obj_ls[locus_index+1]	#the loop breaks after the outside locus has been reached.
						peak_start_data_obj.peak_start = intersection_point.x()	#record the peak starting position
					else:
						locus_index = start_data_obj.index	#forward starting from the start_data_obj
						# start_data_obj could contain something.
						while locus_index<len(data_obj_ls):
							data_obj = data_obj_ls[locus_index]
							start_position = getattr(data_obj, 'position')
							if start_position>intersection_point.x():
								break
							locus_index += 1
						peak_stop_data_obj = data_obj_ls[locus_index-1]
						no_of_loci = peak_stop_data_obj.index - peak_start_data_obj.index + 1
						peak_data_obj = self.findPeakDataObjWithinPeak(data_obj_ls=data_obj_ls, connected_component_node_list=cc, \
												peak_start_data_obj=peak_start_data_obj, peak_stop_data_obj=peak_stop_data_obj)
						result_peak = Stock_250kDB.ResultPeak(result_id = rm.id, chromosome=peak_start_data_obj.chromosome,\
								start=peak_start_data_obj.peak_start, stop=intersection_point.x(), no_of_loci=no_of_loci, \
								peak_score=peak_data_obj.value)
						if rm.cnv_method_id:
							# 2011-4-21 assuming locus might not be in db.
							result_peak.start_locus = self.getLocusBasedOnDataObj(db_250k, peak_start_data_obj)
							result_peak.stop_locus = self.getLocusBasedOnDataObj(db_250k, peak_stop_data_obj)
							result_peak.peak_locus = self.getLocusBasedOnDataObj(db_250k, peak_data_obj)
						elif rm.call_method_id:
							# 2011-4-21 assuming all relevant loci are in db.
							result_peak.start_locus_id = peak_start_data_obj.db_id
							result_peak.stop_locus_id = peak_stop_data_obj.db_id
							result_peak.peak_locus_id = peak_data_obj.db_id
						
						result_peak.result_peak_type = result_peak_type
						
						db_250k.session.add(result_peak)
						peak_start_data_obj = None	#reset for the next peak
						no_of_peaks += 1
			db_250k.session.flush()
			#subgraph = nx.subgraph(locusLandscapeNeighborGraph, cc)
		sys.stderr.write("%s peaks.\n"%(no_of_peaks))
	
	def getLocusBasedOnDataObj(self, db_250k, data_obj):
		"""
		2011-4-21
			based on attributes of data_obj (class SNP.DataObject), fetch/create a locus.
		"""
		return db_250k.getSNP(chromosome=data_obj.chromosome, start=data_obj.position, stop=data_obj.stop_position)
	
	def addAllDataObjsIntoDb(self, db_250k, data_obj_ls): 
		'''
		2011-4-21
			call getLocusBasedOnDataObj() to make sure all loci are in db now.
			not used now since all data_obj's corresponding loci are stored in db.
		'''
		sys.stderr.write("Adding %s data objs into table Snps as loci ..."%(len(data_obj_ls)))
		no_of_newly_added_loci = 0
		for data_obj in data_obj_ls:
			locus = self.getLocusBasedOnDataObj(db_250k, data_obj)
			if not locus.id:	#if it's new , their id should be None.
				no_of_newly_added_loci += 1
		db_250k.session.flush()
		sys.stderr.write("%s loci inserted into db. Done.\n"%(no_of_newly_added_loci))
	
	def getResultPeakType(self, db_250k, min_score=None, neighbor_distance=None, max_neighbor_distance=None):
		"""
		2011-4-20
		"""
		
		result_peak_type = Stock_250kDB.ResultPeakType.query.filter_by(min_score=min_score).filter_by(neighbor_distance=neighbor_distance).\
			filter_by(max_neighbor_distance=max_neighbor_distance).first()
		if result_peak_type is None:
			result_peak_type = Stock_250kDB.ResultPeakType(min_score=self.min_score, neighbor_distance=self.neighbor_distance, \
						max_neighbor_distance=self.max_neighbor_distance)
			db_250k.session.add(result_peak_type)
			db_250k.session.flush()
		return result_peak_type
	
	def run(self):
		"""
		2011-3-28
			Read in the association result
			
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_250k = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user, password=self.db_passwd, \
									hostname=self.hostname, database=self.dbname)
		db_250k.setup(create_tables=False)
		
		genome_db = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.genome_dbname, )
		genome_db.setup(create_tables=False)
		
		oneGenomeData = genome_db.getOneGenomeData(tax_id=self.tax_id, chr_gap=0)
		#oneGenomeData.chr_id2cumu_start
		#cumuSpan2ChrRBDict = oneGenomeData.cumuSpan2ChrRBDict
		
		#param_obj = PassingData(no_of_total_annotations=0, session=db_250k.session, \
		#			cnv_method_id=self.cnv_method_id, no_of_total_contexts=0, no_of_into_db=0, report=self.report,\
		#			no_of_cnv_contexts_already_in_db=0, no_of_cnv_annotations_already_in_db=0)
		
		pd = PassingData(min_MAF=self.min_MAF,\
					results_directory=self.results_directory, \
					need_chr_pos_ls=0,)
		
		for result_id in self.result_id_ls:
			#establish the map from cnv.id from chr_pos
			rm = Stock_250kDB.ResultsMethod.get(result_id)
			if rm.cnv_method_id and not db_250k._cnv_id2chr_pos:
				db_250k.cnv_id2chr_pos = rm.cnv_method_id
				pd.db_id2chr_pos = db_250k.cnv_id2chr_pos
			elif rm.call_method_id:
				pd.db_id2chr_pos = db_250k.snp_id2chr_pos
			
			gwr = db_250k.getResultMethodContent(result_id, pdata=pd)
			gwr.data_obj_ls.sort(cmp=SNP.cmpDataObjByChrPos)
			#if rm.cnv_method_id:	#for cnvs, need to make sure all loci are in db.
			#	self.addAllDataObjsIntoDb(db_250k, gwr.data_obj_ls)
			landscapeData = self.findLandscape(gwr.data_obj_ls, neighbor_distance=self.neighbor_distance,\
										max_neighbor_distance=self.max_neighbor_distance)
			"""
			#2011-4-21 for inspection
			outputFname = '/tmp/result_%s_landscape.png'%(result_id)
			self.drawBridgeLs(bridge_ls=landscapeData.bridge_ls, outputFname=outputFname, oneGenomeData=oneGenomeData)
			"""
			
			result_peak_type = self.getResultPeakType(db_250k, min_score=self.min_score, neighbor_distance=self.neighbor_distance, \
												max_neighbor_distance=self.max_neighbor_distance)
			self.findPeaks(db_250k, landscapeData.locusLandscapeNeighborGraph, bridge_ls=landscapeData.bridge_ls, data_obj_ls=gwr.data_obj_ls, \
						ground_score=self.ground_score, min_score=self.min_score, rm=rm, result_peak_type=result_peak_type)
			"""
			#2011-4-21 to check how far things are from each other.
			#self.drawBridgeChromosomalLengthHist(bridge_ls)
			"""
		db_250k.session.flush()
		if self.commit:
			db_250k.session.commit()
	
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = DefineAssociationLandscape
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()