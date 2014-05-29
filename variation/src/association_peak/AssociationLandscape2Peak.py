#!/usr/bin/env python
"""

Examples:
	%s
	
	# set the minimum score to 4
	%s  -i  association_min_score_5.0/call_75_analysis_1_min_score_5.0_min_overlap_0.5_result_5371_landscape.h5
		-o  association_min_score_5.0/call_75_analysis_1_min_score_5.0_min_overlap_0.5_result_5371_peak.h5  --min_score 5.0 --ground_score 0.0
	
Description:
	2012.11.18 this program derives association peak from association landscape (output of DefineAssociationLandscape.py)

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import numpy
import CGAL	#computational geometry algorithm python binding
import networkx as nx
from pymodule import ProcessOptions, PassingData, PassingDataList
from pymodule import HDF5MatrixFile
from pymodule import getGenomeWideResultFromHDF5MatrixFile, AbstractMapper
from pymodule import yhio, AssociationTableFile, AssociationLandscapeTableFile, AssociationPeakTableFile
from variation.src import AbstractVariationMapper

class AssociationLandscape2Peak(AbstractVariationMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVariationMapper.option_default_dict.copy()
	#get rid of db-related options
	ProcessOptions.removeCertainOptions(option_default_dict=option_default_dict, option_dict_to_remove=AbstractMapper.db_option_dict)
	#option_default_dict.pop(('min_MAF', 0, float))
	my_option_dict = {
					('min_score', 0, float): [4, 'f', 1, 'minimum score to call a peak'],\
					('ground_score', 0, float): [0, 's', 1, 'minimum score possible in this test'],\
					}
	option_default_dict.update(my_option_dict.copy())
	
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		2008-08-19
		"""
		AbstractVariationMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def connectDB(self):
		"""
		2012.11.20 no db connection required
		"""
		pass
	
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
	
	def findPeaks(self, locusLandscapeNeighborGraph=None, bridge_ls=None, data_obj_ls=None, \
				ground_score=0, min_score=4):
		"""
		2012.3.12
			bugfixing, bounds the peak start  locus and peak stop locus within each connected component (one landscape)
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
		association_peak_ls = []
		for cc in cc_list:
			peak_start_data_obj = None
			cc.sort()	#each node is identified by index in data_obj_ls
			cc_start_obj_index = cc[0]
			cc_stop_obj_index = cc[-1]
			cc.insert(0, -1)	#-1 is fake index. standing for the ground point
			cc.append(-1)	#-1 is fake index. standing for the ground point in the end
			no_of_nodes = len(cc)
			no_of_loci_in_cc = cc_stop_obj_index - cc_start_obj_index + 1	#2012.3.13
			
			for i in xrange(no_of_nodes-1):
				start_data_obj_index = cc[i]
				stop_data_obj_index = cc[i+1]
				if start_data_obj_index!=-1:
					start_data_obj = data_obj_ls[start_data_obj_index]
				else:
					start_data_obj = None
				bridge_start_x = getattr(start_data_obj, 'stop', None)	#failsafe if start_data_obj is None.
				bridge_start_y = getattr(start_data_obj, 'value', None)
				if stop_data_obj_index!=-1:
					stop_data_obj = data_obj_ls[stop_data_obj_index]
				else:
					stop_data_obj = None
				bridge_stop_x = getattr(stop_data_obj, 'start', None)
				bridge_stop_y = getattr(stop_data_obj, 'value', None)
				
				if bridge_start_x is None:	#start_data_obj_index is -1
					bridge_start_x = bridge_stop_x
					bridge_start_y = ground_score
				if bridge_stop_x is None:	#stop_data_obj_index is -1
					bridge_stop_x = bridge_start_x
					bridge_stop_y = ground_score
				
				#cast them into float because CGAL takes float values, not numpy.int64 (data is from HDF5 file)
				bridge_start_x = float(bridge_start_x)
				bridge_start_y = float(bridge_start_y)
				bridge_stop_x = float(bridge_stop_x)
				bridge_stop_y = float(bridge_stop_y)
				
				p1 = CGAL.Kernel.Point_2(bridge_start_x, bridge_start_y)
				p2 = CGAL.Kernel.Point_2(bridge_stop_x, bridge_stop_y)
				segment = CGAL.Kernel.Segment_2(p1, p2)
				
				cutoff_segment = CGAL.Kernel.Segment_2(CGAL.Kernel.Point_2(bridge_start_x, min_score),\
										CGAL.Kernel.Point_2(bridge_stop_x, min_score))
				intersection_point = CGAL.intersectionPointOfTwoSegments(segment, cutoff_segment)
				if hasattr(intersection_point, 'is_degenerate'):
					pass	#it's an empty segment. no intersection.
				else:
					if peak_start_data_obj is None:
						#there is an intersection, but no starting obj for this potential peak.
						# so find this starting object first.
						locus_index = stop_data_obj.index	#backwards starting from the stop_data_obj (which could not be None)
						# start_data_obj could be None.
						while locus_index>=0 and locus_index>=cc_start_obj_index:	#2012.3.12 within the cc.
							data_obj = data_obj_ls[locus_index]
							
							#2011-10-11 bugfix
							#this data_obj could be a single-position locus (SNP) which doesn't have stopPosition defined.
							stop_position = data_obj.stopPosition	#same as getattr(data_obj, 'stopPosition', getattr(data_obj, 'position'))
							
							if stop_position<intersection_point.x():
								break
							locus_index -= 1
						peak_start_data_obj = data_obj_ls[locus_index+1]	#the loop breaks after the outside locus has been reached.
						peak_start_data_obj.peak_start = intersection_point.x()	#record the peak starting position
					else:
						#intersection exists again.
						# this should be the ending position for the peak, now find the ending object (peak_stop_data_obj).
						locus_index = start_data_obj.index	#forward starting from the start_data_obj (which could not be None)
						# start_data_obj could contain something.
						while locus_index<=cc_stop_obj_index:	#2012.3.12 bounded by the last object in the cc
							data_obj = data_obj_ls[locus_index]
							start_position = getattr(data_obj, 'position')
							if start_position>intersection_point.x():
								break
							locus_index += 1
						peak_stop_data_obj = data_obj_ls[locus_index-1]
						no_of_loci = peak_stop_data_obj.index - peak_start_data_obj.index + 1
						if no_of_loci>no_of_loci_in_cc:
							sys.stderr.write("Error: no_of_loci of a peak, %s, could not be larger than no_of_loci within its landscape: %s.\n"%\
											(no_of_loci, no_of_loci_in_cc))
							#import pdb
							#pdb.set_trace()
							sys.exit(4)
						peak_data_obj = self.findPeakDataObjWithinPeak(data_obj_ls=data_obj_ls, connected_component_node_list=cc, \
												peak_start_data_obj=peak_start_data_obj, peak_stop_data_obj=peak_stop_data_obj)
						association_peak = PassingDataList()	#assigning each value separately imposes an order in which they appear in the PassingDataList.
						#which would determine the order when sorting a list of PassingDataList.
						association_peak.chromosome=peak_start_data_obj.chromosome
						association_peak.start=peak_start_data_obj.peak_start
						association_peak.stop=intersection_point.x()
						association_peak.no_of_loci=no_of_loci
						association_peak.peak_score=peak_data_obj.value
						association_peak.start_locus_id=peak_start_data_obj.db_id
						association_peak.stop_locus_id=peak_stop_data_obj.db_id
						association_peak.peak_locus_id= peak_data_obj.db_id
						association_peak_ls.append(association_peak)
						peak_start_data_obj = None	#reset for the next peak
						no_of_peaks += 1
			#subgraph = nx.subgraph(locusLandscapeNeighborGraph, cc)
		sys.stderr.write("%s peaks.\n"%(no_of_peaks))
		return association_peak_ls

	
	def run(self):
		"""
		2011-3-28
			Read in the association result
			
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		associationLandscapeFile = AssociationLandscapeTableFile(self.inputFname, openMode='r')
		landscapeTable = associationLandscapeFile.associationLandscapeTable
		#landscapeData = yhio.Association.getAssociationLandscapeDataFromHDF5File(inputFname=self.inputFname, associationTableName='association', \
		#			landscapeTableName='landscape', min_MAF=self.min_MAF)
		genome_wide_result = associationLandscapeFile.genome_wide_result
		association_peak_ls = self.findPeaks(associationLandscapeFile.locusLandscapeNeighborGraph, bridge_ls=associationLandscapeFile.bridge_ls, \
									data_obj_ls=genome_wide_result.data_obj_ls, \
					ground_score=self.ground_score, min_score=self.min_score)
		attributeDict = {'result_id':genome_wide_result.result_id, \
					'call_method_id':landscapeTable.getAttribute('call_method_id', None),\
					'cnv_method_id':landscapeTable.getAttribute('cnv_method_id', None),\
					'phenotype_method_id':landscapeTable.getAttribute('phenotype_method_id', None),\
					'analysis_method_id':landscapeTable.getAttribute('analysis_method_id', None),\
					'min_MAF':landscapeTable.getAttribute('min_MAF'), \
					'neighbor_distance':landscapeTable.getAttribute('neighbor_distance'), \
					'max_neighbor_distance':landscapeTable.getAttribute('max_neighbor_distance'), \
					'min_score':self.min_score,\
					'ground_score':self.ground_score}
		
		peakPyTable = AssociationPeakTableFile(self.outputFname, openMode='w')
		peakPyTable.addAttributeDict(attributeDict)
		peakPyTable.appendAssociationPeak(association_peak_ls=association_peak_ls)
		#yhio.Association.outputAssociationPeakInHDF5(association_peak_ls=association_peak_ls, filename=self.outputFname, \
		#					tableName='association_peak', closeFile=True,\
		#					attributeDict=attributeDict)
		

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = AssociationLandscape2Peak
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()