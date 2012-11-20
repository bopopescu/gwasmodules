#!/usr/bin/env python
"""

Examples:
	%s
	
	# set the minimum score to 4
	%s -i /tmp/5566_association.h5  --min_score 4 -o /tmp/5566_association_peak.h5 --result_id 5566
	
Description:
	2012.11.18 this program derives association peak from association landscape (output of DefineAssociationLandscape.py)

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, getopt
import warnings, traceback, numpy
from pymodule import ProcessOptions, PassingData
from pymodule import HDF5MatrixFile
from pymodule import getGenomeWideResultFromHDF5MatrixFile, AbstractMapper
import networkx as nx
from variation.src import Stock_250kDB
from DefineAssociationLandscape import DefineAssociationLandscape
import CGAL	#computational geometry algorithm python binding

class AssociationLandscape2Peak(DefineAssociationLandscape):
	__doc__ = __doc__
	option_default_dict = DefineAssociationLandscape.option_default_dict.copy()
	#get rid of db-related options
	for option_key in AbstractMapper.db_option_dict:
		option_default_dict.pop(option_key)
	
	option_default_dict.update({
						('min_score', 0, float): [4, 'f', 1, 'minimum score to call a peak'],\
						('ground_score', 0, float): [0, 's', 1, 'minimum score possible in this test'],\
						
						})
	
	def __init__(self,  **keywords):
		"""
		2008-08-19
		"""
		DefineAssociationLandscape.__init__(self, **keywords)
	
	def connectDB(self):
		pass
	
	def readLandscape(self, inputFname=None, neighbor_distance=5000,\
					max_neighbor_distance=20000, min_MAF=0.10):
		"""
		2012.11.19
			the inputFname is output of DefineAssociationLandscape.py
		"""
		pdata = PassingData(min_MAF=min_MAF)
		genome_wide_result = getGenomeWideResultFromHDF5MatrixFile(inputFname=inputFname, \
							min_value_cutoff=None, do_log10_transformation=False, pdata=pdata,\
							construct_chr_pos2index=False, construct_data_obj_id2index=False, \
							construct_locus_db_id2index=True,\
							report=True)
		data_obj_ls = genome_wide_result.data_obj_ls
		
		sys.stderr.write("Reading landscape from %s ..."%(inputFname))
		current_obj = None
		bridge_ls = []
		start_index = 0	#starting index in each step of bridge building
		locusLandscapeNeighborGraph = nx.Graph()
		
		reader = HDF5MatrixFile(inputFname, openMode='r')
		landscapeGroupObject = reader.getGroupObject(groupName='landscape')
		start_locus_id_index = landscapeGroupObject.getColIndexGivenColHeader('start_locus_id')
		stop_locus_id_index = landscapeGroupObject.getColIndexGivenColHeader('stop_locus_id')
		no_of_loci_index = landscapeGroupObject.getColIndexGivenColHeader('no_of_loci')
		deltaX_index = landscapeGroupObject.getColIndexGivenColHeader('deltaX')
		
		for row in landscapeGroupObject:
			start_locus_id = row[start_locus_id_index]
			stop_locus_id = row[stop_locus_id_index]
			no_of_loci = row[no_of_loci_index]
			deltaX = row[deltaX_index]
			
			start_obj = genome_wide_result.get_data_obj_by_locus_db_id(start_locus_id)
			stop_obj = genome_wide_result.get_data_obj_by_locus_db_id(stop_locus_id)
			
			bridge_ls.append([start_obj, stop_obj, no_of_loci, deltaX])
			
			source_index = start_obj.index
			#genome_wide_result.get_data_obj_index_by_locus_db_id(start_locus_id)
			target_index = stop_obj.index
			
			locusLandscapeNeighborGraph.add_edge(source_index, target_index, \
										weight=None)
			locusLandscapeNeighborGraph[source_index][target_index]['no_of_loci'] = no_of_loci
			locusLandscapeNeighborGraph[source_index][target_index]['deltaX'] = deltaX
			
		del reader
		
		sys.stderr.write("%s bridges.\n"%(len(bridge_ls)))
		returnData = PassingData(bridge_ls=bridge_ls, locusLandscapeNeighborGraph=locusLandscapeNeighborGraph,\
								genome_wide_result=genome_wide_result)
		return returnData
		
	
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
				ground_score=0, min_score=4, result_id=None):
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
		result_peak_ls = []
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
						result_peak = Stock_250kDB.ResultPeak(result_id = result_id, chromosome=peak_start_data_obj.chromosome,\
								start=peak_start_data_obj.peak_start, stop=intersection_point.x(), no_of_loci=no_of_loci, \
								peak_score=peak_data_obj.value, \
								start_locus_id=peak_start_data_obj.db_id, \
								stop_locus_id=peak_stop_data_obj.db_id, peak_locus_id= peak_data_obj.db_id)
						result_peak_ls.append(result_peak)
						peak_start_data_obj = None	#reset for the next peak
						no_of_peaks += 1
			#subgraph = nx.subgraph(locusLandscapeNeighborGraph, cc)
		sys.stderr.write("%s peaks.\n"%(no_of_peaks))
		return result_peak_ls
	
	def outputPeakInHDF5(self, result_peak_ls=None, filename=None, writer=None, groupName='association_peak', closeFile=True):
		"""
		2012.11.20
		"""
		sys.stderr.write("Dumping association peaks into %s (HDF5 format) ..."%(filename))
		#each number below is counting bytes, not bits
		dtypeList = [('result_id','i8'), ('start_locus_id','i8'), ('stop_locus_id','i8'), \
					('chromosome', HDF5MatrixFile.varLenStrType), ('start','i8'), ('stop', 'i8'), \
					('no_of_loci', 'i8'), ('peak_locus_id', 'i8'), ('peak_score', 'f8')]
		headerList = [row[0] for row in dtypeList]
		dtype = numpy.dtype(dtypeList)
		if writer is None and filename:
			writer = HDF5MatrixFile(filename, openMode='w', dtype=dtype, firstGroupName=groupName)
			writer.writeHeader(headerList)
			groupObject = writer.getGroupObject(groupName=groupName)
		elif writer:
			groupObject = writer.createNewGroup(groupName=groupName, dtype=dtype)
			groupObject.setColIDList(headerList)
		else:
			sys.stderr.write("Error: no writer(%s) or filename(%s) to dump.\n"%(writer, filename))
			sys.exit(3)
		#add neighbor_distance, max_neighbor_distance, min_MAF, min_score, ground_score as attributes
		groupObject.addAttribute(name='min_MAF', value=self.min_MAF, overwrite=True)
		groupObject.addAttribute(name='neighbor_distance', value=self.neighbor_distance, overwrite=True)
		groupObject.addAttribute(name='max_neighbor_distance', value=self.max_neighbor_distance, overwrite=True)
		groupObject.addAttribute(name='min_score', value=self.min_score, overwrite=True)
		groupObject.addAttribute(name='ground_score', value=self.ground_score, overwrite=True)
		cellList = []
		for result_peak in result_peak_ls:
			dataTuple = (result_peak.result_id, result_peak.start_locus_id, result_peak.stop_locus_id, \
						result_peak.chromosome, result_peak.start, result_peak.stop, result_peak.no_of_loci,\
						result_peak.peak_locus_id, result_peak.peak_score)
			cellList.append(dataTuple)
		
		if groupObject is None:
			sys.stderr.write("Error: groupObject (name=%s) is None. could not write.\n"%(groupName))
			sys.exit(3)
		groupObject.writeCellList(cellList)
		if closeFile:
			writer.close()
		sys.stderr.write("%s objects.\n"%(len(cellList)))
		return writer

	
	def run(self):
		"""
		2011-3-28
			Read in the association result
			
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		landscapeData = self.readLandscape(inputFname=self.inputFname, neighbor_distance=self.neighbor_distance,\
								max_neighbor_distance=self.max_neighbor_distance)
		genome_wide_result = landscapeData.genome_wide_result
		result_peak_ls = self.findPeaks(landscapeData.locusLandscapeNeighborGraph, bridge_ls=landscapeData.bridge_ls, \
									data_obj_ls=genome_wide_result.data_obj_ls, \
					ground_score=self.ground_score, min_score=self.min_score, result_id=self.result_id)
		self.outputPeakInHDF5(result_peak_ls=result_peak_ls, filename=self.outputFname, closeFile=True)
		

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = AssociationLandscape2Peak
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()