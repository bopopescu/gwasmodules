#!/usr/bin/env python
"""

Examples:
	# define landscape for result 5566 and others
	%s  -o  association_min_score_5.0/call_75_analysis_1_min_score_5.0_min_overlap_0.5_result_5371_landscape.h5 
		--result_id 5371 --neighbor_distance 5000 --max_neighbor_distance 20000 --min_MAF 0.1 --tax_id 3702
		--data_dir /Network/Data/250k/db/ --drivername mysql --hostname banyan --dbname stock_250k --db_user yh --db_passwd secret
	
	%s
Description:
	2012.11.20 This program finds the landscape of a genome wide association result and outputs them in HDF5.
	The output contains 2 HDF5 groups. One is "association", the input association result.
		Second is "landscape", the bridge_ls of the landscape.
	This program is upstream of AssociationLandscape2Peak.py
	

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import traceback, numpy
import networkx as nx
from pymodule import ProcessOptions, PassingData, GenomeDB
from pymodule import SNP, yh_matplotlib, yhio
from pymodule import AssociationLandscapeTableFile
from variation.src import Stock_250kDB
from variation.src import AbstractVariationMapper

class DefineAssociationLandscape(AbstractVariationMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVariationMapper.option_default_dict.copy()
	option_default_dict.update({
						('result_id', 1, int): [None, 'j', 1, 'comma or dash-separated list of result ids, i.e. 3431-3438,4041'],\
						('neighbor_distance', 0, int): [5000, '', 1, "within this distance, a locus that increases the association score \
									the fastest is chosen as bridge end. outside this distance, whatever the next point is will be picked."],\
						('max_neighbor_distance', 0, int): [20000, 'm', 1, "beyond this distance, no bridge would be formed."],\
						('tax_id', 0, int): [3702, 'x', 1, 'to get the number of total genes from database, which species.'],\
						})
	#('landscapeLocusIDFname', 1, ): ['', '', 1, 'file that contains one-column, snps.id of the gwas landscape'],\
	#('call_method_id', 0, int):[None, 'l', 1, 'Restrict results based on this call_method. Default is no such restriction.'],\
	#('analysis_method_id_ls', 0, ):['1,7', 'a', 1, 'Restrict results based on these analysis_methods. coma or dash-separated list'],\
	#("phenotype_method_id_ls", 0, ): [None, 'e', 1, 'comma/dash-separated phenotype_method id list, like 1,3-7. Default is all.'],\
	
	def __init__(self,  **keywords):
		"""
		2008-08-19
		"""
		AbstractVariationMapper.__init__(self, **keywords)
	
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
		2012.11.12
			add argument outputFname
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
				bridge_ls.append([current_obj, obj_with_fastest_score_increase, no_of_loci, deltaX, deltaY])
				source_index = current_obj.index
				target_index = obj_with_fastest_score_increase.index
				locusLandscapeNeighborGraph.add_edge(source_index, target_index, \
										weight=obj_with_fastest_score_increase.rate)
				locusLandscapeNeighborGraph[source_index][target_index]['no_of_loci'] = no_of_loci
				deltaX = self.getXDistanceBetweenTwoDataObjs(current_obj, obj_with_fastest_score_increase)
				locusLandscapeNeighborGraph[source_index][target_index]['deltaX'] = deltaX
				start_index = obj_with_fastest_score_increase.index		#new starting point
				
			else:
				start_index = start_index+1	#no bridge is built. at chromosome boundary or beyond the max allowable distance.
		
		sys.stderr.write("%s bridges.\n"%(len(bridge_ls)))
		returnData = PassingData(bridge_ls=bridge_ls, locusLandscapeNeighborGraph=locusLandscapeNeighborGraph)
		return returnData
	
	def drawBridgeChromosomalLengthHist(self, bridge_ls=None):
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
	
	def run(self):
		"""
		2011-3-28
			Read in the association result
			
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_250k = self.db_250k
		
		"""
		#2011-10-12 below commented out because only required for drawing
		genome_db = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.genome_dbname, )
		genome_db.setup(create_tables=False)
		
		oneGenomeData = genome_db.getOneGenomeData(tax_id=self.tax_id, chr_gap=0)
		#oneGenomeData.chr_id2cumu_start
		#cumuSpan2ChrRBDict = oneGenomeData.cumuSpan2ChrRBDict
		"""
		
		
		pd = PassingData(min_MAF=self.min_MAF,\
					data_dir=self.data_dir, \
					need_chr_pos_ls=0,)
		
		rm = Stock_250kDB.ResultsMethod.get(self.result_id)
		
		if rm.cnv_method_id and not db_250k._cnv_id2chr_pos:
			db_250k.cnv_id2chr_pos = rm.cnv_method_id
			pd.db_id2chr_pos = db_250k.cnv_id2chr_pos
		elif rm.call_method_id:
			pd.db_id2chr_pos = db_250k.snp_id2chr_pos
		
		gwr = db_250k.getResultMethodContent(self.result_id, pdata=pd)
		gwr.data_obj_ls.sort(cmp=SNP.cmpDataObjByChrPos)
		#gwr.setResultID(self.result_id)	#already done in db_250k.getResultMethodContent()
		attributeDict = {'result_id':self.result_id, 'min_MAF':self.min_MAF,\
				'call_method_id': rm.call_method_id, 'cnv_method_id': rm.cnv_method_id, \
				'phenotype_method_id':rm.phenotype_method_id, 'analysis_method_id':rm.analysis_method_id,\
				'no_of_accessions': rm.no_of_accessions, 'do_log10_transformation':getattr(gwr,'do_log10_transformation', None)}
		landscapeFile = AssociationLandscapeTableFile(self.outputFname, openMode='w')
		
		writer = gwr.outputInHDF5MatrixFile(tableObject=landscapeFile.associationTable, attributeDict=attributeDict)
		
		landscapeData = self.findLandscape(gwr.data_obj_ls, neighbor_distance=self.neighbor_distance,\
								max_neighbor_distance=self.max_neighbor_distance)
		
		attributeDict.update({'neighbor_distance':self.neighbor_distance, \
							'max_neighbor_distance':self.max_neighbor_distance, \
							})
		
		landscapeFile.associationLandscapeTable.addAttributeDict(attributeDict)
		landscapeFile.appendAssociationLandscapeBridgeList(bridge_ls=landscapeData.bridge_ls)
		#yhio.Association.outputAssociationLandscapeInHDF5(bridge_ls=landscapeData.bridge_ls, writer=writer, closeFile=True, \
		#					tableName='landscape', attributeDict=attributeDict)
		"""
		#2011-4-21 for inspection
		outputFname = '/tmp/result_%s_landscape.png'%(result_id)
		self.drawBridgeLs(bridge_ls=landscapeData.bridge_ls, outputFname=outputFname, oneGenomeData=oneGenomeData)
		"""
		
		self.closeLogF()
		
	
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = DefineAssociationLandscape
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()