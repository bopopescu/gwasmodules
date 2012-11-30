#!/usr/bin/env python
"""

Examples:
	# 2012.11.14 
	%s -i  folder/4498_4498_57_37_1_results.tsv  -o  folderAssociationLandscape/Result4498_LandscapeType1_ReducedGWAS.tsv
		--result_id 4498 --phenotype_method_id 37 --analysis_method_id 1
		--results_method_type_id 3 --neighbor_distance 5000
		--max_neighbor_distance 20000 --min_MAF 0.1
		--call_method_id 57 --data_dir /Network/Data/250k/db/ --commit
		--landscapeLocusIDFname  folderAssociationLandscape/Result4498_LandscapeType1.tsv
		--drivername mysql --hostname banyan --dbname stock_250k --db_user yh
	
Description:
	2012.11.12 This program deduces association loci from peaks within db and then add them into database.
	
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], )
"""
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
"""
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv, stat, getopt, re
import traceback, gc, subprocess
from pymodule import figureOutDelimiter, PassingData
from pymodule import MatrixFile
from pymodule import AbstractDBInteractingJob
from variation.src.common import getOneResultJsonData
from variation.src.db.Results2DB_250k import Results2DB_250k
from variation.src import Stock_250kDB

class AssociationLocus2DB(Results2DB_250k):
	__doc__ = __doc__
	option_default_dict = Results2DB_250k.option_default_dict.copy()
	option_default_dict.update({
						('landscapeLocusIDFname', 1, ): ['', '', 1, 'file that contains one-column, snps.id,'],\
						('result_id', 1, int): [None, '', 1, "ResultsMethod.id of the association result from which the landscape is derived."],\
						('neighbor_distance', 0, int): [5000, '', 1, "within this distance, a locus that increases the association score \
									the fastest is chosen as bridge end. outside this distance, whatever the next point is will be picked."],\
						('max_neighbor_distance', 0, int): [20000, '', 1, "beyond this distance, no bridge would be formed."],\
						('min_MAF', 0, float): [0.1, '', 1, 'minimum Minor Allele Frequency.'],\
						('outputFname', 1, ): [None, 'o', 1, 'this output file stores the reduced gwas (landscape-only).'],\
						
						})
	
	def __init__(self, **keywords):
		"""
		2012.11.12
		"""
		Results2DB_250k.__init__(self, **keywords)
	
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
	
	def add2DB(self, db_250k=None, inputFname=None, result=None, result_landscape_type=None, \
			call_method_id=None, cnv_method_id=None, phenotype_method_id=None, \
			analysis_method_id=None, results_method_type_id=None, \
			no_of_loci = None, \
			data_dir=None, commit=0,\
			comment=None, user=None):
		"""
		2012.11.13
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
		sys.stderr.write("%s association loci added into db.\n"%(counter))
		
		db_250k.session.flush()
		if self.commit:
			db_250k.session.commit()

	def run(self):
		"""
		2012.11.12
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		if not os.path.isfile(self.inputFname):
			sys.stderr.write("Error: file, %s,  is not a file.\n"%(self.inputFname))
			sys.exit(3)
		
		#extract full data from an old association result file (snps.id in landscapeLocusIDFname) and generate a new one
		no_of_loci = self.reduceAssociationResult(inputFname=self.inputFname, landscapeLocusIDFname=self.landscapeLocusIDFname, \
									outputFname=self.outputFname)
		
		result = Stock_250kDB.ResultsMethod.get(self.result_id)
		
		result_peak_type = db_250k.getResultPeakType(min_MAF=self.min_MAF, min_score=self.min_score, \
											neighbor_distance=self.neighbor_distance, \
											max_neighbor_distance=self.max_neighbor_distance)
		#2011-10-12 check to see if ResultPeak contains the peaks from this result already.
		query = Stock_250kDB.ResultPeak.query.filter_by(result_peak_type_id=result_peak_type.id).filter_by(result_id=rm.id)
		if query.first():
			logString = "result_id=%s, result_peak_type_id=%s already exists in ResultPeak. exit.\n"%(result_id, result_peak_type.id)
			self.outputLogMessage(logString)
			sys.stderr.write(logString)
			sys.exit(0)
		
		no_of_peaks = self.findPeaks(db_250k, landscapeData.locusLandscapeNeighborGraph, bridge_ls=landscapeData.bridge_ls, data_obj_ls=gwr.data_obj_ls, \
					ground_score=self.ground_score, min_score=self.min_score, rm=rm, result_peak_type=result_peak_type)
		"""
		#2011-4-21 to check how far things are from each other.
		#self.drawBridgeChromosomalLengthHist(bridge_ls)
		"""
		self.outputLogMessage("Result %s has %s peaks.\n"%(result_id, no_of_peaks))
		
		#add the extracted association result into db
		self.add2DB(db=self.db, inputFname=self.outputFname, result=result, result_landscape_type=result_landscape_type, \
				call_method_id=result.call_method_id, cnv_method_id=result.cnv_method_id, phenotype_method_id=result.phenotype_method_id, \
				analysis_method_id=result.analysis_method_id, results_method_type_id=result.results_method_type_id, \
				no_of_loci=no_of_loci,\
				data_dir=self.data_dir, commit=self.commit, comment=self.comment, user=self.db_user)
		
		#2012.6.5
		self.outputLogMessage("submission done.\n")
		self.closeLogF()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = AssociationLocus2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()