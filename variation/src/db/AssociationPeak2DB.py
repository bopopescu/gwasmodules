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
	2012.12.21 This program submits association peak file from one association result into database (table GenomeWideAssociationPeak).
		not table AssociationTeak.
		The corresponding association landscape has to be in the database already.
	
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], )
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv, stat, getopt, re
import traceback, gc, subprocess
from pymodule import figureOutDelimiter, PassingData
from pymodule import AssociationPeakTableFile
from variation.src import Stock_250kDB
from AssociationLandscape2DB import AssociationLandscape2DB

class AssociationPeak2DB(AssociationLandscape2DB):
	__doc__ = __doc__
	option_default_dict = AssociationLandscape2DB.option_default_dict.copy()
	option_default_dict.update({
						('min_score', 0, float): [3, '', 1, 'minimum score to cut an association landscape into peaks, -log(pvalue).'],\
						})
	
	def __init__(self, **keywords):
		"""
		2012.11.12
		"""
		AssociationLandscape2DB.__init__(self, **keywords)
	
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
	
	def add2DB(self, db=None, inputFname=None, result=None, association_landscape=None, \
			association_peak_type=None, \
			data_dir=None, commit=0,\
			comment=None, db_user=None):
		"""
		2012.12.21
		"""
		
		session = db.session
		session.begin()
		
		
		#2012.11.13 check if it's in db already
		db_entry = db.checkGenomeWideAssociationPeak(result_id=result.id, \
										association_landscape_id=association_landscape.id, \
										association_peak_type_id=association_peak_type.id)
		if db_entry:
			sys.stderr.write("Warning: genome-wide association-peak of (result=%s, landscape %s, peak-type %s) already in db.\n"%\
							(result.id, association_landscape.id, association_peak_type.id))
			sys.exit(3)
		else:
			peakFile = AssociationPeakTableFile(inputFname, openMode='r')
			no_of_peaks = peakFile.nrows
			peakFile.close()
			db_entry = db.getGenomeWideAssociationPeak(result_id=result.id, \
								association_landscape_id=association_landscape.id,\
								association_peak_type_id=association_peak_type.id,\
								no_of_peaks = no_of_peaks,\
								original_path=os.path.abspath(inputFname), comment=comment, created_by=db_user)
		
		
		if commit:
			inputFileBasename = os.path.basename(inputFname)
			#moveFileIntoDBAffiliatedStorage() will also set db_entry.path
			exitCode = db.moveFileIntoDBAffiliatedStorage(db_entry=db_entry, filename=inputFileBasename, \
									inputDir=os.path.split(inputFname)[0], \
									outputDir=data_dir,\
									relativeOutputDir=None, shellCommand='cp -rL', \
									srcFilenameLs=self.srcFilenameLs, dstFilenameLs=self.dstFilenameLs,\
									constructRelativePathFunction=db_entry.constructRelativePath, data_dir=data_dir)
			
			if exitCode!=0:
				sys.stderr.write("Error: moveFileIntoDBAffiliatedStorage() exits with %s code.\n"%(exitCode))
				session.rollback()
				self.cleanUpAndExitOnFailure(exitCode=exitCode)
			
			session.flush()
			session.commit()
		else:	#default is also rollback(). to demonstrate good programming
			session.rollback()
		"""
								association_peak = Stock_250kDB.AssociationPeak(result_id = rm.id, chromosome=peak_start_data_obj.chromosome,\
								start=peak_start_data_obj.peak_start, stop=intersection_point.x(), no_of_loci=no_of_loci, \
								peak_score=peak_data_obj.value)
						if rm.cnv_method_id:
							# 2011-4-21 assuming locus might not be in db.
							association_peak.start_locus = self.getLocusBasedOnDataObj(db_250k, peak_start_data_obj)
							association_peak.stop_locus = self.getLocusBasedOnDataObj(db_250k, peak_stop_data_obj)
							association_peak.peak_locus = self.getLocusBasedOnDataObj(db_250k, peak_data_obj)
						elif rm.call_method_id:
							# 2011-4-21 assuming all relevant loci are in db.
							association_peak.start_locus_id = peak_start_data_obj.db_id
							association_peak.stop_locus_id = peak_stop_data_obj.db_id
							association_peak.peak_locus_id = peak_data_obj.db_id
						
						association_peak.association_peak_type = association_peak_type
						
						db_250k.session.add(association_peak)
						peak_start_data_obj = None	#reset for the next peak
						no_of_peaks += 1
			db_250k.session.flush()
			#subgraph = nx.subgraph(locusLandscapeNeighborGraph, cc)
		sys.stderr.write("%s peaks.\n"%(no_of_peaks))
		"""
	
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
		
		peakFile = AssociationPeakTableFile(self.inputFname, openMode='r')
		result_id = peakFile.getAttribute('result_id')
		peakFile.close()
		result = Stock_250kDB.ResultsMethod.get(result_id)
		
		association_landscape_type = self.db_250k.getAssociationLandscapeType(min_MAF=self.min_MAF, \
											neighbor_distance=self.neighbor_distance, \
											max_neighbor_distance=self.max_neighbor_distance)
		
		association_peak_type = self.db_250k.getAssociationPeakType(association_landscape_type_id=association_landscape_type.id, \
															min_score=self.min_score)
	
		association_landscape = self.db_250k.checkAssociationLandscape(result_id=result.id, \
								association_landscape_type_id=association_landscape_type.id)
		if not association_landscape:
			sys.stderr.write("Error, landscape for result_id=%s, association_landscape_type_id=%s not in db.\n"%\
							(result.id, association_landscape_type.id))
			sys.exit(2)
		
		"""
		#2011-4-21 to check how far things are from each other.
		#self.drawBridgeChromosomalLengthHist(bridge_ls)
		"""
		self.add2DB(db=self.db_250k, inputFname=self.inputFname, result=result, association_landscape=association_landscape, \
			association_peak_type=association_peak_type, \
			data_dir=self.data_dir, commit=self.commit,\
			comment=None, db_user=self.db_user)
		
		#2012.6.5
		self.outputLogMessage("submission done.\n")
		self.closeLogF()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = AssociationPeak2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()