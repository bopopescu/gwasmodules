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
from pymodule import AssociationLocusTableFile, utils
from variation.src import Stock_250kDB
from AssociationPeak2DB import AssociationPeak2DB
 
class AssociationLocus2DB(AssociationPeak2DB):
	__doc__ = __doc__
	option_default_dict = AssociationPeak2DB.option_default_dict.copy()
	option_default_dict.update({
			('min_connectivity', 0, float): [None, '', 1, "minimum connectivity of any peak graph connected component, not used."],\
			('min_overlap_ratio', 0, float): [0.5, '', 1, 'minimum overlap ratio between two peaks for them to merge. overlap length/total'],\
			})
	
	def __init__(self, **keywords):
		"""
		2012.11.12
		"""
		AssociationPeak2DB.__init__(self, **keywords)
	
	def add2DB(self, db=None, inputFname=None, association_locus_type=None, \
			data_dir=None, commit=0,\
			comment=None, db_user=None):
		"""
		2012.12.21
		"""
		
		session = db.session
		session.begin()
		
		associationLocusTableFile = AssociationLocusTableFile(inputFname, openMode='r', constructLocusRBDict=False)	#don't need it
		call_method_id_ls = associationLocusTableFile.getAttribute('call_method_id_ls')
		if len(call_method_id_ls)==1:	#only assign singular call_method_id when there is only one call method
			call_method_id = call_method_id_ls[0]
		else:
			call_method_id = None
		analysis_method_id_ls = associationLocusTableFile.getAttribute('analysis_method_id_ls')
		if len(analysis_method_id_ls)==1:
			analysis_method_id = analysis_method_id_ls[0]
		else:
			analysis_method_id = None
		phenotype_method_id_ls_str = utils.getSuccinctStrOutOfList(associationLocusTableFile.getAttribute('phenotype_method_id_ls'))
		total_no_of_results = associationLocusTableFile.getAttribute('total_no_of_results')
		no_of_loci = associationLocusTableFile.getTableObject().nrows
		
		associationLocusTableFile.close()
		
		call_method_id_ls_str = utils.getSuccinctStrOutOfList(call_method_id_ls)
		analysis_method_id_ls_str = utils.getSuccinctStrOutOfList(analysis_method_id_ls)
		#2012.11.13 check if it's in db already
		db_entry = db.checkGenomeWideAssociationLocus(association_locus_type_id=association_locus_type.id,\
								call_method_id=call_method_id, analysis_method_id=analysis_method_id,\
								call_method_id_ls=call_method_id_ls_str,\
								analysis_method_id_ls=analysis_method_id_ls_str,\
								phenotype_method_id_ls=phenotype_method_id_ls_str)
		if db_entry:
			sys.stderr.write("Warning: genome-wide association-locus of (association_locus_type_id %s, call %s, analysis %s, phenotype %s) already in db.\n"%\
							(association_locus_type.id, call_method_id_ls_str, analysis_method_id_ls_str, phenotype_method_id_ls_str))
			sys.exit(3)
		else:
			
			db_entry = db.getGenomeWideAssociationLocus(association_locus_type_id=association_locus_type.id,\
								no_of_loci=no_of_loci, total_no_of_results=total_no_of_results, \
								call_method_id=call_method_id, analysis_method_id=analysis_method_id,\
								call_method_id_ls=call_method_id_ls_str,\
								analysis_method_id_ls=analysis_method_id_ls_str,\
								phenotype_method_id_ls=phenotype_method_id_ls_str,
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
				
		association_landscape_type = self.db_250k.getAssociationLandscapeType(min_MAF=self.min_MAF, \
											neighbor_distance=self.neighbor_distance, \
											max_neighbor_distance=self.max_neighbor_distance)
		
		association_peak_type = self.db_250k.getAssociationPeakType(association_landscape_type_id=association_landscape_type.id, \
															min_score=self.min_score)
		
		association_locus_type = self.db_250k.getAssociationLocusType(association_peak_type_id=association_peak_type.id, \
									min_overlap_ratio=self.min_overlap_ratio, \
									min_connectivity=self.min_connectivity)
			
		"""
		#2011-4-21 to check how far things are from each other.
		#self.drawBridgeChromosomalLengthHist(bridge_ls)
		"""
		self.add2DB(db=self.db_250k, inputFname=self.inputFname, association_locus_type=association_locus_type, \
			data_dir=self.data_dir, commit=self.commit,\
			comment=None, db_user=self.db_user)
		
		#2012.6.5
		self.outputLogMessage("submission done.\n")
		self.closeLogF()


if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = AssociationLocus2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()