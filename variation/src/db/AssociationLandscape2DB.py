#!/usr/bin/env python
"""

Examples:
	# 2012.11.14 
	%s -i  folder/4498_4498_57_37_1_results.tsv  -o  folderAssociationLandscape/Result4498_LandscapeType1_ReducedGWAS.tsv
		--result_id 4498 --phenotype_method_id 37 --analysis_method_id 1
		--results_method_type_id 3 --neighbor_distance 5000
		--max_neighbor_distance 20000 --min_MAF 0.1
		--call_method_id 57 --data_dir /Network/Data/250k/db/ --commit
		--drivername mysql --hostname banyan --dbname stock_250k --db_user yh
	
Description:
	2012.11.12 This program would submit landscape of association results into database.
	
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
from pymodule import AssociationLandscapeTableFile
from variation.src import Stock_250kDB
from variation.src.mapper.AbstractVariationMapper import AbstractVariationMapper

class AssociationLandscape2DB(AbstractVariationMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVariationMapper.option_default_dict.copy()
	option_default_dict.update({
						('result_id', 1, int): [None, '', 1, "ResultsMethod.id of the association result from which the landscape is derived."],\
						('neighbor_distance', 0, int): [5000, '', 1, "within this distance, a locus that increases the association score \
									the fastest is chosen as bridge end. outside this distance, whatever the next point is will be picked."],\
						('max_neighbor_distance', 0, int): [20000, '', 1, "beyond this distance, no bridge would be formed."],\
						('min_MAF', 0, float): [0.1, '', 1, 'minimum Minor Allele Frequency.'],\
						
						})
	
	def __init__(self, **keywords):
		"""
		2012.11.12
		"""
		AbstractVariationMapper.__init__(self, **keywords)
		
	def add2DB(self, db=None, inputFname=None, result=None, association_landscape_type=None, \
			data_dir=None, commit=0,\
			comment=None, db_user=None):
		"""
		2012.11.13
		"""
		session = db.session
		session.begin()
		
		
		#2012.11.13 check if it's in db already
		db_entry = db.checkAssociationLandscape(result_id=result.id, \
										association_landscape_type_id=association_landscape_type.id)
		if db_entry:
			sys.stderr.write("Warning: landscape of (result=%s, landscape-type %s) already in db.\n"%\
							(result.id, association_landscape_type.id))
			sys.exit(3)
		else:
			landscapeFile = AssociationLandscapeTableFile(inputFname, openMode='r')
			no_of_accessions = landscapeFile.getAttribute('no_of_accessions')
			no_of_bridges = landscapeFile.associationLandscapeTable.nrows 
			landscapeFile.close()
			db_entry = db.getAssociationLandscape(result_id=result.id, \
								association_landscape_type_id=association_landscape_type.id,\
								no_of_accessions = no_of_accessions,\
								no_of_bridges = no_of_bridges,\
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
		
		landscapeFile = AssociationLandscapeTableFile(inputFname, openMode='r')
		result_id = landscapeFile.getAttribute('result_id')
		landscapeFile.close()
		
		result = Stock_250kDB.ResultsMethod.get(result_id)
		
		association_landscape_type = self.db_250k.getAssociationLandscapeType(min_MAF=self.min_MAF, \
											neighbor_distance=self.neighbor_distance, \
											max_neighbor_distance=self.max_neighbor_distance)
		
		
		#add the extracted association result into db
		self.add2DB(db=self.db_250k, inputFname=self.inputFname, result=result, association_landscape_type=association_landscape_type, \
				data_dir=self.data_dir, commit=self.commit, db_user=self.db_user)
		
		#2012.6.5
		self.outputLogMessage("submission done.\n")

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = AssociationLandscape2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()