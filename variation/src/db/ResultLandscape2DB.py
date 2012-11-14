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
from variation.src.Stock_250kDB import Results, ResultsMethod, PhenotypeMethod, CallMethod, \
	ResultsMethodType, AnalysisMethod, ResultsMethodJson
from variation.src import Stock_250kDB
from pymodule import figureOutDelimiter, PassingData
from pymodule.MatrixFile import MatrixFile

from variation.src.common import getOneResultJsonData
from pymodule.pegasus.mapper.AbstractDBInteractingJob import AbstractDBInteractingJob
from variation.src.Results2DB_250k import Results2DB_250k


class ResultLandscape2DB(Results2DB_250k):
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
		
		from pymodule import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def add2DB(self, db=None, inputFname=None, result=None, result_landscape_type=None, \
			call_method_id=None, cnv_method_id=None, phenotype_method_id=None, \
			analysis_method_id=None, results_method_type_id=None, \
			no_of_loci = None, \
			data_dir=None, commit=0,\
			comment=None, user=None):
		"""
		2012.11.13
		"""
		session = db.session
		session.begin()
		
		
		rmt = Stock_250kDB.ResultsMethodType.get(results_method_type_id)
		
		if not rmt:
			sys.stderr.write("No results method type available for results_method_type_id=%s.\n"%results_method_type_id)
			sys.exit(3)
		
		pm = Stock_250kDB.PhenotypeMethod.query.get(phenotype_method_id)
		if not pm:
			sys.stderr.write("No phenotype method available for phenotype_method_id=%s.\n"%phenotype_method_id)
			sys.exit(3)
		if call_method_id:	#2012.6.6
			cm = Stock_250kDB.CallMethod.query.get(call_method_id)
		else:
			cm = None
		if cnv_method_id:	#2012.6.6 this could None
			cnv_m = Stock_250kDB.CNVMethod.get(cnv_method_id)
		else:
			cnv_m = None
		if not cm and not cnv_m:
			sys.stderr.write("Neither call method (%s) nor cnv method (%s) available.\n"%(call_method_id, cnv_method_id))
			sys.exit(3)
		
		am = AnalysisMethod.query.get(analysis_method_id)
		if not am:
			sys.stderr.write("No analysis method available for analysis_method_id=%s.\n"%analysis_method_id)
			sys.exit(3)
		
		#2012.11.13 check if it's in db already
		query = Stock_250kDB.ResultLandscape.query.filter_by(result_id=result.id).\
			filter_by(result_landscape_type_id=result_landscape_type.id)
		if query.count()>0:
			db_entry = query.first()
			sys.stderr.write("There is already an entry in result_landscape (id=%s) with same (result_id, result_landscape_type_id)=(%s, %s).\n"\
							%(db_entry.id, db_entry.result_id, result_landscape_type.id))
			sys.exit(3)
		#2012.11.13 check if it's in db already
		query = Stock_250kDB.ResultLandscape.query.filter_by(phenotype_method_id=pm.id).\
			filter_by(analysis_method_id=am.id).filter_by(results_method_type_id=rmt.id).\
			filter_by(result_landscape_type_id=result_landscape_type.id)
		if cm:	#2011-2-22
			query = query.filter_by(call_method_id=cm.id)
		elif cnv_m:	#2011-2-22
			query = query.filter_by(cnv_method_id=cnv_m.id)
		if query.count()>0:
			db_entry = query.first()
			sys.stderr.write("There is already an entry in result_landscape (id=%s) with same \
	(call_method_id, phenotype_method_id, analysis_method_id, results_method_type_id, result_landscape_type_id)=(%s, %s, %s, %s, %s).\n"\
							%(db_entry.id, call_method_id, phenotype_method_id, analysis_method_id, results_method_type_id, \
							result_landscape_type.id))
			sys.exit(3)
		
		short_name = '%s_%s_%s_%s_%s_%s'%(am.short_name, pm.short_name, getattr(cm, 'id',0), getattr(cnv_m, 'id',0),\
										results_method_type_id, result_landscape_type.id)
		
		db_entry = Stock_250kDB.ResultLandscape(short_name=short_name, \
											locus_type_id=cm.locus_type_id, no_of_loci=no_of_loci,\
											comment=comment, created_by=user)	#2012.3.9
		db_entry.phenotype_method = pm
		db_entry.call_method = cm
		db_entry.analysis_method = am
		db_entry.cnv_method = cnv_m
		db_entry.result = result
		db_entry.result_landscape_type = result_landscape_type
		db_entry.results_method_type = rmt
		session.add(db_entry)
		
		session.flush()	#not necessary as no immediate query on the new results after this and commit() would execute this.
		if commit:
			db_entry.filename = os.path.join(db.data_dir, db_entry.constructRelativePath(data_dir=data_dir))
			localAbsPath = os.path.join(data_dir, db_entry.constructRelativePath(data_dir=data_dir))
			if db_entry.analysis_method_id==13:
				return_value = self.copyResultsFile(db, inputFname, db_entry, user, localAbsPath)
			else:
				return_value = self.submit_results(db, inputFname, db_entry, user, localAbsPath)
			if return_value:
				session.add(db_entry)
			else:	#bad thing happened when getting data out of the file. don't save this results_method.
				session.delete(db_entry)
			session.flush()
			session.commit()
			self.reset_marker_pos2snp_id()
		else:	#default is also rollback(). to demonstrate good programming
			session.rollback()
			self.reset_marker_pos2snp_id()
	
	def reduceAssociationResult(self, inputFname=None, landscapeLocusIDFname=None, outputFname=None):
		"""
		2012.11.13
		"""
		sys.stderr.write("Getting a set of landscape-locus from %s ..."%(landscapeLocusIDFname))
		reader = MatrixFile(inputFname=landscapeLocusIDFname)
		locus_id_set = set()
		for row in reader:
			locus_id = row[0].strip()
			locus_id_set.add(locus_id)
		del reader
		sys.stderr.write("%s loci.\n"%(len(locus_id_set)))
		
		sys.stderr.write("Reducing association result %s into a landscape (%s) using %s ... "%(inputFname, outputFname, landscapeLocusIDFname))
		reader = MatrixFile(inputFname=inputFname)
		reader.constructColName2IndexFromHeader()
		
		outf = open(outputFname, 'w')
		writer = csv.writer(outf, delimiter='\t')
		writer.writerow(reader.header)
		
		counter = 0
		real_counter = 0
		for row in reader:
			counter += 1
			locus_id = row[0]
			if locus_id in locus_id_set:
				writer.writerow(row)
				real_counter += 1
		outf.close()
		del reader, writer
		sys.stderr.write("%s out of %s loci chosen.\n"%(real_counter, counter))
		return real_counter
	
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
		
		result_landscape_type = self.db.getResultLandscapeType(min_MAF=self.min_MAF, \
											neighbor_distance=self.neighbor_distance, \
											max_neighbor_distance=self.max_neighbor_distance)
		
		
		#add the extracted association result into db
		self.add2DB(db=self.db, inputFname=self.outputFname, result=result, result_landscape_type=result_landscape_type, \
				call_method_id=result.call_method_id, cnv_method_id=result.cnv_method_id, phenotype_method_id=result.phenotype_method_id, \
				analysis_method_id=result.analysis_method_id, results_method_type_id=result.results_method_type_id, \
				no_of_loci=no_of_loci,\
				data_dir=self.data_dir, commit=self.commit, comment=self.comment, user=self.db_user)
		
		#2012.6.5
		if self.logFilename:
			logFile = open(self.logFilename, 'w')
			logFile.write("submission done.\n")
			logFile.close()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = ResultLandscape2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()