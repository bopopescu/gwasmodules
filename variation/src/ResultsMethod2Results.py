#!/usr/bin/env python
"""
Examples:
	~/script/variation/src/ResultsMethod2Results.py  -l 29 -u yh -p secret -c
	
Description:
	program to pull one results_method from db, pick a certain number of top SNPs and put them into db.
	It checks table Results to see whether it already contains anything from a particular result.
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, getopt, cPickle
import warnings, traceback
from pymodule import PassingData, figureOutDelimiter, getListOutOfStr
import Stock_250kDB
from GeneListRankTest import GeneListRankTest
from DrawSNPRegion import DrawSNPRegion

class ResultsMethod2Results(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("max_rank", 1, int): [1000, 'm', 1, 'maximum rank of the SNP to be put into db. ties-method is random.'],\
							('min_score', 0, float):[None, 'n', 1, 'minimum association score for a SNP.\
								log transformation is automatically determined based on analysis_method.smaller_score_more_significant in db. \
								This argument has an AND relationship with max_rank.'],\
							('call_method_id', 1, int):[None, 'l', 1, 'Restrict results based on this call_method. Default is no such restriction.'],\
							('analysis_method_id_ls', 0, ):['1,7', 'a', 1, 'Restrict results based on these analysis_methods. coma or dash-separated list'],\
							("phenotype_id_ls", 0, ): [None, 'e', 1, 'comma/dash-separated phenotype_method id list, like 1,3-7. Default is all.'],\
							('results_directory', 0, ):[None, 't', 1, 'The results directory. Default is None, then use the one given by db.'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2010-3-8
			add argument min_score, phenotype_id_ls
		2009-2-16
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		self.analysis_method_id_ls = getListOutOfStr(self.analysis_method_id_ls, data_type=int)
		self.phenotype_id_ls = getListOutOfStr(self.phenotype_id_ls, data_type=int)
		
	@classmethod
	def rm2result(cls, session, rm, chr_pos2db_id, max_rank=1000, commit=False, min_rank=1, results_directory=None, \
				min_score=None,update=True,db_id2chr_pos=None):
		"""
		2010-3-8
			add argument min_score to exclude SNPs whose scores are too low. This argument has an AND relationship with max_rank.
				log transformation is automatically determined based on analysis_method.smaller_score_more_significant in db.
		2009-11-2
			split out of run()
		"""
		
		# 2009-5-1 check whether it's already in db.
		db_entries = Stock_250kDB.Results.query.filter_by(results_id=rm.id)
		result_exists = False
		if db_entries.count()==max_rank-min_rank+1:
			if update:
				db_entries.delete()
				sys.stderr.write("%s already in db. Deleting rows.\n"%rm.id)
			else:
				sys.stderr.write("%s already in db. Ignore.\n"%rm.id)

				
			
		
		param_data = PassingData(min_MAC=0, db_id2chr_pos = db_id2chr_pos)
		genome_wide_result = GeneListRankTest.getResultMethodContent(rm, results_directory=results_directory, min_MAF=0., \
																pdata=param_data, min_value_cutoff=min_score)
		
		counter = 0
		no_of_saved = 0
		if genome_wide_result:
			for rank in range(min_rank, max_rank+1):
				if rank>len(genome_wide_result.data_obj_ls):	# rank has gone past the total number of SNPs. break the for loop.
					break
				data_obj = genome_wide_result.get_data_obj_at_given_rank(rank) 
				if data_obj is not None:
					counter += 1
					snps_id = chr_pos2db_id.get((data_obj.chromosome, data_obj.position))
					if data_obj.extra_col_ls:
						result_obj = cPickle.dumps(data_obj.extra_col_ls)
					else:
						result_obj = None
					# 2010-3-8 check if it's in db now.
					db_entries = Stock_250kDB.Results.query.filter_by(results_id=rm.id).filter_by(snps_id=snps_id)
					if db_entries.count()==0:
						Stock_250kDB.Results(snps_id=snps_id, results_id=rm.id, score=data_obj.value, rank=rank, beta=getattr(data_obj, 'beta1', None),\
											maf=data_obj.maf, mac=data_obj.mac, genotype_var_perc=data_obj.genotype_var_perc,\
											correlation=getattr(data_obj, 'correlations', None), odds_ratio=getattr(data_obj, "odds_ratio_est",None),\
											statistic=getattr(data_obj, "statistics",None), object=result_obj)
						no_of_saved += 1
		if commit:
			session.flush()
			#session.commit()
		else:
			session.rollback()
		sys.stderr.write("%s out of %s saved in db.\n"%(no_of_saved, counter))
	
	def run(self):
		"""
		2009-6-10
			set Results.beta = getattr(data_obj, 'beta1', None)
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
									password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		
		session = db.session
		session.begin()
		snp_data = db.getSNPChrPos2IDLs() 
		db_id2chr_pos = snp_data.get('snp_id2chr_pos')
		chr_pos2db_id = snp_data.get('chr_pos_2snp_id')
		
		
		query = Stock_250kDB.ResultsMethod.query.filter_by(call_method_id=self.call_method_id)
		if self.analysis_method_id_ls:
			query = query.filter(Stock_250kDB.ResultsMethod.analysis_method_id.in_(self.analysis_method_id_ls))
		if self.phenotype_id_ls:
			query = query.filter(Stock_250kDB.ResultsMethod.phenotype_method_id.in_(self.phenotype_id_ls))
		for rm in query:
			self.rm2result(session, rm, chr_pos2db_id, max_rank=self.max_rank, commit=self.commit, \
						results_directory=self.results_directory, min_score=self.min_score,db_id2chr_pos = db_id2chr_pos)
		if self.commit:
			session.commit()
	
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = ResultsMethod2Results
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()