#!/usr/bin/env python
"""

Examples:
	PlotTopSNPTestRm.py -e 1-5 -l 1  -o /tmp/hist_of_results_by_gene_candidate_score_rank
	
Description:
	2008-09-28
	program to check histogram of ranks of candidate genes vs those of non-candidate genes in results_by_gene.
	it also draws histogram of scores of results_by_gene.
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
#	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
#	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib as mpl; mpl.use("Agg")
import time, csv, cPickle
import warnings, traceback
from pymodule import PassingData, figureOutDelimiter, getColName2IndexFromHeader, getListOutOfStr
import Stock_250kDB
from Stock_250kDB import ResultsByGene, ResultsMethod
from sets import Set
from GeneListRankTest import GeneListRankTest, SnpsContextWrapper
#from sqlalchemy.orm import join
import pylab
import StringIO

class PlotTopSNPTestRm(GeneListRankTest):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("results_method_id_ls", 1, ): [None, '', 1, 'comma/dash-separated results_method id list, like 1,3-7'],\
							("min_distance", 1, int): [20000, 'm', 1, 'minimum distance allowed from the SNP to gene'],\
							("get_closest", 0, int): [0, 'g', 0, 'only get genes closest to the SNP within that distance'],\
							("no_of_snps", 0, int): [200, '', 1, 'number of snps in each rank-window'],\
							('min_MAF', 1, float): [0.1, 'n', 1, 'minimum Minor Allele Frequency. deprecated.'],\
							('test_type', 0, int): [1, 'y', 1, 'minimum size for candidate gene sets to draw histogram'],\
							("results_type", 1, int): [1, 'w', 1, 'which type of results. 1; ResultsMethod, 2: ResultsByGene'],\
							("list_type_id_ls", 1, ): [None, 'l', 1, 'comma/dash-separated Gene list type. must be in table gene_list_type beforehand.'],\
							('results_directory', 0, ):[None, 't', 1, 'The results directory. Default is None. use the one given by db.'],\
							("output_dir", 1, ): [None, 'o', 1, 'directory to store output'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	
	def __init__(self,  **keywords):
		"""
		2008-10-20
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		self.results_method_id_ls = getListOutOfStr(self.results_method_id_ls, data_type=int)
		self.list_type_id_ls = getListOutOfStr(self.list_type_id_ls, data_type=int)
		
		if self.output_dir and not os.path.isdir(self.output_dir):
			os.makedirs(self.output_dir)
	
	def plot(self, ResultsClass, results_id, list_type_id, no_of_snps, min_distance, get_closest, min_MAF, test_type,\
			output_fname_prefix=None, commit=0):
		"""
		"""
		sys.stderr.write("Plotting for result %s on list %s ..."%(results_id, list_type_id))
		pylab.clf()
		rows = ResultsClass.query.filter_by(results_id=results_id).filter_by(list_type_id=list_type_id).\
			filter_by(no_of_top_snps=no_of_snps).filter_by(min_distance=min_distance).filter_by(get_closest=get_closest).\
			filter(ResultsClass.min_MAF>=min_MAF-0.0001).filter(ResultsClass.min_MAF<=min_MAF+0.0001).filter_by(test_type=test_type)
		avg_rank_ls = []
		pvalue_ls = []
		zero_pvalue_avg_rank_ls = []
		title = None
		for row in rows:
			avg_rank = (row.no_of_top_snps+row.starting_rank-1 + row.starting_rank)/2.
			if row.pvalue>0:
				pvalue_ls.append(-math.log10(row.pvalue))
				avg_rank_ls.append(avg_rank)
			else:
				zero_pvalue_avg_rank_ls.append(avg_rank)
			
			if title is None:
				title = '%s on %s (rm.id=%s) with candidate gene list=%s'%(row.result.analysis_method.short_name, row.result.phenotype_method.short_name, row.results_id, row.list_type.short_name)
		
		png_data = None
		svg_data = None
		if avg_rank_ls and pvalue_ls:
			pylab.plot(avg_rank_ls, pvalue_ls, '.', alpha=0.3, markersize=4)
			pylab.title(title)
			pylab.xlabel('average rank')
			pylab.ylabel('-log10(pvalue)')
			if zero_pvalue_avg_rank_ls:
				no_of_zero_pvalues = len(zero_pvalue_avg_rank_ls)
				pylab.plot(zero_pvalue_avg_rank_ls, [max(pvalue_ls)]*no_of_zero_pvalues, '.', c='r', alpha=0.3, markersize=4)
		else:
			return png_data, svg_data
		
		if commit:
			png_data = StringIO.StringIO()
			svg_data = StringIO.StringIO()
			pylab.savefig(png_data, format='png', dpi=300)
			pylab.savefig(svg_data, format='svg', dpi=300)
		elif output_fname_prefix:
			pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
			pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		sys.stderr.write("Done.\n")
		return png_data, svg_data
	
	def run(self):
		"""
		2008-10-20
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		#session.begin()
		
		if self.results_type==1:
			ResultsClass = Stock_250kDB.CandidateGeneTopSNPTestRM
		elif self.results_type==2:
			ResultsClass = Stock_250kDB.CandidateGeneTopSNPTest
		else:
			sys.stderr.write("Invalid results type : %s.\n"%self.results_type)
			return None
		
		for results_id in self.results_method_id_ls:
			for list_type_id in self.list_type_id_ls:
				output_fname_prefix = os.path.join(self.output_dir, 'rank_vs_pvalue_%s_%s'%(results_id, list_type_id))
				png_data, svg_data = self.plot(ResultsClass, results_id, list_type_id, self.no_of_snps, self.min_distance, \
											self.get_closest, self.min_MAF, self.test_type,\
											output_fname_prefix, commit=self.commit)
		

if __name__ == '__main__':
	from matplotlib import rcParams
	rcParams['font.size'] = 6
	rcParams['legend.fontsize'] = 6
	#rcParams['text.fontsize'] = 6	#deprecated. use font.size instead
	rcParams['axes.labelsize'] = 4
	rcParams['axes.titlesize'] = 6
	rcParams['xtick.labelsize'] = 6
	rcParams['ytick.labelsize'] = 6

	from pymodule import ProcessOptions
	main_class = PlotTopSNPTestRm
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()