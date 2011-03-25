#!/usr/bin/env python
"""
Examples:

	#2011-3-16 run test on result 4634, list 28, top 500 loci, 20kb max distance from locus to gene 
	TopCNVTest.py -e 4634 -l 28 -f 500 -s genomeDictTaxID3702MaxDist20kb.pickle -m 20000 -y15 -u yh -z banyan
	
Description:
	2011-3-25
		Enrichment test for top CNV associations.
		
		input:
			max distance between a gene and a CNV to declare they are close
			minimum overlap ratio or deleted length
		fetch CNV association test result
			cnv.id	association_score	MAC	MAF
		take the top X number of deletions
		fetch genomeCNVDict
		for each permutation (input: candidate gene set):
			shift base position by a random number (1-genome size)
				translate the new position into chr, start, stop
				if it spans two chromosomes, move it wholly to the left or right chromosome
			find how many candidate genes to which these shifted top CNVs are adjacent (use genomeCNVDict).
		
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, getopt
import warnings, traceback
import Stock_250kDB
from pymodule import getGenomeWideResultFromFile
from pymodule import ProcessOptions, PassingData, GenomeDB, figureOutDelimiter
from GeneListRankTest import GeneListRankTest, SnpsContextWrapper
import rpy, random, numpy
from pymodule.CNV import CNVCompareBySmallOverlapRatio, CNVSegmentBinarySearchTreeKey, CNVCompareByBigOverlapRatio, CNVCompareByOverlapLen
from pymodule.RBTree import RBDict

class TopCNVTest(GeneListRankTest):
	__doc__ = __doc__
	#option_default_dict = GeneListRankTest.option_default_dict.copy()
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['banyan.usc.edu', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
						('genome_dbname', 1, ): ['genome', 'g', 1, 'genome database name', ],\
						('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						("results_id_ls", 1, ): [None, 'e', 1, 'comma/dash-separated results_by_gene id list, like 1,3-7'],\
						("max_distance", 1, int): [20000, 'm', 1, 'maximum distance allowed from the locus to gene'],\
						('min_MAF', 0, float): [0.1, 'n', 1, 'minimum Minor Allele Frequency.'],\
						('min_sample_size', 0, int): [1, 'i', 1, 'minimum size for both candidate and non-candidate sets to do wilcox.test'],\
						("list_type_id", 1, int): [None, 'l', 1, 'Gene list type. must be in table gene_list_type beforehand.'],\
						('results_directory', 0, ):[None, 't', 1, 'The results directory. Default is None. use the one given by db.'],\
						("output_fname", 0, ): [None, 'o', 1, 'To store rank test results into this file as a backup version of db'],\
						("genomeRBDictPickleFname", 0, ): [None, 's', 1, 'The file to contain pickled genomeRBDict. if the file does not exist, it will be loaded from db and pickled into this file'],\
						("results_type", 1, int): [1, 'w', 1, 'which type of results. 1; ResultsMethod, 2: ResultsByGene (forget it), 3: ResultsMethod but counting #distinct genes'],\
						("no_of_permutations", 1, int): [20000, 'N', 1, 'no of permutations to carry out'],\
						("no_of_min_breaks", 1, int): [10, 'B', 1, 'minimum no of times that rank_sum_stat_perm>=rank_sum_stat to break away. if 0, no breaking'],\
						('null_distribution_type_id', 0, int):[1, 'C', 1, 'Type of null distribution. 1=original, 2=permutation, 3=random gene list. check DB table null_distribution_type'],\
						("min_big_overlap", 1, float): [0.3, '', 1, 'minimum of max(overlap1, overlap2) to declare a CNV is close to a gene'],\
						('min_no_of_genes', 1, int):[10, 'G', 1, 'minimum no of genes one candidate gene list should harbor. effective only in MPI version'],\
						('tax_id', 1, int): [3702, 'x', 1, 'to get the number of total genes from database, which species.'],\
						('no_of_top_loci', 1, int): [200, 'f', 1, 'how many number of top snps based on score or -log(pvalue).'],\
						("test_type_id", 1, int): [15, 'y', 1, 'which type of tests. check db table analysis_method. likely be 14,15 etc.'],\
						('min_score', 0, float): [None, 'M', 1, 'alternative way to get top snps. if specified, no_of_top_loci is ignored.'],\
						('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-08-20
		"""
		GeneListRankTest.__init__(self, **keywords)
	
	def translateChrPosIntoCumuPos(self, top_loci, chr_id2cumu_start):
		"""
		2011-3-21
			top_loci has become a list of DataObject of GWR.
		2011-3-16
		
		"""
		top_loci_in_cumu_pos = []
		for top_locus in top_loci:
			start = top_locus.position
			stop = top_locus.stop_position
			chr = str(top_locus.chromosome) 	#chr in chr_id2cumu_start is of type "str"
			cumu_start = chr_id2cumu_start.get(chr) + start -1	#cumu_start is 1-based.
			cumu_stop = chr_id2cumu_start.get(chr) + stop - 1
			top_loci_in_cumu_pos.append([cumu_start, cumu_stop])
		return top_loci_in_cumu_pos
	
	def translateCumuPosIntoChrPos(self, top_loci_in_cumu_pos, cumuSpan2ChrRBDict):
		"""
		2011-3-16
		"""
		top_loci = []
		for span in top_loci_in_cumu_pos:
			cumu_start, cumu_stop = span[:3]
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=0, \
							span_ls=[cumu_start, cumu_stop], \
							min_reciprocal_overlap=0.00000000000001,)
							#2010-8-17 overlapping keys are regarded as separate instances as long as they are identical.
			node = cumuSpan2ChrRBDict.findNode(segmentKey)
			if node is None:
				sys.stderr.write("(%s, %s) not found in cumuSpan2ChrRBDict.\n"%(cumu_start, cumu_stop))
			else:
				chr = str(node.value[0]) 	#chr in chr_id2cumu_start is of type "str"
				start = cumu_start - node.key.span_ls[0] + 1
				stop = cumu_stop - node.key.span_ls[0] + 1
				top_loci.append([chr, start, stop])
		return top_loci
	
	def prepareDataForPermutationRankTest(self, top_loci, genomeRBDict, param_data, report=False):
		"""
		2011-3-16
		"""
		if report:
			sys.stderr.write("Preparing data out of  %s top loci for permutation test ...\n"%\
							(len(top_loci)))
		permData = PassingData(candidate_gene_snp_rank_ls=[],\
							non_candidate_gene_snp_rank_ls=[],\
							captured_candidate_gene_set = set())
		compareIns = CNVCompareByBigOverlapRatio(min_reciprocal_overlap=param_data.min_big_overlap)
		for i in range(len(top_loci)):
			chr, start, stop = top_loci[i][:3]
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=str(chr), \
							span_ls=[start, stop], \
							min_reciprocal_overlap=0.0000001,)
							#min_reciprocal_overlap doesn't matter here.
							# it's decided by compareIns.
			node_ls = []
			genomeRBDict.findNodes(segmentKey, node_ls=node_ls, compareIns=compareIns)
			isNearCandidate = False
			for node in node_ls:
				geneSegKey = node.key
				for oneGeneData in node.value:
					if oneGeneData.ncbi_gene_id in param_data.candidate_gene_set:
						permData.captured_candidate_gene_set.add(oneGeneData.ncbi_gene_id)
						isNearCandidate = True
						break
			if isNearCandidate:
				permData.candidate_gene_snp_rank_ls.append(i+1)
			else:
				permData.non_candidate_gene_snp_rank_ls.append(i+1)
		if report:
			sys.stderr.write("%s loci near %s candidates. Done.\n"%\
							(len(permData.candidate_gene_snp_rank_ls), \
							len(permData.captured_candidate_gene_set)))
		return permData
	
	def applyGWLoopToCumuPos(self, top_loci_in_cumu_pos, cumuSpan2ChrRBDict, shift=None):
		"""
		2011-3-12
		"""
		genome_size = max([max(rbNode.key.span_ls) for rbNode in cumuSpan2ChrRBDict])
		if shift is None:
			shift = random.randint(1, genome_size)
		top_loci_in_cumu_pos_perm = (numpy.array(top_loci_in_cumu_pos)+shift)%genome_size	#modulo to recycle
		
		return self.translateCumuPosIntoChrPos(top_loci_in_cumu_pos_perm, cumuSpan2ChrRBDict)
	
	
	def get_enrichment_pvalue_by_gw_looping(self, candidate_sample_size, top_loci_in_cumu_pos, candidate_gene_set=None, \
							genomeRBDict=None, cumuSpan2ChrRBDict=None, no_of_permutations=20000, \
							no_of_min_breaks=30,\
							param_data=None):
		"""
		2011-3-18
			do the test against permData.captured_candidate_gene_set
		2011-3-12
			get enrichment pvalue by genome-wide looping of SNP positions. a permutation to preserve LD.
		"""
		if self.debug:
			sys.stderr.write("Getting enrichment pvalue by gw-looping ... ")
		i = 0
		no_of_hits = 0
		while i<no_of_permutations:
			permuted_top_loci_in_chr_start_stop = self.applyGWLoopToCumuPos(top_loci_in_cumu_pos, cumuSpan2ChrRBDict)
			
			permData = self.prepareDataForPermutationRankTest(permuted_top_loci_in_chr_start_stop, genomeRBDict, param_data)
			new_candidate_sample_size = len(permData.captured_candidate_gene_set)
			if new_candidate_sample_size>=candidate_sample_size:	#pvalue = Prob(X>=candidate_sample_size)
				no_of_hits += 1
			i+=1
			if no_of_min_breaks>0 and no_of_hits>=no_of_min_breaks:	#if no_of_min_breaks<=0, no smart breaking
				break
		pvalue = no_of_hits/float(i)
		return_data = PassingData(pvalue=pvalue, no_of_tests=i, no_of_tests_passed=no_of_hits)
		if self.debug:
			sys.stderr.write("%s/%s tests in total. Done.\n"%(no_of_hits, i))
		return return_data
	
	def run(self):
		"""
		2010-8-15
			create a RBDict of gene-forms [chr, start, stop] with min_overlap_ratio=1.
				value is a sub-RBDict of the gene structure (UTR, non-UTR-exon, intron)
			
			given any CNV, use RBDict.findNodes() to find all gene-forms.
				WATCH: use an alternative comparison function.
			
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		genome_db = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.genome_dbname, )
		genome_db.setup(create_tables=False)
		
		db_250k = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user, password=self.db_passwd, \
									hostname=self.hostname, database=self.dbname)
		db_250k.setup(create_tables=False)
		
		genome_db.tax_id = self.tax_id
		genome_db.chr_id2size = (self.tax_id,)
		genome_db.chr_id2cumu_start = (self.tax_id, 0)
		cumuSpan2ChrRBDict = genome_db.createCumuSpan2ChrRBDict(self.tax_id, chr_gap=0)
		genomeRBDict = genome_db.dealWithGenomeRBDict(self.genomeRBDictPickleFname, tax_id=self.tax_id, \
									max_distance=self.max_distance, debug=self.debug)
		#genomeRBDict = None
		pd = PassingData(min_MAF=self.min_MAF,\
					min_score=self.min_score, \
					results_directory=self.results_directory, \
					no_of_top_loci=self.no_of_top_loci, \
					starting_rank=0, \
					need_chr_pos_ls=0,\
					need_candidate_association=False,\
					min_big_overlap=self.min_big_overlap,\
					no_of_permutations=self.no_of_permutations,\
					no_of_min_breaks=self.no_of_min_breaks)
		
		for result_id in self.results_id_ls:
			#establish the map from cnv.id from chr_pos
			rm = Stock_250kDB.ResultsMethod.get(result_id)
			if not rm.cnv_method_id:
				sys.stderr.write("ResultsMethod %s doesn't have cnv_method_id. Skip.\n"%(result_id))
				continue
			if not db_250k._cnv_id2chr_pos:
				db_250k.cnv_id2chr_pos = rm.cnv_method_id
			pd.db_id2chr_pos = db_250k.cnv_id2chr_pos
			
			candidate_gene_set = db_250k.dealWithCandidateGeneList(self.list_type_id, return_set=True)	#internal cache
			pd.candidate_gene_set = candidate_gene_set
			
			gwr = db_250k.getResultMethodContent(result_id, pdata=pd)
		
			top_loci = gwr.getTopLoci(no_of_top_loci=self.no_of_top_loci)
			top_loci_in_cumu_pos = self.translateChrPosIntoCumuPos(top_loci, genome_db.chr_id2cumu_start)
			
			param_data = pd
			permData = self.prepareDataForPermutationRankTest(top_loci, genomeRBDict, param_data, report=True)
			
			#m = self.dealWithNoOfSNPsAssociatedWithCandidateGeneList(pd.list_type_id, rm, pd)	#cache is internally going on
			#n = permData.no_of_total_snps - m
			
			candidate_sample_size = len(permData.captured_candidate_gene_set)
			non_candidate_sample_size = len(permData.non_candidate_gene_snp_rank_ls)
			
			return_data = self.get_enrichment_pvalue_by_gw_looping(candidate_sample_size, top_loci_in_cumu_pos, candidate_gene_set, \
							genomeRBDict, cumuSpan2ChrRBDict=cumuSpan2ChrRBDict, \
							no_of_permutations=pd.no_of_permutations, no_of_min_breaks=pd.no_of_min_breaks, param_data=param_data)
			pvalue = return_data.pvalue
			no_of_tests = return_data.no_of_tests
			no_of_tests_passed = return_data.no_of_tests_passed
			sys.stderr.write("%s pvalue: %s.\n"%(result_id, pvalue))
		if self.commit:
			db_250k.session.flush()
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = TopCNVTest
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
