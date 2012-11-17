#!/usr/bin/env python
"""
Examples:

	#2011-10-15 run test on (cnv_20_LD_KW, LDV, SD, SDV), peak type 1(min_score=4),  list 28, 20kb max distance from locus to gene 
	%s.py -e 4634,4635,4636,4637 -t 1 -l 28 -s genomeDictTaxID3702MaxDist20kb.pickle -m 20000 -y15 -u yh -z banyan
	
	#2011-10-15 run test on call_32_LD_KW, LDV, SD, SDV;  peak type 2(min_score=5), 
	%s  -e 3395,3396,3397,3398 -t 2 -l 28 -s genomeDictTaxID3702MaxDist20kb.pickle -m 20000 -y15 -u yh -z banyan

	#2011-10-15 run test on call_80_LD_KW, LDV, SD, SDV;  peak type 2(min_score=5),
	%s -e 4885-4891 -t 2 -l 28 -s genomeDictTaxID3702MaxDist20kb.pickle -m 20000 -y15 -u yh -z banyan
	
	#2011-10-15 run test on call_80_LD_KW, LDV, SD, SDV;  peak type 2(min_score=4),
	%s -e 4885-4891 -t 1 -l 28 -s genomeDictTaxID3702MaxDist20kb.pickle -m 20000 -y15 -u yh -z banyan
	
	
Description:
	2011-10-15
		Enrichment test for ResultPeak, based off an old version of TopCNVTest.py
		
		input:
			max distance between a gene and a peak to declare they are close
			minimum overlap ratio or deleted length
		fetch result peak
			cnv.id	association_score	MAC	MAF
		fetch genomeCNVDict
		for each permutation (input: candidate gene set):
			shift base position by a random number (1-genome size)
				translate the new position into chr, start, stop
				if it spans two chromosomes, move it wholly to the left or right chromosome
			find how many candidate genes to which these shifted top CNVs are adjacent (use genomeCNVDict).
		
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

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
import random, numpy
from pymodule.CNV import CNVCompareBySmallOverlapRatio, CNVSegmentBinarySearchTreeKey, CNVCompareByBigOverlapRatio, CNVCompareByOverlapLen
from pymodule.RBTree import RBDict
from pymodule.SNP import DataObject

class ResultPeakEnrichmentTest(GeneListRankTest):
	__doc__ = __doc__
	#option_default_dict = GeneListRankTest.option_default_dict.copy()
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['banyan.usc.edu', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
						('genome_dbname', 1, ): ['genome', 'g', 1, 'genome database name', ],\
						('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						("result_id_ls", 1, ): [None, 'e', 1, 'comma/dash-separated results_by_gene id list, like 1,3-7'],\
						("result_peak_type_id", 1, int): [None, 't', 1, 'id in ResultPeakType, used to restrict data from ResultPeak'],\
						("max_distance", 1, int): [20000, 'm', 1, 'maximum distance allowed from the locus to gene'],\
						('min_MAF', 0, float): [0.1, 'n', 1, 'minimum Minor Allele Frequency.'],\
						('min_sample_size', 0, int): [1, 'i', 1, 'minimum size for both candidate and non-candidate sets to do wilcox.test'],\
						("list_type_id", 1, int): [None, 'l', 1, 'Gene list type. must be in table gene_list_type beforehand.'],\
						("output_fname", 0, ): [None, 'o', 1, 'To store rank test results into this file as a backup version of db'],\
						("genomeRBDictPickleFname", 0, ): [None, 's', 1, 'The file to contain pickled genomeRBDict. if the file does not exist, it will be loaded from db and pickled into this file'],\
						("no_of_permutations", 1, int): [20000, 'N', 1, 'no of permutations to carry out'],\
						("no_of_min_breaks", 1, int): [10, 'B', 1, 'minimum no of times that rank_sum_stat_perm>=rank_sum_stat to break away. if 0, no breaking'],\
						('null_distribution_type_id', 0, int):[1, 'C', 1, 'Type of null distribution. 1=original, 2=permutation, 3=random gene list. check DB table null_distribution_type'],\
						("min_big_overlap", 1, float): [0.3, 'a', 1, 'minimum of max(overlap1, overlap2) to declare a CNV is close to a gene'],\
						('min_no_of_genes', 1, int):[10, 'G', 1, 'minimum no of genes one candidate gene list should harbor. effective only in MPI version'],\
						('tax_id', 1, int): [3702, 'x', 1, 'to get the number of total genes from database, which species.'],\
						("test_type_id", 1, int): [15, 'y', 1, 'which type of tests. check db table analysis_method. likely be 14,15 etc.'],\
						('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-08-20
		"""
		GeneListRankTest.__init__(self, **keywords)
	
	def translateChrPosDataObjectIntoCumuPos(self, top_loci, chr_id2cumu_start):
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
			cumu_start = chr_id2cumu_start.get(chr) + start	#cumu_start is 0-based.
			cumu_stop = chr_id2cumu_start.get(chr) + stop
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
	
	def prepareDataForPermutationRankTest(self, top_loci_in_chr_pos, genomeRBDict, param_data, report=False):
		"""
		2011-3-16
		"""
		if report:
			sys.stderr.write("Preparing data out of  %s top loci for permutation test ...\n"%\
							(len(top_loci_in_chr_pos)))
		permData = PassingData(candidate_gene_snp_rank_ls=[],\
							non_candidate_gene_snp_rank_ls=[],\
							captured_candidate_gene_set = set())
		compareIns = CNVCompareByBigOverlapRatio(min_reciprocal_overlap=param_data.min_big_overlap)
		for i in range(len(top_loci_in_chr_pos)):
			chr, start, stop = top_loci_in_chr_pos[i][:3]
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
					if oneGeneData.gene_id in param_data.candidate_gene_set:
						permData.captured_candidate_gene_set.add(oneGeneData.gene_id)
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
	
	def getResultPeak(self, result_id, result_peak_type_id=1, param_data=None):
		"""
		2011-10-15
			get a list of peaks in the form of pymodule.SNP.DataObject
		"""
		sys.stderr.write("Getting result peak for result %s, peak type %s ..."%(result_id, result_peak_type_id))
		query = Stock_250kDB.ResultPeak.query.filter_by(result_id=result_id).filter_by(result_peak_type_id=result_peak_type_id)
		counter = 0
		real_counter = 0
		top_loci = []
		for row in query:
			counter += 1
			data_obj = DataObject(db_id=row.id, chromosome=row.chromosome, position=row.start, stop_position=row.stop, \
								value =row.peak_score)
			top_loci.append(data_obj)
		sys.stderr.write(" %s peaks.\n"%(len(top_loci)))
		return top_loci

	
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
		
		oneGenomeData = genome_db.getOneGenomeData(tax_id=self.tax_id, chr_gap=0)
		cumuSpan2ChrRBDict = oneGenomeData.cumuSpan2ChrRBDict
		genomeRBDict = genome_db.dealWithGenomeRBDict(self.genomeRBDictPickleFname, tax_id=self.tax_id, \
									max_distance=self.max_distance, debug=self.debug)
		#genomeRBDict = None
		pd = PassingData(min_MAF=self.min_MAF,\
					starting_rank=0, \
					need_chr_pos_ls=0,\
					need_candidate_association=False,\
					min_big_overlap=self.min_big_overlap,\
					no_of_permutations=self.no_of_permutations,\
					no_of_min_breaks=self.no_of_min_breaks)
		
		for result_id in self.result_id_ls:
			candidate_gene_set = db_250k.dealWithCandidateGeneList(self.list_type_id, return_set=True)	#internal cache
			pd.candidate_gene_set = candidate_gene_set
			
			#gwr = db_250k.getResultMethodContent(result_id, pdata=pd)
		
			#top_loci = gwr.getTopLoci(no_of_top_loci=self.no_of_top_loci)
			top_loci = self.getResultPeak(result_id, self.result_peak_type_id, pd)
			
			top_loci_in_cumu_pos = self.translateChrPosDataObjectIntoCumuPos(top_loci, oneGenomeData.chr_id2cumu_start)
			top_loci_in_chr_pos = self.translateCumuPosIntoChrPos(top_loci_in_cumu_pos, cumuSpan2ChrRBDict)
			permData = self.prepareDataForPermutationRankTest(top_loci_in_chr_pos, genomeRBDict, pd, report=True)
			
			#m = self.dealWithNoOfSNPsAssociatedWithCandidateGeneList(pd.list_type_id, rm, pd)	#cache is internally going on
			#n = permData.no_of_total_snps - m
			
			candidate_sample_size = len(permData.captured_candidate_gene_set)
			non_candidate_sample_size = len(permData.non_candidate_gene_snp_rank_ls)
			
			return_data = self.get_enrichment_pvalue_by_gw_looping(candidate_sample_size, top_loci_in_cumu_pos, candidate_gene_set, \
							genomeRBDict, cumuSpan2ChrRBDict=cumuSpan2ChrRBDict, \
							no_of_permutations=pd.no_of_permutations, no_of_min_breaks=pd.no_of_min_breaks, param_data=pd)
			pvalue = return_data.pvalue
			no_of_tests = return_data.no_of_tests
			no_of_tests_passed = return_data.no_of_tests_passed
			sys.stderr.write("%s pvalue: %s.\n"%(result_id, pvalue))
		if self.commit:
			db_250k.session.flush()
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = ResultPeakEnrichmentTest
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
