#!/usr/bin/env python
"""
Examples:

	#2011-3-16 run test on result 4634, list 28, top 500 loci, 20kb max distance from locus to gene 
	%s -e 4634 -l 28 -f 500 -s genomeDictTaxID3702MaxDist20kb.pickle -m 20000 -y15 -u yh -z banyan
	
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
__doc__ = __doc__%(sys.argv[0])
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
from pymodule.CNV import CNVCompareBySmallOverlapRatio, CNVSegmentBinarySearchTreeKey, CNVCompareByBigOverlapRatio, \
	CNVCompareByOverlapLen, get_overlap_ratio
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
						("result_id_ls", 1, ): [None, 'e', 1, 'comma/dash-separated results_by_gene id list, like 1,3-7'],\
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
						("min_big_overlap", 1, float): [0.3, 'a', 1, 'minimum of max(overlapFraction1, overlapFraction2) to declare a CNV is close to a gene'],\
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
	
	def translateChrPosDataObjectIntoCumuPos(self, top_loci, chrSpan2cumuStartRBDict=None):
		"""
		2011-4-22
			change chr_id2cumu_start to chrSpan2cumuStartRBDict
		2011-3-21
			top_loci has become a list of DataObject of GWR.
		2011-3-16
		
		"""
		sys.stderr.write("Translating %s loci from chr-span coordinates into cumu-span ..."%(len(top_loci)))
		top_loci_in_cumu_pos = []
		no_of_loci_skipped = 0
		compareIns = CNVCompareBySmallOverlapRatio(min_reciprocal_overlap=0.0000001)
		for top_locus in top_loci:
			chr = top_locus.chromosome
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=chr, \
							span_ls=[top_locus.position, top_locus.stop_position], \
							min_reciprocal_overlap=0.00000000000001,)
			node_ls = []
			chrSpan2cumuStartRBDict.findNodes(segmentKey, node_ls=node_ls, compareIns=compareIns)
			if len(node_ls)==0:
				no_of_loci_skipped += 1
			for node in node_ls:
				overlapData = get_overlap_ratio(segmentKey.span_ls, [node.key.start, node.key.stop])
				overlapFraction1 = overlapData.overlapFraction1
				overlapFraction2 = overlapData.overlapFraction2
				overlap_length = overlapData.overlap_length
				overlap_start_pos = overlapData.overlap_start_pos
				overlap_stop_pos = overlapData.overlap_stop_pos
				
				cumu_start = overlap_start_pos - node.key.start + 1 + node.value	#overlap_start_pos is in normal genome coordinates.
				cumu_stop = overlap_stop_pos - node.key.start + 1 + node.value
				top_loci_in_cumu_pos.append([cumu_start, cumu_stop])
		sys.stderr.write("%s loci skipped. now %s loci.\n"%(no_of_loci_skipped, len(top_loci_in_cumu_pos)))
		return top_loci_in_cumu_pos
	
	def translateCumuPosIntoChrPos(self, top_loci_in_cumu_pos, cumuSpan2ChrSpanRBDict=None, compareIns=None):
		"""
		2011-4-22
			adjust because chr_id2cumu_start is now 0-based.
		2011-4-22
			For CNVs, one (cumu_start, cumu_stop) could span multiple keys in cumuSpan2ChrSpanRBDict
		2011-3-16
		"""
		top_loci = []
		compareIns = CNVCompareBySmallOverlapRatio(min_reciprocal_overlap=0.0000001)
		for span in top_loci_in_cumu_pos:
			cumu_start, cumu_stop = span[:2]
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=0, \
							span_ls=[cumu_start, cumu_stop], \
							min_reciprocal_overlap=0.00000000000001,)
							#2010-8-17 overlapping keys are regarded as separate instances as long as they are identical.
			node_ls = []
			cumuSpan2ChrSpanRBDict.findNodes(segmentKey, node_ls=node_ls, compareIns=compareIns)
			if len(node_ls)==0:
				sys.stderr.write("(%s, %s) not found in cumuSpan2ChrSpanRBDict.\n"%(cumu_start, cumu_stop))
			for node in node_ls:
				chr, node_chr_start, node_chr_stop = node.value[:3]
				overlapData = get_overlap_ratio(segmentKey.span_ls, [node.key.start, node.key.stop])
				overlapFraction1 = overlapData.overlapFraction1
				overlapFraction2 = overlapData.overlapFraction2
				overlap_length = overlapData.overlap_length
				overlap_start_pos = overlapData.overlap_start_pos
				overlap_stop_pos = overlapData.overlap_stop_pos
				
				
				start = overlap_start_pos - node.key.span_ls[0] + node_chr_start	#overlap_start_pos is in cumu coordinates.
				stop = overlap_stop_pos - node.key.span_ls[0] + node_chr_start
				if stop>node_chr_stop:	#truncate it. shouldn't happen though
					stop = node_chr_stop
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
	
	def applyGWLoopToCumuPos(self, top_loci_in_cumu_pos, cumuSpan2ChrSpanRBDict, shift=None, compareIns=None):
		"""
		2011-3-12
		"""
		genome_size = max([max(rbNode.key.span_ls) for rbNode in cumuSpan2ChrSpanRBDict])
		if shift is None:
			shift = random.randint(1, genome_size)
		top_loci_in_cumu_pos_perm = (numpy.array(top_loci_in_cumu_pos)+shift)%genome_size	#modulo to recycle
		
		return self.translateCumuPosIntoChrPos(top_loci_in_cumu_pos_perm, cumuSpan2ChrSpanRBDict)
	
	
	def get_enrichment_pvalue_by_gw_looping(self, candidate_sample_size, top_loci_in_cumu_pos, candidate_gene_set=None, \
							genomeRBDict=None, cumuSpan2ChrSpanRBDict=None, no_of_permutations=20000, \
							no_of_min_breaks=30, param_data=None, compareIns=None):
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
			permuted_top_loci_in_chr_start_stop = self.applyGWLoopToCumuPos(top_loci_in_cumu_pos, cumuSpan2ChrSpanRBDict, \
															compareIns=compareIns)
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
	
	def getTranslationDataStructureForBackgroundLoci(self, db_250k, cnv_method_id=None, min_MAF=0.1):
		"""
		2011-4-22
			1. get all loci whose MAF is above min_MAF
			2. construct a (chr,start,stop) 2 cumu_start dictionary
			3. construct a (cumu_start, cumu_stop) 2 (chr, start, stop) RBDict
			
		"""
		sys.stderr.write("Getting translation structures between (chr, start, stop) and (cumu_start, cumu_stop) for cnv method %s ..."%\
						cnv_method_id)
		TableClass = Stock_250kDB.CNV
		query = TableClass.query.filter_by(cnv_method_id=cnv_method_id).order_by(TableClass.chromosome).order_by(TableClass.start)
		
		chrSpan2cumuStartRBDict = RBDict()
		cumuSpan2ChrSpanRBDict = RBDict()
		
		cumu_start = 0
		counter = 0
		real_counter = 0
		for row in query:
			counter += 1
			maf = min(row.frequency, 1-row.frequency)
			if maf<=min_MAF:
				continue
			
			real_counter += 1
			chrSpanKey = CNVSegmentBinarySearchTreeKey(chromosome=row.chromosome, \
							span_ls=[row.start, row.stop], \
							min_reciprocal_overlap=0.00000000000001,)
					#2010-8-17 overlapping keys are regarded as separate instances as long as they are not identical.
			chrSpan2cumuStartRBDict[chrSpanKey] = cumu_start	#cumu_start is 0-based
			
			size = row.stop-row.start+1
			span_ls=[cumu_start+1, cumu_start+size]
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=0, \
							span_ls=span_ls, \
							min_reciprocal_overlap=0.00000000000001,)
					#2010-8-17 overlapping keys are regarded as separate instances as long as they are not identical.
			if segmentKey not in cumuSpan2ChrSpanRBDict:
				cumuSpan2ChrSpanRBDict[segmentKey] = (row.chromosome, row.start, row.stop)
			else:
				sys.stderr.write("Error: %s of chr %s is already in cumuSpan2ChrSpanRBDict.\n"%(segmentKey, row.chromosome))
			
			cumu_start += size
		sys.stderr.write("%s out of %s CNVs are included. Done.\n"%(real_counter, counter))
		return PassingData(cumuSpan2ChrSpanRBDict=cumuSpan2ChrSpanRBDict, chrSpan2cumuStartRBDict=chrSpan2cumuStartRBDict)
	
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
		
		compareIns = CNVCompareByOverlapLen(min_overlap_len=100)	#any overlap is an overlap
		translationData = None
		for result_id in self.result_id_ls:
			#establish the map from cnv.id from chr_pos
			rm = Stock_250kDB.ResultsMethod.get(result_id)
			if not rm.cnv_method_id:
				sys.stderr.write("ResultsMethod %s doesn't have cnv_method_id. Skip.\n"%(result_id))
				continue
			if not db_250k._cnv_id2chr_pos:
				db_250k.cnv_id2chr_pos = rm.cnv_method_id
				translationData = self.getTranslationDataStructureForBackgroundLoci(db_250k, cnv_method_id=rm.cnv_method_id, min_MAF=self.min_MAF)
				if not translationData.chrSpan2cumuStartRBDict:
					sys.stderr.write("Error: translationData.chrSpan2cumuStartRBDict is empty for cnv method %s. exit.\n"%(rm.cnv_method_id))
					sys.exit(3)
			pd.db_id2chr_pos = db_250k.cnv_id2chr_pos
			
			candidate_gene_set = db_250k.dealWithCandidateGeneList(self.list_type_id, return_set=True)	#internal cache
			pd.candidate_gene_set = candidate_gene_set
			
			gwr = db_250k.getResultMethodContent(result_id, pdata=pd, min_value_cutoff=self.min_score)
		
			top_loci = gwr.getTopLoci(no_of_top_loci=self.no_of_top_loci, min_score=self.min_score)
			top_loci_in_cumu_pos = self.translateChrPosDataObjectIntoCumuPos(top_loci, translationData.chrSpan2cumuStartRBDict)
			top_loci_in_chr_pos = self.translateCumuPosIntoChrPos(top_loci_in_cumu_pos, translationData.cumuSpan2ChrSpanRBDict, \
													compareIns=compareIns)
			permData = self.prepareDataForPermutationRankTest(top_loci_in_chr_pos, genomeRBDict, pd, report=True)
			
			#m = self.dealWithNoOfSNPsAssociatedWithCandidateGeneList(pd.list_type_id, rm, pd)	#cache is internally going on
			#n = permData.no_of_total_snps - m
			
			candidate_sample_size = len(permData.captured_candidate_gene_set)
			non_candidate_sample_size = len(permData.non_candidate_gene_snp_rank_ls)
			
			return_data = self.get_enrichment_pvalue_by_gw_looping(candidate_sample_size, top_loci_in_cumu_pos, candidate_gene_set, \
							genomeRBDict, cumuSpan2ChrSpanRBDict=translationData.cumuSpan2ChrSpanRBDict, \
							no_of_permutations=pd.no_of_permutations, no_of_min_breaks=pd.no_of_min_breaks, param_data=pd,\
							compareIns=compareIns)
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
