#!/usr/bin/env python
"""
Examples:
	#output SNP plots by ranks according to results_by_gene's id=260 (one particular phenotype).
	DrawSNPRegion.py -e 260 -l 28 -L /Network/Data/250k/tmp-yh/call_method_17_LD_m0.2.tsv  -o /Network/Data/250k/tmp-yh/snp_region/
		-i phenotype_1_c10_f200
	
	#output SNP plots by ranks according to results_by_gene (analysis_method_id=7-Emma, call_method_id=17) (covering all phenotypes)
	DrawSNPRegion.py -l 28 -L /Network/Data/250k/tmp-yh/call_method_17_LD_m0.2.tsv -o /Network/Data/250k/tmp-yh/snp_region_all/
		-i phenotype_1_c10_f200
	
	#take LD_info and gene_annotation from a file contains the pickled data structure
	DrawSNPRegion.py -i ./banyan_fs/tmp/GWA_res_FT.csv -l 28 -D /Network/Data/250k/tmp-yh/call_method_17_LD_m0.1_n0.1_m40000
		-o /Network/Data/250k/tmp-yh/snp_region/  -j /tmp/at_gene_model_pickelf
	
	#get SNPs from ResultsByGene, which are associated with a top percentage of candidate genes
	DrawSNPRegion.py -q 1 -a 7 -y Top30PercCandidateGene  -l 28 -j /tmp/at_gene_model_pickelf
		-D /Network/Data/250k/tmp-yh/call_method_17_LD_pickle_dummy_n0.1_m20000 -o /tmp/snp_region -u yh 
	
	#replace LD axe with axe_snp_matrix, axe_strain_pca, axe_map, axe_phenotype
	#draw gene labels for every gene, not just candidate genes
	DrawSNPRegion.py -i /tmp/DrawSNPRegion_snps.csv -I /Network/Data/250k/tmp-yh/call_method_17.tsv -N /tmp/phenotype.tsv
		-l 28 -o /tmp/snp_region/ -j /tmp/at_gene_model_pickelf  -u yh -s
	
	#ditto but highlight imputed genotypes with blue color
	DrawSNPRegion.py -i /tmp/DrawSNPRegion_snps.csv -I /Network/Data/250k/tmp-yh/call_method_17.tsv
		-J /Network/Data/250k/tmp-yh/call_method_17_from_18_before_imputation.tsv -N /tmp/phenotype.tsv -l 28
		-o /tmp/snp_region/ -j /tmp/at_gene_model_pickelf -e 17 -u yh -s
	
	#database commit, specify plot_type_short_name by '-y'
	DrawSNPRegion.py -i $tmp_fname -I /Network/Data/250k/tmp-yh/call_method_17.tsv -N /Network/Data/250k/tmp-yh/phenotype.tsv
		-l 28 -o /Network/Data/250k/tmp-yh/snp_region -j /Network/Data/250k/tmp-yh/at_gene_model_pickelf -e 17 -u yh -s
		-y SuziHandPick20081209 -c
	
	#use CNV discrete calls (1/-1) as SNP matrix
	DrawSNPRegion.py -i /tmp/CycB25CycB23Region.tsv -I /Network/Data/250k/tmp-yh/CNV/call_method_17_CNV_array_intensity_norm_chr1_GADA_M10_discrete_m.3125.tsv
		-N /Network/Data/250k/tmp-yh/phenotype.tsv -l 28 -o /Network/Data/250k/tmp-yh/snp_region_CNV_M10_discrete_m.3125
		-j /Network/Data/250k/tmp-yh/at_gene_model_pickelf -e 17 -u yh -s -V 2
	
	#use CNV amplitude (float) output as SNP matrix
	DrawSNPRegion.py -i /tmp/CycB25CycB23Region.tsv -I /Network/Data/250k/tmp-yh/CNV/call_method_17_CNV_array_intensity_norm_chr1_GADA_M10_out_amp.tsv
		-N /Network/Data/250k/tmp-yh/phenotype/phenotype.tsv -l 28 -o /Network/Data/250k/tmp-yh/snp_region_CNV_M10_amp
		-j /Network/Data/250k/tmp-yh/at_gene_model_pickelf -e 17 -u yh -s -V 3
	
	# 2010-2-2 draw SNP regions with gene labels (-s), exclude accession with NA phenotype from the SNP matrix (-E), use KW & RF scores (-a 1,6)
	DrawSNPRegion.py -i /tmp/AdnanePNASRegionsToDraw.csv -I /Network/Data/250k/db/dataset/call_method_32.tsv 
		-N /Network/Data/250k/tmp-yh/phenotype/phenotype.tsv -l 148 -o /Network/Data/250k/tmp-yh/snp_region
		-j /Network/Data/250k/tmp-yh/at_gene_model_pickelf -e 32 -u yh -s -a 1,6 -E
	
Description:
	2008-09-24 program to draw pvalues, gene-models, LD around one SNP.
		Top panel is pvalues from all different methods. Margarita and RF's values will be normalized in range of KW.
		Middle panel is gene model. Displaying CDS, intron, strand.
		Bottom panel is haplotype plot representing all SNPs in that range.
	
	The input file contains at least 3 columns, chromosome/chr, position/pos/start, phenotype_id/phenot_id, stop (optional). tab or comma-delimited.
		1. The 1st row is a header telling which column is which. No particular order is required.
		2. The name could either be full (chromosome) or short(chr).
		3. Output of FindTopSharedGenes.py conforms to this standard.
	2009-3-10 call_method_id (-e) has to match the input file (-I). GWA results would get into trouble in matching SNPs from the input SNP matrix.
	
	argument snp_matrix_data_type (-V)
	if snp_matrix_data_type==1:		#SNP data (snps are subset of what db has).
		need_convert_alleles2binary = True
		useAlleleToDetermineAlpha = False
	elif snp_matrix_data_type==2:	#CNV discrete call data (1 or -1);
		need_convert_alleles2binary = True
		useAlleleToDetermineAlpha = False
	elif snp_matrix_data_type==3:	#2009-3-23 for CNV amplitude data
		need_convert_alleles2binary = False
		useAlleleToDetermineAlpha = True
	elif snp_matrix_data_type==4:	#2009-3-27, for arbitrary non-diallelic SNP matrix
		need_convert_alleles2binary = False
		useAlleleToDetermineAlpha = False
	else:	# all other cases
		need_convert_alleles2binary = True
		useAlleleToDetermineAlpha = False
"""

import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, cPickle
import warnings, traceback
import matplotlib; matplotlib.use("Agg")	#to avoid popup and collapse in X11-disabled environment
import matplotlib as mpl
import pylab
import ImageColor
import numpy, StringIO
from matplotlib.patches import Polygon, CirclePolygon
from pymodule import PassingData, figureOutDelimiter, getColName2IndexFromHeader, getListOutOfStr, read_data, SNPData, SNPInfo, SNP
from pymodule import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
from pymodule import RBDict
from pymodule.Genome import GeneModel
from pymodule.plot.yh_matplotlib_artists import ExonIntronCollection
from pymodule import outputMatrixInLatexTable, outputFigureInLatex
from pymodule import RBTree, RBDict, RBTreeIter, RBList, RBNode
from variation.src import Stock_250kDB
from variation.src import AbstractVariationMapper
from variation.src.enrichment.GeneListRankTest import GeneListRankTest	#GeneListRankTest.getGeneList()
from variation.src.common import getEcotypeInfo
from variation.src.association.Kruskal_Wallis import Kruskal_Wallis
from PlotGroupOfSNPs import PlotGroupOfSNPs

class LD_statistic(object):
	"""
	2008-09-30
		a class to get label/name conveniently given which LD_statistic is chosen
	"""
	which_LD_statistic2name = {1:'r2', 2:'D_prime', 3:'D'}
	which_LD_statistic2label = {1:'r^2', 2:"|D^'|", 3:'D'}
	
	def get_name(cls, which_LD_statistic):
		return cls.which_LD_statistic2name.get(which_LD_statistic)
	get_name = classmethod(get_name)
	
	def get_label(cls, which_LD_statistic):
		return cls.which_LD_statistic2label.get(which_LD_statistic)
	get_label = classmethod(get_label)


class SNPPassingData(PassingData):
	chromosome = None
	position = None
	def __cmp__(self, other):
		"""
		2008-09-30
			define how to compare to itself
		"""
		return cmp((self.chromosome, self.position), (other.chromosome, other.position))


class GeneSNPPhenotypeAssoData(object):
	"""
	2008-10-02
		the SNP could also denote a region
	2008-10-02
		a data structure to hold gene-snp-phenotype association data got from input file
		iterator not fully tested, ...
	"""
	def __init__(self, **keywords):
		self.db_250k = None	#needs this to be connected.
		for argument_key, argument_value in keywords.iteritems():
			setattr(self, argument_key, argument_value)
		self.snp2phenotype_id_ls = {}
		self.gene_id2snp_data = {}
		self.gene_id2phenotype_id_set = {}
		
		self._gene_id_ls = None
		self.gene_index = -1
		self.iter_block = []
		self.phenotype_id2phenotype = {}	#caches for return_matrix_of_snp_descriptions()
	
	def return_phenotype_short_name(self, phenotype_id):
		"""
		2008-10-1
			just for convenience, for DrawSNPRegion.run()
		"""
		phenotype = self.phenotype_id2phenotype.get(phenotype_id)
		if not phenotype:
			phenotype = Stock_250kDB.PhenotypeMethod.get(phenotype_id)
			self.phenotype_id2phenotype[phenotype_id] = phenotype
		return phenotype.short_name
	
	def return_phenotype_id2snp_ls(self):
		phenotype_id2snp_ls = {}
		for snp, phenotype_id_set in self.snp2phenotype_id_ls.iteritems():
			for phenotype_id in phenotype_id_set:
				if phenotype_id not in phenotype_id2snp_ls:
					phenotype_id2snp_ls[phenotype_id] = []
				phenotype_id2snp_ls[phenotype_id].append(snp)
		self.phenotype_id2snp_ls = phenotype_id2snp_ls
		return phenotype_id2snp_ls
	
	def no_of_genes(self):
		return len(self.gene_id2snp_data)
	
	no_of_genes = property(no_of_genes)
	
	def no_of_snps(self):
		return len(self.snp2phenotype_id_ls)
	
	no_of_snps = property(no_of_snps)
	
	def gene_id_ls(self):
		if self._gene_id_ls == None:
			self._gene_id_ls = self.gene_id2snp_data.keys()
			self._gene_id_ls.sort()
		return self._gene_id_ls
	gene_id_ls = property(gene_id_ls)
	
	def addGeneSNPPhenotype(self, gene_id, snp, phenotype_id):
		"""
		2012.3.23
			snp has a new field. snp.fileNamePrefix
		"""
		this_snp = (snp.chromosome, snp.position, snp.snps_id, snp.stop, snp.fileNamePrefix)
		if gene_id is not None:
			if gene_id not in self.gene_id2snp_data:
				self.gene_id2snp_data[gene_id] = {}
			snp_data = self.gene_id2snp_data[gene_id]
			if this_snp not in snp_data:
				snp_data[this_snp] = set()
			snp_data[this_snp].add(phenotype_id)
			
			if gene_id not in self.gene_id2phenotype_id_set:
				self.gene_id2phenotype_id_set[gene_id] = set()
			self.gene_id2phenotype_id_set[gene_id].add(phenotype_id)
		
		if this_snp not in self.snp2phenotype_id_ls:
			self.snp2phenotype_id_ls[this_snp] = set()
		self.snp2phenotype_id_ls[this_snp].add(phenotype_id)
	
	def get_no_of_phenotypes_given_gene_id(self, gene_id):
		"""
		"""
		phenotype_id_set = self.gene_id2phenotype_id_set.get(gene_id)
		if phenotype_id_set:
			return len(phenotype_id_set)
		else:
			return 0
	
	snp_desc_names = ['chromosome', 'position', 'phenotypes']
	def return_matrix_of_snp_descriptions(self, gene_id):
		"""
		2008-10-02
			return a 2-D list
		"""
		snp_data = self.gene_id2snp_data.get(gene_id)
		matrix_of_snp_descriptions = []
		if not snp_data:
			return matrix_of_snp_descriptions
		snp_ls = snp_data.keys()
		snp_ls.sort()
		for this_snp in snp_ls:
			snp_desc_ls = []
			for snp_desc_name in self.snp_desc_names:
				if snp_desc_name=='phenotypes':
					element_ls = []
					for phenotype_id in snp_data[this_snp]:
						phenotype = self.phenotype_id2phenotype.get(phenotype_id)
						if not phenotype:
							phenotype = Stock_250kDB.PhenotypeMethod.get(phenotype_id)
							self.phenotype_id2phenotype[phenotype_id] = phenotype
						element_ls.append('%s(id %s)'%(phenotype.short_name, phenotype_id))
					element = ', '.join(element_ls)
				elif snp_desc_name=='chromosome':
					element = this_snp[0]
				elif snp_desc_name=='position':
					element = this_snp[1]
				else:
					element = ''
				snp_desc_ls.append(element)
			matrix_of_snp_descriptions.append(snp_desc_ls)
		return matrix_of_snp_descriptions
	
	def return_snp_data_given_gene_id(self, gene_id):
		return self.gene_id2snp_data.get(gene_id)
	
	def __iter__(self):
		return self
	
	def next(self):
		"""
		2008-10-02
			iterator not fully tested, ...
		"""
		self.read()
		return self.iter_block
	
	def read(self):
		"""
		2008-10-02
			iterator not fully tested, ...
		"""
		self.gene_index += 1
		if self.gene_index<self.no_of_genes:
			gene_id = self.gene_id_ls[self.gene_index]
			self.iter_block = [gene_id, self.gene_id2snp_data[gene_id]]
		if self.gene_index==len(self.gene_id_ls):
			self.gene_index = -1
			raise StopIteration
		
class DrawSNPRegion(PlotGroupOfSNPs):
	__doc__ = __doc__
	option_default_dict = AbstractVariationMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('LD_fname', 0, ): [None, 'L', 1, 'the file containing LD info, output of MpiLD.py', ],\
						("min_distance", 1, int): [20000, 'm', 1, 'minimum distance allowed from the SNP to gene'],\
						("get_closest", 0, int): [0, 'g', 0, 'only get genes closest to the SNP within that distance'],\
						('min_MAF', 0, float): [0.1, 'n', 1, 'minimum Minor Allele Frequency.'],\
						("list_type_id_list", 1, ): [None, 'l', 1, 'comma/dash separated list of Gene list type. must be in table gene_list_type.'],\

						('rbg_results_directory', 0, ):[os.path.expanduser('~/panfs/db/results_by_gene/'), '', 1, 'The rbg results directory. Default is None. use the one given by db.'],\
						("input_fname", 0, ): [None, 'i', 1, 'Filename which contains at least 3 columns: chromosome, position, phenotype_id. if not given, get from db'],\
						('snp_matrix_fname', 0, ): ['', 'I', 1, 'genotype matrix. Strain X SNP format.', ],\
						('snp_matrix_before_impute_fname', 0, ): ['', 'J', 1, 'genotype matrix, before imputation. Strain X SNP format, specify to check which positions are imputed. (colored in blue)', ],\
						('phenotype_fname', 0, ): [None, 'N', 1, 'phenotype file, if snp_matrix_fname is given, this is needed as well.', ],\
						("output_dir", 1, ): [None, 'o', 1, 'directory to store all images'],\
						('call_method_id', 0, int):[17, 'e', 1, 'Restrict results based on this call_method. Default is no such restriction.'],\
						('analysis_method_id_ls', 0, ):['7', 'a', 1, 'Restrict results based on this list of analysis_method. Default is no such restriction.'],\
						('phenotype_method_id_ls', 0, ):[None, 'q', 1, 'Restrict results based on this list of phenotype_method. Default is no such restriction.'],\
						('no_of_top_hits', 1, int): [1000, 'f', 1, 'how many number of top hits based on score or -log(pvalue).'],\
						("which_LD_statistic", 1, int): [2, 'w', 1, 'which LD_statistic to plot, 1=r2, 2=|D_prime|, 3=|D|'],\
						("LD_info_picklef", 0, ): [None, 'D', 1, 'given the option, If the file does not exist yet, store a pickled LD_info into it (min_MAF and min_distance will be attached to the filename). If the file exists, load LD_info out of it.'],\
						("gene_annotation_picklef", 0, ): [None, 'j', 1, 'given the option, If the file does not exist yet, store a pickled gene_annotation into it. If the file exists, load gene_annotation out of it.'],\
						("latex_output_fname", 0, ): [None, 'x', 1, 'a file to store the latex of gene description tables and figures'],\
						("label_gene", 0, int): [0, 's', 0, 'toggle this to label every gene by gene symbol (candidate genes in blue, others in black). \
	Otherwise only candidate genes are labeled.'],\
						("draw_LD_relative_to_center_SNP", 0, int): [1, '', 0, 'toggle this to draw LD between other SNPs and the center SNP'],\
						("drawMap", 0, int): [0, 'M', 0, 'toggle this to draw a map (require toolkits.basemap) right to the SNP matrix to show geo-location of each accession/row'],\
						("plot_type_short_name", 0, ): [None, 'y', 1, 'A short name (<256 chars) to characterize the type of all the plots. it could have existed in table SNPRegionPlotType'],\
						('snp_matrix_data_type', 1, int):[1, 'V', 1, 'what type of data in snp_matrix_fname? \n\
	1: SNP data (snps are subset of what db has). \n \
	2: CNV discrete call data (1 or -1); \n \
	3: CNV amplitude data; \n\
	4: arbitrary non-diallelic SNP data (SNPs are not required to be included in db). \n\
	input file has to contain start, stop position for 2 & 3. Providing central position in input works only for SNP data.'],\
						('exclude_accessions_with_NA_phenotype', 0, int):[0, 'E', 0, 'toggle to exclude accessions from the SNP matrix to be drawn'],\
						('markSNPInGenes', 0, int):[0, 'P', 0, 'toggle to mark the SNPs in gene models'],\
						('drawGrid', 0, int):[0, 'G', 0, 'toggle to draw a grid for the SNP-matrix/haplotype'],\
						('phenotypeDrawType', 1, int):[1, '', 1, 'phenotype drawing type. 1: same color but the length varies by magnitude. 2: same length but color varies'],\
						('snpInfoPickleFname', 1, ): ['', 'F', 1, 'The file to contain pickled SNPInfo.'],\
						('logFilename', 0, ): [None, 'A', 1, 'file to contain logs. use it only if this program is at the end of pegasus workflow \
					and has no output file'],\
						})
	
	def __init__(self,  **keywords):
		"""
		2008-11-30
			now it's inherited from PlotGroupOfSNPs
		2008-09-24
		"""
		GeneListRankTest.__init__(self, **keywords)
		#from pymodule import ProcessOptions
		#self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		self.analysis_method_id_ls = getListOutOfStr(self.analysis_method_id_ls, data_type=int)
		if self.phenotype_method_id_ls:
			self.phenotype_method_id_ls = getListOutOfStr(self.phenotype_method_id_ls, data_type=int)
		
		if self.list_type_id_list:
			self.list_type_id_list = getListOutOfStr(self.list_type_id_list, data_type=int)
		else:
			self.list_type_id_list = []
		
	def get_LD(self, LD_fname, min_MAF, min_gap=40000):
		"""
		2008-09-29
			add min_MAF, min_gap, which_LD_statistic
			min_gap is the distance allowed between two SNPs
			which_LD_statistic returns appropriate LD statistic
		2008-09-24
		"""
		sys.stderr.write("Reading in LD info from %s ...\n"%(LD_fname))
		snp_pair2r2 = {}
		reader = csv.reader(open(LD_fname), delimiter='\t')
		col_name2index = getColName2IndexFromHeader(reader.next())
		counter = 0
		for row in reader:
			snp1 = row[col_name2index['snp1']].split('_')
			snp1 = map(int, snp1)
			snp2 = row[col_name2index['snp2']].split('_')
			snp2 = map(int, snp2)
			if snp1[0]==snp2[0] and abs(snp1[1]-snp2[1])<=min_gap:	#on the same chromosome, and less than a certain distance
				allele1_freq = float(row[col_name2index['allele1_freq']])
				allele2_freq = float(row[col_name2index['allele2_freq']])
				if allele1_freq>=min_MAF and allele2_freq>=min_MAF:	#meet the minimum minor-allele-frequency
					r2 = float(row[col_name2index['r2']])
					D_prime = float(row[col_name2index['D_prime']])
					D = float(row[col_name2index['D']])
					if snp1<snp2:
						snp_pair = (snp1[0], snp1[1], snp2[0], snp2[1])
					else:
						snp_pair = (snp2[0], snp2[1], snp1[0], snp1[1])
					snp_pair2r2[snp_pair] = PassingData(r2=r2, D_prime=D_prime, D=D)
			counter += 1
			if counter%100000==0:
				sys.stderr.write('%s\t%s'%('\x08'*100, counter))
				if self.debug>0:
					break
					pass
		LD_info = PassingData(snp_pair2r2=snp_pair2r2)
		sys.stderr.write("%s entries. Done.\n"%len(snp_pair2r2))
		return LD_info
	
	def dealLD_info(self, LD_info_picklef, LD_fname=None, min_MAF=0.1, min_distance=20000):
		"""
		2008-09-30
			if the LD_info_picklef does not exist yet, store a pickled LD_info into it. If the file exists, load LD_info out of it.
		"""
		sys.stderr.write("Dealing with LD_info ...")
		if LD_info_picklef:
			if os.path.isfile(LD_info_picklef):	#if this file is already there, suggest to un-pickle it.
				picklef = open(LD_info_picklef)
				LD_info = cPickle.load(picklef)
				del picklef
			else:	#if the file doesn't exist, but the filename is given, pickle into it
				LD_info = self.get_LD(LD_fname, min_MAF, min_distance*2+100)	#min_gap is 2 X min_distance + a few more bases
				picklef = open('%s_n%s_m%s'%(LD_info_picklef, min_MAF, min_distance), 'w')
				cPickle.dump(LD_info, picklef, -1)
				picklef.close()
		else:
			LD_info = self.get_LD(LD_fname, min_MAF, min_distance*2+100)
		sys.stderr.write("Done.\n")
		return LD_info
	
	@classmethod
	def getGeneAnnotation(cls, tax_id=3702, cls_with_db_args=None):
		"""
		2011-1-27
			db.get_gene_id2model() returns 3 objects now. Make sure only the first 2 are needed.
		2009-1-15
			add argument cls_with_db_args to get drivername, username, dbname, etc.
		2008-12-16
			become classmethod
			add option tax_id
		2008-10-01
			substitute GenomeBrowser.get_gene_id2model() with GenomeDB.GenomeDatabase.get_gene_id2model()
			
		2008-09-24
		"""
		from pymodule import GenomeDB
		if cls_with_db_args is None:
			cls_with_db_args = cls
		db = GenomeDB.GenomeDatabase(drivername=cls_with_db_args.drivername, username=cls_with_db_args.db_user,
				   password=cls_with_db_args.db_passwd, hostname=cls_with_db_args.hostname, database='genome', \
				   schema=getattr(cls_with_db_args, 'schema', None))
		db.setup(create_tables=False)
		gene_id2model, chr_id2gene_id_ls, = db.get_gene_id2model(tax_id=tax_id)[:2]
		gene_annotation = PassingData()
		gene_annotation.gene_id2model = gene_id2model
		gene_annotation.chr_id2gene_id_ls = chr_id2gene_id_ls
		return gene_annotation
	
	@classmethod	
	def dealWithGeneAnnotation(cls, gene_annotation_picklef, tax_id=3702, cls_with_db_args=None):
		"""
		2010-4-20
			use os.access to test whether the gene_annotation_picklef is writeable or not before dumping the db-fetched annotaion into it.
		2009-1-15
			add argument cls_with_db_args to be passed to getGeneAnnotation()
		2008-12-16
			become classmethod
			add option tax_id
		2008-10-01
			similar to dealLD_info()
		"""
		sys.stderr.write("Dealing with gene_annotation %s ..." %(gene_annotation_picklef))
		if gene_annotation_picklef:
			if os.path.isfile(gene_annotation_picklef):	#if this file is already there, suggest to un-pickle it.
				from pymodule import utils	#2012.3.26 "utils" has to be present for cPickle.load() to work
				picklef = open(gene_annotation_picklef)
				gene_annotation = cPickle.load(picklef)
				del picklef
			else:	#if the file doesn't exist, but the filename is given, pickle into it
				gene_annotation = cls.getGeneAnnotation(tax_id=tax_id, cls_with_db_args=cls_with_db_args)	#min_gap is 2 X min_distance + a few more bases
				try: #os.access(gene_annotation_picklef, os.W_OK):	# 2010-4-20 make sure the file is writable
					picklef = open(gene_annotation_picklef, 'w')
					cPickle.dump(gene_annotation, picklef, -1)
					picklef.close()
				except:
					sys.stderr.write('Except while trying to pickle gene_annotation to a file: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
		else:
			gene_annotation = cls.getGeneAnnotation(tax_id=tax_id, cls_with_db_args=cls_with_db_args)
		sys.stderr.write(" %s genes.\n"%(len(gene_annotation.gene_id2model)))
		return gene_annotation
	
	@classmethod
	def getSimilarGWResultsGivenResultsByGene(cls, db_250k=None, phenotype_method_id=None, call_method_id=None, results_directory=None, \
											analysis_method_id_ls=[1,7]):
		"""
		2012.3.7
			add argument db_250k
			pass db_id2chr_pos to db_250k.getResultMethodContent()
		2009-3-24
			default value of analysis_method_id_ls is changed from [1,7,17,18,19] to [1,7]
		2008-12-14
			default value of analysis_method_id_ls is changed from [1,5,6,7] to [1,7,17,18,19]
		2008-09-30
			use phenotype_method_id and call_method_id to find out all genome-wide results
		2008-09-24
		"""
		sys.stderr.write("Getting results with phenotype=%s and call_method=%s ..."%(phenotype_method_id, call_method_id))
		if getattr(cls, 'debug', False):
			analysis_method_id_set = set([1,7])
		else:
			analysis_method_id_set = set(analysis_method_id_ls)
		rows = Stock_250kDB.ResultsMethod.query.filter_by(phenotype_method_id=phenotype_method_id, call_method_id=call_method_id)
		analysis_method_id2gwr = {}
		pd = PassingData(min_MAF=0,\
						results_directory=results_directory, \
						need_chr_pos_ls=0,)
		counter = 0
		for rm in rows:
			if rm.analysis_method_id in analysis_method_id_set:
				if rm.cnv_method_id and not db_250k._cnv_id2chr_pos:
					db_250k.cnv_id2chr_pos = rm.cnv_method_id
					pd.db_id2chr_pos = db_250k.cnv_id2chr_pos
				elif rm.call_method_id:
					pd.db_id2chr_pos = db_250k.snp_id2chr_pos
				#genome_wide_result = db_250k.getResultMethodContent(rm.id, pdata=pd)
				genome_wide_result = db_250k.getResultMethodContent(rm.id, results_directory=results_directory, min_MAF=0, \
												construct_chr_pos2index=True, pdata=pd)
				analysis_method_id2gwr[rm.analysis_method_id] = genome_wide_result
				counter += 1
		sys.stderr.write("%s gwas results.\n"%(counter))
		return analysis_method_id2gwr
	
	@classmethod
	def getSNPInfo(cls, db, locus_type_id=1):
		"""
		2012.3.19 moved to Stock_250kDB
		2012.3.7
			add end_position into chr_pos
			todo :
				#. have to somehow separate snp-gwas loci and cnv-gwas loci , because they overlap
				#. generate RB dictionary to store all loci for fast retrieval based on (chr, start, stop)
		2009-2-18
			replace PassingData with class SNPInfo from pymodule.SNP to wrap up all data
		2009-1-22
			order by chromosome, position
			no chr_pos2snps_id, can get snps_id from chr_pos2index => data_ls
		2009-1-5
			add chr_pos2snps_id
		2008-09-24
			in order
		"""
		return db.getSNPInfo(locus_type_id=locus_type_id)
	
	@classmethod
	def add_mid_point(cls, chr_pos_ls, chr_pos2adjacent_window):
		"""
		2012.3.8
			defunct. findSNPsInRegion() has this functionality built in.
		2008-09-24
			called by getSNPsAroundThisSNP()
		"""
		new_chr_pos = chr_pos_ls[-1]
		old_chr_pos = chr_pos_ls[-2]
		if old_chr_pos not in chr_pos2adjacent_window:
			chr_pos2adjacent_window[old_chr_pos] = []
		if new_chr_pos not in chr_pos2adjacent_window:
			chr_pos2adjacent_window[new_chr_pos] = []
		mid_point = (new_chr_pos[1]+old_chr_pos[1])/2.
		chr_pos2adjacent_window[old_chr_pos].append(mid_point)
		chr_pos2adjacent_window[new_chr_pos].append(mid_point)
	
	@classmethod
	def getSNPsAroundThisSNP(cls, this_snp, snp_info, min_distance=20000):
		"""
		2012.3.8
			defunct. replaced by findSNPsInRegion().
		2008-09-24
		"""
		sys.stderr.write("\t Get SNPs around this snp ...")
		#chr_pos = snp_info.chr_pos_ls[snp_info.snps_id2index[this_snp.snps_id]]
		chromosome, position = this_snp.chromosome, this_snp.position
		chr_pos_ls = []
		chr_pos2adjacent_window = {}
		j = 0
		for i in range(min_distance*2):
			new_pos = position - min_distance + i
			new_chr_pos = (chromosome, new_pos)
			if new_chr_pos in snp_info.chr_pos2index:
				chr_pos_ls.append(new_chr_pos)
				if j!=0:
					cls.add_mid_point(chr_pos_ls, chr_pos2adjacent_window)
				j += 1
		#deal with the leftest point of the 1st chr_pos
		chr_pos = chr_pos_ls[0]
		window_size = chr_pos2adjacent_window[chr_pos][0]-chr_pos[1]
		chr_pos2adjacent_window[chr_pos] = [chr_pos[1]-window_size, chr_pos[1]+window_size]
		
		#deal with the rightest point of the 1st chr_pos
		chr_pos = chr_pos_ls[-1]
		window_size = chr_pos[1] - chr_pos2adjacent_window[chr_pos][0]
		chr_pos2adjacent_window[chr_pos] = [chr_pos[1]-window_size, chr_pos[1]+window_size]
		snp_region = PassingData(chr_pos_ls=chr_pos_ls, chr_pos2adjacent_window=chr_pos2adjacent_window, center_snp=this_snp)
		sys.stderr.write("Done.\n")
		return snp_region
	
	analysis_method_id2color = {1:'r',
							2: 'r',
							5:'g',
							6:'c',
							7:'b',
							17:'g',
							18:'c',
							19:'k',
							32:'#87CEFA'}	#lightskyblue: #87CEFA	R=135 G=206	B=250 ACCESS=16436871
	
	# 2010-2-2 cls.analysis_method_id2name is deprecated in this program.
	# 	the analysis method name is now retrieved directly from db (genome_wide_result.rm.analysis_method.short_name)
	analysis_method_id2name = {1:'KW',
							2:'FisherExact',
							5:'Margarita',
							6:'RF',
							7:'Emma',
							17:'LM_with_PC12',
							18:'LM_with_PC1_10',
							19:'LM_with_PC1_34'}
	
	@classmethod
	def getXY(cls, snps_within_this_region=None, genome_wide_result=None,\
			analysis_method_id2gwr=None, analysis_method_id=None, LD_info=None, which_LD_statistic=2):
		"""
		2012.11.14 modernizes it so that PlotAssociationLocus.py could use it.
		2012.3.7
		2008-10-1
			return LD value (with regard to the center SNP)
		2008-09-24
			of GW results, get values for each SNP position, adjust value for analysis_method_id=5,6
		"""
		x_ls = []
		y_ls = []
		if genome_wide_result is None and analysis_method_id and analysis_method_id2gwr:
			genome_wide_result = analysis_method_id2gwr[analysis_method_id]
		if genome_wide_result is None and LD_info is None:	#2008-10-01
			return x_ls, y_ls
		
		if analysis_method_id2gwr and 1 in analysis_method_id2gwr:
			ref_gwr = analysis_method_id2gwr[1]
		elif analysis_method_id2gwr and 7 in analysis_method_id2gwr:
			ref_gwr = analysis_method_id2gwr[7]
		else:
			ref_gwr = None
		for chr_pos in snps_within_this_region.chr_pos_ls:
			if genome_wide_result:
				data_obj = genome_wide_result.get_data_obj_by_chr_pos(chr_pos[0], chr_pos[1], stopPosition=chr_pos[2])	#2012.3.7
			else:	#2008-10-1 fake a data_obj for LD
				chr_pos1 = (snps_within_this_region.center_snp.chromosome, snps_within_this_region.center_snp.position)
				chr_pos2 = chr_pos
				if chr_pos1<chr_pos2:
					snp_pair = (chr_pos1[0], chr_pos1[1], chr_pos2[0], chr_pos2[1])
				else:
					snp_pair = (chr_pos2[0], chr_pos2[1], chr_pos1[0], chr_pos1[1])
				if snp_pair in LD_info.snp_pair2r2:
					LD_stat = getattr(LD_info.snp_pair2r2[snp_pair], LD_statistic.get_name(which_LD_statistic), None)
					LD_stat = abs(LD_stat)	#D_prime, D need abs()
				elif chr_pos1==chr_pos2:	#fake a perfect LD point
					LD_stat = 1.
				else:	#skip this point
					continue
				data_obj = PassingData(value=LD_stat)
				
			if data_obj is not None:
				x_ls.append(chr_pos[1])
				if (analysis_method_id==5 or analysis_method_id==6) and ref_gwr:
					value = (data_obj.value-genome_wide_result.min_value)/(genome_wide_result.max_value-genome_wide_result.min_value)*(ref_gwr.max_value-ref_gwr.min_value)
				else:
					value = data_obj.value
				y_ls.append(value)
		return x_ls, y_ls
	
	@classmethod
	def drawPvalue(cls, axe_pvalue, axe_to_put_pvalue_legend, axe_LD_center_SNP, snps_within_this_region, analysis_method_id2gwr, \
				LD_info=None, \
				which_LD_statistic=2, draw_LD_relative_to_center_SNP=0, legend_loc='upper right'):
		"""
		2010-4-19
			add the ylable for axe_LD_center_SNP only when LD stat is actually drawn
		2010-4-15 set the legend box semi-transparent
		2010-4-14
			add "fancybox=True" to axe_to_put_pvalue_legend.legend()
		2010-4-13
			alpha of scatter in matplotlib 0.99.0 controls both edge and face. set it to 1.
		2010-2-2
			get the analysis method name from db
		2008-10-1
			refine the LD line on axe_LD_center_SNP, with markersize=2 and etc.
		2008-10-1
			draw LD of other SNPs w.r.t the center SNP in a different axe
		2008-09-24
		"""
		sys.stderr.write("\t Drawing pvalues  ...")
		analysis_method_id_ls = analysis_method_id2gwr.keys()
		analysis_method_id_ls.sort()
		pscatter_ls = []
		legend_ls = []
		for analysis_method_id in analysis_method_id_ls:
			genome_wide_result = analysis_method_id2gwr[analysis_method_id]
			x_ls, y_ls = cls.getXY(snps_within_this_region, analysis_method_id2gwr, analysis_method_id)
			if x_ls and y_ls:
				pscatter = axe_pvalue.scatter(x_ls, y_ls, s=4, linewidth=0.6, \
											edgecolor=cls.analysis_method_id2color.get(analysis_method_id, 'b'), facecolor='none', alpha=1)
				
				#legend_ls.append(cls.analysis_method_id2name[analysis_method_id])
				legend_ls.append(genome_wide_result.rm.analysis_method.short_name)	# 2010-2-2 get the analysis method name from db
				pscatter_ls.append(pscatter)
		if LD_info and draw_LD_relative_to_center_SNP:	#draw LD with regard to the center SNP
			x_ls, y_ls = cls.getXY(snps_within_this_region, LD_info=LD_info, which_LD_statistic=which_LD_statistic)
			if x_ls and y_ls:
				apl = axe_LD_center_SNP.plot(x_ls, y_ls, 'k-o', linewidth=0.5, markersize=2, alpha=0.3, markerfacecolor='w')
				legend_ls.append(r'$%s$ with center SNP'%LD_statistic.get_label(which_LD_statistic))
				pscatter_ls.append(apl[0])
				# 2010-4-19 add the ylable for axe_LD_center_SNP only when LD stat is actually drawn
				axe_LD_center_SNP.set_ylabel(r'$%s$ with center SNP'%LD_statistic.get_label(which_LD_statistic))				
		axe_pvalue.set_ylabel(r'-log(pvalue)/normalized score')
		legend = axe_to_put_pvalue_legend.legend(pscatter_ls, legend_ls, loc=legend_loc, fancybox=True, handlelength=0.02)	#cut the legend length to 0.02, default 0.05 (5% of x-axis).
		frame  = legend.get_frame()
		frame.set_alpha(0.60)	#2010-4-15 set the legend box semi-transparent
		
		sys.stderr.write("Done.\n")
		#return legend_ls
	
	gene_desc_names = ['gene_id', 'gene_symbol', 'type_of_gene', 'chr', 'start', 'stop', 'protein_label', 'protein_comment', 'protein_text']
	
	@classmethod
	def returnGeneDescLs(cls, gene_desc_names, gene_model, gene_commentary=None, cutoff_length=200, replaceNoneElemWithEmptyStr=0):
		"""
		2009-2-4
			add argument replaceNoneElemWithEmptyStr, which toggles the option whether to set None element to empty string ('').
		2008-10-02
		"""
		#2008-10-01	get the gene descriptions
		gene_desc_ls = []
		if gene_commentary is None:	#take gene_model
			gene_commentary = gene_model
		for gene_desc_name in gene_desc_names:
			if gene_desc_name=='protein_label':	#if not available, get the rna label
				element = getattr(gene_commentary, gene_desc_name, None)
				if not element:
					element = getattr(gene_commentary, 'label', '')
					if cutoff_length is not None and type(element) ==str:
						element = element[:cutoff_length/2]
				if cutoff_length is not None and type(element) ==str:	#cut label length half
						element = element[:cutoff_length/2]
			elif gene_desc_name=='protein_comment':
				element = getattr(gene_commentary, gene_desc_name, None)
				if not element:
					element = getattr(gene_commentary, 'comment', '')
					if cutoff_length is not None and type(element) ==str:
						element = element[:cutoff_length/2]
			elif gene_desc_name=='gene_symbol' or gene_desc_name=='chromosome':	#have to get it from the gene_model
				element = getattr(gene_model, gene_desc_name, '')
			elif gene_desc_name=='type_of_gene':
				if getattr(gene_commentary, 'protein_label', None) is not None:
					element = getattr(gene_model, gene_desc_name, '')
				else:	#it doesn't have protein, get gene_commentary_type
					element = getattr(gene_commentary, 'gene_commentary_type', '')
			elif gene_desc_name=='chr':
				element = getattr(gene_model, 'chromosome', '')
			else:
				element = getattr(gene_commentary, gene_desc_name, '')
				if not element:	#2010-8-19 if gene_commentary doesn't have anything, try gene_model itself
					element = getattr(gene_model, gene_desc_name, '')
			if cutoff_length is not None and type(element)==str and len(element)>cutoff_length:	#only 20 characters
				element = element[:cutoff_length]
			if replaceNoneElemWithEmptyStr and element is None:
				element = ''
			gene_desc_ls.append(element)
		return gene_desc_ls
	
	@classmethod
	def plot_one_gene(cls, ax, gene_id, param_data, base_y_value=1):
		"""
		2008-12-16
			handle drawing on a genome-wide plot
			param_data has chr_id2cumu_size, chr_id2size, chr_gap
			add artist_obj_id2artist_gene_id_ls & gene_id2artist_object_id param_data
		2008-10-02
			plot all different forms of a gene out.
		2008-10-01
			gene_model is changed. now it's from GenomeDB.py
		2008-09-24
			draw a single gene on the canvas, adapted from GenomeBrowser.plot_one_gene()
		"""
		gene_model = param_data.gene_id2model.get(gene_id)
		chr_id2cumu_size = getattr(param_data, 'chr_id2cumu_size', None)
		chr_id2size = getattr(param_data, 'chr_id2size', None)
		chr_gap = getattr(param_data, 'chr_gap', None)
		artist_obj_id2artist_gene_id_ls = getattr(param_data, 'artist_obj_id2artist_gene_id_ls', None)
		gene_id2artist_object_id = getattr(param_data, 'gene_id2artist_object_id', None)
		#2008-12-16
		if artist_obj_id2artist_gene_id_ls is not None and gene_id2artist_object_id is not None:
			if gene_id in gene_id2artist_object_id:	#drawn already
				return
		if gene_model:
			if len(gene_model.gene_commentaries)==0:
				gene_commentaries = [gene_model]	#fake one here
			else:
				gene_commentaries = gene_model.gene_commentaries
			for gene_commentary in gene_commentaries:	#multiple commentary
				gene_desc_ls = cls.returnGeneDescLs(cls.gene_desc_names, gene_model, gene_commentary)
				param_data.matrix_of_gene_descriptions.append(gene_desc_ls)
				
				if getattr(gene_commentary, 'box_ls', None):
					box_ls = gene_commentary.box_ls
				elif gene_commentary.start and gene_commentary.stop:	#no box_ls, just use start, stop
					box_ls = [(gene_commentary.start, gene_commentary.stop, 'exon')]
				elif gene_model.start and gene_model.stop:	#use gene_model's coordinate
					box_ls = [(gene_model.start, gene_model.stop, 'exon')]
				else:	#no coordinate. skip
					sys.stderr.write("Warning: gene %s has either no start or no stop.\n"%(gene_commentary.gene_id))
					continue
				if param_data.candidate_gene_set and gene_id in param_data.candidate_gene_set:
					gene_symbol_color = 'b'
				else:
					gene_symbol_color = 'k'
				
				#2008-12-16 draw genes on a genome-wide plot
				if chr_id2cumu_size and chr_id2size and chr_gap:
					this_chr_starting_pos_on_plot = chr_id2cumu_size[gene_model.chromosome]-chr_id2size[gene_model.chromosome]-chr_gap
					for i in range(len(box_ls)):
						box = box_ls[i]
						box_ls[i] = [box[0]+this_chr_starting_pos_on_plot, box[1]+this_chr_starting_pos_on_plot]+list(box[2:])
				else:
					this_chr_starting_pos_on_plot = 0
				y_value = base_y_value+param_data.no_of_genes_drawn%param_data.gene_position_cycle	#cycling through the y position to avoid clogging
				
				g_artist = ExonIntronCollection(box_ls, y=y_value, strand=gene_model.strand, width=param_data.gene_width, alpha=0.7, \
											picker=True, linewidths=0.7, box_line_widths=0.3, rotate_xy=param_data.rotate_xy)
				ax.add_artist(g_artist)	# 2010-4-12 new ExonIntronCollection is a pathPatch. add_patch and add_artist are same
				
				#2008-12-16 for GenomeBrowser.py to keep track of objects
				if artist_obj_id2artist_gene_id_ls is not None and gene_id2artist_object_id is not None:
					artist_obj_id = id(g_artist)
					param_data.artist_obj_id2artist_gene_id_ls[artist_obj_id] = [g_artist, gene_id]
					param_data.gene_id2artist_object_id[gene_id] = artist_obj_id
				
				text_start_pos = box_ls[-1][1]
				if this_chr_starting_pos_on_plot:	#2008-12-16
					text_start_pos += this_chr_starting_pos_on_plot
				#mid_point = (c_start_ls[0]+c_end_ls[-1])/2.
				if param_data.label_gene or gene_symbol_color=='b':	#specified to label or it's candidate gene
					if param_data.rotate_xy:
						ax.text(y_value, text_start_pos+param_data.gene_box_text_gap, gene_model.gene_symbol, size=4, \
							color=gene_symbol_color, alpha=0.8, verticalalignment='center', rotation='vertical')
					else:
						ax.text(text_start_pos+param_data.gene_box_text_gap, y_value, gene_model.gene_symbol, size=4, \
							color=gene_symbol_color, alpha=0.8, verticalalignment='center')
				param_data.no_of_genes_drawn += 1
	
	@classmethod
	def drawGeneModel(cls, ax, region=None, gene_annotation=None, candidate_gene_set=None, gene_width=1.0, \
					gene_position_cycle=4, base_y_value=1, gene_box_text_gap=100, label_gene=0, rotate_xy=False,\
					chr_id2cumu_size=None, chr_id2size=None, chr_gap=None, artist_obj_id2artist_gene_id_ls=None,\
					gene_id2artist_object_id=None, drawGeneOnTheBoundary=True):
		"""
		2012.11.14 rename snps_within_this_region to region (3 attributes, chromosome, start, stop)
		2009-4-30
			handle chr,pos,offset in snps_within_this_region.chr_pos_ls by ignoring the offset bit
		2008-12-22
			add argument drawGeneOnTheBoundary, which controls whether to draw genes that are not fully within this region.
				genes that sit on the boundary of snps_within_this_region.
			GenomeBrowser.py set it to False because later adding text to these genes would corrupt the running program.
		2008-12-16
			add options: chr_id2cumu_size, chr_id2size, chr_gap
		2008-10-01
		2008-09-24
		"""
		sys.stderr.write("\t Drawing gene model  ...")
		#2012.3.8 use the user-indicated span, rather than boundary of polymorphism loci
		left_chr, left_pos = region.chromosome, region.start
		right_chr, right_pos = region.chromosome, region.stop
		left_chr = str(left_chr)
		right_chr = str(right_chr)
		param_data = PassingData(gene_id2model=gene_annotation.gene_id2model, candidate_gene_set=candidate_gene_set, \
								gene_width=gene_width, gene_position_cycle=gene_position_cycle, no_of_genes_drawn=0, \
								gene_box_text_gap=gene_box_text_gap, matrix_of_gene_descriptions = [], label_gene=label_gene,\
								rotate_xy=rotate_xy, chr_id2cumu_size=chr_id2cumu_size, chr_id2size=chr_id2size, chr_gap=chr_gap,\
								artist_obj_id2artist_gene_id_ls=artist_obj_id2artist_gene_id_ls, gene_id2artist_object_id=gene_id2artist_object_id)
		for gene_id in gene_annotation.chr_id2gene_id_ls[left_chr]:
			gene_model = gene_annotation.gene_id2model[gene_id]
			if drawGeneOnTheBoundary:
				gene_position = gene_model.stop
			else:
				gene_position = gene_model.start
			if gene_model.start!=None and gene_model.stop!=None and gene_position>left_pos:
				if left_chr==right_chr:	#same chromosome
					if gene_model.start>right_pos:	#totally out of range, skip it
						continue
				cls.plot_one_gene(ax, gene_id, param_data, base_y_value=base_y_value)
		if left_chr!=right_chr:
			for gene_id in gene_annotation.chr_id2gene_id_ls[right_chr]:
				gene_model = gene_annotation.gene_id2model[gene_id]
				if drawGeneOnTheBoundary:
					gene_position = gene_model.start
				else:
					gene_position = gene_model.stop
				if gene_model.start!=None and gene_model.stop!=None and gene_position<right_pos:
					cls.plot_one_gene(ax, gene_id, param_data, base_y_value=base_y_value)
		sys.stderr.write("Done.\n")
		return param_data
	
	@classmethod
	def drawLD(cls, axe_gene_model, ax2, snps_within_this_region, LD_info, gene_model_min_y,\
				gene_model_max_y, which_LD_statistic=1, markSNPInGenes=True):
		"""
		2010-4-23
			add argument markSNPInGenes, whether to draw vertical lines in the axe_gene_model
		2009-4-30
			handle chr,pos,offset in snps_within_this_region.chr_pos_ls by ignoring the offset bit
		2008-09-28
			represent r2 by HSL color, rather than a grayscale intensity. matplotlib doesn't support hsl notation (i'm not aware).
			Use ImageColor.getrgb to convert hsl representation to RGB format.
		2008-09-24
			draw LD in the bottom axe
		"""
		sys.stderr.write("\t Drawing LD info  ...")
		no_of_snps = len(snps_within_this_region.chr_pos_ls)
		left_chr, left_pos = snps_within_this_region.chr_pos_ls[0][:2]
		right_chr, right_pos = snps_within_this_region.chr_pos_ls[-1][:2]
		#axe_pvalue.hlines(y_value, left_pos, right_pos, linewidth=0.3)
		for i in range(no_of_snps):
			chr_pos1 = snps_within_this_region.chr_pos_ls[i]
			if markSNPInGenes:	#2010-4-23
				axe_gene_model.vlines(chr_pos1[1], gene_model_min_y, gene_model_max_y, linestyle='dashed', alpha=0.3, linewidth=0.3)
			if not LD_info:	#skip drawng LD boxes
				continue
			for j in range(i+1, no_of_snps):
				chr_pos2 = snps_within_this_region.chr_pos_ls[j]
				if chr_pos1<chr_pos2:
					snp_pair = (chr_pos1[0], chr_pos1[1], chr_pos2[0], chr_pos2[1])
				else:
					snp_pair = (chr_pos2[0], chr_pos2[1], chr_pos1[0], chr_pos1[1])
				if snp_pair in LD_info.snp_pair2r2:
					LD_stat = getattr(LD_info.snp_pair2r2[snp_pair], LD_statistic.get_name(which_LD_statistic), None)
					LD_stat = abs(LD_stat)	#D_prime, D need abs()
					s11, s12 = snps_within_this_region.chr_pos2adjacent_window[chr_pos1]
					s21, s22 = snps_within_this_region.chr_pos2adjacent_window[chr_pos2]
					#draw a plot to understand this. a parallegram below SNPs, it's the intersection of two adjacent windows 
					x1 = (s12+s21)/2.
					y1 = (s12-s21)/2.
					x2 = (s12+s22)/2.
					y2 = (s12-s22)/2.
					x3 = (s11+s22)/2.
					y3 = (s11-s22)/2.
					x4 = (s11+s21)/2.
					y4 = (s11-s21)/2.
					xs = [x1, x2, x3, x4]
					ys = [y1, y2, y3, y4]
					fc = cls.r2ToRGBColor(LD_stat)
					poly = Polygon(zip(xs, ys), facecolor=fc, linewidth=0)
					ax2.add_patch(poly)
		sys.stderr.write("Done.\n")
	
	@classmethod
	def r2ToRGBColor(cls, r2):
		"""
		2008-09-29
			convert r2 to matplotlib RGB color
		"""
		hue_value = int(round((1-r2)*255))	#r2 is [0,1], 255 is the maximum hue value. lower r2, higher hue. r2=0 -> blue, r2=1 -> red
		fc = ImageColor.getrgb('hsl(%s'%hue_value+',100%,50%)')	#"hsl(hue, saturation%, lightness%)" where hue is the colour given as an
		# angle between 0 and 360 (red=0, green=120, blue=240),
		#saturation is a value between 0% and 100% (gray=0%, full color=100%), and lightness is a value between 0% and 100% (black=0%, normal=50%, white=100%).
		fc = [color_value/255. for color_value in fc]	#matplotlib accepts rgb in [0-1] range
		return fc
	
	@classmethod
	def drawLDLegend(cls, ax, which_LD_statistic):
		"""
		2008-09-29
			left half is 10 color patchs, right half is r2 from 0 to 1
		"""
		sys.stderr.write("\t Drawing LD legend  ...")
		xs = [0,0,0.5,0.5]	#X-axis value for the 4 points of the rectangle starting from upper left corner. 
		ys = numpy.array([0.1,0,0,0.1])
		for i in range(10):
			r2 = ys[0]	#upper bound of ys correspond to r2
			fc = cls.r2ToRGBColor(r2)
			poly = Polygon(zip(xs, ys), facecolor=fc, linewidth=0)
			ax.add_patch(poly)
			if i%2==0:	#every other band
				ax.text(0.6, ys[0], '%.1f'%r2, horizontalalignment ='left', verticalalignment='top', size=4)
			ys += 0.1	#increase y-axis
		ax.text(0.6, 1.1, r"$%s=1.0$"%LD_statistic.get_label(which_LD_statistic), horizontalalignment ='left', verticalalignment='top', size=4)
		sys.stderr.write("Done.\n")
	
	@classmethod
	def construct_chr_pos2index_forSNPData(cls, snpData, snp_info=None):
		"""
		2012.3.7
			handle the situation that snpData.col_id_ls is comprised of db-ids, not chr_pos.
			add argument snp_info
			
			snp_info = SNPInfo()
			snp_info.chr_pos_ls = chr_pos_ls
			snp_info.chr_pos2index = chr_pos2index
			snp_info.snps_id2index = snps_id2index
			snp_info.data_ls = data_ls
		#2008-12-05 fake a snp_info for findSNPsInRegion
			
		"""
		chr_pos2index = {}
		for i in range(len(snpData.col_id_ls)):
			db_id = snpData.col_id_ls[i]
			chr_pos_index = snp_info.snps_id2index.get(db_id)
			if chr_pos_index is not None:
				snp_db_entry = snp_info.data_ls[chr_pos_index]
				chr_pos = (snp_db_entry.chromosome, snp_db_entry.position, snp_db_entry.end_position)
				if len(chr_pos)==3 and chr_pos[2]!='0':	#2009-3-27 ignore insertion right now
					continue
				chr_pos = tuple(map(int, chr_pos))
				chr_pos2index[chr_pos] = i
		snpData.chr_pos2index = chr_pos2index
	
	@classmethod
	def drawRegionAroundThisSNP(cls, phenotype_method_id=None, this_snp=None, candidate_gene_set=None, gene_annotation=None,
							snp_info=None, \
							analysis_method_id2gwr=None, \
							LD_info=None, output_dir=None, which_LD_statistic=None, snp_region=None, min_distance=40000, \
							list_type_id_list=None, label_gene=0,\
							draw_LD_relative_to_center_SNP=0, commit=0, snpData=None, phenData=None, ecotype_info=None,\
							snpData_before_impute=None, snp_matrix_data_type=1, call_method_id=None,\
							drawMap=False, drawStrainPCA=True, markSNPInGenes=True, drawGrid=False, phenotypeDrawType=1):
		"""
		2012.3.19
			if output_dir is not a directory, it's treated as a filename prefix. newly minted filename is attached to it.
			add this_snp.fileNamePrefix into output filename if it's not None.
		2010-4-26
			add argument phenotypeDrawType
				1: the length of each accession's stick varies according to phenotype magnitude but same color
				2: the length of each stick is same but color varies according to the phenotype hsv bar
		2010-4-23
			add argument markSNPInGenes, drawGrid
		2010-4-19
			add argument drawStrainPCA, and adjust axe width if this part is not drawn.
		2010-4-13
			add argument drawMap. controlling whether to draw a geographic map or not.
		2010-1-22 if call_method_id is not specified, try cls.call_method_id. if still no, just use 'Unknown'.
		2009-5-30
			add call_method_id into the output figure filename
		2009-4-30
			deal with argument snp_matrix_data_type =4 (arbitrary non-diallelic SNP matrix)
			could handle both (chr,pos) and (chr,pos,offset) SNP representation
			skip drawing gene models if no SNPs in the region at all.
			return None if no SNPs are found in the region.
		2009-3-23
			add arguments snp_matrix_data_type to allow CNV or CNV amplitude to fill the SNP matrix
		2008-12-01
			add option snpData_before_impute
		2008-11-30
			add code to replace axe_LD with axe_strain_pca, axe_snp_matrix, axe_map to demonstrate the haplotype structure,
				phenotype and geographic source of strains.
		2008-10-24
			handle option commit to return png_data & svg_data
		2008-10-01
			remove the frame of axe_pvalue and add a grid to axe_pvalue
			leave axe_gene_model's xticks there as otherwise axe_pvalue's xticks will go with it as they share xticks.
		2008-10-01
			draw gene models on a separate axe,
			add a twinx axe to draw LD w.r.t the center SNP
			if output_dir is not a directory, it's treated as a filename.
		2008-09-24
		"""
		sys.stderr.write("Drawing region ... \n")
		phenotype = Stock_250kDB.PhenotypeMethod.get(phenotype_method_id)
		
		#list_type = Stock_250kDB.GeneListType.get(list_type_id)
		if call_method_id is None:
			call_method_id = getattr(cls, 'call_method_id', "Unknown")
		fname_basename = 'call_%s_snp_%s_%s_%s_phenotype_%s_%s'%\
						(call_method_id, this_snp.chromosome, this_snp.position, this_snp.stop, \
						phenotype.id, phenotype.getProperShortName())
		if hasattr(this_snp, 'fileNamePrefix'):	#2012.3.19
			maxPrefixLength=150
			if len(this_snp.fileNamePrefix)>maxPrefixLength:	#cut the filename prefix shorter than maxPrefixLength
				# 2012.3.28 length of entire filename has to be <=260 
				fileNamePrefix = this_snp.fileNamePrefix[:maxPrefixLength]
			else:
				fileNamePrefix = this_snp.fileNamePrefix
			fname_basename = "%s_call%s_phenotype_%s_%s"%(fileNamePrefix, call_method_id, phenotype.id, phenotype.getProperShortName())
		fname_basename = fname_basename.replace('/', '_')	#to avoid file system error
		if os.path.isdir(output_dir):
			output_fname_prefix = os.path.join(output_dir, fname_basename)
		else:	#output_dir is a filename prefix
			output_fname_prefix = "%s_%s"%(output_dir, fname_basename)
		# 2010-9-17 create the directory if it doesn't exist.
		output_dir = os.path.split(output_fname_prefix)[0]
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		
		if snp_region:
			snps_within_this_region = snp_region
		elif getattr(this_snp, 'stop', None):
			snps_within_this_region = cls.findSNPsInRegion(snp_info, chromosome=this_snp.chromosome, start=this_snp.position, \
												stop=this_snp.stop)
			snps_within_this_region_snpData = snps_within_this_region
			#snps_within_this_region_snpData = cls.findSNPsInRegion(snpData, chromosome=this_snp.chromosome, start=this_snp.position, \
			#									stop=this_snp.stop)
		else:
			snps_within_this_region = cls.findSNPsInRegion(snp_info, chromosome=this_snp.chromosome, start=this_snp.position-min_distance, \
														stop=this_snp.stop+min_distance, center_snp_position=this_snp.position)
			snps_within_this_region_snpData = snps_within_this_region
		if len(snps_within_this_region.chr_pos_ls)==0 and len(snps_within_this_region_snpData.chr_pos_ls)==0:
			return None
		
		pylab.clf()
		#fig = pylab.figure()
		axe_y_offset1 = 0.05	#y_offset for axe_LD, axe_strain_pca, axe_phenotype, axe_map
		axe_y_offset2 = 0.6
		axe_height1 = axe_y_offset2 - axe_y_offset1	#height of axe_LD or axe_snp_matrix
		axe_y_offset3 = 0.7
		axe_height2 = axe_y_offset3 - axe_y_offset2	#height of axe_gene_model
		axe_y_offset4 = 0.95
		axe_height3 = axe_y_offset4 - axe_y_offset3	#height of axe_pvalue
		
		axe_x_offset1 = 0.02	#
		if drawStrainPCA:	#2010-4-19
			axe_x_offset2 = axe_x_offset1 + 0.2
		else:
			axe_x_offset2 = axe_x_offset1
		axe_width1 = axe_x_offset2 - axe_x_offset1	#width of axe_strain_pca
		axe_x_offset3 = 0.77
		axe_width2 = axe_x_offset3 - axe_x_offset2	#width of axe_pvalue, axe_LD, or axe_snp_matrix
		axe_x_offset4 = 0.79 
		axe_width3 = axe_x_offset4 - axe_x_offset3	#width of axe_phenotype
		axe_x_offset5 = 0.99
		axe_width4 = axe_x_offset5 - axe_x_offset4	#width of axe_map, axe_map_frame
		no_of_axes_drawn = 0
		
		axe_pvalue = pylab.axes([axe_x_offset2, axe_y_offset3, axe_width2, axe_height3], frameon=False)	#left gap, bottom gap, width, height, axes for pvalue, gene models
		axe_pvalue.grid(True, alpha=0.3)
		axe_pvalue.set_xticklabels([])	#remove xtick labels on axe_pvalue because axe_LD's xtick labels cover this.
		axe_LD_center_SNP = pylab.twinx()	#axes for LD with center SNP, copy axe_pvalue's
		axe_LD_center_SNP.set_xticklabels([])
		axe_gene_model = pylab.axes([axe_x_offset2, axe_y_offset2, axe_width2, axe_height2], frameon=False, sharex=axe_pvalue)
		#axe_gene_model.set_xticks([])	#this will set axe_pvalue's xticks off as well because the x-axis is shared.
		axe_gene_model.set_yticks([])
		
		# use actual loci data to restrict the region to be shown
		#snp_region_tup = [snps_within_this_region_snpData.chr_pos_ls[0][0], snps_within_this_region_snpData.chr_pos_ls[0][1],\
		#				snps_within_this_region_snpData.chr_pos_ls[-1][0], snps_within_this_region_snpData.chr_pos_ls[-1][1]]
		# 2012.3.8 region decided by snps_within_this_region_snpData
		snp_region_tup = [snps_within_this_region_snpData.chromosome, snps_within_this_region_snpData.start,\
						snps_within_this_region_snpData.chromosome, snps_within_this_region_snpData.stop]
		axe_snp_matrix_margin = abs(snp_region_tup[3]-snp_region_tup[1])/15.	#offset to push strain labels on even rows further right
		if LD_info:
			axe_LD = pylab.axes([axe_x_offset2, axe_y_offset1, axe_width2, axe_height1], frameon=False)	#axes for LD
			axe_LD_legend = pylab.axes([axe_x_offset3-0.1, axe_y_offset1+0.03, 0.1, 0.13], frameon=False)	#axes for the legend of LD
			axe_LD_legend.set_xticks([])
			axe_LD_legend.set_yticks([])
			axe_to_put_pvalue_legend = axe_LD
			legend_loc = 'lower left'
			axe_pvalue_xlim = [snp_region_tup[1]-axe_snp_matrix_margin, snp_region_tup[3]+axe_snp_matrix_margin]
		elif snpData:
			phenotype_col_index = cls.findOutWhichPhenotypeColumn(phenData, set([phenotype_method_id]))[0]
			genome_wide_result = analysis_method_id2gwr.get(1)
			if not genome_wide_result:
				sys.stderr.write("No genome association results for phenotype_method_id=%s, analysis_method_id=%s. Take a random one out of analysis_method_id2gwr.\n"%\
							(phenotype_method_id, 1))
				genome_wide_result = analysis_method_id2gwr.values()[0]	#take random genome_wide_result
			if snp_matrix_data_type==1:
				chr_pos_ls = []
			elif snp_matrix_data_type==3:
				chr_pos_ls = snpData.chr_pos2index.keys()
				#2008-12-08 for CNV probes. use snpData.chr_pos2index.keys() to locate top_snp_data because here snpData doesn't match genome_wide_result.
			else:	#2012.3.27 for other types, genome_wide_result matches snpData. 
				chr_pos_ls = []
			top_snp_data = cls.getTopSNPData(genome_wide_result, None, snp_region_tup, chr_pos_ls=chr_pos_ls)
			if len(top_snp_data.snp_id_ls)==0:	#2012.3.23 no plot
				return None
			if snp_matrix_data_type==3:	#2009-3-23 for CNV amplitude data, don't convert alleles into binary 0/1=major/minor form and use allele/amplitude to determine alpha
				need_convert_alleles2binary = False
				useAlleleToDetermineAlpha = True
			elif snp_matrix_data_type==4:
				#2009-3-27, for arbitrary non-diallelic SNP matrix
				need_convert_alleles2binary = False
				useAlleleToDetermineAlpha = False
			else:
				need_convert_alleles2binary = True
				useAlleleToDetermineAlpha = False
			
			subSNPData = cls.getSubStrainSNPMatrix(snpData, phenData, phenotype_method_id=phenotype_method_id, \
												phenotype_col_index=phenotype_col_index, \
												snp_id_ls=top_snp_data.snp_id_ls, \
												need_convert_alleles2binary=need_convert_alleles2binary)	#2009-3-23 last argument is for CNV intensity matrix
			snp_value2color = None
			if snp_matrix_data_type==4:
				##2009-3-27 it's for SNP matrix inferred from raw sequences, might have >2 alleles, heterozygous calls, deletions etc.
				from DrawSNPMatrix import DrawSNPMatrix
				subSNPData.data_matrix = DrawSNPMatrix.transformMatrixIntoTwoAllelesAndHetero(subSNPData.data_matrix)
				snp_value2color = cls.snp_value2five_color
			
			#the two offsets below decides where the label of strains/snps should start in axe_snp_matrix
			last_chr_pos = snps_within_this_region_snpData.chr_pos_ls[-1]
			strain_id_label_x_offset=snps_within_this_region_snpData.chr_pos2adjacent_window[last_chr_pos][1]	#right next to the rightmost SNP
			snp_id_label_y_offset=0.95
			
			StrainID2PCAPosInfo = cls.getStrainID2PCAPosInfo(subSNPData, pca_range=[0,1], snp_id_label_y_offset=snp_id_label_y_offset)
			
			#fake one SNPID2PCAPosInfo only for drawSNPMtrix()
			SNPID2PCAPosInfo = PassingData(step=None, snp_id2img_x_pos={})
			for locus_id, adjacent_window in snps_within_this_region_snpData.locus_id2adjacent_window.iteritems():
				#2012.3.8 use locus_id, rather than chr_pos,
				#chr_pos = map(str, chr_pos)
				#snp_id = '_'.join(chr_pos)
				SNPID2PCAPosInfo.snp_id2img_x_pos[str(locus_id)] = adjacent_window
			
			phenotype_cmap = mpl.cm.jet
			max_phenotype = numpy.nanmax(phenData.data_matrix[:,phenotype_col_index])
			min_phenotype = numpy.nanmin(phenData.data_matrix[:,phenotype_col_index])
			phenotype_gap = max_phenotype - min_phenotype
			phenotype_jitter = phenotype_gap/10.
			phenotype_norm = mpl.colors.Normalize(vmin=min_phenotype-phenotype_jitter, vmax=max_phenotype+phenotype_jitter)
			axe_map_phenotype_legend = pylab.axes([axe_x_offset4+0.02, axe_y_offset1, axe_width4-0.02, axe_height1/10.], frameon=False)
			cb = mpl.colorbar.ColorbarBase(axe_map_phenotype_legend, cmap=phenotype_cmap,
										norm=phenotype_norm,
										orientation='horizontal')
			cb.set_label('Phenotype Legend On the Map')
			
			axe_snp_matrix = pylab.axes([axe_x_offset2, axe_y_offset1, axe_width2, axe_height1], frameon=False)
			#axe_snp_matrix.set_xticks([])
			axe_snp_matrix.set_yticks([])
			cls.drawSNPMtrix(axe_snp_matrix, subSNPData, top_snp_data, StrainID2PCAPosInfo, SNPID2PCAPosInfo, \
							ecotype_info, strain_id_label_x_offset, snp_id_label_y_offset, \
							strain_id_label_x_offset_extra=axe_snp_matrix_margin, \
							draw_snp_id_label=False, snpData_before_impute=snpData_before_impute, \
							strain_snp_label_font_size=1, \
							useAlleleToDetermineAlpha=useAlleleToDetermineAlpha, \
							snp_value2color=snp_value2color, drawGrid=drawGrid)	#2008-11-14 turn draw_snp_id_label off
			#axe_snp_matrix.set_xlim([0,1])
			#axe_snp_matrix.set_ylim([0,1])
			no_of_axes_drawn += 1
			#pylab.savefig('%s_%s.png'%(self.output_fname_prefix, no_of_axes_drawn), dpi=400)
			
			if drawStrainPCA:	#2010-4-19
				axe_strain_map = None	#no strain map
				axe_strain_pca = pylab.axes([axe_x_offset1, axe_y_offset1, axe_width1, axe_height1], frameon=False, sharey=axe_snp_matrix)
				axe_strain_map_pca_cover = None	#not used.
				axe_strain_pca_xlim = [-0.05,1.05]
				axe_strain_pca_ylim = [0, 1]
				axe_strain_pca.set_xlim(axe_strain_pca_xlim)
				axe_strain_pca.set_ylim(axe_strain_pca_ylim)
				
				axe_strain_pca.grid(True, alpha=0.3)
				axe_strain_pca.set_xticks([])
				axe_strain_pca.set_yticks([])
				axe_strain_pca_legend =None
				
				cls.drawStrainPCA(axe_strain_pca, axe_strain_map, axe_strain_map_pca_cover, axe_strain_pca_legend, StrainID2PCAPosInfo, \
								ecotype_info, phenData, \
							phenotype_col_index, phenotype_cmap, phenotype_norm, rightmost_x_value=axe_strain_pca_xlim[1],\
							country_order_name='', strain_color_type=2, draw_axe_strain_map=False)
				
				axe_strain_pca.set_xlim(axe_strain_pca_xlim)
				axe_strain_pca.set_ylim(axe_strain_pca_ylim)
				no_of_axes_drawn += 1
				if getattr(cls, "debug", False):
					pylab.savefig('%s_%s.png'%(output_fname_prefix, no_of_axes_drawn), dpi=400)
			
			if drawMap:
				#mark ecotypes on the map colored according to phenotype
				axe_map = pylab.axes([axe_x_offset4, axe_y_offset1, axe_width4, axe_height1], frameon=False)
				#axe_map_frame is used to connect strains from axe_phenotype to dot on the axe_map (another axe due to reasons stated in drawMap())
				axe_map_frame = pylab.axes([axe_x_offset4, axe_y_offset1, axe_width4, axe_height1], frameon=False, sharey=axe_snp_matrix)
				axe_map_frame.set_xticks([])
				axe_map_frame.set_yticks([])
				cls.drawMap(axe_map_frame, axe_map, StrainID2PCAPosInfo, phenData, phenotype_col_index, phenotype_method_id, \
						ecotype_info, phenotype_cmap, phenotype_norm)
			
				#axe_map.set_ylim([0,1])
				no_of_axes_drawn += 1
				if getattr(cls, "debug", False):
					pylab.savefig('%s_%s.png'%(output_fname_prefix, no_of_axes_drawn), dpi=400)
			
			axe_phenotype = pylab.axes([axe_x_offset3, axe_y_offset1, axe_width3, axe_height1], frameon=False, sharey=axe_snp_matrix)
			axe_phenotype.set_yticks([])
			axe_phenotype.set_xticklabels([])	#no tick labels (axe_map_phenotype_legend has it already)
			if phenotypeDrawType==1:	#2010-4-26
				cls.drawPhenotype(axe_phenotype, StrainID2PCAPosInfo, phenData, phenotype_col_index, \
								phenotype_method_id, ecotype_info)
			else:
				cls.drawPhenotype(axe_phenotype, StrainID2PCAPosInfo, phenData, phenotype_col_index, phenotype_method_id, \
								ecotype_info,\
							phenotype_cmap=phenotype_cmap, phenotype_norm=phenotype_norm)	# 2010-2-5 add phenotype_cmap & phenotype_norm
			no_of_axes_drawn += 1
			axe_phenotype.set_ylim([0,1])
			
			
			axe_snp_matrix.set_ylim([0,1])	#without this, ylim of all 3 axes are set to [0,0.9] because axe_map automatically adjust to 0-0.9
			#pylab.savefig('%s_%s.png'%(self.output_fname_prefix, no_of_axes_drawn), dpi=400)
			
			axe_to_put_pvalue_legend = axe_pvalue	#axe_LD is gone. put legend into axe_pvalue itself.
			legend_loc = 'upper right'
			axe_LD = None
			axe_LD_legend = None
			
			axe_pvalue_xlim = [snp_region_tup[1]-axe_snp_matrix_margin, snp_region_tup[3]+axe_snp_matrix_margin*2]
			
		fig_title = 'SNP chr %s. pos %s.'%(this_snp.chromosome, this_snp.position)
		if getattr(this_snp, 'stop', None):
			fig_title += ' - %s. '%this_snp.stop
		fig_title += "Phenotype %s (id=%s)."%(phenotype.short_name, phenotype.id)
		axe_pvalue.title.set_text(fig_title)	#main title using this snp.
		
		#2009-5-29 temporarily snps_within_this_region is replaced by snps_within_this_region_snpData below
		cls.drawPvalue(axe_pvalue, axe_to_put_pvalue_legend, axe_LD_center_SNP, snps_within_this_region_snpData, \
					analysis_method_id2gwr, LD_info, \
					which_LD_statistic, draw_LD_relative_to_center_SNP=draw_LD_relative_to_center_SNP, legend_loc=legend_loc)
		gene_position_cycle = 5
		base_y_value = 1
		gene_width=0.8
		gene_box_text_gap = min_distance*2*0.005
		
		skip_gene_model = False
		if len(snps_within_this_region.chr_pos_ls)>0:
			_snps_within_this_region = snps_within_this_region
		elif len(snps_within_this_region_snpData.chr_pos_ls)>0:
			_snps_within_this_region = snps_within_this_region_snpData
		else:
			skip_gene_model = True
		if not skip_gene_model:
			return_data = cls.drawGeneModel(axe_gene_model, _snps_within_this_region, gene_annotation, \
											candidate_gene_set, gene_width=gene_width, gene_position_cycle=gene_position_cycle, \
											base_y_value=base_y_value, gene_box_text_gap=gene_box_text_gap,\
											label_gene=label_gene)
		matrix_of_gene_descriptions = return_data.matrix_of_gene_descriptions
		gene_model_min_y = base_y_value-gene_width
		gene_model_max_y = gene_position_cycle + base_y_value -1 + gene_width	#"-1" because genes never sit on y=gene_position_cycle + base_y_value
		
		if not skip_gene_model:
			cls.drawLD(axe_gene_model, axe_LD, _snps_within_this_region, LD_info, gene_model_min_y=gene_model_min_y,\
					gene_model_max_y=gene_model_max_y, which_LD_statistic=which_LD_statistic, \
					markSNPInGenes=markSNPInGenes)
		if LD_info:
			cls.drawLDLegend(axe_LD_legend, which_LD_statistic)
		#adjust x, y limits and etc
		axe_pvalue.set_xlim(axe_pvalue_xlim)
		axe_pvalue_ylim = axe_pvalue.get_ylim()
		axe_pvalue.set_ylim((0, axe_pvalue_ylim[1]))	#set axe_pvalue to 0 to sit right above axe_gene_model
		axe_gene_model.set_ylim((gene_model_min_y, gene_model_max_y))	#LD panel right under gene models
		
		if LD_info:
			axe_LD.set_xlim(axe_pvalue.get_xlim())	#make the axe_LD and axe_pvalue within the same X range
			axe_LD_x_span = (axe_LD.get_xlim()[1]-axe_LD.get_xlim()[0])
			axe_LD.set_ylim((-axe_LD_x_span/2., 0))	#has to force here, don't know why. otherwise it's (0,1)
			axe_LD.set_yticks([])	#remove all Y ticks on LD plot
		elif snpData:
			axe_snp_matrix.set_xlim(axe_pvalue.get_xlim())
		
		png_data = None
		svg_data = None
		png_output_fname = None
		if len(snps_within_this_region.chr_pos_ls)>0:
			distance = abs(snps_within_this_region.chr_pos_ls[-1][1] - snps_within_this_region.chr_pos_ls[0][1])
		elif len(snps_within_this_region_snpData.chr_pos_ls)>0:
			distance = abs(snps_within_this_region_snpData.chr_pos_ls[-1][1] - snps_within_this_region_snpData.chr_pos_ls[0][1])
		else:
			distance = 0
		if commit:	#2008-10-24
			png_data = StringIO.StringIO()
			svg_data = StringIO.StringIO()
			pylab.savefig(png_data, format='png', dpi=600)
			if distance<=5000:		#save the svg format if less than 80kb
				pylab.savefig(svg_data, format='svg', dpi=300)
		else:
			png_output_fname = '%s.png'%output_fname_prefix
			pylab.savefig(png_output_fname, dpi=600)
			
			#if distance<=5000:		#save the svg format if less than 80kb
			#	pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
			#pylab.savefig('%s.pdf'%output_fname_prefix, dpi=300)	# 2010-2-7 swedish letters can't be saved in pdf.
		if getattr(cls, "debug", False):
			pylab.show()
		del axe_pvalue, axe_LD_center_SNP, axe_gene_model, axe_LD, axe_LD_legend
		sys.stderr.write("Done.\n")
		after_plot_data = PassingData(png_output_fname=png_output_fname, matrix_of_gene_descriptions=matrix_of_gene_descriptions, \
									png_data=png_data,\
									svg_data=svg_data,\
									snps_within_this_region=snps_within_this_region)
		return after_plot_data
	
	def generate_params(self, param_obj):
		"""
		2008-09-24
			copied from a version of MpiGeneListRankTest.py
		"""
		sys.stderr.write("Generating parameters ...")
		i = 0
		block_size = 5000
		query = Stock_250kDB.ResultsByGene.query
		if param_obj.call_method_id!=0:
			query = query.filter(Stock_250kDB.ResultsByGene.results_method.has(call_method_id=param_obj.call_method_id))
		if param_obj.analysis_method_id!=0 and param_obj.analysis_method_id is not None:
			query = query.filter(Stock_250kDB.ResultsByGene.results_method.has(analysis_method_id=param_obj.analysis_method_id))
		query = query.filter_by(min_distance=param_obj.min_distance).filter_by(get_closest=param_obj.get_closest)
		rows = query.offset(i).limit(block_size)
		results_id_ls = []
		while rows.count()!=0:
			for row in rows:
				results_id_ls.append(row.id)
				i += 1
			rows = query.offset(i).limit(block_size)
		
		sys.stderr.write("%s results. "%(len(results_id_ls)))
		return results_id_ls
	
	def getSNPsFromInputFile(self, input_fname=None, snp_info=None, db_250k=None):
		"""
		2010-2-5
			ignore empty rows
		2008-11-30
			phenotype_id_index could also be under column name 'phenotype'
		2008-10-27
			add code to deal with suzi's input file
			suzi's input '/Network/Data/250k/tmp-Suzi/sig_marker_ft_files/sig_marker_1_LD.csv' use 'phen' and the phenotype identifier, '1_LD', needs split to get id.
		2008-10-04
			input could contain one more column, 'stop', indicating it is a region, not a single SNP. 'start' becomes alternative name for 'position'.
		2008-10-02
			use GeneSNPPhenotypeAssoData to hold input data
		2008-09-30
			read input from a file. check module doc for format.
		"""
		sys.stderr.write("Getting input from %s ..."%input_fname)
		input_data = GeneSNPPhenotypeAssoData(db_250k=db_250k)
		delimiter = figureOutDelimiter(input_fname)
		reader = csv.reader(open(input_fname), delimiter=delimiter)
		col_name2index = getColName2IndexFromHeader(reader.next())
		counter = 0
		for row in reader:
			if not row:	# 2010-2-5 ignore empty rows
				continue
			if row[0].find('#')==0:	#2010-4-28 ignore lines beginned with '#'
				continue
			phenotype_need_split = 0
			phenotype_id_index = col_name2index.get('phenotype_id')
			if phenotype_id_index == None:
				phenotype_id_index = col_name2index.get('phenot_id')
			if phenotype_id_index == None:
				phenotype_id_index = col_name2index.get('phenotype')
			if phenotype_id_index == None:	#2010-4-22
				phenotype_id_index = col_name2index.get('phenotype_method_id')
			if phenotype_id_index==None:	#2008-10-27 suzi's input '/Network/Data/250k/tmp-Suzi/sig_marker_ft_files/sig_marker_1_LD.csv' use 'phen' and the phenotype identifier, '1_LD', needs split to get id.
				phenotype_id_index = col_name2index.get('phen')
				phenotype_need_split = 1
			if phenotype_id_index is not None:
				if phenotype_need_split:	#2008-10-27 suzi's input '/Network/Data/250k/tmp-Suzi/sig_marker_ft_files/sig_marker_1_LD.csv' use 'phen' and the phenotype identifier, '1_LD', needs split to get id.
					phenotype_id = row[phenotype_id_index].split('_')[0]
				else:
					phenotype_id = row[phenotype_id_index]
				phenotype_id = int(phenotype_id)
			else:
				continue
			
			chromosome_index = col_name2index.get('chromosome')
			if chromosome_index ==None:
				chromosome_index = col_name2index.get('chr')
			
			position_index = col_name2index.get('position')
			if position_index==None:
				position_index = col_name2index.get('pos')
				if position_index is None:
					position_index = col_name2index.get('start')
			
			stop_index = col_name2index.get('stop')
			if stop_index is None:
				stop = None
			else:
				stop = int(row[stop_index])
			
			if chromosome_index is not None and position_index is not None:
				chromosome = row[chromosome_index]	#2012.3.8 kept in string format
				position = int(row[position_index])
				if stop:
					snps_id = '%s_%s_%s'%(chromosome, position, stop)
				else:
					snps_id = '%s_%s'%(chromosome, position)
			else:
				try:
					snps_id = int(row[col_name2index['snps_id']])
					snp_db_entry = snp_info.data_ls[snp_info.snps_id2index[snps_id]]
					chromosome = snp_db_entry.chromosome
					position = snp_db_entry.position
					
				except:	#forget it, skip
					continue
			
			gene_id_index = col_name2index.get('gene_id')
			if gene_id_index:
				gene_id = int(row[gene_id_index])
			else:
				gene_id = None
				
			fileNamePrefixIndex = col_name2index.get("fileNamePrefix")
			if fileNamePrefixIndex:
				fileNamePrefix = row[fileNamePrefixIndex]
			else:
				fileNamePrefix = None
			counter += 1
			#snps_id = int(row[col_name2index['snps_id']])
			#disp_pos = int(row[col_name2index['disp_pos']])
			this_snp = SNPPassingData(chromosome=chromosome, position=position, stop=stop, snps_id=snps_id, \
									fileNamePrefix=fileNamePrefix)
			input_data.addGeneSNPPhenotype(gene_id, this_snp, phenotype_id)
		sys.stderr.write("%s genes. %s total snps. Done.\n"%(input_data.no_of_genes, input_data.no_of_snps))
		return input_data
	
	def getSNPsFromRBGCandidateGenes(self, phenotype_method_id_ls, call_method_id, min_distance, \
									get_closest, min_MAF, candidate_gene_set, snp_info, results_directory=None, analysis_method_id_ls=None,\
									perc_of_top_candidate_genes=0.3):
		"""
		2008-10-23
			equivalent to getSNPsFromInputFile()
			
			iterate through a rbg file from top to bottom, if a gene belongs to candidate_gene_set, include its associated SNP
		"""
		sys.stderr.write("Get SNPs From RBG Candidate Genes ...")
		query = Stock_250kDB.ResultsByGene.query
		if call_method_id!=0:
			query = query.filter(Stock_250kDB.ResultsByGene.results_method.has(call_method_id=call_method_id))
		if analysis_method_id_ls:
			query = query.filter(Stock_250kDB.ResultsByGene.results_method.has(Stock_250kDB.ResultsMethod.analysis_method_id.in_(analysis_method_id_ls)))
		if phenotype_method_id_ls:
			query = query.filter(Stock_250kDB.ResultsByGene.results_method.has(Stock_250kDB.ResultsMethod.phenotype_method_id.in_(phenotype_method_id_ls)))
		
		query = query.filter_by(min_distance=min_distance).filter_by(get_closest=get_closest).\
			filter(Stock_250kDB.ResultsByGene.min_MAF>=min_MAF-0.0001).filter(Stock_250kDB.ResultsByGene.min_MAF<=min_MAF+0.0001)
		
		no_of_top_hits = len(candidate_gene_set)*perc_of_top_candidate_genes
		param_data = PassingData(results_directory=results_directory, no_of_top_lines=no_of_top_hits, \
								candidate_gene_set=candidate_gene_set)
		
		input_data = GeneSNPPhenotypeAssoData()
		counter = 0
		i = 0
		for rbg in query:
			i += 1
			pdata_ls = self.getTopResultsByGene(rbg, param_data)
			for pdata in pdata_ls:
				snp_db_entry = snp_info.data_ls[snp_info.snps_id2index[pdata.snps_id]]
				chromosome = snp_db_entry.chromosome
				position = snp_db_entry.position
				#chromosome, position = chr_pos
				this_snp = SNPPassingData(chromosome=chromosome, position=position, snps_id=pdata.snps_id, stop=None)
				input_data.addGeneSNPPhenotype(pdata.gene_id, this_snp, rbg.results_method.phenotype_method_id)
				counter += 1
		
		sys.stderr.write(" %s entries from %s rbgs. Done.\n"%(counter, i))
		return input_data
	
	@classmethod
	def findSNPsInRegion(cls, snp_info, chromosome=None, start=None, stop=None, center_snp_position=None):
		"""
		2009-4-30
			decide whether to use (chr,pos) or (chr,pos,offset) to represent a SNP based on snp_info.chr_pos2index
				?? use db-id to represent snp in snp_id2adjacent_window
			if snp_info has a RBDict storing all loci, then base-by-base searching is not necessary anymore (which is good for snp only)
		2008-10-1
			called by plotSNPRegion()
			find SNPs in this region, if center_snp_position is not given, find one.
			similar to getSNPsAroundThisSNP()
		"""
		sys.stderr.write("Get SNPs in this region (chr=%s, start=%s, stop=%s)..."%(chromosome, start, stop))
		chr_pos_ls = []
		chr_pos2adjacent_window = {}
		locus_id2adjacent_window = {}	#2012.3.7 because SNPID2PCAPosInfo.snp_id2img_x_pos have to be in db-id cuz drawSNPMtrix() uses db-id
		j = 0
		midpoint = (start+stop)/2.
		if center_snp_position is None:
			_center_snp_position = start
		else:
			_center_snp_position = center_snp_position
		center_snp = SNPPassingData(chromosome=chromosome, position=_center_snp_position, stopPosition=_center_snp_position, snps_id=None)
		
		segmentQueryKey = CNVSegmentBinarySearchTreeKey(chromosome=str(chromosome), \
							span_ls=[start, stop], min_reciprocal_overlap=0.0000001, )
		compareIns = CNVCompare(min_reciprocal_overlap=0.0000001)	#any overlap is an overlap
		
		node_ls = []
		snp_info.locusRBDict.findNodes(segmentQueryKey, node_ls=node_ls, compareIns=compareIns)
		snp_db_entry_ls = []
		for node in node_ls:
			segKey = node.key
			for snp_db_entry_index in node.value:
				snp_db_entry = snp_info.data_ls[snp_db_entry_index]
				snp_db_entry_ls.append(snp_db_entry)
		snp_db_entry_ls.sort(cmp=SNP.cmpDataObjByChrPos)
		adjacentWindowLs = [[obj.position, obj.end_position] for obj in snp_db_entry_ls]
		
		#2012.3.8 fill up the adjacentWindowLs
		for i in xrange(1, len(snp_db_entry_ls)):
			snp_db_entry = snp_db_entry_ls[i]
			prev_snp_db_entry = snp_db_entry_ls[i-1]
			mid_point = (prev_snp_db_entry.end_position + snp_db_entry.position)/2.
			if prev_snp_db_entry.end_position==prev_snp_db_entry.position:	#it's a SNP, rather than a segment
				adjacentWindowLs[i-1][1] = mid_point
			if snp_db_entry.end_position == snp_db_entry.position:
				adjacentWindowLs[i][0] = mid_point
			if center_snp_position is None and abs(snp_db_entry.position-midpoint)<abs(center_snp.position-midpoint):	#this SNP is closer to the center
					center_snp.position = snp_db_entry.position
					center_snp.stopPosition = snp_db_entry.end_position
		
		fst_chr_pos = snp_info.chr_pos2index.keys()[0]
		snp_representation_type = len(fst_chr_pos)
		"""
		#2012.3.8 old way of finding out adjacent windows
		#2009-3-27 get a SNP representation from snp_info to see whether the SNP is represented by chr_pos or chr_pos_offset
		
		for i in range(start-1, stop+2):
			new_pos = i
			if snp_representation_type==3:	#2009-3-27, 3-number representation
				new_chr_pos = (chromosome, new_pos, 0)
			else:
				new_chr_pos = (chromosome, new_pos, new_pos)
			if new_chr_pos in snp_info.chr_pos2index:
				if center_snp_position is None and abs(new_pos-midpoint)<abs(center_snp.position-midpoint):	#this SNP is closer to the center
					center_snp.position = new_pos
					center_snp.stopPosition = new_pos
				chr_pos_ls.append(new_chr_pos)
				if j!=0:
					cls.add_mid_point(chr_pos_ls, chr_pos2adjacent_window)
				j += 1
		"""
		if len(snp_db_entry_ls)>1:
			#deal with the leftest point of the 1st chr_pos
			window_size = adjacentWindowLs[0][1] - adjacentWindowLs[0][0]
			adjacentWindowLs[0][0] -= window_size
		
			#deal with the rightest point of the 1st chr_pos
			window_size = adjacentWindowLs[-1][1] - adjacentWindowLs[-1][0]
			adjacentWindowLs[0][1] += window_size
		
		for i in xrange(len(snp_db_entry_ls)):
			snp_db_entry = snp_db_entry_ls[i]
			chr_pos = (snp_db_entry.chromosome, snp_db_entry.position, snp_db_entry.end_position)
			chr_pos_ls.append(chr_pos)
			chr_pos2adjacent_window[chr_pos] = adjacentWindowLs[i]
			locus_id2adjacent_window[snp_db_entry.id] = adjacentWindowLs[i]
		center_snp.snps_id = '%s_%s_%s'%(center_snp.chromosome, center_snp.position, center_snp.stopPosition)
		snp_region = PassingData(chr_pos_ls=chr_pos_ls, chr_pos2adjacent_window=chr_pos2adjacent_window, center_snp=center_snp,\
								snp_representation_type=snp_representation_type, locus_id2adjacent_window=locus_id2adjacent_window,\
								chromosome=chromosome, start=start, stop=stop)
		sys.stderr.write("%s loci.\n"%(len(chr_pos2adjacent_window)))
		return snp_region
	
	def loadDataStructure(self, gene_annotation_picklef, LD_info_picklef, LD_fname=None, min_MAF=0.1, min_distance=20000, \
						list_type_id_list=None, snp_matrix_fname=None, phenotype_fname=None, snp_matrix_before_impute_fname=None,\
						snp_matrix_data_type=1, locus_type_id=1, call_method_id=None):
		"""
		2012.3.8
			add argument call_method_id to fetch a sample ResultsMethod to determine locus_type_id
		2009-3-23
			add argument snp_matrix_data_type to allow CNV or CNV amplitude data
		2008-12-04
			add snp_matrix_before_impute_fname
		2008-11-25
			if snp_matrix_fname is provided
				1. no LD_info will be loaded.
				2. load snpData from snp_matrix_fname
				3. load phenData from phenotype_fname
				4. get ecotype_info
				
		2008-10-01
			wrap a few functions up, convenient for both run() and drawSNPRegion()
		"""
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		self.db = db
		#2012.3.8 get locus_type_id from a sample ResultsMethod
		rm = Stock_250kDB.ResultsMethod.query.filter_by(call_method_id=call_method_id).first()
		locus_type_id = rm.locus_type_id
		gene_annotation = self.dealWithGeneAnnotation(gene_annotation_picklef)
		snp_info = db.dealWithSNPInfo(self.snpInfoPickleFname, locus_type_id=locus_type_id)	#2012.3.8
		if snp_matrix_fname:
			LD_info = None
		else:
			LD_info = self.dealLD_info(LD_info_picklef, LD_fname, min_MAF, min_distance)
		candidate_gene_set = set()
		if list_type_id_list:
			for list_type_id in list_type_id_list:
				candidate_gene_list = db.getGeneList(list_type_id)
				candidate_gene_set |= set(candidate_gene_list)
		
		if snp_matrix_fname:
			if snp_matrix_data_type==3:
				matrix_data_type=float		#2009-3-23 for CNV amplitude file
			else:
				matrix_data_type=int
			snpData = SNPData(input_fname=snp_matrix_fname, turn_into_integer=1, turn_into_array=1, ignore_2nd_column=1,\
							matrix_data_type=matrix_data_type)
			if snpData.data_matrix is None:
				sys.stderr.write("Error. snpData.data_matrix is None.\n")
				sys.exit(3)
			header_phen, strain_acc_list_phen, category_list_phen, data_matrix_phen = read_data(phenotype_fname, turn_into_integer=0)
			
			phenData = SNPData(header=header_phen, strain_acc_list=snpData.strain_acc_list, data_matrix=data_matrix_phen)
			#row label is that of the SNP matrix, because the phenotype matrix is gonna be re-ordered in that way
			
			phenData.data_matrix = Kruskal_Wallis.get_phenotype_matrix_in_data_matrix_order(snpData.row_id_ls, \
																		strain_acc_list_phen, phenData.data_matrix)
			#tricky, using strain_acc_list_phen
			
			#2008-12-05 fake a snp_info for findSNPsInRegion
			self.construct_chr_pos2index_forSNPData(snpData, snp_info=snp_info)
			"""
			chr_pos2index = {}
			for i in range(len(snpData.col_id_ls)):
				col_id = snpData.col_id_ls[i]
				chr_pos = col_id.split('_')
				if len(chr_pos)==3 and chr_pos[2]!='0':	#2009-3-27 ignore insertion right now
					continue
				chr_pos = tuple(map(int, chr_pos))
				chr_pos2index[chr_pos] = i
			snpData.chr_pos2index = chr_pos2index
			"""
			ecotype_info = getEcotypeInfo(db)
		else:
			snpData = None
			phenData = None
			ecotype_info = None
		
		if self.snp_matrix_before_impute_fname:
			snpData_before_impute = SNPData(input_fname=self.snp_matrix_before_impute_fname, \
							turn_into_integer=1, turn_into_array=1, ignore_2nd_column=1)
			#snpData_before_impute, allele2index_ls = snpData_before_impute.convertSNPAllele2Index(self.report)
		else:
			snpData_before_impute = None
		return_data = PassingData(gene_annotation=gene_annotation, snp_info=snp_info, LD_info=LD_info, \
								candidate_gene_set=candidate_gene_set, snpData=snpData, phenData=phenData,\
								ecotype_info=ecotype_info, snpData_before_impute=snpData_before_impute)
		return return_data
	
	@classmethod
	def getNonNAPhenotypeRowIndexLs(cls, phenotype_ls):
		"""
		2010-2-2
			get the indices that are not NA
		"""
		non_NA_phenotype_row_index_ls = []
		for i in range(len(phenotype_ls)):
			if not numpy.isnan(phenotype_ls[i]):
				non_NA_phenotype_row_index_ls.append(i)
		return non_NA_phenotype_row_index_ls
	
	@classmethod
	def handleOneGeneOnePhenotype(cls, db, phenotype_id2analysis_method_id2gwr, phenotype_id, this_snp, input_data, \
								param_data, latex_f=None, delete_gwr=False, plot_type=None, commit=0):
		"""
		2010-4-26
			pass param_data.drawMap, param_data.markSNPInGenes, param_data.drawGrid, param_data.phenotypeDrawType
				to drawRegionAroundThisSNP()
		2010-2-2
			deal with param_data.exclude_accessions_with_NA_phenotype, an argument controlling whether the SNP matrix
				would contain accessions with NA phenotype or not
			pass getattr(param_data, 'analysis_method_id_ls', [1,7]) to getSimilarGWResultsGivenResultsByGene()
		2008-10-26
			fix a bug in the stop position when submitting data into db.
			previously "stop = after_plot_data.snps_within_this_region.chr_pos_ls[1][1]", the 1st "1" should be "-1".
		2008-10-24
			pass commit to drawRegionAroundThisSNP
		2008-10-02
			wrap up towards latex output
		"""
		analysis_method_id2gwr = phenotype_id2analysis_method_id2gwr.get(phenotype_id)
		if not analysis_method_id2gwr:
			analysis_method_id2gwr = cls.getSimilarGWResultsGivenResultsByGene(db_250k=db, phenotype_method_id=phenotype_id, \
											call_method_id=param_data.call_method_id, \
											results_directory = param_data.results_directory,\
											analysis_method_id_ls=getattr(param_data, 'analysis_method_id_ls', [1,7]))
			if not delete_gwr:
				phenotype_id2analysis_method_id2gwr[phenotype_id] = analysis_method_id2gwr
		if not analysis_method_id2gwr:	#still nothing, skip
			sys.stderr.write("No gwas results available. skip.\n")
			return
		
		# 2010-2-2
		exclude_accessions_with_NA_phenotype = getattr(param_data, 'exclude_accessions_with_NA_phenotype', False)
		if exclude_accessions_with_NA_phenotype:
			phenotype_index, = cls.findOutWhichPhenotypeColumn(param_data.phenData, set([phenotype_id]))
			non_NA_phenotype_row_index_ls = cls.getNonNAPhenotypeRowIndexLs(param_data.phenData.data_matrix[:,phenotype_index])
			phenData = SNPData.keepRowsByRowIndex(param_data.phenData, non_NA_phenotype_row_index_ls)
			snpData = SNPData.keepRowsByRowIndex(param_data.snpData, non_NA_phenotype_row_index_ls)
			cls.construct_chr_pos2index_forSNPData(snpData, snp_info=param_data.snp_info)
			if param_data.snpData_before_impute:
				snpData_before_impute = SNPData.keepRowsByRowIndex(param_data.snpData_before_impute, non_NA_phenotype_row_index_ls)
			else:
				snpData_before_impute = None
		else:
			phenData = param_data.phenData
			snpData = param_data.snpData
			snpData_before_impute = param_data.snpData_before_impute
		after_plot_data = cls.drawRegionAroundThisSNP(phenotype_id, this_snp, param_data.candidate_gene_set, param_data.gene_annotation, \
								param_data.snp_info, \
								analysis_method_id2gwr, param_data.LD_info, param_data.output_dir, param_data.which_LD_statistic, \
								min_distance=param_data.min_distance, list_type_id_list=param_data.list_type_id_list,
								label_gene=param_data.label_gene, \
								draw_LD_relative_to_center_SNP=param_data.draw_LD_relative_to_center_SNP,\
								commit=commit, snpData=snpData, phenData=phenData, \
								ecotype_info=param_data.ecotype_info, snpData_before_impute=snpData_before_impute,\
								snp_matrix_data_type=param_data.snp_matrix_data_type, call_method_id=param_data.call_method_id,\
								drawMap=param_data.drawMap, markSNPInGenes=param_data.markSNPInGenes, \
								drawGrid=param_data.drawGrid, phenotypeDrawType=param_data.phenotypeDrawType)
		param_data.no_of_snps_drawn += 1
		if commit and after_plot_data:
			#2008-10-24 first check if it's in db or not
			start = after_plot_data.snps_within_this_region.chr_pos_ls[0][1]
			stop = after_plot_data.snps_within_this_region.chr_pos_ls[-1][1]
			rows = Stock_250kDB.SNPRegionPlot.query.filter_by(chromosome=this_snp.chromosome).\
													filter_by(start=start).\
													filter_by(stop=stop).\
													filter_by(phenotype_method_id=phenotype_id).\
													filter_by(plot_type_id=plot_type.id)
			if rows.count()>0:
				row = rows.first()
				sys.stderr.write("region plot (chr=%s, start=%s, stop=%s, phenotype_method_id=%s, plot_type_id=%s) already in db(id=%s). skip.\n"%\
								(row.chromosome, row.start, row.stop, row.phenotype_method_id, row.plot_type_id, row.id))
				return
			snp_region_plot = Stock_250kDB.SNPRegionPlot(chromosome=this_snp.chromosome, \
														start=start,\
														stop=stop,\
													center_snp_position=after_plot_data.snps_within_this_region.center_snp.position,\
													phenotype_method_id=phenotype_id)
			snp_region_plot.plot_type = plot_type
			snp_region_plot.png_data = after_plot_data.png_data.getvalue()
			snp_region_plot.svg_data = after_plot_data.svg_data.getvalue()
			db.session.save(snp_region_plot)
			
			#save related genes into db as well			
			gene_id2count = {}
			for gene_desc_ls in after_plot_data.matrix_of_gene_descriptions:
				gene_id = gene_desc_ls[0]
				if gene_id not in gene_id2count:
					gene_id2count[gene_id] = 0
					plot2gene = Stock_250kDB.SNPRegionPlotToGene(gene_id=gene_id)
					plot2gene.snp_region_plot = snp_region_plot
					db.session.save(plot2gene)
				gene_id2count[gene_id] += 1
			db.session.flush()
		
		if delete_gwr:
			del analysis_method_id2gwr
		
		if latex_f:
			annot_table_label = 'annottable%s'%param_data.no_of_snps_drawn
			snp_figure_label = 'snpfig%s'%param_data.no_of_snps_drawn
			no_of_genes_in_annot = len(after_plot_data.matrix_of_gene_descriptions)
			snp_label_region = 'SNP at %s, %s'%(this_snp.chromosome, this_snp.position)
			if getattr(this_snp, 'stop', None):
				snp_label_region += ' to %s'%this_snp.stop
			
			caption = 'Annotation of %s genes for Figure~\\ref{%s}. %s. Phenotype is %s(id=%s).'%\
				(no_of_genes_in_annot, snp_figure_label, snp_label_region, \
				input_data.return_phenotype_short_name(phenotype_id), phenotype_id)
			
			if after_plot_data.matrix_of_gene_descriptions:	#if matrix is empty, don't output
				latex_f.write(outputMatrixInLatexTable(after_plot_data.matrix_of_gene_descriptions, caption, annot_table_label, cls.gene_desc_names))
				annot_table_ref = 'Annotation of %s genes is in Table~\\ref{%s}.'%\
					(no_of_genes_in_annot, annot_table_label)
			else:	#no mentioning this table if it doesn't exist
				annot_table_ref = ''
			caption = '%s. Phenotype is %s(id=%s). %s'%(snp_label_region, \
						input_data.return_phenotype_short_name(phenotype_id), phenotype_id, annot_table_ref)
			latex_f.write(outputFigureInLatex(after_plot_data.png_output_fname, caption, snp_figure_label))
			latex_f.flush()
	
	def run(self):
		"""
		2008-10-02
			modify to enable latex output of tables and figures
		2008-09-24
		"""
		if self.debug==1:
			import pdb
			pdb.set_trace()
		import gc
		gc.set_threshold(250,10,10)
		
		grand_dataStructure = self.loadDataStructure(self.gene_annotation_picklef, self.LD_info_picklef, self.LD_fname, \
							self.min_MAF, self.min_distance, self.list_type_id_list, snp_matrix_fname=self.snp_matrix_fname,\
							phenotype_fname=self.phenotype_fname, snp_matrix_before_impute_fname=self.snp_matrix_before_impute_fname, \
							snp_matrix_data_type=self.snp_matrix_data_type, call_method_id=self.call_method_id)
		gene_annotation, snp_info, LD_info = grand_dataStructure.gene_annotation, grand_dataStructure.snp_info, grand_dataStructure.LD_info
		if not os.path.isdir(self.output_dir):
			os.makedirs(self.output_dir)
		"""
		if not self.results_id_ls:
			param_obj = PassingData(call_method_id=self.call_method_id, analysis_method_id=self.analysis_method_id,\
								min_distance=self.min_distance, get_closest=self.get_closest)
			self.results_id_ls = self.generate_params(param_obj)
		"""
		which_LD_statistic = self.which_LD_statistic
		input_fname = self.input_fname
		latex_fname = self.latex_output_fname
		
		phenotype_id2analysis_method_id2gwr = {}
		
		param_data = grand_dataStructure
		param_data.call_method_id = self.call_method_id
		param_data.results_directory = self.data_dir
		param_data.output_dir = self.output_dir
		param_data.which_LD_statistic = which_LD_statistic
		param_data.min_distance = self.min_distance
		param_data.list_type_id_list = self.list_type_id_list
		param_data.no_of_snps_drawn = 0
		param_data.label_gene = self.label_gene
		param_data.draw_LD_relative_to_center_SNP = self.draw_LD_relative_to_center_SNP
		param_data.snp_matrix_data_type = self.snp_matrix_data_type
		param_data.analysis_method_id_ls = self.analysis_method_id_ls
		param_data.exclude_accessions_with_NA_phenotype = self.exclude_accessions_with_NA_phenotype
		param_data.drawMap = self.drawMap	#2010-4-13
		param_data.markSNPInGenes = self.markSNPInGenes	#2010-4-23
		param_data.drawGrid = self.drawGrid	#2010-4-23
		param_data.phenotypeDrawType = self.phenotypeDrawType	# 2010-4-26
		
		no_of_genes_handled = 0
		
		#2008-10-23 handle plot_type
		if self.commit:
			plot_type = Stock_250kDB.SNPRegionPlotType.query.filter_by(short_name=self.plot_type_short_name).first()
			if not plot_type:
				plot_type = Stock_250kDB.SNPRegionPlotType(short_name=self.plot_type_short_name)
				self.db.session.save(plot_type)
				self.db.session.flush()
		else:
			plot_type = None
		
		while 1:
			#try:
			if input_fname:
				input_data = self.getSNPsFromInputFile(input_fname, snp_info=snp_info, db_250k=self.db)
			else:
				input_data = self.getSNPsFromRBGCandidateGenes(self.phenotype_method_id_ls, self.call_method_id, self.min_distance, \
								self.get_closest, self.min_MAF, grand_dataStructure.candidate_gene_set, snp_info, \
								self.rbg_results_directory, \
								self.analysis_method_id_ls, perc_of_top_candidate_genes=0.3)
			
			param_data.which_LD_statistic = which_LD_statistic
			"""
			for gene_id in input_data.gene_id_ls:
				sys.stderr.write("Handling Gene %s ...\n"%gene_id)
				snp2phenotype_id_ls = input_data.return_snp_data_given_gene_id(gene_id)
				gene_model = gene_annotation.gene_id2model.get(gene_id)
				if not gene_model:
					continue
				no_of_genes_handled += 1
				snp_ls = snp2phenotype_id_ls.keys()
				snp_ls.sort()
				if latex_fname:
					latex_fname_appx = 'gene_%s_%s'%(gene_id, getattr(gene_model, 'gene_symbol', ''))
					latex_fname_appx = latex_fname_appx.replace('/', '_')
					_latex_fname = '%s_%s'%(latex_fname, latex_fname_appx)
					latex_f = open(_latex_fname, 'w')
					matrix_of_gene_descriptions = []
					if len(gene_model.gene_commentaries)>0:
						for gene_commentary in gene_model.gene_commentaries:
							gene_desc_ls = self.returnGeneDescLs(self.gene_desc_names, gene_model, gene_commentary)
							matrix_of_gene_descriptions.append(gene_desc_ls)
					else:
						gene_desc_ls = self.returnGeneDescLs(self.gene_desc_names, gene_model)
						matrix_of_gene_descriptions.append(gene_desc_ls)
					gene_table_label = 'gene%s'%no_of_genes_handled
					snps_table_label = 'snps%s'%no_of_genes_handled
					caption = 'Gene %s (id=%s) has %s SNPs (check Tabel~\\ref{%s}) and %s phenotypes associated with.'%\
						(gene_model.gene_symbol, gene_id, len(snp2phenotype_id_ls), snps_table_label, input_data.get_no_of_phenotypes_given_gene_id(gene_id))
					latex_f.write(outputMatrixInLatexTable(matrix_of_gene_descriptions, caption, gene_table_label, self.gene_desc_names))
					
					matrix_of_snp_descriptions = input_data.return_matrix_of_snp_descriptions(gene_id)
					caption = '%s SNPs are related to Gene %s, id=%s, (check Tabel~\\ref{%s}).'%\
						(len(snp_ls), gene_model.gene_symbol, gene_id, gene_table_label)
					latex_f.write(outputMatrixInLatexTable(matrix_of_snp_descriptions, caption, snps_table_label, input_data.snp_desc_names))
				else:
					latex_f = None
				for snp in snp_ls:
					phenotype_id_ls = snp2phenotype_id_ls[snp]
					this_snp = SNPPassingData(chromosome=snp[0], position=snp[1], snps_id=snp[2])
					for phenotype_id in phenotype_id_ls:
						self.handleOneGeneOnePhenotype(self.db, phenotype_id2analysis_method_id2gwr, phenotype_id, this_snp, \
													input_data, param_data, latex_f)
				if latex_f:
					del latex_f
			"""
			snps_id2latex_f = {}
			if input_data.no_of_genes>=0:	#no gene in input_fname, iterate over snps
				phenotype_id2snp_ls = input_data.return_phenotype_id2snp_ls()
				for phenotype_id, snp_ls in phenotype_id2snp_ls.iteritems():
					snp_ls.sort()
					for snp in snp_ls:
						
						this_snp = SNPPassingData(chromosome=snp[0], position=snp[1], snps_id=snp[2], stop=snp[3],
												fileNamePrefix=snp[4])
						latex_f = snps_id2latex_f.get(this_snp.snps_id)
						if this_snp.snps_id not in snps_id2latex_f and latex_fname:
							latex_f = open('%s_%s'%(latex_fname, this_snp.snps_id), 'w')
							snps_id2latex_f[this_snp.snps_id] = latex_f
						self.handleOneGeneOnePhenotype(self.db, phenotype_id2analysis_method_id2gwr, phenotype_id, this_snp, \
													input_data, param_data, latex_f, plot_type=plot_type, commit=self.commit)
				del phenotype_id2analysis_method_id2gwr	#clean up the memory
				phenotype_id2analysis_method_id2gwr = {}
				"""
				
				for snp, phenotype_id_ls in input_data.snp2phenotype_id_ls.iteritems():
					this_snp = SNPPassingData(chromosome=snp[0], position=snp[1], snps_id=snp[2])
					if latex_fname:
						latex_f = open('%s_%s'%(latex_fname, this_snp.snps_id), 'w')
					else:
						latex_f = None
					for phenotype_id in phenotype_id_ls:
						self.handleOneGeneOnePhenotype(self.db, phenotype_id2analysis_method_id2gwr, phenotype_id, this_snp, \
													input_data, param_data, latex_f)
				"""
			del snps_id2latex_f
			#except:
			#	sys.stderr.write('Except: %s\n'%repr(sys.exc_info()))
			#	traceback.print_exc()
				#raise
			if self.logFilename:
				logF = open(self.logFilename, 'w')
				logF.write("%s phenotype(s).\n"%(len(input_data.return_phenotype_id2snp_ls())))
				del logF
				
			sys.exit(0)	#2008-10-24 don't wanna the while loop
			
			to_continue = raw_input("Continue running?([y]/n): ")
			if to_continue.upper()=='N':
				break
			_input_fname = raw_input("File containing phenotype_id, chromosome, position, default is %s: "%input_fname)
			_which_LD_statistic = raw_input("which LD statistic (1=r2, 2=|D_prime|), default is %s: "%(which_LD_statistic))
			_latex_fname = raw_input("Latex Filename (or prefix) to store gene description tables and figures, default is %s: "%(latex_fname))
			if _input_fname:
				input_fname = _input_fname
			if _which_LD_statistic:
				which_LD_statistic = _which_LD_statistic
			if _latex_fname:
				latex_fname = _latex_fname

if __name__ == '__main__':
	from pymodule import yh_matplotlib
	yh_matplotlib.setFontAndLabelSize(6)
	
	from pymodule import ProcessOptions
	main_class = DrawSNPRegion
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
