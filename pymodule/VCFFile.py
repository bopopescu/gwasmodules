#!/usr/bin/env python
"""
2011-9-27
	a module to handle vcf file http://www.1000genomes.org/node/101
	
	need to add row-based iterator
"""
import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from ProcessOptions import  ProcessOptions
from utils import dict_map, importNumericArray, figureOutDelimiter, PassingData, getColName2IndexFromHeader
import copy, csv



class VCFFile(object):
	__doc__ = __doc__
	option_default_dict = {('inputFname', 1, ): [None, 'o', 1, 'a VCF input file.'],\
						('openMode', 1, ): ['rb', '', 1, 'rb: bam file. r: sam file.'],\
						('minorAlleleDepthLowerBoundCoeff', 1, float): [1/4., 'M', 1, 'minimum read depth multiplier for an allele to be called (heterozygous or homozygous)', ],\
						('minorAlleleDepthUpperBoundCoeff', 1, float): [3/4., 'A', 1, 'maximum read depth multiplier for the minor allele of a heterozygous call', ],\
						('majorAlleleDepthUpperBoundCoeff', 1, float): [7/8., 'a', 1, 'maximum read depth multiplier for the major allele of het call'],\
						('maxNoOfReadsForGenotypingError', 1, float): [1, 'x', 1, 'if read depth for one allele is below or equal to this number, regarded as genotyping error ', ],\
						('depthUpperBoundCoeff', 1, float): [2, 'm', 1, 'depthUpperBoundCoeff*coverage = maximum read depth for one base'],\
						('depthLowerBoundCoeff', 1, float): [1/4., 'O', 1, 'depthLowerBoundCoeff*coverage = minimum read depth multiplier for one base'],\
						('depthUpperBoundMultiSampleCoeff', 1, float): [3, 'N', 1, 'across n samples, ignore bases where read depth > coverageSum*depthUpperBoundCoeff*multiplier.'],\
						('seqCoverageFname', 0, ): ['', 'q', 1, 'The sequence coverage file. tab/comma-delimited: individual_sequence.id coverage'],\
						('defaultCoverage', 1, float): [5, 'd', 1, 'default coverage when coverage is not available for a read group'],\
						("site_type", 1, int): [2, 's', 1, '1: all sites, 2: variants only'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self, **keywords):
		"""
		2011-9-27
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		
		self.header = None
		self.sample_id_ls = []
		self.sample_id2index = {}
		self.locus_id_ls = []
		self.locus_id2row_index = {}
		self.locus_id2data = {}
		self.genotype_call_matrix = []
		self.parseFile(self.inputFname)
	
	@classmethod
	def getIndividual2ColIndex(cls, header, col_name2index, sampleStartingColumn=9):
		"""
		2011-9-27
			called by parseFile
		"""
		sys.stderr.write("\t Finding all individuals ...")
		no_of_cols = len(header)
		individual_name2col_index = {}	#individual's column name -> an opened file handler to store genetic data
		counter = 0
		for i in xrange(sampleStartingColumn, no_of_cols):
			individualName = header[i]
			col_index = col_name2index.get(individualName)
			if not individualName:	#ignore empty column
				continue
			if individualName[:-4]=='.bam':
				individualCode = individualName[:-4]	#get rid of .bam
			else:
				individualCode = individualName
			individual_name2col_index[individualCode] = col_index
			counter += 1
		sys.stderr.write("%s individuals added. Done.\n"%(counter))
		return individual_name2col_index
	
	
	def parseFile(self, inputFname=None):
		"""
		2011-9-27
			based on discoverFromVCF() from vervet/src/GenotypeCallByCoverage.py
		"""
		sys.stderr.write("Parsing %s ... \n"%(os.path.basename(inputFname),))
		if inputFname[-3:]=='.gz':
			import gzip
			inf = gzip.open(inputFname)
		else:
			inf = open(inputFname)
		reader =csv.reader(inf, delimiter='\t')
		
		sample_id2index = self.sample_id2index
		locus_id2row_index = self.locus_id2row_index
		sample_id2index['ref'] = 0	#ref is at column 0. "ref" must not be equal to any read_group.
		self.sample_id_ls.append('ref')
		read_group2coverage = {}	#2011-9-2
		
		
		tid2refName = {}	#dictionary storing the target references which have SNP calls
		refNameSet = set()
		"""
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		header = ['sample', 'snp_id', 'chr', 'pos', 'qual', 'DP', 'minDP4', 'DP4_ratio', 'MQ']
		moreHeader = ['GQ', 'GL', 'SB', 'QD', 'sndHighestGL', 'deltaGL']
		#['AF', 'AC','AN', 'Dels', 'HRun', 'HaplotypeScore','MQ0', 'QD']	#2011-3-4 useless
		if VCFOutputType==2:
			header += moreHeader
		chr_pure_number_pattern = re.compile(r'[a-z_A-Z]+(\d+)')
		chr_number_pattern = re.compile(r'chr(\d+)')
		"""
		
		individual_name2col_index = None
		col_name2index = None
		counter = 0
		real_counter = 0
		
		for row in reader:
			if row[0] =='#CHROM':
				row[0] = 'CHROM'	#discard the #
				header = row
				col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
				individual_name2col_index = self.getIndividual2ColIndex(header, col_name2index)
				continue
			elif row[0][0]=='#':	#2011-3-4
				continue
			"""
			if chr_number_pattern.search(row[0]):
				chr = chr_number_pattern.search(row[0]).group(1)
			elif chr_pure_number_pattern.search(row[0]):
				chr = chr_pure_number_pattern.search(row[0]).group(1)
			else:
				sys.stderr.write("Couldn't parse the chromosome number/character from %s.\n Exit.\n"%(row[0]))
				sys.exit(4)
			"""
			chr = row[0]
			refNameSet.add(chr)
			
			pos = row[1]
			quality = row[5]
			
			info = row[7]
			info_ls = info.split(';')
			info_tag2value = {}
			for info in info_ls:
				try:
					tag, value = info.split('=')
				except:
					#sys.stderr.write("Error in splitting %s by =.\n"%info)	###Error in splitting DS by =.
					continue
				info_tag2value[tag] = value
			
			current_locus = '%s_%s'%(chr, pos)
			refBase = row[col_name2index['REF']]
			altBase = row[col_name2index['ALT']]
			if current_locus not in self.locus_id2data:
				self.locus_id2data[current_locus] = info_tag2value
				self.locus_id2data[current_locus]['REF'] = refBase
				self.locus_id2data[current_locus]['ALT'] = altBase
			
			format_column = row[col_name2index['FORMAT']]
			format_column_ls = format_column.split(':')
			format_column_name2index = getColName2IndexFromHeader(format_column_ls)
			data_row = ['NA']*(len(individual_name2col_index)+1)	# extra 1 for the ref
			data_row[0] = refBase
			allele2count = {}
			for individual_name, individual_col_index in individual_name2col_index.iteritems():
				read_group = individual_name
				if read_group not in sample_id2index:
					sample_id2index[read_group] = len(sample_id2index)
					self.sample_id_ls.append(read_group)
				
				#coverage = read_group2coverage[read_group]
				genotype_data = row[individual_col_index]
				genotype_data_ls = genotype_data.split(':')
				genotype_call_index = format_column_name2index.get('GT')
				genotype_quality_index = format_column_name2index.get('GQ')
				if genotype_quality_index is None:
					genotype_quality_index = format_column_name2index.get('DP')
				depth_index = format_column_name2index.get("DP")
				#GL_index = format_column_name2index.get('GL')
				if len(genotype_data_ls)<len(format_column_name2index):	#this genotype call is probably empty "./." due to no reads
					continue
				if depth_index is None or genotype_call_index is None:
					sys.stderr.write("")
					continue
				#genotype_quality = genotype_data_ls[genotype_quality_index]
				genotype_call = genotype_data_ls[genotype_call_index]
				depth = int(genotype_data_ls[depth_index])
				if depth<1:	#no read. samtools would still assign ref/ref to this individual
					continue
				#if depth>maxNoOfReads*coverage or depth<minNoOfReads*coverage:	#2011-3-29 skip. coverage too high or too low
				#	continue
				allele = 'NA'
				if genotype_call=='0/1' or genotype_call =='1/0':	#heterozygous, the latter notation is never used though.
					allele = '%s%s'%(refBase, altBase)
					"""
					GL_list = genotype_data_ls[GL_index]
					GL_list = GL_list.split(',')
					GL_list = map(float, GL_list)
					GL = GL_list[1]
					sndHighestGL = max([GL_list[0], GL_list[2]])
					deltaGL = GL-sndHighestGL
					"""
					"""
					
					AD = genotype_data_ls[format_column_name2index.get('AD')]
					AD = map(int, AD.split(','))
					minorAlleleCoverage = min(AD)
					majorAlleleCoverage = max(AD)
					
					if minorAlleleCoverage<=self.minorAlleleDepthUpperBoundCoeff*coverage and \
							minorAlleleCoverage>=self.minorAlleleDepthLowerBoundCoeff*coverage and \
							majorAlleleCoverage<=self.majorAlleleDepthUpperBoundCoeff*coverage:
						DP4_ratio = float(AD[0])/AD[1]
						allele = '%s%s'%(refBase, altBase)
					"""
				elif genotype_call=='./.':	#missing
					pass
				elif genotype_call =='1/1':
					allele = '%s%s'%(altBase, altBase)
				elif genotype_call =='0/0':
					allele = '%s%s'%(refBase, refBase)
				col_index = sample_id2index.get(read_group)
				data_row[col_index] = allele
				if allele!='NA':
					if allele not in allele2count:
						allele2count[allele] = 0
					allele2count[allele] += 1
				
			if len(allele2count)>self.site_type-1:	#whether polymorphic across samples or all sites in vcf
				real_counter += 1
				locus_id2row_index[current_locus] = len(locus_id2row_index)
				self.locus_id_ls.append(current_locus)
				self.genotype_call_matrix.append(data_row)
			
			counter += 1
			if counter%2000==0 and self.report:
				sys.stderr.write("%s\t%s\t%s"%("\x08"*80, counter, real_counter))
		del reader
		sys.stderr.write("%s loci X %s samples.\n"%(len(self.locus_id_ls), len(self.sample_id_ls)))
		