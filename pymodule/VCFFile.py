#!/usr/bin/env python
"""
Example:

	vcfFile = VCFFile(inputFname=self.inputFname)
	trio_col_index_data = self.findTrioIndex(vcfFile.sample_id2index, self.trio_isq_id_ls)
	for vcfRecord in vcfFile.parseIter():
		locus_id = vcfRecord.locus_id
		chr = vcfRecord.chr
		pos = vcfRecord.pos
		pos = int(pos)
	
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
import re
#".|.", "0|1" for phased 
#"./.", "0/1" for unphased
diploidGenotypePattern = re.compile(r'([\d.])[|/]([\d.])')

def parseOneVCFRow(row, col_name2index, col_index_individual_name_ls, sample_id2index, minDepth=1,\
				dataEntryType=1):
	"""
	2012.9.6 turn pos into integer
	2012.5.10
		complete representation of one locus
	2012.1.17
		common snippet split out of VCFFile & VCFRecord
		row is a list of input columns from one VCF file line
		dataEntryType
			1: each cell is base call
			2: each cell is a dictionary {'GT': base-call, 'DP': depth}
	"""
	chr = row[0]
	pos = int(row[1])	#2012.9.6 turn pos into integer
	vcf_locus_id=row[2]
	quality = row[5]
	filter=row[6]
	info = row[7]
	format = row[8]
	info_ls = info.split(';')
	info_tag2value = {}
	for info_entry in info_ls:
		try:
			tag, value = info_entry.split('=')
		except:
			#sys.stderr.write("Error in splitting %s by =.\n"%info)	###Error in splitting DS by =.
			continue
		info_tag2value[tag] = value
	
	locus_id = (chr, pos)
	refBase = row[col_name2index['REF']]
	altBase = row[col_name2index['ALT']]
	
	altBaseLs = altBase.split(',')	#altBase could be just "C" or "C,G" (multi-nucleotide)
	alleleLs = [refBase] + altBaseLs
	alleleNumber2Base = {'.':'NA'}
	for i in xrange(len(alleleLs)):
		alleleNumber2Base[repr(i)] = alleleLs[i]
	
	format_column = row[col_name2index['FORMAT']]
	format_column_ls = format_column.split(':')
	format_column_name2index = getColName2IndexFromHeader(format_column_ls)
		
	if dataEntryType==1:
		data_row = ['NA']*(len(col_index_individual_name_ls)+1)	# extra 1 for the ref
		data_row[0] = refBase
	else:
		data_row = [None]*(len(col_index_individual_name_ls)+1)	# extra 1 for the ref
		data_row[0] = {'GT':refBase, 'DP':-1}
	genotypeCall2Count = {}
	for individual_col_index, individual_name in col_index_individual_name_ls:
		individual_name = individual_name
		if individual_name not in sample_id2index:
			sample_id2index[individual_name] = len(sample_id2index)
		
		#coverage = read_group2coverage[individual_name]
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
		if genotype_call_index is None:
			sys.stderr.write("")
			continue
		#genotype_quality = genotype_data_ls[genotype_quality_index]
		genotype_call = genotype_data_ls[genotype_call_index]
		callData = {}
		if genotype_call=='./.':
			continue
		if minDepth>0 and depth_index is None:
			continue
		elif minDepth>0 and depth_index is not None:
			depth = genotype_data_ls[depth_index]
			if depth=='.':	#this means depth=0
				depth = 0
			else:
				depth = int(depth)
			if depth<minDepth:	#no read. samtools would still assign ref/ref to this individual
				continue
			#if depth>maxNoOfReads*coverage or depth<minNoOfReads*coverage:	#2011-3-29 skip. coverage too high or too low
			#	continue
			callData['DP'] = depth
		patternSearchResult = diploidGenotypePattern.search(genotype_call)
		genotypeCallInBase = 'NA'
		if patternSearchResult:
			allele1 = alleleNumber2Base[patternSearchResult.group(1)]
			allele2 = alleleNumber2Base[patternSearchResult.group(2)]
			if allele1!='N' and allele2!='N':
				genotypeCallInBase = '%s%s'%(allele1, allele2)
		"""
		if genotype_call=='0/1' or genotype_call =='1/0':	#heterozygous, the latter notation is never used though.
			allele = '%s%s'%(refBase, altBase)
			GL_list = genotype_data_ls[GL_index]
			GL_list = GL_list.split(',')
			GL_list = map(float, GL_list)
			GL = GL_list[1]
			sndHighestGL = max([GL_list[0], GL_list[2]])
			deltaGL = GL-sndHighestGL
			
			AD = genotype_data_ls[format_column_name2index.get('AD')]
			AD = map(int, AD.split(','))
			minorAlleleCoverage = min(AD)
			majorAlleleCoverage = max(AD)
			
			if minorAlleleCoverage<=minorAlleleDepthUpperBoundCoeff*coverage and \
					minorAlleleCoverage>=minorAlleleDepthLowerBoundCoeff*coverage and \
					majorAlleleCoverage<=majorAlleleDepthUpperBoundCoeff*coverage:
				DP4_ratio = float(AD[0])/AD[1]
				allele = '%s%s'%(refBase, altBase)

		elif genotype_call=='./.' or genotype_call=='.|.':	#missing
			allele = 'NA'
		elif genotype_call =='1/1' or genotype_call =='1|1':
			allele = '%s%s'%(altBase, altBase)
		elif genotype_call =='0/0' or genotype_call=='0|0':
			allele = '%s%s'%(refBase, refBase)
		"""
		col_index = sample_id2index.get(individual_name)
		if dataEntryType==1:
			data_row[col_index] = genotypeCallInBase
		else:
			callData['GT'] = genotypeCallInBase
			data_row[col_index] = callData
		if genotypeCallInBase!='NA':
			if genotypeCallInBase not in genotypeCall2Count:
				genotypeCall2Count[genotypeCallInBase] = 0
			genotypeCall2Count[genotypeCallInBase] += 1
	return PassingData(chr=chr, pos=pos, locus_id=locus_id, quality=quality, info_tag2value=info_tag2value, \
					refBase=refBase, altBase=altBase, \
					alleleLs=alleleLs, alleleNumber2Base=alleleNumber2Base, genotypeCall2Count=genotypeCall2Count, data_row=data_row,\
					info=info, format=format, filter=filter, vcf_locus_id=vcf_locus_id,)

class VCFRecord(object):
	"""
	2011-10-20
		
	"""
	def __init__(self, row, col_name2index=None, individual_name2col_index=None, sample_id2index=None,\
				col_index_individual_name_ls=None, minDepth=1, **keywords):
		"""
		2012.1.17
			add argument minDepth
		"""
		self.row = row
		self.col_name2index = col_name2index
		self.individual_name2col_index = individual_name2col_index
		self.sample_id2index = sample_id2index
		self.col_index_individual_name_ls = col_index_individual_name_ls	#list of [column-index, sample-id], column-index >=9
		self.minDepth = minDepth
		self.data_row = []
		self.alleleLs = []	#index 0 is refBase, 1 is first altBase, 2 is 2nd altBase .. 
		self.alleleNumber2Base = {}	#map string type of allele 0,1,2 to the actual base
		self._parse(row)
		
	def _parse(self, row):
		"""
		2012.5.2
			add refBase, altBase to self
		"""
		returnData = parseOneVCFRow(row, self.col_name2index, self.col_index_individual_name_ls, self.sample_id2index, \
					minDepth=self.minDepth, dataEntryType=2)
		
		self.chr = returnData.chr
		self.pos = returnData.pos
		self.locus_id = returnData.locus_id
		self.vcf_locus_id = returnData.vcf_locus_id
		self.info_tag2value = returnData.info_tag2value
		self.data_row = returnData.data_row	#a list of dictionary {'GT': base-call, 'DP': depth}, including one extra sample, the ref
		self.alleleNumber2Base = returnData.alleleNumber2Base
		self.alleleLs = returnData.alleleLs
		
		self.refBase = returnData.refBase
		self.altBase = returnData.altBase
		self.quality = returnData.quality
		self.filter = returnData.filter
		self.info = returnData.info
		self.format = returnData.format
	
	def getAAF(self, useMax=True, defaultValue=0):
		"""
		2012.10.5
			add argument useMax, which determines to return max or min when there are >1 alternative alleles.
			add argument defaultValue
		2011.12.16
		"""
		frequencyLs = self.info_tag2value.get("AF", self.info_tag2value.get("AF1", None))
		if frequencyLs is not None:
			frequencyLs = frequencyLs.split(',')	#if more than one alt alleles, AF is a string separated by ",".
			frequencyLs = map(float ,frequencyLs)
			if useMax:
				frequency = max(frequencyLs)
			else:
				frequency = min(frequencyLs)
		else:
			frequency = defaultValue
		return frequency
	
class VCFFile(object):
	__doc__ = __doc__
	option_default_dict = {
						('inputFname', 0, ): [None, 'i', 1, 'a VCF input file to read in the VCF content.'],\
						('outputFname', 0, ): [None, 'o', 1, 'a VCF output file to hold the the VCF content.'],\
						('openMode', 1, ): ['rb', '', 1, 'rb: bam file. r: sam file.'],\
						('minorAlleleDepthLowerBoundCoeff', 1, float): [1/4., 'M', 1, 'minimum read depth multiplier for an allele to be called (heterozygous or homozygous)', ],\
						('minorAlleleDepthUpperBoundCoeff', 1, float): [3/4., 'A', 1, 'maximum read depth multiplier for the minor allele of a heterozygous call', ],\
						('majorAlleleDepthUpperBoundCoeff', 1, float): [7/8., 'a', 1, 'maximum read depth multiplier for the major allele of het call'],\
						('maxNoOfReadsForGenotypingError', 1, float): [1, 'x', 1, 'if read depth for one allele is below or equal to this number, regarded as genotyping error ', ],\
						('depthUpperBoundCoeff', 1, float): [2, 'm', 1, 'depthUpperBoundCoeff*coverage = maximum read depth for one base'],\
						('depthLowerBoundCoeff', 1, float): [1/4., 'O', 1, 'depthLowerBoundCoeff*coverage = minimum read depth multiplier for one base'],\
						('depthUpperBoundMultiSampleCoeff', 1, float): [3, 'N', 1, 'across n samples, ignore bases where read depth > coverageSum*depthUpperBoundCoeff*multiplier.'],\
						('seqCoverageFname', 0, ): ['', 'q', 1, 'The sequence coverage file. tab/comma-delimited: individual_sequence.id coverage'],\
						('minDepth', 0, float): [0, '', 1, 'minimum depth for a call to regarded as non-missing', ],\
						('defaultCoverage', 1, float): [5, 'd', 1, 'default coverage when coverage is not available for a read group'],\
						("site_type", 1, int): [1, 's', 1, '1: all sites, 2: variants only'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						("sampleStartingColumn", 1, int): [9, '', 1, 'The column index (starting from 0) of the 1st sample genotype'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']
						}
	
	def __init__(self, **keywords):
		"""
		2012.9.5 set default minDepth=0
		2011-9-27
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		
		self.header = None	#2012.5.10 the list of column headers (the header line starting by #CHROM)
		self.headerWithoutHash = None	#2012.5.10 same as self.header except, instead of "#CHROM", it is "CHROM".
		self.sample_id_ls = []
		self.sample_id2index = {}	#the index is the index of its column in the genotype_call_matrix
		self.locus_id_ls = []
		self.locus_id2row_index = {}
		self.locus_id2data = {}
		self.genotype_call_matrix = []
		self.col_name2index = {}	#column index in file
		self.col_index_individual_name_ls = None
		self.individual_name2col_index = {}	#not the matrix column, the column in input file
		self.metaInfoLs = []	#2012.3.28 anything before the "#CHROM" line. each entry is a raw line content, including '\n'
		self.sampleIDHeader = []	#2012.3.20 a list of sample column headers (from sampleStartingColumn)
		
		self.inf = None
		self.reader = None
		self._initializeInput(self.inputFname)
		
		self.outf = None
		self.writer = None
		self._initializeOutput(self.outputFname)
	
	def _initializeInput(self, inputFname=None):
		"""
		2012.5.10
			split out of __init__()
		"""
		if inputFname:
			if inputFname[-3:]=='.gz':
				import gzip
				self.inf = gzip.open(inputFname, 'rb')
			else:
				self.inf = open(inputFname)
			self.reader =csv.reader(self.inf, delimiter='\t')
			self._parseHeader()
	
	def _initializeOutput(self, outputFname=None):
		"""
		2012.9.6 handle gzipped file as well.
		2012.5.10
			split out of __init__()
		"""
		if outputFname:
			if outputFname[-3:]=='.gz':
				import gzip
				self.outf = gzip.open(outputFname, 'wb')
			else:
				self.outf = open(outputFname, 'w')
			self.writer = csv.writer(self.outf, delimiter='\t')
	
	def _getIndividual2ColIndex(self, header, col_name2index=None, sampleStartingColumn=9):
		"""
		2012.10.5 add "_" in front of the function name
		2011-9-27
			called by other function, not meant for public
		"""
		if self.report:
			sys.stderr.write("\t Finding all individuals ...")
		no_of_cols = len(header)
		individual_name2col_index = {}	#individual's column name -> an opened file handler to store genetic data
		col_index_individual_name_ls = []
		
		counter = 0
		for i in xrange(sampleStartingColumn, no_of_cols):
			individualName = header[i].strip()
			col_index = col_name2index.get(individualName)
			if not individualName:	#ignore empty column
				continue
			if individualName[:-4]=='.bam':
				individualCode = individualName[:-4]	#get rid of .bam
			else:
				individualCode = individualName
			individual_name2col_index[individualCode] = col_index
			col_index_individual_name_ls.append([col_index, individualCode])
			counter += 1
		self.individual_name2col_index = individual_name2col_index
		
		col_index_individual_name_ls.sort()	#sorted by column index
		if self.report:
			sys.stderr.write("%s individuals added. Done.\n"%(counter))
		return col_index_individual_name_ls
	
	def get_col_index_individual_name_ls(self):
		"""
		2012.10.5 accessor of self.col_index_individual_name_ls
			a list of [column-index, individual_name].
				The individual_name excludes the trailing ".bam" if there is one.
			The column-index starts from the sampleStartingColumn (=9 usually). 
		"""
		return self.col_index_individual_name_ls
	
	def getIndividualName2ColIndex(self):
		"""
		2012.10.5 public version of _getIndividual2ColIndex()
			a dictionary between samplID and column-index (which column in the VCF file)
			
			difference between sampleID and individual_name
				the latter excludes the trailing ".bam" if there is one. 
			
			The column-index starts from the sampleStartingColumn (=9 usually).
		"""
		return self.individual_name2col_index
	
	def getSampleID2DataIndex(self):
		"""
		2012.10.5 accessor of self.sample_id2index
			It includes "ref" as index=0.
			The index is the index of each sample in the variant data_matrix or data_row of VCFRecord.
			
		"""
		
		return self.sample_id2index
	
	def parseFile(self):
		"""
		#2012.8.8 bug fix, same locus appearing >1 times
		2011.11.30
			bugfix: genotype_call could be ./.(missing) and depth=. (missing) in "./.:.:.:.:.".
		2011-9-27
			based on discoverFromVCF() from vervet/src/GenotypeCallByCoverage.py
			
			It reads all genotype calls into the matrix so its memory footprint might explode.
		"""
		sys.stderr.write("Parsing %s ... \n"%(os.path.basename(self.inputFname),))
		
		locus_id2row_index = self.locus_id2row_index
		
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
		
		individual_name2col_index = self.individual_name2col_index
		col_name2index = self.col_name2index
		counter = 0
		real_counter = 0
		
		for row in self.reader:
			counter += 1
			returnData = parseOneVCFRow(row, self.col_name2index, self.col_index_individual_name_ls, self.sample_id2index, \
					minDepth=self.minDepth, dataEntryType=1)
			
			genotypeCall2Count = returnData.genotypeCall2Count
			data_row = returnData.data_row
			
			current_locus = returnData.locus_id
			if current_locus not in self.locus_id2data:
				self.locus_id2data[current_locus] = returnData.info_tag2value
				self.locus_id2data[current_locus]['REF'] = returnData.refBase
				self.locus_id2data[current_locus]['ALT'] = returnData.altBase
			
			if len(genotypeCall2Count)>self.site_type-1:	#whether polymorphic across samples or all sites in vcf
				real_counter += 1
				if current_locus in locus_id2row_index:	#2012.8.8 bug fix, same locus appearing >1 times
					sys.stderr.write("WARNING: Locus %s at line %s already in locus_id2row_index with row_index = %s. Skip\n"%\
									(repr(current_locus), counter, locus_id2row_index.get(current_locus)))
				else:
					locus_id2row_index[current_locus] = len(locus_id2row_index)
					self.locus_id_ls.append(current_locus)
					self.genotype_call_matrix.append(data_row)
					real_counter += 1
			
			if counter%2000==0 and self.report:
				sys.stderr.write("%s\t%s\t%s"%("\x08"*80, counter, real_counter))
		sys.stderr.write("%s loci X %s samples.\n"%(len(self.locus_id_ls), len(self.sample_id_ls)))
	
	def _parseHeader(self):
		"""
		2012.3.28
			add all header content into self.metaInfoLs
				except the last header line, which goes into self.sampleIDHeader
		2011-11-2
			this function is run inside __init__()
		"""
		self.metaInfoLs = []	#2012.3.28 anything before the "#CHROM" line. each entry is a raw line content, including '\n'
		self.sampleIDHeader = []	#2012.3.20 a list of column headers (#CHROM)
		self.sample_id_ls = []
		
		self.sample_id2index['ref'] = 0	#ref is at column 0. "ref" must not be equal to any read_group.
		self.sample_id_ls.append('ref')
		read_group2coverage = {}	#2011-9-2
		
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
		
		counter = 0
		real_counter = 0
		
		for line in self.inf:
			if line[:6]=='#CHROM':
				line = line.strip()	#get rid of the trailing \n
				row = line.split('\t')
				self.sampleIDHeader = row[self.sampleStartingColumn:]
				self.header = row[:]
				self.headerWithoutHash= row[:]
				self.headerWithoutHash[0] = 'CHROM'	#discard the #
				self.col_name2index = getColName2IndexFromHeader(self.headerWithoutHash, skipEmptyColumn=True)
				self.col_index_individual_name_ls = self._getIndividual2ColIndex(self.headerWithoutHash, self.col_name2index)
				for individual_col_index, individual_name in self.col_index_individual_name_ls:
					read_group = individual_name.strip()
					if read_group not in self.sample_id2index:
						self.sample_id2index[read_group] = len(self.sample_id2index)
						self.sample_id_ls.append(read_group)
				break	# "#CHROM" is the last line of the self.headerWithoutHash
			elif line[0]=='#':	#2011-3-4
				self.metaInfoLs.append(line)
				#continue
			else:	#leave everything for parseFile or parseIter
				break
	
	def parseIter(self):
		"""
		2011-11-2
			an iterator over each line (call data) in VCF
		"""
		for row in self.reader:
			vcfRecord = VCFRecord(row, col_name2index=self.col_name2index, \
						individual_name2col_index=self.individual_name2col_index, sample_id2index=self.sample_id2index,\
						col_index_individual_name_ls=self.col_index_individual_name_ls,\
						minDepth=self.minDepth)
			yield vcfRecord
	
	def __iter__(self):
		"""
		2012.8.25 make itself an iterator
		"""
		return self
	
	def next(self):
		"""
		2012.8.25 make itself an iterator
		"""
		try:
			row = self.reader.next()
			if not row:
				raise StopIteration
			vcfRecord = VCFRecord(row, col_name2index=self.col_name2index, \
						individual_name2col_index=self.individual_name2col_index, sample_id2index=self.sample_id2index,\
						col_index_individual_name_ls=self.col_index_individual_name_ls,\
						minDepth=self.minDepth)
			return vcfRecord
		except:
			raise StopIteration
	
	def _countHomoHetCallsForEachSampleFromVCF(self):
		"""
		2011-11-2
			given a VCF file, count the number of homo-ref, homo-alt, het calls
			
		"""
		sys.stderr.write("Count the number of homozygous-ref/alt & het.\n")
		sampleID2data = {}	#key is sampleID, value is a list of 3 numbers. 'NoOfHomoRef', 'NoOfHomoAlt', 'NoOfHet'
		
		no_of_total = 0.
		minStart = None
		for vcfRecord in self.parseIter():
			locus_id = vcfRecord.locus_id
			chr = vcfRecord.chr
			pos = vcfRecord.pos
			pos = int(pos)
			refBase = vcfRecord.data_row[0].get("GT")[0]
			
			for sample_id, sample_index in self.sample_id2index.iteritems():
				if sample_id=='ref':	#ignore the reference
					continue
				if sample_id not in sampleID2data:
					sampleID2data[sample_id] = [0, 0, 0]
				if not vcfRecord.data_row[sample_index]:	#None for this sample
					continue
				callForThisSample = vcfRecord.data_row[sample_index].get('GT')
				if not callForThisSample or callForThisSample=='NA':
					continue
				if callForThisSample[0]==refBase and callForThisSample[1]==refBase:
					#homozygous reference allele
					sampleID2data[sample_id][0]+=1
				elif callForThisSample[0]==callForThisSample[1] and callForThisSample[0]!=refBase:
					#homozygous alternative allele
					sampleID2data[sample_id][1]+=1
				elif callForThisSample[0]!=callForThisSample[1]:
					sampleID2data[sample_id][2]+=1
		
		sampleIDLs = sampleID2data.keys()
		sampleIDLs.sort()
		for sampleID in sampleIDLs:
			count_data = sampleID2data.get(sampleID)
			noOfHomoRef, noOfHomoAlt, noOfHet = count_data[:3]
			no_of_calls = float(sum(count_data))
			if no_of_calls>0:
				fractionOfHomoAlt = noOfHomoAlt/no_of_calls
				fractionOfHet = noOfHet/no_of_calls
			else:
				fractionOfHomoAlt = -1
				fractionOfHet = -1
		sys.stderr.write("Done.\n")
	
	def getSampleIDList(self):
		"""
		2012.7.18
			self.sample_id_ls includes a "ref" in index=0 cell. this function helps to reduce confusion by returning
			a list of sample IDs without "ref"
		"""
		return self.sample_id_ls[1:]
	
	def getNoOfLoci(self):
		"""
		2012.9.25
			fast way to get the number of loci, no VCF parsing.
		"""
		if self.inf is None:
			return None
		
		if len(self.sample_id_ls)==0 and len(self.sampleIDHeader)==0:
			self._parseHeader()
		no_of_loci = 0
		for line in self.inf:
			no_of_loci += 1
		self._resetInput()
		
		return no_of_loci
	
	def _resetInput(self):
		"""
		2012.10.5
			run only to reset this VCFFile object and the read cursor is at the beginning of all VCF records.
			split out of getNoOfLoci()
		"""
		#reset the inf
		self.inf.seek(0)
		#run this to move the cursor to where the VCF records start
		self._parseHeader()
	
	def writeMetaAndHeader(self, outf=None):
		"""
		2012.5.10
			vcfFile = VCFFile(outputFname=...)
			vcfFile.metaInfoLs = oldVCFFile.metaInfoLs
			vcfFile.header = [...]
			vcfFile.writeMetaAndHeader()
			
			write out the self.header and self.metaInfoLs
		"""
		if outf is None:
			outf = self.outf
			writer = self.writer
		else:
			writer = csv.writer(outf, delimiter='\t')
		if outf is None:
			sys.stderr.write("Warning: outf is None. nothing is written out.\n")
		else:
			for metaInfo in self.metaInfoLs:
				outf.write(metaInfo)
			writer.writerow(self.header)
	
	def writeVCFRecord(self, vcfRecord=None, outf=None):
		"""
		2012.5.10
		"""
		if outf is None:
			outf = self.outf
			writer = self.writer
		else:
			writer = csv.writer(outf, delimiter='\t')
		if outf is None:
			sys.stderr.write("Warning: outf is None. nothing is written out.\n")
		else:
			self.locus_id_ls.append(vcfRecord.vcf_locus_id)#2012.9.6
			data_row = [vcfRecord.chr, vcfRecord.pos, vcfRecord.vcf_locus_id, vcfRecord.refBase, vcfRecord.altBase, \
					vcfRecord.quality, vcfRecord.filter,\
					vcfRecord.info, vcfRecord.format] + vcfRecord.row[self.sampleStartingColumn:]
			writer.writerow(data_row)
	
	def close(self,):
		"""
		2012.5.10
		"""
		if self.inf:
			del self.inf, self.reader
		if self.outf:
			del self.outf, self.writer
	
	def getLocus2AlternativeAlleleFrequency(self):
		"""
		2012.10.5
			
			adapted from VariantDiscovery.getContig2Locus2Frequency()
		
		"""
		locus2frequency = {}
		counter = 0
		real_counter = 0
		for vcfRecord in self:
			current_locus = (vcfRecord.chr, vcfRecord.pos)
			AAF = vcfRecord.getAAF()
			locus2frequency[current_locus] = AAF
		
		self._resetInput()
		return locus2frequency