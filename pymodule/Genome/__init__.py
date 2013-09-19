"""
2008-10-01
	module related to Genome
"""

import os,sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from pymodule.utils import PassingData

class GeneModel(object):
	def __init__(self, **keywords):
		"""
		2008-10-01
			moved from transfac.src.GenomeDB in order for the return of GenomeDB.get_gene_id2model() to be pickled independent of transfac.src
		2008-10-01
			a class to hold all stuff related to a gene (Gene+EntrezgeneMapping)
			it's hierarchical. Its gene_commentaries contains also GeneModel.
		"""
		for argument_key, argument_value in keywords.iteritems():
			setattr(self, argument_key, argument_value)
		if not hasattr(self, 'gene_commentaries'):
			self.gene_commentaries = []

class fasta_block_iterator(object):
	'''
	2011-1-5 moved from Transfac.src.transfacdb
	09-13-05
		fasta format iterator
		a little bit tricky, '>', the block starter is used as a tokenizer
	2006-09-01
		it seems 'for line in self.inf' doesn't work on hpc-cmb.
		check https://dl403k-1.cmb.usc.edu/log/hpc-cmb
	'''
	def __init__(self, inf):
		self.inf = inf
		self.block = ''
		self.previous_line = ''
	def __iter__(self):
		return self
	def next(self):
		self.read()
		return self.block
	def read(self):
		self.block = self.previous_line	#don't forget the starting line
		line = self.inf.readline()
		while(line):
			if line[0]=='>':
				self.previous_line = line
				if self.block:	#not the first time
					break
				else:	#first time to read the file, block is still empty
					self.block += line
			else:
				self.block += line
			line = self.inf.readline()
		if self.block==self.previous_line:	#nothing new into the block
			raise StopIteration

class LargeFastaFileTraverse:
	"""
	2010-1-5
	"""
	def __init__(self):
		pass
	
	def traverse(self, input_dir, headerFunctor=None, seqlineFunctor=None):
		"""
		2011-1-5
			orignal idea is in transfac/src/chromosome_fasta2db.py.
			Each fasta file in input_dir could contain many fasta blocks.
		"""
		import re
		p_chromosome = re.compile(r'chromosome (\w+)[,\n\r]?')	#the last ? means [,\n\r] is optional
		files = os.listdir(input_dir)
		no_of_total_files = len(files)
		
		for i in xrange(no_of_total_files):
			fname = files[i]
			input_fname = os.path.join(input_dir, fname)
			sys.stderr.write("\t %s/%s %s ..."%(i+1, no_of_total_files, fname))
			if fname[-2:]=='gz':
				import gzip
				inf = gzip.open(input_fname)
			else:
				inf = open(input_fname)
			line = inf.readline()
			new_fasta_block = 1	#'line' is not enough to stop the 'while' loop. after the file reading is exhausted by "for line in inf:", 'line' still contains the stuff from the last line.
			no_of_fasta_blocks = 0
			while line and new_fasta_block:
				new_fasta_block = 0	#set it to 0, assuming only one fasta block, change upon new fasta block
				if line[0]!='>':	#not fasta block header
					for line in inf:	#exhaust this fasta block as it's not what's wanted.
						if line[0]=='>':
							new_fasta_block = 1
							break	#start from while again
					continue
				
				headerFunctor(line)
				"""
				# 2010-1-5 an example header function which parses the title of the fasta block 
				#possible header lines:
				#>gi|51511461|ref|NC_000001.8|NC_000001 Homo sapiens chromosome 1, complete sequence
				#>gi|186497660|ref|NC_003070.6| Arabidopsis thaliana chromosome 1, complete sequence
				#>gi|26556996|ref|NC_001284.2| Arabidopsis thaliana mitochondrion, complete genome
				#>gi|115442598|ref|NC_008394.1| Oryza sativa (japonica cultivar-group) genomic DNA, chromosome 1
				header = line[1:-1]	#discard '>' and '\n'
				header = header.split('|')
				
				if tax_id is None:
					_tax_id = FigureOutTaxID_ins.returnTaxIDGivenSentence(header[4])
					if not _tax_id:
						_tax_id = 'null'
				else:
					_tax_id = tax_id
				
				if p_chromosome.search(header[4]) is not None:
					chromosome = p_chromosome.search(header[4]).groups()[0]
				elif header[4].find('mitochondrion')!=-1:
					chromosome = 'mitochondrion'
				elif header[4].find('chloroplast')!=-1:
					chromosome = 'chloroplast'
				else:	#something else, take the whole before ','
					chromosome = header[4].split(',')[0]
				
				outf.write(">chr%s\n"%(chromosome))
				"""
				for line in inf:
					if line[0]=='>':
						
						new_fasta_block = 1
						break	#start from while again
					else:
						seqlineFunctor(line)
			sys.stderr.write("\n")

import re
chr_pattern = re.compile(r'([a-zA-Z]+[\dXY]+)[._\-:]*')	#the last - has special meaning in [] when it's not the last character. 
#contig_id_pattern = re.compile(r'Contig(\d+)[._\-:]*')	#2013.3.19
#contig_id_pattern = re.compile(r'Scaffold(\d+)[._\-:]*')	#2013.3.19 new vervet ref is scaffold-based.
contig_id_pattern = re.compile(r'[CcS][a-zA-Z]+([\dXY]+)[._\-:]*')	#2013.07.16 new chromosome (Chlorocebus aethiops) vervet ref, add X, Y
chr_start_stop_pattern = re.compile(r'([a-zA-Z]+[\dXY]+)_(\d+)_(\d+)[._\-:]*')	#the last - has special meaning in [] when it's not the last character. 

def getContigIDFromFname(filename):
	"""
	2012.7.14 copied from  pymodule.pegasus.AbstractNGSWorkflow
	2011-10-20
		
		If filename is like .../Contig0.filter_by_vcftools.recode.vcf.gz,
			It returns "0", excluding the "Contig".
			If you want "Contig" included, use getChrIDFromFname().
		If search fails, it returns the prefix in the basename of filename.
	"""
	contig_id_pattern_sr = contig_id_pattern.search(filename)
	if contig_id_pattern_sr:
		contig_id = contig_id_pattern_sr.group(1)
	else:
		contig_id = os.path.splitext(os.path.split(filename)[1])[0]
	return contig_id

def getChrFromFname(filename):
	"""
	2012.7.14 copied from  pymodule.pegasus.AbstractNGSWorkflow
	2011-10-20
		filename example: Contig0.filter_by_vcftools.recode.vcf.gz
			It returns "Contig0".
			If you want just "0", use getContigIDFromFname().
		If search fails, it returns the prefix in the basename of filename.
	"""
	searchResult = chr_pattern.search(filename)
	if searchResult:
		chromosome = searchResult.group(1)
	else:
		chromosome = os.path.splitext(os.path.split(filename)[1])[0]
	return chromosome

def parseChrStartStopFromFilename(filename=None, chr2size=None, defaultChromosomeSize=10000000000):
	"""
	2013.09.18
		#10000000000 is used when filename contains data from a whole chromosome and chr2size is not available or not containing chromosome
		make it very big so that it could be intersected with any interval from any chromosome.
	"""
	searchResult = chr_start_stop_pattern.search(filename)
	if searchResult:
		chromosome = searchResult.group(1)
		start = int(searchResult.group(2))
		stop = int(searchResult.group(3))
	else:	#try
		chromosome = getChrFromFname(filename=filename)
		start =1
		if chr2size is not None:
			stop = chr2size.get(chromosome, defaultChromosomeSize)
		else:
			stop = defaultChromosomeSize
	return PassingData(chromosome=chromosome, start=start, stop=stop)

class IntervalData(PassingData):
	"""
	2013.09.05 a more dynamic way to store intervals
		required keywords: chr, chromosome, chromosomeSize, start, stop
		overlapStart and overlapStop are optional.
	"""
	def __init__(self, chr=None, chromosome=None, chromosomeSize=None, start=None, stop=None, 
				overlapStart=None, overlapStop=None, **keywords):
		PassingData.__init__(self, chr=chr, chromosome=chromosome, chromosomeSize=chromosomeSize, \
							start=start, stop=stop, 
				overlapStart=overlapStart, overlapStop=overlapStop, **keywords)
		if not hasattr(self, 'file'):
			self.file = None
		if not hasattr(self, 'jobLs'):
			self.jobLs = []
		if self.chr is None and self.chromosome:
			self.chr = self.chromosome
		elif self.chr and self.chromosome is None:
			self.chromosome = self.chr
		
		if self.overlapStart is None:
			self.overlapStart = self.start
		
		if self.overlapStop is None:
			self.overlapStop = self.stop
		
		self.subIntervalLs = []
		self.subIntervalLs.append((self.overlapStart, self.overlapStop))
	
	def unionOneIntervalData(self, otherIntervalData=None):
		"""
		"""
		if self.chromosome!=otherIntervalData.chromosome:
			sys.stderr.write("Error: this interval %s is trying to union an interval from a different chromosome %s.\n"%\
							(self.interval, otherIntervalData.interval))
			raise
		self.subIntervalLs.append((otherIntervalData.overlapStart, otherIntervalData.overlapStop))
		self.subIntervalLs.sort()
		self.start = min(self.start, otherIntervalData.start)
		self.stop = max(self.stop, otherIntervalData.stop)
		self.overlapStart = min(self.overlapStart, otherIntervalData.overlapStart)
		self.overlapStop = max(self.overlapStop, otherIntervalData.overlapStop)
		
		
	@property
	def span(self):
		deltaLs = [a[1]-a[0]+1 for a in self.subIntervalLs]
		return sum(deltaLs)
	
	@property
	def overlapIntervalFileBasenameSignature(self):
		return '%s_%s_%s'%(self.chromosome, self.overlapStart, self.overlapStop)
	
	@property
	def intervalFileBasenameSignature(self):
		return '%s_%s_%s'%(self.chromosome, self.start, self.stop)
		
		
	@property
	def overlapInterval(self):
		return "%s:%s-%s"%(self.chromosome, self.overlapStart, self.overlapStop)
	@property
	def interval(self):
		return "%s:%s-%s"%(self.chromosome, self.start, self.stop)
	@property
	def parentInterval(self):
		return self.interval