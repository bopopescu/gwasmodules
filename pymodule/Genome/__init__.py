"""
2008-10-01
	module related to Genome
"""
import os,sys

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

class fasta_block_iterator:
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
chr_pattern = re.compile(r'([a-zA-Z]+\d+)[._\-:]*')	#the last - has special meaning in [] when it's not the last character. 
contig_id_pattern = re.compile(r'Contig(\d+)[._\-:]*')

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
	chr_pattern_sr = chr_pattern.search(filename)
	if chr_pattern_sr:
		chr = chr_pattern_sr.group(1)
	else:
		chr = os.path.splitext(os.path.split(filename)[1])[0]
	return chr