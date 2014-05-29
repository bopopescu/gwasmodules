#!/usr/bin/env python
"""
2012.10.8
	module related to Fasta file
"""
import os,sys
from Bio import SeqIO
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from pymodule import ProcessOptions
from pymodule.utils import dict_map, importNumericArray, figureOutDelimiter, PassingData, getColName2IndexFromHeader

class FastaFile(object):
	__doc__ = __doc__
	option_default_dict = {
						('inputFname', 0, ): [None, 'i', 1, 'a VCF input file to read in the VCF content.'],\
						('outputFname', 0, ): [None, 'o', 1, 'a VCF output file to hold the the VCF content.'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']
						}
	
	def __init__(self, **keywords):
		"""
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		self.inf = None
		self._initializeInput(inputFname=self.inputFname)
		self._seqTitle2sequence = None
		self._seqTitleList = None	#2013.07.08
		
		self.outf = None
	
	def _initializeInput(self, inputFname=None):
		"""
		"""
		if inputFname:
			if inputFname[-3:]=='.gz':
				import gzip
				self.inf = gzip.open(inputFname, 'rb')
			else:
				self.inf = open(inputFname)
		
	def close(self,):
		"""
		2012.5.10
		"""
		if self.inf:
			del self.inf
			
	def __iter__(self):
		"""
		2012.10.8 make itself an iterator
		"""
		return self 	#or should it be SeqIO.parse(self.inf, "fasta") ?
	
	def next(self):
		"""
		2012.10.8 make itself an iterator
		"""
		try:
			return SeqIO.parse(self.inf, "fasta").next()
		except:
			raise StopIteration
	
	@property
	def seqTitleList(self):
		"""
		2013.07.08
		"""
		if self._seqTitleList is None:
			self.seqTitleList = ()
		return self._seqTitleList
	
	@seqTitleList.setter
	def seqTitleList(self, argument_ls=[]):
		"""
		2013.07.08
		"""
		sys.stderr.write("Getting seqTitleList from %s ..."%(self.inputFname))
		if self.inf.tell()!=0:
			self.inf.seek(0)
		self._seqTitleList = []
		for record in SeqIO.parse(self.inf, "fasta"):
			seqTitle = record.id
			self._seqTitleList.append(seqTitle)
		sys.stderr.write("%s sequences.\n"%(len(self._seqTitleList)))
		
	@property
	def seqTitle2sequence(self):
		"""
		2012.10.8
		"""
		if self._seqTitle2sequence is None:
			self.seqTitle2sequence = ()
		return self._seqTitle2sequence
	
	@seqTitle2sequence.setter
	def seqTitle2sequence(self, argument_ls=[]):
		"""
		2012.10.8
		"""
		sys.stderr.write("Getting seqTitle2sequence from %s ..."%(self.inputFname))
		if self.inf.tell()!=0:
			self.inf.seek(0)
		self._seqTitle2sequence = {}
		for record in SeqIO.parse(self.inf, "fasta"):
			seqTitle = record.id
			seq = record.seq
			self._seqTitle2sequence[seqTitle] = seq
		sys.stderr.write("%s sequences.\n"%(len(self._seqTitle2sequence)))

	
	def getSequence(self, seqTitle=None, start=None, stop=None):
		"""
		2012.10.8
			coordinates are 1-based and full closure.
			example: start=1, stop=10 => [1,10]
		"""
		sequence = self.seqTitle2sequence.get(seqTitle)
		if sequence:
			return sequence[start-1:stop]
		else:
			return None