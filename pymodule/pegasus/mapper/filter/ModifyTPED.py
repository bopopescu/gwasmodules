#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i input.tped.gz -o /tmp/output.tped

Description:
	2012.7.20
		all input files could be gzipped or not.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, MatrixFile
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper

class ModifyTPED(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.update({
						('run_type', 1, int): [1, 'y', 1, 'which run type. 1: modify snp_id (2nd-column) = chr_phyiscalPosition,\n\
		2: snp_id=chr_physicalPosition (original data), set chromosome = X (chromosome X, for sex check by plink), pos += positionStartBase. \n\
			assuming input is on X already,\n\
		3: snp_id=chr_physicalPosition (original data), chromosome (1st column) = newChr, pos += positionStartBase,\n\
		4: mark genotype calls involved in mendel inconsistency as missing, requiring mendelErrorFname.\n', ],\
						('mendelErrorFname', 0, ): ["", '', 1, 'the plink mendel error output file (the full version). looks like\
 FID       KID  CHR                 SNP   CODE                 ERROR\n\
   1   1996027    0   Contig1149_104181      2      T/T x T/T -> G/T'],\
						('tfamFname', 0, ): ["", '', 1, 'the plink tfam file used to figure out which individual is on which column'],\
						('newChr', 0, ): ["1", 'n', 1, 'the new chromosome for the TPED file (replace the old one)'],\
						('positionStartBase', 0, int): [0, 'p', 1, 'the number to be added to position of every SNP'],\
						})
	def __init__(self,  **keywords):
		"""
		"""
		AbstractMapper.__init__(self, **keywords)
	
	def processRow(self, row):
		"""
		2012.8.9
		"""
		chromosome, snp_id, genetic_distace, physical_distance = row[:4]
		#chromosome = Genome.getContigIDFromFname(chromosome)	# 2012.8.16 getting rid of the string part of chromosome ID doesn't help.
			#   non-human chromosome numbers would still be regarded as 0 by plink.
		snp_id = '%s_%s'%(chromosome, physical_distance)
		new_row = [chromosome, snp_id, genetic_distace, physical_distance] + row[4:]
		return new_row
	
	def processRow_ChangeChromosomeIDToX(self, row):
		"""
		2012.8.9
		"""
		chromosome, snp_id, genetic_distace, physical_distance = row[:4]
		snp_id = '%s_%s'%(chromosome, physical_distance)	#the snp_id is the original contig & position
		#chromosome = self.newChr	#new chromosome, new position
		physical_distance = int(physical_distance) + self.positionStartBase
		chromosome = "X"
		new_row = [chromosome, snp_id, genetic_distace, physical_distance] + row[4:]
		return new_row
	
	def processRow_addPositionStartBase(self, row):
		"""
		2012.8.9
		"""
		chromosome, snp_id, genetic_distace, physical_distance = row[:4]
		snp_id = '%s_%s'%(chromosome, physical_distance)	#the snp_id is the original contig & position
		chromosome = self.newChr	#new chromosome, new position
		physical_distance = int(physical_distance) + self.positionStartBase
		#chromosome = Genome.getContigIDFromFname(chromosome)
		new_row = [chromosome, snp_id, genetic_distace, physical_distance] + row[4:]
		return new_row
	
	def getIndividualID2IndexFromTFAMFile(self, tfamFname=None):
		"""
		2013.1.29
		"""
		sys.stderr.write("Getting individualID2Index from tfam file %s ..."%(tfamFname))
		individualID2Index = {}
		reader = MatrixFile(inputFname=tfamFname)
		index = 0
		for row in reader:
			individualID = row[1]
			individualID2Index[individualID] = len(individualID2Index)
			index += 1
		del reader
		sys.stderr.write(" %s individuals.\n"%(len(individualID2Index)))
		return individualID2Index
	
	def getMendelErrorIndividualLocusData(self, mendelErrorFname=None, individualID2Index=None):
		"""
		2013.1.29
		
		"""
		sys.stderr.write("Getting data on loci involved in mendel-errors from %s ..."%(mendelErrorFname))
		locus_id2individual_index_ls = {}
		#inf = utils.openGzipFile(mendelErrorFname, 'r')
		reader = MatrixFile(inputFname=mendelErrorFname)
		#header = reader.next()
		reader.constructColName2IndexFromHeader()
		counter = 0
		for row in reader:
			individual_id = row[reader.getColIndexGivenColHeader('KID')]
			if individual_id in individualID2Index:
				index =individualID2Index.get(individual_id)
			else:
				sys.stderr.write("Individual %s not in individualID2Index.\n"%(individual_id))
				sys.exit(3)
			snp_id = row[3]
			if snp_id not in locus_id2individual_index_ls:
				locus_id2individual_index_ls[snp_id] = []
			locus_id2individual_index_ls[snp_id].append(index)
			counter += 1
		del reader
		sys.stderr.write(" %s calls of %s loci, involved in mendel errors.\n"%\
						(counter, len(locus_id2individual_index_ls)))
		return locus_id2individual_index_ls
	
	def markGenotypeMissingIfInvolvedInMendelError(self, row=None, locus_id2individual_index_ls=None):
		"""
		2013.1.29 
			need to read the tfam file to figure out which column one individual is on.
			Starting from 5th column, data is genotype of individuals (order is recorded in tfam).
				Each individual's genotype occupies two columns (diploid). 
		"""
		chromosome, snp_id, genetic_distace, physical_distance = row[:4]
		individual_index_ls = locus_id2individual_index_ls.get(snp_id)
		if individual_index_ls:
			for index in individual_index_ls:
				row[index*2] = 0	#mark it missing
				row[index*2+1] = 0
		return row
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		#inf = utils.openGzipFile(self.inputFname)
		reader = MatrixFile(inputFname=self.inputFname)
		writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		counter = 0
		if self.run_type==4:	#2013.2.1
			individualID2Index = self.getIndividualID2IndexFromTFAMFile(tfamFname=self.tfamFname)
			locus_id2individual_index_ls = self.getMendelErrorIndividualLocusData(mendelErrorFname=self.mendelErrorFname, \
												individualID2Index=individualID2Index)
		else:
			individualID2Index = None
			locus_id2individual_index_ls = None
		for row in reader:
			if self.run_type==2:
				new_row = self.processRow_ChangeChromosomeIDToX(row)
			elif self.run_type==3:
				new_row = self.processRow_addPositionStartBase(row)
			elif self.run_type==4:
				new_row = self.markGenotypeMissingIfInvolvedInMendelError(row=row, \
														locus_id2individual_index_ls=locus_id2individual_index_ls)
			else:
				new_row = self.processRow(row)
			writer.writerow(new_row)
			counter += 1
			
		del reader
		del writer
		sys.stderr.write("%s lines modified.\n"%(counter))

if __name__ == '__main__':
	main_class = ModifyTPED
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()