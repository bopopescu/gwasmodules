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
						('markMissingStatFname', 0, ): ["", '', 1, 'only required for run_type 4, used to store genotypes that have been marked missing'],\
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
		2013.07.24 return individualIDList as well
		2013.1.29
		"""
		sys.stderr.write("Getting individualID2Index from tfam file %s ..."%(tfamFname))
		individualID2Index = {}
		individualIDList = []
		reader = MatrixFile(inputFname=tfamFname)
		counter = 0
		for row in reader:
			individualID = row[1]
			individualID2Index[individualID] = len(individualID2Index)
			individualIDList.append(individualID)
			counter += 1
		del reader
		sys.stderr.write(" %s individuals.\n"%(len(individualID2Index)))
		return PassingData(individualID2Index=individualID2Index, individualIDList=individualIDList)
	
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
	
	def markGenotypeMissingIfInvolvedInMendelError(self, row=None, locus_id2individual_index_ls=None, \
												individual_index2no_of_genotype_marked_missing=None):
		"""
		2013.07.24 bugfix
		2013.1.29 
			need to read the tfam file to figure out which column one individual is on.
			Starting from 5th column, data is genotype of individuals (order is recorded in tfam).
				Each individual's genotype occupies two columns (diploid). 
		"""
		chromosome, snp_id, genetic_distace, physical_distance = row[:4]
		individual_index_ls = locus_id2individual_index_ls.get(snp_id)
		
		if individual_index_ls:
			for index in individual_index_ls:
				row[4+index*2] = 0	#mark it missing
				row[4+index*2+1] = 0
				if individual_index2no_of_genotype_marked_missing is not None:
					if index not in individual_index2no_of_genotype_marked_missing:
						individual_index2no_of_genotype_marked_missing[index] = 0
					individual_index2no_of_genotype_marked_missing[index] += 1
		return row
	
	def outputGenotypeMarkedMissingStat(self, outputFname=None, \
									individual_index2no_of_genotype_marked_missing=None,\
									individualIDList=None):
		"""
		2013.07.24
		"""
		if outputFname and individual_index2no_of_genotype_marked_missing is not None:
			writer = MatrixFile(inputFname=outputFname, openMode='w', delimiter='\t')
			header = ["individualID", "noOfGenotypesMarkedMissing"]
			writer.writeHeader(header)
			for individual_index, no_of_genotype_marked_missing in individual_index2no_of_genotype_marked_missing.iteritems():
				individual_id = individualIDList[individual_index]
				writer.writerow([individual_id, no_of_genotype_marked_missing])
			writer.close()
		
	
	def run(self):
		"""
		2013.07.24
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		#inf = utils.openGzipFile(self.inputFname)
		reader = MatrixFile(inputFname=self.inputFname)
		writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		counter = 0
		if self.run_type==4:	#2013.2.1
			tfamIndividualData = self.getIndividualID2IndexFromTFAMFile(tfamFname=self.tfamFname)
			individualID2Index = tfamIndividualData.individualID2Index
			individualIDList = tfamIndividualData.individualIDList
			locus_id2individual_index_ls = self.getMendelErrorIndividualLocusData(mendelErrorFname=self.mendelErrorFname, \
												individualID2Index=individualID2Index)
			individual_index2no_of_genotype_marked_missing = {}
		else:
			individualID2Index = None
			individualIDList = None
			locus_id2individual_index_ls = None
			individual_index2no_of_genotype_marked_missing = None
		for row in reader:
			if self.run_type==2:
				new_row = self.processRow_ChangeChromosomeIDToX(row)
			elif self.run_type==3:
				new_row = self.processRow_addPositionStartBase(row)
			elif self.run_type==4:
				new_row = self.markGenotypeMissingIfInvolvedInMendelError(row=row, \
											locus_id2individual_index_ls=locus_id2individual_index_ls,\
											individual_index2no_of_genotype_marked_missing=individual_index2no_of_genotype_marked_missing)
				
			else:
				new_row = self.processRow(row)
			writer.writerow(new_row)
			counter += 1
		sys.stderr.write("%s lines modified.\n"%(counter))
		
		del reader
		del writer
		self.outputGenotypeMarkedMissingStat(outputFname=self.markMissingStatFname, \
								individual_index2no_of_genotype_marked_missing=individual_index2no_of_genotype_marked_missing, \
								individualIDList=individualIDList)
	
	
if __name__ == '__main__':
	main_class = ModifyTPED
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()