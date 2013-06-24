#!/usr/bin/env python
"""
Examples:
	%s -s 1 -o /tmp/Contig685_largeSiteGap.tsv -V 5000  -W distanceToNextSite  /tmp/Contig685_siteGap.tsv
	
	# 2012.8.6 select rows based on the fraction of heterozygotes per individual [0.2, 0.8].
	%s -s 1 -o /tmp/hetPerMonkey_hist.png
		-W NoOfHet_by_NoOfTotal -V 0.2 -x 0.8 /tmp/homoHetCountPerSample.tsv
	

Description:
	2013.06.06
		This program converts Beagle genotype file (could be >1) into VCF.
		If one sample's genotype appears in >1 beagle genotype files, genotype in the last file will be used.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])


#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:	   #64bit
#	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
#	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
#else:   #32bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import random
import networkx as nx
from pymodule import ProcessOptions
from pymodule.utils import PassingData
from pymodule.yhio.AbstractMatrixFileWalker import AbstractMatrixFileWalker
from pymodule.yhio.BeagleGenotypeFile import BeagleGenotypeFile
from pymodule.yhio.VCFFile import VCFFile
from pymodule.statistics import NumberContainer, DiscreteProbabilityMassContainer

class CombinePhasedBeagleOutputsIntoVCF(AbstractMatrixFileWalker):
	__doc__ = __doc__
	
	option_default_dict = AbstractMatrixFileWalker.option_default_dict
	option_default_dict.update({
			('replicateIndividualTag', 0, ): ['copy', '', 1, 'the tag that separates the true ID and its replicate count'],\
			('originalVCFFname', 1, ): ['', '', 1, 'original VCF file on which both Beagle phased output and output VCF will be based. \n\
	The output VCF will be same as originalVCFFname, except GT field, to be replaced by phased genotypes from Beagle-phased files'],\
			})
	
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractMatrixFileWalker.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		#a map from one sample to specific beagle file
		self.sampleID2BeagleFile = None
	
	def setup(self, **keywords):
		"""
		2012.10.15
			run before anything is run
		"""
		#2013.05.30 comment out AbstractMatrixFileWalker.setup() to open the output file differently
		#AbstractMatrixFileWalker.setup(self, **keywords)
		self.writer = VCFFile(inputFname=self.outputFname, openMode='w')
		self.reader = VCFFile(inputFname=self.originalVCFFname, openMode='r')
		self.writer.metaInfoLs = self.reader.metaInfoLs
		self.writer.header = self.reader.header
		self.writer.writeMetaAndHeader()
		
		# read all the Beagle files
		sampleID2BeagleFile = {}
		for inputFname in self.inputFnameLs:
			beagleFile = BeagleGenotypeFile(inputFname=inputFname)
			beagleFile.readInAllHaplotypes()
			for individualID in beagleFile.sampleIDList:
				sampleID2BeagleFile[individualID] = beagleFile
			# get all haplotypes , etc.
			# get all sample IDs
		self.sampleID2BeagleFile = sampleID2BeagleFile
	
	def reduce(self, **keywords):
		"""
		2012.10.15
			run after all files have been walked through
		"""
		#sample the data
		
		real_counter = 0
		counter = 0
		no_of_loci = 0
		for vcfRecord in self.reader:
			for sampleID, sample_index in vcfRecord.sample_id2index.iteritems():
				beagleFile = self.sampleID2BeagleFile.get(sampleID)
				"""
				if beagleFile is None:
					sys.stderr.write("Warning: sampleID %s is not affiliated with any Beagle file.\n"%(sampleID)
					raise
				"""
				beagleGenotype = beagleFile.getGenotypeOfOneSampleOneLocus(sampleID=sampleID, locusID=None)
				vcfRecord.setGenotypeCallForOneSample(sampleID=sampleID, genotype='%s|%s'%(beagleGenotype[0], beagleGenotype[1]))
				counter += 1
			self.writer.writeVCFRecord(vcfRecord)
			no_of_loci += 1
		sys.stderr.write("%s genotypes, %s loci.\n"%(counter, no_of_loci))
		
		#close the self.invariantPData.writer and self.writer
		AbstractMatrixFileWalker.reduce(self, **keywords)


if __name__ == '__main__':
	main_class = CombinePhasedBeagleOutputsIntoVCF
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
