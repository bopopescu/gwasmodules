#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i ~/script/lh3_foreign/msHOT-lite/msHOT-lite.output -o /tmp/msHOT-lite.polymorphism.h5

Description:
	2013.3.26 this turns msHOT-lite (Heng Li's custom ms software git@github.com:lh3/foreign.git)
		output into PolymorphismTableFile format.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from pymodule import ProcessOptions, utils
from pymodule import AbstractMapper
from pymodule.yhio.PolymorphismTableFile import PolymorphismTableFile

class msOutput2PolymorphismTableFile(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	#option_default_dict.pop(('inputFname', 1, ))
	option_default_dict.update({
						('outputChromosomeSequenceFname', 0, ): ['', '', 1, 'the fasta file to store the chromosome sequences\n\
	(different iteration is regarded as different species, different locus is regarded as different chromosome)'],\
						('ploidy', 0, int):[2, '', 1, 'how many sets of chromosomes one individual carries (in input file=SFS_CODE output, output will be identical).\n\
	ploidy=2 means diploid.'],\
							})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
	
	def _convert(self, inputFile=None, outputPolymorphismFile=None, ploidy=2):
		"""
		2013.3.26
			
			Heng Li msHOT-lite output looks like:
			./msHOT-lite 2 1 -t 84989.8346003745 -r 34490.1412746802 30000000 -l -en 0.0013 1 0.0670 -en 0.0022 1 0.3866 -en 0.0032 1 0.3446 -en 0.0044 1 0.21
				79 -en 0.0059 1 0.1513 -en 0.0076 1 0.1144 -en 0.0096 1 0.0910 -en 0.0121 1 0.0757 -en 0.0150 1 0.0662 -en 0.0184 1 0.0609 -en 0.0226 1 0.0583 -en
				 0.0275 1 0.0572 -en 0.0333 1 0.0571 -en 0.0402 1 0.0577 -en 0.0485 1 0.0589 -en 0.0583 1 0.0603 -en 0.0700 1 0.0615 -en 0.0839 1 0.0624 -en 0.100
				5 1 0.0632 -en 0.1202 1 0.0641 -en 0.1437 1 0.0651 -en 0.1716 1 0.0663 -en 0.2048 1 0.0678 -en 0.2444 1 0.0696 -en 0.2914 1 0.0719 -en 0.3475 1 0.
				0752 -en 0.4935 1 0.0794 
				//
				@begin 6422
				
				30000000
				1100    01
				6074    10
				
				29966899        10
				29971027        01
				29973740        01
				29982767        01
				29985696        10
				@end
		this does not work on standard-ms output, which looks like:
				msHOT-lite 2 1 -t 4781.50413187402 -r 790.4466018 ...
				//
				segsites: 40567
				
				positions: 0.0002 0.0003
				001001101011011001...
				101001010100101111...
				...
		"""
		sys.stderr.write("Converting ...")
		#add species
		speciesNo = 0
		speciesName = repr(speciesNo)
		speciesEntry = outputPolymorphismFile.addSpecies(name=speciesName, ploidy=ploidy)
		
		
		populationSize = None
		populationIndex = 0
		populationName = '%s.%s'%(speciesName, populationIndex)
		outputPolymorphismFile.addPopulation(name=populationName, size=populationSize, speciesName=speciesName)
		
		isSampleBegun =False	#True when the sample data is up in the next line
		individualsAdded = False
		chromosomePureName = 0
		chromosomeName = None
		previousPolymorhicSitePosition = 0
		no_of_segsites_encountered = 0
		no_of_individuals_added = 0
		for line in inputFile:
			content = line.strip()
			if content=="//":
				segsitesLineContent = inputFile.next().strip()
				no_of_segsites = int(segsitesLineContent.split()[-1])
				chromosomeLength = int(inputFile.next().strip())
				
				chromosomeName = '%s.%s'%(speciesName, chromosomePureName)
				#add chromosome
				outputPolymorphismFile.addChromosome(name=chromosomeName, length=chromosomeLength, \
											speciesName=speciesName, ploidy=ploidy, \
											path=None)
				isSampleBegun =True
			elif content =="@end":
				break
			elif isSampleBegun:
				polymorphicPosition, haplotypeAlleleList = content.split()
				polymorphicPosition = int(polymorphicPosition)
				
				# add the locus
				locusName = '%s.%s.%s'%(speciesName, chromosomePureName, polymorphicPosition)
				locus = outputPolymorphismFile.getLocus(name=locusName, chromosomeName=chromosomeName,\
					start = polymorphicPosition, stop = polymorphicPosition, \
					ref_allele ="0", ref_allele_length=1,\
					ref_allele_frequency=None, \
					alt_allele="1", alt_allele_length=1, alt_allele_frequency=None,\
					generation_mutation_arose=None, generation_mutation_fixed=None,\
					mutation_type ="0", fitness = None, ancestral_amino_acid =None, \
					derived_amino_acid =None)
				# mutation_type 0 means non-coding
				
				# add the individuals
				if not individualsAdded:
					for individualIndex in xrange(0, len(haplotypeAlleleList), self.ploidy):	#go through each indivdiual at a time
						populationName = '%s.%s'%(speciesName, populationIndex)
						individualName = '%s.%s.%s'%(speciesName, populationIndex, individualIndex)
						outputPolymorphismFile.addIndividual(name=individualName, family_id = None, father_name = None, \
							mother_name = None, sex = 0, phenotype = None, \
							populationName=populationName, speciesName=speciesName, ploidy=ploidy)
						no_of_individuals_added += 1
					individualsAdded = True
				
				#add allele for all individuals	
				for individualIndex in xrange(0, len(haplotypeAlleleList), self.ploidy):
					individualGenotype = haplotypeAlleleList[individualIndex*self.ploidy:(individualIndex+1)*self.ploidy]
					
					for chromosome_copy in xrange(len(individualGenotype)):
						individualName = '%s.%s.%s'%(speciesName, populationIndex, individualIndex)
						polymorphismName = 'ind%s.locus%s.c%s'%(individualName, locus.name, chromosome_copy)
						outputPolymorphismFile.getPolymorphism(name=polymorphismName, individualName=individualName, \
									locusName=locus.name, chromosome_copy = chromosome_copy,\
									allele_sequence=individualGenotype[chromosome_copy], allele_sequence_length=1, \
									allele_type =None)
				no_of_segsites_encountered  += 1
				previousPolymorhicSitePosition = polymorphicPosition
		sys.stderr.write(" added %s polymorphic loci, %s individuals.\n"%(no_of_segsites_encountered, no_of_individuals_added))
	
	def run(self):
		"""
		input looks like (inputFileFormat=1)
				msHOT-lite 2 1 -t 4781.50413187402 -r 790.4466018 ...
				//
				segsites: 40567
				
				positions: 0.0002 0.0003
				001001101011011001...
				101001010100101111...
				...
			
			./msHOT-lite 2 1 -t 84989.8346003745 -r 34490.1412746802 30000000 -l -en 0.0013 1 0.0670 -en 0.0022 1 0.3866 -en 0.0032 1 0.3446 -en 0.0044 1 0.21
				79 -en 0.0059 1 0.1513 -en 0.0076 1 0.1144 -en 0.0096 1 0.0910 -en 0.0121 1 0.0757 -en 0.0150 1 0.0662 -en 0.0184 1 0.0609 -en 0.0226 1 0.0583 -en
				 0.0275 1 0.0572 -en 0.0333 1 0.0571 -en 0.0402 1 0.0577 -en 0.0485 1 0.0589 -en 0.0583 1 0.0603 -en 0.0700 1 0.0615 -en 0.0839 1 0.0624 -en 0.100
				5 1 0.0632 -en 0.1202 1 0.0641 -en 0.1437 1 0.0651 -en 0.1716 1 0.0663 -en 0.2048 1 0.0678 -en 0.2444 1 0.0696 -en 0.2914 1 0.0719 -en 0.3475 1 0.
				0752 -en 0.4935 1 0.0794 
				//
				@begin 6422
				
				30000000
				1100    01
				6074    10
				
				29966899        10
				29971027        01
				29973740        01
				29982767        01
				29985696        10
				@end
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		if not os.path.isfile(self.inputFname):
			sys.stderr.write("Error: file, %s,  is not a file.\n"%(self.inputFname))
			sys.exit(3)
			
		inputFile = utils.openGzipFile(self.inputFname, 'r')
		outputPolymorphismFile = PolymorphismTableFile(self.outputFname, openMode='w', isPhased=1, \
														ploidy=self.ploidy)
		
		commandline = inputFile.next().strip()
		outputPolymorphismFile.addAttribute('commandline', value=commandline, overwrite=True, tableName='polymorphism')
		
		self._convert(inputFile=inputFile, outputPolymorphismFile=outputPolymorphismFile, ploidy=self.ploidy)
		
		inputFile.close()
		outputPolymorphismFile.close()

	
if __name__ == '__main__':
	main_class = msOutput2PolymorphismTableFile
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()