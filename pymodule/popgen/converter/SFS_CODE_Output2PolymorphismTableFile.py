#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i ./sfs_code.ex1.output -u sfs_code.ex1.fasta -o sfs_code.ex1.polymorphism.h5

Description:
	2013.3.8 this turns SFS_CODE (Hernandez 2008) output into PolymorphismTableFile format.
	Warning: this program is coded with the capability to handle indel/inversion mutations. but un-tested. 
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
from pymodule.yhio.PolymorphismTableFile import PolymorphismTableFile, OneIndividualPolymorphismData


class SFS_CODE_Output2PolymorphismTableFile(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	#option_default_dict.pop(('inputFname', 1, ))
	option_default_dict.update({
						('outputChromosomeSequenceFname', 1, ): ['', '', 1, 'the fasta file to store the chromosome sequences\n\
		(different iteration is regarded as different species, different locus is regarded as different chromosome)'],\
						('ploidy', 0, int):[2, '', 1, 'how many sets of chromosomes one individual carries (in input file=SFS_CODE output, output will be identical).\n\
		ploidy=2 means diploid.'],\
							})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
		self.iterationPattern = re.compile('^//iteration:(\d+)/(\d+)')
		self.chromosomeIndexPattern = re.compile('^locus_(\d+)')
	
	def outputOneIteration(self, inputFile=None, iterationLine=None, outputPolymorphismFile=None, \
						outputChromosomeSequenceFile=None, ploidy=2):
		"""
		2013.3.8
		
		./sfs_code 1 2 -L 2 66
		SEED = -382034166
		//iteration:1/2
		>locus_0
		GTTCCAGGAAGCTGGACAGTCTCTTATGGCGACATGGTAAATAAATTTGCGGTCCTGAAATCGGGT
		>locus_1
		AATGGGTTAGATTATGATTATGTGCATCGTCTTTACACGAGTGGAGTCATTCGACTTTTCGTACTA
		Nc:500;
		MALES:3;
		0,A,24,-237,0,TTA,A,1,Y,N,0.0,2,0.1,0.8;
		//iteration:2/2
		>locus_0
		TTTTCAGTCTGTTTTGTCAAAGATTATTCTTTTGGGCTCCTCACGCACCTTAAGAAGTGTATATAC
		>locus_1
		TACGCTATCAACTACAATATACATAGTGTGGTTTTCGATGGCCTTAGGTCAGTTGACCTACGTAAC
		Nc:500;
		MALES:3;
		1,A,28,-523,0,GTG,A,1,V,E,0.0,2,0.1,0.3;

		"""
		iterationPatternSearchResult = self.iterationPattern.search(iterationLine)
		#add species
		speciesNo = iterationPatternSearchResult.group(1)
		noOfSpecies = iterationPatternSearchResult.group(2)
		speciesName = speciesNo
		speciesEntry = outputPolymorphismFile.addSpecies(name=speciesName, ploidy=ploidy)
		
		# output reference sequence, form chromosome names, add into table
		populationSizeLine = None
		for line in inputFile:
			line = line.strip()	#to get rid of \n or \r
			if line[:3]=='Nc:':	#break as the population size line has been reached
				populationSizeLine = line
				break
			else:
				if line[0]=='>':
					chromosomePureName = self.chromosomeIndexPattern.search(line[1:]).group(1)	#take 0 from locus_0
					chromosomeName = '%s.%s'%(speciesName, chromosomePureName)
					chromosomeSequence = inputFile.next().strip()
					#add chromosome
					outputPolymorphismFile.addChromosome(name=chromosomeName, length=len(chromosomeSequence), \
												speciesName=speciesName, ploidy=ploidy, \
												path=os.path.abspath(outputChromosomeSequenceFile.name))
					# output the chromosome sequence
					record=SeqRecord(Seq(chromosomeSequence, IUPAC.unambiguous_dna), \
									id=chromosomeName, name='', description='')
					SeqIO.write([record], outputChromosomeSequenceFile, "fasta")
		
		# add population + size
		if populationSizeLine is None:
			sys.stderr.write("Error: population size line is None after parsing the locus (chromosome) sequences.\n")
			raise
		#Nc:500,300;
		populationSizeList = map(int, populationSizeLine[3:-1].split(','))
		totalPopulationSize = sum(populationSizeList)
		populationPureName2populationSize = {}
		for populationIndex in xrange(len(populationSizeList)):
			#add population + size
			populationSize = populationSizeList[populationIndex]
			populationName = '%s.%s'%(speciesName, populationIndex)
			populationPureName2populationSize[repr(populationIndex)] = populationSize
			outputPolymorphismFile.addPopulation(name=populationName, size=populationSize, speciesName=speciesName)
		
		#. add individuals (from "MALES:3;" )
		maleStartIndexLine = inputFile.next().strip()
		maleStartIndexList = map(int, maleStartIndexLine[6:-1].split(','))
		for populationIndex in xrange(len(maleStartIndexList)):
			maleStartIndex = maleStartIndexList[populationIndex]
			populationPureName = repr(populationIndex)
			populationName = '%s.%s'%(speciesName, populationPureName)
			for individualIndex in xrange(maleStartIndex):	#females first
				individualName = '%s.%s.%s'%(speciesName, populationPureName, individualIndex)
				outputPolymorphismFile.addIndividual(name=individualName, family_id = None, father_name = None, \
					mother_name = None, sex = 2, phenotype = None, \
					populationName=populationName, speciesName=speciesName, ploidy=ploidy)
			for individualIndex in xrange(maleStartIndex, 2*maleStartIndex):	#males
				individualName = '%s.%s.%s'%(speciesName, populationPureName, individualIndex)
				outputPolymorphismFile.addIndividual(name=individualName, family_id = None, father_name = None, \
					mother_name = None, sex = 1, phenotype = None, \
					populationName=populationName, speciesName=speciesName, ploidy=ploidy)
		
		#. add polymorphism data from last mutation line
		for line in inputFile:
			line = line.strip()
			if self.iterationPattern.search(line):	#a new iteration, a new species
				self.outputOneIteration(inputFile=inputFile, iterationLine=line, outputPolymorphismFile=outputPolymorphismFile,\
									outputChromosomeSequenceFile=outputChromosomeSequenceFile, ploidy=ploidy)
			else:
				"""
				#. add polymorphism data from last mutation line
				"""
				mutationList = line.strip()[:-1].split(';')	#[:-1] gets rid of the last ;, which would result in an empty mutation string
				for mutation in mutationList:
					mutation = mutation.split(',')
					try:
						chromosomePureName, autosomeOrSexChromosome, positionIndex, generationMutationArose, generationMutationFixed,\
							ancestralSequence, derivedSequence, mutation_type = mutation[:8]
					except:
						import pdb
						pdb.set_trace()
					position = int(positionIndex) + 1
					generation_mutation_arose = int(generationMutationArose)
					generation_mutation_fixed = int(generationMutationFixed)
					ancestralSequence = ancestralSequence[1:-1]	#the first and last are the adjacent non-changing nucleotide
					stop = None
					ref_allele_length = None
					alt_allele, alt_allele_length=None,None
					ancestral_amino_acid, derived_amino_acid =None, None 
					if mutation_type =='0' or mutation_type=='1':	#0: synonymous & 1: non-synonymous & 0: noncoding
						ancestral_amino_acid, derived_amino_acid, fitness,\
							noOfChromosomes = mutation[8:12]
						populationHaplotypeIndexList = mutation[12:]
						stop = position
						ref_allele_length = 1
						alt_allele = derivedSequence
						alt_allele_length = 1
					else:
						if mutation_type=='i':
							noOfBasesInvolved, derivedSequence, fitness, noOfChromosomes = mutation[8:12]
							populationHaplotypeIndexList = mutation[12:]
						else:
							noOfBasesInvolved, fitness, noOfChromosomes = mutation[8:11]
							populationHaplotypeIndexList = mutation[11:]
						noOfBasesInvolved = int(noOfBasesInvolved)
						if mutation_type=='i':
							stop = position
							ref_allele_length = 0
							alt_allele = derivedSequence
							alt_allele_length = len(derivedSequence)
						elif mutation_type=='d':
							stop = position + noOfBasesInvolved
							ref_allele_length = noOfBasesInvolved
							alt_allele = None
							alt_allele_length = 0
						elif mutation_type=='v':
							stop = position + noOfBasesInvolved
							ref_allele_length =  noOfBasesInvolved
							alt_allele = derivedSequence
							alt_allele_length = ref_allele_length
					fitness = float(fitness)
					noOfChromosomes = int(noOfChromosomes)
					
					#find the populations in which the mutation is fixed, and then calculate the alt_allele_frequency
						#if mutation is fixed in any population, noOfChromosomes counts only chromosomes in the mutation-segregating populations
					noOfChromosomes = 0
					for i in xrange(len(populationHaplotypeIndexList)):
						populationHaplotypeIndex = populationHaplotypeIndexList[i]
						populationPureName, haplotypeIndex = populationHaplotypeIndex.split('.')
						haplotypeIndex = int(haplotypeIndex)
						if haplotypeIndex ==-1:
							#this mutation is fixed in this population
							noOfChromosomes += populationPureName2populationSize.get(populationPureName)
						else:	#just one chromosome
							noOfChromosomes += 1
						#put it back
						populationHaplotypeIndexList[i] = [populationPureName, haplotypeIndex]
					
					alt_allele_frequency = noOfChromosomes/float(totalPopulationSize)
					
					#add the locus
					chromosomeName = '%s.%s'%(speciesName, chromosomePureName)
					locusName = '%s.%s.%s'%(speciesName, chromosomePureName, position)
					locus = outputPolymorphismFile.getLocus(name=locusName, chromosomeName=chromosomeName,\
						start = position, stop = stop, \
						ref_allele = ancestralSequence, ref_allele_length=ref_allele_length,\
						ref_allele_frequency=1-alt_allele_frequency, \
						alt_allele=alt_allele, alt_allele_length=alt_allele_length, alt_allele_frequency=alt_allele_frequency,\
						generation_mutation_arose=generation_mutation_arose, generation_mutation_fixed=generation_mutation_fixed,\
						mutation_type =mutation_type, fitness = fitness, ancestral_amino_acid =ancestral_amino_acid, \
						derived_amino_acid =derived_amino_acid)
					#add allele for each individual	
					
					
					for populationPureName, haplotypeIndex in populationHaplotypeIndexList:
						if haplotypeIndex==-1:	#no individuals in this population as it's fixed in this population
							continue
						individualIndex = haplotypeIndex/ploidy
						individualName = '%s.%s.%s'%(speciesName, populationPureName, individualIndex)
						chromosome_copy = haplotypeIndex%ploidy
						polymorphismName = 'ind%s.locus%s.c%s'%(individualName, locus.name, chromosome_copy)
						outputPolymorphismFile.getPolymorphism(name=polymorphismName, individualName=individualName, \
									locusName=locus.name, chromosome_copy = chromosome_copy,\
									allele_sequence=alt_allele, allele_sequence_length=alt_allele_length, \
									allele_type =None)
						
						
	
	def run(self):
		"""
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
		outputChromosomeSequenceFile = open(self.outputChromosomeSequenceFname, "w")
		
		commandline = inputFile.next().strip()
		outputPolymorphismFile.addAttribute('commandline', value=commandline, overwrite=True, tableName='polymorphism')
		
		for line in inputFile:
			if self.iterationPattern.search(line):	#one iteration is regarded as one species
				self.outputOneIteration(inputFile=inputFile, iterationLine=line, outputPolymorphismFile=outputPolymorphismFile,\
									outputChromosomeSequenceFile=outputChromosomeSequenceFile, ploidy=self.ploidy)
			
		inputFile.close()
		outputPolymorphismFile.close()
		outputChromosomeSequenceFile.close()

	
if __name__ == '__main__':
	main_class = SFS_CODE_Output2PolymorphismTableFile
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()