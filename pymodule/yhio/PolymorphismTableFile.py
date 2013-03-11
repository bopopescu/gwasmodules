#!/usr/bin/env python
"""
2013.3.6 data structure for polymorphism data in YHPyTables (custom of PyTables) format
	testing examples could be found in ~/script/pymodule/popgen/examples
"""
import sys, os, math
__doc__ = __doc__%()

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import csv
import tables
import numpy
import networkx as nx
from tables import UInt64Col, Float64Col, StringCol, UInt16Col, Int64Col
from pymodule.utils import PassingData, PassingDataList
from pymodule.ProcessOptions import  ProcessOptions
from pymodule.yhio.YHPyTables import YHTable, YHFile, castPyTablesEntryIntoPassingData
from pymodule.yhio.SNP import SNPData

class SpeciesTable(tables.IsDescription):
	"""
	2013.3.8
	"""
	id = UInt64Col(pos=0)
	name = StringCol(512, pos=1)
	scientific_name = StringCol(512, pos=2)
	ploidy = UInt16Col(pos=3)

class PopulationTable(tables.IsDescription):
	"""
	2013.3.8
	"""
	id = UInt64Col(pos=0)
	name = StringCol(512, pos=1)	#name should like 'species_name.population_name'
	size = UInt64Col(pos=2)
	species_id = UInt64Col(pos=3)

class ChromosomeTable(tables.IsDescription):
	"""
	2013.3.8
	"""
	id = UInt64Col(pos=0)
	name = StringCol(512, pos=1)	#should be unique, but not 100% enforced
		#name should look like species_name.chromosome_name, if species_name is available
	length = UInt64Col(pos=2)
	sex_chromosome = StringCol(4, pos=3)	#A=autosome, X=chr X, Y=chr Y
	species_id = UInt64Col(pos=4)
	path = StringCol(4512, pos=5)	#path to the file containing the chromosome sequences, should be unique, 

class LocusTable(tables.IsDescription):
	id = UInt64Col(pos=0)
	name = StringCol(512, pos=1)	#'species_name.chromosome_name.locus_name'
	chromosome_id = UInt64Col(pos=2)	#64 byte-long
	start = UInt64Col(pos=3)
	stop = UInt64Col(pos=4)
	ref_allele = StringCol(64, pos=5)
	ref_allele_length = UInt64Col(pos=6)	#this could be larger than 64, good when ref_allele is simply too long
	ref_allele_frequency = Float64Col(pos=7)
	
	alt_allele = StringCol(64, pos=8)
	alt_allele_length = UInt64Col(pos=9)	#this could be larger than 64, good when ref_allele is simply too long
	alt_allele_frequency = Float64Col(pos=10)
	
	generation_mutation_arose = Int64Col(pos=11)
	generation_mutation_fixed = Int64Col(pos=12)
	mutation_type = StringCol(4, pos=13)	#0=synonymous/non-coding, 1=non-synonymous, i=insertion, d=deletion, v=inversion
	fitness = Float64Col(pos=14)
	ancestral_amino_acid = StringCol(4, pos=15)	#only for mutation_type 0/1
	derived_amino_acid = StringCol(4, pos=16)
	

class IndividualTable(tables.IsDescription):
	id = UInt64Col(pos=0)
	family_id = StringCol(512, pos=1)	#64 byte-long
	name = StringCol(512, pos=2)	# name should look like 'species_name.population_name.individual_name' to ensure uniqueness
	father_name = StringCol(512, pos=3)
	mother_name = StringCol(512, pos=4)
	sex = UInt64Col(pos=5)	#0 is unknown, 1=male, 2=female
	phenotype = Float64Col(pos=6)
	species_id = UInt64Col(pos=7)
	population_id = UInt64Col(pos=8)
	

class PolymorphismTable(tables.IsDescription):
	"""
	2013.3.6 table that records the polymorphism alleles of individuals at relevant loci
	"""
	id = UInt64Col(pos=0)
	name = StringCol(512, pos=1)	#name should be "individualName.locusName.chromosome_copy" to ensure uniqueness
	individual_id = UInt64Col(pos=2)
	locus_id = UInt64Col(pos=3)
	chromosome_copy = UInt64Col(pos=4, dflt=0)
	allele_sequence = StringCol(512, pos=5)
	allele_sequence_length = UInt64Col(pos=6)
	allele_type = StringCol(512, pos=7)

class OneIndividualPolymorphismData(PassingData):
	"""
	2013.3.8 to store an individual's polymorphism data, i.e. haplotypes or genotypes
		used in SimulatePedigreeHaplotype.py
	"""
	def __init__(self, **keywords):
		self.isPhased = None
		self.ploidy = None
		self.locusIDList = []
		self.haplotypeList = []
		self.locusPositionList = []
		PassingData.__init__(self, **keywords)
	
	def addHaplotype(self, haplotype=None):
		self.haplotypeList.append(haplotype)


class RecombinationTable(tables.IsDescription):
	"""
	2013.3.7 in case a pedigree is stored in IndividualTable, this table records the recombination events in meiosis.
		Used in SimulatePedigreeHaplotype.py
	"""
	id = UInt64Col(pos=0)
	parent_id = UInt64Col(pos=1)
	child_id = UInt64Col(pos=2)
	position = UInt64Col(pos=3, dflt=0)

class PolymorphismTableFile(YHFile):
	"""
	2013.3.6 usage examples:
		
		#for writing
		landscapeFile = AssociationLandscapePyTable(self.outputFname, openMode='w', isPhased=1, ploidy=2)
		landscapeFile.addAttributeDict(attributeDict)
		
		#for read-only,
		landscapeTable = AssociationLandscapePyTable(self.inputFname, openMode='r')
	"""
	def __init__(self, inputFname=None, openMode='r', \
				tableName='polymorphism', groupNamePrefix='group', tableNamePrefix='table',\
				filters=None, autoRead=True, autoWrite=True, \
				isPhased=None, ploidy=None, constructSNPData=True, **keywords):
		
		
		self.bridge_ls = None
		self.locusLandscapeNeighborGraph = None
		
		YHFile.__init__(self, inputFname=inputFname, openMode=openMode, \
				tableName=tableName, groupNamePrefix=groupNamePrefix, tableNamePrefix=tableNamePrefix,\
				rowDefinition=None, filters=filters, \
				debug=0, report=0, autoRead=False, autoWrite=False)
		
		self.speciesTableName = 'species'
		self.populationTableName = 'population'
		self.individualTableName = "individual"
		self.chromosomeTableName = 'chromosome'
		self.locusTableName = 'locus'
		self.recombinationTableName = 'recombination'
		
		self.isPhased = isPhased
		self.ploidy = ploidy
		self.constructSNPData = constructSNPData
		
		#to overwrite self.autoRead that is set by YHFile.__init__
		self.autoRead = autoRead
		self.autoWrite = autoWrite
		
		self.snpData = None	#the SNPData structure that holds all polymorphism, locus, individual info
		
		if self.autoRead and (self.openMode=='r' or self.openMode=='a'):
			self.speciesTable = self.getTableObject(tableName=self.speciesTableName)
			self.populationTable = self.getTableObject(tableName=self.populationTableName)
			self.individualTable = self.getTableObject(tableName=self.individualTableName)
			self.chromosomeTable = self.getTableObject(tableName=self.chromosomeTableName)
			self.locusTable = self.getTableObject(tableName=self.locusTableName)
			self.recombinationTable = self.getTableObject(tableName=self.recombinationTableName)
			self.polymorphismTable = self.getTableObject(tableName=self.tableName)
			
			#read the isPhased, ploidy from pytables attributes, overwrites the arguments
			self.isPhased = self.polymorphismTable.getAttribute(name='isPhased', defaultValue=0)
			self.ploidy = self.polymorphismTable.getAttribute(name='ploidy', defaultValue=2)
			
			self._readInData(tableName=self.tableName, tableObject=self.associationLandscapeTable)
		if self.autoWrite and self.openMode=='w':
			self.speciesTable = self.createNewTable(tableName=self.speciesTableName, rowDefinition=SpeciesTable,\
													expectedrows=500)
			self.populationTable = self.createNewTable(tableName=self.populationTableName, rowDefinition=PopulationTable,\
													expectedrows=500)
			self.individualTable = self.createNewTable(tableName=self.individualTableName, rowDefinition=IndividualTable,\
													expectedrows=30000)
			self.chromosomeTable = self.createNewTable(tableName=self.chromosomeTableName, rowDefinition=ChromosomeTable,\
													expectedrows=500)
			self.locusTable = self.createNewTable(tableName=self.locusTableName, rowDefinition=LocusTable,\
												expectedrows=300000)
			self.recombinationTable = self.createNewTable(tableName=self.recombinationTableName, rowDefinition=RecombinationTable,\
												expectedrows=300000)
			self.polymorphismTable = self.createNewTable(tableName=self.tableName, rowDefinition=PolymorphismTable,\
												expectedrows=500000)
			#set the attributes of isPhased, ploidy
			self.polymorphismTable.addAttribute(name='isPhased', value=self.isPhased, overwrite=True)
			self.polymorphismTable.addAttribute(name='ploidy', value=self.ploidy, overwrite=True)
		
		#2013.3.8 these dictionaries are for outputting purposes 
		self._individualName2ID = {}
		self._locus_index2id = {}
		
		#2013.3.8 helper structures
		self._locusStartPositionList = None
		self._locusChrStartStopList = None
	
	def _readInData(self, tableName=None, tableObject=None, bugfixType=None):
		"""
		2013.3.6
		"""
		if tableName is None:
			tableName = self.tableName
		YHFile._readInData(self, tableName=tableName, tableObject=tableObject)
		if not self.constructSNPData:
			return
		
		sys.stderr.write("Reading everything into a SNPData structure ...")
		row_id_list = []
		row_id_number2row_index = {}
		col_id_list = []
		col_id_number2col_index = {}
		for row in self.individualTable:
			row_id_list.append(row.name)
			row_id_number2row_index[row.id] = len(row_id_list)-1
		for row in self.locusTable:
			#col_id_list.append(row.id)
			col_id_list.append((row.chromosome_id, row.start, row.stop))
			col_id_number2col_index[row.id] = len(col_id_list)-1
		
		allele_sequence2allele_number = {}
		allele_number2allele_sequence = {}
		
		#each cell in data_matrix is an array of alleles for one individual at one locus, but different chromosomes
		# alleles are encoded in numbers starting from 1. 0 is missing.
		data_matrix = numpy.zeros([len(row_id_list), len(col_id_list), self.ploidy], dtype=numpy.int16)
		
		if self.ploidy>1:
			#chromosome_copy_matrix is used to keep track of the chromosomes for particular individual & locus
			chromosome_copy_matrix = numpy.zeros([len(row_id_list), len(col_id_list)], dtype=numpy.int8)
		else:
			chromosome_copy_matrix = None
		
		for row in self.polymorphismTable:
			row_index = row_id_number2row_index.get(row.individual_id)
			col_index = col_id_number2col_index.get(row.locus_id)
			
			#figure out which chromosome to hold this allele
			if self.ploidy>1:
				chromosome_copy_matrix[row_index][col_index] = chromosome_copy_matrix[row_index][col_index]+1
				if row.chromosome_copy == 0:	#unphased genotype
					chromosome_copy_index = chromosome_copy_matrix[row_index][col_index] -1
				else:
					chromosome_copy_index = row.chromosome_copy-1
			else:
				chromosome_copy_index = 0
				if row.chromosome_copy>1:
					sys.stderr.write("Warning: ploidy=%s, but encounter chromosome_copy (%s) >1.\n"%\
									(self.ploidy, row.chromosome_copy))
			
			#allele_number starts from 1. 0 is reserved for missing.
			if row.allele_sequence not in allele_sequence2allele_number:
				allele_sequence2allele_number[row.allele_sequence] = len(allele_sequence2allele_number)+1
				allele_number = allele_sequence2allele_number.get(row.allele_sequence)
				allele_number2allele_sequence[allele_number] = row.allele_sequence
			
			allele_number = allele_sequence2allele_number.get(row.allele_sequence)
			data_matrix[row_index][col_index][chromosome_copy_index] = allele_number
		self.snpData = SNPData(row_id_list=row_id_list, col_id_list=col_id_list, data_matrix=data_matrix)
		
		self.snpData.allele_sequence2allele_number = allele_sequence2allele_number
		self.snpData.allele_number2allele_sequence = allele_number2allele_sequence
		sys.stderr.write(" %s individuals, %s loci, ploidy=%s, isPhased=%s.\n"%(len(self.snpData.row_id_ls),\
																	len(self.snpData.col_id_ls), \
																	self.ploidy, self.isPhased))
		
		return self.snpData
	
	def getOneHaplotype(self, individual_id=None, chromosome_copy=None):
		"""
		2013.3.7
			if data is unphased and multi-ploidy, the allele combination is entirely random.
			if data is phased, and multi-ploidy, it will randomly select one chromosome copy in advance.
		"""
		if chromosome_copy is None and self.isPhased:	#randomly select one copy if data is phased
			if self.ploidy>=2:	#
				chromosome_copy  = numpy.random.randint(1, self.ploidy+1, size=1)	#sample one integer from [1, self.ploidy]
				
			else:
				chromosome_copy = None	#do not add it to query restriction condition
		where_condition_ls = ["individual_id==%s"%(individual_id)]
		if chromosome_copy is not None:
			where_condition_ls.append(("chromosome_copy==%s"%(chromosome_copy)))
		query = self.polymorphismTable.query(""" %s """%(" & ".join(where_condition_ls)))
		haplotype = ['']*self.locusTable.nrows
		for row in query:
			locus_index =row['locus_id']-1
			haplotype[locus_index] = row['allele_sequence']
		return haplotype

	def getOneIndividualGenotypeList(self, individual_id=None):
		"""
		2013.3.7
			to get unphased genotype data
		"""
		where_condition_ls = ["individual_id==%s"%(individual_id)]
		query = self.polymorphismTable.query(""" %s """%(" & ".join(where_condition_ls)))
		genotypeList = ['']*self.locusTable.nrows
		for row in query:
			locus_index =row['locus_id']-1
			genotypeList[locus_index] = genotypeList[locus_index] + row['allele_sequence']
		return genotypeList
		
	def sampleOneRandomHaplotypeWithReplacement(self):
		"""
		2013.3.7 randomly select a haplotype and return
		"""
		no_of_individuals = self.individualTable.no_of_rows
		random_individual_id = numpy.random.randint(1, no_of_individuals+1, size=1)
		return self.getOneHaplotype(individual_id=random_individual_id, chromosome_copy=None)
	
	def sampleOneIndividualPolymorphismWithReplacement(self):
		"""
		2013.3.7 randomly select a haplotype and return
		"""
		no_of_individuals = self.individualTable.no_of_rows
		random_individual_id = numpy.random.randint(1, no_of_individuals+1, size=1)
		polymorphismData = OneIndividualPolymorphismData(isPhased=self.isPhased, ploid=self.ploidy)
		if self.isPhased and self.ploidy>0:
			for i in xrange(1,self.ploidy+1):
				haplotype = self.getOneHaplotype(individual_id=random_individual_id, chromosome_copy=i)
				polymorphismData.addHaplotype(haplotype)
		else:
			genotypeList = self.getOneIndividualGenotypeList(individual_id=random_individual_id)
			polymorphismData.addHaplotype(genotypeList)
		return polymorphismData
	
	@property
	def locusStartPositionList(self):
		"""
		2013.3.7 return a list of locus.start
		"""
		if self._locusStartPositionList is None:
			self._locusStartPositionList = []
			for row in self.locusTable:
				self._locusStartPositionList.append(row['start'])
			self._locusStartPositionList.sort()
		return self._locusStartPositionList
	
	@property
	def locusChrStartStopList(self):
		"""
		2013.3.8 return a list of (locus.chromosome, locus.start, locus.stop)
		"""
		if self._locusChrStartStopList is None:
			self._locusChrStartStopList = []
			for row in self.locusTable:
				chromosomeEntry = self.getChromosome(id=row['chromosome_id'])
				self._locusChrStartStopList.append((chromosomeEntry.name, row['start'], row['stop']))
			self._locusChrStartStopList.sort()
		return self._locusChrStartStopList
	
	def checkIfEntryInTable(self, tableObject=None, name=None, id=None):
		"""
		2013.3.8
			for Population, Species, Chromosome, similar structure
		"""
		returnValue = None
		if tableObject and (name or id):
			where_condition_ls = []
			if name:
				where_condition_ls.append("name=='%s'"%(name))
			if id:
				where_condition_ls.append("id==%s"%(id))
			rows = tableObject.readWhere("""%s"""%(' & '.join(where_condition_ls)))
			if len(rows)==1:
				returnValue = rows[0]
			elif len(rows)>1:
				sys.stderr.write("Error: %s (>1) %s entries with same name.\n"%(len(rows), tableObject.name))
				raise
			else:
				returnValue = False
		if returnValue:
			returnValue = castPyTablesEntryIntoPassingData(returnValue)
		return returnValue
	
	def checkChromosome(self, name=None, id=None):
		"""
		2013.3.8
		"""
		return self.checkIfEntryInTable(tableObject=self.chromosomeTable, name=name, id=id)
	
	def addChromosome(self, name=None, length=None, speciesName=None, ploidy=None, path=None):
		"""
		2013.3.8
		"""
		if name:
			species_id = None
			if speciesName:
				speciesEntry = self.getSpecies(name=speciesName, ploidy=ploidy)
			if speciesEntry:
				species_id = speciesEntry.id
			oneCell = PassingData(name=name, length=length, species_id=species_id, path=path)
			self.chromosomeTable.writeOneCell(oneCell, cellType=2)
			self.flush()
		return self.checkChromosome(name=name)	#would this work without flush()?
	
	def getChromosome(self, name=None, id=None, length=None, speciesName=None, path=None, \
					ploidy=None, getRandomChromosome=False):
		"""
		2013.3.8
			#. if name is None, try get 1st entry in chromosome table
			#. ploidy is for species
		"""
		entry = None
		if name or id:
			entry= self.checkChromosome(name=name, id=id)
			if not entry:
				entry = self.addChromosome(name=name, length=length, speciesName=speciesName, ploidy=ploidy, path=path)
		elif getRandomChromosome:
			entry = self.getFirstChromosome()
		return entry
	
	def getFirstChromosome(self):
		"""
		2013.3.8
		"""
		entry = None
		rows = self.chromosomeTable.read(start=0, stop=1)
		if len(rows)==1:
			entry= rows[0]
		return entry
	
	def checkPopulation(self, name=None, id=None):
		"""
		2013.3.8
		"""
		return self.checkIfEntryInTable(tableObject=self.populationTable, name=name, id=id)
	
	def addPopulation(self, name=None, size=None, speciesName=None):
		"""
		2013.3.8
		"""
		if name:
			species_id = None
			if speciesName:
				species = self.getSpecies(name=speciesName)
			if species:
				species_id = species.id
			oneCell = PassingData(name=name, size=size, species_id=species_id)
			self.populationTable.writeOneCell(oneCell, cellType=2)
			self.flush()
		return self.checkPopulation(name=name)	#would this work without flush()?
	
	def getPopulation(self, name=None, id=None, size=None, speciesName=None):
		"""
		2013.3.8
		"""
		entry = None
		if name or id:
			entry= self.checkPopulation(name=name, id=id)
			if not entry:
				entry = self.addPopulation(name=name, size=size, speciesName=speciesName)
		return entry
	
	def checkSpecies(self, name=None, id=None):
		"""
		2013.3.8
		"""
		return self.checkIfEntryInTable(tableObject=self.speciesTable, name=name, id=id)
	
	def addSpecies(self, name=None, scientific_name=None, ploidy=None):
		"""
		2013.3.8
		"""
		if name:
			oneCell = PassingData(name=name, scientific_name=scientific_name, ploidy=ploidy)
			self.speciesTable.writeOneCell(oneCell, cellType=2)
			self.flush()
		return self.checkSpecies(name=name)	#would this work without flush()?
	
	def getSpecies(self, name=None, id=None, scientific_name=None, ploidy=None):
		"""
		2013.3.8
		"""
		entry = None
		if name or id:
			entry= self.checkSpecies(name=name, id=id)
			if not entry:
				entry = self.addSpecies(name=name, scientific_name=scientific_name, ploidy=ploidy)
		return entry
	
	def checkIndividual(self, name=None, id=None):
		"""
		2013.3.8
		"""
		return self.checkIfEntryInTable(tableObject=self.individualTable, name=name, id=id)
	
	def addIndividual(self, name=None, family_id = None, father_name = None, \
					mother_name = None, sex = None, phenotype = None, \
					populationName=None, speciesName=None, ploidy=None):
		"""
		2013.3.8
		"""
		if name:
			population_id = None
			species_id = None
			if speciesName:
				species = self.getSpecies(name=speciesName, ploidy=ploidy)
				if species:
					species_id = species.id
			if populationName:
				population = self.getPopulation(name=populationName, speciesName=speciesName)
				if population:
					population_id = population.id
			oneCell = PassingData(name=name,family_id =family_id, father_name=father_name,\
								mother_name=mother_name, sex=sex, phenotype=phenotype,\
								species_id=species_id, population_id=population_id)
			self.individualTable.writeOneCell(oneCell, cellType=2)
			self.flush()
			if name in self._individualName2ID:
				sys.stderr.write("Error: individual %s is not unique, already in _individualName2ID with id=%s.\n"%\
								(name, self._individualName2ID.get(name)))
				raise
			else:
				self._individualName2ID[name] = self.individualTable.no_of_rows	#nrows is not updated until flush()
		
		return self.checkIndividual(name=name)	#would this work without flush()?
	
	def getIndividual(self, name=None, id=None, family_id = None, father_name = None, \
					mother_name = None, sex = None, phenotype = None, \
					populationName=None, speciesName=None, ploidy=None):
		"""
		2013.3.8 query the individual table
		"""
		entry = None
		if name or id:
			entry= self.checkIndividual(name=name, id=id)
			if not entry:
				entry = self.addIndividual(name=name, family_id=family_id, father_name=father_name, \
								mother_name=mother_name, sex=sex, phenotype=phenotype, \
								populationName=populationName, speciesName=speciesName, ploidy=ploidy)
		return entry
	
	def checkLocus(self, name=None, id=None):
		"""
		2013.3.10
		"""
		return self.checkIfEntryInTable(tableObject=self.locusTable, name=name, id=id)
	
	def addLocus(self, name=None, chromosomeName=None,\
				start = None, stop = None, ref_allele = None, ref_allele_length=None,\
				ref_allele_frequency =None, alt_allele=None, alt_allele_length=None,\
				alt_allele_frequency=None, generation_mutation_arose=None, generation_mutation_fixed=None,\
				mutation_type =None, fitness = None, ancestral_amino_acid =None, \
				derived_amino_acid =None, **keywords):
		"""
		2013.3.8
		"""
		if name:
			chromosome_id = None
			if chromosomeName:
				chromosome_id = self.getChromosome(name=chromosomeName).id
			oneCell = PassingData(name=name, chromosome_id=chromosome_id, start = start, stop = stop, \
					ref_allele = ref_allele, ref_allele_length=ref_allele_length,\
					ref_allele_frequency =ref_allele_frequency, alt_allele=alt_allele, alt_allele_length=alt_allele_length,\
					alt_allele_frequency=alt_allele_frequency, generation_mutation_arose=generation_mutation_arose, \
					generation_mutation_fixed=generation_mutation_fixed,\
					mutation_type =mutation_type, fitness = fitness, ancestral_amino_acid =ancestral_amino_acid, \
					derived_amino_acid =derived_amino_acid, **keywords)
			self.locusTable.writeOneCell(oneCell, cellType=2)
			self.flush()
		return self.checkLocus(name=name)	#would this work without flush()?
	
	def getLocus(self, name=None, id=None, chromosomeName=None,\
				start = None, stop = None, ref_allele = None, ref_allele_length=None,\
				ref_allele_frequency =None, alt_allele=None, alt_allele_length=None,\
				alt_allele_frequency=None, generation_mutation_arose=None, generation_mutation_fixed=None,\
				mutation_type =None, fitness = None, ancestral_amino_acid =None, \
				derived_amino_acid =None, **keywords):
		"""
		2013.3.8
		"""
		entry = None
		if name or id:
			entry= self.checkLocus(name=name, id=id)
			if not entry:
				entry = self.addLocus(name=name, chromosomeName=chromosomeName,\
					start = start, stop = stop, ref_allele = ref_allele, ref_allele_length=ref_allele_length,\
					ref_allele_frequency =ref_allele_frequency, alt_allele=alt_allele, alt_allele_length=alt_allele_length,\
					alt_allele_frequency=alt_allele_frequency, generation_mutation_arose=generation_mutation_arose, \
					generation_mutation_fixed=generation_mutation_fixed,\
					mutation_type =mutation_type, fitness = fitness, ancestral_amino_acid =ancestral_amino_acid, \
					derived_amino_acid =derived_amino_acid, **keywords)
		return entry
	
	def checkPolymorphism(self, name=None, id=None):
		"""
		2013.3.10
		"""
		return self.checkIfEntryInTable(tableObject=self.polymorphismTable, name=name, id=id)
	
	def addPolymorphism(self, name=None, individualName=None, locusName=None, chromosome_copy = None,\
					allele_sequence=None, allele_sequence_length=None, allele_type =None, **keywords):
		"""
		2013.3.10
		"""
		if name:
			individual_id = self.getIndividual(name=individualName).id
			locus_id = self.getLocus(name=locusName).id
			oneCell = PassingData(name=name, individual_id=individual_id, locus_id = locus_id, \
					chromosome_copy=chromosome_copy,\
					allele_sequence = allele_sequence, allele_sequence_length=allele_sequence_length,\
					allele_type=allele_type, **keywords)
			self.polymorphismTable.writeOneCell(oneCell, cellType=2)
			self.flush()
		return self.checkPolymorphism(name=name)	#would this work without flush()?
	
	def getPolymorphism(self, name=None, id=None, individualName=None, locusName=None, chromosome_copy = None,\
					allele_sequence=None, allele_sequence_length=None, allele_type =None, **keywords):
		"""
		2013.3.10
		"""
		entry = None
		if name or id:
			entry= self.checkPolymorphism(name=name, id=id)
			if not entry:
				entry = self.addPolymorphism(name=name, individualName=individualName, locusName =locusName, \
					chromosome_copy=chromosome_copy,\
					allele_sequence = allele_sequence, allele_sequence_length=allele_sequence_length,\
					allele_type=allele_type, **keywords)
		return entry
	
	def writePedigreeDiGraph2IndividualTable(self, diGraph=None, populationName=None, speciesName=None, \
											ploidy=None, **keywords):
		"""
		2013.3.7 diGraph is an instance of pymodule.algorithm.graph.DiGraphWrapper
			#. make sure all parents have an entry in the table (so that they could be referenced in recombination-table)
			#. establish _individualName2ID map so that it could be used in outputting recombination or polymorphism data
			#. sex is set to 0 (unknown). 1=male, 2=female
			#. no good for appending to an individual table with entries in it already
		"""
		sys.stderr.write("Writing the pedigree directed graph out ...")
		
		family_id = "1"
		for individualName in diGraph:
			parents = diGraph.predecessors(individualName)
			father_name = None
			mother_name = None
			if len(parents)>=1:
				father_name = parents[0]
			if len(parents)>=2:
				mother_name = parents[1]
			self.addIndividual(name=individualName, family_id=family_id, father_name=father_name, \
							mother_name=mother_name, sex=0, phenotype=None, populationName=populationName, \
							speciesName=speciesName, ploidy=ploidy)
			
		sys.stderr.write(" %s individuals.\n"%(self.individualTable.no_of_rows))
	
	def writeChrStartStopTupleList2LocusTable(self, chr_start_stop_list=None, chromosomeLength=None,\
											speciesName=None, ploidy=None):
		"""
		2013.3.7
			#. establish _locus_index2id, to be used in writeIndividualName2PolymorphismData()
			#. make sure chr_start_stop_list is in the same order as the haplotype in writeIndividualName2PolymorphismData()
		"""
		sys.stderr.write("Writing a %s-element list of (chr, start,stop) out ..."%(len(chr_start_stop_list)))
		chr_start_stop_list.sort()	#make sure it's sorted
		if ploidy is None:
			ploidy=self.ploidy
		for i in xrange(len(chr_start_stop_list)):
			chromosomeName, start, stop = chr_start_stop_list[i][:3]
			if chromosomeName:
				chromosomeEntry = self.getChromosome(name=chromosomeName, length=chromosomeLength, speciesName=speciesName,\
													ploidy=ploidy)
			else:
				chromosomeEntry = None
			name = '%s_%s_%s'%(chromosomeName, start, stop)
			oneCell = PassingData(name=name, chromosome_id=getattr(chromosomeEntry, 'id', None), start=start, stop=stop)
			self.locusTable.writeOneCell(oneCell, cellType=2)
			self._locus_index2id[i] = self.locusTable.no_of_rows
		sys.stderr.write("%s loci \n")
		return self._locus_index2id
	
	def writeIndividualName2PolymorphismData(self, individualName2polymorphismData=None, \
											locus_index2id=None, speciesName=None, ploidy=None):
		"""
		2013.3.7
			if locus_index2id is not available, raise error
		"""
		sys.stderr.write("Writing individualName2polymorphismData (%s individuals) out ..."%\
						(len(individualName2polymorphismData)))
		if locus_index2id is None:
			locus_index2id = self._locus_index2id
		counter = 0
		for individualName, polymorphismData in individualName2polymorphismData.iteritems():
			individual_id = self.getIndividual(individualName, speciesName=speciesName, ploidy=ploidy).id
			for i in xrange(len(polymorphismData.haplotypeList)):
				haplotype = polymorphismData.haplotypeList[i]
				for j in xrange(len(haplotype)):
					locus_id = locus_index2id.get(j)
					if locus_id is None:
						sys.stderr.write("Error: no locus_id for locus index %s.\n"%(j))
						raise
					oneCell = PassingData(individual_id=individual_id, locus_id=locus_id,\
										chromosome_copy=i, allele_sequence=haplotype[j],\
										allele_sequence_length=len(haplotype[j]), allele_type=1)
					self.polymorphismTable.writeOneCell(oneCell, cellType=2)
					counter += 1
		sys.stderr.write(" %s alleles outputted.\n"%(counter))
		
	
	def writeRecombinationEvents(self, parentName=None, childName=None, recombinationLocationList=None):
		"""
		2013.3.7
		"""
		parent_id = self.getIndividual(parentName).id
		child_id = self.getIndividual(childName).id
		for position in recombinationLocationList:
			oneCell = PassingData(parent_id=parent_id, child_id=child_id, position=position)
			self.recombinationTable.writeOneCell(oneCell, cellType=2)
		
		
	