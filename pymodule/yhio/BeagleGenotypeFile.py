#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s

Description:
	2013.05.03 a child class of MatrixFile. used to describe Beagle likelihood file format, which looks like:
		marker alleleA alleleB 1000_709_1996093_GA_vs_524 1000_709_1996093_GA_vs_524 1000_709_1996093_GA_vs_524 1001_710_1995025_GA_vs_524 1001_710_1995025_GA_vs_524 1001_710_1995025_GA_vs_524 1002_711_2001039_GA_vs_524
		Contig791:1086 C A 1 0 0 0.9997 0.0003 0 0
		Contig791:1649 G C 1 0 0 1 0 0 0
		Contig791:4084 A C 1 0 0 1 0 0 0
		Contig791:4118 A G 1 0 0 1 0 0 0
		Contig791:4143 C A 1 0 0 1 0 0 0
		Contig791:4168 G C 1 0 0 0.9999 0.0001 0 0
		Contig791:4203 C G 1 0 0 1 0 0 0

	Example:
		beagleFile = BeagleGenotypeFile(inputFname='/tmp/input.bgl.gz')
		beagleFile.readInAllHaplotypes()
		
		for individualID, firstHaplotypeIndex in beagleFile.snpData.col_id2col_index.iteritems():
			haplotypeList = []
			for j in xrange(firstHaplotypeIndex, firstHaplotypeIndex+self.ploidy):
				haplotypeList.append(beagleFile.snpData.data_matrix[:,j]) 
			# another way
			#haplotypeList = beagleFile.getHaplotypeListOfOneSample(individualID)
		
		reader = MatrixFile(inputFname='/tmp/input.txt', openMode='r')
		reader = MatrixFile('/tmp/input.txt', openMode='r')
		reader.constructColName2IndexFromHeader()
		for row in reader:
			row[reader.getColName2IndexFromHeader('KID')]
		
		inf = utils.openGzipFile(inputFname, openMode='r')
		reader = MatrixFile(inputFile=inf)
		
		#2013.2.1 writing
		writer = MatrixFile('/tmp/output.txt', openMode='w', delimiter='\t')
		writer.writeHeader(...)
		writer.writerow(row)
		writer.close()

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy, numpy
from pymodule import utils, PassingData
from pymodule.ProcessOptions import ProcessOptions
from pymodule.yhio.SNP import SNPData
from MatrixFile import MatrixFile

class BeagleGenotypeFile(MatrixFile):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(MatrixFile.option_default_dict)
	option_default_dict.update({
						('delimiter', 0, ): [' ', '', 1, 'delimiter for Beagle likelihood format is single-space'],\
						('phased', 0, int): [0, '', 1, 'whether the genotypes are phased or not'],\
						('ploidy', 0, int): [2, '', 1, 'ploidy of the organism. = #genotype-columns for one individual'],\
						})
	def __init__(self, inputFname=None, **keywords):
		
		MatrixFile.__init__(self, inputFname=inputFname, **keywords)
		
		self.header = None
		self.col_name2index = None	#key is sampleID, value is index of first haplotype
		
		self.sampleIDList = []
		self.sampleID2index = {}	#same as col_name2index
		self.locusIDList = []
		self.haplotypeMatrix = []
		
		self.snpData = None	#to store everything above . SNPData type
	
	def _completeSNPData(self):
		"""
		2013.06.05
		"""
		if self.snpData is None:
			if self.openMode=='r':
				self.readInAllHaplotypes()
			else:	#writing mode
				self._consolidateHaplotypeMatrixIntoSNPData()
	
	def getSampleIndexGivenSampleID(self, sampleID=None):
		"""
		2013.06.05
		"""
		return self.getColIndexGivenColHeader(colHeader=sampleID)
	
	def constructColName2IndexFromHeader(self):
		"""
		2013.05.03
			First two columns are for designators.
			one diploid sample occupies two columns.
			
		I id 1466_1867_2001024_GA_vs_524copy3 1466_1867_2001024_GA_vs_524copy3 1243_1338_2005057_GA_vs_524copy3 1511_639_1987079_GA_vs_524copy3 1511_639_1987079_GA_vs_524copy3 1016_725_1995116_GA_vs_524copy7 1521_767_1993041_GA_vs_524copy4 1521_767_1993041_GA_vs_524copy4
		M Contig791:1086 C C C A A C C C
		M Contig791:1649 G G G C C G G G
		
		"""
		self.col_name2index = {}	#col_name2index does not include the first two columns, starting only from samples
			#same as sampleID2index
		headerLineData = self.next()
		if headerLineData.lineIndicator=='I':
			self.header = headerLineData.genotypeList
			for i in xrange(2,len(self.header), self.ploidy):
				self.sampleIDList.append(self.header[i])
		else:
			if headerLineData.lineIndicator=='M':	#there is no header line. roll back to beginning
				self.csvFile.seek(0)
			#make up a header
			self.header = []
			for i in xrange(0, len(headerLineData.genotypeList),self.ploidy):
				sampleID = 'sampleColumn%s'%(i)
				self.sampleIDList.append(sampleID)
				for j in xrange(self.ploidy):	#same samples in multiple columns
					self.header.append(sampleID)
		for i in xrange(len(self.header)):
			sampleID = self.header[i]
			if sampleID not in self.col_name2index:
				self.col_name2index[sampleID] = i	# only first haplotype is recorded if ploidy >1
			#self.col_name2index[sampleID].append(i)
			#the index corresponds to genotypeLikelihoodList in next().
			#so it starts from 0
		return self.col_name2index
	
	def getHaplotypeListOfOneSample(self, sampleID=None):
		"""
		2013.05.23
			return a list of haplotypes for one individual, from self.snpData
			the number of haplotypes fetched depends on self.ploidy
		"""
		self._completeSNPData()
		haplotypeList = []
		firstHaplotypeIndex = self.snpData.col_id2col_index.get(sampleID)
		for j in xrange(firstHaplotypeIndex, firstHaplotypeIndex+self.ploidy):
			#
			haplotypeList.append(self.snpData.data_matrix[:,j])
		return haplotypeList
	
	def getGenotypeOfOneSampleFromOneLocus(self, oneLocus=None, sampleID=None):
		"""
		2013.05.06
			oneLocus is output of next(), upon iterating through the file
		"""
		try:
			sampleStartIndex = self.getSampleIndexGivenSampleID(sampleID)
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
			sys.stderr.write("  sampleID %s not in file.\n"%(sampleID))
			raise
		genotypeList = oneLocus.genotypeList[sampleStartIndex:sampleStartIndex+self.ploidy]
		return genotypeList
	
	def getGenotypeOfOneSampleOneLocus(self, sampleID=None, locusID=None):
		"""
		2013.06.05
			get genotype from self.snpData
		"""
		self._completeSNPData()
		
		sampleStartIndex = self.getSampleIndexGivenSampleID(sampleID)
		locusIndex = self.snpData.row_id2row_index.get(locusID)
		sampleIndexList = []
		for i in xrange(sampleStartIndex, sampleStartIndex+self.ploidy):
			sampleIndexList.append(i)
		return self.snpData.data_matrix[locusIndex, sampleIndexList]
	
	def next(self):
		try:
			row = self.csvFile.next()
		except:
			raise StopIteration
		if not self.isRealCSV:
			row = row.strip().split()
		lineIndicator, markerID = row[0:2]
		
		return PassingData(markerID=markerID, lineIndicator=lineIndicator, genotypeList=row[2:])
	
	def readInAllHaplotypes(self, ):
		"""
		2013.05.23
			this function reads in the entire file and fills in these data structures
				self.sampleIDList
				self.locusIDList
				self.snpData
		"""
		sys.stderr.write("Reading in all haplotypes from %s ..."%(self.inputFname))
		if self.inputFile:
			if self.inputFile.tell()!=0:
				self.inputFile.seek(0)
		self.constructColName2IndexFromHeader()
		dataMatrix = []
		for oneLocusData in self:
			self.locusIDList.append(oneLocusData.markerID)
			dataMatrix.append(oneLocusData.genotypeList)
		
		self.snpData = SNPData(row_id_ls=self.locusIDList, col_id_ls=self.sampleIDList, \
							data_matrix=dataMatrix, ploidy=self.ploidy, matrix_orientation=2)
			#matrix_orientation=2, rows are loci, columns are samples
		self.sampleID2index = self.snpData.col_id2col_index
		sys.stderr.write("%s samples, %s loci.\n"%(len(self.sampleIDList), len(self.locusIDList)))
		
	def addOneIndividual(self, sampleID=None, locusIDList=None, haplotypeList=None):
		"""
		2013.05.30
			this is for writing mode
		"""
		self.sampleIDList.append(sampleID)
		self.haplotypeMatrix.append(haplotypeList)
		if locusIDList is not None:
			if not self.locusIDList:
				self.locusIDList = locusIDList
			elif self.locusIDList!=locusIDList:
				sys.stderr.write("Error: this individual %s has a different list of locus ID than this file (%s) already has %s.\n"%\
								(sampleID, self.inputFname, self.locusIDList))
				raise
			
	def _consolidateHaplotypeMatrixIntoSNPData(self, ):
		"""
		2013.06.04
			if self.sampleIDList, self.locusIDList, self.haplotypeMatrix are all available (from self.addOneIndividual()),
				but self.snpData is not
			the key step is to turn a list of lists into a numpy 2D array
		"""
		sys.stderr.write("Consolidating sampleIDList, locusIDList, haplotypeMatrix into SNPData ...")
		dataMatrix = numpy.array(self.haplotypeMatrix).transpose()	#transpose to make a column corresponding to haplotype 
		self.snpData = SNPData(row_id_ls=self.locusIDList, col_id_ls=self.sampleIDList, \
							data_matrix=dataMatrix, ploidy=self.ploidy, matrix_orientation=2)
			#matrix_orientation=2, rows are loci, columns are samples
		self.sampleID2index = self.snpData.col_id2col_index
		sys.stderr.write("%s samples, %s loci.\n"%(len(self.sampleIDList), len(self.locusIDList)))
	
	def writeDataToDisk(self):
		"""
		2013.05.30
		"""
		sys.stderr.write("Writing data to file %s ..."%(self.inputFname))
		if self.snpData is None:
			if self.sampleIDList and self.locusIDList and self.haplotypeMatrix:
				self._consolidateHaplotypeMatrixIntoSNPData()
			else:
				sys.stderr.write("Error: either sampleIDList (%s) or locusIDList (%s) or haplotypeMatrix is None/empty.\n"%\
								(repr(self.sampleIDList[:3]), repr(self.locusIDList[:3]),  ))
				raise
		self.writerow(['I', 'id'] + self.sampleIDList)
		for i in xrange(len(self.locusIDList)):
			locusID = self.locusIDList[i]
			genotypeList = self.haplotypeMatrix[i:]
			oneRow = ['M', locusID] + genotypeList
			self.writerow(oneRow)
		sys.stderr.write(" .\n")
	
if __name__ == '__main__':
	main_class = BeagleGenotypeFile
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()