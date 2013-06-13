#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -o ... input1.bgl input2.bgl input3.bgl
	

Description:
	2013.06.05 select distant members from 
		input Beagle phased files should be added to the end of commandline, as many as you want.
			Whether it is Singleton/Trio/Duo beagle files, it does not matter. Pedigree is extracted from plinkIBDCheckOutputFname.
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
from pymodule.yhio import SNP
from pymodule.yhio.MatrixFile import MatrixFile
from pymodule.yhio.PlinkPedigreeFile import PlinkPedigreeFile
from pymodule.statistics import NumberContainer, DiscreteProbabilityMassContainer

class SelectDistantMembersFromGenotypeFile(AbstractMatrixFileWalker):
	__doc__ = __doc__
	
	option_default_dict = AbstractMatrixFileWalker.option_default_dict
	option_default_dict.update({
			('sampleSize', 0, int): [40, '', 1, 'max number of individuals to be selected into output file.\n\
	It could be low this number because maxPairwiseKinship is another threshold that is required.'],\
			('maxPairwiseKinship', 0, float): [0.2, '', 1, 'maximum pairwise kinship allowed among selected individuals.'],\
			('plinkIBDCheckOutputFname', 1, ): [None, '', 1, 'file that contains IBD check result, PI_HAT=relatedness.\n\
	at least 3-columns with header: IID1, IID2, PI_HAT. IID1 and IID2 should match the whichColumn (whichColumnHeader) of inputFname.\n\
	The sampling will try to avoid sampling close pairs, PI_HAT(i,j)<=maxIBDSharing'],\
			('replicateIndividualTag', 0, ): ['copy', '', 1, 'the tag that separates the true ID and its replicate count'],\
			('individualAlignmentCoverageFname', 1, ): ['', '', 1, 'file contains two columns, individual-alignment-ID, coverage.'],\
			('pedigreeFname', 1, ): ['', '', 1, 'pedigree file that covers IDs in all beagle input files, but not more, in plink format.\n\
	This is used to figure out the family context for different replicates of the same individual => used to select one final replicate to represent the whole'],\
			
			})
	
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractMatrixFileWalker.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def getOriginalIndividualID(self, individualID=None, replicateIndividualTag='copy'):
		"""
		2013.05.24
		"""
		originalIndividualID = individualID.split(replicateIndividualTag)[0]
		return originalIndividualID
	
	def getIndividualCoverage(self, individualID=None, alignmentIDTuple2coverageLs=None):
		"""
		2013.05.24
		"""
		alignmentID = individualID.split('_')[0]
		alignmentIDTuple = (alignmentID, )
		coverageLs = alignmentIDTuple2coverageLs.get(alignmentIDTuple)
		if not coverageLs:
			sys.stderr.write("Error: coverage for %s is %s.\n"%(individualID, repr(coverageLs)))
			sys.exit(3)
		if len(coverageLs)>1:
			sys.stderr.write("Warning: coverage for %s has more than 1 entries %s. Take first one.\n "%\
							(individualID, repr(coverageLs)))
		return coverageLs[0]
	
	def setup(self, **keywords):
		"""
		2012.10.15
			run before anything is run
		"""
		#2013.05.30 comment out AbstractMatrixFileWalker.setup() to open the output file differently
		#AbstractMatrixFileWalker.setup(self, **keywords)
		self.writer = BeagleGenotypeFile(inputFname=self.outputFname, openMode='w')
		
		#read in the IBD check result
		self.ibdData = SNP.readAdjacencyListDataIntoMatrix(inputFname=self.plinkIBDCheckOutputFname, \
								rowIDHeader="IID1", colIDHeader="IID2", \
								rowIDIndex=None, colIDIndex=None, \
								dataHeader="PI_HAT", dataIndex=None, hasHeader=True)
		
		#. read in the alignment coverage data
		alignmentCoverageFile = MatrixFile(inputFname=self.individualAlignmentCoverageFname)
		alignmentCoverageFile.constructColName2IndexFromHeader()
		alignmentIDTuple2coverageLs = alignmentCoverageFile.constructDictionary(keyColumnIndexList=[0], valueColumnIndexList=[1])
		alignmentCoverageFile.close()
		
		
		#. read in the pedigree or deduce it from Beagle Trio/Duo genotype file (columns)
		#. construct individualID2pedigreeContext, context: familySize=1/2/3, familyPosition=1/2 (parent/child)
		plinkPedigreeFile = PlinkPedigreeFile(inputFname=self.pedigreeFname)
		pGraphStructure = plinkPedigreeFile.getPedigreeGraph()
		sys.stderr.write("Constructing individualID2pedigreeContext ...")
		pGraph = pGraphStructure.DG
		cc_subgraph_list = nx.connected_component_subgraphs(pGraph)
		individualID2familyContext = {}
		familySizeContainer = NumberContainer(minValue=0)
		individualCoverageContainer = NumberContainer(minValue=0)
		familyCoverageContainer = NumberContainer(minValue=0)
		for cc_subgraph in cc_subgraph_list:
			familySize= len(cc_subgraph)
			familySizeContainer.addOneValue(familySize)
			
			familyCoverage = None
			for n in cc_subgraph:	#assuming each family is a two-generation trio/nuclear family
				individualCoverage = self.getIndividualCoverage(individualID=n, alignmentIDTuple2coverageLs=alignmentIDTuple2coverageLs)
				individualCoverageContainer.addOneValue(individualCoverage)
				familyCoverage += individualCoverage
				in_degree = cc_subgraph.in_degree(n)
				out_degree = cc_subgraph.out_degree(n)
				if in_degree==2:
					familyPosition = 3
				elif out_degree==1:
					familyPosition=1	#parent
				familyContext = PassingData(familySize=familySize, familyPosition=familyPosition, \
										individualCoverage=individualCoverage,\
										familyCoverage=None)
				if n not in individualID2familyContext:
					individualID2familyContext[n] = familyContext
				else:
					sys.stderr.write("Node %s already in individualID2familyContext.\n"%(n))
			familyCoverageContainer.addOneValue(familyCoverage)
			#set the family coverage for each member, used in weighing the individual. better covered family => better haplotype
			for n in cc_subgraph:
				individualID2familyContext[n].familyCoverage = familyCoverage
		plinkPedigreeFile.close()
		sys.stderr.write("%s individuals.\n"%(len(individualID2familyContext)))
		
		# read all the Beagle files
		individualID2HaplotypeData = {}
		for inputFname in self.inputFnameLs:
			beagleFile = BeagleGenotypeFile(inputFname=inputFname)
			beagleFile.readInAllHaplotypes()
			for individualID in beagleFile.sampleIDList:
				haplotypeList = beagleFile.getHaplotypeListOfOneSample(individualID)
				individualID2HaplotypeData[individualID] = PassingData(haplotypeList=haplotypeList,
																	locusIDList=beagleFile.locusIDList)
			# get all haplotypes , etc.
			# get all sample IDs
		
		
		sys.stderr.write("Getting originalIndividualID2individualIDList ...")
		originalIndividualID2individualIDList = {}
		for individualID in individualID2HaplotypeData:
			originalIndividualID = self.getOriginalIndividualID(individualID, self.replicateIndividualTag)
			if originalIndividualID not in originalIndividualID2individualIDList:
				originalIndividualID2individualIDList[originalIndividualID] = []
			originalIndividualID2individualIDList[originalIndividualID].append(individualID)
		
		sys.stderr.write("%s original IDs from %s individuals.\n"%(len(originalIndividualID2individualIDList),\
																len(individualID2HaplotypeData)))
		
		sys.stderr.write("Selecting representative individual for each replicate group, Weighing each individual , assigning probability mass  ...")
		# select representative for each replicate group, trio > duo > singleton
		# Note: if an individual is replicated, it will not show up as singleton at all.
		#  			so the candidates lie in trios, duos
		#  importance score = familySize + individual-coverage + coverage of all members
		#		equal weight
		# weigh each unique individual based on its sequencing coverage => probability mass for each individual
		originalIndividualID2representativeData = {}
		for originalIndividualID, individualIDList in originalIndividualID2individualIDList.iteritems():
			representativeID = None
			representativeImportanceScore = None
			for individualID in individualIDList:
				familyContext = individualID2familyContext.get(individualID)
				familySizeQuotient = familySizeContainer.normalizeValue(familyContext.familySize)
				individualCoverageQuotient = individualCoverageContainer.normalizeValue(familyContext.individualCoverage)
				familyCoverageQuotient = familyCoverageContainer.normalizeValue(familyContext.familyCoverage)
				importanceScore = familySizeQuotient + individualCoverageQuotient + familyCoverageQuotient
				if representativeImportanceScore is None or representativeImportanceScore < importanceScore:
					representativeImportanceScore = importanceScore
					representativeID = individualID
			originalIndividualID2representativeData[originalIndividualID] = PassingData(\
																representativeImportanceScore=representativeImportanceScore,\
																representativeID=representativeID)
		sys.stderr.write(" \n")
		
		self.originalIndividualID2representativeData = originalIndividualID2representativeData
		self.individualID2HaplotypeData = individualID2HaplotypeData
		self.originalIndividualID2individualIDList = originalIndividualID2individualIDList
	
	def processRow(self, row=None, pdata=None):
		"""
		2012.10.7
		returnValue = 1
		self.data_matrix.append(row)
		self.col_name2index = getattr(pdata, 'col_name2index', None)
		
		col_name2index = getattr(pdata, 'col_name2index', None)
		y_ls = getattr(pdata, 'y_ls', None)	#don't add anything to it, then self.plot won't be called in self.fileWalker
		if col_name2index and y_ls is not None:
			if self.whichColumnHeader:
				whichColumn = col_name2index.get(self.whichColumnHeader, None)
			else:
				whichColumn = self.whichColumn
			
			yValue = self.handleYValue(row[whichColumn])
			if self.minWhichColumnValue is not None and yValue<self.minWhichColumnValue:
				return
			if self.maxWhichColumnValue is not None and yValue>self.maxWhichColumnValue:
				return
			self.invariantPData.writer.writerow(row)
		"""
		
		return 0
	
	def detectSampledSetSizeHistoryChangeInLastRounds(self, sampledSetSizeHistoryData=None):
		"""
		2013.05.30
		"""
		historyList = sampledSetSizeHistoryData.historyList
		noOfLastRounds = sampledSetSizeHistoryData.noOfLastRounds
		if len(historyList)<=1:	#can't do anything. too short.
			return 1
		
		#1. add the last step difference into the sum 
		sampledSetSizeHistoryData.sumOfAbsStepDifference += abs(historyList[-1]-historyList[-2])
		if len(historyList)<noOfLastRounds:	#hasn't reached that the minimal number of rounds
			return 1	#True
		else:
			#shorten the list to of length (noOfLastRounds+3)
			historyList = historyList[-noOfLastRounds-3:]
			if len(historyList)>noOfLastRounds:	#need to kick out the rounds before the last few rounds
				#this is to make sure the sum only contains the step difference within the last few rounds
				sampledSetSizeHistoryData.sumOfAbsStepDifference = sampledSetSizeHistoryData.sumOfAbsStepDifference - \
					abs(historyList[-noOfLastRounds]-historyList[-noOfLastRounds-1])
			if sampledSetSizeHistoryData.sumOfAbsStepDifference==0:
				return 0
			else:
				return 1
	
	def reduce(self, **keywords):
		"""
		2012.10.15
			run after all files have been walked through
		"""
		#sample the data
		probabilityMassContainer = DiscreteProbabilityMassContainer(object2proabilityMassDict=self.originalIndividualID2representativeData)
		
		self.originalIndividualID2individualIDList.keys()
		
		noOfTotalRows = len(self.originalIndividualID2representativeData)
		real_counter = 0
		if self.sampleSize<noOfTotalRows:
			if self.ibdData:
				#complicated sampling starts here
				#
				sampledSetSizeHistoryData = PassingData(historyList= [], sumOfAbsStepDifference = 0, \
													noOfLastRounds=20)
					#a metre about whether sampledIndividualIDSet has stopped growing
				sampledIndividualIDSet = set()
				while len(sampledIndividualIDSet)<self.sampleSize and \
						self.detectSampledSetSizeHistoryChangeInLastRounds(sampledSetSizeHistoryData=sampledSetSizeHistoryData):
					sampledIndividualID = probabilityMassContainer.sampleObject()
					if sampledIndividualID:
						includeInTheSampling = True
						for alreadySampledIndividualID in sampledIndividualIDSet:	#not too close to anyone previously sampled
							#getting the relatedness
							relatedness = self.ibdData.getCellDataGivenRowColID(sampledIndividualID, alreadySampledIndividualID)
							if relatedness>=self.maxIBDSharing:
								includeInTheSampling = False
						if includeInTheSampling:
							sampledIndividualIDSet.add(sampledIndividualID)
					sampledSetSizeHistoryData.historyList.append(len(sampledIndividualIDSet))
				#turn into list
				sampledIndividualIDList = list(sampledIndividualIDSet)
			else:
				sampledIndividualIDList = random.sample(self.originalIndividualID2representativeData.keys(), self.sampleSize)
		else:	#take all
			sampledIndividualIDList = self.originalIndividualID2representativeData.keys()
		
		#output the sampled individuals
		for originalIndividualID in sampledIndividualIDList:
			representativeData = self.originalIndividualID2representativeData.get(originalIndividualID)
			haplotypeData = self.individualID2HaplotypeData.get(representativeData.representativeID)
			self.writer.addOneIndividual(sampleID=originalIndividualID, haplotypeList=haplotypeData.haplotypeList, \
										locusIDList=haplotypeData.locusIDList)
			real_counter += 1
		
		self.writer.writeDataToDisk()	#first write out header (sampleID), then one locus a line
		fraction = float(real_counter)/float(noOfTotalRows)
		sys.stderr.write("%s/%s (%.3f) selected.\n"%(real_counter, noOfTotalRows, fraction))
		
		#close the self.invariantPData.writer and self.writer
		AbstractMatrixFileWalker.reduce(self, **keywords)


if __name__ == '__main__':
	main_class = SelectDistantMembersFromGenotypeFile
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
