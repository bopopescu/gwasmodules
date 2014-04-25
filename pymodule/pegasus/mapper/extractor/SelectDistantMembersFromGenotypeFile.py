#!/usr/bin/env python
"""
Examples:
	#2013.06.25
	m=0.2;
	%s  -i  folderHighCoveragePanel/Scaffold96_88778_VCF_87875_VCF_Scaffold96_splitVCF_u1.minCoverage8.beagled.vcf.gz 
		-o  folderHighCoveragePanel/Scaffold96_88778_VCF_87875_VCF_Scaffold96_splitVCF_u1.phasedRefPanel.m$m.tsv
		--maxPairwiseKinship $m --sampleSize 40
		--pedigreeKinshipFilePath ../Kinx2Sept2012.txt.gz
		--replicateIndividualTag copy
		--individualAlignmentCoverageFname  folderAuxilliary/88179_VCF_86965_VCF_Scaffold99.filterByMaxSNPMissingRate.recode.alignmentDepth.tsv
		--pedigreeFname  folderAuxilliary/pedigree.88179_VCF_86965_VCF_Scaffold99.filterByMaxSNPMissingRate.recode.format1.txt
	
	%s -o ... input1.vcf input2.vcf input3.vcf
	

Description:
	2013.06.05 select haplotypes of distant pedigree members. IBD distance of all pairs must be < --maxPairwiseKinship. 
		Input VCF files could be passed through "-i" and/or appended to the end of commandline.
			Unlimited number of them.
		--pedigreeFname provides the pedigree linking all samples.
			Pedigree file could include more samples (identical ID system) than those in input Beagle files.
				The extra will be removed from the graph before proceeding.
			The pedigree is better not to be split into trios/duos because an individual's number of offspring affects
				its probability of being sampled.
		Sample IDs in input Beagle and --pedigreeFname files could be in replicate form, containing --replicateIndividualTag, or not.
			first number is sample ID is alignment ID, matching samples in --individualAlignmentCoverageFname.
		Sample IDs in --pedigreeKinshipFilePath are individual.code.
		Sample IDs in outputFname would be in the same format as those of input.
		
		Two factors (linearly) affect a sample's probability of being chosen.
			1. sequence coverage
			2. number of offspring
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
from pymodule.yhio.VCFFile import VCFFile
from pymodule.statistics import NumberContainer, DiscreteProbabilityMassContainer

class SelectDistantMembersFromGenotypeFile(AbstractMatrixFileWalker):
	__doc__ = __doc__
	
	option_default_dict = AbstractMatrixFileWalker.option_default_dict
	option_default_dict.update({
			('sampleSize', 0, int): [40, '', 1, 'max number of individuals to be selected into output file.\n\
	The final number could be lower because maxPairwiseKinship is another threshold to meet.'],\
			('maxPairwiseKinship', 0, float): [0.2, '', 1, 'maximum pairwise kinship allowed among selected individuals.'],\
			('pedigreeKinshipFilePath', 1, ): [None, '', 1, 'file that contains pairwise kinship between individuals (ID: ucla_id/code).\n\
	no header. coma-delimited 3-column file: individual1, individual2, kinship\n\
	The sampling will try to avoid sampling close pairs, kinship(i,j)<=maxPairwiseKinship'],\
			('replicateIndividualTag', 0, ): ['copy', '', 1, 'the tag that separates the true ID and its replicate count'],\
			('individualAlignmentCoverageFname', 1, ): ['', '', 1, 'file contains two columns, individual-alignment.read_group, coverage.'],\
			('pedigreeFname', 1, ): ['', '', 1, 'pedigree (trios/duos/singletons) file that covers IDs in all beagle input files, but not more, in plink format.\n\
	This is used to figure out the family context of different replicates of the same individual => select one to represent all replicates'],\
			
			})
	option_default_dict[('outputFileFormat', 0, int)][0] = 4	#a no-header matrix
	
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
	
	def getIndividualCoverage(self, individualID=None, alignmentReadGroup2coverageLs=None):
		"""
		2013.05.24
		"""
		alignmentID = individualID.split('_')[0]
		#alignmentIDTuple = (alignmentID, )
		coverageLs = alignmentReadGroup2coverageLs.get(individualID)
		if not coverageLs:
			sys.stderr.write("Error: coverage for alignmentID=%s, %s is %s.\n"%(alignmentID, individualID, repr(coverageLs)))
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
		AbstractMatrixFileWalker.setup(self, **keywords)
		#self.writer = BeagleGenotypeFile(inputFname=self.outputFname, openMode='w')
		
		#read in the IBD check result
		self.ibdData = SNP.readAdjacencyListDataIntoMatrix(inputFname=self.pedigreeKinshipFilePath, \
								rowIDHeader=None, colIDHeader=None, \
								rowIDIndex=0, colIDIndex=1, \
								dataHeader=None, dataIndex=2, hasHeader=False)
		
		#. read in the alignment coverage data
		alignmentCoverageFile = MatrixFile(inputFname=self.individualAlignmentCoverageFname)
		alignmentCoverageFile.constructColName2IndexFromHeader()
		alignmentReadGroup2coverageLs = alignmentCoverageFile.constructDictionary(keyColumnIndexList=[0], valueColumnIndexList=[1])
		alignmentCoverageFile.close()
		
		sys.stderr.write("Reading in all samples from %s VCF input files ... \n"%(len(self.inputFnameLs)))
		# read all the Beagle files
		individualID2HaplotypeData = {}
		for inputFname in self.inputFnameLs:
			vcfFile = VCFFile(inputFname=inputFname)
			#vcfFile.readInAllHaplotypes()
			for individualID in vcfFile.getSampleIDList():
				individualID2HaplotypeData[individualID] = None
				#haplotypeList = vcfFile.getHaplotypeListOfOneSample(individualID)
				#individualID2HaplotypeData[individualID] = PassingData(haplotypeList=haplotypeList,
				#													locusIDList=vcfFile.locusIDList)
			# get all haplotypes , etc.
			# get all sample IDs
		sys.stderr.write("%s individuals total.\n"%(len(individualID2HaplotypeData)))
		
		#. read in the pedigree or deduce it from Beagle Trio/Duo genotype file (columns)
		#. construct individualID2pedigreeContext, context: familySize=1/2/3, familyPosition=1/2 (parent/child)
		sys.stderr.write("Constructing individualID2pedigreeContext ...")
		plinkPedigreeFile = PlinkPedigreeFile(inputFname=self.pedigreeFname)
		pGraph = plinkPedigreeFile.pedigreeGraph
		#shrink the graph to only individuals with data
		pGraph = nx.subgraph(pGraph, individualID2HaplotypeData.keys())
		
		cc_subgraph_list = nx.connected_component_subgraphs(pGraph.to_undirected())
		individualID2familyContext = {}
		outDegreeContainer = NumberContainer(minValue=0)
		familySizeContainer = NumberContainer(minValue=0)
		individualCoverageContainer = NumberContainer(minValue=0)
		familyCoverageContainer = NumberContainer(minValue=0)
		for cc_subgraph in cc_subgraph_list:
			familySize= len(cc_subgraph)
			familySizeContainer.addOneValue(familySize)
			
			familyCoverage = 0
			for n in cc_subgraph:	#assuming each family is a two-generation trio/nuclear family
				individualCoverage = self.getIndividualCoverage(individualID=n, alignmentReadGroup2coverageLs=alignmentReadGroup2coverageLs)
				individualCoverage = float(individualCoverage)
				individualCoverageContainer.addOneValue(individualCoverage)
				familyCoverage += individualCoverage
				in_degree = pGraph.in_degree(n)
				out_degree = pGraph.out_degree(n)
				outDegreeContainer.addOneValue(out_degree)
				familyContext = PassingData(familySize=familySize, in_degree=in_degree, out_degree=out_degree, \
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
		
		
		# weigh each unique individual based on its sequencing coverage + no of offspring => probability mass for each individual
		sys.stderr.write("Weighing each individual , assigning probability mass  ...")
		individualID2probabilityMass = {}
		for individualID, familyContext in individualID2familyContext.iteritems():
			outDegreeQuotient = outDegreeContainer.normalizeValue(familyContext.familySize)
			individualCoverageQuotient = individualCoverageContainer.normalizeValue(familyContext.individualCoverage)
			#familyCoverageQuotient = familyCoverageContainer.normalizeValue(familyContext.familyCoverage)
			importanceScore = outDegreeQuotient + individualCoverageQuotient
			representativeImportanceScore = importanceScore
			individualID2probabilityMass[individualID] = representativeImportanceScore
		sys.stderr.write(" %s IDs with probability mass assigned.\n"%(len(individualID2probabilityMass)))
		
		self.individualID2probabilityMass = individualID2probabilityMass
		self.individualID2HaplotypeData = individualID2HaplotypeData
	
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
	
	def mapSampleIDToIDInIBDFile(self, genotypeSampleIDList=None, ibdFileSampleIDList=None):
		"""
		2013.06.20
			sample IDs are alignment read_group, with or without replicateIndividualTag,
			ID in IBD file are individual.code (or maybe ucla_id)
			do pairwise matching , not relying on DB
		"""
		sys.stderr.write("Finding a map between %s genotype sample IDs and %s IBD sample IDs ...\n"%\
						(len(genotypeSampleIDList), len(ibdFileSampleIDList)))
		genotypeSampleID2IBDSampleID = {}
		for genotypeSampleID in genotypeSampleIDList:
			for ibdSampleID in ibdFileSampleIDList:
				if genotypeSampleID.find(ibdSampleID)>=0:
					if genotypeSampleID in genotypeSampleID2IBDSampleID:
						sys.stderr.write("Warning: %s already in genotypeSampleID2IBDSampleID with value=%s, overwritten with ibdSampleID=%s.\n"%\
										(genotypeSampleID, genotypeSampleID2IBDSampleID.get(genotypeSampleID), ibdSampleID))
					genotypeSampleID2IBDSampleID[genotypeSampleID] = ibdSampleID
		sys.stderr.write("\t %s pairs.\n"%(len(genotypeSampleID2IBDSampleID)))
		return genotypeSampleID2IBDSampleID
	
	def reduce(self, **keywords):
		"""
		2012.10.15
			run after all files have been walked through
		"""
		#sample the data
		probabilityMassContainer = DiscreteProbabilityMassContainer(object2proabilityMassDict=self.individualID2probabilityMass)
				
		noOfTotalRows = len(self.individualID2probabilityMass)
		genotypeSampleID2IBDSampleID = self.mapSampleIDToIDInIBDFile(genotypeSampleIDList=self.individualID2probabilityMass.keys(), \
															ibdFileSampleIDList=self.ibdData.row_id_ls)
		counter = 0
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
					counter += 1
					if sampledIndividualID:
						includeInTheSampling = True
						for alreadySampledIndividualID in sampledIndividualIDSet:	#not too close to anyone previously sampled
							#getting the relatedness
							relatedness = self.ibdData.getCellDataGivenRowColID(genotypeSampleID2IBDSampleID.get(sampledIndividualID), \
																			genotypeSampleID2IBDSampleID.get(alreadySampledIndividualID))
							if relatedness is not None and relatedness>=self.maxPairwiseKinship:
								includeInTheSampling = False
						if includeInTheSampling:
							sampledIndividualIDSet.add(sampledIndividualID)
					sampledSetSizeHistoryData.historyList.append(len(sampledIndividualIDSet))
				#turn into list
				sampledIndividualIDList = list(sampledIndividualIDSet)
			else:
				sampledIndividualIDList = random.sample(self.individualID2probabilityMass.keys(), self.sampleSize)
		else:	#take all
			sampledIndividualIDList = self.individualID2probabilityMass.keys()
		
		#output the sampled individuals
		for individualID in sampledIndividualIDList:
			self.writer.writerow([individualID])
			real_counter += 1
		
		fraction = float(real_counter)/float(noOfTotalRows)
		sys.stderr.write("%s/%s (%.3f) selected out of %s samplings.\n"%(real_counter, noOfTotalRows, fraction, counter))
		
		#close the self.invariantPData.writer and self.writer
		AbstractMatrixFileWalker.reduce(self, **keywords)


if __name__ == '__main__':
	main_class = SelectDistantMembersFromGenotypeFile
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
