#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s  -i --gatkPrintBeagleFname method_36_Contig791_replicated.bgl 
		--plinkPedigreeFname trioCaller_VRC.merlin -O method_36_Contig791_replicated_unphased
	
Description:
	2013.05.03
		This program splits a VCF file from pedigree members into trios/duos beagle input files.
			Pedigree is provided by --plinkPedigreeFname and contains individuals from --gatkPrintBeagleFname.
				If more (a pedigree that encompasses more individuals from --gatkPrintBeagleFname), the extra will be ignored.
				If less (individuals in --gatkPrintBeagleFname but not in --plinkPedigreeFname), would be ignored.
				Output is the intersection between the two. 
			The --gatkPrintBeagleFname file (beagle likelihood file) is output of running GATK PrintBeagleInput on VCF.
		Output includes
			1 or 2 or 3 (depends on composition of pedigree) beagle files 'outputFnamePrefix_familySize??.bgl
				If a pedigree is composed entirely of one type of family (trio/duo/singleton), then only one beagle file.
				The Beagle file for singletons is in likelihood format.
			One additional output file is outputFnamePrefix.markers, which contains the markers.
				Contig791:1086 1086 C A
	
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy, numpy
import networkx as nx
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.yhio.VCFFile import VCFFile
from pymodule.yhio.MatrixFile import MatrixFile
from pymodule.yhio.BeagleLikelihoodFile import BeagleLikelihoodFile
from pymodule.yhio.PlinkPedigreeFile import PlinkPedigreeFile
from pymodule.pegasus.mapper.AbstractVCFMapper import AbstractVCFMapper

class SplitPedigreeVCFIntoBeagleTriosDuosFiles(AbstractVCFMapper):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractVCFMapper.option_default_dict)
	#option_default_dict.pop(("chromosome", 0, ))
	#option_default_dict.pop(("chrLength", 1, int))
	option_default_dict.pop(('outputFname', 0, ))
	#option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						
						('gatkPrintBeagleFname', 0, ): [None, '', 1, 'beagle likelihood file, output of running GATK PrintBeagleInput on the VCF (inputFname).', ],\
						('plinkPedigreeFname', 0, ): [None, '', 1, 'pedigree filename in plink format. tab or space delimiter is fine', ],\
						('minProbForValidCall', 0, float): [0.7, '', 1, 'minimum probability (from gatkPrintBeagleFname) for the most likely genotype (among 3) to be valid,\n\
	not marked as missing ?', ],\
						('dummyIndividualNamePrefix', 1, ): ['dummy', '', 1, 'the prefix to name a dummy parent (TrioCaller format). The suffix is its order among all dummies.'],\
						#('noOfTotalSites', 1, int): [None, 'n', 1, 'The total number of sites in input VCF', ],\
						})
	def __init__(self,  **keywords):
		"""
		"""
		AbstractVCFMapper.__init__(self, **keywords)
	
	def constructPedigreeFamilyData(self, plinkPedigreeFname=None, gatkPrintBeagleFname=None, dummyIndividualNamePrefix='dummy'):
		"""
		
			
			gatkPrintBeagleFname is the output of this command:
				java -Xmx2g -jar ~/script/gatk2/GenomeAnalysisTK.jar
					-R /Network/Data/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta
					-T ProduceBeagleInput
					-V /Network/Data/vervet/db/genotype_file/method_36/36860_VCF_33518_VCF_26869_VCF_24449_VCF_Contig791.subset.vcf.gz
					-o method_36_Contig791.bgl
				
				.bgl format:
					marker alleleA alleleB 1000_709_1996093_GA_vs_524 1000_709_1996093_GA_vs_524 1000_709_1996093_GA_vs_524 1001_710_1995025_GA_vs_524 1001_710_1995025_GA_vs_524 1001_710_1995025_GA_vs_524 1002_711_2001039_GA_vs_524
					Contig791:1086 C A 0.9693 0.0307 0.0000 0.6660 0.3338 0.0003 0.0000
					Contig791:1649 G C 0.9406 0.0594 0.0000 0.9693 0.0307 0.0000 0.0000
					Contig791:4084 A C 0.9980 0.0020 0.0000 0.9844 0.0156 0.0000 0.0000
		"""
		sys.stderr.write("Constructing pedigree-family data from pedigree file %s, vcf file %s ...."%\
						(plinkPedigreeFname, gatkPrintBeagleFname))
		
		beagleLikelihoodFile = BeagleLikelihoodFile(inputFname=gatkPrintBeagleFname)
		sampleID2ColIndexList = beagleLikelihoodFile.constructColName2IndexFromHeader()
		#read in the pedigree, make sure samples exist in gatkPrintBeagleFname	
		counter = 0
		real_counter = 0
		plinkPedigreeFile = PlinkPedigreeFile(inputFname=plinkPedigreeFname, dummyIndividualNamePrefix=dummyIndividualNamePrefix)
		pedigreeGraph = plinkPedigreeFile.getPedigreeGraph().DG
		# shrink the pedigree to only individuals that exist in beagleLikelihoodFile
		pedigreeGraph = nx.subgraph(pedigreeGraph, sampleID2ColIndexList.keys())
		
		familySize2SampleIDList = {}
		#use pedigreeGraph above to figure out
		for node in pedigreeGraph:
			if node in sampleID2ColIndexList:
				no_of_incoming_edges = pedigreeGraph.in_degree(node)
				no_of_outgoing_edges = pedigreeGraph.out_degree(node)
				memberList = []
				familySize = None
				if no_of_incoming_edges==0 and no_of_outgoing_edges==0:	#check the outgoing as well because founders have no incoming edges.
					familySize = 1
				elif no_of_incoming_edges==1:
					familySize = 2
					memberList.extend(pedigreeGraph.predecessors(node))
				elif no_of_incoming_edges==2:
					familySize = 3
					memberList.extend(pedigreeGraph.predecessors(node))
				if familySize is not None:
					memberList.append(node)
					if familySize not in familySize2SampleIDList:
						familySize2SampleIDList[familySize] = []
					
					familySize2SampleIDList[familySize].extend(memberList)
		
		sys.stderr.write("\tFamilySize\tNoOfIndividuals\n")
		for familySize, sampleIDList in familySize2SampleIDList.iteritems():
			sys.stderr.write("\t%s\t%s\n"%(familySize, len(sampleIDList)))
		return PassingData(familyID2MemberList=None, familySize2SampleIDList=familySize2SampleIDList,\
						beagleLikelihoodFile=beagleLikelihoodFile)
	
	def openWriteBeagleFiles(self, pedigreeFamilyData=None, outputFnamePrefix=None):
		"""
		2013.05.02
			
		The non-likelihood (unphased, trios, pairs) Beagle format:
			I id sample1 sample1 sample2 sample2
			A diabetes 1 1 2 2
			M rs12082861 C C C C
			M rs4912233 T C C C
			M rs12732823 G A A A
			M rs17451521 C C C C
			M rs12033358 C T T T
		
		The likelihood version is
			marker alleleA alleleB 1000_709_1996093_GA_vs_524 1000_709_1996093_GA_vs_524 1000_709_1996093_GA_vs_524 1001_710_1995025_GA_vs_524 1001_710_1995025_GA_vs_524 1001_710_1995025_GA_vs_524 1002_711_2001039_GA_vs_524
			Contig791:1086 C A 0.9693 0.0307 0.0000 0.6660 0.3338 0.0003 0.0000
			Contig791:1649 G C 0.9406 0.0594 0.0000 0.9693 0.0307 0.0000 0.0000
			Contig791:4084 A C 0.9980 0.0020 0.0000 0.9844 0.0156 0.0000 0.0000
		
		The markers file has this format (markerID, position, alleleA, alleleB)
			Contig791:1086 1086 C A
		"""
		sys.stderr.write("Opening beagle files (outputFnamePrefix =%s) to write ..."%(outputFnamePrefix))
		familySize2BeagleFileHandler = {}
		familySize2SampleIDList =  pedigreeFamilyData.familySize2SampleIDList
		counter = 0
		for familySize, sampleIDList in familySize2SampleIDList.iteritems():
			if familySize not in familySize2BeagleFileHandler:
				tmpOutputFnamePrefix = '%s_familySize%s'%(outputFnamePrefix, familySize)
				writer = MatrixFile(inputFname='%s.bgl'%(tmpOutputFnamePrefix), openMode='w', delimiter=' ')
				familySize2BeagleFileHandler[familySize] = writer
				if familySize==1:
					headerRow = ['marker', 'alleleA', 'alleleB']
				else:
					headerRow = ['I', 'id'] 
				for sampleID in sampleIDList:
					if familySize==1:	#likelihood format has sample name replicated three times, rather than 2 times
						headerRow.extend([sampleID]*3)
					else:
						headerRow.extend([sampleID]*2)
				writer.writeHeader(headerRow)
				counter += 1
		markersFile = MatrixFile(inputFname='%s.markers'%(outputFnamePrefix), openMode='w', delimiter=' ')
		
		counter += 1
		sys.stderr.write("%s files outputted.\n"%(counter))
		
		return PassingData(familySize2BeagleFileHandler=familySize2BeagleFileHandler, markersFile=markersFile)
	
	def checkConcordanceBetweenBeagleAndVCFCall(self, vcfGenotypeCall=None, diploidCallFromBeagle=None):
		"""
		2013.05.06
			vcfGenotypeCall is a string, like 'NA', 'AA,'AG', etc.
			diploidCallFromBeagle is a list of two calls, ['?', '?'], ['A', 'A'], ['A', 'G']
		"""
		beagleCallInOneString1 = '%s%s'%(diploidCallFromBeagle[0], diploidCallFromBeagle[1])
		beagleCallInOneString2 = '%s%s'%(diploidCallFromBeagle[1], diploidCallFromBeagle[0])
		if vcfGenotypeCall!=beagleCallInOneString1 and vcfGenotypeCall!=beagleCallInOneString2:
			return False
		else:
			return True
	
	def splitVCFIntoBeagleInputs(self, inputFname=None, beagleLikelihoodFile=None, \
						familySize2BeagleFileHandler=None, pedigreeFamilyData=None, \
						minProbForValidCall=0.9, markersFile=None):
		"""
		2013.05.03
		
		The non-likelihood (unphased, trios, pairs) Beagle format:
			I id sample1 sample1 sample2 sample2
			A diabetes 1 1 2 2
			M Contig791:1086 C C C C
			M Contig791:1649 T C C C
			M Contig791:4084 G A A A
		"""
		sys.stderr.write("Splitting VCFFile %s (+ one beagle Likelihood file %s) into Beagle trios/duos files, minProbForValidCall=%s ... \n"%\
						(inputFname, beagleLikelihoodFile.inputFname, minProbForValidCall))
		counter = 0
		no_of_trios = 0
		no_of_duos = 0
		no_of_singletons = 0
		totalNoOfCalls = 0
		noOfCallsMarkedMissing = 0
		vcfFile = VCFFile(inputFname=inputFname)
		familySize2SampleIDList = pedigreeFamilyData.familySize2SampleIDList
		
		for vcfRecord in vcfFile:
			oneLocus = beagleLikelihoodFile.next()
			counter += 1
			familySize2CallList = {}
			genotypeLikelihoodList = oneLocus.genotypeLikelihoodList
			for familySize, sampleIDList in familySize2SampleIDList.iteritems():
				if familySize not in familySize2CallList:
					familySize2CallList[familySize] = []
				for sampleID in sampleIDList:
					totalNoOfCalls += 1
					vcfGenotypeCallData = vcfRecord.getGenotypeCallForOneSample(sampleID)
					tripleLikelihood = beagleLikelihoodFile.getLikelihoodListOfOneGenotypeOneSample(oneLocus=oneLocus, sampleID=sampleID)
					if familySize==1:
						no_of_singletons += 1
						familySize2CallList[familySize].extend(tripleLikelihood)
					else:
						if familySize==2:
							no_of_duos += 1
						elif familySize==3:
							no_of_trios += 1
						tripleLikelihood = map(float, tripleLikelihood)
						maxLikelihoodIndex = numpy.argmax(tripleLikelihood)
						maxLikelihood = tripleLikelihood[maxLikelihoodIndex]
						if maxLikelihood>=minProbForValidCall:
							if maxLikelihoodIndex==0:
								diploidCallFromBeagle = [oneLocus.alleleA, oneLocus.alleleA]
							elif maxLikelihoodIndex==1:
								diploidCallFromBeagle = [oneLocus.alleleA, oneLocus.alleleB]
							else:
								diploidCallFromBeagle = [oneLocus.alleleB, oneLocus.alleleB]
						else:
							noOfCallsMarkedMissing += 1
							diploidCallFromBeagle = ['?', '?']
						#if vcfGenotypeCallData is None:	#DP is zero
						#	sys.stderr.write("vcfGenotypeCallData for sample %s at locus %s, %s is None.\n"%\
						#					(sampleID, vcfRecord.chr, vcfRecord.pos))
						#	import pdb
						#	pdb.set_trace()
						if vcfGenotypeCallData and self.checkConcordanceBetweenBeagleAndVCFCall(vcfGenotypeCallData['GT'], diploidCallFromBeagle):
								diploidCall = [vcfGenotypeCallData['GT'][0], vcfGenotypeCallData['GT'][1]]
						else:
							diploidCall = ['?', '?']
						familySize2CallList[familySize].extend(diploidCall)
			
			for familySize, callList in familySize2CallList.iteritems():
				if familySize==1:
					rowHeaderList = [oneLocus.markerID, oneLocus.alleleA, oneLocus.alleleB]
				else:
					rowHeaderList = ['M', oneLocus.markerID]
				beagleFileHandler = familySize2BeagleFileHandler[familySize]
				
				beagleFileHandler.writerow(rowHeaderList+callList)
			if markersFile is not None:
				markersFile.writerow([oneLocus.markerID, oneLocus.markerID.split(':')[1], oneLocus.alleleA, oneLocus.alleleB])
		vcfFile.close()
		sys.stderr.write("%s loci, total %s calls, %s calls for singletons, %s calls for duos, %s calls for trios. %s calls marked missing.\n"%\
						(counter, totalNoOfCalls, no_of_singletons, no_of_duos, no_of_trios, noOfCallsMarkedMissing))
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		pedigreeFamilyData = self.constructPedigreeFamilyData(plinkPedigreeFname=self.plinkPedigreeFname, gatkPrintBeagleFname=self.gatkPrintBeagleFname)
		beagleFileData = self.openWriteBeagleFiles(pedigreeFamilyData=pedigreeFamilyData, outputFnamePrefix=self.outputFnamePrefix)
		
		self.splitVCFIntoBeagleInputs(inputFname=self.inputFname, beagleLikelihoodFile=pedigreeFamilyData.beagleLikelihoodFile, \
									familySize2BeagleFileHandler=beagleFileData.familySize2BeagleFileHandler, \
									pedigreeFamilyData=pedigreeFamilyData, minProbForValidCall=self.minProbForValidCall,\
									markersFile=beagleFileData.markersFile)
		
		#close all output files
		for familySize, beagleFileHandler in beagleFileData.familySize2BeagleFileHandler.iteritems():
			beagleFileHandler.close()
		beagleFileData.markersFile.close()
		pedigreeFamilyData.beagleLikelihoodFile.close()
		

if __name__ == '__main__':
	main_class = SplitPedigreeVCFIntoBeagleTriosDuosFiles
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()