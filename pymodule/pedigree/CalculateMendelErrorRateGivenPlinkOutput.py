#!/usr/bin/env python
"""
Examples:
	%s  -i mendelRuntype3_s1Gzip/meanMedianModePerLocusMendelError.tsv.gz
		--pedigreeFname vcf2plinkRuntype3_s1VCF2PlinkMerged/pedigree.tfam
		-o mendelRuntype3_s1Gzip/meanMendelErrorPerLocusPerFamily.tsv
		
	%s
	

Description:
	2013.07.19 a program that calculates the number of nuclear families (plink definition) and divide above meanMendelError with that
		# => meanMendelErrorRate per site
		# the number of nuclear families => the number of unique parents pairs in the tfam file (both parents are not missing)
		
input 1 is output of CalculateMedianMeanOfInputColumn
input 2 is .tfam file
input 3 is column stat name from input 1. should be meanMendelError.
output: tsv file with the divided rate (meanMendelError/#nuclear families)
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
from pymodule.yhio.MatrixFile import MatrixFile
from pymodule.yhio.PlinkPedigreeFile import PlinkPedigreeFile
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper

parentClass = AbstractMapper
class CalculateMendelErrorRateGivenPlinkOutput(parentClass):
	__doc__ = __doc__
	
	option_default_dict = parentClass.option_default_dict
	option_default_dict.update({
			('inputMendelErrorColumnHeader', 0, ): ["meanMendelError", '', 1, 'stat column header in inputFname, \n\
	to be divided by the number of nuclear families with both parents in pedigree'],\
			('pedigreeFname', 1, ): ['', '', 1, 'pedigree (trios/duos/singletons) file that covers IDs in all beagle input files, but not more, in plink format.\n\
	This is used to figure out the family context of different replicates of the same individual => select one to represent all replicates'],\
			
			})
	#option_default_dict[('outputFileFormat', 0, int)][0] = 4	#a no-header matrix
	
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		parentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def getNoOfFamiliesAndKidsGivenParentSetSize(self, noOfParents2FamilyData=None, parentSetSize=2):
		"""
		2013.07.19
		"""
		familyData = noOfParents2FamilyData.get(parentSetSize, None)
		
		if familyData:
			noOfFamilies = len(familyData.parentTupleSet)
			noOfParents = len(familyData.parentIDSet)
			noOfKids = len(familyData.childIDSet)
		else:
			noOfFamilies = 0
			noOfParents = 0
			noOfKids = 0
		return PassingData(noOfFamilies=noOfFamilies, noOfParents=noOfParents, noOfKids=noOfKids)

	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		reader = MatrixFile(inputFname=self.inputFname)
		reader.constructColName2IndexFromHeader()
		
		meanMendelErrorIndex = reader.getColIndexGivenColHeader("meanMendelError")
		noOfLociIndex = reader.getColIndexGivenColHeader("sampled_base_count")
		sumOfMendelErrorIndex = reader.getColIndexGivenColHeader("sumOfMendelError")
		
		plinkPedigreeFile = PlinkPedigreeFile(inputFname=self.pedigreeFname)
		familyStructureData = plinkPedigreeFile.getFamilyStructurePlinkWay()
		
		twoParentFamilyCountData = self.getNoOfFamiliesAndKidsGivenParentSetSize(noOfParents2FamilyData=familyStructureData.noOfParents2FamilyData, \
																		parentSetSize=2)
		singleParentFamilyCountData = self.getNoOfFamiliesAndKidsGivenParentSetSize(noOfParents2FamilyData=familyStructureData.noOfParents2FamilyData, \
																		parentSetSize=1)
		zeroParentFamilyCountData = self.getNoOfFamiliesAndKidsGivenParentSetSize(noOfParents2FamilyData=familyStructureData.noOfParents2FamilyData, \
																		parentSetSize=0)
		
		writer = MatrixFile(self.outputFname, openMode='w', delimiter='\t')
		header = ["ID", "noOfTotalLoci", "noOfTwoParentFamilies", "noOfParentsInTwoParentFamilies", "noOfKidsInTwoParentFamilies", \
				"noOfSingleParentFamilies", "noOfParentsInSingleParentFamilies", "noOfKidsInSingleParentFamilies", \
				"noOfZeroParentFamilies", "noOfParentsInZeroParentFamilies", "noOfKidsInZeroParentFamilies", \
				"noOfTotalMendelErrors", \
				"mendelErrorRatePerLocusPerNuclearFamily", "mendelErrorRatePerNuclearFamily"]
		writer.writeHeader(header)
		for row in reader:
			meanMendelError = float(row[meanMendelErrorIndex])
			noOfLoci = int(row[noOfLociIndex])
			sumOfMendelError = int(row[sumOfMendelErrorIndex])
			noOfNuclearFamilies = twoParentFamilyCountData.noOfFamilies
			data_row = [row[0], noOfLoci, noOfNuclearFamilies, twoParentFamilyCountData.noOfParents, twoParentFamilyCountData.noOfKids, \
					singleParentFamilyCountData.noOfFamilies,  singleParentFamilyCountData.noOfParents,  singleParentFamilyCountData.noOfKids,\
					zeroParentFamilyCountData.noOfFamilies, zeroParentFamilyCountData.noOfParents,  zeroParentFamilyCountData.noOfKids,\
					sumOfMendelError, \
					meanMendelError/float(noOfNuclearFamilies), sumOfMendelError/float(noOfNuclearFamilies)]
			writer.writerow(data_row)
		
		plinkPedigreeFile.close()
		reader.close()
		writer.close()

if __name__ == '__main__':
	main_class = CalculateMendelErrorRateGivenPlinkOutput
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
