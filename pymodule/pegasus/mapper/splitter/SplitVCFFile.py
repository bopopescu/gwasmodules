#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i ~namtran/panasas/Experiment/RNA-seq/Freimer/Developmental/ASE/Variant/MultipleSampleCalling/genome.algn.split.part17/samtools.var.filt.vcf.gz 
		-o /tmp/ -m 2
	
	#run in a for loop. "ls -d" only lists the folders and doesn't list contents of those folders.
	for i in `ls -d ~namtran/panasas/Experiment/RNA-seq/Freimer/Developmental/ASE/Variant/MultipleSampleCalling/genome*`; 
		do echo $i; 
		%s -i  $i/samtools.var.filt.vcf.gz  -o VariantsOf36RNA-SeqMonkeysFromNam_minDepth5/ -m 5;
	done
	
Description:
	2012.8.25
		This program splits a VCF file into multiple small ones.
	
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.yhio.VCFFile import VCFFile
from pymodule.pegasus.mapper.AbstractVCFMapper import AbstractVCFMapper

class SplitVCFFile(AbstractVCFMapper):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractVCFMapper.option_default_dict)
	option_default_dict.pop(("chromosome", 0, ))
	option_default_dict.pop(("chrLength", 1, int))
	option_default_dict.pop(('outputFname', 0, ))
	#option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('noOfOverlappingSites', 0, int): [1000, 'o', 1, 'The number of overlapping sites between two adjacent split output.', ],\
						('noOfSitesPerUnit', 1, int): [5000, 's', 1, 'The number of sites in each split output', ],\
						('noOfTotalSites', 1, int): [None, 'n', 1, 'The total number of sites in input VCF', ],\
						})
	def __init__(self,  **keywords):
		"""
		"""
		AbstractVCFMapper.__init__(self, **keywords)
		#turn them into nonnegative
		self.noOfOverlappingSites = abs(self.noOfOverlappingSites)
		self.noOfSitesPerUnit = abs(self.noOfSitesPerUnit)
		self.noOfTotalSites = abs(self.noOfTotalSites)
		
		if self.noOfSitesPerUnit > self.noOfTotalSites:
			sys.stderr.write("Error: noOfTotalSites, %s, could not be smaller than noOfSitesPerUnit %s.\n"%\
							(self.noOfTotalSites, self.noOfSitesPerUnit))
			sys.exit(3)
		if self.noOfOverlappingSites>self.noOfSitesPerUnit:
			sys.stderr.write("Error: noOfSitesPerUnit, %s, could not be smaller than noOfOverlappingSites %s.\n"%\
							(self.noOfSitesPerUnit, self.noOfOverlappingSites))
			sys.exit(2)
	
	def splitVCF(self, inputFname, outputFnamePrefix=None, noOfOverlappingSites=1000, noOfSitesPerUnit=5000,\
				noOfTotalSites=None):
		"""
		2012.8.25
			
		"""
		sys.stderr.write("Splitting VCF %s into files each with %s sites and %s overlapping ... \n"%(inputFname, noOfSitesPerUnit,\
																		noOfOverlappingSites))
		
		vcfFile = VCFFile(inputFname=inputFname)
		
		unitNumber2OutVCFFile = {}
		counter = 0
		real_counter = 0
		#make it 1 less than total so the last unit is >=s
		noOfUnits = max(1, utils.getNoOfUnitsNeededToCoverN(N=noOfTotalSites, s=noOfSitesPerUnit, o=noOfOverlappingSites)-1)
		sys.stderr.write(" will be split into %s units ... "%(noOfUnits))
		overlappingRecordLs = []
		for vcfRecord in vcfFile:
			counter += 1
			#below the maximum: noOfUnits.
			unitNumber = min(noOfUnits, max(1, utils.getNoOfUnitsNeededToCoverN(N=counter, s=noOfSitesPerUnit, o=noOfOverlappingSites)))
			if unitNumber not in unitNumber2OutVCFFile:
				outputFname = '%s_unit%s.vcf'%(outputFnamePrefix, unitNumber)
				outVCFFile = VCFFile(outputFname=outputFname)
				outVCFFile.metaInfoLs = vcfFile.metaInfoLs
				outVCFFile.header = vcfFile.header
				outVCFFile.writeMetaAndHeader()
				outVCFFile.noOfLoci =0
				#output the overlapping vcf records (from previous unit
				if overlappingRecordLs:
					for overlappingVCFRecord in overlappingRecordLs:
						outVCFFile.writeVCFRecord(overlappingVCFRecord)
						outVCFFile.noOfLoci += 1
					overlappingRecordLs = []	#reset it
				unitNumber2OutVCFFile[unitNumber] = outVCFFile
			outVCFFile = unitNumber2OutVCFFile[unitNumber]
			outVCFFile.writeVCFRecord(vcfRecord)
			outVCFFile.noOfLoci += 1
			#store the overlapping records
			if unitNumber<noOfUnits:
				if outVCFFile.noOfLoci>(noOfSitesPerUnit-noOfOverlappingSites):
					overlappingRecordLs.append(vcfRecord)
			
		
		vcfFile.close()
		#close all output files
		for unitNumber, outVCFFile in unitNumber2OutVCFFile.iteritems():
			outVCFFile.close()
		
		sys.stderr.write("%s loci split into %s files.\n"%(counter, len(unitNumber2OutVCFFile)))
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		self.splitVCF(inputFname=self.inputFname, outputFnamePrefix=self.outputFnamePrefix, \
					noOfOverlappingSites=self.noOfOverlappingSites, noOfSitesPerUnit=self.noOfSitesPerUnit,\
					noOfTotalSites=self.noOfTotalSites)

if __name__ == '__main__':
	main_class = SplitVCFFile
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()