#!/usr/bin/env python
"""
Examples:
	%s -i folderRun/Scaffold301_splitVCF_unit1.vcf -o folderRun/Scaffold301_splitVCF_unit1.info.vcf
	
	%s 
	
	%s 
	
Description:
	2013.07.10 program that adds description of info fields/tags that appear in sites but not in VCF header, into VCF header.
		At this moment, it just addes description of these fields into VCF header, and does NOT check missing info tags.
	LDAF
	ERATE
	AVGPOST
	AF_Orig
	AC_Orig
	AN_Orig

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import cStringIO, re, csv
from pymodule import ProcessOptions, figureOutDelimiter
from pymodule.utils import sortCMPBySecondTupleValue
from pymodule.yhio.VCFFile import VCFFile
from pymodule.pegasus.mapper.AbstractVCFMapper import AbstractVCFMapper

parentClass = AbstractVCFMapper
class AddMissingInfoDescriptionToVCFHeader(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update({
						}
						)
	knownInfoTag2DescriptionLine = {"LDAF": """##INFO=<ID=LDAF,Number=1,Type=Float,Description="MLE Allele Frequency Accounting for LD. Range: 0 - 1">\n""",\
				"ERATE": """##INFO=<ID=ERATE,Number=1,Type=Float,Description="Per-marker Mutation rate from MaCH/Thunder. Range: 0.0001 - 0.2051">\n""",\
				"AVGPOST": """##INFO=<ID=AVGPOST,Number=1,Type=Float,Description="Average posterior probability from MaCH/Thunder. Range: 0.5242 - 1">\n""",\
				"RSQ": """##INFO=<ID=RSQ,Number=1,Type=Float,Description="Genotype imputation quality from MaCH/Thunder. Range:0 - 1">\n""",\
				"THETA": """##INFO=<ID=THETA,Number=1,Type=Float,Description="Per-marker Transition rate from MaCH/Thunder. Range:0 - 0.1493">\n""",\
				"AC_Orig": """##INFO=<ID=AC_Orig,Number=1,Type=Integer,Description="Original AC">\n""",\
				"AF_Orig": """##INFO=<ID=AF_Orig,Number=1,Type=Float,Description="Original AF">\n""",\
				"AN_Orig": """##INFO=<ID=AN_Orig,Number=1,Type=Integer,Description="Original AN">\n""",\
				}
	
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		parentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	
	def getAllInfoTags(self, inputFname=None, **keywords):
		"""
		2013.07.10
			not used right now.
		"""
		sys.stderr.write("Extracting info tags from  VCF %s ..."%(inputFname))
		vcfFile = VCFFile(inputFname=inputFname)
		
		info_tag_set = set()
		counter = 0
		real_counter = 0
		for vcfRecord in vcfFile:
			for info_tag in vcfRecord.info_tag2value:
				info_tag_set.add(info_tag)
			counter += 1
		vcfFile.close()
		
		sys.stderr.write("%s unique info tags.\n"%(len(info_tag_set)))
		return info_tag_set
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
			debug = True
		else:
			debug =False
		
		
		
		outputDir = os.path.split(self.outputFname)[0]
		if outputDir and not os.path.isdir(outputDir):
			os.makedirs(outputDir)
		
		
		self.reader = VCFFile(inputFname=self.inputFname)
		
		self.writer = VCFFile(outputFname=self.outputFname, openMode='w')
		self.writer.metaInfoLs = self.reader.metaInfoLs
		for info_tag, description in self.knownInfoTag2DescriptionLine.iteritems():
			self.writer.metaInfoLs.append(description)
		self.writer.header = self.reader.header
		self.writer.writeMetaAndHeader()
		
		counter = 0
		for vcfRecord in self.reader:
			counter += 1
			self.writer.writeVCFRecord(vcfRecord)
		
		self.reader.close()
		self.writer.close()
		

if __name__ == '__main__':
	main_class = AddMissingInfoDescriptionToVCFHeader
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()