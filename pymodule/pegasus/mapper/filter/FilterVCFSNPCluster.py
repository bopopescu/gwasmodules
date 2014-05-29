#!/usr/bin/env python
"""
Examples:
	%s -i ~/NetworkData/vervet/db/genotype_file/method_27/*Contig0.vcf.gz -o /tmp/Contig0.filter.vcf.gz -m 10
	
	%s 
	
	%s 
	
Description:
	2012.9.6
		program that filters out SNPs that are too close to each other.
		If two SNPs are within minNeighborDistance of each other, both will be removed.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import cStringIO, re, csv
from pymodule import ProcessOptions, figureOutDelimiter
from pymodule.utils import sortCMPBySecondTupleValue
from pymodule.yhio.VCFFile import VCFFile
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper

class FilterVCFSNPCluster(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.update({
						('minNeighborDistance', 1, int): [10, 'm', 1, 'minimum distance between two adjacent SNPs'],\
						}
						)

	def __init__(self,  inputFnameLs=None, **keywords):
		"""
		2011-7-12
		"""
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	
	def filterVCFSNPCluster(self, inputFname=None, outputFname=None, minNeighborDistance=10, **keywords):
		"""
		#2012.8.20 locus_id2row_index from VCFFile is using (chr, pos) as key, not chr_pos
			need a conversion in between
		2012.5.8
		"""
		sys.stderr.write("Filtering VCF %s to get rid of SNPs that are %s distance apart ..."%(inputFname, minNeighborDistance))
		vcfFile = VCFFile(inputFname=inputFname)
		
		outVCFFile = VCFFile(outputFname=outputFname)
		outVCFFile.metaInfoLs = vcfFile.metaInfoLs
		outVCFFile.header = vcfFile.header
		outVCFFile.writeMetaAndHeader()
		
		previousVCFRecord = None
		previousVCFRecordIsBad = False #indicator whether previous record is bad or not. based on distance to the previous-previous record
		counter = 0
		for vcfRecord in vcfFile:
			if previousVCFRecord is not None:
				if previousVCFRecord.chr == vcfRecord.chr:
					distanceToPreviousRecord = abs(vcfRecord.pos - previousVCFRecord.pos)
					if distanceToPreviousRecord <minNeighborDistance:
						previousVCFRecordIsBad = True
					else:
						if not previousVCFRecordIsBad:	#distance to current & previous-previous record is >=minNeighborDistance
							outVCFFile.writeVCFRecord(previousVCFRecord)
						previousVCFRecordIsBad = False
				else:
					#handle the last record from the previous chromosome (assuming loci are in chromosomal order)
					if not previousVCFRecordIsBad:	#distance to previous-previous record is >=minNeighborDistance
						outVCFFile.writeVCFRecord(previousVCFRecord)
					
					previousVCFRecordIsBad = False#reset
					
			previousVCFRecord = vcfRecord
			counter += 1
		vcfFile.close()
		
		#handle the last record
		if previousVCFRecord is not None and not previousVCFRecordIsBad:	#distance to previous-previous record is >=minNeighborDistance
			outVCFFile.writeVCFRecord(previousVCFRecord)
		outVCFFile.close()
		
		noOfLociAfterFilter = len(outVCFFile.locus_id_ls)
		delta = counter-noOfLociAfterFilter
		if counter>0:
			fraction = delta/float(counter)
		else:
			fraction = -0.0
		sys.stderr.write(" %s (%s -> %s) or %.2f%% loci filtered out.\n"%(delta, counter, noOfLociAfterFilter, fraction*100))
		
	
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
		
		self.filterVCFSNPCluster(inputFname=self.inputFname, outputFname=self.outputFname, minNeighborDistance=self.minNeighborDistance)
		

if __name__ == '__main__':
	main_class = FilterVCFSNPCluster
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()