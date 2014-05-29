#!/usr/bin/env python
"""
Examples:
	%s -i tmp/testVCF/CAEY.sorted.vcf.gz  -o /tmp/VCF_CAEY.unique.vcf
	
	%s 
	
	%s 
	
Description:
	2013.07.12 program that removes any redundant SNPs. By redundant, (chromosome, position) appears >1.

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from Bio.Seq import Seq
from pymodule import ProcessOptions, MatrixFile, PassingData
from pymodule.yhio.VCFFile import VCFFile
from pymodule.pegasus.mapper.AbstractVCFMapper import AbstractVCFMapper
from pymodule import SNP
from pymodule.yhio.SNP import nt2number
from CalculateSameSiteConcordanceInVCF import CalculateSameSiteConcordanceInVCF

parentClass = CalculateSameSiteConcordanceInVCF
class RemoveRedundantLociFromVCF(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update({
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		parentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		
		outputDir = os.path.split(self.outputFname)[0]
		if outputDir and not os.path.isdir(outputDir):
			os.makedirs(outputDir)
		
		snp_pos2count =self.readInSNPID2GenotypeVectorLs(self.inputFname, returnType=2).snp_pos2returnData
		
		
		reader = VCFFile(inputFname=self.inputFname)
		
		writer = VCFFile(outputFname=self.outputFname, openMode='w')
		writer.metaInfoLs = reader.metaInfoLs
		writer.header = reader.header
		writer.writeMetaAndHeader()
		
		
		counter = 0
		real_counter = 0
		for vcfRecord in reader:	#assuming input VCF is sorted
			counter += 1
			key = (vcfRecord.chromosome, vcfRecord.position)
			frequency = snp_pos2count.get(key)
			if frequency==1:
				writer.writeVCFRecord(vcfRecord)
				real_counter += 1
			
		reader.close()
		writer.close()
		if counter >0:
			fraction=real_counter/float(counter)
		else:
			fraction = 0
		sys.stderr.write("%s (out of %s, %s) snps are unique.\n"%(real_counter, counter, fraction))

if __name__ == '__main__':
	main_class = RemoveRedundantLociFromVCF
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()