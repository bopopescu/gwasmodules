#!/usr/bin/env python
"""
Examples:
	%s -i ~/NetworkData/vervet/db/genotype_file/method_27/*Contig0.vcf.gz -o /tmp/Contig0_HaplotypeScore.tsv -k HaplotypeScore
	
	%s 
	
	%s 
	
Description:
	2012.9.17 program that extracts certain key in the INFO column of VCF and outputs them in tsv format.
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

class ExtractInfoFromVCF(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.update({
						('infoKey', 1, ): ['HaplotypeScore', 'k', 1, 'the key of the INFO to be extracted'],\
						}
						)

	def __init__(self,  **keywords):
		"""
		2011-7-12
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
	
	
	def extract(self, inputFname=None, outputFname=None, infoKey='HaplotypeScore', **keywords):
		"""
		2012.9.17
		"""
		sys.stderr.write("Extracting %s from  VCF %s ..."%(infoKey, inputFname))
		vcfFile = VCFFile(inputFname=inputFname)
		
		writer= csv.writer(open(outputFname, 'w'), delimiter='\t')
		header = ['CHROM', 'POS', infoKey]
		writer.writerow(header)
		
		counter = 0
		real_counter = 0
		for vcfRecord in vcfFile:
			value = vcfRecord.info_tag2value.get(infoKey)
			if value is not None:
				data_row = [vcfRecord.chr, vcfRecord.pos, value]
				writer.writerow(data_row)
				real_counter += 1
			counter += 1
		vcfFile.close()
		
		del writer
		if counter>0:
			fraction = real_counter/float(counter)
		else:
			fraction = -0.0
		sys.stderr.write("%.2f%% (%s/%s) loci have %s.\n"%( fraction*100, real_counter, counter, infoKey))
	
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
		
		self.extract(inputFname=self.inputFname, outputFname=self.outputFname, infoKey=self.infoKey)
		

if __name__ == '__main__':
	main_class = ExtractInfoFromVCF
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()