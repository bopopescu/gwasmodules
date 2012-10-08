#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s --inputHeaderLs StKitts,Nevis  -o /tmp/Contig315_StKitts_vs_Nevis.tsv /tmp/Contig316_country144.vcf /tmp/Contig316_country148.vcf

Description:
	2012.10.6 program that juxtaposes alternative allele frequencies of same loci from all input VCF files.
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

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.VCFFile import VCFFile
from pymodule import SNP
from AbstractVCFMapper import AbstractVCFMapper

class JuxtaposeAlleleFrequencyFromMultiVCFInput(AbstractVCFMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVCFMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 1, ))
	option_default_dict.update({
						('inputHeaderLs', 1, ): ['AF1,AF2', 'H', 1, 'header of the output column to represent 1st input.', ],\
						})

	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractVCFMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		
		listArgumentName_data_type_ls = [("inputHeaderLs", str)]
		listArgumentName2hasContent = ProcessOptions.processListArguments(listArgumentName_data_type_ls,\
												emptyContent=[], class_to_have_attr=self)
	
	def _juxtaposeAlleleFrequencyFromMultiVCFInput(self, inputFnameLs=None, inputHeaderLs=None, outputFname=None, \
										defaultNullFrequency=-0, **keywords):
		"""
		2012.10.5
		
		"""
		sys.stderr.write("Getting allele frequency from %s input ..."%(len(inputFnameLs)))
		
		#get locus2AF from inputFname
		locus2frequencyList = []
		
		locus_id_set = set()
		for inputFname in inputFnameLs:
			vcfFile = VCFFile(inputFname=inputFname)
			locus2frequency = vcfFile.getLocus2AlternativeAlleleFrequency()
			vcfFile.close()
			locus2frequencyList.append(locus2frequency)
			locus_id_set = locus_id_set.union(set(locus2frequency.keys()))
		sys.stderr.write("%s loci.\n"%(len(locus_id_set)))
		
		sys.stderr.write("Outputting frequency collected from all input to %s ..."%(outputFname))
		#output them in juxtaposition
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		header = ['locusID'] + inputHeaderLs + ['count']
		writer.writerow(header)
		
		locus_id_list = list(locus_id_set)
		locus_id_list.sort()
		
		for locus_id in locus_id_list:
			locus_id_str_ls = map(str, locus_id)
			data_row = ['_'.join(locus_id_str_ls)]
			for i in xrange(len(locus2frequencyList)):
				locus2frequency = locus2frequencyList[i]
				frequency = locus2frequency.get(locus_id, defaultNullFrequency)
				data_row.append(frequency)
			data_row.append(1)
			writer.writerow(data_row)
		del writer
		sys.stderr.write("\n")
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		self._juxtaposeAlleleFrequencyFromMultiVCFInput(inputFnameLs=self.inputFnameLs, inputHeaderLs=self.inputHeaderLs, \
													outputFname=self.outputFname)
		
		

if __name__ == '__main__':
	main_class = JuxtaposeAlleleFrequencyFromMultiVCFInput
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()