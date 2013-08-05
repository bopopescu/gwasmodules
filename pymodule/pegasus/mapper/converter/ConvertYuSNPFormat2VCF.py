#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i /Network/Data/250k/db/dataset/call_method_32.tsv -o /tmp/call_32.vcf

Description:
	2013.03

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions
from pymodule import SNPData, VCFFile
from pymodule import utils
from pymodule.yhio.SNP import number2di_nt
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper

class ConvertYuSNPFormat2VCF(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.update({
							('min_MAF', 0, float): [None, 'n', 1, 'minimum minor allele frequency', ],\
							})
	def __init__(self,  **keywords):
		"""
		"""
		AbstractMapper.__init__(self, **keywords)
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		snpData = SNPData(input_fname=self.inputFname, turn_into_array=1, ignore_2nd_column=1)
		snpData = SNPData.removeMonomorphicCols(snpData, NA_set=set([]))
		if self.min_MAF and self.min_MAF>0:
			snpData = SNPData.removeColsByMAF(snpData,min_MAF=self.min_MAF, NA_set=set([]))
		
		self.writer = VCFFile(outputFname=self.outputFname, openMode='w')
		self.writer.makeupHeaderFromSampleIDList(sampleIDList=snpData.row_id_ls)
		self.writer.writeMetaAndHeader()
		
		counter = 0
		for j in xrange(len(snpData.col_id_ls)):
			snp_id = snpData.col_id_ls[j]
			chromosome, start = snp_id.split('_')[:2]
			genotype_ls = snpData.data_matrix[:,j]
			genotype_ls = utils.dict_map(number2di_nt, genotype_ls)
			genotype_ls_vcf = []
			alleleNucleotide2Number = {}
			alleleNumber2Nucleotide = {}
			for genotype in genotype_ls:
				if genotype=='NA':
					genotype_ls_vcf.append("./.")
				elif len(genotype)==2:
					for allele in genotype:
						if allele not in alleleNucleotide2Number:
							alleleNumber = len(alleleNucleotide2Number)
							alleleNucleotide2Number[allele] = alleleNumber
							alleleNumber2Nucleotide[alleleNumber] = allele
					genotype_ls_vcf.append("%s/%s"%(alleleNucleotide2Number[genotype[0]], alleleNucleotide2Number[genotype[1]]))
					
				else:
					genotype_ls_vcf.append("./.")
			refAllele = alleleNumber2Nucleotide[0]
			if 1 not in alleleNumber2Nucleotide:
				altAllele = refAllele
			else:
				altAllele = alleleNumber2Nucleotide[1]
			row=[chromosome, start, ".", refAllele, altAllele, 999, 'PASS', "DP=100", "GT"] + genotype_ls_vcf
			self.writer.writerow(row)
			counter += 1
		sys.stderr.write("  %s records.\n"%(counter))
		self.writer.close()
		

if __name__ == '__main__':
	main_class = ConvertYuSNPFormat2VCF
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()