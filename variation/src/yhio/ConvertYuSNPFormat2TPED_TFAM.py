#!/usr/bin/env python
"""

Examples:
	ConvertYuSNPFormat2TPED_TFAM.py -i /Network/Data/250k/tmp-yh/call_method_17.tsv -o /tmp/call_method_17_eigenstrat
	
	# get the TFAM & EMMAX-bin phenotype files (-y1 = take phenotype id 1) as well.
	# The StrainXSNP matrix is a watered-down version of original to speed up.
	# -d = recode12 for Hyun Min Kang's EMMAX binary.
	ConvertYuSNPFormat2TPED_TFAM.py -i /Network/Data/250k/tmp-yh/250k_data/call_method_17_test.tsv -o /tmp/call_method_17 -p /Network/Data/250k/tmp-yh/phenotype/phenotype.tsv -y 1 -d
	# after above run emmax
	./emmax-kin -v -h -s -d 10 /tmp/call_method_17
	./emmax -v -d 10 -t /tmp/call_method_17 -p /tmp/call_method_17_1_LD.pheno -k /tmp/call_method_17.hIBS.kinf -o /tmp/call_method_17_EMMAX
	
Description:
	Convert Yu's SNP format (input) to TPED, TFAM, EMMAX's phenotype format.
	
	Input format is strain X snp. 2nd column ignored by default (change by array_id_2nd_column).
		delimiter (either tab or comma) is automatically detected. 	header is 'chromosome_position'.
	
	tped (transposed PED fortmat), tfam are used in Plink. http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped
	
	There are 3 output files.
		<--------- xxxx.tped (genotype) -------->      <- xxxx.tfam (phenotype) ->
	     1 snp1 0 5000650 A A A C C C A C C C C C      1  1  0  0  1  1
	     1 snp2 0 5000830 G T G T G G T T G T T T      2  1  0  0  1  1
	                                                   3  1  0  0  1  1
	                                                   4  1  0  0  1  2
	                                                   5  1  0  0  1  2
	                                                   6  1  0  0  1  2
	<---- xxxx.pheno (phenotype for EMMAX) ---->
	1  1  1
	2  1  1
	3  1  1
	4  1  2
	5  1  2
	6  1  2
	
	First 4 columns of xxxx.tped are:
		chromosome (1-22, X, Y or 0 if unplaced)
		rs# or snp identifier
		Genetic distance (morgans)
		Base-pair position (bp units)
		
		After which the two genotypes for each individual are in the order of xxxx.tfam. Two columns (diploid) for one individual.
	
	First 6 columns of xxxx.tfam are:
		Family ID
		Individual ID
		Paternal ID
		Maternal ID
		Sex (1=male; 2=female; other=unknown)
		Phenotype
	
	3 columns of xxxx.pheno are in the same order of xxxx.tfam.
		Family ID
		Individual ID
		Phenotype
"""
import sys, os, math
import csv, numpy
bit_number = math.log(sys.maxint)/math.log(2)
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script/')))
from pymodule.yhio.SNP import transposeSNPData, number2nt
from variation.src.association.Association import Association

class ConvertYuSNPFormat2TPED_TFAM(object):
	__doc__ = __doc__
	option_default_dict = {('input_fname', 1, ): [None, 'i', 1, ],\
						('output_fname_prefix', 1, ): [None, 'o', 1, 'Output Filename prefix'],\
						('array_id_2nd_column', 0, int): [0, 'a', 0, 'whether 2nd column in input_fname is array id or not'],\
						('phenotype_fname', 0, ): [None, 'p', 1, 'phenotype file, same number and order of individuals as input_fname', ],\
						('phenotype_method_id', 0, int): [0, 'y', 1, 'which phenotype',],\
						('phenotype_is_binary', 0, int): [0, 'e', 1, 'is phenotype binary?',],\
						('recode12', 0, int):[0, 'd', 0, 'argument from plink. recode the alleles as 1 and 2 (and the missing genotype is 0). '],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-12-02
			modelled after ConvertYuSNPFormat2Bjarni.py
		"""
		from pymodule import ProcessOptions
		self.ad=ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def run(self):
		"""
		2008-12-02
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		initData = Association.readInData(self.phenotype_fname, self.input_fname, phenotype_method_id_ls=[self.phenotype_method_id],\
								report=self.report, snpAlleleOrdinalConversion=False)	# no ordinal (binary) conversion
		
		if self.recode12:
			# 2010-2-28 encode allele as 1,2,... (0 is missing).
			snpData, allele_index2allele_ls= initData.snpData.convert2Binary(in_major_minor_order=False, NA_notation=0, alleleStartingIndex=1)
		else:
			snpData = initData.snpData
		no_of_rows = len(snpData.data_matrix)
		
		tfam_w = csv.writer(open('%s.tfam'%(self.output_fname_prefix), 'w'), delimiter=' ')
		if initData.phenData:
			phenData = initData.phenData
			phenotype_col_index = initData.which_phenotype_ls[0]
			phenotype_label = phenData.col_id_ls[phenotype_col_index]
			phenotype_w = csv.writer(open('%s_%s.pheno'%(self.output_fname_prefix, phenotype_label.replace('/', '_')), 'w'), delimiter=' ')
			#no_of_rows = len(phenData.data_matrix) 
		for i in range(no_of_rows):
			if initData.phenData:
				phenotype_value = phenData.data_matrix[i,phenotype_col_index]
			else:
				phenotype_value = 1
			ecotype_id = initData.snpData.row_id_ls[i][0]
			if self.array_id_2nd_column:
				array_id = initData.snpData.row_id_ls[i][1]
			else:
				array_id = ecotype_id
			tfam_w.writerow([1, ecotype_id, 0, 0, 1, phenotype_value])	#1 to make everyone in the same family
			if numpy.isnan(phenotype_value):
				phenotype_value = 'NA'
			if initData.phenData:
				phenotype_w.writerow([ecotype_id, array_id, phenotype_value])
		del tfam_w
		
		genotype_w = csv.writer(open('%s.tped'%self.output_fname_prefix, 'w'), delimiter=' ')
		
		#transpose it
		newSnpData = transposeSNPData(snpData)
		no_of_rows = len(newSnpData.data_matrix)
		no_of_cols = len(newSnpData.data_matrix[0])
		for i in range(no_of_rows):
			snp_id = newSnpData.row_id_ls[i]
			chromosome, pos = snp_id.split('_')[:2]
			row = [chromosome, "%s_%s"%(chromosome, pos), 0, pos,]
			for j in range(no_of_cols):
				allele = newSnpData.data_matrix[i][j]
				if not self.recode12:	# it's in my nucleotide number format.
					allele = number2nt[allele]
				if allele !='NA':
					row.append(allele[0])
					if len(allele)>1:	#2012.8.20 diploid, 2nd allele
						row.append(allele[1])
					else:
						row.append(allele[0])
				else:
					row.append("N")	# "N" is plink for missing.
					row.append("N")
			genotype_w.writerow(row)
		del genotype_w

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = ConvertYuSNPFormat2TPED_TFAM
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
