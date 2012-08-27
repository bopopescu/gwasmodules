#!/usr/bin/env python
"""
Examples:
	%s -i /tmp/vs_top150Contigs_by_aln_Contig0.vcf -o /tmp/vs_top150Contigs_by_aln_Contig0.tsv
		-u ,
	
	%s 
	
	%s 
	
Description:
	2011-9-1
		Two functions:
		1. A multi-sample genotype caller based entirely on coverage of reads.
			sam/bam file has to be indexed beforehand.
		2. A coverage-based GATK-generated VCF file filter.
		
		For multi-read-group input, seqCoverageFname will provide coverage for each individual sequence.
			Each read group has the individual_sequence.id embedded as 2nd value if it's split by "_".
			Read group formula: alnID_isqID_individual.code_sequencer_vs_alignment.ref_ind_seq_id
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import cStringIO, re, csv
from pymodule import ProcessOptions, figureOutDelimiter
from pymodule.utils import sortCMPBySecondTupleValue
from pymodule.VCFFile import VCFFile

class ConvertVCF2BjarniFormat(object):
	__doc__ = __doc__
	option_default_dict = {('inputFname', 1, ): ['', 'i', 1, 'The input Bam file.', ],\
						('outputFname', 1, ): [None, 'o', 1, 'output the SNP data.'],\
						('outputDelimiter', 1, ): ['\t', 'u', 1, 'delimiter in the csv output file.'],\
						("run_type", 1, int): [1, 'y', 1, ''],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self,  **keywords):
		"""
		2011-7-12
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
	
	
	def outputCallMatrix(self, data_matrix=None, refFastaFname=None, outputFname=None, refNameSet=None, read_group2col_index=None, \
						locus_id2row_index=None, outputDelimiter='\t'):
		"""
		2012.8.20
			add argument outputDelimiter
		"""
		
		# output the matrix in the end
		#read_group2col_index.pop('ref', None)	#remove the "ref" item if "ref" is in read_group2col_index. None is for failsafe when "ref" is not present.
		read_group_col_index_ls = read_group2col_index.items()
		read_group_col_index_ls.sort(cmp=sortCMPBySecondTupleValue)
		header = ['Chromosome', 'Positions']+[row[0] for row in read_group_col_index_ls]
		writer = csv.writer(open(outputFname, 'w'), delimiter=outputDelimiter)
		writer.writerow(header)
		
		locus_id_and_row_index_ls = locus_id2row_index.items()
		locus_id_and_row_index_ls.sort(cmp=sortCMPBySecondTupleValue)
		for i in xrange(len(locus_id_and_row_index_ls)):
			locus_id, row_index = locus_id_and_row_index_ls[i]
			data_row = data_matrix[i]
			refName, pos = locus_id[:2]
			# if data_row is shorter than read_group_col_index_ls, add "NA" to fill it up
			for j in xrange(len(data_row), len(read_group_col_index_ls)):
				data_row.append('NA')
			writer.writerow([refName, pos] + data_row)
		del writer
	
	def convertVCF2BjarniFormat(self, inputFname, outputFname, **keywords):
		"""
		#2012.8.20 locus_id2row_index from VCFFile is using (chr, pos) as key, not chr_pos
			need a conversion in between
		2012.5.8
		"""
		vcfFile = VCFFile(inputFname=inputFname)
		vcfFile.parseFile()
		
		read_group2col_index = vcfFile.sample_id2index
		locus_id2row_index = vcfFile.locus_id2row_index	
		
		data_matrix = vcfFile.genotype_call_matrix
		
		self.outputCallMatrix(data_matrix, refFastaFname=None, outputFname=outputFname, refNameSet=None, \
					read_group2col_index=read_group2col_index, \
					locus_id2row_index=locus_id2row_index, outputDelimiter=self.outputDelimiter)
	
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
		
		self.convertVCF2BjarniFormat(self.inputFname, self.outputFname, \
					outputDelimiter=self.outputDelimiter,\
					report=self.report)
		
		

if __name__ == '__main__':
	main_class = ConvertVCF2BjarniFormat
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()