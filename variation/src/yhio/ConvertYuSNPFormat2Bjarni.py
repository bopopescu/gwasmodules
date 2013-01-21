#!/usr/bin/env python
"""

Examples:
	./src/ConvertYuSNPFormat2Bjarni.py -i data/2010/data_2010_ecotype_id_y0002_n1c1d110_mergedup.tsv
		-o data/2010/data_2010_ecotype_id_y0002_n1c1d110_mergedup_bjarni.csv -r
	
	./src/ConvertYuSNPFormat2Bjarni.py -i genotyping/149snp/stock_149SNP_y0000110101_mergedup.csv
		-o genotyping/149snp/stock_149SNP_y0000110101_mergedup_bjarni.csv  -r
	
	./src/ConvertYuSNPFormat2Bjarni.py  -i genotyping/384-illumina_y0000110101_mergedup.csv
		-o genotyping/384-illumina_y0000110101_mergedup_bjarni.csv -r
	
	./src/ConvertYuSNPFormat2Bjarni.py -i /mnt/nfs/NPUTE_data/input/250k_l3_y0.85.tsv
		-o /mnt/nfs/NPUTE_data/input/250k_l3_y0.85_bjarni.tsv -a -b
		
	# 2011-2-28 convert call method 75 into Bjarni format with TAIR9 coordinates.
	# it doesn't matter if call_method doesn't exist in banyan db. It only needs the Snps table.
	ConvertYuSNPFormat2Bjarni.py -i /Network/Data/250k/db/dataset/call_method_75.tsv
		-o /Network/Data/250k/db/dataset/call_method_75_TAIR9.csv -a -z banyan

Description:
	Convert Yu's SNP format (input) to Bjarni's.

	Input format is strain X snp. 2nd column ignored by default (change by array_id_2nd_column).
		delimiter (either tab or comma) is automatically detected.
		Column/SNP ID is 'chromosome_position' or Snps.id.
			If it's the latter, program uses db to translate it into chr_pos.
		Different db connection might result in different chr_pos due to different TAIR versions.
	
	Output format is snp X strain, bjarni's format. The addition of an extra first row for array id depends
		on the array_id_2nd_column option. If it's toggled, yes.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/')))
import numpy
from pymodule import process_function_arguments, write_data_matrix, figureOutDelimiter, read_data, SNPData,\
	SNPData2RawSnpsData_ls
from variation.src import Stock_250kDB
from variation.src import AbstractVariationMapper
import snpsdata

class ConvertYuSNPFormat2Bjarni(AbstractVariationMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVariationMapper.option_default_dict.copy()
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('ecotype_table', 1, ): ['stock.ecotype', '', 1, 'ecotype Table to get ecotypeid2nativename'],\
						('array_id_2nd_column', 0, int): [0, 'a', 0, 'whether 2nd column in inputFname is array id or not'],\
						('snp_id_type', 1, int): [1, 'y', 1, 'type 1: db_id (need db translation); type 2: chr_pos'],\
						})
	def __init__(self, **keywords):
		"""
		2008-05-18
			add argument array_id_2nd_column
		2008-5-12
		"""
		AbstractVariationMapper.__init__(self, inputFnameLs=None, **keywords)
	
	def run(self):
		"""
		2008-5-12
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		#database connection and etc
		db = self.db_250k
		
		session = db.session
		session.begin()
		
		delimiter = figureOutDelimiter(self.inputFname, report=self.report)
		header, strain_acc_list, category_list, data_matrix = read_data(self.inputFname, delimiter=delimiter)
		
		if self.snp_id_type==1:
			#2011-2-27 translate the db_id into chr_pos because the new StrainXSNP dataset uses db_id to identify SNPs.
			# but if col-id is already chr_pos, it's fine.
			new_header = header[:2]
			data_matrix_col_index_to_be_kept = []
			for i in xrange(2, len(header)):
				snp_id = header[i]
				chr_pos = db.get_chr_pos_given_db_id2chr_pos(snp_id,)
				if chr_pos is not None:
					data_matrix_col_index_to_be_kept.append(i-2)
					new_header.append(chr_pos)
			# to remove no-db_id columns from data matrix
			data_matrix = numpy.array(data_matrix)
			data_matrix = data_matrix[:, data_matrix_col_index_to_be_kept]
			header = new_header
		
		if self.array_id_2nd_column:
			snpData = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list,\
							data_matrix=data_matrix)
		else:
			snpData = SNPData(header=header, strain_acc_list=strain_acc_list,\
							data_matrix=data_matrix)	#ignore category_list
		
		rawSnpsData_ls = SNPData2RawSnpsData_ls(snpData, need_transposeSNPData=1, report=self.report)
		chromosomes = [rawSnpsData.chromosome for rawSnpsData in rawSnpsData_ls]
		snpsdata.writeRawSnpsDatasToFile(self.outputFname, rawSnpsData_ls, chromosomes=chromosomes, deliminator=',', withArrayIds=self.array_id_2nd_column)
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = ConvertYuSNPFormat2Bjarni
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
