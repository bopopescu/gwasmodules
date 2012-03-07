#!/usr/bin/env python
"""
Examples:
	%s  -i /Network/Data/250k/db/dataset/call_method_32.tsv 
		-j /Network/Data/250k/db/dataset/call_method_80.tsv -o /tmp/call_method_80_in_32_order.tsv  -m1
	
	%s
	
Description:
	2012.3.2
		Program to order the rows of the 2nd dataset (-j) in the order of the 1st datset (-i).
		Rows that are not present in 1st dataset will be discarded.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule import SNPData, TwoSNPData


class Order2ndSNPDataRowsSameAs1stSNPData(object):
	__doc__ = __doc__
	option_default_dict = {('input_fname1',1, ): [None, 'i', 1, 'to form SNPData1'],\
						('input_fname2',1, ): [None, 'j', 1, 'to form SNPData2'],\
						('output_fname', 1, ): [None, 'o', 1, 'Final output.', ],\
						('row_matching_by_which_value', 0, int):[0, 'm', 1, 'which column in the input_fname1 should be used to establish row-id linking to input_fname2. \
							0=both inputs discard the 2nd column and use the 1st column;\
							1=1st column(input_fname1 keeps both columns);\
							2=2nd column(ditto).'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']
							}
							#('input_fname1_format',1,int): [1, 'k', 1, 'Format of input_fname1. 1=strain X snp (Yu). 2=snp X strain (Bjarni) without arrayId. 3=snp X strain with arrayId.'],\
							#('input_fname2_format',1,int): [1, 'l', 1, 'Format of input_fname2. 1=strain X snp (Yu). 2=snp X strain (Bjarni) without arrayId'],\
	def __init__(self, **keywords):
		"""
		2008-09-17
			allow priority=3 or 4
		2008-06-02
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def run(self):
		"""
		2008-06-02
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		if self.row_matching_by_which_value==0:
			snpData1 = SNPData(input_fname=self.input_fname1, turn_into_array=1, ignore_2nd_column=1)
		else:
			snpData1 = SNPData(input_fname=self.input_fname1, turn_into_array=1)
		snpData2 = SNPData(input_fname=self.input_fname2, turn_into_array=1)
		
		if self.row_matching_by_which_value==1 or self.row_matching_by_which_value==2:
			row_matching_by_which_value = self.row_matching_by_which_value-1
		else:
			row_matching_by_which_value = None
		twoSNPData = TwoSNPData(SNPData1=snpData1, SNPData2=snpData2, debug=self.debug, row_matching_by_which_value=row_matching_by_which_value)
		newSnpData= twoSNPData.order2ndSNPDataRowsSameAs1stSNPData()
		newSnpData.tofile(self.output_fname)

if __name__ == '__main__':
	#do simple intersectSNPData1_and_SNPData2_row_wise()
	main_class = Order2ndSNPDataRowsSameAs1stSNPData
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()