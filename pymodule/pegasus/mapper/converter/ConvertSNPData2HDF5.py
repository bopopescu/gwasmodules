#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i /Network/Data/250k/db/dataset/call_method_32.tsv -o /tmp/call_32.hdf5

Description:
	2012.3.1
		assuming row and column ids are of type integer.

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
import h5py
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule import SNPData
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper

class ConvertSNPData2HDF5(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.update({
							('min_MAF', 1, float): [0.1, 'n', 1, 'minimum minor allele frequency', ],\
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
		if self.min_MAF>0:
			snpData = SNPData.removeColsByMAF(snpData,min_MAF=self.min_MAF, NA_set=set([]))
		snpData.col_id_ls = map(int, snpData.col_id_ls)
		snpData.row_id_ls = map(int, snpData.row_id_ls)
		f = h5py.File(self.outputFname, 'w')
		import numpy
		#snpData.data_matrix.dtype = numpy.int16
		dset = f.create_dataset("data_matrix", data=snpData.data_matrix, maxshape=(None, None))
		#numpy.array(snpData.data_matrix, dtype=numpy.int64)
		col_id_ls_dset = f.create_dataset('col_id_ls', data=snpData.col_id_ls, maxshape=(None,))
		row_id_ls_dset = f.create_dataset('row_id_ls', data=snpData.row_id_ls, maxshape=(None,))
		f.close()
		

if __name__ == '__main__':
	main_class = ConvertSNPData2HDF5
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()