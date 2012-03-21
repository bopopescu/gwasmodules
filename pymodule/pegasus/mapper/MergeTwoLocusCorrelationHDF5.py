#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -d correlation  -o /tmp/output2 /tmp/out /tmp/out

Description:
	2012.3.1
	This program merges multiple output of CalculateColCorBetweenTwoHDF5/FindMaxLDBetweenPeakAndEachLocus into one.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule import SNPData
from AbstractMapper import AbstractMapper
import h5py, numpy

class MergeTwoLocusCorrelationHDF5(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 1, ))
	option_default_dict.update({
							('datasetName', 1, ): ['correlation', 'd', 1, 'name of the dataset in both the HDF5 input & output', ],\
							('datasetExpansionUnit', 1, int): [10000000, '', 1, 'How many rows to expand every time the HDF5 dataset is too small.', ],\
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		AbstractMapper.__init__(self, **keywords)
		self.inputFnameLs = inputFnameLs
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		outputF = h5py.File(self.outputFname, 'w')
		dtype = numpy.dtype([('input1LocusID', 'i'), ('input2LocusID', 'i'), ('correlation', 'f')])
		shape = (self.datasetExpansionUnit,)	#initial shape
		ds = outputF.create_dataset(self.datasetName, shape, dtype, compression='gzip', compression_opts=4, \
								maxshape=(None,), chunks=True)
		#compression seems to have little impact on the running speed.
		no_of_records = 0
		for inputFname in self.inputFnameLs:
			if not os.path.isfile(inputFname):
				continue
			f1 = h5py.File(inputFname, 'r')
			d1 = f1[self.datasetName]
			d1_length = d1.shape[0]
			no_of_records += d1_length
			if ds.shape[0]<no_of_records:
				ds.resize((no_of_records+self.datasetExpansionUnit,))
			#2012.3.21 d1 and d1[:] show different performance. the latter is a lot faster.
			ds[no_of_records-d1_length:no_of_records] = d1[:]
			f1.close() 
		
		if ds.shape[0]!=no_of_records:
			ds.resize((no_of_records,))
		outputF.close()
		

if __name__ == '__main__':
	main_class = MergeTwoLocusCorrelationHDF5
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()