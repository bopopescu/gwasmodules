#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s  -i 59444 -u yh -z banyan -o /tmp/peak_59444.h5

Description:
	2012.3.7
		program to output id of loci within numerous result peaks.

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule import SNPData, SNP
#from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
from AbstractVariationMapper import AbstractVariationMapper
import h5py, numpy
from variation.src import Stock_250kDB

class OutputLociIDOfResultPeakInHDF5(AbstractVariationMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVariationMapper.option_default_dict.copy()
	option_default_dict.update({
							('peak_id_ls', 1, ): ['', 'i', 1, 'loci from these peaks will be outputted', ],\
							('datasetName', 1, ): ['locus_id_ls', '', 1, 'name of the dataset in both the HDF5 output', ],\
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		AbstractVariationMapper.__init__(self, inputFnameLs, **keywords)
		self.peak_id_ls = getListOutOfStr(self.peak_id_ls, data_type=int)
		self.peak_id_ls.sort()
	
	def run(self):
		"""
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_250k = self.db_250k
		
		pd = PassingData(min_MAF=self.min_MAF,\
							data_dir=self.data_dir, \
							need_chr_pos_ls=0,)
		for peak_id in self.peak_id_ls:
			result_peak = Stock_250kDB.ResultPeak.get(peak_id)
			
			rm = result_peak.result
			
			if rm.cnv_method_id and not db_250k._cnv_id2chr_pos:
				db_250k.cnv_id2chr_pos = rm.cnv_method_id
				pd.db_id2chr_pos = db_250k.cnv_id2chr_pos
			elif rm.call_method_id:
				pd.db_id2chr_pos = db_250k.snp_id2chr_pos
			
			gwr = db_250k.getResultMethodContent(rm.id, pdata=pd)
			gwr.data_obj_ls.sort(cmp=SNP.cmpDataObjByChrPos)
			
			locus_id_ls = []
			for data_obj in gwr.data_obj_ls:
				if data_obj.chromosome>result_peak.chromosome:
					break
				if data_obj.chromosome==result_peak.chromosome:
					if data_obj.position>result_peak.stop:
						break
					
					if data_obj.stopPosition>=result_peak.start and data_obj.position<=result_peak.stop:
						locus_id_ls.append(data_obj.db_id)
			sys.stderr.write("%s loci within peak %s.\n"%(len(locus_id_ls), peak_id))
		
		outputF = h5py.File(self.outputFname, 'w')
		shape = (len(locus_id_ls),)	#initial shape
		ds = outputF.create_dataset(self.datasetName, shape, numpy.int, compression='gzip', compression_opts=4, maxshape=(None,), chunks=True)
		for i in xrange(len(locus_id_ls)):
			ds[i] = locus_id_ls[i]
		outputF.close()
		

if __name__ == '__main__':
	main_class = OutputLociIDOfResultPeakInHDF5
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()