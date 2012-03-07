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
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
import h5py, numpy
from variation.src import Stock_250kDB

class OutputLociIDOfResultPeakInHDF5(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 1, ))
	option_default_dict.update({
							('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'stock_250k database name', ],\
							('schema', 0, ): ['', 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('port', 0, ):[None, '', 1, 'database port number'],\
							('results_directory', 0, ):[None, 't', 1, 'The results directory. Default is None. use the one given by db.'],\
							('min_MAF', 0, float): [0.1, 'n', 1, 'minimum Minor Allele Frequency.'],\
							
							('peak_id_ls', 1, ): ['', 'i', 1, 'loci from these peaks will be outputted', ],\
							('datasetName', 1, ): ['locus_id_ls', '', 1, 'name of the dataset in both the HDF5 output', ],\
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		AbstractMapper.__init__(self, **keywords)
		self.inputFnameLs = inputFnameLs
		self.peak_id_ls = getListOutOfStr(self.peak_id_ls, data_type=int)
		self.peak_id_ls.sort()
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_250k = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user, password=self.db_passwd, \
									hostname=self.hostname, database=self.dbname)
		db_250k.setup(create_tables=False)
		
		for peak_id in self.peak_id_ls:
			result_peak = Stock_250kDB.ResultPeak.get(peak_id)
			
			rm = result_peak.result
			
			pd = PassingData(min_MAF=self.min_MAF,\
							results_directory=self.results_directory, \
							need_chr_pos_ls=0,)
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