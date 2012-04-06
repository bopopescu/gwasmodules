#!/usr/bin/env python
"""
Examples:
	%s 
	
	
	%s  -w GenomeDBPassword -i /tmp/out -l 59444 -O /tmp/gw_LD_pattern_between_snp_and_peak_59444 -u yh -p stock_250k_db_password

Description:
	2012.3.7
		program to output id of loci within numerous result peaks.
		GenomeDB is needed only for the chromosome sizes.
		If inputFname doesn't exist, the program still exits 0 to not disrupt a cluster of jobs in FindGenomeWideLDPatternBetweenSNPsAndPeakWorkflow.py
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])
import matplotlib; matplotlib.use("Agg")

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from AbstractVariationMapper import AbstractVariationMapper
import h5py, numpy
from pymodule import GenomeDB
from pymodule import SNP 
from variation.src import Stock_250kDB

class DrawManhattanPlotForLDInHDF5(AbstractVariationMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVariationMapper.option_default_dict.copy()
	option_default_dict.update({
							('genome_drivername', 1,):['postgresql', '', 1, 'which type of database is the genome database? mysql or postgresql', ],\
							('genome_hostname', 1, ): ['uclaOffice', '', 1, 'hostname of the genome db server', ],\
							('genome_dbname', 1, ): ['vervetdb', '', 1, 'genome database name', ],\
							('genome_schema', 0, ): ['genome', '', 1, 'genome database schema name', ],\
							('genome_db_user', 1, ): ['yh', 'g', 1, 'genome database username', ],\
							('genome_db_passwd', 1, ): [None, 'w', 1, 'genome database password', ],\
							
							('datasetName', 1, ): ['correlation', 'N', 1, 'name of the dataset in the HDF5 input', ],\
							('peak_id_ls', 1, ): ['', 'l', 1, 'peaks to be highlighted. dash/comma-separated ResultPeak.id', ],\
							('inputFname', 1, ): ['', 'i', 1, 'the HDF5 file which contains pairwise correlation between two types of loci', ],\
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		AbstractVariationMapper.__init__(self, inputFnameLs, **keywords)
		self.peak_id_ls = getListOutOfStr(self.peak_id_ls, data_type=int)
		self.peak_id_ls.sort()
	
	def getLocusTypeIDFromInput(self, inputFname=None, datasetName=None):
		"""
		2012.3.9
			guess the locus_type_id of the 1st locus in inputFname by checking the first data entry against db
		"""
		import h5py
		f1 = h5py.File(inputFname, 'r')
		d1 = f1[datasetName]
		d1_length = d1.shape[0]
		locus_type_id=None
		if d1_length>1:
			fstLocusID = d1[0][0]
			locus = Stock_250kDB.Snps.get(fstLocusID)
			locus_type_id = locus.locus_type_id
		del f1, d1
		return locus_type_id
	
	def run(self):
		"""
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_250k = self.db_250k
		
		#construct the bands to be highlighted in manhattan plot
		highlightBandLs = []
		rm = None
		for peak_id in self.peak_id_ls:
			result_peak = Stock_250kDB.ResultPeak.get(peak_id)
			highlightBandLs.append([result_peak.chromosome, result_peak.start, result_peak.stop])
			#take the 1st result_peak's result as the association result to get locus_type_id
			if rm is None:
				rm = result_peak.result
		if not rm:
			sys.stderr.write("Error: no results_method (association result) fetched from db.\n")
			sys.exit(1)
		if self.inputFname and os.path.isfile(self.inputFname):
			locus_type_id = self.getLocusTypeIDFromInput(self.inputFname, datasetName=self.datasetName)
			
			pd = PassingData()
			if rm.cnv_method_id and not db_250k._cnv_id2chr_pos:
				db_250k.cnv_id2chr_pos = rm.cnv_method_id
				pd.db_id2chr_pos = db_250k.cnv_id2chr_pos
			elif rm.call_method_id:
				db_250k.snp_id2chr_pos = (False, locus_type_id)	#correspond to priorTAIRVersion, locus_type_id
				pd.db_id2chr_pos = db_250k.snp_id2chr_pos
			
			#need to setup a different db setting
			db_genome = GenomeDB.GenomeDatabase(drivername=self.genome_drivername, username=self.genome_db_user,
							password=self.genome_db_passwd, hostname=self.genome_hostname, database=self.genome_dbname, \
							schema=self.genome_schema)
			db_genome.setup(create_tables=False)
			
			gwr_name = ''
			gwr = SNP.GenomeWideResult(name=gwr_name, construct_chr_pos2index=False, \
					construct_data_obj_id2index=False)
			gwr.fillGenomeWideResultFromHDF5CorrelationFile(self.inputFname, datasetName=self.datasetName, pdata=pd)
			gwr.drawManhattanPlot(db_genome, outputFnamePrefix=self.outputFnamePrefix,\
								min_value=None, need_svg=False, ylim_type=2,\
								drawBonferroni=False, highlightBandLs=highlightBandLs)
		else:	#2012.3.28 input is invalid.
			sys.stderr.write("inputFname %s is not valid (inexistent).\n"%(self.inputFname))
			sys.exit(0)	#fake ok as I want pegasus workflow to keep running.
	
if __name__ == '__main__':
	main_class = DrawManhattanPlotForLDInHDF5
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()