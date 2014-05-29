#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s  -o /tmp/count_association_locus.h5 /Network/Data/250k/db/genome_wide_association_locus/type_1/15_251results_call75_ana1_type1.h5
		/Network/Data/250k/db/genome_wide_association_locus/type_3/11_251results_call75_ana1_type3.h5

Description:
	2012.11.22
	If "-i ..." is given, it is regarded as one of the input files (plus the ones in trailing arguments). 
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import copy, tables
from tables import UInt64Col, Float64Col, StringCol
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, getColName2IndexFromHeader, figureOutDelimiter,\
	yh_matplotlib, YHFile
from pymodule import AbstractMapper
from pymodule import AbstractMatrixFileWalker

class CountAssociationLocusTable(tables.IsDescription):
	"""
	2012.12.24 new PyTables-based table definition
	"""
	id = UInt64Col(pos=0)
	min_score = Float64Col(pos=1)
	min_overlap_ratio = Float64Col(pos=2)
	total_no_of_results = UInt64Col(pos=3)
	no_of_association_loci = UInt64Col(pos=4)
	call_method_id_ls = StringCol(1000, pos=5)
	cnv_method_id_ls = StringCol(1000, pos=6)
	phenotype_method_id_ls = StringCol(1000, pos=7)
	analysis_method_id_ls  = StringCol(1000, pos=8)
	
class CountAssociationLocus(AbstractMatrixFileWalker):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractMatrixFileWalker.option_default_dict)
	for key in [('samplingRate', 1, float), ('whichColumnHeader', 0, ), ('whichColumn', 0, int),\
			('minNoOfTotal', 1, int), ('maxNoOfTotal', 0, int), ('logY', 0, int), \
			('valueForNonPositiveYValue', 1, float), ('missingDataNotation', 0, )]:
		option_default_dict.pop(key)
	# change default of file-format
	option_default_dict[('inputFileFormat', 0, int)][0] = 2
	option_default_dict[('outputFileFormat', 0, int)][0] = 2
	#option_default_dict.update(AbstractMapper.db_option_dict.copy())
	option_default_dict.update({
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractMatrixFileWalker.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		#self.connectDB() called within its __init__()
	
	def processHeader(self, pdata=None, **keywords):
		"""
		2012.12.24 use YHFile
		2012.11.22 called right after the header of an input file is derived in fileWalker()
		
		"""
		AbstractMatrixFileWalker.processHeader(self, pdata=pdata, rowDefinition=CountAssociationLocusTable)
	
	def fileWalker(self, inputFname=None, afterFileFunction=None, processRowFunction=None , run_type=1, **keywords):
		"""
		2012.11.22
		"""
		sys.stderr.write("walking through %s ..."%(inputFname))
		counter = 0
		real_counter = 0
		noOfSampled = 0
		pdata = self.initiatePassingData()
		if afterFileFunction is None:
			afterFileFunction = self.afterFileFunction
		try:
			reader = YHFile(inputFname, openMode='r')
			#output the header to the output file if necessary 
			self.processHeader(pdata=pdata) #2012.8.13
			
			firstTableObject = reader.getTableObject()
			no_of_association_loci = firstTableObject.nrows
			
			dataTuple = (firstTableObject.getAttribute('min_score', -1.0), \
						firstTableObject.getAttribute('min_overlap_ratio', -1.0),\
						firstTableObject.getAttribute('total_no_of_results', 1),\
						no_of_association_loci,
						firstTableObject.getListAttributeInStr('call_method_id_ls'),\
						firstTableObject.getListAttributeInStr('cnv_method_id_ls'),\
						firstTableObject.getListAttributeInStr('phenotype_method_id_ls'),\
						firstTableObject.getListAttributeInStr('analysis_method_id_ls'),\
						)
			self.invariantPData.writer.writeCellList([dataTuple])
			del reader
			counter += 1
			real_counter += 1
			noOfSampled += 1
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
			sys.exit(3)
		if counter>0:
			fraction = float(noOfSampled)/float(counter)
		else:
			fraction = 0
		sys.stderr.write("%s/%s (%.3f) data sampled. real_counter=%s.\n"%(noOfSampled, counter, fraction, real_counter))
		

if __name__ == '__main__':
	main_class = CountAssociationLocus
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()