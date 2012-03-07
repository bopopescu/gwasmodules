#!/usr/bin/env python
"""

Examples:
	# compare cnv20 LD (min_score=4) vs call 32 LD (min_score=5) 
	TwoGWASPeakOverlap.py -z banyan -u yh -i 4634 -j 3395 -x 1 -y2
	
	# compare cnv80 LD (min_score=4) vs call 32 LD (min_score=5) 
	TwoGWASPeakOverlap.py -z banyan -u yh -i 4885 -j 3395 -x 1 -y2 -o /tmp/overlapStat.tsv
		
Description:
	filter the set of ResultPeak into two types:
		one covered by SNP array (snp probe within) -> calculate how many one GWAS's peaks fall into the others and vice versa
		one not covered by SNP-array

"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, getopt
import warnings, traceback
from pymodule import ProcessOptions, PassingData, GenomeDB
import Stock_250kDB
from pymodule import SNP, getListOutOfStr, yh_matplotlib, PassingData
from pymodule.db import formReadmeObj
from pymodule.CNV import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
from pymodule.RBTree import RBDict

class TwoGWASPeakOverlap(object):
	__doc__ = __doc__
	option_default_dict = {
						('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['stock_250k', 'd', 1, 'stock_250k database name', ],\
						('genome_dbname', 1, ): ['genome', 'g', 1, 'genome database name', ],\
						('schema', 0, ): ['', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('result1_id', 1, int): [None, 'i', 1, 'ResultMethod.id for result1'],\
						('result2_id', 1, int): [None, 'j', 1, ''],\
						('result1_peak_type_id', 1, int): [None, 'x', 1, 'peak type id for result1'],\
						('result2_peak_type_id', 1, int): [None, 'y', 1, ''],\
						('peakPadding', 1, int): [10000, 'e', 1, 'the extension for each peak on both sides. Rationale is if two peaks are '],\
						('logFilename', 0, ): [None, 'l', 1, 'file to contain logs. use it only if this program is at the end of pegasus workflow'],\
						('outputFname', 0, ): [None, 'o', 1, 'if given, output the statistics.'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.'],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
						}
	#('call_method_id', 0, int):[None, 'l', 1, 'Restrict results based on this call_method. Default is no such restriction.'],\
	#('analysis_method_id_ls', 0, ):['1,7', 'a', 1, 'Restrict results based on these analysis_methods. coma or dash-separated list'],\
	#("phenotype_method_id_ls", 0, ): [None, 'e', 1, 'comma/dash-separated phenotype_method id list, like 1,3-7. Default is all.'],\
	
	def __init__(self,  **keywords):
		"""
		2008-08-19
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
	
	def constructRBDictFromResultPeak(self, result_id, result_peak_type_id, peakPadding=10000):
		"""
		2011-10-16
		"""
		sys.stderr.write("Constructing RBDict for peaks from result %s, (peak type %s) ..."%(result_id, result_peak_type_id))
		result_peakRBDict = RBDict()
		query = Stock_250kDB.ResultPeak.query.filter_by(result_id=result_id).filter_by(result_peak_type_id=result_peak_type_id)
		counter = 0
		real_counter = 0
		for row in query:
			counter += 1
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=row.chromosome, \
							span_ls=[max(1, row.start - peakPadding), row.stop + peakPadding], \
							min_reciprocal_overlap=1,)
							#2010-8-17 overlapping keys are regarded as separate instances as long as they are not identical.
			if segmentKey not in result_peakRBDict:
				result_peakRBDict[segmentKey] = []
			result_peakRBDict[segmentKey].append(row)
		sys.stderr.write("%s peaks. Done.\n"%counter)
		return result_peakRBDict
	
	def calculateNoOfResult1PeaksCoveredByResult2(self, result1_peakRBDict, result2_peakRBDict):
		"""
		2011-10-16
		"""
		compareIns = CNVCompare(min_reciprocal_overlap=0.0000001)	#any overlap is an overlap
		no_of_peaks_not_in_result2 = 0
		overlap_ls = []
		counter = 0
		for node in result1_peakRBDict:
			counter += 1
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=node.key.chromosome, \
							span_ls=[node.key.start, node.key.stop], \
							min_reciprocal_overlap=0.0000001, )	#min_reciprocal_overlap doesn't matter here.
				# it's decided by compareIns.
			node_ls = []
			result2_peakRBDict.findNodes(segmentKey, node_ls=node_ls, compareIns=compareIns)
			total_perc_overlapped_by_result2 = 0.
			for node in node_ls:
				overlap1, overlap2, overlap_length, overlap_start_pos, overlap_stop_pos = get_overlap_ratio(segmentKey.span_ls, \
										[node.key.start, node.key.stop])[:5]
				total_perc_overlapped_by_result2 += overlap1
			if total_perc_overlapped_by_result2==0:
				no_of_peaks_not_in_result2 += 1
				overlap_ls.append(-0.5)
			else:
				overlap_ls.append(total_perc_overlapped_by_result2)
		return counter - no_of_peaks_not_in_result2
	
	def comparePeaksFromTwoAssociationResults(self, db_250k, result1_id=None, result1_peak_type_id=None,\
											result2_id=None, result2_peak_type_id=None, peakPadding=10000,\
											outputDir='/Network/Data/250k/tmp-yh/CNV/', outputFname=None):
		"""
		2011-4-22
			peak data is from table ResultPeak.
			The function tells how many peaks from result1 are not found, and how many are found
				in result2.
		"""
		sys.stderr.write("Comparing peaks from result %s against result %s ...\n"%(result1_id, result2_id))
		
		result1_peakRBDict = self.constructRBDictFromResultPeak(result1_id, result1_peak_type_id, peakPadding=peakPadding)
		result2_peakRBDict = self.constructRBDictFromResultPeak(result2_id, result2_peak_type_id, peakPadding=peakPadding)
		
		no_of_result1_peaks = len(result1_peakRBDict)
		no_of_result2_peaks = len(result2_peakRBDict)
		no_of_result1_peaks_in_result2 = self.calculateNoOfResult1PeaksCoveredByResult2(result1_peakRBDict, result2_peakRBDict)
		no_of_result2_peaks_in_result1 = self.calculateNoOfResult1PeaksCoveredByResult2(result2_peakRBDict, result1_peakRBDict)
		
		
		sys.stderr.write("%s/%s peaks in result %s found in result %s.\n"%(no_of_result1_peaks_in_result2,\
								no_of_result1_peaks, result1_id, result2_id))
		sys.stderr.write("%s/%s peaks in result %s found in result %s.\n"%(no_of_result2_peaks_in_result1,\
								no_of_result2_peaks, result2_id, result1_id))
		if outputFname:
			writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
			overlapFraction = (no_of_result1_peaks_in_result2 + no_of_result2_peaks_in_result1)/float(no_of_result1_peaks + no_of_result2_peaks)
			result1 = Stock_250kDB.ResultsMethod.get(result1_id)
			result2 = Stock_250kDB.ResultsMethod.get(result2_id)
			
			writer.writerow(['result1_id', 'result2_id|string', \
							'result1_category_id', 'result2_category_id', \
							"no_of_result1_peaks", "no_of_result2_peaks", \
							"no_of_result1_peaks_in_result2", "no_of_result2_peaks_in_result1", \
							"overlapFraction"])
			writer.writerow(['%s_%s'%(result1_id, result1.phenotype_method.short_name), \
							'%s_%s'%(result2_id, result2.phenotype_method.short_name), \
							result1.phenotype_method.biology_category_id,\
							result2.phenotype_method.biology_category_id,\
							no_of_result1_peaks, no_of_result2_peaks, \
							no_of_result1_peaks_in_result2, no_of_result2_peaks_in_result1, \
							overlapFraction])
		
		"""
		from pymodule import yh_matplotlib
		outputFname  = os.path.join(outputDir, 'hist_of_result%s_overlap_rate_in_result%s.png'%(result1_id, result2_id))
		yh_matplotlib.drawHist(overlap_ls, title='Histogram of overlap rate of %s peaks'%(len(overlap_ls)), \
							xlabel_1D='overlap rate',\
							outputFname=outputFname, \
							min_no_of_data_points=10, needLog=False)
		
		outputFname  = os.path.join(outputDir, 'hist_of_result%s_overlap_rate_in_result%s_2D.png'%(result1_id, result2_id))
		#outputFname = 'hist_of_result%s_overlap_rate_in_result%s_2D.png'%(result1_id, result2_id)
		C_ls = [1]*len(overlap_ls)
		colorBarLabel='log(count)'
		reduce_C_function = yh_matplotlib.logSum
		yh_matplotlib.drawHexbin(overlap_ls, score_ls, C_ls, fig_fname=outputFname, gridsize=10, \
								title='%s peaks in result %s vs %s'%(counter, result1_id, result2_id), \
								xlabel = 'overlap rate', \
								ylabel = 'association-score',\
								colorBarLabel=colorBarLabel, reduce_C_function= reduce_C_function)
		"""
	"""
		# 2011-4-22
		result1_id = 4634	#cnv_20_LD_KW
		result1_peak_type_id = 1
		result2_id = 3395	#call_32_LD_KW
		result2_peak_type_id = 2
		CNV.comparePeaksFromTwoAssociationResults(db_250k, result1_id=result1_id, result1_peak_type_id=result1_peak_type_id,\
										result2_id=result2_id, result2_peak_type_id=result2_peak_type_id, peakPadding=0)
		sys.exit(0)
		
		# 2011-4-22
		result1_id = 3395	#call_32_LD_KW
		result1_peak_type_id = 2 #min_score=5
		result2_id = 3396	#call_32_LDV_KW
		result2_peak_type_id = 2	#min_score=5
		CNV.comparePeaksFromTwoAssociationResults(db_250k, result1_id=result1_id, result1_peak_type_id=result1_peak_type_id,\
										result2_id=result2_id, result2_peak_type_id=result2_peak_type_id, peakPadding=0)
		sys.exit(0)
		
		result1_id_ls = [4634, 4635, 4636, 4637]	#cnv_20_LD_KW, LDV, SD, SDV
		result1_peak_type_id = 1 #min_score=4
		result2_id_ls = [3395, 3396, 3397, 3398]	#call_32_LD_KW, LDV, SD, SDV
		result2_id = 3396	#call_32_LDV_KW
		result2_peak_type_id = 2	#min_score=5
		for i in xrange(len(result1_id_ls)):
			result1_id = result1_id_ls[i]
			result2_id = result2_id_ls[i]
		CNV.comparePeaksFromTwoAssociationResults(db_250k, result1_id=result1_id, result1_peak_type_id=result1_peak_type_id,\
			result2_id=result2_id, result2_peak_type_id=result2_peak_type_id, peakPadding=0)
		sys.exit(0)
	"""
	def run(self):
		"""
		2011-3-28
			Read in the association result
			
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_250k = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user, password=self.db_passwd, \
									hostname=self.hostname, database=self.dbname)
		db_250k.setup(create_tables=False)
		db_250k.session.begin()
		
		self.comparePeaksFromTwoAssociationResults(db_250k, result1_id=self.result1_id, result1_peak_type_id=self.result1_peak_type_id,\
			result2_id=self.result2_id, result2_peak_type_id=self.result2_peak_type_id, peakPadding=self.peakPadding,\
			outputFname=self.outputFname)
		
		db_250k.session.flush()
		if self.commit:
			db_250k.session.commit()
		

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = TwoGWASPeakOverlap
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()