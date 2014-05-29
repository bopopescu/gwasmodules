#!/usr/bin/env python
"""
Examples:
	# 2011-9-29
	%s -j 32 -l 80 -a 1 -A 1 -x 2 -X 1 -i 1 -o cm32am1peakType2_vs_cm80am1peakType1_peakOverlap.xml -u yh -z banyan
		-B condorpool -D condorpool
	%s -j 32 -l 80 -a 7 -A 32 -x 3 -X 3 -i 1 -o cm32am7peakType3_vs_cm80am32peakType3_peakOverlap.xml -u yh -z banyan
		-B condorpool -D condorpool
		
	# 2012.2.22
	%s -j 32 -l 80 -a 1 -A 1 -x 1 -X 3 -o workflow/cm32am1peakType1_vs_cm80am1peakType3_peakOverlap.xml -u yh -z banyan -s 1
		-B condorpool -D condorpool
	
Description:
	2011-10-12
		output a workflow that runs the PairwiseGWASPeakOverlap.py on all specified gwas results
		
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

from sqlalchemy.types import LargeBinary

bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:	   #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *
import Stock_250kDB
from PairwiseGWASPeakOverlapPipeline import PairwiseGWASPeakOverlapPipeline
from AbstractVariationWorkflow import AbstractVariationWorkflow

class PairwiseSamePhenotypeGWASPeakOverlapPipeline(PairwiseGWASPeakOverlapPipeline):
	__doc__ = __doc__
	option_default_dict = AbstractVariationWorkflow.option_default_dict.copy()
	option_default_dict.update(PairwiseGWASPeakOverlapPipeline.common_option_dict)

	option_default_dict.update({
						('call_method1_id', 0, int): [None, 'j', 1, ''],\
						('cnv_method1_id', 0, int):[None, 'c', 1, 'Restrict results to this cnv_method.'],\
						('call_method2_id', 0, int):[None, 'l', 1, ''],\
						('cnv_method2_id', 0, int):[None, 'q', 1, 'Restrict result 2 to this cnv_method.'],\
						('analysis_method1_id', 1, int):[1, 'a', 1, 'Restrict result 1 to this analysis_method'],\
						('analysis_method2_id', 1, int):[1, 'A', 1, 'Restrict result 2 to this analysis_method'],\
						("association1_peak_type_id", 1, int): [None, 'x', 1, 'PeakType for result1'],\
						("association2_peak_type_id", 1, int): [None, 'X', 1, 'PeakType for result2'],\
						})
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AbstractVariationWorkflow.__init__(self, **keywords)
		
		if getattr(self, 'phenotype_method_id_ls', None):
			self.phenotype_method_id_ls = getListOutOfStr(self.phenotype_method_id_ls, data_type=int)
			self.phenotype_method_id_ls.sort()
		else:
			self.phenotype_method_id_ls = []
		
	
	def getPhenotypeMethodId2ResultID(self, result_id_ls):
		"""
		"""
		sys.stderr.write("Getting phenotype_method_id2result_id given %s entries in result_id_ls ... "%(len(result_id_ls)))
		phenotype_method_id2result_id = {}
		for result_id in result_id_ls:
			rm = Stock_250kDB.ResultsMethod.get(result_id)
			phenotype_method_id2result_id[rm.phenotype_method_id] = rm.id
		sys.stderr.write("Done.\n")
		return phenotype_method_id2result_id
	
	def run(self):
		"""
		2011-9-28
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_250k = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user, password=self.db_passwd, \
									hostname=self.hostname, database=self.dbname)
		db_250k.setup(create_tables=False)
		
		
		sameCategoryPhenotypeMethodLs = db_250k.getPhenotypeMethodLsGivenBiologyCategoryID(self.biology_category_id, access=self.access)
		sameCategoryPhenotypeMethodIDLs = [pm.id for pm in sameCategoryPhenotypeMethodLs]
		#merge the two lists of phenotype method id together
		phenotype_method_id_ls = list(set(self.phenotype_method_id_ls + sameCategoryPhenotypeMethodIDLs))
		
		result1_query = db_250k.getResultLs(call_method_id=self.call_method1_id, analysis_method_id_ls=[self.analysis_method1_id], \
						phenotype_method_id_ls=phenotype_method_id_ls, cnv_method_id=self.cnv_method1_id)
		result2_query = db_250k.getResultLs(call_method_id=self.call_method2_id, analysis_method_id_ls=[self.analysis_method2_id], \
						phenotype_method_id_ls=phenotype_method_id_ls, cnv_method_id=self.cnv_method2_id)
		
		result1_id_ls = []
		for result in result1_query:
			result1_id_ls.append(result.id)
		
		result2_id_ls = []
		for result in result2_query:
			result2_id_ls.append(result.id)
		
		#make sure the entries with (result_id, self.association_peak_type_id) exists in AssociationPeak
		result1_id_ls = db_250k.filterResultIDLsBasedOnAssociationPeak(result1_id_ls, self.association1_peak_type_id)
		result2_id_ls = db_250k.filterResultIDLsBasedOnAssociationPeak(result2_id_ls, self.association2_peak_type_id)
		
		phenotype_method_id2result1_id = self.getPhenotypeMethodId2ResultID(result1_id_ls)
		phenotype_method_id2result2_id = self.getPhenotypeMethodId2ResultID(result2_id_ls)
		phenotype_method_id2result_id_pair = {}
		for phenotype_method_id , result1_id in phenotype_method_id2result1_id.iteritems():
			if phenotype_method_id in phenotype_method_id2result2_id:
				result2_id = phenotype_method_id2result2_id.get(phenotype_method_id)
				phenotype_method_id2result_id_pair[phenotype_method_id] = [result1_id, result2_id]
		
		# Create a abstract dag
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = self.initiateWorkflow(workflowName)
		
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		counter = 0
		
		overlapStatDir = "overlapStat"
		overlapStatDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=overlapStatDir)
		
		counter += 1
		"""
		overlapPlotDir = "overlapPlot"
		overlapPlotDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=overlapPlotDir)
		
		outputFnamePrefix = os.path.join(overlapPlotDir, 'cm%s_cnvM%s_am%s_peakType%s_vs_cm%s_cnvM%s_am%s_peakType%s_biologyCategory%s_overlapPeak'%\
								(self.call_method1_id, self.cnv_method1_id, self.analysis_method1_id, self.association1_peak_type_id, \
								self.call_method2_id, self.cnv_method2_id, self.analysis_method2_id, self.association2_peak_type_id, \
								self.biology_category_id))
		plotAssociationPeakOverlapJob = slef.addPlotPeakOverlapJob(workflow, executable=workflow.plotAssociationPeakOverlapJob, \
							outputFnamePrefix=outputFnamePrefix, \
							parentJobLs=[overlapPlotDirJob], job_max_memory=100, walltime = 60, \
							extraDependentInputLs=[], \
							transferOutput=True)
		"""
		#each contig in each trio gets a summary.
		peakOverlapStatMergeFile = File('peak_overlap_stat.tsv')
		peakOverlapStatMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=peakOverlapStatMergeFile, transferOutput=True, parentJobLs=[])
		counter += 1
		no_of_input =0
		for phenotype_method_id, result_id_pair in phenotype_method_id2result_id_pair.iteritems():
			result1_id = result_id_pair[0]
			result2_id = result_id_pair[1]
			outputFnamePrefix = 'result_%s_vs_%s_peak_type_%s_vs_%s'%(result1_id, result2_id, self.association1_peak_type_id, self.association2_peak_type_id)
			outputF = File(os.path.join(overlapStatDir, '%s.tsv'%(outputFnamePrefix)))
			
			if no_of_input==0:	#add one random input, otherwise replica catalog error occurs
				rm1 = Stock_250kDB.ResultsMethod.get(result1_id)
				inputFile1 = self.registerOneInputFile(workflow, rm1.filename)
				extraDependentInputLs=[inputFile1]
			else:
				extraDependentInputLs=[]
			
			no_of_input += 1
			gwasPeakOverlapJob = self.addGWASPeakOverlapJob(workflow, executable=workflow.twoGWASPeakOverlap, \
							result1_id=result1_id, result2_id=result2_id, association1_peak_type_id=self.association1_peak_type_id, \
							association2_peak_type_id=self.association2_peak_type_id, peak_padding=self.peak_padding, \
							outputF=outputF, \
							commit=1, results_directory=None, logFile=None, \
							parentJobLs=[overlapStatDirJob], job_max_memory=100, walltime = 60, \
							extraDependentInputLs=extraDependentInputLs, \
							transferOutput=True)
			counter += 1
			
			self.addInputToStatMergeJob(workflow, statMergeJob=peakOverlapStatMergeJob, \
									inputF=outputF, parentJobLs=[gwasPeakOverlapJob])
			"""
			self.addInputToStatMergeJob(workflow, statMergeJob=plotAssociationPeakOverlapJob, inputF=outputF, \
							parentJobLs=[gwasPeakOverlapJob])
			"""
			
		sys.stderr.write("%s jobs.\n"%(counter))
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
	
if __name__ == '__main__':
	main_class = PairwiseSamePhenotypeGWASPeakOverlapPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
