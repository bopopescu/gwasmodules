#!/usr/bin/env python
"""
Examples:
	# 2011-10-17
	%s -i 1 -x 1 -c 20 -a 1 -o cnv20Analysis1BiologyCategory1_peakType1_pairwisePeakOverlap.xml -B condorpool -D condorpool
		-u yh -z banyan
	
	%s -x 1 -E 80 -a 1 -o call80Analysis1_peakType1_pairwisePeakOverlap.xml -B condorpool -D condorpool
		-u yh -z banyan
	
	%s -x 3 -E 80 -a 32 -o call80Analysis32_peakType3_pairwisePeakOverlap.xml -B condorpool -D condorpool
		-u yh -z banyan
	
	%s -x 2 -E 32 -a 1 -o call32Analysis1_peakType2_pairwisePeakOverlap.xml -B condorpool -D condorpool
		-u yh -z banyan
	
	#2012.2.28
	%s -x 2 -E 32 -a 1 -i 1 -o call32Analysis1_peakType2_biologyCategory1_publicPhenotype_pairwisePeakOverlap.xml
		-B condorpool -D condorpool -u yh -z banyan -s1
	
	#2012.2.28
	%s -x 1 -E 80 -a 1 -i 1 -o call80Analysis1_peakType1_biologyCategory1_publicPhenotype_pairwisePeakOverlap.xml
		-B condorpool -D condorpool -u yh -z banyan -s1
	
Description:
	2011-10-12
		output a workflow that runs the PairwiseGWASPeakOverlap.py on all specified gwas results
		
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *
import Stock_250kDB
from AbstractVariationWorkflow import AbstractVariationWorkflow

class PairwiseGWASPeakOverlapPipeline(AbstractVariationWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVariationWorkflow.option_default_dict.copy()
	common_option_dict = {
						("phenotype_method_id_ls", 0, ): [None, 'y', 1, 'comma/dash-separated phenotype_method id list, like 1,3-7. Default is all.'],\
						('biology_category_id', 0, int): [None, 'i', 1, 'phenotype biology category id. Default is no filter'],\
						('access', 0, int): [None, 's', 1, 'Restrict phenotype via access field in db. 1: public phenotypes, 2: restricted. Default is no filter'],\
						('peak_padding', 1, int): [10000, '', 1, 'the extension for each peak on both sides. Rationale is if two peaks are '],\
						}
	option_default_dict.update(common_option_dict)
	option_default_dict.update({
						('result_id_ls', 0, ): [None, 'w', 1, 'comma or dash-separated list of result ids, i.e. 3431-3438,4041'],\
						('call_method_id', 0, int):[None, 'E', 1, 'Restrict results based on this call_method. Default is no such restriction.'],\
						('cnv_method_id', 0, int):[None, 'c', 1, 'Restrict results based on this cnv_method. Default is no such restriction.'],\
						('analysis_method_id_ls', 0, ):[1, 'a', 1, 'Restrict results based on these analysis_methods. coma or dash-separated list'],\
						('association_peak_type_id', 1, int): [None, 'x', 1, 'peak type id for result peaks'],\
						})
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AbstractVariationWorkflow.__init__(self, **keywords)
		
		if getattr(self, 'result_id_ls', None):
			self.result_id_ls = getListOutOfStr(self.result_id_ls, data_type=int)
			self.result_id_ls.sort()
		else:
			self.result_id_ls = []
		
		if getattr(self, 'analysis_method_id_ls', None):
			self.analysis_method_id_ls = getListOutOfStr(self.analysis_method_id_ls, data_type=int)
			self.analysis_method_id_ls.sort()
		else:
			self.analysis_method_id_ls = []
		if getattr(self, 'phenotype_method_id_ls', None):
			self.phenotype_method_id_ls = getListOutOfStr(self.phenotype_method_id_ls, data_type=int)
			self.phenotype_method_id_ls.sort()
		else:
			self.phenotype_method_id_ls = []
	
	
	def registerCustomExecutables(self, workflow):
		"""
		2012.2.15
		"""
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		
		twoGWASPeakOverlap = Executable(namespace=namespace, name="TwoGWASPeakOverlap", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		twoGWASPeakOverlap.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "TwoGWASPeakOverlap.py"), site_handler))
		#twoGWASPeakOverlap.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(twoGWASPeakOverlap)
		workflow.twoGWASPeakOverlap = twoGWASPeakOverlap
		
		PlotAssociationPeakOverlap = Executable(namespace=namespace, name="PlotAssociationPeakOverlap", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		PlotAssociationPeakOverlap.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "plot/PlotAssociationPeakOverlap.py"), site_handler))
		PlotAssociationPeakOverlap.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(PlotAssociationPeakOverlap)
		workflow.plotAssociationPeakOverlap = PlotAssociationPeakOverlap
	
	def addGWASPeakOverlapJob(self, workflow, executable=None, \
							result1_id=None, result2_id=None, association1_peak_type_id=None, \
							association2_peak_type_id=None, peak_padding=None, outputF=None, \
							commit=0, results_directory=None, logFile=None, \
							parentJobLs=[], job_max_memory=100, job_max_walltime = 60, \
							extraDependentInputLs=[], \
							transferOutput=False, **keywords):
		"""
		2012.2.22
			job_max_walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		#apply int because result1_id is of long int. and casting "434L" to integer will raise exception
		job.addArguments("-v", self.drivername, "-z", self.hostname, "-d", self.dbname, \
						"-u", self.db_user, "-p", self.db_passwd,\
						"-i", repr(int(result1_id)), "-j", repr(int(result2_id)), \
						"-x", repr(association1_peak_type_id), "-y", repr(association2_peak_type_id), \
						"-e", repr(peak_padding), "-o", outputF)
		job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		if commit:
			job.addArguments("-c")
		if results_directory:
			job.addArguments("-t", results_directory,)
		if self.schema:
			job.addArguments("-k", self.schema,)
		if logFile:
			job.addArguments("--logFilename=%s"%(logFile.name))
			job.uses(logFile, transfer=transferOutput, register=transferOutput, link=Link.OUTPUT)
			job.output = logFile
		workflow.addJob(job)
		
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory, max_walltime=job_max_walltime)
		return job
	
	def addPlotPeakOverlapJob(self, workflow, executable=None, \
							outputFnamePrefix=None, \
							parentJobLs=[], job_max_memory=100, job_max_walltime = 60, \
							extraDependentInputLs=[], \
							transferOutput=False, **keywords):
		"""
		2012.2.22
			job_max_walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		
		job.addArguments("-o", outputFnamePrefix)
		
		job.uses(File("%s_hist_of_fraction_of_association1_peaks_in_result2.png"%outputFnamePrefix), \
							transfer=True, register=True, link=Link.OUTPUT)
		job.uses(File("%s_hist_of_fraction_of_association2_peaks_in_result1.png"%outputFnamePrefix), \
							transfer=True, register=True, link=Link.OUTPUT)
		job.uses(File("%s_hist_of_fraction_of_recurrent_peaks_in_combined.png"%outputFnamePrefix), \
							transfer=True, register=True, link=Link.OUTPUT)
		job.uses(File("%s_no_of_peaks_result1_vs_result2.png"%outputFnamePrefix), \
							transfer=True, register=True, link=Link.OUTPUT)
		job.uses(File("%s_result1_no_of_peak_vs_fraction.png"%outputFnamePrefix), \
							transfer=True, register=True, link=Link.OUTPUT)
		job.uses(File("%s_result2_no_of_peak_vs_fraction.png"%outputFnamePrefix), \
							transfer=True, register=True, link=Link.OUTPUT)
		job.uses(File("%s_1_fraction_in2_vs_2_fraction_in1.png"%outputFnamePrefix), \
							transfer=True, register=True, link=Link.OUTPUT)
		job.uses(File("%s_combined_no_of_peak_vs_fraction.png"%outputFnamePrefix), \
							transfer=True, register=True, link=Link.OUTPUT)
		
		workflow.addJob(job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory, max_walltime=job_max_walltime)
		return job

	
	def run(self):
		"""
		2011-9-28
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_250k = self.db_250k
		
		sameCategoryPhenotypeMethodLs = db_250k.getPhenotypeMethodLsGivenBiologyCategoryID(self.biology_category_id, access=self.access)
		sameCategoryPhenotypeMethodIDLs = [pm.id for pm in sameCategoryPhenotypeMethodLs]
		
		#merge the two lists of phenotype method id together
		phenotype_method_id_ls = list(set(self.phenotype_method_id_ls + sameCategoryPhenotypeMethodIDLs))
		phenotype_method_id_ls.sort()
		
		result_query = db_250k.getResultLs(call_method_id=self.call_method_id, analysis_method_id_ls=self.analysis_method_id_ls, \
						phenotype_method_id_ls=phenotype_method_id_ls, cnv_method_id=self.cnv_method_id)
		result_id_ls = self.result_id_ls
		for result in result_query:
			result_id_ls.append(result.id)
		
		#make sure the entries with (result_id, self.association_peak_type_id) exists in AssociationPeak
		result_id_ls = db_250k.filterResultIDLsBasedOnAssociationPeak(result_id_ls, self.association_peak_type_id)
		
		# Create a abstract dag
		workflow = self.initiateWorkflow()
		
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		overlapStatDir = "overlapStat"
		overlapStatDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=overlapStatDir)
		overlapPlotDir = "overlapPlot"
		overlapPlotDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=overlapPlotDir)
		
		analysis_method_id_ls = map(str, self.analysis_method_id_ls)
		outputFnamePrefix = os.path.join(overlapPlotDir, 'callMethod%s_cnvMethod%s_analysisMethod%s_biologyCategory%s_peakType%s_overlapPeak'%\
								(self.call_method_id, self.cnv_method_id, '_'.join(analysis_method_id_ls), self.biology_category_id,\
								self.association_peak_type_id))
		plotAssociationPeakOverlapJob = self.addPlotPeakOverlapJob(workflow, executable=workflow.plotAssociationPeakOverlap, \
							outputFnamePrefix=outputFnamePrefix, \
							parentJobLs=[overlapPlotDirJob], job_max_memory=100, job_max_walltime = 60, \
							extraDependentInputLs=[], \
							transferOutput=True)
		
		counter = 0
		no_of_input =0
		for i in xrange(len(result_id_ls)):
			for j in range(i+1, len(result_id_ls)):
				result1_id = result_id_ls[i]
				result2_id = result_id_ls[j]
				outputFnamePrefix = 'result_%s_vs_%s_peak_type_%s'%(result1_id, result2_id, self.association_peak_type_id)
				outputF = File(os.path.join(overlapStatDir, '%s.tsv'%(outputFnamePrefix)))
				if no_of_input==0:	#add one random input, otherwise replica catalog error occurs
					rm1 = Stock_250kDB.ResultsMethod.get(result1_id)
					inputFile1 = self.registerOneInputFile(workflow, rm1.filename)
					extraDependentInputLs=[inputFile1]
				else:
					extraDependentInputLs=[]
				
				no_of_input += 1
				gwasPeakOverlapJob = self.addGWASPeakOverlapJob(workflow, executable=workflow.twoGWASPeakOverlap, \
							result1_id=result1_id, result2_id=result2_id, association1_peak_type_id=self.association_peak_type_id, \
							association2_peak_type_id=self.association_peak_type_id, peak_padding=self.peak_padding, \
							outputF=outputF, \
							commit=1, results_directory=None, logFile=None, \
							parentJobLs=[overlapStatDirJob], job_max_memory=100, job_max_walltime = 60, \
							extraDependentInputLs=extraDependentInputLs, \
							transferOutput=True)
				counter += 1
				
				self.addInputToStatMergeJob(workflow, statMergeJob=plotAssociationPeakOverlapJob, inputF=outputF, \
							parentJobLs=[gwasPeakOverlapJob])
		
		sys.stderr.write("%s gwas peak overlap jobs.\n"%(counter))
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		self.writeXML(outf)


if __name__ == '__main__':
	main_class = PairwiseGWASPeakOverlapPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
