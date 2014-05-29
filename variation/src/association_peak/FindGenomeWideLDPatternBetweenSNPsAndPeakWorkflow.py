#!/usr/bin/env python
"""
Examples:
	# 2012.3.12 only phenotype 1 (-y1, "-j -1" is used to match no biology category. default is all categories).
	%s -i LD_between_call_32_and_80.2012.3.6T1846/LD/ -o GW_LD_pattern_between_call32_call80Type1Peaks_phenotype1LD.xml 
		-p STOCK_250K_PASSWD -u yh -z banyan -w GENOME_DB_PASSWD -x 1 -l 80 -s 1 -j -1 -y1
		-F GWLD -B condorpool -D condorpool
	
	# 2012.3.12 all public phenotypes ( -s 1,), use 50M as chunkSize, only kruskal wallis results, type 1(-x1 min-score=4)
	%s -i LD_between_call_32_and_80.2012.3.6T1846/LD/ -o GW_LD_pattern_between_call32_call80Type1Peaks_KW.xml 
		-p STOCK_250K_PASSWD -u yh -z banyan -w GENOME_DB_PASSWD -x 1 -l 80 -s 1 -a1
		-F GWLD -B condorpool -D condorpool -f 50000000
	
	# 2012.3.12 all public phenotypes ( -s 1,), use 50M as chunkSize, EMMAX results, type 3 (-x 3, min-score=3)
	%s-i LD_between_call_32_and_80.2012.3.6T1846/LD/ -o GW_LD_pattern_between_call32_call80Type3Peaks_EMMAX.xml 
		-p STOCK_250K_PASSWD -u yh -z banyan -w GENOME_DB_PASSWD -x 3 -l 80 -s 1 -a32
		-F GWLD -B condorpool -D condorpool -f 50000000
	
Description:
	2012.3.8
		1. note: OR between argument biology_category_id and phenotype_method_id_ls 
	GenomeDB is needed only for the chromosome sizes.


"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], )

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
import h5py
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, figureOutDelimiter
from variation.src import Stock_250kDB
from variation.src.pegasus.AbstractVariationWorkflow import AbstractVariationWorkflow

class FindGenomeWideLDPatternBetweenSNPsAndPeakWorkflow(AbstractVariationWorkflow):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractVariationWorkflow.option_default_dict)
	common_option_dict = {
						('genome_drivername', 1,):['postgresql', '', 1, 'which type of database is the genome database? mysql or postgresql', ],\
						('genome_hostname', 1, ): ['uclaOffice', '', 1, 'hostname of the genome db server', ],\
						('genome_dbname', 1, ): ['vervetdb', '', 1, 'genome database name', ],\
						('genome_schema', 0, ): ['genome', '', 1, 'genome database schema name', ],\
						('genome_db_user', 1, ): ['yh', 'g', 1, 'genome database username', ],\
						('genome_db_passwd', 1, ): [None, 'w', 1, 'genome database password', ],\
						
						('datasetName', 1, ): ['correlation', 'N', 1, 'name of the dataset in the HDF5 input', ],\
						
						("inputFolder", 1, ): [None, 'i', 1, 'the input folder which contains all the hdf5 correlation between 2 types of loci'],\
						('chunkSize', 0, int): [10000000, 'f', 1, 'Loci from each will be partitioned into chunks of 5000'],\
						('biology_category_id', 0, int): [None, 'j', 1, 'phenotype biology category id. Default is no filter.'],\
						('call_method_id', 0, int):[None, 'l', 1, 'peaks are based on this call method and \
				also 2nd type of loci are also from this call method'],\
						('cnv_method_id', 0, int):[None, 'A', 1, 'same functionality as call_method_id'],\
						('analysis_method_id_ls', 0, ):[None, 'a', 1, 'Restrict results based on these analysis_methods. coma or dash-separated list'],\
						("phenotype_method_id_ls", 0, ): [None, 'y', 1, 'comma/dash-separated phenotype_method id list, like 1,3-7. Default is all.'],\
						('pegasusFolderName', 0, ): ['GWLD', 'F', 1, 'the folder relative to pegasus workflow root to contain input & output.\
								It will be created during the pegasus staging process. It is useful to separate multiple workflows.\
								If empty, everything is in the pegasus root.', ],\
						('result_peak_type_id', 1, int): [None, 'x', 1, 'peak type id for result peaks'],\
						('access', 0, int): [1, 's', 1, 'Restrict phenotype via access field in db. 1: public phenotypes, 2: restricted. Default is no filter'],\
						}
	
	option_default_dict.update(common_option_dict)
	
	def __init__(self,  **keywords):
		"""
		2012.2.28
		"""
		AbstractVariationWorkflow.__init__(self, **keywords)
		self.inputFolder = os.path.abspath(self.inputFolder)
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
		2012.3.2
		"""
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		
		#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		executableClusterSizeMultiplierList = []
		
		FindMaxLDBetweenPeakAndEachLocus = Executable(namespace=namespace, name="FindMaxLDBetweenPeakAndEachLocus", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		FindMaxLDBetweenPeakAndEachLocus.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/FindMaxLDBetweenPeakAndEachLocus"), \
									site_handler))
		executableClusterSizeMultiplierList.append((FindMaxLDBetweenPeakAndEachLocus, 0.3))
		
		
		MergeTwoLocusCorrelationHDF5 = Executable(namespace=namespace, name="MergeTwoLocusCorrelationHDF5", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		MergeTwoLocusCorrelationHDF5.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/MergeTwoLocusCorrelationHDF5.py"), \
									site_handler))
		executableClusterSizeMultiplierList.append((MergeTwoLocusCorrelationHDF5, 0.3))
		
		OutputLociIDOfResultPeakInHDF5 = Executable(namespace=namespace, name="OutputLociIDOfResultPeakInHDF5", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		OutputLociIDOfResultPeakInHDF5.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "association_peak/OutputLociIDOfResultPeakInHDF5.py"), \
									site_handler))
		executableClusterSizeMultiplierList.append((OutputLociIDOfResultPeakInHDF5, 0.3))
		
		DrawManhattanPlotForLDInHDF5 = Executable(namespace=namespace, name="DrawManhattanPlotForLDInHDF5", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		DrawManhattanPlotForLDInHDF5.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "plot/DrawManhattanPlotForLDInHDF5.py"), \
									site_handler))
		executableClusterSizeMultiplierList.append((DrawManhattanPlotForLDInHDF5, 0.4))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	
	def addOutputLociIDOfResultPeakInHDF5Job(self, workflow, executable=None, peak_id=None, outputFile=None,\
					parentJobLs=[], extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
					job_max_memory=10, **keywords):
		"""
		2012.3.10
			-i 59444 -u yh -z banyan -o /tmp/peak_59444.h5
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments("-v", self.drivername, "-z", self.hostname, "-d", self.dbname, \
						"-u", self.db_user, "-p", self.db_passwd,\
						"-i", repr(int(peak_id)), "-o", outputFile)
		
		if extraArguments:
			job.addArguments(extraArguments)
		job.uses(outputFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = outputFile
		workflow.addJob(job)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job
	
	def getNoOfRowsFromHDF5Data(self, inputFname, datasetName="correlation"):
		"""
		2012.3.10
		"""
		f1 = h5py.File(inputFname, 'r')
		d1 = f1[datasetName]
		d1_length = d1.shape[0]
		del f1
		return d1_length
	
	def addDrawManhattanPlotForLDInHDF5Job(self, workflow, executable=None, correlationFile=None, peak_id=None, \
										datasetName=None, outputFile=None,\
					outputFnamePrefix=None, parentJobLs=[], extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
					job_max_memory=10, **keywords):
		"""
		2012.3.10
			DrawManhattanPlotForLDInHDF5.py -w secret -i /tmp/output.2.h5 -l 59444 -N correlation
				-O /tmp/gw_LD_pattern_between_snp_and_peak_59444 -u yh -p secret
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments("-v", self.drivername, "-z", self.hostname, "-d", self.dbname, \
						"-u", self.db_user, "-p", self.db_passwd, "-w", self.genome_db_passwd,\
						"-l %s"%(peak_id),"-N", datasetName, "-i", correlationFile, "-O", outputFnamePrefix)
		if extraArguments:
			job.addArguments(extraArguments)
		job.uses(correlationFile, transfer=True, register=True, link=Link.INPUT)
		job.uses(outputFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = outputFile
		workflow.addJob(job)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job
	
	def addFindMaxLDBetweenPeakAndEachLocusJob(self, workflow, executable=None, correlationFile=None, peakLociH5File=None, \
											outputFile=None, row_start=None, row_stop=None, \
					parentJobLs=[], extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
					job_max_memory=100, **keywords):
		"""
		2012.3.10
			FindMaxLDBetweenPeakAndEachLocus -j /tmp/peak_59444.h5
				-s 0 -t 1000
				-i /Network/Data/250k/tmp-yh/pegasus/LD_between_call_32_and_80.2012.3.9T2005/LD/cor_i1_0_4999_i2_0_4999.h5
				-o /tmp/output.3.h5
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments("-i", correlationFile, "-j", peakLociH5File, "-o", outputFile)
		if row_start is not None:
			job.addArguments("-s %s"%(row_start))
		if row_stop is not None:
			job.addArguments("-t %s"%(row_stop))
		if extraArguments:
			job.addArguments(extraArguments)
		job.uses(peakLociH5File, transfer=False, register=True, link=Link.INPUT)
		job.uses(correlationFile, transfer=True, register=True, link=Link.INPUT)
		job.uses(outputFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = outputFile
		workflow.addJob(job)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job
		
	
	def mapMaxLDJobsGivenInputData(self, workflow, inputData=None, datasetName=None, peak_id=None, \
								peakLociH5File=None, outputFile=None, outputDirJob=None, \
						chunkSize=None, parentJobLs=[], extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
						job_max_memory=10, **keywords):
		"""
		2012.3.10
			for each input hdf5, figure out how many data it has, split by chunks
		"""
		redundantPeakCorrelationFile = File(os.path.join(outputDirJob.folder, 'redundantMaxCorrelationBetweenFstLociAndPeak%s.h5'%(peak_id)))
		
		correlationMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.MergeTwoLocusCorrelationHDF5, \
				outputF=redundantPeakCorrelationFile, transferOutput=False, extraArguments='-d %s'%(datasetName), parentJobLs=[outputDirJob])
		no_of_jobs = 1
		for jobData in inputData.jobDataLs:
			inputFile = jobData.output
			no_of_rows = self.getNoOfRowsFromHDF5Data(inputFile.abspath, datasetName=datasetName)
			for i1_start in range(0, no_of_rows, chunkSize):
				i1_stop = min(i1_start + chunkSize -1, no_of_rows-1)
				subOutputFile = File(os.path.join(outputDirJob.folder, 'maxCorrelationBetweenFstLociAndPeak%s_%s_%s.h5'%(peak_id, i1_start, i1_stop)))
				maxLDJob = self.addFindMaxLDBetweenPeakAndEachLocusJob(workflow, executable=workflow.FindMaxLDBetweenPeakAndEachLocus, \
											correlationFile=inputFile, \
											peakLociH5File=peakLociH5File, \
											outputFile=subOutputFile, row_start=i1_start, row_stop=i1_stop, \
					parentJobLs=jobData.jobLs + [outputDirJob] + parentJobLs, extraDependentInputLs=[], transferOutput=False, extraArguments=None, \
					job_max_memory=job_max_memory)
				no_of_jobs += 1
				self.addInputToStatMergeJob(workflow, statMergeJob=correlationMergeJob, \
								inputF=maxLDJob.output, parentJobLs=[maxLDJob])
		
		maxLDJob = self.addFindMaxLDBetweenPeakAndEachLocusJob(workflow, executable=workflow.FindMaxLDBetweenPeakAndEachLocus, \
											correlationFile=correlationMergeJob.output, \
											peakLociH5File=peakLociH5File, \
											outputFile=outputFile, row_start=None, row_stop=None, \
					parentJobLs=[correlationMergeJob], extraDependentInputLs=[], transferOutput=transferOutput, extraArguments=None, \
					job_max_memory=job_max_memory)
		no_of_jobs += 1
		maxLDJob.no_of_jobs = no_of_jobs
		return maxLDJob
		
	
	def addJobs(self, workflow=None, result_peak_ls=None, inputData=None, datasetName=None, chunkSize=None, pegasusFolderName=""):
		"""
		2012.3.3
		"""
		if workflow is None:
			workflow = self
		
		returnJobData = PassingData()
		
		no_of_jobs = 0
		
		topOutputDir = pegasusFolderName
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		no_of_jobs += 1
		
		biology_category_id2peak_data = {}
		no_of_peaks = 0
		for row in result_peak_ls:
			biology_category_id = row.result.phenotype_method.biology_category_id
			if biology_category_id not in biology_category_id2peak_data:
				biology_category_id2peak_data[biology_category_id] = PassingData(job=None, result_peak_ls=[])
				#add a mkdirJob
				folderName = os.path.join(topOutputDir, 'biology_category_%s'%(biology_category_id))
				folderJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=folderName,\
												parentJobLs=[topOutputDirJob])
				biology_category_id2peak_data[biology_category_id].job = folderJob
				no_of_jobs += 1
			biology_category_id2peak_data[biology_category_id].result_peak_ls.append(row)
			no_of_peaks +=1
		
		sys.stderr.write("%s peaks. %s biology categories.\n"%(no_of_peaks, len(biology_category_id2peak_data)))
		
		sys.stderr.write("Finding max LD between one type of loci and peaks on %s input correlation files ... "%(len(inputData.jobDataLs)))
		
		prevReportedStr = ""
		for biology_category_id, peak_data in biology_category_id2peak_data.iteritems():
			outputDirJob = peak_data.job
			for peak in peak_data.result_peak_ls:
				
				#identify the proper output folder & its creation job
				outputFile = File(os.path.join(outputDirJob.folder, 'peak_%s_loci.h5'%(peak.id)))
				peakLociOutputJob = self.addOutputLociIDOfResultPeakInHDF5Job(workflow, executable=workflow.OutputLociIDOfResultPeakInHDF5, \
														peak_id=peak.id, outputFile=outputFile,\
														parentJobLs=[outputDirJob], extraDependentInputLs=[], \
														transferOutput=False, extraArguments=None, \
														job_max_memory=200)
				no_of_jobs += 1
				peakCorrelationFile = File(os.path.join(outputDirJob.folder, 'maxCorrelationBetweenFirstLociAndPeak%s.h5'%(peak.id)))
				maxLDJob = self.mapMaxLDJobsGivenInputData(workflow, inputData=inputData, datasetName=datasetName, peak_id=peak.id, \
								peakLociH5File=peakLociOutputJob.output, outputFile=peakCorrelationFile, outputDirJob=outputDirJob, \
								chunkSize=chunkSize, parentJobLs=[peakLociOutputJob], extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
								job_max_memory=200)
				no_of_jobs += maxLDJob.no_of_jobs
				
				#final output filename is call_method_biology_category_phenotype_analysis_chr_start_stop_peak_id
				outputFnamePrefix = os.path.join(outputDirJob.folder, 'call_%s_category_%s_phenotype_%s_%s_analysis_%s_chr_%s_%s_%s_%s'%\
												(peak.result.call_method_id, peak.result.phenotype_method.biology_category_id, \
												peak.result.phenotype_method.id, \
												peak.result.phenotype_method.getProperShortName(), peak.result.analysis_method.id,\
												peak.chromosome, peak.start, peak.stop, peak.id))
				outputFile = File("%s.png"%(outputFnamePrefix))
				plotJob = self.addDrawManhattanPlotForLDInHDF5Job(workflow, executable=workflow.DrawManhattanPlotForLDInHDF5, \
										correlationFile=maxLDJob.output, peak_id=peak.id, \
										datasetName=datasetName, outputFile=outputFile,\
										outputFnamePrefix=outputFnamePrefix, parentJobLs=[maxLDJob], \
										extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
					job_max_memory=300)
				no_of_jobs += 1
				if no_of_jobs%2==0:
					sys.stderr.write("%s%s"%("\x08"*len(prevReportedStr), no_of_jobs))
					prevReportedStr = str(no_of_jobs)
		sys.stderr.write("  %s jobs.\n"%(no_of_jobs))
		return no_of_jobs
	
	
	def run(self):
		"""
		2012.3.2
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_250k = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user, password=self.db_passwd, \
									hostname=self.hostname, database=self.dbname)
		db_250k.setup(create_tables=False)
		
		# Create a abstract dag
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = self.initiateWorkflow(workflowName)
		
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		#find all hdf5 correlation files
		inputFnameLs = self.getFilesWithProperSuffixFromFolder(self.inputFolder, suffix='.h5')
		inputData = self.registerAllInputFiles(workflow, inputFnameLs=inputFnameLs, input_site_handler=self.input_site_handler, \
								pegasusFolderName=self.pegasusFolderName)
		#organize final output plots by biology_category, biology_category_id2outputfolder 
		sameCategoryPhenotypeMethodLs = db_250k.getPhenotypeMethodLsGivenBiologyCategoryID(self.biology_category_id, access=self.access)
		sameCategoryPhenotypeMethodIDLs = [pm.id for pm in sameCategoryPhenotypeMethodLs]
		phenotype_method_id_ls = self.phenotype_method_id_ls + sameCategoryPhenotypeMethodIDLs
		result_list = db_250k.getResultLs(call_method_id=self.call_method_id, analysis_method_id_ls=self.analysis_method_id_ls, \
										phenotype_method_id_ls=phenotype_method_id_ls, cnv_method_id=self.cnv_method_id)
		result_id_ls=[result.id for result in result_list]
		sys.stderr.write("%s results.\n"%(len(result_id_ls)))
		result_peak_ls = db_250k.getResultPeakList(result_id_ls=result_id_ls, \
												result_peak_type_id=self.result_peak_type_id)
		
		self.addJobs(workflow, result_peak_ls=result_peak_ls, inputData=inputData, datasetName=self.datasetName, chunkSize=self.chunkSize, \
					pegasusFolderName=self.pegasusFolderName)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)


if __name__ == '__main__':
	main_class = FindGenomeWideLDPatternBetweenSNPsAndPeakWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()