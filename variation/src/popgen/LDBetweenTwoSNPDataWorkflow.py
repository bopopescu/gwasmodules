#!/usr/bin/env python
"""
Examples:
	# 2012.2.28
	%s -i /Network/Data/250k/db/dataset/call_method_32.tsv -j /Network/Data/250k/db/dataset/call_method_80.tsv
		-o LD_between_call_32_and_80.xml -F LD -B condorpool -D condorpool -u yh -z banyan
	
	%s 
	
	%s
	
Description:
	2012.2.28
		output a workflow that calculates LD between two Yu-format (Strain X Locus) SNP datasets.
		
		0. figure out how many columns each dataset has
		1. (defunct, superceded by step 2) remove non-intersecting rows via TwoSNPData.py for each dataset
			~/script/pymodule/TwoSNPData.py -i /Network/Data/250k/db/dataset/call_method_57.tsv
				-j /Network/Data/250k/db/dataset/call_method_32.tsv
				-o /tmp/call_method_57_only_accessions_in_call_32.tsv -c0 -m0 -p3 -w2
		
		2. order one dataset's rows in the same order as the other
				~/script/pymodule/pegasus/mapper/Order2ndSNPDataRowsSameAs1stSNPData.py -i /Network/Data/250k/db/dataset/call_method_32.tsv
					-j /Network/Data/250k/db/dataset/call_method_80.2011.5.2.same.as.57.tsv -o /tmp/call_method_57_in_32_order.tsv  -m1
			
		3. convert each dataset into hdf5 format with min-MAF filter
			~/script/pymodule/pegasus/mapper/ConvertSNPData2HDF5.py
		4. call a C++ program to calculate correlation between two datasets given column-index-span
			~/script/pymodule/pegasus/mapper/CalculateColCorBetweenTwoHDF5 -o /tmp/out
				-i /tmp/call_80.hdf5 -j /tmp/call_80.hdf5 -s 0 -t 5000 -u 20000 -v 25885 -c 0.4
		

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], )

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO, csv
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, figureOutDelimiter
from variation.src import Stock_250kDB
from variation.src import AbstractVariationWorkflow

class LDBetweenTwoSNPDataWorkflow(AbstractVariationWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVariationWorkflow.option_default_dict.copy()
	common_option_dict = {
						("input1Fname", 1, ): [None, 'i', 1, '1st input dataset'],\
						('input2Fname', 1, ): [None, '', 1, '2nd input dataset'],\
						('chunkSize', 0, int): [5000, '', 1, 'Loci from each will be partitioned into chunks of 5000'],\
						('min_MAF', 1, float): [0.1, '', 1, 'minimum minor allele frequency to filter the input dataset'],\
						('min_cor', 1, float): [0.1, '', 1, 'minimum correlation for output'],\
						}
	
	option_default_dict.update(common_option_dict)
	
	def __init__(self,  **keywords):
		"""
		2012.2.28
		"""
		AbstractVariationWorkflow.__init__(self, **keywords)
	
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
		
		Order2ndSNPDataRowsSameAs1stSNPData = Executable(namespace=namespace, name="Order2ndSNPDataRowsSameAs1stSNPData", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		Order2ndSNPDataRowsSameAs1stSNPData.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/Order2ndSNPDataRowsSameAs1stSNPData.py"), \
													site_handler))
		#twoGWASPeakOverlap.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(Order2ndSNPDataRowsSameAs1stSNPData)
		workflow.Order2ndSNPDataRowsSameAs1stSNPData = Order2ndSNPDataRowsSameAs1stSNPData
		
		CalculateColCorBetweenTwoHDF5 = Executable(namespace=namespace, name="CalculateColCorBetweenTwoHDF5", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		CalculateColCorBetweenTwoHDF5.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/CalculateColCorBetweenTwoHDF5"), \
									site_handler))
		#CalculateColCorBetweenTwoHDF5.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(CalculateColCorBetweenTwoHDF5)
		workflow.CalculateColCorBetweenTwoHDF5 = CalculateColCorBetweenTwoHDF5
	
		ConvertSNPData2HDF5 = Executable(namespace=namespace, name="ConvertSNPData2HDF5", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		ConvertSNPData2HDF5.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/ConvertSNPData2HDF5.py"), \
									site_handler))
		#ConvertSNPData2HDF5.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(ConvertSNPData2HDF5)
		workflow.ConvertSNPData2HDF5 = ConvertSNPData2HDF5
		
		
		MergeTwoLocusCorrelationHDF5 = Executable(namespace=namespace, name="MergeTwoLocusCorrelationHDF5", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		MergeTwoLocusCorrelationHDF5.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/MergeTwoLocusCorrelationHDF5.py"), \
									site_handler))
		#MergeTwoLocusCorrelationHDF5.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(MergeTwoLocusCorrelationHDF5)
		workflow.MergeTwoLocusCorrelationHDF5 = MergeTwoLocusCorrelationHDF5
			
	def addCalculateColCorBetweenTwoHDF5Job(self, workflow, executable=None, inputFile1=None, inputFile2=None, outputFile=None, \
					i1_start=None, i1_stop=None, i2_start=None, i2_stop=None, min_cor=None, \
					parentJobLs=[], extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
					job_max_memory=100, **keywords):
		"""
		2012.3.2
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments('-o', outputFile,'-i', inputFile1, '-j', inputFile2, '-s %s'%i1_start, '-t %s'%i1_stop, \
						'-u %s'%i2_start, '-v %s'%i2_stop, '-c %s'%min_cor)
		job.uses(inputFile1, transfer=True, register=True, link=Link.INPUT)
		job.uses(inputFile2, transfer=True, register=True, link=Link.INPUT)
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
	
	
	
	def getNoOfLociFromSNPData(self, inputFname):
		"""
		2012.3.2
		"""
		sys.stderr.write("Getting no of loci from %s ..."%(os.path.basename(inputFname)))
		reader = csv.reader(open(inputFname), delimiter=figureOutDelimiter(inputFname))
		header = reader.next()
		first_data_row = reader.next()
		no_of_cols = len(first_data_row) - 2
		del reader
		sys.stderr.write("%s columns.\n"%(no_of_cols))
		return no_of_cols
	
	def addOrderDatasetRowJob(self, workflow, executable=None, inputFile1=None, inputFile2=None, outputFile=None,\
					parentJobLs=[], extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
					job_max_memory=100, **keywords):
		"""
		2012.3.2
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments('-i', inputFile1, '-j', inputFile2, '-o', outputFile, '-m1')
		job.uses(inputFile1, transfer=True, register=True, link=Link.INPUT)
		job.uses(inputFile2, transfer=True, register=True, link=Link.INPUT)
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
	
	def addConvertSNPData2HDF5Job(self, workflow, executable=None, inputFile=None, outputFile=None, min_MAF=None, \
					parentJobLs=[], extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
					job_max_memory=100, **keywords):
		"""
		2012.3.2
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments('-i', inputFile, '-o', outputFile, '-n %s'%(min_MAF))
		job.uses(inputFile, transfer=True, register=True, link=Link.INPUT)
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
	
	def addJobs(self, workflow, inputData=None, min_MAF=None, min_cor=None, chunkSize=None, pegasusFolderName=""):
		"""
		2012.3.3
		"""
		
		sys.stderr.write("Adding LD-calculating jobs on %s input datasets ..."%(len(inputData.jobDataLs)))
		returnJobData = PassingData()
		
		no_of_jobs = 0
		
		topOutputDir = pegasusFolderName
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		
		finalCorrelationOutputFile = File(os.path.join(topOutputDir, 'correlation.h5'))
		
		correlationMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.MergeTwoLocusCorrelationHDF5, \
						outputF=finalCorrelationOutputFile, transferOutput=True, extraArguments='-d correlation', parentJobLs=[topOutputDirJob])
		
		inputJobData1, inputJobData2 = inputData.jobDataLs[:2]
		inputFile1 = inputJobData1.output
		inputFile2 = inputJobData2.output
		outputFile = File(os.path.join(topOutputDir, 'input1_in_input2_order.tsv'))
		orderDatasetRowJob1 = self.addOrderDatasetRowJob(workflow, executable=workflow.Order2ndSNPDataRowsSameAs1stSNPData, inputFile1=inputFile2, \
								inputFile2=inputFile1, outputFile=outputFile,\
								parentJobLs=[topOutputDirJob]+inputJobData1.jobLs, extraDependentInputLs=[], \
								transferOutput=False, extraArguments=None, \
								job_max_memory=1000)
		outputFile = File(os.path.join(topOutputDir, 'input2_in_same_order.tsv'))
		orderDatasetRowJob2 = self.addOrderDatasetRowJob(workflow, executable=workflow.Order2ndSNPDataRowsSameAs1stSNPData, inputFile1=orderDatasetRowJob1.output, \
								inputFile2=inputFile2, outputFile=outputFile,\
								parentJobLs=[orderDatasetRowJob1] + inputJobData2.jobLs, extraDependentInputLs=[], transferOutput=False, extraArguments=None, \
								job_max_memory=1000)
		
		outputFile = File(os.path.join(topOutputDir, 'input1.hdf5'))
		convertDataset2HDF5Job1 = self.addConvertSNPData2HDF5Job(workflow, executable=workflow.ConvertSNPData2HDF5, \
																inputFile=orderDatasetRowJob1.output, \
								outputFile=outputFile, min_MAF=min_MAF, \
								parentJobLs=[orderDatasetRowJob1], extraDependentInputLs=[], transferOutput=False, extraArguments=None, \
								job_max_memory=100)
		
		outputFile = File(os.path.join(topOutputDir, 'input2.hdf5'))
		convertDataset2HDF5Job2 = self.addConvertSNPData2HDF5Job(workflow, executable=workflow.ConvertSNPData2HDF5, \
																inputFile=orderDatasetRowJob2.output, \
									outputFile=outputFile, min_MAF=min_MAF, \
									parentJobLs=[orderDatasetRowJob2], extraDependentInputLs=[], transferOutput=False, extraArguments=None, \
									job_max_memory=100)
		no_of_jobs += 5
		
		no_of_cols_input1 = self.getNoOfLociFromSNPData(inputFile1.abspath)
		no_of_cols_input2 = self.getNoOfLociFromSNPData(inputFile2.abspath)
		for i1_start in range(0, no_of_cols_input1, chunkSize):
			i1_stop = min(i1_start + chunkSize -1, no_of_cols_input1-1)
			for i2_start in range(0, no_of_cols_input2, chunkSize):
				i2_stop = min(i2_start + chunkSize -1, no_of_cols_input2-1)
				outputFile = os.path.join(topOutputDir, 'cor_i1_%s_%s_i2_%s_%s.h5'%(i1_start, i1_stop, i2_start, i2_stop))
				corCalulationJob = self.addCalculateColCorBetweenTwoHDF5Job(workflow, executable=workflow.CalculateColCorBetweenTwoHDF5, \
												inputFile1=convertDataset2HDF5Job1.output, inputFile2=convertDataset2HDF5Job2.output, \
												outputFile=outputFile, i1_start=i1_start, i1_stop=i1_stop, i2_start=i2_start, i2_stop=i2_stop,\
												min_cor=min_cor, parentJobLs=[convertDataset2HDF5Job1, convertDataset2HDF5Job2], \
												extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
												job_max_memory=50)
				no_of_jobs += 1
				self.addInputToStatMergeJob(workflow, statMergeJob=correlationMergeJob, \
								inputF=corCalulationJob.output, parentJobLs=[corCalulationJob])
		
		sys.stderr.write("%s jobs.\n"%(no_of_jobs))
		return correlationMergeJob
	
	def run(self):
		"""
		2012.3.2
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_250k = self.db_250k
		
		# Create a abstract dag
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = self.initiateWorkflow(workflowName)
		
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		inputData = self.registerAllInputFiles(workflow, inputFnameLs=[self.input1Fname, self.input2Fname], input_site_handler=self.input_site_handler, \
								pegasusFolderName=self.pegasusFolderName)
		self.addJobs(workflow, inputData=inputData, min_MAF=self.min_MAF, min_cor=self.min_cor, chunkSize=self.chunkSize, \
					pegasusFolderName=self.pegasusFolderName)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)


if __name__ == '__main__':
	main_class = LDBetweenTwoSNPDataWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()