#!/usr/bin/env python
"""
Examples:
	# 2011-10-10 call method 80, analysis method 1, min_score 5
	%s  --call_method_id_ls 80 --analysis_method_id_ls 1 --min_score 5
		-o call80analysis1MinScore5.xml -u yh -z banyan -c -l condorpool
	
	# 2011-10-16 call method 57 or cnv method 20, analysis method 1, min score 4
	%s --call_method_id_ls 57 --cnv_method_id 20 --analysis_method_id_ls 1 --min_score 4
		-o call57cnv20analysis1MinScore4.xml -u yh -z banyan -c -l condorpool
	
	#ditto but analysis method 7 and min score 3 (make sure -c is added, otherwise nothing will be stored in db)
	%s --call_method_id_ls 57 --cnv_method_id 20 --analysis_method_id_ls 7 --min_score 3
		-o call57cnv20analysis7MinScore3.xml -u yh -z banyan -c -j condorpool -l condorpool
	#ditto but analysis method 32
	%s --call_method_id_ls 57 --cnv_method_id 20 --analysis_method_id_ls 32 --min_score 3
		-o call57cnv20analysis32MinScore3.xml -u yh  -z banyan -c -j condorpool -l condorpool

Description:
	2012.11.21 change it to save the landscape data into db
	2011-10-12
		output a workflow that runs the DefineAssociationLandscape.py on specified gwas results.
		
		OR among three arguments, call_method_id, call_method_id_ls, cnv_method_id in fetching ResultsMethod.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])


#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:	   #64bit
#	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
#	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
#else:   #32bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *
from variation.src import Stock_250kDB
from variation.src import AbstractVariationWorkflow

class DefineAssociationLandscapePipeline(AbstractVariationWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVariationWorkflow.option_default_dict.copy()
	my_option_dict = {
					('min_score', 0, float): [4, 'f', 1, 'minimum score to call a peak'],\
					}
	option_default_dict.update({
						('result_id_ls', 0, ): [None, 'w', 1, 'comma or dash-separated list of result ids, i.e. 3431-3438,4041'],\
						('call_method_id_ls', 0, ):[None, 'q', 1, 'Restrict results based on list of call_method_id. Default is no such restriction.'],\
						('cnv_method_id', 0, int):[None, 'A', 1, 'Restrict results based on this cnv_method. Default is no such restriction.'],\
						('analysis_method_id_ls', 0, ):['1,7', 'a', 1, 'Restrict results based on these analysis_methods. coma or dash-separated list'],\
						("phenotype_method_id_ls", 0, ): [None, 'y', 1, 'comma/dash-separated phenotype_method id list, like 1,3-7. Default is all.'],\
						('neighbor_distance', 0, int): [5000, 'i', 1, "within this distance, a locus that increases the association score \
									the fastest is chosen as bridge end. outside this distance, whatever the next point is will be picked."],\
						('max_neighbor_distance', 0, int): [20000, 'm', 1, "beyond this distance, no bridge would be formed."],\
						('ground_score', 0, float): [0, 's', 1, 'minimum score possible in this test'],\
						('min_MAF', 0, float): [0.1, 'n', 1, 'minimum Minor Allele Frequency.'],\
						('tax_id', 1, int): [3702, 'x', 1, 'to get the number of total genes from database, which species.'],\
						})
	option_default_dict.update(my_option_dict)
	def __init__(self,  inputFnameLs=None, **keywords):
		"""
		2011-10
		"""
		AbstractVariationWorkflow.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		
		
		listArgumentName_data_type_ls = [('result_id_ls', int), ("call_method_id_ls", int), \
								("analysis_method_id_ls", int), ('phenotype_method_id_ls', int)]
		listArgumentName2hasContent = ProcessOptions.processListArguments(listArgumentName_data_type_ls, emptyContent=[],\
																class_to_have_attr=self)
		
		self.result_id_ls.sort()
		self.call_method_id_ls.sort()
		self.analysis_method_id_ls.sort()
		self.phenotype_method_id_ls.sort()
	
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2012.11.13
			add more
		2012.2.15
		"""
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		
		DefineAssociationLandscape = Executable(namespace=namespace, name="DefineAssociationLandscape", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		DefineAssociationLandscape.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "association_peak/DefineAssociationLandscape.py"), site_handler))
		executableClusterSizeMultiplierList.append((DefineAssociationLandscape, 0.8))
		
		AssociationLandscape2DB = Executable(namespace=namespace, name="AssociationLandscape2DB", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		AssociationLandscape2DB.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "db/AssociationLandscape2DB.py"), site_handler))
		executableClusterSizeMultiplierList.append((AssociationLandscape2DB, 0.5))
		
		AssociationLandscape2Peak = Executable(namespace=namespace, name="AssociationLandscape2Peak", version=version, \
					os=operatingSystem, arch=architecture, installed=True)
		AssociationLandscape2Peak.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "association_peak/AssociationLandscape2Peak.py"), site_handler))
		executableClusterSizeMultiplierList.append((AssociationLandscape2Peak, 0.2))
		
		AssociationPeak2DB = Executable(namespace=namespace, name="AssociationPeak2DB", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		AssociationPeak2DB.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "db/AssociationPeak2DB.py"), site_handler))
		executableClusterSizeMultiplierList.append((AssociationPeak2DB, 0.5))
		
		AssociationPeak2AssociationLocus = Executable(namespace=namespace, name="AssociationPeak2AssociationLocus", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		AssociationPeak2AssociationLocus.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "association_peak/AssociationPeak2AssociationLocus.py"), site_handler))
		executableClusterSizeMultiplierList.append((AssociationPeak2AssociationLocus, 0))
		
		AssociationLocus2DB = Executable(namespace=namespace, name="AssociationLocus2DB", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		AssociationLocus2DB.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "db/AssociationLocus2DB.py"), site_handler))
		executableClusterSizeMultiplierList.append((AssociationLocus2DB, 0.4))
		
		CountAssociationLocus = Executable(namespace=namespace, name="CountAssociationLocus", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		CountAssociationLocus.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "association_peak/CountAssociationLocus.py"), site_handler))
		executableClusterSizeMultiplierList.append((CountAssociationLocus, 0.8))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
	
	def addDefineLandscapeJob(self, workflow=None, executable=None, \
							result_id=None, neighbor_distance=None, max_neighbor_distance=None,\
							min_MAF=None, tax_id=None, \
							data_dir=None, logFile=None, landscapeOutputFile=None, \
							parentJobLs=None, job_max_memory=100, walltime = 60, sshDBTunnel=None,\
							extraDependentInputLs=None, \
							transferOutput=False, **keywords):
		"""
		2012.11.21 renamed
		2012.2.15
			walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
		"""
		extraArgumentList = ["--result_id %s"%(result_id), "--neighbor_distance %s"%(neighbor_distance), \
						"--max_neighbor_distance %s"%(max_neighbor_distance), "--min_MAF %s"%(min_MAF), "--tax_id %s"%(tax_id)]
		extraOutputLs = []
		if data_dir:
			extraArgumentList.extend(["--data_dir", data_dir])
		if logFile:
			extraArgumentList.extend(["--logFilename", logFile])
			extraOutputLs.append(logFile)
		key2ObjectForJob = {}
		
		job = self.addGenericDBJob(workflow=workflow, executable=executable, inputFile=None, \
					outputFile=landscapeOutputFile, \
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, extraOutputLs=extraOutputLs, \
					transferOutput=transferOutput, \
					extraArguments=None, extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, \
					sshDBTunnel=sshDBTunnel, walltime=walltime,\
					key2ObjectForJob=key2ObjectForJob, objectWithDBArguments=self, **keywords)
		return job
	
	def addAssociationLandscape2PeakJob(self, workflow=None, executable=None, \
				inputFile=None, outputFile=None,  min_score=None, ground_score=None, \
				parentJobLs=None, job_max_memory=100, walltime = 60, \
				extraDependentInputLs=None, extraArguments=None, \
				transferOutput=False, **keywords):
		"""
		2012.11.21
			walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
		"""
		extraArgumentList = []
		key2ObjectForJob = {}
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		if min_score is not None:
			extraArgumentList.append("--min_score %s"%(min_score))
		if ground_score is not None:
			extraArgumentList.append("--ground_score %s"%(ground_score))
		
		job = self.addGenericJob(workflow=workflow, executable=executable, inputFile=inputFile, \
					outputFile=outputFile, \
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
					transferOutput=transferOutput, \
					extraArguments=extraArguments, extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, \
					walltime=walltime,\
					key2ObjectForJob=key2ObjectForJob, **keywords)
		return job
	
	def addAssociationPeak2LocusJob(self, workflow=None, executable=None, \
				inputFile=None, outputFile=None, min_overlap_ratio=None, \
				parentJobLs=None, job_max_memory=100, walltime = 60, \
				extraDependentInputLs=None, \
				transferOutput=False, **keywords):
		"""
		2012.11.21
			walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
		"""
		extraArgumentList = []
		if min_overlap_ratio is not None:
			extraArgumentList.append("--min_overlap_ratio %s"%(min_overlap_ratio))
		job = self.addGenericJob(workflow=workflow, executable=executable, inputFile=inputFile, \
					outputFile=outputFile, \
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
					transferOutput=transferOutput, \
					extraArguments=None, extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, \
					walltime=walltime,\
					key2ObjectForJob=None, **keywords)
		return job
	
	def addCountAssociationLocusJob(self, workflow=None, executable=None, \
				inputFileList=None, inputFile=None, outputFile=None, \
				parentJobLs=None, job_max_memory=100, walltime = 60, \
				extraDependentInputLs=None, \
				transferOutput=False, **keywords):
		"""
		2012.11.21
			walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
		"""
		extraArgumentList = []
		job = self.addAbstractMatrixFileWalkerJob(workflow=workflow, executable=executable, inputFileList=inputFileList, \
								inputFile=inputFile, outputFile=outputFile, \
					outputFnamePrefix=None, whichColumn=None, whichColumnHeader=None, \
					logY=False, valueForNonPositiveYValue=None, \
					minNoOfTotal=None,\
					samplingRate=None, \
					inputFileFormat=2, outputFileFormat=2,\
					parentJobLs=parentJobLs, \
					extraDependentInputLs=extraDependentInputLs, extraArgumentList=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=job_max_memory, \
					walltime=walltime,\
					sshDBTunnel=False, \
					objectWithDBArguments=None, **keywords)
		return job
	
	def addAssociationLandscape2DBJob(self, executable=None, inputFile=None, result_id=None, \
					outputFile=None, data_dir=None, logFile=None, commit=False,\
					neighbor_distance=None, max_neighbor_distance=None, min_MAF=None, \
					parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					extraArguments=None, job_max_memory=2000, sshDBTunnel=None, **keywords):
		"""
		2012.11.13
		"""
		extraArgumentList = ['--neighbor_distance %s'%(neighbor_distance),\
							'--max_neighbor_distance %s'%(max_neighbor_distance), \
							'--min_MAF %s'%(min_MAF)]
		if result_id:
			extraArgumentList.append('--result_id %s'%(result_id))
		extraOutputLs = []
		
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		
		job = self.addGenericFile2DBJob(executable=executable, inputFile=inputFile, outputFile=outputFile,\
					data_dir=data_dir, logFile=logFile, commit=commit,\
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
					extraOutputLs=extraOutputLs, transferOutput=transferOutput, \
					extraArgumentList=extraArgumentList, extraArguments=extraArguments, \
					job_max_memory=job_max_memory, sshDBTunnel=sshDBTunnel, objectWithDBArguments=self,\
					**keywords)
		return job
	
	def run(self):
		"""
		2011-10
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		pd = PassingData(min_MAF=self.min_MAF,\
					data_dir=self.data_dir, \
					need_chr_pos_ls=0,)
		
		result_query = self.db_250k.getResultLs(analysis_method_id_ls=self.analysis_method_id_ls, \
						phenotype_method_id_ls=self.phenotype_method_id_ls, call_method_id_ls=self.call_method_id_ls,\
						cnv_method_id=self.cnv_method_id)
		result_id_ls = self.result_id_ls
		for result in result_query:
			result_id_ls.append(result.id)
			
		workflow = self.initiateWorkflow()
		
		self.registerExecutables()
		self.registerCustomExecutables(workflow)
		
		counter = 0
		
		topOutputDir = "%sAssociationLandscape"%(self.pegasusFolderName)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow=self, mkdir=self.mkdirWrap, outputDir=topOutputDir)
		
		resultLandscapeType = self.db_250k.getResultLandscapeType(min_MAF=self.min_MAF, \
									neighbor_distance=self.neighbor_distance, \
									max_neighbor_distance=self.max_neighbor_distance)
		
		for result_id in result_id_ls:
			result = Stock_250kDB.ResultsMethod.get(result_id)
			
			associationResultFile = self.registerOneInputFile(inputFname=result.getFileAbsPath(oldDataDir=self.db_250k.data_dir, newDataDir=self.data_dir), \
												folderName=self.pegasusFolderName)
			logFile = File(os.path.join(topOutputDirJob.output, 'Result%s_LandscapeType%s.log'%\
									(result_id, resultLandscapeType.id)))
			landscapeOutputFile = File(os.path.join(topOutputDirJob.output, 'Result%s_LandscapeType%s.h5'%\
									(result_id, resultLandscapeType.id)))
			
			defineLandscapeJob = self.addDefineLandscapeJob(workflow, executable=workflow.DefineAssociationLandscape, \
							result_id=result_id, neighbor_distance=self.neighbor_distance, \
							max_neighbor_distance=self.max_neighbor_distance,\
							min_MAF=self.min_MAF, tax_id=self.tax_id, \
							data_dir=self.data_dir, logFile=logFile,\
							landscapeOutputFile=landscapeOutputFile,\
							extraDependentInputLs=[associationResultFile], \
							parentJobLs=[topOutputDirJob], sshDBTunnel=self.needSSHDBTunnel,\
							transferOutput=False)
			
			logFile = File(os.path.join(topOutputDirJob.output, 'Result%s_LandscapeType%s_log.tsv'%\
										(result_id, resultLandscapeType.id)))
			landscape2DBJob = self.addAssociationLandscape2DBJob(executable=self.AssociationLandscape2DB, inputFile=defineLandscapeJob.output, \
						result_id=result_id, \
						data_dir=self.data_dir, logFile=logFile, commit=self.commit, \
						min_MAF=self.min_MAF, \
						neighbor_distance=self.neighbor_distance, max_neighbor_distance=self.max_neighbor_distance, \
						parentJobLs=[topOutputDirJob, defineLandscapeJob], \
						extraDependentInputLs=None, transferOutput=True, extraArguments=None, job_max_memory=1000, sshDBTunnel=self.needSSHDBTunnel)
			
			#add landscape -> peak job
			outputFile = File('%s_peak.h5'%(outputFnamePrefix))
			self.addAssociationLandscape2PeakJob(executable=self.AssociationLandscape2Peak, \
				inputFile=defineLandscapeJob.output, outputFile=outputFile, min_score=min_score, ground_score=ground_score, \
				data_dir=data_dir, \
				parentJobLs=[defineLandscapeJob], job_max_memory=100, walltime = 60, \
				extraDependentInputLs=None, \
				transferOutput=False)
			
			counter += 1
		sys.stderr.write("%s total jobs.\n"%(self.no_of_jobs))
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		self.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = DefineAssociationLandscapePipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
