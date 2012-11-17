#!/usr/bin/env python
"""
Examples:
	# 2011-10-10 call method 80, analysis method 1, min_score 5
	%s  -E 80 -a 1 -f 5 -o call80analysis1MinScore5.xml -u yh -z banyan -c -l condorpool
	
	# 2011-10-16 call method 57 or cnv method 20, analysis method 1, min score 4
	%s -q 57 -A 20 -a 1 -f 4 -o call57cnv20analysis1MinScore4.xml -u yh -z banyan -c -l condorpool
	
	#ditto but analysis method 7 and min score 3 (make sure -c is added, otherwise nothing will be stored in db)
	%s -q 57 -A 20 -a 7 -f 3 -o call57cnv20analysis7MinScore3.xml -u yh -z banyan -c -j condorpool -l condorpool
	#ditto but analysis method 32
	%s -q 57 -A 20 -a 32 -f 3 -o call57cnv20analysis32MinScore3.xml -u yh  -z banyan -c -j condorpool -l condorpool

Description:
	2012.11.21 change it to save the landscape data into db
	2011-10-12
		output a workflow that runs the DefineAssociationLandscape.py on specified gwas results.
		
		OR among three arguments, call_method_id, call_method_id_ls, cnv_method_id in fetching ResultsMethod.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

from sqlalchemy.types import LargeBinary

#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:	   #64bit
#	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
#	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
#else:   #32bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *
import Stock_250kDB
from AbstractVariationWorkflow import AbstractVariationWorkflow

class DefineAssociationLandscapePipeline(AbstractVariationWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVariationWorkflow.option_default_dict
	option_default_dict.update({
						('result_id_ls', 0, ): [None, 'w', 1, 'comma or dash-separated list of result ids, i.e. 3431-3438,4041'],\
						('call_method_id', 0, int):[None, 'E', 1, 'Restrict results based on this call_method. Default is no such restriction.'],\
						('call_method_id_ls', 0, ):[None, 'q', 1, 'Restrict results based on list of call_method_id. Default is no such restriction.'],\
						('cnv_method_id', 0, int):[None, 'A', 1, 'Restrict results based on this cnv_method. Default is no such restriction.'],\
						('analysis_method_id_ls', 0, ):['1,7', 'a', 1, 'Restrict results based on these analysis_methods. coma or dash-separated list'],\
						("phenotype_method_id_ls", 0, ): [None, 'y', 1, 'comma/dash-separated phenotype_method id list, like 1,3-7. Default is all.'],\
						('neighbor_distance', 0, int): [5000, 'i', 1, "within this distance, a locus that increases the association score \
									the fastest is chosen as bridge end. outside this distance, whatever the next point is will be picked."],\
						('max_neighbor_distance', 0, int): [20000, 'm', 1, "beyond this distance, no bridge would be formed."],\
						('min_MAF', 0, float): [0.1, 'n', 1, 'minimum Minor Allele Frequency.'],\
						('min_score', 0, float): [4, 'f', 1, 'minimum score to call a peak'],\
						('ground_score', 0, float): [0, 's', 1, 'minimum score possible in this test'],\
						('tax_id', 1, int): [3702, 'x', 1, 'to get the number of total genes from database, which species.'],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction (individual_alignment and/or individual_alignment.path'],\
						
						})
	
	def __init__(self,  **keywords):
		"""
		2011-10
		"""
		AbstractVariationWorkflow.__init__(self, **keywords)
		
		if getattr(self, 'result_id_ls', None):
			self.result_id_ls = getListOutOfStr(self.result_id_ls, data_type=int)
			self.result_id_ls.sort()
		else:
			self.result_id_ls = []
		
		if getattr(self, 'call_method_id_ls', None):
			self.call_method_id_ls = getListOutOfStr(self.call_method_id_ls, data_type=int)
			self.call_method_id_ls.sort()
		else:
			self.call_method_id_ls = []
		
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
		
		ResultLandscape2DB = Executable(namespace=namespace, name="ResultLandscape2DB", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		ResultLandscape2DB.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "db/ResultLandscape2DB.py"), site_handler))
		executableClusterSizeMultiplierList.append((ResultLandscape2DB, 0.5))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
	
	def addPeakFindingJob(self, workflow=None, executable=None, \
							result_id=None, neighbor_distance=None, max_neighbor_distance=None,\
							min_MAF=None, min_score=None, ground_score=None, tax_id=None, \
							commit=0, data_dir=None, logFile=None, landscapeOutputFile=None,\
							parentJobLs=None, job_max_memory=100, job_max_walltime = 60, sshDBTunnel=None,\
							extraDependentInputLs=None, \
							transferOutput=False, **keywords):
		"""
		2012.2.15
			job_max_walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
		"""
		extraArgumentList = ["--result_id_ls", repr(int(result_id)), "--neighbor_distance", repr(neighbor_distance), \
						"--max_neighbor_distance", repr(max_neighbor_distance), "--min_MAF", repr(min_MAF), "--min_score", repr(min_score),\
						"--ground_score", repr(ground_score), "--tax_id", repr(tax_id)]
		extraOutputLs = []
		if commit:
			extraArgumentList.append("--commit")
		if data_dir:
			extraArgumentList.extend(["--data_dir", data_dir])
		if logFile:
			extraArgumentList.extend(["--logFilename", logFile])
			extraOutputLs.append(logFile)
		if landscapeOutputFile:
			extraArgumentList.extend(["--landscapeLocusIDFname", landscapeOutputFile])
			extraOutputLs.append(landscapeOutputFile)
		job = self.addGenericDBJob(workflow=workflow, executable=executable, inputFile=None, \
					outputFile=None, \
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, extraOutputLs=extraOutputLs, \
					transferOutput=transferOutput, \
					extraArguments=None, extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, \
					sshDBTunnel=sshDBTunnel, max_walltime=job_max_walltime,\
					key2ObjectForJob=None, objectWithDBArguments=self, **keywords)
		return job
	
	def addResultLandscape2DBJob(self, executable=None, inputFile=None, result_id=None, \
					call_method_id=None, phenotype_method_id=None, \
					analysis_method_id=None, results_method_type_id=1, data_dir=None, \
					neighbor_distance=None, max_neighbor_distance=None, min_MAF=None, \
					landscapeLocusIDFile=None, logFile=None, outputFile=None, commit=False,\
					parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					extraArguments=None, job_max_memory=2000, sshDBTunnel=None, **keywords):
		"""
		2012.11.13
		"""
		extraArgumentList = ['--result_id %s'%(result_id), \
							'--phenotype_method_id %s'%(phenotype_method_id), \
							'--analysis_method_id %s'%(analysis_method_id), \
							'--results_method_type_id %s'%(results_method_type_id), \
							'--neighbor_distance %s'%(neighbor_distance),\
							'--max_neighbor_distance %s'%(max_neighbor_distance), \
							'--min_MAF %s'%(min_MAF)]
		extraOutputLs = []
		if call_method_id:
			extraArgumentList.append('--call_method_id %s'%(call_method_id))
		if data_dir:
			extraArgumentList.append('--data_dir %s'%(data_dir))
		if commit:
			extraArgumentList.append('--commit')
		if logFile:
			extraArgumentList.extend(["--logFilename", logFile])
			extraOutputLs.append(logFile)
		
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		if landscapeLocusIDFile:
			extraArgumentList.extend(['--landscapeLocusIDFname', landscapeLocusIDFile])
			extraDependentInputLs.append(landscapeLocusIDFile)
		
		job = self.addGenericDBJob(executable=executable, inputFile=inputFile, outputFile=outputFile,\
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=extraOutputLs, transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, extraArguments=extraArguments, \
						job_max_memory=job_max_memory, sshDBTunnel=sshDBTunnel,\
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
		
		result_query = self.db.getResultLs(call_method_id=self.call_method_id, analysis_method_id_ls=self.analysis_method_id_ls, \
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
		
		result_landscape_type = self.db.getResultLandscapeType(min_MAF=self.min_MAF, \
									neighbor_distance=self.neighbor_distance, \
									max_neighbor_distance=self.max_neighbor_distance)
		
		for result_id in result_id_ls:
			result = Stock_250kDB.ResultsMethod.get(result_id)
			
			inputFile = self.registerOneInputFile(inputFname=result.getFilePath(oldDataDir=self.db.data_dir, newDataDir=self.data_dir), \
												folderName=self.pegasusFolderName)
			logFile = File(os.path.join(topOutputDirJob.output, 'Result%s_LandscapeType%s.log'%\
									(result_id, result_landscape_type.id)))
			landscapeOutputFile = File(os.path.join(topOutputDirJob.output, 'Result%s_LandscapeType%s.tsv'%\
									(result_id, result_landscape_type.id)))
			
			peakFindingJob = self.addPeakFindingJob(workflow, executable=workflow.DefineAssociationLandscape, \
							result_id=result_id, neighbor_distance=self.neighbor_distance, \
							max_neighbor_distance=self.max_neighbor_distance,\
							min_MAF=self.min_MAF, min_score=self.min_score, ground_score=self.ground_score, tax_id=self.tax_id, \
							commit=self.commit, data_dir=self.data_dir, logFile=logFile,\
							landscapeOutputFile=landscapeOutputFile,\
							extraDependentInputLs=[inputFile], \
							parentJobLs=[topOutputDirJob], sshDBTunnel=self.needSSHDBTunnel,\
							transferOutput=False)
			
			outputFile = File(os.path.join(topOutputDirJob.output, 'Result%s_LandscapeType%s_ReducedGWAS.tsv'%\
										(result_id, result_landscape_type.id)))
			landscape2DBJob = self.addResultLandscape2DBJob(executable=self.ResultLandscape2DB, inputFile=inputFile, \
						result_id=result_id, call_method_id=result.call_method_id, \
						phenotype_method_id=result.phenotype_method_id, analysis_method_id=result.analysis_method_id, \
						results_method_type_id=result.results_method_type_id, \
						data_dir=self.data_dir,\
						neighbor_distance=self.neighbor_distance, max_neighbor_distance=self.max_neighbor_distance, \
						min_MAF=self.min_MAF, landscapeLocusIDFile=landscapeOutputFile, \
						outputFile=outputFile, commit=self.commit, parentJobLs=[topOutputDirJob, peakFindingJob], \
						extraDependentInputLs=None, transferOutput=True, job_max_memory=1000, sshDBTunnel=self.needSSHDBTunnel)
			
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
