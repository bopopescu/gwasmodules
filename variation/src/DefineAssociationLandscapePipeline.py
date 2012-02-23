#!/usr/bin/env python
"""
Examples:
	# 2011-10-10 call method 80, analysis method 1, min_score 5
	%s  -l 80 -a 1 -f 5 -o call80analysis1MinScore5.xml -u yh -z banyan -c -D condorpool
	
	# 2011-10-16 call method 57 or cnv method 20, analysis method 1, min score 4
	%s -q 57 -A 20 -a 1 -f 4 -o call57cnv20analysis1MinScore4.xml -u yh -z banyan -c -D condorpool
	
	#ditto but analysis method 7 and min score 3
	%s -q 57 -A 20 -a 7 -f 3 -o call57cnv20analysis7MinScore3.xml -u yh -z banyan -c -B condorpool -D condorpool
	#ditto but analysis method 32
	%s -q 57 -A 20 -a 32 -f 3 -o call57cnv20analysis32MinScore3.xml -u yh  -z banyan -c -B condorpool -D condorpool

Description:
	2011-10-12
		output a workflow that runs the DefineAssociationLandscape.py on specified gwas results.
		
		OR among three arguments, call_method_id, call_method_id_ls, cnv_method_id in fetching ResultsMethod.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

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
from AbstractVariationWorkflow import AbstractVariationWorkflow

class DefineAssociationLandscapePipeline(AbstractVariationWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVariationWorkflow.option_default_dict
	option_default_dict.update({
						('result_id_ls', 0, ): [None, 'j', 1, 'comma or dash-separated list of result ids, i.e. 3431-3438,4041'],\
						('call_method_id', 0, int):[None, 'l', 1, 'Restrict results based on this call_method. Default is no such restriction.'],\
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
		
		defineAssociationLandscape = Executable(namespace=namespace, name="DefineAssociationLandscape", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		defineAssociationLandscape.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "DefineAssociationLandscape.py"), site_handler))
		#nucmer.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(defineAssociationLandscape)
		workflow.defineAssociationLandscape = defineAssociationLandscape
	
	def addPeakFindingJob(self, workflow, executable=None, \
							result_id=None, neighbor_distance=None, max_neighbor_distance=None,\
							min_MAF=None, min_score=None, ground_score=None, tax_id=None, \
							commit=0, results_directory=None, logFile=None,\
							parentJobLs=[], job_max_memory=100, job_max_walltime = 60, \
							extraDependentInputLs=[], \
							transferOutput=False, **keywords):
		"""
		2012.2.15
			job_max_walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments("-v", self.drivername, "-z", self.hostname, "-d", self.dbname, \
						"-u", self.db_user, "-p", self.db_passwd,\
						"-j", repr(int(result_id)), "-i", repr(neighbor_distance), \
						"-m", repr(max_neighbor_distance), "-n", repr(min_MAF), "-f", repr(min_score),\
						"-s", repr(ground_score), "-x", repr(tax_id))
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
		
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory, max_walltime=job_max_walltime)
		workflow.addJob(job)
		return job
	
	def run(self):
		"""
		2011-10
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_250k = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user, password=self.db_passwd, \
									hostname=self.hostname, database=self.dbname)
		db_250k.setup(create_tables=False)
		
		pd = PassingData(min_MAF=self.min_MAF,\
					results_directory=self.results_directory, \
					need_chr_pos_ls=0,)
		
		result_query = db_250k.getResultLs(call_method_id=self.call_method_id, analysis_method_id_ls=self.analysis_method_id_ls, \
						phenotype_method_id_ls=self.phenotype_method_id_ls, call_method_id_ls=self.call_method_id_ls,\
						cnv_method_id=self.cnv_method_id)
		result_id_ls = self.result_id_ls
		for result in result_query:
			result_id_ls.append(result.id)
			
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = self.initiateWorkflow(workflowName)
		
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		counter = 0
		for result_id in result_id_ls:
			logFile = File('DefineAssociationPeak_%s.log'%(result_id))
			rm = Stock_250kDB.ResultsMethod.get(result_id)
			inputFile = self.registerOneInputFile(workflow, rm.filename)
			peakFindingJob = self.addPeakFindingJob(workflow, executable=workflow.defineAssociationLandscape, \
							result_id=result_id, neighbor_distance=self.neighbor_distance, max_neighbor_distance=self.max_neighbor_distance,\
							min_MAF=self.min_MAF, min_score=self.min_score, ground_score=self.ground_score, tax_id=self.tax_id, \
							commit=self.commit, results_directory=self.results_directory, logFile=logFile,\
							parentJobLs=[], \
							extraDependentInputLs=[inputFile], \
							transferOutput=True)
			counter += 1
		sys.stderr.write("%s DefineAssociationLandscape jobs.\n"%(counter))
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = DefineAssociationLandscapePipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
