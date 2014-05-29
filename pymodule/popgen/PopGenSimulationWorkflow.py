#!/usr/bin/env python
"""
Examples:
	# scaffold (120) is used as reference in alignment. get genome sequence id from 1 to 8.
	%s -o workflow.xml -a 120 --ind_seq_id_ls 1-8 -y2 -s2
	
	# 1Mb-BAC (9) is used as reference.
	%s -o workflow.xml -a 9 --ind_seq_id_ls 1-4 -y2 -s2
	
	# 8 genomes versus top 156 contigs
	%s -o workflow_8GenomeVsTop156Contigs.xml -u yh  -a 128 --ind_seq_id_ls 1-8 -N 156 -y2 -s2
	
	
Description:
	2013.2.26
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])


#bit_number = math.log(sys.maxint)/math.log(2)
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from Pegasus.DAX3 import *
from pymodule import ProcessOptions, PassingData, AbstractWorkflow, utils
from pymodule.pegasus import yh_pegasus

parentClass = AbstractWorkflow
class PopGenSimulationWorkflow(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	#option_default_dict.pop(('inputDir', 0, ))
	
	
	option_default_dict.update({
						#('seqCoverageFname', 0, ): ['', 'q', 1, 'The sequence coverage file. tab/comma-delimited: individual_sequence.id coverage'],\
						("recombinationRate", 1, float): [1e-9, '', 1, 'recombination rate per base per meiosis.'],\
						("mutationRate", 1, float): [1e-8, '', 1, 'mutation rate per base per generation.'],\
						
						("initialEffectivePopulationSize", 1, int ): [50000, '', 1, 'the effective population size for pop-gen simulation to start with'],\
						("otherParametersPassedToPopGenSimulator", 0, ): [None, '', 1, ''],\
						
						("sampleSize", 1, int ): [50, '', 1, 'how many individuals from the pop-gen simulation to be sampled out to constitute the final output.\n\
			measured in terms of the number of chromosomes (haploid).'],\
						("noOfLociToSimulate", 1, int): [1, '', 1, 'how many loci to simulate'],\
						("simulateLocusLengthList", 1, ): ['10000000', '', 1, 'a coma-separated list of locus length to simulate. if not a list, assuming just one locus'],\
						("sfs_code_path", 1, ): ['%s/script/repository/sfscode/bin/sfs_code', '', 1, 'path to the sfs_code pop-gen forward simulator'],\
						("slim_path", 1, ): ['%s/script/repository/slim/slim', '', 1, 'path to the slim pop-gen forward simulator'],\
						("msHOT_lite_path", 1, ): ['%s/script/lh3_foreign/msHOT-lite/msHOT-lite', '', 1, 'path to the ms pop-gen backwards simulator'],\
						
						
						("simulateIndels", 0, int): [0, '', 1, 'toggle to simulate insertion/deletion polymorphism as well in pop-gen simulation'],\
						("noOfReplicates", 1, int): [5, '', 1, 'how many replicates for this simulation setup'],\
						
						})

	def __init__(self,  **keywords):
		"""
		2012.9.17 call parentClass.__init__() directly
		2011-7-11
		"""
		self.pathToInsertHomePathList.extend(['sfs_code_path'])
		
		parentClass.__init__(self, **keywords)
		
		if hasattr(self, 'simulateLocusLengthList', None):
			self.simulateLocusLengthList = utils.getListOutOfStr(self.simulateLocusLengthList, data_type=int)
		else:
			self.simulateLocusLengthList = []
	
	def addMSSimulationJob(self, executableFile=None, outputFile=None, \
					recombinationRate=None, mutationRate=None, \
					initialEffectivePopulationSize=50000,\
					otherParametersPassedToPopGenSimulator="",\
					sampleSize=50, noOfLociToSimulate=1,\
					simulateLocusLengthList=None, lh3SpecificOutput=True, outputGeneTree=False, \
					parentJobLs=None, extraDependentInputLs=None, extraOutputLs=None, transferOutput=False, \
					extraArguments=None, extraArgumentList=None,\
					job_max_memory=2000, walltime=120, **keywords):
		"""
		2013.08.04
			sampleSize is measured in terms of the number of chromosomes (haploid)
			
			rho (Effective population recombination rate) =  4N*r*simulateLocusLength
			theta (Effective population mutation rate) = 4N*mu*simulateLocusLength
		
		"""
		if executableFile is None:
			executableFile = self.msHOT_liteExecutableFile
		if extraArgumentList is None:
			extraArgumentList = []
		if simulateLocusLengthList:
			locusLength = simulateLocusLengthList[0]
		else:
			locusLength = None
		if recombinationRate and locusLength and initialEffectivePopulationSize:
			rho = 4 * initialEffectivePopulationSize * recombinationRate * locusLength
		else:
			rho=None
		if mutationRate and locusLength and initialEffectivePopulationSize:
			theta = 4 * initialEffectivePopulationSize * mutationRate * locusLength
		else:
			theta=None
		
		noOfSamplings = 1	#the number of independent samples(output) from simulation
		extraArgumentList.extend([repr(sampleSize), noOfSamplings])
		if theta:
			extraArgumentList.append("-t %s"%(theta))
		else:
			sys.stderr.write("ERROR: theta (%s) is not valid, can not add ms simulation job.\n"%(theta))
			sys.exit(4)
		
		if outputGeneTree:
			extraArgumentList.append("-T")	#-T is to output the gene tree
		
		
		if rho is not None and locusLength is not None:
			extraArgumentList.extend(["-r %s %s"%(rho, locusLength)])
		if initialEffectivePopulationSize is not None:
			extraArgumentList.append()
		if otherParametersPassedToPopGenSimulator:
			extraArgumentList.append(otherParametersPassedToPopGenSimulator)
		if lh3SpecificOutput:	#only for msHOT-lite from lh3_foreign
			extraArgumentList.append("-l")
		commandline = " ".join(extraArgumentList)
		job = self.addGenericPipeCommandOutput2FileJob(executable=self.msShellPipe, \
					executableFile=executableFile, \
					outputFile=outputFile, \
					parentJobLs=parentJobLs, extraDependentInputLs=None,\
					extraOutputLs=extraOutputLs, transferOutput=transferOutput, \
					extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
					job_max_memory=job_max_memory, walltime=walltime, **keywords)
		job.rho = rho
		job.theta = theta
		job.locusLength = locusLength
		job.commandline = commandline
		return job
	
	
	def addAllJobs(self, workflow=None, \
				data_dir=None, \
				outputDirPrefix="", transferOutput=True, **keywords):
		"""
		2013.2.27
			run ms
			estimate parameters from ms
			ms2SLiM
			SLiM forward simulator with estimated ms-parameters or take the output of ms as input
			SLiM2PolymorphismTableFile
			
			AddPopGenSimulation2DB.py
			
		"""
		if workflow is None:
			workflow = self
		
		sys.stderr.write("Adding jobs for pop-gen simulation #jobs=%s... \n"%\
							(self.no_of_jobs))
		
		returnData = PassingData()
		returnData.jobDataLs = []
		
		passingData = PassingData(fileBasenamePrefix=None, \
					outputDirPrefix=outputDirPrefix, \
					jobData=None,\
					preReduceReturnData=None,\
					association_group_key2orderIndex = {},\
					association_group_key2resultList = {},\
					association_group_key2reduceAssociationPeakJobMatrix = {},\
					association_group_key2countAssociationLocusJobList = {},\
					resultID2defineLandscapeJobData = {},
					)
		
		preReduceReturnData = self.preReduce(workflow=workflow, outputDirPrefix=outputDirPrefix, \
									passingData=passingData, transferOutput=False,\
									**keywords)
		
		mapDirJob = preReduceReturnData.mapDirJob
		plotOutputDirJob = preReduceReturnData.plotOutputDirJob
		countAssociationLocusOutputDirJob = preReduceReturnData.countAssociationLocusOutputDirJob
		reduceOutputDirJob = preReduceReturnData.reduceOutputDirJob
		
		passingData.preReduceReturnData = preReduceReturnData
		
		#add output pedigree job
		
		for i in xrange(self.noOfReplicates):
			popGenSimulationFolderJob = self.addMkDirJob(outputDir=os.path.join(mapDirJob.output, 'popGenSim%s'%(i)), \
														parentJobLs=[mapDirJob])
			#pending user choice, use ms/sfs-code/slim/ms & slim combination 
			msOutputFile = File(os.path.join(popGenSimulationFolderJob.output, \
									'sim%s_msOutput.txt.gz'%(i)))
			popSimulationJob = self.addMSSimulationJob(outputFile=msOutputFile, \
								recombinationRate=self.recombinationRate, mutationRate=self.mutationRate, \
								initialEffectivePopulationSize=self.initialEffectivePopulationSize, \
								otherParametersPassedToPopGenSimulator=self.otherParametersPassedToPopGenSimulator, \
								sampleSize=self.sampleSize, noOfLociToSimulate=self.noOfLociToSimulate, \
								simulateLocusLengthList=self.simulateLocusLengthList, \
								parentJobLs=[popGenSimulationFolderJob], \
								extraDependentInputLs=None, extraOutputLs=None, \
								transferOutput=False, extraArguments=None, extraArgumentList=None, \
								job_max_memory=2000, walltime=180)
			
			#. convert ms pop-gen output 2 polymorphism-table file
			msOutputHDF5File = File(os.path.join(popGenSimulationFolderJob.output, \
									'sim%s_msOutput.h5'%(i)))
			msOutput2PolymorphismTableFileJob = self.addGenericJob(executable=self.msOutput2PolymorphismTableFile, \
					inputFile=popSimulationJob.output, \
					outputFile=msOutputHDF5File,\
					parentJob=None, parentJobLs=[popGenSimulationFolderJob, popSimulationJob], \
					extraDependentInputLs=None, extraOutputLs=None, transferOutput=False, \
					frontArgumentList=None, \
					extraArguments=None, \
					extraArgumentList=None, job_max_memory=2000,  \
					no_of_cpus=None, walltime=None)
			
			#. add polymorphism-table file to db
			logFile = File(os.path.join(popGenSimulationFolderJob.output, "sim%s_2DB.log"%(i)))
			extraArgumentList = ["--r %s"%self.recombinationRate, "--rho %s"%popSimulationJob.rho, "--mu %s"%self.mutationRate,\
								"--theta %s"%popSimulationJob.theta, "--n0 %s"%self.initialEffectivePopulationSize,\
								"--no_of_populations 1", "--no_of_chromosomes %s"%self.sampleSize,\
								"--chromosome_length %s"%popSimulationJob.locusLength,\
								"--replicate_index %s"%(i)]
			"""
			extraArgumentList.append("--parent_pop_gen_simulation_type_id %s"%self.parent_pop_gen_simulation_type_id)
			"""
			simulation2DBJob = self.addPutStuffIntoDBJob(executable=self.AddPopGenSimulation2DB, \
					inputFileList=[msOutput2PolymorphismTableFileJob.output], \
					logFile=logFile, commit=True, \
					parentJobLs=[popGenSimulationFolderJob, msOutput2PolymorphismTableFileJob], \
					extraDependentInputLs=None, transferOutput=True, extraArguments=None, \
					extraArgumentList=extraArgumentList,\
					job_max_memory=10, sshDBTunnel=self.needSSHDBTunnel)
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2013.2.26
		"""
		parentClass.registerCustomExecutables(self, workflow=workflow)
		
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		clusters_size = self.clusters_size
		site_handler = self.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		#2013.3.8
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.expanduser(self.sfs_code_path), \
												name="sfs_code", clusterSizeMultipler=0.05)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.expanduser(self.slim_path), \
												name="slim", clusterSizeMultipler=0.3)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.expanduser(self.msHOT_lite_path), \
												name="msHOT_lite", clusterSizeMultipler=0.5)
		self.registerOneExecutableAsFile(pythonVariableName="msHOT_liteExecutableFile", path=self.msHOT_lite_path)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, 'shell/pipeCommandOutput2File.sh'), \
										name='msShellPipe', clusterSizeMultipler=1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(\
										path=os.path.join(self.vervetSrcPath, 'db/input/AddPopGenSimulation2DB.py'), \
										name="AddPopGenSimulation2DB", clusterSizeMultipler=0.2)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(\
										path=os.path.join(self.pymodulePath, 'popgen/converter/SFS_CODE_Output2PolymorphismTableFile.py'), \
										name="SFS_CODE_Output2PolymorphismTableFile", clusterSizeMultipler=0.2)
		self.addOneExecutableFromPathAndAssignProperClusterSize(\
										path=os.path.join(self.pymodulePath, 'popgen/converter/msOutput2PolymorphismTableFile.py'), \
										name="msOutput2PolymorphismTableFile", clusterSizeMultipler=0.2)
		
	def run(self):
		"""
		2011-7-11
		"""
		
		pdata = self.setup_run()
		
		self.addAllJobs(data_dir=self.data_dir, outputDirPrefix=self.pegasusFolderName, transferOutput=True)
		
		self.end_run()


if __name__ == '__main__':
	main_class = PopGenSimulationWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()