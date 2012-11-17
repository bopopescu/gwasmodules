#!/usr/bin/env python
"""
Examples:
	# 2011-10-10 call method 80, analysis method 1, min_score 5
	%s  -E 80 -a 1 -f 5 -o call80analysis1MinScore5.xml -u yh -z banyan -c -D condorpool
	
	# 2011-10-16 call method 57 or cnv method 20, analysis method 1, min score 4
	%s -q 57 -A 20 -a 1 -f 4 -o call57cnv20analysis1MinScore4.xml -u yh -z banyan -c -D condorpool
	
	#ditto but analysis method 7 and min score 3 (make sure -c is added, otherwise nothing will be stored in db)
	%s -q 57 -A 20 -a 7 -f 3 -o call57cnv20analysis7MinScore3.xml -u yh -z banyan -c -B condorpool -D condorpool
	#ditto but analysis method 32
	%s -q 57 -A 20 -a 32 -f 3 -o call57cnv20analysis32MinScore3.xml -u yh  -z banyan -c -B condorpool -D condorpool

	mpiexec ~/script/variation/src/MpiAssociation.py -i ~/panfs/250k/dataset/call_method_17_test.tsv
		-p ~/panfs/250k/phenotype.tsv -o ~/panfs/250k/association_results/call_method_17_test_y3.tsv -y3 -w 187-190,210-221


	#2010-8-8 association on CNV. watch the "-n" argument
	mpiexec ~/script/variation/src/MpiAssociation.py -i ~/panfs/250k/CNV/NonOverlapCNVAsSNP_cnvMethod20.tsv -p ~/panfs/250k/phenotype/phenotype.tsv
		-o ~/panfs/250k/association_results/cnvMethod20/cnvMethod20_y3_pheno.tsv
		-y3 -w 39-59,80-82,210-212,32-38,65-74,183-186,191-194,260-263,362-380,161-176,264-267,849-947 -n
	
	#2012.9.28 add --getPublicPhenotype if u want to run association only on published phenotype and its values
	# -n: --noSNPAlleleOrdinalConversion 
	# -c: commit, (make sure -c is added, otherwise nothing will be stored in db)
	%s -A 57,75 -n -c -w 314-351 -a 1,32 -K /Network/Data/250k/db/dataset/call_method_75_K.tsv
		-o workflow/association/association/Call57_75Phenotype314_351Analysis1_32.xml
		-u yh -d stock_250k -z banyan  -j condorpool -l condorpool
		#--getPublicPhenotype

Description:
	2012.6.5
		output a workflow that runs Association.py on specified input and imports the result into db.
		associations already in db will be skipped.
		
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

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
from Association import Association
from AbstractVariationWorkflow import AbstractVariationWorkflow

class AssociationWorkflow(AbstractVariationWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVariationWorkflow.option_default_dict
	option_default_dict.update(Association.common_option_dict.copy())
	option_default_dict.pop(('input_fname', 1, ))
	option_default_dict.pop(('phenotype_fname', 1, ))
	
	option_default_dict.update({
			('getPublicPhenotype', 0, int):[0, '', 0, 'toggle to get public phenotype only, phenotype_method.access==public and phenotype_avg.ready_for_publication=1'],\
			('call_method_id_ls', 1, ):[None, 'A', 1, 'Restrict results based on list of call_method_id. Default is no such restriction.'],\
			('analysis_method_id_ls', 1, ):['1,7', 'a', 1, 'Restrict results based on these analysis_methods. coma or dash-separated list'],\
			('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
			})
	"""
						{
						('cnv_method_id', 0, int):[None, 'A', 1, 'Restrict results based on this cnv_method. Default is no such restriction.'],\
						("phenotype_method_id_ls", 0, ): [None, 'y', 1, 'comma/dash-separated phenotype_method id list, like 1,3-7. Default is all.'],\
						
						})
	"""
	
	def __init__(self,  **keywords):
		"""
		2012.6.5
			to replace MpiAssociation.py
		"""
		AbstractVariationWorkflow.__init__(self, **keywords)
		
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
		2012.6.5
		"""
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		clusters_size = self.clusters_size
		site_handler = self.site_handler
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		
		Association = Executable(namespace=namespace, name="Association", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		Association.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "Association.py"), site_handler))
		executableClusterSizeMultiplierList.append((Association, 0))
		
		Results2DB_250k = Executable(namespace=namespace, name="Results2DB_250k", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		Results2DB_250k.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "Results2DB_250k.py"), site_handler))
		executableClusterSizeMultiplierList.append((Results2DB_250k, 0))
		
		OutputPhenotype = Executable(namespace=namespace, name="OutputPhenotype", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		OutputPhenotype.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "OutputPhenotype.py"), site_handler))
		executableClusterSizeMultiplierList.append((OutputPhenotype, 0))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	
	def addAssociationJob(self, executable=None, datasetFile=None, phenotypeFile=None, phenotype_method_id=None, \
						outputFile=None, kinshipFile=None, eigenVectorFile=None, genotypeFileToGenerateKinship=None, test_type=1,\
						min_data_point=3, noSNPAlleleOrdinalConversion=False, which_PC_index_ls=None,\
						parentJobLs=None, job_max_memory=100, job_max_walltime = 600, \
						extraArguments=None, extraDependentInputLs=None, \
						transferOutput=False, **keywords):
		"""
		2012.6.5
			job_max_walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
		"""
		if not extraDependentInputLs:
			extraDependentInputLs = []
		extraDependentInputLs.append(phenotypeFile)
		
		extraArgumentList = ['-P', phenotypeFile, '-y %s'%(test_type), '-w %s'%(phenotype_method_id),\
							'-m %s'%(min_data_point), ]
		
		if kinshipFile:
			extraArgumentList.extend(['-K', kinshipFile])
			extraDependentInputLs.append(kinshipFile)
		if eigenVectorFile:
			extraDependentInputLs.append(eigenVectorFile)
			extraArgumentList.extend(["-f", eigenVectorFile])
		if genotypeFileToGenerateKinship:
			extraArgumentList.extend(["-G", genotypeFileToGenerateKinship])
			extraDependentInputLs.append(genotypeFileToGenerateKinship)
		if noSNPAlleleOrdinalConversion:
			extraArgumentList.append("-n")
		if which_PC_index_ls:
			extraArgumentList.append("-W %s"%(which_PC_index_ls))
		if extraArguments:
			extraArgumentList.append(extraArguments)
		
		return self.addGenericJob(executable=executable, inputFile=datasetFile, outputFile=outputFile, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, job_max_walltime=job_max_walltime,
						**keywords)
	
	def addOutputPhenotypeJob(self, executable=None, outputFile=None, getRawPhenotypeData=False,\
						ecotype_table='stock.ecotype', phenotype_method_table='stock_250k.phenotype_method', \
						phenotype_avg_table='stock_250k.phenotype_avg', getPublicPhenotype=False,\
						parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=2000, sshDBTunnel=None, **keywords):
		"""
		2012.9.28
			add argument getPublicPhenotype
		2012.6.5
		"""
		extraArgumentList = ['-e %s'%(ecotype_table), '-m %s'%(phenotype_method_table), '-q %s'%(phenotype_avg_table)]
		if getRawPhenotypeData:
			extraArgumentList.append('-g')
		if getPublicPhenotype:
			extraArgumentList.append("--getPublicPhenotype")
		if extraArguments:
			extraArgumentList.append(extraArguments)
		
		job = self.addGenericJob(executable=executable, outputFile=outputFile, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, sshDBTunnel=sshDBTunnel,\
						**keywords)
		self.addDBArgumentsToOneJob(job, objectWithDBArguments=self)
		return job
	
	def addResult2DBJob(self, executable=None, inputFile=None, call_method_id=None, phenotype_method_id=None, \
					analysis_method_id=None, data_dir=None, results_method_type_id=1,\
					logFile=None, commit=False,\
					parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					extraArguments=None, job_max_memory=2000, sshDBTunnel=None, **keywords):
		"""
		2012.11.13 expand all short-arguments to become long ones
		2012.6.5
		"""
		extraArgumentList = ['--phenotype_method_id %s'%(phenotype_method_id), \
							'--analysis_method_id %s'%(analysis_method_id),\
							'--results_method_type_id %s'%(results_method_type_id)]
		if call_method_id:
			extraArgumentList.append('--call_method_id %s'%(call_method_id))
		if data_dir:
			extraArgumentList.append('--data_dir %s'%(data_dir))
		if commit:
			extraArgumentList.append('--commit')
		if logFile:
			extraArgumentList.extend(["--logFilename", logFile])
		
		job = self.addGenericDBJob(executable=executable, inputFile=inputFile, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, extraArguments=extraArguments,\
						job_max_memory=job_max_memory, sshDBTunnel=sshDBTunnel,\
						**keywords)
		return job
	
	
	def addJobs(self, db_250k=None, callMethodID2FileObject=None, kinshipFile=None, eigenVectorFile=None, phenotype_method_id_ls=[],\
			analysis_method_id_ls=[], genotypeFileToGenerateKinship=None, data_dir=None, getPublicPhenotype=False, commit=False, \
			transferOutput=True, needSSHDBTunnel=False, outputDirPrefix=""):
		"""
		2012.9.28
			add argument getPublicPhenotype
		2012.6.5
		"""
		sys.stderr.write("Adding association jobs for %s polymorphism data ... "%(len(callMethodID2FileObject)))
		no_of_jobs= 0
		
		
		topOutputDir = "%sAssociation"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow=self, mkdir=self.mkdirWrap, outputDir=topOutputDir)
		no_of_jobs += 1
		
		phenotypeFile = File(os.path.join(topOutputDir, 'phenotype.tsv'))
		outputPhenotypeJob = self.addOutputPhenotypeJob(executable=self.OutputPhenotype, outputFile=phenotypeFile, \
									getRawPhenotypeData=False, getPublicPhenotype=getPublicPhenotype,\
									parentJobLs=[topOutputDirJob], transferOutput=True, job_max_memory=2000,\
									sshDBTunnel=needSSHDBTunnel)
		no_of_jobs += 1
		
		returnData = PassingData()
		returnData.jobDataLs = []
		for analysis_method_id in analysis_method_id_ls:
			analysisMethod = Stock_250kDB.AnalysisMethod.get(analysis_method_id)
			if not analysisMethod:
				sys.stderr.write("Warning: analysis_method_id %s not in db. Skip.\n"%(analysis_method_id))
				continue
			for phenotype_method_id in phenotype_method_id_ls:
				phenotypeMethod = Stock_250kDB.PhenotypeMethod.get(phenotype_method_id)
				if not phenotypeMethod:
					sys.stderr.write("Warning: phenotype_method_id %s not in db. Skip.\n"%(phenotype_method_id))
					continue
				for callMethodID, fileObject in callMethodID2FileObject.iteritems():
					test_type = analysisMethod.association_test_type
					if not test_type:
						sys.stderr.write("Warning: analysisMethod %s has non-None test_type %s. Skip.\n"%(test_type))
						continue
					#2012.9.28 not in db
					rm = db_250k.getResultsMethod(call_method_id=callMethodID, phenotype_method_id=phenotype_method_id, \
											analysis_method_id=analysis_method_id, \
											cnv_method_id=None)
					if rm:
						sys.stderr.write("Warning: skip association for c=%s, p=%s, a=%s.\n"%(callMethodID, phenotype_method_id, analysis_method_id))
						continue
					outputFile = File(os.path.join(topOutputDir, '%s_%s_%s.tsv'%(callMethodID, phenotype_method_id, \
																				analysis_method_id)))
					associationJob = self.addAssociationJob(executable=self.Association, datasetFile=fileObject, \
						phenotypeFile=outputPhenotypeJob.output, phenotype_method_id=phenotype_method_id, \
						outputFile=outputFile, kinshipFile=kinshipFile, eigenVectorFile=eigenVectorFile, \
						genotypeFileToGenerateKinship=genotypeFileToGenerateKinship, test_type=test_type,\
						min_data_point=self.min_data_point, noSNPAlleleOrdinalConversion=self.noSNPAlleleOrdinalConversion, \
						which_PC_index_ls=self.which_PC_index_ls,\
						parentJobLs=[outputPhenotypeJob], job_max_memory=3500, job_max_walltime =200, \
						extraDependentInputLs=None, transferOutput=False)
					logFile = File(os.path.join(topOutputDir, '%s_%s_%s_2DB.log'%(callMethodID, phenotype_method_id, \
																				analysis_method_id)))
					result2DBJob = self.addResult2DBJob(executable=self.Results2DB_250k, inputFile=associationJob.output, \
									call_method_id=callMethodID, phenotype_method_id=phenotype_method_id, \
									analysis_method_id=analysis_method_id, data_dir=data_dir, \
									results_method_type_id=1,\
									logFile=logFile, commit=commit,\
									parentJobLs=[associationJob], transferOutput=transferOutput, \
									job_max_memory=500, sshDBTunnel=needSSHDBTunnel)
					returnData.jobDataLs.append(PassingData(jobLs=[result2DBJob], file=logFile, \
											fileList=[logFile]))
					no_of_jobs += 2
		sys.stderr.write("%s jobs. Done.\n"%(no_of_jobs))

		return returnData
	
	def run(self):
		"""
		2011-10
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		workflow = self.initiateWorkflow()
		
		self.registerExecutables()
		self.registerCustomExecutables()
		
		callMethodID2FileObject = {}
		for call_method_id in self.call_method_id_ls:
			callMethod = Stock_250kDB.CallMethod.get(call_method_id)
			if callMethod and callMethod.filename:
				datasetFile = self.registerOneInputFile(inputFname=callMethod.filename, folderName=self.pegasusFolderName)
				callMethodID2FileObject[callMethod.id] = datasetFile
			else:
				sys.stderr.write("WARNING: call method %s has no db entry or its .filename is empty.\n"%(call_method_id))
		if self.kinship_fname:
			kinshipFile = self.registerOneInputFile(inputFname=self.kinship_fname, folderName=self.pegasusFolderName)
		else:
			kinshipFile = None
		if self.eigen_vector_fname:
			eigenVectorFile = self.registerOneInputFile(inputFname=self.eigen_vector_fname, folderName=self.pegasusFolderName)
		else:
			eigenVectorFile = None
		if self.genotype_fname_to_generate_kinship:
			genotypeFileToGenerateKinship = self.registerOneInputFile(inputFname=self.genotype_fname_to_generate_kinship, \
																folderName=self.pegasusFolderName)
		else:
			genotypeFileToGenerateKinship = None
		
		self.addJobs(db_250k=self.db, callMethodID2FileObject=callMethodID2FileObject, kinshipFile=kinshipFile, \
					eigenVectorFile=eigenVectorFile, phenotype_method_id_ls=self.phenotype_method_id_ls,\
					analysis_method_id_ls=self.analysis_method_id_ls, \
					genotypeFileToGenerateKinship=genotypeFileToGenerateKinship, \
					data_dir=self.data_dir, \
					getPublicPhenotype=self.getPublicPhenotype,\
					commit=self.commit, \
					transferOutput=True, needSSHDBTunnel=self.needSSHDBTunnel, outputDirPrefix="")
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		self.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = AssociationWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
