#!/usr/bin/env python
"""
Examples:
	# 2012.2.28
	%s -o workflow/DrawResultPeakWorkflowPublicPhenotype.xml -p DB_PASSWORD -u yh -w GENOME_DB_PASSWORD
		 -i 32:1:2,32:7:3,80:1:3,80:32:3
		 -N /Network/Data/250k/tmp-yh/phenotype/phenotype.tsv -s1 -B condorpool -D condorpool
	
	%s 
	
	%s
	
Description:
	2012.3.26
		1. genome db and stock_250k db have to be sync'ed. The candidate genes in stock_250k db refers to genome.gene.id.
		2. one biology category has one preferred gene list or a list of preferred gene list (if table BiologyCategory doesn't specify).
		
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], )

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, figureOutDelimiter
from Pegasus.DAX3 import *
import Stock_250kDB, csv
from AbstractVariationWorkflow import AbstractVariationWorkflow

class DrawResultPeakWorkflow(AbstractVariationWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVariationWorkflow.option_default_dict.copy()
	common_option_dict = {
						('genome_drivername', 1,):['mysql', '', 1, 'which type of database is the genome database? mysql or postgresql', ],\
						('genome_hostname', 1, ): ['banyan', '', 1, 'hostname of the genome db server', ],\
						('genome_dbname', 1, ): ['genome', '', 1, 'genome database name', ],\
						('genome_schema', 0, ): ['', '', 1, 'genome database schema name', ],\
						('genome_db_user', 1, ): ['yh', 'g', 1, 'genome database username', ],\
						('genome_db_passwd', 1, ): [None, 'w', 1, 'genome database password', ],\
						
						('call_analysis_peak_type_id_list', 1, ):[None, 'i', 1, 'a comma separated list of call_method_id:analysis_method_id:peak_type_id. \
			i.e.  32:1:2,32:7:3,80:1:3,80:32:3'],\
						("phenotype_method_id_ls", 0, ): [None, 'y', 1, 'comma/dash-separated phenotype_method id list, like 1,3-7. Default is all.'],\
						
						('peakPadding', 1, int): [10000, '', 1, 'the extension for each peak on both sides. Rationale is if two peaks are ...'],\
						('genePadding', 0, int): [20000, 'x', 1, "the extension around a gene on both sides to allow association between a locus and a gene. Proper distance is LD-decay."],\
						('tax_id', 0, int): [3702, '', 1, 'Taxonomy ID to get gene position and coordinates.'],\
						
						('biology_category_id', 0, int): [None, '', 1, 'filter phenotype based on biology category id. Default is no filter'],\
						('access', 0, int): [None, 's', 1, 'Restrict phenotype via access field in db. 1: public phenotypes, 2: restricted ones. Default is no filter'],\
						('phenotype_fname', 1, ): [None, 'N', 1, 'phenotype file, if snp_matrix_fname is given, this is needed as well.', ],\
						
						}
	
	option_default_dict.update(common_option_dict)
	
	def __init__(self, inputLs, **keywords):
		"""
		2012.2.28
		"""
		AbstractVariationWorkflow.__init__(self, **keywords)
		
		if getattr(self, 'phenotype_method_id_ls', None):
			self.phenotype_method_id_ls = getListOutOfStr(self.phenotype_method_id_ls, data_type=int)
			self.phenotype_method_id_ls.sort()
		else:
			self.phenotype_method_id_ls = []
		
		self.call_analysis_peak_type_id_list = getListOutOfStr(self.call_analysis_peak_type_id_list, data_type=str)
		call_analysis_peak_type_id_list = []
		for call_analysis_peak_type_id in self.call_analysis_peak_type_id_list:
			call_analysis_peak_type_id = call_analysis_peak_type_id.split(':')
			call_analysis_peak_type_id = map(int, call_analysis_peak_type_id)
			call_analysis_peak_type_id_list.append(call_analysis_peak_type_id)
		self.call_analysis_peak_type_id_list = call_analysis_peak_type_id_list
	
	def addPickleGenomeRBDictJob(self, workflow, executable=None, \
							outputF=None, genePadding=None, tax_id=3702,\
							parentJobLs=[], job_max_memory=100, extraDependentInputLs=[], \
							transferOutput=False, **keywords):
		"""
		2012.3.22
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		
		job.addArguments("-v", self.genome_drivername, "-z", self.genome_hostname, "-d", self.genome_dbname, \
						"-u", self.genome_db_user, "-p", self.genome_db_passwd,\
						"--genePadding=%s"%(genePadding), "--tax_id=%s"%(tax_id), "-o", outputF)
		if self.genome_schema:
			job.addArguments("--schema=%s"%self.genome_schema)
		job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = outputF
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		workflow.addJob(job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		return job
	
	def addPickleGeneAnnotationJob(self, workflow, executable=None, \
							outputF=None, tax_id=3702,\
							parentJobLs=[], job_max_memory=100, extraDependentInputLs=[], \
							transferOutput=False, **keywords):
		"""
		2012.3.26
			
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		
		job.addArguments("-v", self.genome_drivername, "-z", self.genome_hostname, "-d", self.genome_dbname, \
						"-u", self.genome_db_user, "-p", self.genome_db_passwd,\
						"--tax_id=%s"%(tax_id), "-o", outputF)
		if self.genome_schema:
			job.addArguments("--schema=%s"%self.genome_schema)
		job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = outputF
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		workflow.addJob(job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		return job
	
	
	def addPickleSNPInfoJob(self, workflow, executable=None, \
						outputF=None, call_method_id=None, \
						parentJobLs=[], job_max_memory=100, extraDependentInputLs=[], \
						transferOutput=False, **keywords):
		"""
		2012.3.22
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments("-v", self.drivername, "-z", self.hostname, "-d", self.dbname, \
						"-u", self.db_user, "-p", self.db_passwd, \
						"--call_method_id=%s"%(call_method_id), "-F", outputF)
		job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = outputF
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		workflow.addJob(job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		return job
	
	def addOutputMultiGWASOverlapPeakSpanJob(self, workflow, executable=None, \
						outputF=None, peakPadding=None,list_type_id_list=None, result_id_peak_type_id_ls=None,\
						genePadding=None, tax_id=3702, genomeRBDictPickleFile=None, \
						parentJobLs=[], job_max_memory=100, extraDependentInputLs=[], \
						transferOutput=False, **keywords):
		"""
		2012.3.22
			argument list_type_id_list is in string format. "129,137"
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments("-v", self.drivername, "-z", self.hostname, "-d", self.dbname, \
						"-u", self.db_user, "-p", self.db_passwd, \
						"--genome_drivername=%s"%self.genome_drivername, "--genome_hostname=%s"%self.genome_hostname, \
						"--genome_dbname=%s"%self.genome_dbname, \
						"--genome_db_user=%s"%self.genome_db_user, "--genome_schema=%s"%self.genome_schema, \
						"--genome_db_passwd=%s"%self.genome_db_passwd,\
						"--peakPadding=%s"%(peakPadding), \
						"--result_id_peak_type_id_ls=%s"%(result_id_peak_type_id_ls), \
						"--genePadding=%s"%(genePadding), \
						"--tax_id=%s"%(tax_id), "-o", outputF)
		if list_type_id_list:
			job.addArguments("--list_type_id_list=%s"%(list_type_id_list))
		if genomeRBDictPickleFile:
			job.addArguments("-m", genomeRBDictPickleFile)
			job.uses(genomeRBDictPickleFile, transfer=True, register=True, link=Link.INPUT)
		job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = outputF
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		workflow.addJob(job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		return job
	
	def addDrawSNPRegionJob(self, workflow, executable=None, \
				inputF=None, call_method_id=None, snpMatrixFile=None, phenotypeFile=None, output_dir=None,\
				results_directory=None,analysis_method_id_ls=None, geneAnnotationPickleFile=None,\
				list_type_id_list=None, snp_matrix_data_type=1, exclude_accessions_with_NA_phenotype=0,\
				snpInfoPickleFile=None, label_gene=1, min_MAF=0.1, min_distance=20000,\
				logFile=None,\
				parentJobLs=[], job_max_memory=2000, extraDependentInputLs=[], \
				transferOutput=False, **keywords):
		"""
		2012.3.22
			argument analysis_method_id_ls, list_type_id_list are in string format. "129,137"
		
			DrawSNPRegion.py -I /Network/Data/250k/db/dataset/call_method_80.tsv -N /Network/Data/250k/tmp-yh/phenotype/phenotype.tsv 
				-l 129 -o /Network/Data/250k/tmp-yh/snp_region -j /Network/Data/250k/tmp-yh/at_gene_model_pickelf 
				-e 80 -u yh -s -a 1,32 -z banyan -u yh
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments("-v", self.drivername, "-z", self.hostname, "-d", self.dbname, \
						"-u", self.db_user, "-p", self.db_passwd, \
						"-i", inputF,\
						"--call_method_id=%s"%(call_method_id), \
						"-I", snpMatrixFile, \
						"-N", phenotypeFile, "--output_dir=%s"%(output_dir), \
						"--analysis_method_id_ls=%s"%(analysis_method_id_ls),\
						"-j", geneAnnotationPickleFile, "--list_type_id_list=%s"%(list_type_id_list), \
						"--snp_matrix_data_type=%s"%(snp_matrix_data_type), \
						"--min_MAF=%s"%(min_MAF), "--min_distance=%s"%(min_distance),
						)
		job.uses(inputF, transfer=True, register=True, link=Link.INPUT)
		job.uses(snpMatrixFile, transfer=True, register=True, link=Link.INPUT)
		job.uses(phenotypeFile, transfer=True, register=True, link=Link.INPUT)
		job.uses(geneAnnotationPickleFile, transfer=True, register=True, link=Link.INPUT)
		if exclude_accessions_with_NA_phenotype:
			job.addArguments("--exclude_accessions_with_NA_phenotype")
		if label_gene:
			job.addArguments("--label_gene")
		if results_directory:
			job.addArguments("--results_directory=%s"%(results_directory))
		if snpInfoPickleFile:
			job.addArguments("-F", snpInfoPickleFile)
			job.uses(snpInfoPickleFile, transfer=True, register=True, link=Link.INPUT)
		if logFile:
			job.addArguments("-A", logFile)
			job.uses(logFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
			job.output = logFile
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		workflow.addJob(job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		return job
	
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
	
		PickleGenomeRBDict = Executable(namespace=namespace, name="PickleGenomeRBDict", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		PickleGenomeRBDict.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/PickleGenomeRBDict.py"), site_handler))
		#PickleGenomeRBDict.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(PickleGenomeRBDict)
		workflow.PickleGenomeRBDict = PickleGenomeRBDict
		
		#used to generate the gene_annotation pickle for DrawSNPRegion.py
		GenomeDB = Executable(namespace=namespace, name="GenomeDB", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		GenomeDB.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "GenomeDB.py"), site_handler))
		#GenomeDB.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(GenomeDB)
		workflow.GenomeDB = GenomeDB
		
		
		OutputMultiGWASOverlapPeakSpan = Executable(namespace=namespace, name="OutputMultiGWASOverlapPeakSpan", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		OutputMultiGWASOverlapPeakSpan.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "mapper/OutputMultiGWASOverlapPeakSpan.py"), site_handler))
		OutputMultiGWASOverlapPeakSpan.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(OutputMultiGWASOverlapPeakSpan)
		workflow.OutputMultiGWASOverlapPeakSpan = OutputMultiGWASOverlapPeakSpan
		
		PickleSNPInfo = Executable(namespace=namespace, name="PickleSNPInfo", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		PickleSNPInfo.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "mapper/PickleSNPInfo.py"), site_handler))
		#PickleSNPInfo.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(PickleSNPInfo)
		workflow.PickleSNPInfo = PickleSNPInfo
		
		
		DrawSNPRegion = Executable(namespace=namespace, name="DrawSNPRegion", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		DrawSNPRegion.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "DrawSNPRegion.py"), site_handler))
		#DrawSNPRegion.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(DrawSNPRegion)
		workflow.DrawSNPRegion = DrawSNPRegion
		
		
	def addJobs(self, workflow, db_250k=None, inputData=None, \
			biologyCategoryID2PhenotypeID2Data=None, pegasusFolderName="", \
			genePadding=20000, tax_id=3702, peakPadding=10000, phenotypeFile=None, call_method_id_set=None,\
			results_directory=None):
		"""
		2012.3.21
		"""
		sys.stderr.write("Adding SNPRegion drawing jobs on %s biology categories ..."%(len(biologyCategoryID2PhenotypeID2Data)))
		returnJobData = PassingData(jobDataLs = [])
		
		topOutputDir = pegasusFolderName
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		
		#add the PickleGenomeRBDict job
		genomeRBDictPickleFile = os.path.join(topOutputDir, 'genomeRBDict_tax%s_padding%s.pickle'%(tax_id, genePadding))
		pickleGenomeRBDictJob = self.addPickleGenomeRBDictJob(workflow, workflow.PickleGenomeRBDict, outputF=genomeRBDictPickleFile, \
									genePadding=genePadding, tax_id=tax_id, \
									parentJobLs=[topOutputDirJob], job_max_memory=200, \
									extraDependentInputLs=[], transferOutput=True)
		
		#add the PickleGeneAnnotation job
		geneAnnotationPickleFile = os.path.join(topOutputDir, 'geneAnnotation_tax%s.pickle'%(tax_id))
		geneAnnotationPickleJob = self.addPickleGeneAnnotationJob(workflow, workflow.GenomeDB, outputF=geneAnnotationPickleFile, \
									tax_id=tax_id, \
									parentJobLs=[topOutputDirJob], job_max_memory=200, \
									extraDependentInputLs=[], transferOutput=True)
		
		no_of_jobs = 3
		
		#add PickleSNPInfo job for each call method
		call_method_id2JobData = {}
		for call_method_id in call_method_id_set:
			call_method = Stock_250kDB.CallMethod.get(call_method_id)
			snpMatrixFile = self.registerOneInputFile(workflow, \
											db_250k.reScalePathByNewDataDir(filePath=call_method.filename, newDataDir=results_directory),\
											folderName=pegasusFolderName)
			outputF = File(os.path.join(topOutputDir, 'SNPInfo_LocusType%s.pickle'%(call_method.locus_type_id)))
			pickleSNPInfoJob = self.addPickleSNPInfoJob(workflow, workflow.PickleSNPInfo, \
									outputF=outputF, call_method_id=call_method_id, \
									parentJobLs=[topOutputDirJob], job_max_memory=100, \
									extraDependentInputLs=[], transferOutput=True)
			call_method_id2JobData[call_method_id] = PassingData(job=pickleSNPInfoJob, snpMatrixFile=snpMatrixFile)
			no_of_jobs += 1
		
		#one folder for each biology category
		for biology_category_id, phenotype_id2data in biologyCategoryID2PhenotypeID2Data.iteritems():
			#add a mkdirJob
			folderName = os.path.join(topOutputDir, 'biology_category_%s'%(biology_category_id))
			folderJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=folderName,\
											parentJobLs=[topOutputDirJob])
			no_of_jobs += 1
			biologyCategory = Stock_250kDB.BiologyCategory.get(biology_category_id)
			
			list_type_id_list = biologyCategory.returnGeneListIDList()
			list_type_id_in_str_list = map(str, list_type_id_list)
			list_type_id_list_str = ','.join(list_type_id_in_str_list)
			
			for phenotype_id, result_data in phenotype_id2data.iteritems():
				result_peak_type_id_ls = result_data.result_peak_type_id_ls
				call_method_id_set = result_data.call_method_id_set
				result_id_peak_type_id_ls = []
				analysis_method_id_in_str_ls = []
				for result, peak_type_id in result_peak_type_id_ls:
					result_id_peak_type_id_ls.append('%s:%s'%(result.id, peak_type_id))
					analysis_method_id_in_str_ls.append(str(result.analysis_method_id))
				result_id_peak_type_id_ls_str = ','.join(result_id_peak_type_id_ls)
				peakSpanOutputFile = File(os.path.join(folderName, 'phenotype_%s_result_%s.tsv'%
													(phenotype_id, result_id_peak_type_id_ls_str)))
				multiPeakSpanJob = self.addOutputMultiGWASOverlapPeakSpanJob(workflow, workflow.OutputMultiGWASOverlapPeakSpan, \
								outputF=peakSpanOutputFile, peakPadding=peakPadding, \
								list_type_id_list=list_type_id_list_str, result_id_peak_type_id_ls=result_id_peak_type_id_ls_str, \
								genePadding=genePadding, tax_id=tax_id, genomeRBDictPickleFile=genomeRBDictPickleFile, \
								parentJobLs=[folderJob, pickleGenomeRBDictJob], job_max_memory=500, \
								extraDependentInputLs=[], transferOutput=True)
				no_of_jobs += 1
				for call_method_id in call_method_id_set:
					call_method = Stock_250kDB.CallMethod.get(call_method_id)
					if call_method.locus_type_id==2:	#2012.3.27 CNV locus type. no need to convert alleles into binary form.
						snp_matrix_data_type=4
						#2012.3.26 these CNV-derived SNP dataset doesn't need its alleles to be converted to binary form as it's already binary.
						#need_convert_alleles2binary = False
						#useAlleleToDetermineAlpha = False
					else:
						snp_matrix_data_type=1
					callMethodJobData = call_method_id2JobData[call_method_id]
					pickleSNPInfoJob = callMethodJobData.job
					snpMatrixFile = callMethodJobData.snpMatrixFile
					output_dir = folderName 	#go to the biology category
					logFile = File(os.path.join(folderName, 'call_%s_phenotype_%s_result_%s_drawSNPRegion.log')%\
								(call_method_id, phenotype_id, result_id_peak_type_id_ls_str))
					analysis_method_id_ls_str = ','.join(analysis_method_id_in_str_ls)
					drawSNPRegionJob = self.addDrawSNPRegionJob(workflow, executable=workflow.DrawSNPRegion, \
								inputF=peakSpanOutputFile, call_method_id=call_method_id, snpMatrixFile=snpMatrixFile, \
								phenotypeFile=phenotypeFile, output_dir=output_dir, results_directory=results_directory,\
								analysis_method_id_ls=analysis_method_id_ls_str, \
								geneAnnotationPickleFile=geneAnnotationPickleJob.output, \
								list_type_id_list=list_type_id_list_str, \
								snp_matrix_data_type=snp_matrix_data_type, exclude_accessions_with_NA_phenotype=0,\
								snpInfoPickleFile=pickleSNPInfoJob.output, label_gene=1, min_MAF=0.1, min_distance=20000,\
								logFile=logFile,\
								parentJobLs=[geneAnnotationPickleJob, pickleSNPInfoJob, multiPeakSpanJob], \
								job_max_memory=3500, extraDependentInputLs=[], \
								transferOutput=True)
					no_of_jobs += 1
		sys.stderr.write("%s jobs.\n"%(no_of_jobs))
		return returnJobData
	
	def getBiologyCategoryID2PhenotypeMethodID2Data(self, result_id_ls_peak_type_list=None):
		"""
		2012.3.22
		
		"""
		sys.stderr.write("Getting biologyCategoryID2PhenotypeID2Data ...")
		call_method_id_set = set()
		phenotype_method_id_set = set()
		
		biologyCategoryID2PhenotypeID2Data = {}
		no_of_results = 0
		for result_id_ls, peak_type_id in result_id_ls_peak_type_list:
			for result_id in result_id_ls:
				no_of_results += 1
				result = Stock_250kDB.ResultsMethod.get(result_id)
				biologyCategoryID = result.phenotype_method.biology_category_id
				call_method_id_set.add(result.call_method_id)
				if biologyCategoryID not in biologyCategoryID2PhenotypeID2Data:
					biologyCategoryID2PhenotypeID2Data[biologyCategoryID] = {}
				PhenotypeID2Data = biologyCategoryID2PhenotypeID2Data[biologyCategoryID]
				phenotype_method_id = result.phenotype_method_id
				if phenotype_method_id not in PhenotypeID2Data:
					PhenotypeID2Data[phenotype_method_id] = PassingData(result_peak_type_id_ls=[], call_method_id_set=set())
				PhenotypeID2Data[phenotype_method_id].result_peak_type_id_ls.append((result, peak_type_id))
				PhenotypeID2Data[phenotype_method_id].call_method_id_set.add(result.call_method_id)
				phenotype_method_id_set.add(phenotype_method_id)
		
		sys.stderr.write("%s categories, %s results, %s phenotypes.\n"%(len(biologyCategoryID2PhenotypeID2Data),\
								no_of_results, len(phenotype_method_id_set)))
		return biologyCategoryID2PhenotypeID2Data, call_method_id_set
	
	def run(self):
		"""
		2012.3.22
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
		phenotype_method_id_ls.sort()
		
		result_id_ls_peak_type_list = []
		for call_analysis_peak_type_id in self.call_analysis_peak_type_id_list:
			call_method_id, analysis_method_id, peak_type_id = call_analysis_peak_type_id
			
			result_query = db_250k.getResultLs(call_method_id=call_method_id, analysis_method_id_ls=[analysis_method_id], \
						phenotype_method_id_ls=phenotype_method_id_ls, cnv_method_id=None)
			result_id_ls = []
			for result in result_query:
				result_id_ls.append(result.id)
			
			#make sure the entries with (result_id, self.result_peak_type_id) exists in ResultPeak
			result_id_ls = db_250k.filterResultIDLsBasedOnResultPeak(result_id_ls, peak_type_id)
			result_id_ls_peak_type_list.append((result_id_ls, peak_type_id))
		biologyCategoryID2PhenotypeID2Data, call_method_id_set = self.getBiologyCategoryID2PhenotypeMethodID2Data(result_id_ls_peak_type_list)
		
		
		# Create a abstract dag
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = self.initiateWorkflow(workflowName)
		
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		#get list of result_id,peak_type_id 
		phenotypeFile = self.registerOneInputFile(workflow, self.phenotype_fname, folderName=self.pegasusFolderName)
		#2012.3.26 generated by GenomeDB.py on the fly
		#geneAnnotationPickleFile =  self.registerOneInputFile(workflow, self.gene_annotation_pickleFname, folderName=self.pegasusFolderName)
		
		self.addJobs(workflow, db_250k=db_250k, inputData=None, \
			biologyCategoryID2PhenotypeID2Data=biologyCategoryID2PhenotypeID2Data, pegasusFolderName=self.pegasusFolderName, \
			genePadding=self.genePadding, tax_id=self.tax_id, peakPadding=self.peakPadding, \
			phenotypeFile=phenotypeFile, call_method_id_set=call_method_id_set,\
			results_directory=self.results_directory)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		del outf


if __name__ == '__main__':
	main_class = DrawResultPeakWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()