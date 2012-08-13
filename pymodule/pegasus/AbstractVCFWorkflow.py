#!/usr/bin/env python
"""
2012.1.17
	a common class for pegasus workflows that work on VCF variant files
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq
from Pegasus.DAX3 import *
from AbstractNGSWorkflow import AbstractNGSWorkflow

class AbstractVCFWorkflow(AbstractNGSWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractNGSWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('inputDir', 0, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('minDepth', 0, float): [0, 'm', 1, 'minimum depth for a call to regarded as non-missing', ],\
						})
	def __init__(self,  **keywords):
		"""
		2012.1.17
		"""
		AbstractNGSWorkflow.__init__(self, **keywords)
		if getattr(self, "inputDir", None):
			self.inputDir = os.path.abspath(self.inputDir)
		
	def registerExecutables(self, workflow=None):
		"""
		"""
		AbstractNGSWorkflow.registerExecutables(self, workflow)
		
		if not workflow:
			workflow = self
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		clusters_size = self.clusters_size
		site_handler = self.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		executableList = []
		noClusteringExecutableSet = set()	#2012.8.2 you don't want to cluster for some jobs.
		
		PlotVCFtoolsStat= Executable(namespace=namespace, name="PlotVCFtoolsStat", version=version, os=operatingSystem, arch=architecture, installed=True)
		PlotVCFtoolsStat.addPFN(PFN("file://" +  os.path.join(self.vervetSrcPath, "plot/PlotVCFtoolsStat.py"), site_handler))
		executableList.append(PlotVCFtoolsStat)
		noClusteringExecutableSet.add(PlotVCFtoolsStat)
		
		for executable in executableList:
			if executable not in noClusteringExecutableSet:
				executable.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.clusters_size))
			workflow.addExecutable(executable)
			setattr(workflow, executable.name, executable)
		
	def registerCommonExecutables(self, workflow=None):
		"""
		"""
		AbstractNGSWorkflow.registerCommonExecutables(self, workflow)
	
	def registerAllInputFiles(self, workflow=None, inputDir=None, input_site_handler=None, checkEmptyVCFByReading=False, pegasusFolderName='',\
							maxContigID=None, minContigID=None):
		"""
		2012.8.23 add maxContigID and minContigID to restrict input
		2012.7.27 add attribute file to each object in returnData.jobDataLs
		2012.5.9
			register the tbi file if it exists
		2012.3.1
			moved from CalculateTrioInconsistencyPipeline.py
		2012.1.9
			the returning data structure is changed to conform to some standard used across several functions
		2011-9-29
			vcf files only
		"""
		sys.stderr.write("Registering input files from %s ..."%(inputDir))
		returnData = PassingData(jobDataLs = [])
		fnameLs = os.listdir(inputDir)
		counter = 0
		real_counter = 0
		previous_reported_real_counter = ''
		for fname in fnameLs:
			counter += 1
			inputFname = os.path.join(inputDir, fname)
			if (maxContigID is not None and maxContigID!=0) or (minContigID is not None and minContigID!=0):
				try:
					contigID = int(self.getContigIDFromFname(os.path.basename(fname)))
					if (maxContigID is not None and maxContigID!=0) and contigID>maxContigID:
						continue
					if (minContigID is not None and minContigID!=0) and contigID<minContigID:
						continue
				except:
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
			if NextGenSeq.isFileNameVCF(fname, includeIndelVCF=False) and \
					not NextGenSeq.isVCFFileEmpty(inputFname, checkContent=checkEmptyVCFByReading):
				real_counter += 1
				inputBaseFname = os.path.basename(inputFname)
				inputF = File(os.path.join(pegasusFolderName, inputBaseFname))
				inputF.addPFN(PFN("file://" + inputFname, input_site_handler))
				inputF.abspath = inputFname
				workflow.addFile(inputF)
				
				tbi_F_absPath = "%s.tbi"%inputFname
				if os.path.isfile(tbi_F_absPath):	#it exists
					tbi_F = File(os.path.join(pegasusFolderName, "%s.tbi"%inputBaseFname))
					tbi_F.addPFN(PFN("file://" + tbi_F_absPath, input_site_handler))
					tbi_F.abspath = tbi_F_absPath
					workflow.addFile(tbi_F)
				else:
					tbi_F = None
				
				returnData.jobDataLs.append(PassingData(vcfFile=inputF, jobLs=[], tbi_F=tbi_F, file=inputF))
				if real_counter%200==0:
					sys.stderr.write("%s%s"%('\x08'*len(previous_reported_real_counter), real_counter))
					previous_reported_real_counter = repr(real_counter)
		sys.stderr.write("  %s non-empty VCF out of %s files.\n"%(len(returnData.jobDataLs), counter))
		return returnData

	def addPlotVCFtoolsStatJob(self, workflow=None, executable=None, inputFileList=None, outputFnamePrefix=None, \
							whichColumn=None, whichColumnLabel=None, whichColumnPlotLabel=None, need_svg=False, \
							logWhichColumn=True,\
							posColumnPlotLabel=None, chrLengthColumnLabel=None, chrColumnLabel=None, \
							minChrLength=1000000, posColumnLabel=None, minNoOfTotal=100,\
							figureDPI=300, ylim_type=2, samplingRate=0.0001,\
							parentJobLs=None, \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=True, job_max_memory=2000, sshDBTunnel=False, **keywords):
		"""
		# whichColumnPlotLabel and posColumnPlotLabel should not contain spaces or ( or ). because they will disrupt shell commandline
		
		2012.8.2 moved from vervet/src/CalculateVCFStatPipeline.py
		2012.8.1
			
			('whichColumn', 0, int): [3, 'w', 1, 'data from this column (index starting from 0) is plotted as y-axis value'],\
			('whichColumnLabel', 0, ): ["", 'W', 1, 'column label (in the header) for the data to be plotted as y-axis value, substitute whichColumn'],\
			('logWhichColumn', 0, int): [0, 'g', 0, 'whether to take log of whichColumn'],\
			('need_svg', 0, ): [0, 'n', 0, 'whether need svg output', ],\
			('whichColumnPlotLabel', 1, ): ['#SNPs in 100kb window', 'D', 1, 'plot label for data of the whichColumn', ],\
			('posColumnPlotLabel', 1, ): ['position', 'x', 1, 'x-axis label (posColumn) in manhattan plot', ],\
			('chrLengthColumnLabel', 1, ): ['chrLength', 'c', 1, 'label of the chromosome length column', ],\
			('chrColumnLabel', 1, ): ['CHR', 'C', 1, 'label of the chromosome column', ],\
			('minChrLength', 1, int): [1000000, 'm', 1, 'minimum chromosome length for one chromosome to be included', ],\
			('posColumnLabel', 1, ): ['BIN_START', 'l', 1, 'label of the position column, BIN_START for binned vcftools output. POS for others.', ],\
			('outputFnamePrefix', 0, ): [None, 'O', 1, 'output filename prefix (optional).'],\
			
				('minNoOfTotal', 1, int): [100, 'i', 1, 'minimum no of total variants (denominator of inconsistent rate)'],\
				('title', 1, ): [None, 't', 1, 'title for the figure.'],\
				('figureDPI', 1, int): [200, 'f', 1, 'dpi for the output figures (png)'],\
				('formatString', 1, ): ['-', '', 1, 'formatString passed to matplotlib plot'],\
				('ylim_type', 1, int): [1, 'y', 1, 'y-axis limit type, 1: 0 to max. 2: min to max'],\
				('samplingRate', 1, float): [0.001, 's', 1, 'how often you include the data'],\
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		if inputFileList:
			extraDependentInputLs.extend(inputFileList)
		extraArgumentList = ["-O %s"%outputFnamePrefix, '-m %s'%(minChrLength), '-i %s'%(minNoOfTotal), \
							'-f %s'%(figureDPI), '-y %s'%(ylim_type), '-s %s'%(samplingRate), '-l %s'%(posColumnLabel)]
		extraOutputLs = [File('%s.png'%(outputFnamePrefix)), File('%s_hist.png'%(outputFnamePrefix))]
		if need_svg:
			extraOutputLs.append(File('%s.svg'%(outputFnamePrefix)))
		key2ObjectForJob = {}
		if whichColumnLabel:
			extraArgumentList.append("-W %s"%(whichColumnLabel))
		if whichColumn:
			extraArgumentList.append("-w %s"%(whichColumn))
		if logWhichColumn:
			extraArgumentList.append('-g')
		if whichColumnPlotLabel:
			extraArgumentList.append("-D %s"%(whichColumnPlotLabel))
		if posColumnPlotLabel:
			extraArgumentList.append("-x %s"%(posColumnPlotLabel))
		if chrLengthColumnLabel:
			extraArgumentList.append("-c %s"%(chrLengthColumnLabel))
		if chrColumnLabel:
			extraArgumentList.append("-C %s"%(chrColumnLabel))
		if extraArguments:
			extraArgumentList.append(extraArguments)
		job= self.addGenericJob(executable=executable, inputFile=None, outputFile=None, \
				parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, \
				sshDBTunnel=sshDBTunnel, **keywords)
		self.addDBArgumentsToOneJob(job=job, objectWithDBArguments=self)
		#add all input files to the last (after db arguments,) otherwise, it'll mask others (cuz these don't have options).
		if inputFileList:
			for inputFile in inputFileList:
				if inputFile:
					job.addArguments(inputFile)
		return job