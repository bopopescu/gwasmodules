#!/usr/bin/env python
"""
2011-10-11
	module that deals with pegasus-related stuff
"""
import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from Pegasus.DAX3 import Executable, File, Job, Link, PFN, Profile, Namespace
from pymodule.utils import PassingData

def addMkDirJob(workflow=None, mkdir=None, outputDir=None, namespace=None, version=None,\
			parentJobLs=None, extraDependentInputLs=None):
	"""
	2012.10.2, increment workflow.no_of_jobs
	2012.9.11
		make sure that parentJobLs and extraDependentInputLs are not None.
	2012.3.10
		add argument parentJobLs, extraDependentInputLs
	2011-11-28
		get namespace and version from workflow first
	2011-9-14
	"""
	# Add a mkdir job for any directory.
	job = Job(namespace=getattr(workflow, 'namespace', namespace), name=mkdir.name, \
				version=getattr(workflow, 'version', version))
	job.addArguments(outputDir)
	job.folder = outputDir	#custom attribute
	job.output = outputDir	#custom attribute
	workflow.addJob(job)
	if parentJobLs:
		for parentJob in parentJobLs:
			if parentJob:
				workflow.depends(parent=parentJob, child=job)
	if extraDependentInputLs:
		for input in extraDependentInputLs:
			if input is not None:
				job.uses(input, transfer=True, register=True, link=Link.INPUT)
	if hasattr(workflow, 'no_of_jobs'):	#2012.10.2
		workflow.no_of_jobs += 1
	return job

def registerRefFastaFile(workflow=None, refFastaFname=None, registerAffiliateFiles=True, input_site_handler='local',\
						checkAffiliateFileExistence=True, addPicardDictFile=True,\
						affiliateFilenameSuffixLs=['fai', 'amb', 'ann', 'bwt', 'pac', 'sa', 'rbwt', 'rpac', 'rsa', \
						'stidx', 'sthash'], folderName="reference"):
	"""
	suffix here doesn't include ".".
	
	2013.08.23 bugfix, check if workflow has a file registered before adding it
	2013.3.26 added refSAMtoolsFastaIndexF, refPicardFastaDictF into returnData
	2013.3.20 deduce needBWARefIndexJob, needSAMtoolsFastaIndexJob, needPicardFastaDictJob, needStampyRefIndexJob from missing suffixes
	2010.10.10 added argument folderName
	2012.5.23
		add an argument "addPicardDictFile" to offer user option to exclude this file (i.e. in registerBlastNucleotideDatabaseFile)
	2012.2.24
		dict is via picard, also required for GATK
		fai is via "samtools faidx" (index reference). also required for GATK
		amb', 'ann', 'bwt', 'pac', 'sa', 'rbwt', 'rpac', 'rsa' are all bwa index.
		stidx is stampy index.
		sthash is stampy hash.
	2012.2.23
		add two suffixes, stidx (stampy index) and sthash (stampy hash)
	2011-11-11
		if needAffiliatedFiles,
			all other files, with suffix in affiliateFilenameSuffixLs, will be registered (symlinked or copied) as well.
	"""
	returnData = PassingData(refFastaFList = [], needBWARefIndexJob=False, needSAMtoolsFastaIndexJob=False, \
							needPicardFastaDictJob=False, needStampyRefIndexJob=False, needBlastMakeDBJob=False,\
							refPicardFastaDictF=None, refSAMtoolsFastaIndexF=None)
	missingSuffixSet = set()	#2013.3.20
	
	if registerAffiliateFiles:
		refFastaF = File(os.path.join(folderName, os.path.basename(refFastaFname)))	#use relative path, otherwise, it'll go to absolute path
		# Add it into replica only when needed.
		refFastaF.addPFN(PFN("file://" + refFastaFname, input_site_handler))
		if not workflow.hasFile(refFastaF):	#2013.08.12
			workflow.addFile(refFastaF)
		returnData.refFastaFList.append(refFastaF)
		# If it's not needed, assume the index is done and all relevant files are in absolute path.
		# and no replica transfer
		
		#add extra affiliated files
		suffix2PathToFileLs = {}
		if addPicardDictFile:	#2012.5.23
			picardDictSuffix = 'dict'
			pathToFile = '%s.%s'%(os.path.splitext(refFastaFname)[0], picardDictSuffix)	#remove ".fasta" from refFastaFname
			if checkAffiliateFileExistence and not os.path.isfile(pathToFile):
				sys.stderr.write("Warning: %s don't exist or not a file on file system. skip registration.\n"%(pathToFile))
				missingSuffixSet.add(picardDictSuffix)
				#suffix2PathToFileLs.append(pathToFile)
			else:
				suffix2PathToFileLs[picardDictSuffix] = pathToFile
		for suffix in affiliateFilenameSuffixLs:
			pathToFile = '%s.%s'%(refFastaFname, suffix)
			if checkAffiliateFileExistence and not os.path.isfile(pathToFile):
				sys.stderr.write("Warning: %s don't exist or not a file on file system. skip registration.\n"%(pathToFile))
				missingSuffixSet.add(suffix)
				continue
			suffix2PathToFileLs[suffix]= pathToFile
		for suffix, pathToFile in suffix2PathToFileLs.iteritems():
			if checkAffiliateFileExistence and not os.path.isfile(pathToFile):
				sys.stderr.write("Warning: %s don't exist or not a file on file system. skip registration.\n"%(pathToFile))
				continue
			affiliateF = File(os.path.join(folderName, os.path.basename(pathToFile)))
			#use relative path, otherwise, it'll go to absolute path
			affiliateF.addPFN(PFN("file://" + pathToFile, input_site_handler))
			if not workflow.hasFile(affiliateF):	#2013.08.12
				workflow.addFile(affiliateF)
			returnData.refFastaFList.append(affiliateF)
			
			if suffix=='dict':	#2013.3.26
				returnData.refPicardFastaDictF = affiliateF
			elif suffix=='fai':
				returnData.refSAMtoolsFastaIndexF = affiliateF
	else:
		refFastaF = File(os.path.join(folderName, os.path.basename(refFastaFname)))
		returnData.refFastaFList.append(refFastaF)
	if 'bwt' in missingSuffixSet or 'pac' in missingSuffixSet:
		returnData.needBWARefIndexJob = True
	if 'fai' in missingSuffixSet:
		returnData.needSAMtoolsFastaIndexJob = True
		returnData.needPicardFastaDictJob = True
	if 'stidx' in missingSuffixSet or 'sthash' in missingSuffixSet:
		returnData.needStampyRefIndexJob = True
	if 'dict' in missingSuffixSet:
		returnData.needPicardFastaDictJob = True
	if 'nin' in missingSuffixSet or 'nhr' in missingSuffixSet or 'nsq' in missingSuffixSet:
		returnData.needBlastMakeDBJob = True
	return returnData

def setJobToProperMemoryRequirement(job=None, job_max_memory=500, no_of_cpus=1, walltime=180, sshDBTunnel=0):
	"""
	2013.06.22 if job_max_memory is None, then skip setting memory requirement
		if job_max_memory is "" or 0 or "0", then assign 500 (mb) to it.
	2012.8.15
		increase default walltime to 180
	2012.4.16
		add argument sshDBTunnel.
			=1: this job needs a ssh tunnel to access psql db on dl324b-1.
			=anything else: no need for that.
	2011-11-23
		set walltime default to 120 minutes (2 hours)
	2011-11-16
		add more requirements
	2011-11-11
		job_max_memory is in MB.
		walltime is in minutes.
	"""
	condorJobRequirementLs = []
	if job_max_memory == "" or job_max_memory == 0 or job_max_memory =="0":
		job_max_memory = 500
	if job_max_memory is not None: 
		job.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%(job_max_memory)))
		job.addProfile(Profile(Namespace.CONDOR, key="request_memory", value="%s"%(job_max_memory)))	#for dynamic slots
		condorJobRequirementLs.append("(memory>=%s)"%(job_max_memory))
	#2012.4.16
	if sshDBTunnel==1:
		condorJobRequirementLs.append("(sshDBTunnel==%s)"%(sshDBTunnel))	#use ==, not =.
	
	if no_of_cpus is not None:
		job.addProfile(Profile(Namespace.CONDOR, key="request_cpus", value="%s"%(no_of_cpus)) )	#for dynamic slots
	
	if walltime is not None:
		#2013.3.21 scale walltime according to clusters_size
		job.addProfile(Profile(Namespace.GLOBUS, key="maxwalltime", value="%s"%(walltime)) )
		#TimeToLive is in seconds
		condorJobRequirementLs.append("(Target.TimeToLive>=%s)"%(int(walltime)*60) )
	#key='requirements' could only be added once for the condor profile
	job.addProfile(Profile(Namespace.CONDOR, key="requirements", value=" && ".join(condorJobRequirementLs) ))

setJobProperRequirement = setJobToProperMemoryRequirement

def registerFile(workflow, filename):
	"""
	2011.12.13
		function to register any file to the workflow.input_site_handler, 
	"""
	file = File(os.path.basename(filename))
	file.addPFN(PFN("file://" + os.path.abspath(filename), \
								workflow.input_site_handler))
	workflow.addFile(file)
	return file

def getAbsPathOutOfExecutable(executable):
	"""
	2012.2.24
		This function extracts path out of a registered executable.
			executable is a registered pegasus executable with PFNs.
	"""
	pfn = (list(executable.pfns)[0])
	#the url looks like "file:///home/crocea/bin/bwa"
	return pfn.url[7:]


def getAbsPathOutOfFile(file):
	"""
	2012.7.25
		call getAbsPathOutOfExecutable
	"""
	return getAbsPathOutOfExecutable(file)

def getExecutableClustersSize(executable=None):
	"""
	2013.03.21, default is None
	"""
	clusters_size = None
	clusteringProf = Profile(Namespace.PEGASUS, key="clusters.size", value="1")
	for profile in executable.profiles:
		if clusteringProf.__hash__() == profile.__hash__():	#__hash__ only involves namespace + key 
			clusters_size = profile.value
	return clusters_size
