#!/usr/bin/env python
"""
2011-10-11
	module that deals with pegasus-related stuff
"""
import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from ProcessOptions import  ProcessOptions
from Pegasus.DAX3 import *

def addMkDirJob(workflow=None, mkdir=None, outputDir=None, namespace=None, version=None,\
			parentJobLs=[], extraDependentInputLs=[]):
	"""
	2012.3.10
		add argument parentJobLs, extraDependentInputLs
	2011-11-28
		get namespace and version from workflow first
	2011-9-14
	"""
	# Add a mkdir job for any directory.
	mkDirJob = Job(namespace=getattr(workflow, 'namespace', namespace), name=mkdir.name, \
				version=getattr(workflow, 'version', version))
	mkDirJob.addArguments(outputDir)
	mkDirJob.folder = outputDir	#custom attribute
	workflow.addJob(mkDirJob)
	for parentJob in parentJobLs:
		workflow.depends(parent=parentJob, child=mkDirJob)
	for input in extraDependentInputLs:
		mkDirJob.uses(input, transfer=True, register=True, link=Link.INPUT)
	return mkDirJob

def registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, input_site_handler='local',\
						checkAffiliateFileExistence=True, addPicardDictFile=True,\
						affiliateFilenameSuffixLs=['dict', 'fai', 'amb', 'ann', 'bwt', 'pac', 'sa', 'rbwt', 'rpac', 'rsa', \
												'stidx', 'sthash']):
	"""
	2012.5.23
		add an argument "addPicardDictFile" to offer user option to exclude this file (i.e. registerBlastNucleotideDatabaseFile)
	2012.2.24
		dict is via picard
		fai is via "samtools faidx" (index reference)
		amb', 'ann', 'bwt', 'pac', 'sa', 'rbwt', 'rpac', 'rsa' are all bwa index.
		stidx is stampy index.
		sthash is stampy hash.
	2012.2.23
		add two suffixes, stidx (stampy index) and sthash (stampy hash)
	2011-11-11
		if needAffiliatedFiles,
			all other files, with suffix in affiliateFilenameSuffixLs, will be registered (symlinked or copied) as well.
	"""
	refFastaFList = []
	if registerAffiliateFiles:
		refFastaF = File(os.path.basename(refFastaFname))	#use relative path, otherwise, it'll go to absolute path
		# Add it into replica only when needed.
		refFastaF.addPFN(PFN("file://" + refFastaFname, input_site_handler))
		workflow.addFile(refFastaF)
		refFastaFList.append(refFastaF)
		# If it's not needed, assume the index is done and all relevant files are in absolute path.
		# and no replica transfer
		
		#add extra affiliated files
		pathToFileLs = []
		if addPicardDictFile:	#2012.5.23
			dictFname = '%s.dict'%(os.path.splitext(refFastaFname)[0])	#remove ".fasta" from refFastaFname
			pathToFileLs.append(dictFname)
		for suffix in affiliateFilenameSuffixLs:
			pathToFile = '%s.%s'%(refFastaFname, suffix)
			pathToFileLs.append(pathToFile)
		for pathToFile in pathToFileLs:
			if checkAffiliateFileExistence and not os.path.isfile(pathToFile):
				sys.stderr.write("Warning: %s don't exist or not a file on file system. skip registration.\n"%(pathToFile))
				continue
			affiliateF = File(os.path.basename(pathToFile))
			#use relative path, otherwise, it'll go to absolute path
			affiliateF.addPFN(PFN("file://" + pathToFile, input_site_handler))
			workflow.addFile(affiliateF)
			refFastaFList.append(affiliateF)
	else:
		refFastaF = File(refFastaFname)
		refFastaFList.append(refFastaF)
	return refFastaFList

def setJobToProperMemoryRequirement(job, job_max_memory=500, no_of_cpus=1, max_walltime=120, sshDBTunnel=0):
	"""
	2012.4.16
		add argument sshDBTunnel.
			=1: this job needs a ssh tunnel to access psql db on dl324b-1.
			=anything else: no need for that.
	2011-11-23
		set max_walltime default to 120 minutes (2 hours)
	2011-11-16
		add more requirements
	2011-11-11
		job_max_memory is in MB.
		max_walltime is in minutes.
	"""
	job.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%(job_max_memory)))
	job.addProfile(Profile(Namespace.CONDOR, key="request_memory", value="%s"%(job_max_memory)))	#for dynamic slots
	condorJobRequirementLs = ["(memory>=%s)"%(job_max_memory)]
	#2012.4.16
	if sshDBTunnel==1:
		condorJobRequirementLs.append("(sshDBTunnel==%s)"%(sshDBTunnel))	#use ==, not =.
	
	if no_of_cpus is not None:
		job.addProfile(Profile(Namespace.CONDOR, key="request_cpus", value="%s"%(no_of_cpus)) )	#for dynamic slots
	
	if max_walltime is not None:
		job.addProfile(Profile(Namespace.GLOBUS, key="maxwalltime", value="%s"%(max_walltime)) )
		#TimeToLive is in seconds
		condorJobRequirementLs.append("(Target.TimeToLive>=%s)"%(int(max_walltime)*60) )
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