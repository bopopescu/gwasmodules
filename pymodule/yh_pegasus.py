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

def addMkDirJob(workflow, mkdir=None, outputDir=None, namespace=None, version=None):
	"""
	2011-9-14
	"""
	# Add a mkdir job for any directory.
	mkDirJob = Job(namespace=namespace, name=mkdir.name, version=version)
	mkDirJob.addArguments(outputDir)
	mkDirJob.folder = outputDir	#custom attribute
	workflow.addJob(mkDirJob)
	return mkDirJob

def registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, input_site_handler='local',\
						checkAffiliateFileExistence=True,\
						affiliateFilenameSuffixLs=['amb', 'ann', 'bwt', 'dict', 'fai', 'pac', 'rbwt', 'rpac', 'rsa', 'sa']):
	"""
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

def setJobToProperMemoryRequirement(job, job_max_memory=500, no_of_cpus=1, max_walltime=None):
	"""
	2011-11-16
		add more requirements
	2011-11-11
		job_max_memory is in MB.
		max_walltime is in minutes.
	"""
	job.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
	job.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
	job.addProfile(Profile(Namespace.CONDOR, key="request_memory", value="%s"%job_max_memory))	#for dynamic slots
	
	if no_of_cpus is not None:
		job.addProfile(Profile(Namespace.CONDOR, key="request_cpus", value="%s"%no_of_cpus))	#for dynamic slots
	
	if max_walltime is not None:
		job.addProfile(Profile(Namespace.GLOBUS, key="maxwalltime", value="%s"%max_walltime))

setJobProperRequirement = setJobToProperMemoryRequirement
