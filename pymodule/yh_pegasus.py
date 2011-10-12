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