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
						('minDepth', 1, float): [1, 'm', 1, 'minimum depth for a call to regarded as non-missing', ],\
						})
	def __init__(self,  **keywords):
		"""
		2012.1.17
		"""
		AbstractNGSWorkflow.__init__(self, **keywords)
		
	def registerExecutables(self, workflow):
		"""
		"""
		AbstractNGSWorkflow.registerExecutables(self, workflow)
	
	def registerCommonExecutables(self, workflow):
		"""
		"""
		AbstractNGSWorkflow.registerCommonExecutables(self, workflow)