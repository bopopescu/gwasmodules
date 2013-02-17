#!/usr/bin/env python
"""

Description:
	2013.2.12 an abstract class for reducer
"""

import sys, os, math
__doc__ = __doc__

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper

class AbstractReducer(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractMapper.option_default_dict)
	option_default_dict.update({
						('exitNonZeroIfAnyInputFileInexistent', 0, ): ['', '', 0, "by default, it skips files that don't exist. Toggling this to exit 3.", ],\
						('noHeader', 0, int): [0, 'n', 0, 'all input has no header'],\
						('inputEmptyType', 0, int): [0, '', 1, 'when to call the input is empty.\n\
	0: completely empty; 1: fastq empty, type 0 + a "@" only line followed by "+"-only line; 2: fasta empty, type 0+ a sequence ID line with only >.'],\
		})
	def __init__(self, inputFnameLs=None, **keywords):
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def isInputLineEmpty(self, inputLine=None, inputFile=None, inputEmptyType=0):
		"""
		2013.2.16
			only inputEmptyType 1 needs to access inputFile
		"""
		if inputLine is None or inputLine=='':
			return True
		content = inputLine.strip()
		if content=='':
			return True
		
		answer = False
		if inputEmptyType==1:
			if content=='@':
				currentPosition = inputFile.tell()
				try:
					nextLine = inputFile.next().strip()
					afterNextLinePosition = inputFile.tell()
					if nextLine.strip()=='+':
						answer = True
					else:
						inputFile.seek(currentPosition-afterNextLinePosition, 1)	#rewind it back from the current position
						"""
						use f.seek(offset, from_what). the reference point is selected by the from_what argument.
						A from_what value of 0 measures from the beginning of the file, 1 uses the current file position,
							 and 2 uses the end of the file as the reference point.
						"""
				except StopIteration:
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
					sys.stderr.write("End of file reached in trying to get '+'.\n")
				
		elif inputEmptyType==2:
			if content=='>':
				answer = True
		return answer