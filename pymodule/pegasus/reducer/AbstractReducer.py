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
		})
	def __init__(self, inputFnameLs=None, **keywords):
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)