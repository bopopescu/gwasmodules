#!/usr/bin/env python
"""

Description:
	2012.1.17 an abstract class for vcf mapper
"""

import sys, os, math
__doc__ = __doc__

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.yhio.VCFFile import VCFFile
from pymodule.yhio import SNP
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper

class AbstractVCFMapper(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractMapper.option_default_dict)
	option_default_dict.update({
						('inputFname', 0, ): ['', 'i', 1, 'VCF input file. either plain vcf or gzipped is ok. could be unsorted.', ],\
						("chromosome", 0, ): [None, 'c', 1, 'chromosome name for these two VCF.'],\
						("chrLength", 1, int): [1, 'l', 1, 'length of the reference used for the input VCF file.'],\
						('minDepth', 0, float): [0, 'm', 1, 'minimum depth for a call to regarded as non-missing', ],\
						})

	def __init__(self,  inputFnameLs=None, **keywords):
		"""
		"""
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)


if __name__ == '__main__':
	main_class = AbstractVCFMapper
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()