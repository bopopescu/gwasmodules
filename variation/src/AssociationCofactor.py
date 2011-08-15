#!/usr/bin/env python
"""

Examples:
	SNP 4_268809 (FRI deletion 1) and phenotype 44 (FRI) are cofactor. EMMAX association on phenotype 43 FLC.
	 
	~/script/variation/src/AssociationCofactor.py -i /Network/Data/250k/tmp-yh/250k_data/call_method_32.tsv
		-p /Network/Data/250k/tmp-yh/phenotype/phenotype.tsv
		-o /tmp/call_method_32.tsv  -y4 -w 43 -f 4_268809 -n 44
	

Description:
	a program which runs GWA.cofactorLM() from misc.py
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import csv, numpy, traceback
from pymodule import read_data, ProcessOptions, PassingData, SNPData, getListOutOfStr
from numpy import linalg

from Kruskal_Wallis import Kruskal_Wallis
import rpy
from PlotGroupOfSNPs import PlotGroupOfSNPs
from sets import Set
#from DrawEcotypeOnMap import DrawEcotypeOnMap

from misc import GWA
class AssociationCofactor(Kruskal_Wallis):
	debug = 0
	report = 0
	__doc__ = __doc__
	option_default_dict = Kruskal_Wallis.option_default_dict.copy()
	option_default_dict.pop(("which_phenotype", 1, int))
	option_default_dict.update({('phenotype_method_id_ls', 0, ): ['1', 'w', 1, 'which phenotypes to work on. a comma-dash-separated list of phenotype_method ids in the phenotype file. Check db Table phenotype_method. \
		if not available, take all phenotypes in the phenotype_fname.',]})
	option_default_dict.update({('run_type', 1, int): [1, 'y', 1, '1: pure linear model by python, 2: EMMA, 3: pure linear model by R (same as 1), 4: EMMAX']})
	option_default_dict.update({('cofactor_phenotype_id_ls', 0, ): [None, 'n', 1, 'comma-dash-separated list of phenotype_method ids in the phenotype file, to be included as cofactor. like 43-44']})
	option_default_dict.update({('cofactors', 0, ): [None, 'f', 1, 'comma-separated list of SNP IDs to be included as cofactors. like 4_268809,4_269962']})
	def __init__(self, **keywords):
		"""
		2008-11-10
		"""
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		self.phenotype_method_id_ls = getListOutOfStr(self.phenotype_method_id_ls, data_type=int)
		self.cofactors = getListOutOfStr(self.cofactors, data_type=str)
		self.cofactor_phenotype_id_ls = getListOutOfStr(self.cofactor_phenotype_id_ls, data_type=int)
		
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		start_snp = '1_1'
		stop_snp = '5_60000000'
		output_fname_tag_ls = []
		if self.cofactors:
			output_fname_tag_ls.append('cofactor_%s'%('_'.join(self.cofactors)))
		if self.cofactor_phenotype_id_ls:
			cofactor_phenotype_id_str_ls = map(str, self.cofactor_phenotype_id_ls)
			output_fname_tag_ls.append('cofactor_phenotype_%s'%('_'.join(cofactor_phenotype_id_str_ls)))
		output_fname_prefix = '%s_Runtype_%s_%s'%(os.path.splitext(self.output_fname)[0],\
												self.run_type, '_'.join(output_fname_tag_ls))
		GWA.cofactorLM(self.input_fname, self.phenotype_fname, self.phenotype_method_id_ls, output_fname_prefix, start_snp,\
					stop_snp, cofactors=self.cofactors, cofactor_phenotype_id_ls=self.cofactor_phenotype_id_ls, \
					run_type=self.run_type)
	
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = AssociationCofactor
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()