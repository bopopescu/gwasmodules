#!/usr/bin/env python
"""
Examples:
	%s 
	
	# 
	%s -o genomeAnnotationPerBaseForTax60711.tsv.gz
		--genome_drivername=postgresql --genome_hostname=localhost --genome_dbname=vervetdb --genome_schema=genome
		--genome_db_user=yh --genome_db_passwd ...
		--tax_id 60711 --sequence_type_id 1 --chrOrder 1
		--chromosomeHeader chromosome --positionHeader position --whichColumnHeader is_annotated
		--annotation_type_id_list 1,5

Description:
	2013.08.28 a program that outputs genome annotations for each base.
		Output is 4-column matrix file:
			self.chromosomeHeader, self.positionHeader, "genome_annotation_type_id", self.whichColumnHeader
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, utils
from pymodule.yhio.AbstractGenomeFileWalker import AbstractGenomeFileWalker

parentClass = AbstractGenomeFileWalker
class OutputGenomeAnnotation(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	#option_default_dict.update(AbstractMapper.db_option_dict.copy())
	option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.update({
					('annotation_type_id_list', 1, ): [None, '', 1, 'list of annotation type ids'],\
					
					})
	option_default_dict[('whichColumnHeader', 0, )][0] = "is_annotated"
	
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		parentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
		self.annotation_type_id_list = utils.getListOutOfStr(list_in_str=self.annotation_type_id_list, data_type=int)
		self.annotation_type_id_set = set(self.annotation_type_id_list)
	
	def setup(self, **keywords):
		"""
		2013.08.28
		"""
		parentClass.setup(self, **keywords)
		
		header = [self.chromosomeHeader, self.positionHeader, "genome_annotation_type_id", self.whichColumnHeader]
		self.writer.writerow(header)
		
		sys.stderr.write("Fetching all genome annotations with type id in %s ..."%(repr(self.annotation_type_id_list)))
		
		counter = 0
		real_counter = 0
		for chr_id, annot_assembly in self.oneGenomeData._chr_id2annot_assembly.iteritems():
			for genome_annotation in annot_assembly.genome_annotation_list:
				if genome_annotation.genome_annotation_type_id in self.annotation_type_id_set:
					real_counter += 1
					for i in xrange(genome_annotation.start, genome_annotation.stop+1):
						annotation_row = [annot_assembly.chromosome, i, genome_annotation.genome_annotation_type_id, 1]
						self.writer.writerow(annotation_row)
						counter += 1
		sys.stderr.write("%s features affecting %s bases.\n"%(real_counter, counter))
		
	


if __name__ == '__main__':
	main_class = OutputGenomeAnnotation
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()