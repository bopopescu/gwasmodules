#!/usr/bin/env python
"""
Examples:
	#setup database in postgresql
	TaxonomyDB.py -u crocea -k genome
	
	# 2010-12-15 setup genome schema in vervetdb.
	TaxonomyDB.py -u yh -k genome -d vervetdb -v postgresql
	
	#setup database in mysql
	TaxonomyDB.py -v mysql -u yh -z papaya -d taxonomy -k ""
	
	
Description:
	2011-1-20
	This is a wrapper for the taxonomy database (copy of ftp://ftp.ncbi.nih.gov/pub/taxonomy/), built on top of elixir.
"""
import sys, os
from sqlalchemy.engine.url import URL
from elixir import Unicode, DateTime, String, Integer, UnicodeText, Text
from elixir import Entity, Field, using_options, using_table_options
from elixir import OneToMany, ManyToOne, ManyToMany, OneToOne
from elixir import setup_all, session, metadata, entities
from elixir.options import using_table_options_handler	#using_table_options() can only work inside Entity-inherited class.
from datetime import datetime
from sqlalchemy.schema import ThreadLocalMetaData, MetaData
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy import UniqueConstraint, create_engine
from sqlalchemy import and_, or_, not_


from db import ElixirDB

__session__ = scoped_session(sessionmaker(autoflush=False, autocommit=True))
#__metadata__ = ThreadLocalMetaData() #2008-11-04 not good for pylon

__metadata__ = MetaData()

class TaxonomyDB(ElixirDB):
	__doc__ = __doc__
	option_default_dict = ElixirDB.option_default_dict.copy()
	option_default_dict[('drivername', 1,)][0] = 'mysql'
	option_default_dict[('database', 1,)][0] = 'genome'
	def __init__(self, **keywords):
		"""
		2008-10-08
			simplified further by moving db-common lines to ElixirDB
		2008-07-09
		"""
		from ProcessOptions import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		if self.debug:
			import pdb
			pdb.set_trace()
		self.setup_engine(metadata=__metadata__, session=__session__, entities=entities)


if __name__ == '__main__':
	import sys, os, math
	bit_number = math.log(sys.maxint)/math.log(2)
	if bit_number>40:       #64bit
		sys.path.insert(0, os.path.expanduser('~/lib64/python'))
		sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
	else:   #32bit
		sys.path.insert(0, os.path.expanduser('~/lib/python'))
		sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

	from pymodule import ProcessOptions
	main_class = TaxonomyDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.setup()
	
	if instance.debug:
		import pdb
		pdb.set_trace()