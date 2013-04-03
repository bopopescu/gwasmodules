#!/usr/bin/env python
"""
2008-04-28
A wrapper on top of sqlalchemy around a database. Mostly copied from collective.lead.Database. Can't directly use it because
of trouble in understanding how to use adapter involved in TreadlocalDatabaseTransactions.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import sqlalchemy, threading
from sqlalchemy.engine.url import URL
from sqlalchemy import Table, create_engine
from sqlalchemy.orm import mapper, relation
from sqlalchemy.orm import scoped_session, sessionmaker
from elixir import using_table_options
from pymodule import utils

def supplantFilePathWithNewDataDir(filePath="", oldDataDir=None, newDataDir=None):
	"""
	2012.12.15 argument filePath could be relative path and return join(newDataDir, filePath). oldDataDir is ignored.
		So it could handle both Stock_250kDB (absolute path) and VervetDB (relative path)'s path
		i.e. 
		supplantFilePathWithNewDataDir(filePath='genotype_file/method_31/32093_VCF_26356_VCF_24680_VCF_Contig900.subset.vcf.gz',\
			newDataDir='/Network/Data/vervet/db')
		is same as:
			supplantFilePathWithNewDataDir(filePath='genotype_file/method_31/32093_VCF_26356_VCF_24680_VCF_Contig900.subset.vcf.gz',\
				oldDataDir='/Network/Data/vervet/db', newDataDir='/Network/Data/vervet/db')
		for absolute path:
			supplantFilePathWithNewDataDir(filePath='/Network/Data/250k/db/results/type_1/16_results.tsv',\
				oldDataDir='/Network/Data/250k/db', newDataDir='~/NetworkData/250k/db') 
	2012.11.13
		expose the oldDataDir argument
	2012.3.23
		in case that the whole /Network/Data/250k/db is stored in a different place (=newDataDir)
			how to rescale the filePath ( stored in the database tables) to reflect its new path.
		Stock_250kDB stores absolute path for each file.
		This function helps to adjust that path.
	"""
	if newDataDir:
		newDataDir = os.path.expanduser(newDataDir)
	
	newFilePath = filePath
	if filePath:
		if filePath[0]=='/':	#input is in absolute path.
			if oldDataDir and newDataDir and oldDataDir!=newDataDir:
				if filePath.find(oldDataDir)==0:
					relativePath = filePath[len(oldDataDir):]
					if relativePath[0]=='/':	#2012.12.28 have to remove the initial "/" otherwise os.path.join() won't work
						relativePath = relativePath[1:]
					newFilePath = os.path.join(newDataDir, relativePath)
				else:
					sys.stderr.write("Warning: %s doesn't include old data dir %s. Return Nothing.\n"%(filePath, oldDataDir))
					newFilePath = None
		else:	#relative path
			if newDataDir:
				newFilePath = os.path.join(newDataDir, filePath)
	return newFilePath

class TableClass(object):
	using_table_options(mysql_engine='InnoDB')
	"""
	2008-05-03
		a base class for any class to hold a db table.
		to assign (key, value) to a class corresponding to a table
	"""
	def __init__(self, **keywords):
		for key, value in keywords.iteritems():
			setattr(self, key, value)
	
	def __str__(self):
		"""
		2010-6-17
			a string-formatting function
		"""
		return_ls = []
		for attribute_name in dir(self):
			if attribute_name.find('__')==0:	#ignore the 
				continue
			return_ls.append("%s = %s"%(attribute_name, getattr(self, attribute_name, None)))
		
		return ", ".join(return_ls)

class AbstractTableWithFilename(TableClass):
	"""
	2012.11.13 ancestor of ResultsMethod and AssociationLandscape, and other tables that store paths to files on harddisk.
	"""
	id = None
	short_name = None
	
	path = None
	filename = None
	original_path = None
	md5sum = None	# unique=True
	file_size = None	#2012.7.12
	
	folderName = ''
	
	def getDateStampedFilename(self, oldDataDir=None, newDataDir=None):
		"""
		2012.3.21
			xxx.tsv => xxx.2012_3_21.tsv
		"""
		_filename = self.getFileAbsPath(oldDataDir=oldDataDir, newDataDir=newDataDir)
		
		from datetime import datetime
		lastModDatetime = datetime.fromtimestamp(os.stat(_filename).st_mtime)
		prefix, suffix = os.path.splitext(_filename)
		newFilename = '%s.%s_%s_%s%s'%(prefix, lastModDatetime.year, lastModDatetime.month,\
									lastModDatetime.day, suffix)
		return newFilename
	
	def getFileAbsPath(self, oldDataDir=None, newDataDir=None):
		"""
		2013.1.10 renamed from getFilePath() (which is now another function that returns either self.filename or self.path)
		2012.11.13
			in case that the whole /Network/Data/250k/db is stored in a different place (=data_dir)
				how to modify self.filename (stored in the database tables) to reflect its new path.
		"""
		if self.filename:
			filePath = self.filename
		elif self.path:
			filePath = self.path
		else:
			filePath = None
		return supplantFilePathWithNewDataDir(filePath=filePath, oldDataDir=oldDataDir, newDataDir=newDataDir)
	
	def getFilePath(self):
		"""
		2013.1.10
			returns either self.filename or self.path (old-version tables use self.filename, new version use self.path).
		"""
		if self.filename:
			filePath = self.filename
		elif self.path:
			filePath = self.path
		else:
			filePath = None
		return filePath
	
	def setFilePath(self, newPath=None):
		"""
		2013.1.10
			either filename or path column, depending on which one is a real table column
		"""
		if isinstance(self.__class__.filename, sqlalchemy.orm.attributes.InstrumentedAttribute):
			#type(None) is not good. since null-value column is None as well.
			self.filename = newPath
		elif isinstance(self.__class__.path, sqlalchemy.orm.attributes.InstrumentedAttribute):
			self.path = newPath
		else:
			self.path = newPath
	
	def constructRelativePath(self, data_dir=None, subFolder=None, sourceFilename=None, **keywords):
		"""
		2012.11.13
		"""
		if not subFolder:
			subFolder = self.folderName
		outputDirRelativePath = subFolder
		if data_dir and outputDirRelativePath:
			#make sure the final output folder is created. 
			outputDirAbsPath = os.path.join(data_dir, outputDirRelativePath)
			if not os.path.isdir(outputDirAbsPath):
				os.makedirs(outputDirAbsPath)
		
		filename_part_ls = []
		if self.id:
			filename_part_ls.append(self.id)
		filename_part_ls = map(str, filename_part_ls)
		fileRelativePath = os.path.join(outputDirRelativePath, '%s.h5'%('_'.join(filename_part_ls)))
		return fileRelativePath
	
	def getShortName(self):
		"""
		2012.12.15
		"""
		return None
	
	def setShortName(self):
		"""
		2012.12.28
		"""
		self.short_name = self.getShortName()
		

class Database(object):
	__doc__ = __doc__
	option_default_dict = {('v', 'drivername', 1, '', 1, ):'mysql',\
							('z', 'hostname', 1, '', 1, ):'papaya.usc.edu',\
							('d', 'dbname',1, '', 1, ):None,\
							('k', 'schema',0, '', 1, ):None,\
							('u', 'db_user',1, '', 1, ):None,\
							('p', 'db_passwd',1, '', 1, ):None,\
							('o', 'port', 1, '', 0, ):None,\
							('c', 'commit', 0, '', 0, int):0,\
							('b', 'debug', 0, '', 0, int):0,\
							('r', 'report', 0, '', 0, int):0}
	"""
	2008-02-28
		argument_default_dict is a dictionary of default arguments, the key is a tuple, ('argument_name', is_argument_required, argument_type)
		argument_type is optional
	"""
	def __init__(self, **keywords):
		from pymodule import process_function_arguments, turn_option_default_dict2argument_default_dict
		argument_default_dict = turn_option_default_dict2argument_default_dict(self.option_default_dict)
		self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		self._threadlocal = threading.local()
		self.tables = {}
		self.mappers = {}
		self._engine = None
		
	#@property
	def _url(self):
		return URL(drivername=self.drivername, db_passwd=self.db_user,
				   password=self.db_passwd, host=self.hostname,
				   port=self.port, database=self.dbname)
	_url = property(_url)
	
	def _setup_tables(self, metadata, tables):
		"""Map the database structure to SQLAlchemy Table objects
		"""
		"""
		tables['call_QC'] = Table('call_QC', metadata, autoload=True)
		"""
		pass
	
	def _setup_mappers(self, tables, mappers):
		"""Map the database Tables to SQLAlchemy Mapper objects
		"""
		"""
		mappers['call_QC'] = mapper(Call_QC, tables['call_QC'], properties={'call_info_obj': relation(Call_Info),\
																		'qc_method_obj':relation(QCMethod)})
		"""
		pass
	
	#@property
	def _engine_properties(self):
		return {}
	_engine_properties = property(_engine_properties)
	
	def invalidate(self):
		self._initialize_engine()
		
	# IDatabase implementation - code using (not setting up) the database
	# uses this
	
	#@property
	def session(self):
		"""
		2008-07-09
			use the new sessionmaker() in version 0.4 to create Session
			use scoped_session to create a thread-local context
		"""
		if getattr(self._threadlocal, 'session', None) is None:
			# Without this, we may not have mapped things properly, nor
			# will we necessarily start a transaction when the client
			# code begins to use the session.
			ignore = self.engine
			Session = scoped_session(sessionmaker(autoflush=True, transactional=True, bind=self.engine))
			self._threadlocal.session = Session()
		return self._threadlocal.session
	session = property(session)
	
	#@property
	def connection(self):
		return self.engine.contextual_connect()
	connection = property(connection)
	
	#@property
	def engine(self):
		if self._engine is None:
			self._initialize_engine()
		
		return self._engine
	engine = property(engine)
	
	# Helper methods
	
	def _initialize_engine(self):
		"""
		metadata.create_all()	#2008-07-09 create all tables
		2008-07-09
			close and reset the old session if self._threadlocal.session is not None
			change  self._metadata to metadata in "self.tables[name] = table.tometadata(metadata)"
		2008-07-08
			for postgres, set the schema	(doesn't work)
		"""
		kwargs = dict(self._engine_properties).copy()
		if 'strategy' not in kwargs:
			kwargs['strategy'] = 'threadlocal'
		if 'convert_unicode' not in kwargs:
			kwargs['convert_unicode'] = True
		
		engine = sqlalchemy.create_engine(self._url, **kwargs)
		metadata = sqlalchemy.MetaData(engine)
		
		if getattr(self._threadlocal, 'session', None) is not None:	#2008-07-09 close and reset the old session
			self._threadlocal.session.close()
			self._threadlocal.session = None
		"""
		if getattr(self, 'schema', None):	#2008-07-08 for postgres, set the schema, doesn't help. specify "schema=" in Table()
			con = engine.connect()
			con.execute("set search_path to %s"%self.schema)
			#sys.stderr.write('set schema')
			metadata.bind = con	#necessary. otherwise, schema is still not set.
		"""
		
		# We will only initialize once, but we may rebind metadata if
		# necessary

		if not self.tables:
			self._setup_tables(metadata, self.tables)
			metadata.create_all()	#2008-07-09 create all tables
			self._setup_mappers(self.tables, self.mappers)
		else:
			for name, table in self.tables.items():
				self.tables[name] = table.tometadata(metadata)	#2008-07-09 change self._metadata to metadata
		
		self._engine = engine
		self._metadata = metadata

class ElixirDB(object):
	"""
	2008-11-07
		add option sql_echo
	2008-10-07 add option pool_recycle
	2008-08-07
		expose metadata from elixir
	2008-07-11
		elixir db base class
	"""
	
	option_default_dict = {('drivername', 1,):['postgres', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ):['localhost', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ):[None, 'd', 1, 'database name',],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ):[None, 'u', 1, 'database username',],\
							('db_passwd', 1, ):[None, 'p', 1, 'database password', ],\
							('port', 0, ):[None, '', 1, 'database port number'],\
							('pool_recycle', 0, int):[3600, '', 1, 'the length of time to keep connections open before recycling them.'],\
							('echo_pool', 0, bool):[False, 'e', 0, 'if True, the connection pool will log all checkouts/checkins to the logging stream, which defaults to sys.stdout.'],\
							('sql_echo', 0, ):[False, '', 1, 'wanna echo the underlying sql of every sql query'],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	#from elixir import metadata	#2008-08-07
	#metadata = metadata
	def __init__(self, **keywords):
		"""
		2008-07-09
		"""
		#2013.1.10 backwards compatibility
		keywords = self._setInputArgumentsEquivalentValue(keywords=keywords, argument1Name='username', argument2Name='db_user')
		keywords = self._setInputArgumentsEquivalentValue(keywords=keywords, argument1Name='password', argument2Name='db_passwd')
		keywords = self._setInputArgumentsEquivalentValue(keywords=keywords, argument1Name='database', argument2Name='dbname')
		
		from pymodule import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		if self.echo_pool:	#2010-9-19 passing echo_pool to create_engine() causes error. all pool log disappeared.
			#2010-9-19 Set up a specific logger with our desired output level
			import logging
			#import logging.handlers
			logging.basicConfig()
			#LOG_FILENAME = '/tmp/sqlalchemy_pool_log.out'
			my_logger = logging.getLogger('sqlalchemy.pool')
	
			# Add the log message handler to the logger
			#handler = logging.handlers.RotatingFileHandler(LOG_FILENAME, maxBytes=1000000, backupCount=5)
			# create formatter
			#formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
			# add formatter to handler
			#handler.setFormatter(formatter)
			#my_logger.addHandler(handler)
			my_logger.setLevel(logging.DEBUG)
		
		#self.setup_engine()	#2012.12.18 it needs __metadata__ , __session__ from each db-definition file itself
		self.READMEClass = None	#2012.12.18 required to figure out data_dir
		self._data_dir = None	#2012.11.13
		
		if hasattr(self, 'debug') and self.debug:
			import pdb
			pdb.set_trace()
	
	def _setInputArgumentsEquivalentValue(self, keywords=None, argument1Name=None, argument2Name=None):
		"""
		2013.1.10 to deal with transition that 'username' is now changed to 'db_user' and etc.
		"""
		if argument1Name in keywords and argument2Name not in keywords:
			keywords[argument2Name] = keywords[argument1Name]
		return keywords
		
	def setup_engine(self, metadata=None, session=None, entities=[]):
		"""
		2008-10-09
			a common theme for all other databases
		"""
		from elixir.options import using_table_options_handler
		if getattr(self, 'schema', None):	#for postgres
			for entity in entities:
				if entity.__module__==self.__module__:	#entity in the same module
					using_table_options_handler(entity, schema=self.schema)	#2012.5.16 this changes the entity._descriptor.table_fullname
					#if self.schema!='public':
					#	entity._descriptor.tablename = '%s.%s'%(self.schema, entity._descriptor.tablename)	#2012.5.16 change the tablename
		#2008-10-05 MySQL typically close connections after 8 hours resulting in a "MySQL server has gone away" error.
		metadata.bind = create_engine(self._url, pool_recycle=self.pool_recycle, echo=self.sql_echo)	#, convert_unicode=True, encoding="utf8")
		self.metadata = metadata
		self.session = session
		#2012.12.28
		self.cleanUpMetadatas()
	
	def cleanUpMetadatas(self):
		"""
		2012.12.28 because GenomeDB, TaxonomyDB are by default imported into the namespace through pymodule all the time
			their unbound metadata is automatically added into elixir's metadatas variable.
			which caused this error if create_tables=True:
			
				  File "/home/crocea/script/variation/src/db/Stock_250kDB.py", line 2582, in setup
				    setup_all(create_tables=create_tables)      #create_tables=True causes setup_all to call elixir.create_all(), which in turn calls me
				tadata.create_all()
				  File "/usr/local/lib/python2.7/dist-packages/elixir/__init__.py", line 98, in setup_all
				    create_all(*args, **kwargs)
				  File "/usr/local/lib/python2.7/dist-packages/elixir/__init__.py", line 76, in create_all
				    md.create_all(*args, **kwargs)
				  File "/usr/local/lib/python2.7/dist-packages/sqlalchemy/schema.py", line 2560, in create_all
				    bind = _bind_or_error(self)
				  File "/usr/local/lib/python2.7/dist-packages/sqlalchemy/schema.py", line 3176, in _bind_or_error
				    raise exc.UnboundExecutionError(msg)
				sqlalchemy.exc.UnboundExecutionError: The MetaData is not bound to an Engine or Connection.  Execution can not proceed without a databas
				e to execute against.  Either execute with an explicit connection or assign the MetaData's .bind to enable implicit execution.
		"""
		import elixir
		newMetadatas = set()
		for md in elixir.metadatas:
			if md.bind is not None:
				newMetadatas.add(md)
		elixir.metadatas = newMetadatas
	
	def setup(self, create_tables=True):
		"""
		2008-10-09
			add option create_tables
		2008-08-26
		"""
		from elixir import setup_all
		setup_all(create_tables=create_tables)	#create_tables=True causes setup_all to call elixir.create_all(),
				#which in turn calls metadata.create_all()
	
	def _url(self):
		return URL(drivername=self.drivername, username=self.db_user,
				password=self.db_passwd, host=self.hostname,
				port=self.port, database=self.dbname)
	_url = property(_url)
	
	@property
	def data_dir(self):
		"""
		2012.3.23
			(learnt from VervetDB)
			get the master directory in which all files attached to this db are stored.
		"""
		if not self._data_dir:
			if self.READMEClass:
				data_dirEntry = self.READMEClass.query.filter_by(title='data_dir').first()
				if not data_dirEntry or not data_dirEntry.description:
					# todo: need to test data_dirEntry.description is writable to the user
					sys.stderr.write("data_dir not available in db or not accessible on the harddisk. Raise exception.\n")
					raise
					self._data_dir = None
				else:
					self._data_dir = data_dirEntry.description
		return self._data_dir
	
	def checkIfEntryInTable(self, TableClass=None, short_name=None, id=None):
		"""
		2013.04.03 bugfix. query = query.filter_by...
			and check how many entries return from db query. should be one. otherwise raise exception.
		2013.3.14
			this could be used as generic way to query tables with short_name (unique) & id columns
		"""
		query = TableClass.query
		if short_name or id:
			if short_name:
				query = query.filter_by(short_name=short_name)
			if id is not None:
				db_entry = TableClass.get(id)
				return db_entry
		else:
			sys.stderr.write("Either short_name (%s) or id (%s) have to be non-None.\n"%(short_name, id))
			raise
		db_entry = query.first()
		no_of_entries = query.count()
		if no_of_entries>1:
			sys.stderr.write("Error, query table %s by short_name=%s, id=%s returns %s entries (>1).\n"%\
							(TableClass, short_name, id, no_of_entries))
			raise
		if db_entry:
			return db_entry
		else:
			return None
	
	def getProperTableName(self, tableClass=None):
		"""
		2012.12.31 helper function to add schema in front of tablename.
		"""
		tableName = tableClass.table.name
		return self.appendSchemaToTableName(tableName)
	
	def appendSchemaToTableName(self, tableName=None):
		"""
		2012.12.31 helper function to add schema in front of tablename.
		"""
		if self.drivername.find('postgres')>=0 and self.schema:
			tableName = '%s.%s'%(self.schema, tableName)
		return tableName
	
	def updateDBEntryMD5SUM(self, db_entry=None, data_dir=None, absPath=None):
		"""
		2012.12.15 moved from VervetDB
		2012.7.13
			if absPath is given, take that , rather than construct it from data_dir and db_entry.path
		"""
		from pymodule import utils
		if data_dir is None:
			data_dir = self.data_dir
		
		if hasattr(db_entry, 'path') and db_entry.path:
			db_entry_path = db_entry.path
		elif hasattr(db_entry, 'filename') and db_entry.filename:
			db_entry_path = db_entry.filename
		else:
			db_entry_path = None
		if not absPath and db_entry_path:
			absPath = supplantFilePathWithNewDataDir(filePath=db_entry_path, oldDataDir=self.data_dir,\
													newDataDir=data_dir)
		
		if absPath and not os.path.isfile(absPath):
			sys.stderr.write("Warning: file %s doesn't exist.\n"%(absPath))
			return
		md5sum = utils.get_md5sum(absPath)
		if db_entry.md5sum is  not None and db_entry.md5sum!=md5sum:
			sys.stderr.write("WARNING: The new md5sum %s is not same as the existing md5sum %s.\n"%(md5sum, db_entry.md5sum))
		db_entry.md5sum = md5sum
		self.session.add(db_entry)
		self.session.flush()
	"""
	#2012.7.14 update the md5sum for the existing db entries
		TableClass = VervetDB.IndividualAlignment
		TableClass = VervetDB.IndividualSequenceFile
		for db_entry in TableClass.query:
			if db_entry.md5sum is None:
				absPath = os.path.join(db_vervet.data_dir, db_entry.path)
				if db_entry.path and os.path.isfile(absPath):
					sys.stderr.write("md5sum on %s ... "%(db_entry.path))
					db_vervet.updateDBEntryMD5SUM(db_entry=db_entry, absPath=absPath)
					sys.stderr.write("\n")
		db_vervet.session.flush()
		db_vervet.session.commit()
		sys.exit(0)
	"""
	def updateDBEntryPathFileSize(self, db_entry=None, data_dir=None, absPath=None):
		"""
		2012.12.15 moved from VervetDB
		2012.7.13
			if absPath is given, take that , rather than construct it from data_dir and db_entry.path
		"""
		from pymodule import utils
		if data_dir is None:
			data_dir = self.data_dir
		
		if hasattr(db_entry, 'path') and db_entry.path:
			db_entry_path = db_entry.path
		elif hasattr(db_entry, 'filename') and db_entry.filename:
			db_entry_path = db_entry.filename
		else:
			db_entry_path = None
		if not absPath and db_entry_path:
			absPath = supplantFilePathWithNewDataDir(filePath=db_entry_path, oldDataDir=self.data_dir,\
													newDataDir=data_dir)
			#absPath = os.path.join(data_dir, db_entry.path)
		if absPath and not os.path.isfile(absPath):
			sys.stderr.write("Warning: file %s doesn't exist.\n"%(absPath))
			return
		file_size = utils.getFileOrFolderSize(absPath)
		if db_entry.file_size is not None and file_size!=db_entry.file_size:
			sys.stderr.write("Warning: the new file size %s doesn't match the old one %s.\n"%(file_size, db_entry.file_size))
		db_entry.file_size = file_size
		self.session.add(db_entry)
		self.session.flush()
	
	def copyFileWithAnotherFilePrefix(self, inputFname=None, filenameWithPrefix=None, outputDir=None,\
									logMessage=None, srcFilenameLs=None, dstFilenameLs=None):
		"""
		2013.3.18 bugfix in filename. there was extra . between prefix and suffix.
			moved from vervet/src/VervetDB.py
		2012.9.20
		"""
		srcFilename = inputFname
		prefix, suffix = os.path.splitext(os.path.basename(inputFname))
		newPrefix = os.path.splitext(filenameWithPrefix)[0]
		dstFilename = os.path.join(outputDir, '%s%s'%(newPrefix, suffix))
		utils.copyFile(srcFilename=srcFilename, dstFilename=dstFilename)
		if logMessage:
			logMessage += "file %s has been copied to %s.\n"%(srcFilename, dstFilename)
		if srcFilenameLs:
			srcFilenameLs.append(srcFilename)
		if dstFilenameLs:
			dstFilenameLs.append(dstFilename)
		return logMessage
	
	def moveFileIntoDBAffiliatedStorage(self, db_entry=None, filename=None, inputDir=None, outputDir=None, \
									dstFilename=None,\
								relativeOutputDir=None, shellCommand='cp -rL', srcFilenameLs=None, dstFilenameLs=None,\
								constructRelativePathFunction=None, data_dir=None):
		"""
			filename (required): relative path of input file
			inputDir (required): where 'filename' is from
			outputDir (required): where the output file will be 
			dstFilename: the absolute path of where the output file will be.
				if set to None (usually), then it'll be constructed on the fly. First 
					either through constructRelativePathFunction()
					or use join(relativeOutputDir, '%s_%s'%(db_entry.id, filename))
					or '%s_%s'%(db_entry.id, filename)
			
			relativeOutputDir: used for construct dstFilename if constructRelativePathFunction() is not there.
			constructRelativePathFunction: similar function of relativeOutputDir.
			 	used to construct relative path of output file.
			if neither relativeOutputDir nor constructRelativePathFunction is available, relative path is ='%s_%s'%(db_entry.id, filename).
				relative path is used to set db_entry.path when the latter is None.
			
			srcFilenameLs, dstFilenameLs: optional. two lists used to store the absolute path of input and output files.
				used in case rollback is needed.
			
			data_dir: the top-level folder where all the db-affiliated file storage is. for constructRelativePathFunction 
			 	
		2013.1.31 bugfix: if relativeOutputDir is included in both outputDir and newPath, use newfilename to avoid double usage. 
		2012.12.15 moved from VervetDB. i.e.:
			inputFileBasename = os.path.basename(self.inputFname)
			relativePath = genotypeFile.constructRelativePath(sourceFilename=inputFileBasename)
			exitCode = self.db_vervet.moveFileIntoDBAffiliatedStorage(db_entry=genotypeFile, filename=inputFileBasename, \
									inputDir=os.path.split(self.inputFname)[0], dstFilename=os.path.join(self.data_dir, relativePath), \
									relativeOutputDir=None, shellCommand='cp -rL', \
									srcFilenameLs=self.srcFilenameLs, dstFilenameLs=self.dstFilenameLs,\
									constructRelativePathFunction=genotypeFile.constructRelativePath, data_dir=self.data_dir)
			#same as this
			exitCode = self.db_vervet.moveFileIntoDBAffiliatedStorage(db_entry=genotypeFile, filename=inputFileBasename, \
									inputDir=os.path.split(self.inputFname)[0], \
									outputDir=self.data_dir, \
									relativeOutputDir=None, shellCommand='cp -rL', \
									srcFilenameLs=self.srcFilenameLs, dstFilenameLs=self.dstFilenameLs,\
									constructRelativePathFunction=genotypeFile.constructRelativePath, data_dir=self.data_dir)
									
			if exitCode!=0:
				sys.stderr.write("Error: moveFileIntoDBAffiliatedStorage() exits with %s code.\n"%(exitCode))
				session.rollback()
				self.cleanUpAndExitOnFailure(exitCode=exitCode)
		
		2012.8.30 add argument dstFilename, which if given , overwrites outputDir
		2012.7.18 -L of cp meant "always follow symbolic links in SOURCE".
		2012.7.13 copied from RegisterAndMoveSplitSequenceFiles.moveNewISQFileIntoDBStorage()
			filename could be a folder.
		2012.7.4
			add srcFilename and dstFilename into given arguments (srcFilenameLs, dstFilenameLs) for later undo
		2012.6.8
			return non-zero if failure in move or destination file already exists
		2012.2.10
			this function moves a file to a db-affiliated storage path
			relativeOutputDir is the path part (in relative path) of db_entry.path = os.path.split(db_entry.path)[0]
		"""
		exitCode = 0
		if constructRelativePathFunction is not None:
			newPath = constructRelativePathFunction(db_entry=db_entry, sourceFilename=filename, data_dir=data_dir)
			newfilename = os.path.basename(newPath)
		elif relativeOutputDir:
			newfilename = '%s_%s'%(db_entry.id, filename)
			newPath = os.path.join(relativeOutputDir, newfilename)
		else:
			newfilename = '%s_%s'%(db_entry.id, filename)
			newPath = newfilename
		
		if db_entry.getFilePath()!=newPath:
			db_entry.setFilePath(newPath)
			try:
				self.session.add(db_entry)
				self.session.flush()
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
				exitCode = 4
				return exitCode
		
		srcFilename = os.path.join(inputDir, filename)
		if dstFilename is None:	#2012.8.30
			if relativeOutputDir:
				relativePathIndex = outputDir.find(relativeOutputDir)
				noOfCharsInRelativeOutputDir = len(relativeOutputDir)
				if outputDir[relativePathIndex:relativePathIndex+noOfCharsInRelativeOutputDir]==relativeOutputDir and newPath.find(relativeOutputDir)>=0:
					#2013.1.31 bugfix: if relativeOutputDir is included in both outputDir and newPath, use newfilename to avoid double usage. 
					dstFilename = os.path.join(outputDir, newfilename)
			if dstFilename is None:	#still nothing , use newPath instead
				dstFilename = os.path.join(outputDir, newPath)
		if os.path.isfile(dstFilename):
			sys.stderr.write("Error: destination %s already exists.\n"%(dstFilename))
			exitCode = 2
		else:
			#21012.12.15 create folder if not existent
			dstFolder = os.path.split(dstFilename)[0]
			if not os.path.isdir(dstFolder):
				os.makedirs(dstFolder)
			#move the file
			commandline = '%s %s %s'%(shellCommand, srcFilename, dstFilename)
			return_data = utils.runLocalCommand(commandline, report_stderr=True, report_stdout=True)
			if srcFilenameLs is not None:
				srcFilenameLs.append(srcFilename)
			if dstFilenameLs is not None:
				dstFilenameLs.append(dstFilename)
			if hasattr(db_entry, 'md5sum') and getattr(db_entry, 'md5sum', None) is None:	#2012.7.14 has this attribute but it's None
				try:
					self.updateDBEntryMD5SUM(db_entry=db_entry, absPath=dstFilename)
				except:
					self.session.delete(db_entry)
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
					exitCode = 2
					return exitCode
			if return_data.stderr_content:
				#something wrong. abort
				sys.stderr.write("commandline %s failed: %s\n"%(commandline, return_data.stderr_content))
				#remove the db entry
				self.session.delete(db_entry)
				self.session.flush()
				exitCode = 3
				return exitCode
			if hasattr(db_entry, 'file_size') and db_entry.file_size is None:
				try:
					self.updateDBEntryPathFileSize(db_entry=db_entry, absPath=dstFilename)
				except:
					self.session.delete(db_entry)
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
					exitCode = 2
					return exitCode
			else:
				exitCode = 0
		return exitCode
	
	def reScalePathByNewDataDir(self, filePath=None, newDataDir=None):
		"""
		2013.1.11 moved from Stock_250kDB
			if newDataDir is None, take self.data_dir
		2012.11.13
		"""
		oldDataDir=self.data_dir
		if newDataDir is None:
			newDataDir=self.data_dir
		return supplantFilePathWithNewDataDir(filePath=filePath, oldDataDir=oldDataDir, newDataDir=newDataDir)
	
	"""
	def session(self):
		#2008-07-09
		from elixir import session
		#session = scoped_session(sessionmaker(autoflush=True, transactional=True, bind=self._url))
		return session
	
	session = property(session)
	
	def connection(self):
		return self.session.connection
	connection = property(connection)
	"""

"""
from elixir import Entity, Field, using_options, using_table_options
from elixir import DateTime, String
from datetime import datetime
class README(Entity):
	#2008-08-07
	title = Field(String(2000))
	description = Field(String(60000))
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='readme')
	using_table_options(mysql_engine='InnoDB')
"""

def formReadmeObj(argv, ad, READMEClass):
	"""
	2008-08-26
		if argument contains 'password' or 'passwd', mark its value as '*'
	2008-05-06
		create a readme instance (like a log) based on the program's sys.argv and argument dictionary
	"""
	readme_description_ls = []
	argument_ls = ad.keys()
	argument_ls.sort()
	for argument in argument_ls:
		if argument.find('passwd')!=-1 or argument.find('password')!=-1:	#password info not into db
			value = '*'
		else:
			value = ad[argument]
		readme_description_ls.append('%s=%s'%(argument, value))
	readme = READMEClass(title=' '.join(argv), description='; '.join(readme_description_ls))
	return readme

def db_connect(hostname, dbname, schema=None, password=None, user=None):
	"""
	2008-07-29
		copied from annot.bin.codense.common
		add the code to deal with importing psycopg
	2007-03-07
		add password and user two options
	02-28-05
		establish database connection, return (conn, curs).
		copied from CrackSplat.py
	03-08-05
		parameter schema is optional
	"""
	connection_string = 'host=%s dbname=%s'%(hostname, dbname)
	if password:
		connection_string += ' password=%s'%password
	if user:
		connection_string += ' user=%s'%user
	
	try:
		import psycopg
	except ImportError:
		try:
			import psycopg2 as psycopg
		except ImportError:
			sys.stderr.write("Neither psycopg nor psycopg2 is installed.\n")
			raise
	conn = psycopg.connect(connection_string)
	curs = conn.cursor()
	if schema:
		curs.execute("set search_path to %s"%schema)
	return (conn, curs)


def get_sequence_segment(curs, gi=None, start=None, stop=None, annot_assembly_id=None, \
						annot_assembly_table='sequence.annot_assembly', \
						raw_sequence_table='sequence.raw_sequence', chunk_size=10000):
	"""
	2011-10-25
		add argument annot_assembly_id
		Primary key of table annot_assembly is a self-contained "id". "gi" is no longer a primary key,
			although it's still kept. For some entries, it could be null.
		start and stop are 1-based.
	2010-10-05 if this AnnotAssembly is not associated with any raw sequence (raw_sequence_start_id is None).
		return ''
	2009-01-03
		moved from annot/bin/codense/common.py
		curs could be elixirdb.metadata.bind other than the raw curs from psycopg
		raw_sequence_table replaced column acc_ver with annot_assembly_gi
		make sure no_of_chunks_before >=0. old expression becomes negative for genes that are near the tip of a chromosome, within chunk_size.
	11-12-05 get a specified segment from chromosome sequence table
		reverse is handled but complement(strand) is not. Upper level function should take care of this.
	11-15-05 improve it to be more robust, add acc_ver and report if not found in raw_sequence_table
	2006-08-28 fix a bug when the start%chunk_size =0 (sits on the edge of the previous chunk, two more chunks
		are required)
	"""
	need_reverse = int(start>stop)
	if need_reverse:
		start, stop = stop, start
	if annot_assembly_id:	#2011-10-25
		annot_assembly_key_name = "id"
		annot_assembly_key_value = annot_assembly_id
	else:
		annot_assembly_key_name = "gi"
		annot_assembly_key_value = gi
	rows = curs.execute("select acc_ver, start, stop, raw_sequence_start_id from %s where %s=%s"%(annot_assembly_table, \
															annot_assembly_key_name, annot_assembly_key_value))
	is_elixirdb = 1
	if hasattr(curs, 'fetchall'):	#2009-01-03 this curs is not elixirdb.metadata.bind
		rows = curs.fetchall()
		is_elixirdb = 0
	if is_elixirdb:	#2009-01-03
		row = rows.fetchone()
		acc_ver = row.acc_ver
		orig_start = row.start
		orig_stop = row.stop
		raw_sequence_start_id = row.raw_sequence_start_id
	else:
		acc_ver, orig_start, orig_stop, raw_sequence_start_id = rows[0]
	if stop>orig_stop:	#11-14-05 to avoid exceeding the boundary
		stop = orig_stop
	
	if raw_sequence_start_id is None:	# 2010-10-05 this AnnotAssembly is not associated with any raw sequence.
		return ''
	
	no_of_chunks_before = max(0, start/chunk_size-1)	#how many chunks are before this segment (2006-08-28) -1 ensures the edge.
		#2008-01-03 max(0, ...) to make sure it >=0. old expression becomes negative for genes that are near the tip of a chromosome, within chunk_size.
	
	segment_size = stop - start +1
	no_of_chunks_segment = segment_size/chunk_size + 1	#how many chunks could be included in this segment
	raw_sequence_start_id += no_of_chunks_before	#the first chunk which contains teh segment
	offset = no_of_chunks_segment + 2	#add two more chunks to ensure the segment is enclosed(2006-08-28)
	#get the sequence from raw_sequence_table
	seq = ''
	for i in range(offset):
		rows = curs.execute("select sequence from %s where annot_assembly_%s=%s and id=%s"%\
				(raw_sequence_table, annot_assembly_key_name, annot_assembly_key_value, raw_sequence_start_id+i))
		if is_elixirdb:	#2009-01-03
			rows = rows.fetchone()
		else:
			rows = curs.fetchall()
		
		if rows:	#11/14/05 it's possible to exceed the whole raw_sequence table because the offset adds one more chunk
			if is_elixirdb:
				seq += rows.sequence
			else:
				seq += rows[0][0]
		else:
			sys.stderr.write("id %s missing in raw_sequence table.\n"%(raw_sequence_start_id+i))
			sys.stderr.write("%s: %s, start: %s, stop: %s, raw_sequence_start_id: %s\n"%\
							(annot_assembly_key_name, \
							annot_assembly_key_value, start, stop, raw_sequence_start_id))
	relative_start = start - no_of_chunks_before*chunk_size
	segment = seq[relative_start-1:relative_start-1+segment_size]	#WATCH index needs -1
	if need_reverse:
		segment = list(segment)
		segment.reverse()	#only 
		segment = ''.join(segment)
	return segment

if __name__ == '__main__':
	from pymodule import process_options, generate_program_doc
	main_class = Database
	opts_dict = process_options(sys.argv, main_class.option_default_dict, error_doc=generate_program_doc(sys.argv[0], main_class.option_default_dict)+main_class.__doc__)
	
	"""
	
	#from /usr/lib/python2.5/site-packages/zope/interface/adapter.txt
	from zope.interface.adapter import AdapterRegistry
	registry = AdapterRegistry()
	registry.register([Database], ITransactionAware, '', 12)
	
	from zope.component import registry
	from zope.component import tests
	components = registry.Components('comps')
	
	#from /usr/lib/zope2.10/lib/python/zope/, interface/ component/registry.txt
	print components.registerAdapter(TreadlocalDatabaseTransactions)
	print components.getAdapter(ITransactionAware, Database)
	"""
	instance = main_class(**opts_dict)
	if instance.debug:
		import pdb
		pdb.set_trace()
	session = instance.session
	"""
	print dir(session.query(Results))
	
	for row in session.query(ResultsMethod).list():
		print row.id
		print row.short_name
	
	i = 0
	while i <10:
		row = session.query(Results).offset(i).limit(1).list()	#all() = list() returns a list of objects. first() returns the 1st object. one() woud raise error because 'Multiple rows returned for one()'
		print len(row)
		row = row[0]
		i += 1
		print row.id
		print row.chr
		print row.start_pos
		print row.score
		print row.method_id
		print row.results_method_obj.short_name
		print row.phenotype_method_id
		print row.phenotype_method_obj.short_name
	"""