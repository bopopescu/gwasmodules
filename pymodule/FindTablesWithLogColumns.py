#!/usr/bin/env python
"""
Examples:
	# for postgres tables
	FindTablesWithLogColumns.py -d vervetdb -u yh -k public
	
	# for mysql tables
	# "2>/dev/null" is used to direct stderr to /dev/null so it won't mix with the table names outputed to stdout.
	FindTablesWithLogColumns.py -v mysql -d vervetdb -u yh 2>/dev/null

Description:
	2011-4-14
		It finds all tables in a database (&schema) and check whether each has all
			these columns: created_by, updated_by, date_created, date_updated.
		
		If a table does have all columns, its name will be outputted to stdout;
			otherwise not.
		This program is supposed to feed table names to OutputSQLTrigger.py
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
#	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
#	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
#else:   #32bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

class FindTablesWithLogColumns(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgresql', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['graphdb', 'd', 1, 'stock_250k database name', ],\
						('schema', 0, ): [None, 'k', 1, 'database schema name (postgresql only)', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-07-27
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
	def establishConnection(self):
		"""
		2011-4-14
		"""
		if self.drivername=='mysql':
			import MySQLdb
			conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.db_user, passwd = self.db_passwd)
			cursor = conn.cursor()
		elif self.drivername=='postgresql':
			import psycopg2
			conn_string = "host='%s' dbname='%s' user='%s' password='%s'"%(self.hostname, self.dbname, \
																	self.db_user, self.db_passwd)
			conn = psycopg2.connect(conn_string)
			cursor = conn.cursor()
			if self.schema:
				cursor.execute("set search_path to %s;"%self.schema)
		
		self.conn  = conn
		self.cursor = cursor
	
	def run(self):
		"""
		2011-4-14
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		self.establishConnection()
		
		if self.drivername=='postgresql':
			self.cursor.execute("select tablename from pg_tables where schemaname = '%s'"%(self.schema))
		else:
			self.cursor.execute("show tables")
		rows = self.cursor.fetchall()
		for row in rows:
			tablename = row[0]
			try:
				self.cursor.execute("select created_by, updated_by, date_created, date_updated from %s limit 1"%tablename)
				#self.cursor.fetchall()
				print tablename
			except:	#this table doesn't have these columns.
				self.establishConnection()	# connection/transaction is gone due to exception 
				if self.debug:
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
			

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = FindTablesWithLogColumns
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()