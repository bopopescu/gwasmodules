#!/usr/bin/env python
"""

Examples:
	OutputPhenotype.py -o /tmp/phenotype.tsv -z localhost -u yh
	
	OutputPhenotype.py -o /tmp/phenotype.tsv -e stock.ecotype_usc --getPublicPhenotype
	
	#get raw (un-transformed) phenotype
	OutputPhenotype.py -o /tmp/phenotype_g.tsv -g

Description:
	program to output phenotype_avg table.
	The output format is roughly a ecotype_id X phenotype(shortname) matrix.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, getopt
import traceback
from pymodule import process_function_arguments, write_data_matrix, PassingData, SNPData
from variation.src import AbstractVariationMapper

class OutputPhenotype(AbstractVariationMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVariationMapper.option_default_dict.copy()
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update( {
					('phenotype_avg_table',1, ):['stock_250k.phenotype_avg', '', 1,  ],\
					('phenotype_method_table',1, ):['stock_250k.phenotype_method', '', 1, ],\
					('ecotype_table', 1, ): ['stock.ecotype', '', 1, 'ecotype table to get name related to each ecotype', ],\
					('getPublicPhenotype', 0, int):[0, '', 0, 'toggle to get public phenotype only, phenotype_method.access==public and phenotype_avg.ready_for_publication=1'],\
					('get_raw_data', 0, int):[0, '', 0, 'whether to output raw phenotype data from db or transform according to column transformation_description in phenotype_method table'],\
					})
	def __init__(self,  **keywords):
		"""
		2008-11-10
			upgrade option handling to ProcessOptions
		2008-4-2
		2008-02-28
			argument_default_dict is a dictionary of default arguments, the key is a tuple, ('argument_name', is_argument_required, argument_type)
			argument_type is optional
		"""
		#argument dictionary
		#self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=__doc__, class_to_have_attr=self)
		AbstractVariationMapper.__init__(self, inputFnameLs=None, **keywords)
		#from pymodule import ProcessOptions
		#self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
	@classmethod
	def get_phenotype_method_id_info(cls, curs=None, phenotype_avg_table=None, phenotype_method_table=None, getPublicPhenotype=False):
		"""
		2012.9.28
			add argument getPublicPhenotype
		2009-2-2
			curs could be either MySQLdb cursor or elixirdb.metadata.bind.
			do two selects in one
			
		2008-4-2
		"""
		sys.stderr.write("Getting phenotype_method_id info ... " )
		phenotype_method_id2index = {}	#index of the matrix
		method_id_name_ls = []	#as header for each phenotype
		phenotype_id_ls = []
		subQueryToGetDistinctPhenotypeMethodID = "select distinct method_id from %s"%(phenotype_avg_table)
		if getPublicPhenotype:
			subQueryToGetDistinctPhenotypeMethodID += " where ready_for_publication=1"
		phenotypeMethodQuery = "select m.id, m.short_name, m.transformation_description from %s m, (%s) p where m.id=p.method_id "%\
				(phenotype_method_table, subQueryToGetDistinctPhenotypeMethodID)
		if getPublicPhenotype:
			extraConditionSQL = "and m.access='public'"
		else:
			extraConditionSQL = ""
		orderPhenotypeSQL = " order by id "
		rows = curs.execute(" %s %s %s"%(phenotypeMethodQuery, extraConditionSQL,  orderPhenotypeSQL))
		
		is_elixirdb = 1
		if hasattr(curs, 'fetchall'):	#2009-2-2 this curs is not elixirdb.metadata.bind
			rows = curs.fetchall()
			is_elixirdb = 0
		phenotype_method_id2transformation_description = {}
		for row in rows:
			if is_elixirdb:
				method_id = row.id
				method_short_name = row.short_name
				transformation_description = row.transformation_description
			else:
				method_id, method_short_name, transformation_description = row[:3]
			"""
			curs.execute("select short_name, transformation_description from %s where id=%s"%(phenotype_method_table, method_id))
			pm_rows = curs.fetchall()
			method_short_name = pm_rows[0][0]
			transformation_description = pm_rows[0][1]
			"""
			phenotype_id_ls.append(method_id)
			method_id_name_ls.append('%s_%s'%(method_id, method_short_name))
			phenotype_method_id2index[method_id] = len(phenotype_method_id2index)
			if transformation_description=='None':
				transformation_description = None
			phenotype_method_id2transformation_description[method_id] = transformation_description
		return_data = PassingData(phenotype_method_id2index=phenotype_method_id2index, method_id_name_ls=method_id_name_ls,\
								phenotype_id_ls=phenotype_id_ls,\
								phenotype_method_id2transformation_description=phenotype_method_id2transformation_description)
		sys.stderr.write("Done\n")
		return return_data
	
	@classmethod
	def extractPhenotypeIDFromMethodIDName(cls, method_id_name=None):
		"""
		2010-4-21
			split out of PlotGroupOfSNPs.findOutWhichPhenotypeColumn()
			
			extract an integer phenotype method id out of method_id_name
		"""
		if type(method_id_name)==int or type(method_id_name)==long:
			# 2009-4-30 phenData from HaplotypeView (pylons web framework use integer as id)
			phenotype_method_id = method_id_name
		elif type(method_id_name)==str:
			method_id_name = method_id_name.split('_')
			phenotype_method_id=int(method_id_name[0])
		else:	# extreme case
			phenotype_method_id = method_id_name
		return phenotype_method_id
	
	@classmethod
	def get_ecotype_id2info(cls, curs=None, phenotype_avg_table=None, ecotype_table=None, getPublicPhenotype=False):
		"""
		2013.1.3 #in the mysql stock_250k, nativename is of str type.
			#in psql, nativename is of 'unicode' type. In output, (write_data_matrix), python uses the default encoder ascii
			# to convert nativename to str type. Error occurs. So here encodes manually through 'utf-8'
		2012.9.28
			add argument getPublicPhenotype, 
		2009-2-2
			curs could be either MySQLdb cursor or elixirdb.metadata.bind.
			do two selects in one
		
		2008-4-2
		"""
		sys.stderr.write("Getting ecotype id info ... " )
		ecotype_id2index = {}	#index of the matrix
		ecotype_id_ls = []
		ecotype_name_ls = []
		subQueryToGetDistinctPhenotypeMethodID = "select distinct ecotype_id from %s"%(phenotype_avg_table)
		if getPublicPhenotype:
			subQueryToGetDistinctPhenotypeMethodID += " where ready_for_publication=1 "
		rows = curs.execute("select e.id, e.nativename from %s e, ( %s) p where e.id=p.ecotype_id order by id"%\
					(ecotype_table, subQueryToGetDistinctPhenotypeMethodID))
		is_elixirdb = 1
		if hasattr(curs, 'fetchall'):	#2009-2-2 this curs is not elixirdb.metadata.bind
			rows = curs.fetchall()
			is_elixirdb = 0
		for row in rows:
			if is_elixirdb:
				ecotype_id = row.id
				nativename = row.nativename
			else:
				ecotype_id, nativename = row[:2]
			"""
			curs.execute("select nativename from %s where id=%s"%(ecotype_table, ecotype_id))
			nativename = curs.fetchall()[0][0]
			"""
			ecotype_name_ls.append(nativename.encode('utf-8'))
			#in the mysql stock_250k, nativename is of str type.
			#in psql, nativename is of 'unicode' type. In output, (write_data_matrix), python uses the default encoder ascii
			# to convert nativename to str type. Error occurs. So here encodes manually through 'utf-8'
			ecotype_id_ls.append(ecotype_id)
			ecotype_id2index[ecotype_id] = len(ecotype_id2index)
		sys.stderr.write("Done\n")
		return ecotype_id2index, ecotype_id_ls, ecotype_name_ls
	
	@classmethod
	def get_matrix(cls, curs=None, phenotype_avg_table=None, ecotype_id2index=None, phenotype_info=None, get_raw_data=0, \
				phenotype_method_table='phenotype_method', getPublicPhenotype=False):
		"""
		2012.9.28
			add argument getPublicPhenotype
		2010-1-26
			Comment out the code inserted on 2009-9-2 (right below), since its purpose is to avoid truncation and 
				the truncation wasn't due to the value being too small.
			It was that the conversion of a 2D list containing character 'NA' into a numpy array renders the whole 2D list
				being converted into a character numpy array. The conversion was default in write_data_matrix() if
				transform_to_numpy=False is not passed on.  
		2009-9-2
			if value>-5e-7 and value<+5e-7:	#beyond float resolution by a python float
				value = 0
			
			without condition above, values like -5.32907e-15 would be taken as -5.32907e, -3.76545e-12 as -3.76545
			
		2009-9-2
			add phenotype_method_table to get stddev, min_value to do certain transformation involving these two variables
		2009-2-2
			curs could be either MySQLdb cursor or elixirdb.metadata.bind.
			average phenotype values among replicates in the same phenotype method
		2008-11-10
			add code to transform phenotype according to phenotype_info.phenotype_method_id2transformation_description
			add option get_raw_data, if True/1, no transformation.
		2008-04-23
			#some db entries (phenotype_avg.value) have nothing there. convert None to 'NA'
		2008-04-09
			no longer uses numpy matrix. just simple 2-d list.
		2008-4-2
		"""
		sys.stderr.write("Getting matrix ... " )
		#data_matrix = numpy.zeros([len(ecotype_id2index), len(phenotype_method_id2index)], numpy.float)
		data_matrix = [[]]*len(ecotype_id2index)
		for i in range(len(ecotype_id2index)):
			data_matrix[i] = ['NA']*len(phenotype_info.phenotype_method_id2index)
		#data_matrix[:] = numpy.nan
		query = "select pa.ecotype_id, pa.method_id, pa.value, pm.min_value, pm.stddev from %s pa, %s pm where pm.id=pa.method_id"%\
				(phenotype_avg_table, phenotype_method_table)
		if getPublicPhenotype:
			query += " and pm.access='public' and pa.ready_for_publication=1 "
		rows = curs.execute(query)
		is_elixirdb = 1
		if hasattr(curs, 'fetchall'):	#2009-2-2 this curs is not elixirdb.metadata.bind
			rows = curs.fetchall()
			is_elixirdb = 0
		
		for row in rows:
			if is_elixirdb:
				ecotype_id = row.ecotype_id
				phenotype_method_id = row.method_id
				value = row.value
				min_value = row.min_value
				stddev = row.stddev
			else:
				ecotype_id, phenotype_method_id, value, min_value, stddev = row
			if value==None:	#some db entries have nothing there. convert None to 'NA'
				value = 'NA'
			elif not get_raw_data:	#2008-11-10
				transformation_description = phenotype_info.phenotype_method_id2transformation_description.get(phenotype_method_id)
				#if value>-5e-7 and value<+5e-7:	#beyond float resolution by a python float
				#	value = 0
				
				if not transformation_description:
					pass
				elif transformation_description.find('Log(x)')!=-1:
					try:
						value = math.log10(value)
					except:
						sys.stderr.write("Ecotype ID %s, phenotype_method_id %s, value %s.\n"%(ecotype_id, phenotype_method_id, value))
						sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()[0]))
						traceback.print_exc()
						print sys.exc_info()
						#raise sys.exc_info()[0]
						sys.exit(2)
				elif transformation_description=="Log(SD/10+x-minVal)":	#2009-9-1 new transformation
					if min_value is not None and stddev is not None:
						value = math.log10(stddev/10. + value-min_value)
					else:
						value = value
				elif transformation_description=='Log(5+x)':
					value = math.log10(5+value)
				elif transformation_description=='Log(0.5+x)':
					value = math.log10(0.5+value)
				elif transformation_description=='(x-3)':
					value = value-3
			col_index = phenotype_info.phenotype_method_id2index[phenotype_method_id]
			data_matrix[ecotype_id2index[ecotype_id]][col_index] = value
		sys.stderr.write("Done\n")
		return data_matrix
	
	@classmethod
	def getPhenotypeData(cls, curs, phenotype_avg_table=None, phenotype_method_table=None, ecotype_table='stock.ecotype', get_raw_data=1,\
						getPublicPhenotype=False):
		"""
		2012.9.28
			add argument getPublicPhenotype
		2009-2-2
			wrap up all other 3 methods
		"""
		phenotype_info = cls.get_phenotype_method_id_info(curs, phenotype_avg_table=phenotype_avg_table, \
										phenotype_method_table=phenotype_method_table, getPublicPhenotype=getPublicPhenotype)
		ecotype_id2index, ecotype_id_ls, ecotype_name_ls = cls.get_ecotype_id2info(curs, phenotype_avg_table=phenotype_avg_table,\
																ecotype_table=ecotype_table, getPublicPhenotype=getPublicPhenotype)
		data_matrix = cls.get_matrix(curs, phenotype_avg_table, ecotype_id2index=ecotype_id2index, phenotype_info=phenotype_info, \
								get_raw_data=get_raw_data, phenotype_method_table=phenotype_method_table,\
								getPublicPhenotype=getPublicPhenotype)
		pheno_data = SNPData(col_id_ls=phenotype_info.phenotype_id_ls, row_id_ls=ecotype_id_ls, data_matrix=data_matrix)
		pheno_data.row_label_ls = ecotype_name_ls
		pheno_data.col_label_ls = phenotype_info.method_id_name_ls
		return pheno_data
	
	
	@classmethod
	def transform_log(cls,data,minValue=None,maxValue=None):
		import numpy
		phenotype_values = filter(lambda a: a != 'NA', data.values())
		constant = numpy.std(a = phenotype_values,ddof = 1)*0.1 - minValue
		for ecotype,phenotype_value in data.items():
			if data[ecotype] != "NA":
				data[ecotype] = math.log(phenotype_value + constant)
		return data
		
		
	
	def run(self):
		if self.debug==1:
			import pdb
			pdb.set_trace()
		
		#import MySQLdb
		#conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.db_user, passwd = self.db_passwd)
		#curs = conn.cursor()
		
		pheno_data = self.getPhenotypeData(self.db_250k.metadata.bind, self.phenotype_avg_table, self.phenotype_method_table, \
										self.ecotype_table, get_raw_data=self.get_raw_data,\
										getPublicPhenotype=self.getPublicPhenotype)
		header = ['ecotype id', 'nativename'] + pheno_data.col_label_ls
		write_data_matrix(pheno_data.data_matrix, self.outputFname, header, pheno_data.row_id_ls, pheno_data.row_label_ls, \
						transform_to_numpy=False)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = OutputPhenotype
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
	"""
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "user=", "passwd=", "help", "type=", "debug", "report"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:u:p:o:e:q:m:br", long_options_list)
	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
	hostname = None
	dbname = None
	user = None
	passwd = None
	outputFname = None
	ecotype_table = None
	phenotype_avg_table = None
	phenotype_method_table = None
	debug = None
	report = None
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
			sys.exit(2)
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-u", "--user"):
			user = arg
		elif opt in ("-p", "--passwd"):
			passwd = arg
		elif opt in ("-o",):
			outputFname = arg
		elif opt in ("-e",):
			ecotype_table = arg
		elif opt in ("-q",):
			phenotype_avg_table = arg
		elif opt in ("-m",):
			phenotype_method_table = arg
		
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	instance = OutputPhenotype(hostname=hostname, dbname=dbname, user=user, passwd=passwd, outputFname=outputFname,
					ecotype_table=ecotype_table, phenotype_avg_table=phenotype_avg_table, \
					phenotype_method_table = phenotype_method_table, debug=debug, report=report)
	instance.run()
	"""