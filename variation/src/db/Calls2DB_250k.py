#!/usr/bin/env python

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import sys, getopt, csv
import getopt, csv, subprocess
import traceback, gc
from pymodule import runLocalCommand
from pymodule import process_function_arguments
from pymodule import number2nt
from pymodule import ProcessOptions
from variation.src import Stock_250kDB
from variation.src.mapper.AbstractVariationMapper import AbstractVariationMapper

class Calls2DB_250k_old:
	"""
	2007-12-11
	Usage: Calls2DB_250k.py [OPTIONS] -i input_fname -n strain_name
	
	Option:
		-z ..., --hostname=...	the hostname, localhost(default)
		-d ..., --dbname=...	the database name, stock20071008(default)
		-k ..., --schema=...	which schema in the database, dbsnp(default)IGNORE
		-i ...,	input file, snp.RData outputted by write.table() in R
		-l ...,	calls_table, 'calls_250k'(default)
		-a ...,	calls_comment_table, 'calls_250k_duplicate_comment'(default)
		-s ...,	250k snp table, to find out which snpid based on chr+position, 'snps_250k'(default)
		-n ...,	strain name, it'll be used to match nativename and name to get ecotype id.
		-e ...,	comment for this loading
		-c	commit the database submission
		-b, --debug	enable debug
		-r, --report	enable more progress-related output
		-h, --help	show this help
	
	Examples:
		Calls2DB_250k.py -i /Network/Data/250k/yanli9-11-07/Tamm2B_base-calls.txt -n Tamm  -e "yanli9-11-07/Tamm2B_base-calls.txt" -c
	Description:
	"""
	def __init__(self, hostname='localhost', dbname='stock', schema='dbsnp', \
		input_fname='', calls_table='calls_250k', calls_comment_table='calls_250k_duplicate_comment', snps_table='snps_250k', strain_name=None, \
		comment='', commit=0, debug=0, report=0):
		"""
		2007-12-11
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname = input_fname
		self.calls_table = calls_table
		self.calls_comment_table = 'calls_250k_duplicate_comment'
		self.snps_table = snps_table
		self.strain_name = strain_name
		self.comment = comment
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def create_calls_table(self, curs, calls_table):
		"""
		2007-12-11
			won't be used. 2 tables, calls_250k_duplicate_comment and calls_250k will
			be created separately by mysql.
		"""
		sys.stderr.write("Creating %s ..."%calls_table)
		curs.execute("create table %s(\
			id	integer auto_increment primary key,\
			ecotypeid	integer,\
			snpid	integer,\
			snpcall	varchar(2),\
			duplicate	integer,\
			date_created	timestamp default current_timestamp,\
			date_modified	timestamp)"%calls_table)
		sys.stderr.write("Done.\n")
	
	def create_calls_comment_table(self, curs, calls_comment_table):
		"""
		2007-12-11
		"""
		sys.stderr.write("Creating %s ..."%calls_comment_table)
		curs.execute("create table %s(\
			id	integer auto_increment primary key,\
			ecotypeid	integer,\
			duplicate	integer,\
			comment	varchar(2000),\
			date_created	timestamp default current_timestamp,\
			date_modified	timestamp)"%calls_comment_table)
		sys.stderr.write("Done.\n")
	
	def get_chr_pos2snpid(self, curs, snps_table):
		"""
		2007-12-11
		2007-12-13
			'chr' in snps_table changed to 'chromosome'
		"""
		sys.stderr.write("Getting chr_pos2snpid from %s ..."%(snps_table))
		chr_pos2snpid = {}
		curs.execute("select id, chromosome, position from %s"%(snps_table))
		rows = curs.fetchall()
		for row in rows:
			snpid, chromosome, position = row
			chr_pos_key = (chromosome, position)
			chr_pos2snpid[chr_pos_key] = snpid
		sys.stderr.write("Done.\n")
		return chr_pos2snpid
	
	def find_out_ecotypeid_given_strain_name(self, curs, strain_name, ecotype_table='ecotype'):
		sys.stderr.write("Trying to find ecotypeid based on name=%s ...\n"%(strain_name))
		curs.execute("select e.id, e.name, e.stockparent, e.nativename, s.name, c.abbr from %s e, site s, address a, country c where e.siteid=s.id and s.addressid=a.id and a.countryid=c.id and (e.nativename rlike '%s'or e.name rlike '%s')"%(ecotype_table, strain_name, strain_name))
		rows = curs.fetchall()
		choice_id2ecotypeid = {}
		choice_id = 0
		header_ls = ['\t', 'ecotypeid', 'name', 'stockparent', 'nativename', 'site_name', 'country']
		sys.stderr.write("%s\n"%('\t'.join(header_ls)))
		for row in rows:
			ecotypeid, name, stockparent, nativename, site_name, country = row
			choice_id += 1
			choice_id2ecotypeid[repr(choice_id)] = ecotypeid
			ls = ['%s:'%choice_id]
			ls.append(list(row))
			ls = map(repr, ls)
			sys.stderr.write("%s\n"%'\t'.join(ls))
		if len(choice_id2ecotypeid)==0:
			sys.stderr.write("Found no match for %s.\n"%strain_name)
			sys.exit(3)
		choice_id = raw_input("Enter the row number corresponding to the ecotypeid(q for exit):")
		while 1:
			if choice_id =='q':
				sys.exit(2)
			elif choice_id in choice_id2ecotypeid:
				break
			else:
				sys.stderr.write("%s is not a choice.\n"%(choice_id))
				choice_id = raw_input("Enter the row number corresponding to the ecotypeid(q for exit):")
		return choice_id2ecotypeid[choice_id]
	
	def find_out_duplicate_id_given_ecotypeid(self, curs, ecotypeid, calls_table):
		"""
		2007-12-11
		"""
		sys.stderr.write("Finding out the duplicate_id for ecotypeid=%s ..."%(ecotypeid))
		curs.execute("select distinct duplicate from %s where ecotypeid=%s"%(calls_table, ecotypeid))
		rows = curs.fetchall()
		duplicate_ls = []
		for row in rows:
			duplicate_ls.append(row[0])
		if duplicate_ls:
			duplicate_id = max(duplicate_ls) + 1	#1 plus the maximum
		else:
			duplicate_id = 1
		sys.stderr.write("is %s.\n"%duplicate_id)
		return duplicate_id
	
	def get_calls(self, input_fname, chr_pos2snpid):
		"""
		2007-12-11
		"""
		sys.stderr.write("Getting calls ...")
		reader = csv.reader(open(input_fname), delimiter='\t')
		calls_ls = []
		reader.next()	#toss out the 1st header row
		for row in reader:
			chr, position, allele1, allele2, antisense1, sense1, antisense2, sense2, snpcall = row
			chr = int(chr)
			position = int(position)
			chr_pos_key = (chr, position)
			if snpcall=='?':	#2007-12-11 '?' was used as NA
				snpcall = 'N'
			snpid = chr_pos2snpid[chr_pos_key]
			calls_ls.append([snpid, snpcall])
		del reader
		sys.stderr.write("Done.\n")
		return calls_ls
	
	def submit_comment(self, curs, ecotypeid, duplicate_id, comment, calls_comment_table):
		"""
		2007-12-11
		"""
		sys.stderr.write("Submitting comment ...")
		curs.execute("insert into %s(ecotypeid, duplicate, comment) values(%s, %s, '%s')"%(calls_comment_table, ecotypeid, duplicate_id, comment))
		sys.stderr.write("Done.\n")
	
	def submit_calls_ls(self, curs, calls_ls, ecotypeid, duplicate_id, calls_table):
		sys.stderr.write("Submitting calls_ls ...")
		for snpid, snpcall in calls_ls:
			curs.execute("insert into %s(ecotypeid, snpid, snpcall, duplicate) values(%s, %s, '%s', %s)"%(calls_table, ecotypeid, snpid, snpcall, duplicate_id))
		sys.stderr.write("Done.\n")
	
	def run(self):
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname,host=self.hostname)
		curs = conn.cursor()
		ecotypeid = self.find_out_ecotypeid_given_strain_name(curs, self.strain_name)
		duplicate_id = self.find_out_duplicate_id_given_ecotypeid(curs, ecotypeid, self.calls_table)
		chr_pos2snpid = self.get_chr_pos2snpid(curs, self.snps_table)
		calls_ls = self.get_calls(self.input_fname, chr_pos2snpid)
		if self.commit:
			self.submit_comment(curs, ecotypeid, duplicate_id, self.comment, self.calls_comment_table)
			self.submit_calls_ls(curs, calls_ls, ecotypeid, duplicate_id, self.calls_table)

class Calls2DB_250k(AbstractVariationMapper):
	"""

	Examples:
		Calls2DB_250k.py -i /tmp/simplecalls -m 1 --commit
		
		Calls2DB_250k.py -i /Network/Data/250k/finalData_051808/250K_method_5_after_imputation_noRedundant_051908.csv
			-m 6 -u yh -y 2 --commit -l 262_unique
		
		#2011-5-2 add a CNV dataset (column id is Snps.id, 0=normal, 1=deletion) into db
		# The program will use pymodule.SNP.number2nt to translate 0 & 1 into nucleotides (NA and A).
		Calls2DB_250k.py -i script/variation/data/CNV/NonOverlapCNVAsSNP_cnvMethod20_snpID.tsv -m 80 -u yh -y3
			-l cnvMethod20InCall32 --hostname uclaOffice --commit
		
	Description:
		Turn calling algorithm's results into db and associated filesystem directory.
		
		input_type 1:
			Each file in input_dir shall be named like 'array_id'_call.tsv.
			The file would be ignored if a call with same array_id and same method_id exists in database.
			The format is 2-column and tab-delimited. SNP_ID could be chr_pos or Snps.id. example:
				SNP_ID	'array_id'
				1_657_C_T	C
				3	T
		2011-5-3 Input type 2 is outdated.
		output: Apart from entries inserted into db table call_info, files would be created in the output_dir to store
			actual calls.
			Format is 2-column and tab-delimited. First column is id in table snps. example:
				SNP_DB_ID 'array_id'
				213444	C
	"""
	option_default_dict = AbstractVariationMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
							('input_dir', 1, ): [None, 'i', 1, "directory containing output files of any calling algorithm. it could aslso be a file, assuming it's in bjarni's format with arrayId."],\
							('method_id', 1, int): [None, 'm', 1, 'the id of the calling method. If not in db, new call method is to be created with id=method_id, short_name=(user supplied).'],\
							('call_method_short_name', 0, ): ['', 'l', 1, 'call method shortname, given to create a new call method if method_id is not in db.'],\
							('output_dir', 1, ):['/Network/Data/250k/db/calls/', 'o', 1, 'file system storage for the call files. call_info_table would point each entry to this.'],\
							('call_method_table', 1, ): ['call_method', 'a', 1, 'table storing the calling methods'],\
							('call_info_table', 1, ): ['call_info', 't', 1, 'table to store final call file entries'],\
							('input_type', 1, int): [1, 'y', 1, 'The input type. 1: directory. \n\
		2: SNP X strain format (bjarni, 1st row is array id. 2nd row is ecotype id, not used though.). \n\
		3: Strain X SNP format (Yu)'],\
						})
	"""
	2008-04-40
		option_default_dict is a dictionary for option handling, including argument_default_dict info
		the key is a tuple, ('short_option', 'long_option', has_argument, description_for_option, is_option_required, argument_type)
		argument_type is optional
	"""
	def __init__(self, **keywords):
		"""
		2008-04-08
		"""
		AbstractVariationMapper.__init__(self, inputFnameLs=None, **keywords)
		self.cur_max_call_id = None
		
		call2db_func_dict = {1: self.submit_call_dir2db,\
							2: self.submit_SNPxStrain_file2db,\
							3: self.submit_StrainxSNP_file2db}
		self.submit_call2db = call2db_func_dict.get(self.input_type)
		if not self.submit_call2db:
			sys.stderr.write("Error: Input type %s is not available.\n"%self.input_type)
	
	def get_cur_max_call_id(self, curs, call_info_table):
		"""
		get current maximum call id in db
		"""
		sys.stderr.write("Getting current maximum call id in db.\n")
		curs.execute("select max(id) from %s"%(call_info_table))
		rows = curs.fetchall()
		cur_max_call_id = rows[0][0]
		if cur_max_call_id!=None:
			return cur_max_call_id
		else:
			return 0
	
	def check_method_id_exists(self, db, method_id):
		"""
		2011-5-3
			directly query Stock_250kDB.CallMethod.
		"""
		cm = Stock_250kDB.CallMethod.get(method_id)
		if cm is None:
			return 0
		else:
			return 1
		"""
		curs.execute("select id from %s where id=%s"%(call_method_table, method_id))
		rows = curs.fetchall()
		if len(rows)==0:
			return 0
		else:
			return 1
		"""
	
	def get_new_call_id(self, curs, call_info_table, array_id, method_id):
		"""
		"""
		if self.cur_max_call_id==None:
			self.cur_max_call_id = self.get_cur_max_call_id(curs, call_info_table)
		curs.execute("select id from %s where array_id=%s and method_id=%s"%(call_info_table, array_id, method_id))
		rows = curs.fetchall()
		if len(rows)>0:
			sys.stderr.write("\tarray_id=%s and method_id=%s already exists in %s. Ignored.\n"%(array_id, method_id, call_info_table))
			return -1
		else:
			self.cur_max_call_id += 1
			return self.cur_max_call_id
	
	def submit_one_call_entry(self, curs, call_info_table, call_id, filename, array_id, method_id):
		"""
		2008-04-08
			not used right now.
		"""
		pass
	
	def submit_call_dir2db(self, curs, input_dir, call_info_table, output_dir, method_id, user, chr_pos2db_id=None, db=None):
		"""
		2011-5-2
			curs is completely useless. keywords must contain the "db" key now.
		2011-2-27
			input file could use either db_id or chr_pos to identify locus.
		2011-1-24
			bugfix
				Call files to be added into method 3 have the 3rd column, oligo_probability.
		2010-10-13
			add argument chr_pos2db_id and db
			it replaces the old snp ID (chr_pos_...) with id from table Stock_250kDB.Snps
		2008-04-11
			check if output_fname exists already or not. if yes, ignore.
			use subprocess.Popen to do cp
			add method_ in front of the method_id sub-directory
		2008-04-09
			add method_id as sub-directory
			submit user into table as 'created_by'
		2008-04-08
		"""
		output_dir = os.path.join(output_dir, 'method_%s'%method_id)
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		file_ls = os.listdir(input_dir)
		sys.stderr.write("\n\tTotally, %d files to be processed.\n"%len(file_ls))
		file_ls.sort()
		for i in range(len(file_ls)):
			filename = file_ls[i]
			array_id = filename.split('_')[0]
			sys.stderr.write("%d/%d:\t%s\n"%(i+1,len(file_ls),filename))
			
			call_info = db.getCallInfo(array_id=array_id, call_method_id=method_id, created_by=user)
			output_fname = os.path.join(output_dir, '%s_call.tsv'%call_info.id)
			call_info.filename = output_fname
			db.session.add(call_info)
			db.session.flush()
			
			#new_call_id = self.get_new_call_id(curs, call_info_table, array_id, method_id)
			#if new_call_id!=-1:
				#output_fname = os.path.join(output_dir, '%s_call.tsv'%new_call_id)
			if os.path.isfile(output_fname):
				sys.stderr.write("%s already exists. Ignore.\n"%output_fname)
				continue
			input_fname = os.path.join(input_dir, filename)
			reader = csv.reader(open(input_fname), delimiter='\t')
			writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
			header_row = reader.next()
			writer.writerow(header_row)
			for row in reader:
				db_id = db.get_db_id_given_chr_pos2db_id(row[0],)
				if db_id is not None:
					new_row = [db_id] + row[1:]
					writer.writerow(new_row)
			del reader, writer
			"""
			call_info = Stock_250kDB.CallInfo(id=new_call_id, filename=output_fname, array_id=array_id,\
											method_id=method_id, created_by=user)
			db.session.add(call_info)
			db.session.flush()
			
			
			cp_p = subprocess.Popen(['cp', input_fname, output_fname], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
			cp_p_stdout_out = cp_p.stdout.read()
			cp_p_stderr_out = cp_p.stderr.read()
			if cp_p_stdout_out:
				sys.stderr.write("\tcp stdout: %s\n"%cp_p_stdout_out)
			if cp_p_stderr_out:
				sys.stderr.write("\tcp stderr: %s\n"%cp_p_stderr_out)
				continue	#error in cp. skip the db insertion.
			curs.execute("insert into %s(id, filename, array_id, method_id, created_by) values (%s, '%s', %s, %s, '%s')"%\
					(call_info_table, new_call_id, output_fname, array_id, method_id, user))
			"""
	
	def submit_SNPxStrain_file2db(self, curs, input_fname, call_info_table, output_dir, method_id, user, chr_pos2db_id=None,\
								**keywords):
		"""
		2010-10-13
			add argument chr_pos2db_id and **keywords
			it replaces the old snp ID (chr_pos_...) with id from table Stock_250kDB.Snps
		2008-05-17
			submit the calls from a matrix file to db
		"""
		sys.stderr.write("Submitting %s to db ...\n"%(input_fname))
		output_dir = os.path.join(output_dir, 'method_%s'%method_id)
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		
		reader = csv.reader(open(input_fname))
		array_id_ls = reader.next()
		column_index2writer = {}
		for i in range(2, len(array_id_ls)):
			array_id = array_id_ls[i]
			sys.stderr.write("%s\tAssign new call info id to array id=%s ."%('\x08'*80, array_id))
			new_call_id = self.get_new_call_id(curs, call_info_table, array_id, method_id)
			if new_call_id!=-1:
				output_fname = os.path.join(output_dir, '%s_call.tsv'%new_call_id)
				if os.path.isfile(output_fname):
					sys.stderr.write("%s already exists. Ignore.\n"%output_fname)
					continue
				writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
				writer.writerow(['SNP_ID', array_id])
				column_index2writer[i] = writer
				curs.execute("insert into %s(id, filename, array_id, method_id, created_by) values (%s, '%s', %s, %s, '%s')"%\
						(call_info_table, new_call_id, output_fname, array_id, method_id, user))
		reader.next()	#ignore the ecotype id line
		
		sys.stderr.write("Moving real data to file system storage ...\n")
		counter = 0
		for row in reader:
			chromosome = int(row[0])
			position = int(row[1])
			counter += 1
			chr_pos = (chromosome, position)
			for i in range(2, len(row)):
				if i in column_index2writer:
					if chr_pos2db_id:	#2010-10-13
						db_id = chr_pos2db_id.get(chr_pos)
					else:
						db_id = chr_pos
					if db_id is not None:
						column_index2writer[i].writerow([db_id, row[i]])
			if counter%5000==0:
				sys.stderr.write("%s\t%s"%('\x08'*20, counter))
		sys.stderr.write("%s\t%s"%('\x08'*20, counter))
		del reader
		for column_index, writer in column_index2writer.iteritems():
			del writer
		sys.stderr.write(" %s arrays. Done.\n"%(len(column_index2writer)))
	
	def submit_StrainxSNP_file2db(self, curs, input_fname, call_info_table, output_dir, method_id, user, chr_pos2db_id=None,
								**keywords):
		"""
		2011-5-2
			curs is completely useless. keywords must contain the "db" key now.
		2011-2-27
			input file could use either db_id or chr_pos to identify locus.
		2010-10-13
			add argument chr_pos2db_id and **keywords
			it replaces the old snp ID (chr_pos_...) with id from table Stock_250kDB.Snps
		2008-1-5
			if the output_fname already exists, exit the program.
			if db insertion fails, delete the file written out and exit the program.
		2008-05-19
			submit the calls from a matrix file (Strain X SNP format, tsv, nucleotides in numbers) to db
		"""
		sys.stderr.write("Submitting %s to db ...\n"%(input_fname))
		db = keywords.get("db")	#2011-5-2
		output_dir = os.path.join(output_dir, 'method_%s'%method_id)
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		
		reader = csv.reader(open(input_fname), delimiter='\t')
		header = reader.next()
		counter = 0
		for row in reader:
			ecotype_id, array_id = row[:2]
			sys.stderr.write("%s\tAssign new call info id to array id=%s ."%('\x08'*80, array_id))
			
			call_info = Stock_250kDB.CallInfo(array_id=array_id, method_id=method_id, created_by=user)
			db.session.add(call_info)
			db.session.flush()
			output_fname = os.path.join(output_dir, '%s_call.tsv'%call_info.id)
			call_info.filename = output_fname
			db.session.add(call_info)
			db.session.flush()
			
			#new_call_id = self.get_new_call_id(curs, call_info_table, array_id, method_id)
			#if new_call_id!=-1:
			
			if os.path.isfile(output_fname):
				sys.stderr.write("Error: %s already exists. Check why the file exists while db has no record.\n"%output_fname)
				sys.exit(2)
			writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
			writer.writerow(['SNP_ID', array_id])
			for i in range(2, len(row)):
				snp_id = header[i]
				if chr_pos2db_id:	#2010-10-13
					db_id = db.get_db_id_given_chr_pos2db_id(snp_id)
				else:
					db_id = snp_id
				if db_id is not None:
					writer.writerow([db_id, number2nt[int(row[i])]])	#translate 
			del writer
			"""
			try:
				
				#curs.execute("insert into %s(id, filename, array_id, method_id, created_by) values (%s, '%s', %s, %s, '%s')"%\
				#			(call_info_table, new_call_id, output_fname, array_id, method_id, user))
			except:
				traceback.print_exc()
				sys.stderr.write('%s.\n'%repr(sys.exc_info()))
				sys.stderr.write('Error encountered while inserting record into %s. Delete the file written.\n'%call_info_table)
				commandline = 'rm %s'%(output_fname)
				command_out = runLocalCommand(commandline)
				sys.exit(3)
			"""
			counter += 1
		del reader
		sys.stderr.write(" %s arrays. Done.\n"%counter)
	
	def run(self):
		"""
		2008-05-17
			-check_method_id_exists()
			if input_dir is dir:
				-submit_call_dir2db()
					-get_new_call_id()
						-get_cur_max_call_id()
			elif input_dir is file:
				-submit_call_file2db()
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		#database connection and etc
		db = self.db_250k
		session = db.session
		session.begin()
		
		chr_pos2db_id = db.getSNPChrPos2ID()
		
		#import MySQLdb
		#conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.user, passwd = self.passwd)
		#curs = conn.cursor()
		curs = None
		if not self.check_method_id_exists(db, self.method_id):
			sys.stderr.write("Warning: method_id=%s not in %s. A new entry to be created.\n"%\
							(self.method_id, self.call_method_table))
			cm = Stock_250kDB.CallMethod(short_name=self.call_method_short_name, id=self.method_id)
			session.add(cm)
			session.flush()
			
			self.method_id = cm.id
			
		if self.commit:
			self.submit_call2db(curs, self.input_dir, self.call_info_table, self.output_dir, self.method_id, self.db_user, \
							chr_pos2db_id=chr_pos2db_id, db=db)
			#curs.execute("commit")
			session.flush()
			session.commit()
	
if __name__ == '__main__':
	main_class = Calls2DB_250k
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()