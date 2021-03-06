#!/usr/bin/env python
"""

Examples:
	#test run without commiting database (no records in database in the end)
	%s -i /home/nordborglab/pvalue.log -a 1 -E 1 -f kw_test_96_LD -m kw -n best_96_by_tina -u yh
	
	#commit transaction
	%s -i /home/nordborglab/pvalue.log -a 1 -E 1 -f kw_test_96_LD -m kw -n best_96_by_tina -u yh -c
	
	#omit short_name
	%s -a 17 -E 186 -i /Network/KW_newDataset_186_Bact_titer.pvals -l 1 -u yh -c
	
	# 2011-2-24 submit cnv association results (locus identified by CNV.ID) into db.
	# set the call_method_id (-a) to -1 so that no call method will be fetched from db.
	%s -i cnvMethod20_cnvID/cnvMethod20_cnvID_y1_pheno_pheno_183_Trichome\ avg\ C.tsv 
		-g 20 -a -1 -E 183 -l 1 -u yh -z banyan -s 3 -c
	
Description:
	This program would submit simple association results into database.

	The input file format could be 3-column, tab-delimited or 4-column with MAF (minor allele frequency) or more.
		snp.id	2nd-column-to-ignore	score/pvalue	MAF
		
	If analysis_method_id is 13 (boolean SNP pair), format is totally different. This program just copies it over.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

import csv, stat, getopt, re
import traceback, gc, subprocess
import sqlalchemy as sql
from pymodule import figureOutDelimiter, PassingData, utils
from pymodule import AssociationTableFile, addAttributeDictToYHTableInHDF5Group
from variation.src.common import getOneResultJsonData
from variation.src import Stock_250kDB
from variation.src import AbstractVariationMapper
"""
2008-04-16 temporarily put here
	-s ...,	short_name*	give a short name of what you did. try to incorporate method, genotype data and phenotype data
	-m ...,	method_description*	a longer description of your method
	-a ...,	data_description*	which data you used
"""

class Results2DB_250k(AbstractVariationMapper):
	__doc__ = __doc__	#use documentation in the beginning of the file as this class's doc
	"""
	('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['stock_250k', 'd', 1, '', ],\
						('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
	"""
	option_default_dict = AbstractVariationMapper.option_default_dict.copy()
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('short_name', 0, ): [None, 'f', 1, 'short name for this result. Must be unique from previous ones. \
							combining phenotype, data, method is a good one. If not given, will be automatically generated.' ],\
						('phenotype_method_id',1,int): [None, 'E', 1, 'which phenotype you used, check table phenotype_method'],\
						('call_method_id', 1, int ): [None, 'a', 1, 'data from which call_method, field id in table call_method'],\
						('data_description', 0, ): [None, 'n', 1, 'Describe how your data is derived from that call method. like non-redundant set, 1st 96, etc.'],\
						('method_description', 0, ): [None, 'm', 1, 'Describe your method and what type of score, association (-log or not), recombination etc.'],\
						('results_method_type_id', 1, int): [1, 's', 1, 'which type of method. field id in table results_method_type. 1="association"',],\
						('analysis_method_id', 1, int): [None, 'l', 1, ''],\
						('cnv_method_id', 0, int): [None, 'g', 1, 'for CNV association results, need this cnv_method_id'],\
						('comment',0, ): [None, 't', 1, 'Anything more worth for other people to know?'],\
						})
	
	pa_has_characters = re.compile(r'[a-zA-Z_]')	#2008-05-30 check if a string has character in it, used to judge whether the 1st line is header or not.
	"""
	04/28/08 no longer needed
							('results_table',1, ): 'results',\
							('results_method_table',1, ):'results_method',\
							('phenotype_method_table',1, ):'phenotype_method',\	
	"""
	def __init__(self, **keywords):
		"""
		2008-04-28
			use ProcessOptions, newer option handling class
		2008-04-16
		"""
		AbstractVariationMapper.__init__(self, inputFnameLs=None, **keywords)
		
	
	def check_if_phenotype_method_id_in_db(self, curs, phenotype_method_table, phenotype_method_id):
		"""
		"""
		curs.execute("select id from %s where id=%s"%(phenotype_method_table, phenotype_method_id))
		rows = curs.fetchall()
		if len(rows)>0:
			return 1
		else:
			return 0
	
	def submit_results_method(self, curs, results_method_table, short_name, method_description, data_description):
		"""
		2008-04-16
			submit the method part into db and return the id
		"""
		sys.stderr.write("Submitting results method ...")
		curs.execute("insert into %s(short_name, method_description, data_description) values ('%s', '%s', '%s')"%\
					(results_method_table, short_name, method_description, data_description))
		curs.execute("select id from %s where short_name='%s'"%\
					(results_method_table, short_name))
		rows = curs.fetchall()
		sys.stderr.write("Done.\n")
		return rows[0][0]
		
	marker_pos2snp_id = None
	is_new_marker_added = False	#2008-05-26 flag whether new markers were generated. to check after commit/rollback
	
	def reset_marker_pos2snp_id(self):
		"""
		2008-05-27
			set "self.is_new_marker_added = False"
		2008-05-26
			after commit or rollback in plone, session is closed and those new marker objects are gone. need to reset everything.
		"""
		if self.is_new_marker_added:
			del self.marker_pos2snp_id
			self.marker_pos2snp_id = self.get_marker_pos2snp_id(db)
			self.is_new_marker_added = False
	
	def get_marker_pos2snp_id(self, db):
		"""
		2008-05-24
		"""
		sys.stderr.write("Getting marker_pos2snp_id ...")
		marker_pos2snp_id = {}
		snps_table = db.tables['snps'].alias()
		conn = db.connection
		results = conn.execute(sql.select([snps_table.c.id, snps_table.c.chromosome, snps_table.c.position, snps_table.c.end_position]))
		for row in results:
			key = (row.chromosome, row.position, row.end_position)
			marker_pos2snp_id[key] = row.id
		sys.stderr.write("Done.\n")
		return marker_pos2snp_id
	
	def copyAndReformatResultFile(self, db, inputFname=None, db_entry=None, user=None, output_fname=None):
		"""
		2011-2-22
			Locus are now identified as Snps.id / CNV.id in association result files. (chr, pos) before.
		2009-1-7
			insert float into the middle below
				column_5th=int(float(row[4]))	#int('89.0') would raise an exception
		2008-11-12
			parse lines with column_6(genotype_var_perc) and more (comment)
		2008-09-30
			deal with 5-column file. The 5-th column is minor allele count.
			also return True in the end. return False if error in the middle.
		2008-08-19
			add original_filename to ResultsMethod
		2008-07-16
			if inputFname is neither file name nor file object, exit the program
			better handling of the column_4th and its header
		2008-07-16
			if it's 4-column, the last one is MAF.
			can't deal with segment score anymore.
		2008-05-30
			merged with store_file()
				dump the file onto file system storage if output_fname is given
				db submission is too slow
		2008-05-26
			inputFname from plone is not file object although it has file object interface.
		2008-05-26
			csv.Sniffer() can't figure out delimiter if '\n' is in the string, use own dumb function figureOutDelimiter()
		2008-05-25
			save marker(snps) in database if it's not there.
			use marker id in results table
		2008-05-24
			figure out delimiter automatically
			inputFname could be a file object (from plone)
			phenotype method doesn't go with results anymore. it goes with results_method
		2008-04-28
			changed to use Stock_250kDatabase (SQLAlchemy) to do db submission
		"""
		if (isinstance(inputFname, str) or isinstance(inputFname, unicode)) and os.path.isfile(inputFname):
			sys.stderr.write("Submitting results from %s ..."%(os.path.basename(inputFname)))
			delimiter = figureOutDelimiter(inputFname)
			reader = csv.reader(open(inputFname), delimiter=delimiter)
			db_entry.original_filename = inputFname
		elif hasattr(inputFname, 'readline') or hasattr(inputFname, 'read'):	#inputFname is not a file name, but direct file object. it could also be <ZPublisher.HTTPRequest.FileUpload instance at 0xa1774f4c>
			sys.stderr.write("Submitting results from %s on plone ..."%inputFname.filename)
			cs = csv.Sniffer()
			inputFname.seek(0)	#it's already read by plone to put int data['inputFname'], check results2db_250k.py
			if getattr(inputFname, 'readline', None) is not None:
				test_line = inputFname.readline()
				delimiter = cs.sniff(test_line).delimiter
			else:
				test_line = inputFname.read(200)
				delimiter = figureOutDelimiter(test_line)	#counting is a safer solution. if test_line include '\n', cs.sniff() won't figure it out.
			inputFname.seek(0)
			reader = csv.reader(inputFname, delimiter=delimiter)
			if getattr(inputFname, 'filename', None):
				db_entry.original_filename = getattr(inputFname, 'filename', None)
			else:
				db_entry.original_filename = getattr(inputFname, 'name', None)
		else:
			sys.stderr.write("Error: %s is neither a file name nor a file object.\n"%inputFname)
			sys.exit(4)
		
		if output_fname:
			if os.path.isfile(output_fname):
				sys.stderr.write("Error: file %s already exists. Skip.\n"%output_fname)
				return False
			writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		elif self.marker_pos2snp_id is None:
			self.marker_pos2snp_id = self.get_marker_pos2snp_id(db)
		
		header_outputted = 0
		no_of_lines = 0
		
		session = db.session
		for row in reader:
			#check if 1st line is header or not
			if no_of_lines ==0 and self.pa_has_characters.search(row[1]):	#check the 2nd one, which is strict digits. while the 1st column, chromosome could be 'X' or something
				continue
			snp_id = int(row[0])
			if row[1] and row[1]!='0':	#2011-2-22 something on 2nd column. wrong format.
				chr = int(row[0])
				start_pos = int(row[1])
				sys.stderr.write("Error: current version doesn't take chr,pos as marker ID anymore. Has to be one id (either Snps.id or CNV.id).\n")
				sys.exit(4)
			score = row[2]
			stop_pos = None
			column_4th = None
			column_5th = None
			column_6 = None
			rest_of_row = []
			rest_of_header = []
			
			#marker_name = '%s_%s'%(chr, start_pos)	#2011-2-22
			if len(row)>=4:
				column_4th=row[3]
				#stop_pos = int(row[2])
				#score = row[3]
			if len(row)>=5:
				#column_4th=row[3]
				column_5th=int(float(row[4]))	#2009-1-7 int('89.0') would raise an exception
			if len(row)>=6:
				column_6 = row[5]
			if len(row)>=7:
				rest_of_row = row[6:]
				rest_of_header = ['beta%s'%i for i in range(len(rest_of_row))]
				#sys.stderr.write("ERROR: Found %s columns.\n"%(len(row)))
				#return False
			
			if output_fname:	#go to file system
				if not header_outputted:	#3-column or 4-column header
					#if stop_pos is not None:	#2011-2-22
					#	position_header = ['start_position', 'stop_position']
					#else:
					#	position_header = ['position']
					header = ['snp_id', 'none', 'score']	#2011-2-22
					if column_4th is not None:
						header.append('MAF')
					if column_5th is not None:
						header.append('MAC')	#Minor Allele Count
					if column_6 is not None:
						header.append('genotype_var_perc')	#genotype variance percentage
					if rest_of_row:
						header += rest_of_header
					writer.writerow(header)
					header_outputted = 1
				#data_row = [chr, start_pos]	#2011-2-22
				data_row = [snp_id, '']	#2011-2-22
				#if stop_pos is not None:	#2011-2-22
				#	data_row.append(stop_pos)
				data_row.append(score)
				if column_4th is not None:
					data_row.append(column_4th)
				if column_5th is not None:
					data_row.append(column_5th)
					if db_entry.no_of_accessions is None:	#calculate the no_of_accessions based on MAC/MAF
						db_entry.no_of_accessions = int(round(float(column_5th)/float(column_4th)))
				if column_6 is not None:
					data_row.append(column_6)
				if rest_of_row:
					data_row += rest_of_row
				writer.writerow(data_row)
			"""
			# 2011-2-22 store the results directly into db. only for old SNP association results.
			else:
				
				key = (chr, start_pos, stop_pos)
				if key in self.marker_pos2snp_id:
					snps_id = self.marker_pos2snp_id[key]
					if isinstance(snps_id, SNPs):	#it's a new marker object
						r = Results(score=score)
						r.snps = snps_id
					else:	#others are all integer ids
						r = Results(snps_id=snps_id, score=score)
				else:
					#construct a new marker
					marker = SNPs(name=marker_name, chromosome=chr, position=start_pos, end_position=stop_pos, created_by=user)
					#save it in database to get id
					session.add(marker)
					self.marker_pos2snp_id[key] = marker	#for the next time to encounter same marker
					self.is_new_marker_added = True	#set this flag as new marker was inputted into the dict
					r = Results(score=score)
					r.snps = marker
					del marker
				r.results_method = db_entry
				session.add(r)
				del r
			"""
			no_of_lines += 1
		
		del reader
		if output_fname:
			del writer
		sys.stderr.write("Done.\n")
		return 0
	
	
	def come_up_new_results_filename(self, output_dir=None, db_entry_id=None, results_method_type_id=None,
							call_method_id=None, phenotype_method_id=None,\
							analysis_method_id=None, cnv_method_id=None, extraBitLs=None, suffix='result'):
		"""
		2012.11.13 deprecated, use ResultsMethod.constructRelativePath() instead
		2012.11.13 adjust to handle result_landscape filename as well
			added extraBitLs
		2011-2-22
			add cnv_method_id
		2010-9-21
			add argument call_method_id, phenotype_method_id, analysis_method_id
				and incorporate the non-None ones into the filename 
		2008-05-30
			to save the filename into db
		"""
		if results_method_type_id:
			output_dir = os.path.join(output_dir, 'type_%s'%results_method_type_id)
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		
		filename_part_ls = [db_entry_id]
		if call_method_id is not None:
			filename_part_ls.append(call_method_id)
		if cnv_method_id is not None:
			filename_part_ls.append(cnv_method_id)
		if phenotype_method_id is not None:
			filename_part_ls.append(phenotype_method_id)
		if analysis_method_id is not None:
			filename_part_ls.append(analysis_method_id)
		if extraBitLs:
			filename_part_ls.extend(extraBitLs)
		filename_part_ls = map(str, filename_part_ls)
		if suffix:
			filename_part_ls.append(suffix)
		output_fname = os.path.join(output_dir, '%s.tsv'%('_'.join(filename_part_ls)))
		return output_fname
	
	def copyResultsFile(self, db, inputFname, db_entry, user, output_fname=None):
		"""
		2008-09-30
			return True
		2008-09-09
			similar task to copyAndReformatResultFile, but not look into the file, just copy the file
		"""
		sys.stderr.write("Copying results from %s ..."%(os.path.basename(inputFname)))
		db_entry.original_filename = inputFname
		pipe_f = os.popen('cp %s %s'%(inputFname, output_fname))
		pipe_f_out = pipe_f.read()
		if pipe_f_out:
			sys.stderr.write("\t cp output: %s\n"%pipe_f_out)
		sys.stderr.write("Done.\n")
		return 0
		
	
	def add2DB(self, db=None, short_name=None, phenotype_method_id=None, call_method_id=None, data_description=None, \
				method_description=None, comment=None, inputFname=None, user=None, results_method_type_id=None, \
				analysis_method_id=None, results_method_type_short_name=None, data_dir=None, commit=0,\
				cnv_method_id=None):
		"""
		2012.12.28 overhaul
		2012.6.6
			pass db to getOneResultJsonData()
		2012.3.9
			add locus_type_id to ResultsMethod
		2011-2-22
			add argument cnv_method_id
			deal with the association file format change. locus is now identified by Snps.id or CNV.id
		2010-5-3
			becomes classmethod
			store the json structure of top 10000 SNPs from the db_entry into db
		2008-09-30
			don't save results_method into database if bad thing happend when getting data out of the file.
		2008-09-09
			directly copy the result file if analysis_method_id==13
		2008-08-19
			automatically generate short_name if it's NULL
		2008-07-16
			adjust to new Elixir-based db api.
			new analysis_method_id is added to results_method.
		2008-05-30
			go to output_dir
			drop copyAndReformatResultFile()
			use store_file()
		2008-05-26
			add results_method_type_id and results_method_type_short_name
		2008-05-24
			to conveniently wrap up all codes so that both this program and plone can call
		"""
		session = db.session
		session.begin()
		
		
		rmt = Stock_250kDB.ResultsMethodType.get(results_method_type_id)
		if not rmt and results_method_type_short_name is not None:	#create a new results method type
			rmt = Stock_250kDB.ResultsMethodType(short_name=results_method_type_short_name)
			session.add(rmt)
		
		if not rmt:
			sys.stderr.write("No results method type available for results_method_type_id=%s.\n"%results_method_type_id)
			sys.exit(3)
		
		if call_method_id:	#2012.6.6
			cm = Stock_250kDB.CallMethod.query.get(call_method_id)
			locus_type_id = cm.locus_type_id
		else:
			cm = None
			locus_type_id = None
		
		db_entry = db.checkResultsMethod(call_method_id=call_method_id, phenotype_method_id=phenotype_method_id, \
										analysis_method_id=analysis_method_id, \
			cnv_method_id=cnv_method_id, accession_set_id=None, results_method_type_id=results_method_type_id)
		if db_entry:
			sys.stderr.write("There is already an entry in results_method (id=%s) with same (call_method_id, phenotype_method_id, analysis_method_id, results_method_type_id)=(%s, %s, %s, %s).\n"\
							%(db_entry.id, call_method_id, phenotype_method_id, analysis_method_id, results_method_type_id))
			sys.exit(2)
		
		db_entry = db.getResultsMethod(data_dir=data_dir, call_method_id=call_method_id, phenotype_method_id=phenotype_method_id, \
									analysis_method_id=analysis_method_id, \
			cnv_method_id=cnv_method_id, accession_set_id=None, results_method_type_id=results_method_type_id,\
			method_description=method_description, no_of_accessions=None, \
			no_of_loci=None, filename=None, original_filename=inputFname, \
			data_description=data_description, comment=comment, created_by=user, locus_type_id=locus_type_id)	#2012.3.9
		
		if commit:
			db_entry.filename = os.path.join(db.data_dir, db_entry.constructRelativePath(data_dir=data_dir))
			localAbsPath = os.path.join(data_dir, db_entry.constructRelativePath(data_dir=data_dir))
			
			if db_entry.analysis_method_id==13:
				self.srcFilenameLs.append(inputFname)
				self.dstFilenameLs.append(localAbsPath)
				exit_code = self.copyResultsFile(db, inputFname, db_entry, user=user, output_fname=localAbsPath)
			else:
				#2013.1.10 add some db_entry attributes to the hdf5 file
				db.addAttributesToResultFile(db_entry=db_entry, inputFname=inputFname)
				inputFileBasename = os.path.basename(inputFname)
				#moveFileIntoDBAffiliatedStorage() will also set db_entry.path
				exit_code = db.moveFileIntoDBAffiliatedStorage(db_entry=db_entry, filename=inputFileBasename, \
									inputDir=os.path.split(inputFname)[0], \
									outputDir=data_dir,\
									relativeOutputDir=None, shellCommand='cp -rL', \
									srcFilenameLs=self.srcFilenameLs, dstFilenameLs=self.dstFilenameLs,\
									constructRelativePathFunction=db_entry.constructRelativePath, data_dir=data_dir)
				
				#exit_code = self.copyAndReformatResultFile(db, inputFname, db_entry, user=user, output_fname=localAbsPath)
			if exit_code==0:
				session.add(db_entry)
				if db_entry.file_size is None:
					db.updateDBEntryPathFileSize(db_entry=db_entry, data_dir=data_dir)
				if db_entry.md5sum is None:
					db.updateDBEntryMD5SUM(db_entry=db_entry, data_dir=data_dir)
				
				# 2010-5-3 store the json structure of top 10000 SNPs from the db_entry into db
				no_of_top_snps = 10000
				if db_entry.analysis_method.min_maf is not None:
					min_MAF = db_entry.analysis_method.min_maf
				else:
					min_MAF = 0
				try:
					#2011-2-24
					if call_method_id:	#call method, snp dataset
						db_id2chr_pos = db.snp_id2chr_pos
					elif cnv_method_id:
						if db._cnv_method_id!=cnv_method_id:
							db.cnv_id2chr_pos = cnv_method_id
						db_id2chr_pos = db.cnv_id2chr_pos
					pdata = PassingData(db_id2chr_pos=db_id2chr_pos)
					json_data = db.getOneResultJsonData(result_id=db_entry.id, min_MAF=min_MAF, no_of_top_snps=no_of_top_snps, \
													pdata=pdata, data_dir=data_dir)	#2011-2-24 pass pdata to getOneResultJsonData()
					rm_json = Stock_250kDB.ResultsMethodJson(min_maf=min_MAF, no_of_top_snps=no_of_top_snps)
					rm_json.result = db_entry
					rm_json.json_data = json_data
					session.add(rm_json)
				except:
					sys.stderr.write('Except in saving results_method_json (aborted): %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
					session.rollback()
					self.cleanUpAndExitOnFailure(exitCode=3)
			else:	#bad thing happend when getting data out of the file. don't save this results_method.
				session.delete(db_entry)
				sys.stderr.write("Error: copy file from %s to %s failed.\n"%(inputFname, localAbsPath ))
				session.rollback()
				self.cleanUpAndExitOnFailure(exitCode=3)
			session.flush()
			session.commit()
		else:	#default is also rollback(). to demonstrate good programming
			session.rollback()
		self.reset_marker_pos2snp_id()
	
	def run(self):
		"""
		2008-07-15
			adjust to new Elixir-based db api.
		2008-04-28
			use Stock_250kDatabase to do database stuff
		2008-04-16
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
				
		#db = Stock_250kDatabase(username=self.user,
		#		   password=self.passwd, hostname=self.hostname, database=self.dbname)
		if not os.path.isfile(self.inputFname):
			sys.stderr.write("Error: file, %s,  is not a file.\n"%(self.inputFname))
			sys.exit(3)
		self.add2DB(db=self.db_250k, short_name=self.short_name, phenotype_method_id=self.phenotype_method_id, \
				call_method_id=self.call_method_id, data_description=self.data_description, \
				method_description=self.method_description, comment=self.comment, inputFname=self.inputFname, user=self.db_user,
				results_method_type_id=self.results_method_type_id,\
				analysis_method_id=self.analysis_method_id, data_dir=self.data_dir, commit=self.commit,
				cnv_method_id=self.cnv_method_id)
		#2012.6.5
		self.outputLogMessage("submission done.\n")

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = Results2DB_250k
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()