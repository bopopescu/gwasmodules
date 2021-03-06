#!/usr/bin/env python
"""

Examples:
	#output data with accession mismatch rate<=0.20, snp mismatch rate <=0.15, snp NA rate <= default
	DB_250k2data.py -z localhost -l 3 -y 0.85 -w 0.15 -x 0.20 -o /tmp/250k_l3_w0.15_x0.20_y0.85.tsv 
	
	#output matrix with no SNP filtering -w 1 -m 1
	DB_250k2data.py -l 3 -y 0.85 -w 1 -m 1 -x 0.20 -o /tmp/250k_l3_v1_w1_x0.20_y0.85.tsv
	
	DB_250k2data.py -l 17 -o /mnt/nfs/250k/call_method_17.tsv
	
	#no filtering for either SNP or array, but take unique ecotype, the input data directory is changed.
	DB_250k2data.py -y 0.85 -w 1 -m1 -x1 -l 3 -t -o 250k_l3_y.85_uniq_ecotype_20080919.tsv -i ~/banyan_fs/Network/Data/250k/db/calls/method_3/
	
Description:
	2008-05-06
		program to output/filter 250k data based on QC recorded in database in Strain X SNP format.
	
	2011-2-15
	Input: 1st column is Snps.id, 2nd column is ATCG call. 1st row is header. Data starts from 2nd row.
	Output format: 1st column is array id. 2nd column is ecotype id. 1st row is Snps.id. Data starts from 2nd row.

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
import warnings, traceback
import sqlalchemy
from pymodule import write_data_matrix
from pymodule import ProcessOptions, number2nt, nt2number
from variation.src import Stock_250kDB
from variation.src.qc.QualityControl import QualityControl
from variation.src.qc.QC_250k import QC_250k
from variation.src.mapper.AbstractVariationMapper import AbstractVariationMapper

class DB_250k2Data(AbstractVariationMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVariationMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 0, ))
	#option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
							('call_method_id', 1, int): [None, 'l', 1, 'id in table call_method', ],\
							('input_dir', 0, ): [None, 'i', 1, 'If given is directory, call_info.filename is assumed to be in this directory.'],\
							('min_probability', 0, float): [-1, 'y', 1, 'minimum probability for a call to be non-NA if there IS a 3rd column for probability.', ],\
							('max_array_mismatch_rate', 0, float): [1, 'x', 1, 'maximum mismatch rate of an array call_info entry. used to exclude bad arrays.'],\
							('max_snp_mismatch_rate', 0, float): [1, 'w', 1, 'maximum snp error rate, used to exclude bad SNPs', ],\
							('max_snp_NA_rate', 1, float): [1, 'm', 1, 'maximum snp NA rate, used to exclude SNPs with too many NAs', ],\
							('take_unique_ecotype', 0, int):[0, 't', 0, 'toggle this to keep only one call_info with lowest mismatch_rate for one ecotype'],\
							})
	"""
	2008-09-19
		add option take_unique_ecotype
	"""
	def __init__(self, **keywords):
		"""
		2008-05-06
		"""
		AbstractVariationMapper.__init__(self, inputFnameLs=None, **keywords)
		
		
	def get_snps_with_best_QC_ls(cls, db, call_method_id):
		"""
		2008-05-06
			deprecated, too much memory
		"""
		sys.stderr.write("Getting snps_with_best_QC_ls ... ")
		snps_ls = Snps.query.options(sqlalchemy.orm.eagerload('snps_qc')).join('snps_qc').filter(db.tables['snps_qc'].c.call_method_id==call_method_id).all()
		
		snps_with_best_QC_ls = []
		no_of_entries = len(snps_ls)
		for i in range(no_of_entries):
			snps = snps_ls[i]
			if i%5000==0:
				sys.stderr.write("%s%d/%d:\t\t%s"%('\x08'*100, i+1, no_of_entries, snps.name))
			#choose the snps_QC with maximum no of non-NA pairs to get mismatch_rate
			if snps.snps_QC:
				snps_QC_with_max_no_of_non_NA_pairs = snps.snps_QC[0]
				for snps_QC in snps.snps_QC:
					if snps_QC.no_of_non_NA_pairs>snps_QC_with_max_no_of_non_NA_pairs.no_of_non_NA_pairs:
						snps_QC_with_max_no_of_non_NA_pairs = snps_QC
				snps.snps_QC_with_max_no_of_non_NA_pairs = snps_QC_with_max_no_of_non_NA_pairs				
			else:
				snps.snps_QC_with_max_no_of_non_NA_pairs = None
			snps_with_best_QC_ls.append(snps)
		
		sys.stderr.write("%s SNPs. Done.\n"%(len(snps_with_best_QC_ls)))
		return snps_with_best_QC_ls
	
	get_snps_with_best_QC_ls= classmethod(get_snps_with_best_QC_ls)
	
	def get_snps_name_set_given_criteria(cls, db, call_method_id, max_snp_mismatch_rate, max_snp_NA_rate):
		"""
		2008-05-17
			add chromosome_position to snps_name_set as well due to some call files in file system storage only have that as SNP id
		2008-05-07
			direct sql select, less memory
		2008-05-06
			scan thru all QCs related to one SNP and judge based on the QC with most no of non-NA pairs
			
			cost too much memory
		"""
		sys.stderr.write("Getting snps_name_set ... \n")
		snps_name_set = set()
		s = db.tables['snps'].alias()
		q = db.tables['snps_qc'].alias()
		sql_sentence = sqlalchemy.sql.select([s.c.id, s.c.name, q.c.NA_rate, q.c.mismatch_rate, q.c.no_of_non_NA_pairs], s.c.id==q.c.snps_id, order_by=[s.c.id])
		
		block_size = 5000
		counter = 0
		
		results = db.connection.execute(sql_sentence)
		rows = results.fetchmany(block_size)	#2008-05-07 don't fetchmany() is more memory-efficient than fetchall(). 
		old_row = None
		while rows:
			for row in rows:
				counter += 1
				if old_row==None:
					old_row = row
				elif row.id == old_row.id:
					if row.no_of_non_NA_pairs>old_row.no_of_non_NA_pairs:
						old_row = row
				elif row.id != old_row.id:
					if old_row.NA_rate<=max_snp_NA_rate and old_row.mismatch_rate<=max_snp_mismatch_rate:
						snps_name_set.add(old_row.name)
						snps_name_set.add(old_row.name[:-4])	#2008-05-17 add chromosome_position to it as well
					old_row = row
			sys.stderr.write("%s%d"%('\x08'*100, counter))
			rows = results.fetchmany(block_size)
		
		#take care of the last one
		if old_row.NA_rate<=max_snp_NA_rate and old_row.mismatch_rate<=max_snp_mismatch_rate:
			snps_name_set.add(old_row.name)
			snps_name_set.add(old_row.name[:-4])
		"""
		#2008-05-07 cost too much memory
		snps_ls = db.session.query(SNPs).options(sqlalchemy.orm.eagerload('snps_QC')).offset(0).limit(block_size).list()
		no_of_entries = len(snps_ls)
		i = 0
		while no_of_entries>0:
			for j in range(no_of_entries):
				counter += 1
				snps = snps_ls[j]
				#choose the snps_QC with maximum no of non-NA pairs to get mismatch_rate
				if snps.snps_QC:
					snps_QC_with_max_no_of_non_NA_pairs = snps.snps_QC[0]
					for k in range(1,len(snps.snps_QC)):
						snps_QC = snps.snps_QC[k]
						if snps_QC.no_of_non_NA_pairs>snps_QC_with_max_no_of_non_NA_pairs.no_of_non_NA_pairs:
							snps_QC_with_max_no_of_non_NA_pairs = snps_QC
					if snps_QC_with_max_no_of_non_NA_pairs.NA_rate<=max_snp_NA_rate and snps_QC_with_max_no_of_non_NA_pairs.mismatch_rate<=max_snp_mismatch_rate:
						snps_name_set.add(snps.name)
			del snps_ls	#release memory
			sys.stderr.write("%s%d"%('\x08'*100, counter))
			i += 1
			snps_ls = db.session.query(SNPs).options(sqlalchemy.orm.eagerload('snps_QC')).offset(i*block_size).limit(block_size).list()
			no_of_entries = len(snps_ls)
		del snps_ls
		"""
		sys.stderr.write("\n %s snps names Done.\n"%(len(snps_name_set)))
		return snps_name_set
	
	get_snps_name_set_given_criteria= classmethod(get_snps_name_set_given_criteria)
	
	
	@classmethod
	def getSNPID2index(cls, call_info_fname, db_id2chr_pos):
		"""
		2010-10-13
			based on Snps.id listed in first column of call_info_fname and db_id2chr_pos, return a dictionary of Snps.id:index
				and the snps are in (chr,pos) order.
			
		"""
		sys.stderr.write("Getting snp_id2index ...")
		reader = csv.reader(open(call_info_fname), delimiter='\t')
		header = reader.next()
		chr_pos_db_id_ls = []
		for row in reader:
			db_id = row[0]
			try:
				db_id = int(db_id)
				chr_pos = db_id2chr_pos.get(db_id)
				if chr_pos is not None:	#sometime this could be missing. ignore this snp then.
					chr_pos_db_id_ls.append([chr_pos, db_id])
			except:
				continue
		chr_pos_db_id_ls.sort()	#sort in chr, pos order
		
		db_id2index = {}
		for i in xrange(len(chr_pos_db_id_ls)):
			chr_pos, db_id = chr_pos_db_id_ls[i]
			db_id2index[db_id] = i
		sys.stderr.write(" %s entries. Done.\n"%(len(db_id2index)))
		return db_id2index
		
	def run(self):
		"""
		2008-05-20 read_call_matrix returns PassingData object
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db = self.db_250k
		session = db.session
		QC_method_id = 0 	#just for QC_250k.get_call_info_id2fname()
		call_data = QC_250k.get_call_info_id2fname(db, QC_method_id, self.call_method_id, filter_calls_QCed=0, \
												max_call_info_mismatch_rate=self.max_array_mismatch_rate, input_dir=self.input_dir,\
												take_unique_ecotype=self.take_unique_ecotype)
		#snps_with_best_QC_ls = self.get_snps_with_best_QC_ls(db, self.call_method_id)
		if self.max_snp_mismatch_rate<1 or self.max_snp_NA_rate<1:	#2008-05-18 only do this when it's necessary
			snps_name_set = self.get_snps_name_set_given_criteria(db, self.call_method_id, self.max_snp_mismatch_rate, self.max_snp_NA_rate)
		else:
			snps_name_set = None
		db_id2chr_pos = db.getSNPID2ChrPos()
		if len(call_data.call_info_id2fname)>0:
			db_id2index = self.getSNPID2index(call_data.call_info_id2fname.values()[0][1], db_id2chr_pos)
			pdata = QC_250k.read_call_matrix(call_data.call_info_id2fname, self.min_probability, snps_name_set, \
											db_id2chr_pos=db_id2chr_pos, db_id2index=db_id2index)	#2008-05-20 read_call_matrix returns PassingData object
			strain_acc_list, category_list = pdata.ecotype_id_ls, pdata.array_id_ls
			write_data_matrix(pdata.data_matrix, self.outputFname, pdata.header, strain_acc_list, category_list)

if __name__ == '__main__':
	main_class = DB_250k2Data
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()	
