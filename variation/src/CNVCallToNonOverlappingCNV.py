#!/usr/bin/env python
"""

Examples:
	# 2010-8-2
	CNVCallToNonOverlappingCNV.py -u yh -a 16 -m 29 -z banyan.usc.edu -c -r
	
	# 2010-10-11 work on chromosome 2
	CNVCallToNonOverlappingCNV.py -u yh -a 12 -m 30 -z banyan.usc.edu -c -r -o 2
	
Description:
	If two CNVs share much of the overlap (>=min_reciprocal_overlap), they will be split into
		non-overlapping new CNVs. This algorithm works iteratively until no overlapping CNVs exist in the RBTree.
	
	This program works on cnv data from table CNVCall.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import getopt, numpy
from variation.src.common import nt2number, number2nt
from pymodule import ProcessOptions, PassingData
import Stock_250kDB
from pymodule.CNV import get_overlap_ratio
from pymodule.CNV import CNVSegmentBinarySearchTreeKey
from pymodule.RBTree import RBDict
from CNVToNonOverlapping import CNVToNonOverlapping

class CNVCallToNonOverlappingCNV(CNVToNonOverlapping):
	__doc__ = __doc__
	"""
	2009-2-12
	"""
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
						('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('raw_cnv_method_id', 1, int): [8, 'a', 1, 'CNV method id of the old CNVs to be split'],\
						('chromosome', 1, int): [0, 'o', 1, 'which chromosome to work on. if =0, all chromosomes.'],\
						('cnv_method_id', 1, int): [10, 'm', 1, 'CNV method id for the non-overlapping CNVs.'],\
						('cnv_type_id', 1, int): [1, 'y', 1, 'CNV type id. table CNVType', ],\
						('min_reciprocal_overlap', 1, float): [0.000000001, 'n', 1, 'upper bound of overlap ratios for new CNVs. \
							0.00000001 means on overlapping is tolerated.'],\
						('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self, **keywords):
		"""
		2009-2-12
		"""
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
												class_to_have_attr=self)
		if self.chromosome ==0 or self.chromosome=='0':
			self.chromosome = None
	
	
	def assignNewCNVArrayCall(self, nonOverlappingCNVRBDict, raw_cnv_method_id=None, \
								cnv_method_id=None, param_obj=None):
		"""
		2010-8-2
		"""
		if not hasattr(param_obj, 'no_of_total'):
			setattr(param_obj, 'no_of_total', 0)
		if not hasattr(param_obj, 'no_of_into_db'):
			setattr(param_obj, 'no_of_into_db', 0)
		if not hasattr(param_obj, 'no_of_duplicate_cnv_array_calls'):
			setattr(param_obj, 'no_of_duplicate_cnv_array_calls', 0)
		cnv_type_id = getattr(param_obj, 'cnv_type_id', 1)
		cnv_method_id = getattr(param_obj, 'cnv_method_id', 1)
		
		session = getattr(param_obj, 'session', None)
		if session is None:
			sys.stderr.write("Error: db session is not available.\n")
			return
		
		for cnvNode in nonOverlappingCNVRBDict:
			param_obj.no_of_total += 1
			chromosome = cnvNode.key.chromosome
			start = cnvNode.key.span_ls[0]
			stop = cnvNode.key.span_ls[1]
			rows = Stock_250kDB.CNV.query.filter_by(chromosome=cnvNode.key.chromosome).filter_by(start=start).\
					filter_by(stop = stop).filter_by(cnv_type_id = cnv_type_id).\
					filter_by(cnv_method_id=cnv_method_id)
			parent_cnv_id_ls = cnvNode.key.parent_cnv_id_ls
			
			array_id2probability_ls = {}
			
			if rows.count()==0:	#make sure it's not in db yet.
				size_affected = stop - start + 1
				frequency = cnvNode.key.frequency
				#score = numpy.median(segment_matrix[:,3])
				
				cnv = Stock_250kDB.CNV(chromosome = chromosome, start=start, stop=stop, size_affected=size_affected,\
									frequency=frequency, cnv_method_id=cnv_method_id, cnv_type_id=cnv_type_id,)
				
				#record the sources for this CNV
				all_probability_ls = []
				for parent_cnv_call_id in parent_cnv_id_ls:
					parent_cnv_call = Stock_250kDB.CNVCall.get(parent_cnv_call_id)
					cnv.cnv_call_ls.append(parent_cnv_call)
					array_id = parent_cnv_call.array_id
					if array_id not in array_id2probability_ls:
						array_id2probability_ls[array_id] = []
					array_id2probability_ls[array_id].append(parent_cnv_call.probability)
					all_probability_ls.append(parent_cnv_call.probability)
				score = numpy.median(all_probability_ls)
				score_std = numpy.std(all_probability_ls)
				cnv.score = score
				cnv.score_std = score_std
				session.add(cnv)
				#param_obj.no_of_into_db += 1
				
				"""
				if param_obj.no_of_into_db%10000==0:
					try:
						session.flush()
						#session.expunge_all()	#don't expunge now because cnv object needs to be used next.
					except:
						sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
						import traceback
						traceback.print_exc()
						#import pdb
						#pdb.set_trace()
				"""
			else:
				cnv = rows.first()
			
			array_id2cnv_id = {}	#2010-8-5 make sure no array gets >1 entries into CNVArrayCall
			for array_id, probability_ls in array_id2probability_ls.iteritems():
				# 2010-8-4 move all the old CNVArrayCall to the new cnv_method_id
					if array_id in array_id2cnv_id:	#skip this array
						sys.stderr.write("Array %s already saved for this non-overlapping CNV (%s, %s-%s).\n"%\
										(array_id, cnv.chromosome, cnv.start, cnv.stop))
						#import pdb
						#pdb.set_trace()
						param_obj.no_of_duplicate_cnv_array_calls += 1
						continue
					else:
						array_id2cnv_id[array_id] = cnv.id
					saveNewObj = False
					if cnv.id:	#check whether it's in db already or not.
						rows = Stock_250kDB.CNVArrayCall.query.filter_by(array_id=array_id).filter_by(cnv_id=cnv.id).\
							filter_by(cnv_method_id=cnv_method_id)
						if rows.count()==0:
							saveNewObj = True
					else:
						saveNewObj = True
					if saveNewObj:
						probability_ls = array_id2probability_ls.get(array_id)
						if probability_ls:
							probability = numpy.median(probability_ls)
						else:
							probability = None
						cnv_array_call = Stock_250kDB.CNVArrayCall(array_id=array_id, cnv_method_id=cnv_method_id,\
																score=probability)
						cnv_array_call.cnv = cnv
						session.add(cnv_array_call)
						param_obj.no_of_into_db += 1
						if param_obj.no_of_into_db%10000==0:
							try:
								session.flush()
								#session.expunge_all()
							except:
								sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
								import traceback
								traceback.print_exc()
								#import pdb
								#pdb.set_trace()
			
			report = getattr(param_obj, 'report', False)
			if report and param_obj.no_of_total%10000==0:
				sys.stderr.write('%s%s into db out of %s. %s duplicate CNVArrayCalls.\n'%('\x08'*100, param_obj.no_of_into_db, \
														param_obj.no_of_total, param_obj.no_of_duplicate_cnv_array_calls))
		sys.stderr.write('%s%s into db out of %s. %s duplicate CNVArrayCalls.\n'%('\x08'*100, param_obj.no_of_into_db, \
														param_obj.no_of_total, param_obj.no_of_duplicate_cnv_array_calls))
	
	def run(self):
		"""
		2010-7-29
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_250k = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				  			password=self.db_passwd, hostname=self.hostname, database=self.dbname, 
				   			schema=self.schema)
		db_250k.setup(create_tables=False)
		session = db_250k.session
		session.expire_on_commit = False	#otherwise a flush() would cause:	
		"""
		sqlalchemy.orm.exc.DetachedInstanceError: Instance <...x770eb6d0> is not bound to a Session;
		attribute refresh operation cannot proceed
		
		# 2010-7-31 	session.expire_on_commit = False doesn't help
		"""
		from CNVMergeAcrossArrays import CNVMergeAcrossArrays
		non_duplicate_array_id_ls = CNVMergeAcrossArrays.getNonDuplicateArraysWithHighestMedianIntensity(db_250k, \
																					self.raw_cnv_method_id)
		no_of_total_arrays = len(non_duplicate_array_id_ls)
		frequency = 1.0/no_of_total_arrays
		nonOverlappingCNVRBDict = self.partitionCNVIntoNonOverlapping(db_250k, cnv_method_id=self.raw_cnv_method_id,
									table_name=Stock_250kDB.CNVCall.table.name, frequency=frequency, chromosome=self.chromosome)
		
		param_obj = PassingData(no_of_valid_deletions=0, session=db_250k.session, cnv_type_id=self.cnv_type_id, \
					cnv_method_id=self.cnv_method_id, no_of_total=0, no_of_into_db=0, report=self.report,)
		
		self.assignNewCNVArrayCall(nonOverlappingCNVRBDict, raw_cnv_method_id=self.raw_cnv_method_id, \
								cnv_method_id=self.cnv_method_id, param_obj=param_obj)
		
		
		session.flush()
		session.expunge_all()
		if self.commit:
			session.commit()
	
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = CNVCallToNonOverlappingCNV
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
