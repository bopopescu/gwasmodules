#!/usr/bin/env mpipython
"""

Examples:
	#run it on hpc-cmb cluster
	mpiexec ~/script/variation/src/MpiLDCNVvsSNP.py

	#test parallel run on desktop
	mpirun -np 5 -machinefile  /tmp/hostfile /usr/bin/mpipython ~/script/...
	mpirun -np 5 -machinefile ~/hostfile /usr/bin/mpipython ./script/variation/src/MpiLDCNVvsSNP.py -m 20 -c 0.0
		-i /Network/Data/250k/db/dataset/call_method_32.tsv  -o ~/cnvMethod20_vs_callMethod32_LD.tsv -u yh -p xxx -z banyan
	
	# debug with one node
	MpiLDCNVvsSNP.py -m 20 -c 0.0 -i /Network/Data/250k/tmp-yh/250k_data/call_method_17_test.tsv 
		-o /tmp/cnvMethod20_snp_LD.tsv -u yh -p xxx -z banyan -b
Description:
	2010-9-30 a MPI program to calculate LD between CNV and its nearby SNPs.
		The input SNP format is StrainXSNP matrix, either tab or comma-delimited.
			1st row as header. 1st column as strain id, 2nd column is affiliating information for each strain.
		CNV data is fetched from db.
	
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import getopt, csv, math, random
import cPickle
from pymodule import PassingData, importNumericArray, SNPData, read_data
from Scientific import MPI
from pymodule.MPIwrapper import MPIwrapper
from pymodule.CNV import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
from pymodule.RBTree import RBDict
import Stock_250kDB
from pymodule.algorithm import LD
import zlib
from CNVMergeAcrossArrays import CNVMergeAcrossArrays

num = importNumericArray()

class MpiLDCNVvsSNP(MPIwrapper):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
						('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('input_fname',1, ): [None, 'i', 1, 'a file containing StrainXSNP matrix.'],\
						("output_fname", 1, ): [None, 'o', 1, 'Filename to store data matrix'],\
						("cnv_method_id", 1, ): [None, 'm', 1, 'CNV method id of the CNV data (the non-overlapping ones)'],\
						('message_size', 1, int):[1000, 's', 1, 'number of CNVs to computer its LD with nearby SNPs'],\
						("min_LD_to_output", 1, float): [0, 'L', 1, "output LD data whose |D'| is above this cutoff"],\
						('min_MAF', 1, float): [0.0, 'n', 1, 'minimum Minor Allele Frequency for both SNPs, filter after LD is calculated'],\
						('max_CNV_SNP_dist', 1, float): [20000, 'X', 1, 'maximum distance between a SNP and a CNV'],\
						('discard_perc', 1, float): [0, 'c', 1, 'percentage of SNP Pairs to be randomly discarded. \
							this filter is applied before min_LD_to_output and min_MAF and LD calculation'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-09-29
			add option min_LD_to_output and min_MAF
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def createCNVRBDict(self, db_250k, cnv_method_id=None, max_CNV_SNP_dist=None, array_id2row_index = None, snp_id_ls = []):
		"""
		2010-9-30
			This function is to get 1. CNVs from cnv_method_id, 2 the nearby SNPs for each CNV.
				create a RBDict based on CNV segments (add max_CNV_SNP_dist on each side).
				for each SNP, find out CNV segments which contain it.
		"""
		sys.stderr.write("Creating CNV RBDict ... \n")
		query = Stock_250kDB.CNV.query.filter_by(cnv_method_id=cnv_method_id)
		CNVRBDict = RBDict()
		count = 0
		real_count = 0
		for row in query:
			count += 1
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=row.chromosome, \
							span_ls=[max(1, row.start - max_CNV_SNP_dist), row.stop + max_CNV_SNP_dist], \
							min_reciprocal_overlap=1, cnv_id=row.id, cnv_start=row.start, cnv_stop=row.stop)
							#2010-8-17 any overlap short of identity is tolerated.
			if segmentKey not in CNVRBDict:
				CNVRBDict[segmentKey] = PassingData(snp_id_ls=[], deletionDataLs = [0]*len(array_id2row_index))
			"""
			# 2010-9-30 too much memory
			for cnv_array_call in row.cnv_array_call_ls:
				array_row_index = array_id2row_index.get(cnv_array_call.array_id)
				if array_row_index is not None:	#ignore arrays not in SNPs
					CNVRBDict[segmentKey].deletionDataLs[array_row_index] = 1
					real_count += 1
			"""
			if count%200==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, count, real_count))
				if self.debug:
					break
		sys.stderr.write("%s%s\t%s\n %s Done.\n"%('\x08'*80, count, real_count, repr(CNVRBDict)))
		
		sys.stderr.write("Finding nearby SNPs for CNVs ... \n")
		compareIns = CNVCompare(min_reciprocal_overlap=0.0000001)	#any overlap is an overlap
		count = 0
		real_count = 0
		for snp_id in snp_id_ls:
			chromosome, start = snp_id.split('_')[:2]
			chromosome = int(chromosome)
			start = int(start)
			snpSegmentKey = CNVSegmentBinarySearchTreeKey(chromosome=chromosome, \
							span_ls=[start], \
							min_reciprocal_overlap=0.0000001, )	#min_reciprocal_overlap doesn't matter here.
			node_ls = []
			CNVRBDict.findNodes(snpSegmentKey, node_ls=node_ls, compareIns=compareIns)
			for node in node_ls:
				cnvSegKey = node.key
				node.value.snp_id_ls.append(snp_id)
				real_count += 1
			count += 1
			if count%1000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, count, real_count))
		sys.stderr.write("%s%s\t%s\n Done.\n"%('\x08'*80, count, real_count))
		return CNVRBDict
	
	def generate_params(self, CNVRBDict):
		"""
		"""
		for node in CNVRBDict:
			yield node
		

	
	def computing_node_handler(self, communicator, data, param_obj):
		"""
		2009-3-21
			1. data = [(min_index1,stop1), (min_index2,stop2)]. remove one for loop.
			2. add another filter, discard_perc. discard a certain percentage of SNP Pairs randomly
		2008-09-29
			filter LD using |D'| >= min_LD_to_output
			both SNPs have to pass min_MAF
		2008-09-12
			filter r2 below min_r2_to_output
		2008-09-06
			data from input_node is changed
		2008-09-05
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		node_ls = cPickle.loads(data)
		result_ls = []
		snpData = param_obj.snpData
		array_id2row_index = param_obj.array_id2row_index
		
		for node in node_ls:
			cnvSegKey = node.key
			query = Stock_250kDB.CNVArrayCall.query.filter_by(cnv_id=cnvSegKey.cnv_id)
			deletionDataLs = node.value.deletionDataLs
			for cnv_array_call in query:
				array_row_index = array_id2row_index.get(cnv_array_call.array_id)
				if array_row_index is not None:	#ignore arrays not in SNPs
					deletionDataLs[array_row_index] = 1
			snp_id_ls = node.value.snp_id_ls
			for snp_id in snp_id_ls:
				if param_obj.discard_perc==0:
					u = 1
				else:	#sample a uniform unless discard_perc is not 0
					u = random.random()
				if u>=param_obj.discard_perc:
					col_index = snpData.col_id2col_index[snp_id]
					snpAlleleLs = snpData.data_matrix[:, col_index]
					
					LD_data = LD.calLD(snpAlleleLs, deletionDataLs, locus1_id=snp_id, locus2_id = cnvSegKey.cnv_id)
					if LD_data is not None and LD_data.allele_freq[0]>=param_obj.min_MAF and LD_data.allele_freq[1]>=param_obj.min_MAF \
							and abs(LD_data.r2)>=param_obj.min_LD_to_output:
						result_ls.append(LD_data)
			
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return result_ls
	
	def output_node_handler(self, communicator, param_obj, data):
		"""
		2008-09-05
		"""
		writer = param_obj.writer
		result_ls = cPickle.loads(data)
		for LD_data in result_ls:
			if writer:
				row = []
				if not param_obj.is_header_written:
					header_row = []
					for column_name in LD_data.__dict__:
						if column_name=='allele_freq':
							header_row.append('allele1_freq')
							header_row.append('allele2_freq')
						elif column_name=='snp_pair_ls':
							header_row.append('snp1')
							header_row.append('snp2')
						else:
							header_row.append(column_name)
					writer.writerow(header_row)
					param_obj.is_header_written = True
				
				for key,value in LD_data.__dict__.iteritems():
					if key=='allele_freq':
						row.append(value[0])
						row.append(value[1])
					elif key=='snp_pair_ls':
						row.append(value[0])
						row.append(value[1])
					else:
						row.append(value)
				writer.writerow(row)
	
	def run(self):
		self.communicator = MPI.world.duplicate()
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		free_computing_node_set = set(free_computing_nodes)
		output_node_rank = self.communicator.size-1
		
		if node_rank == 0:
			if self.debug:
				#for one-node testing purpose
				import pdb
				pdb.set_trace()
			
			db_250k = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, 
						schema=self.schema)
			db_250k.setup(create_tables=False)
			session = db_250k.session
			
			# 2010-9-30 get total number of arrays in this CNV method 
			non_duplicate_array_id_ls = CNVMergeAcrossArrays.getNonDuplicateArraysWithHighestMedianIntensity(db_250k, \
										self.cnv_method_id, table_name=Stock_250kDB.CNVArrayCall.table.name)
			non_duplicate_array_id_set = set(non_duplicate_array_id_ls)
			no_of_total_arrays = len(non_duplicate_array_id_ls)
			
			# read in the SNP set with only arrays in the CNV method set
			snpData = SNPData(input_fname=self.input_fname, turn_into_array=1)
			row_index_to_be_kept = []
			for row_id, row_index in snpData.row_id2row_index.iteritems():
				array_id = int(row_id[1])
				if array_id in non_duplicate_array_id_set:
					row_index_to_be_kept.append(row_index)
			snpData = snpData.keepRowsByRowIndex(snpData, row_index_to_be_kept)
			# a map between array_id and its row index in the SNP dataset
			array_id2row_index = {}
			for row_id, row_index in snpData.row_id2row_index.iteritems():
				array_id = int(row_id[1])
				array_id2row_index[array_id] = row_index
			
			# create a map (RBDict) between each CNV and its nearby SNPs
			# get all CNVs from db
			CNVRBdict = self.createCNVRBDict(db_250k, self.cnv_method_id, self.max_CNV_SNP_dist, array_id2row_index = array_id2row_index, \
											snp_id_ls = snpData.col_id_ls)
			snpData.array_id2row_index = array_id2row_index	# passed to computer node later
			
			snpData_pickle = cPickle.dumps(snpData, -1)
			snpData_pickle = zlib.compress(snpData_pickle)
			for node in free_computing_nodes:	#send it to the computing_node
				sys.stderr.write("passing initial data to nodes from %s to %s ... "%(node_rank, node))
				self.communicator.send(snpData_pickle, node, 0)
				sys.stderr.write(".\n")
			del snpData_pickle
			del snpData
			params_ls = self.generate_params(CNVRBdict,)
		elif node_rank in free_computing_node_set:
			db_250k = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, 
						schema=self.schema)
			db_250k.setup(create_tables=False)
			session = db_250k.session
			
			data, source, tag = self.communicator.receiveString(0, 0)
			data = zlib.decompress(data)	# 2010-10-1 decompress
			snpData =  cPickle.loads(data)
			del data
		else:
			pass
		
		self.synchronize()
		if node_rank == 0:
			param_obj = PassingData(params_ls=params_ls, output_node_rank=output_node_rank, report=self.report, counter=0)
			self.inputNode(param_obj, free_computing_nodes, param_generator = params_ls, message_size=self.message_size)
		elif node_rank in free_computing_node_set:
			computing_parameter_obj = PassingData(snpData=snpData, min_LD_to_output=self.min_LD_to_output, \
									min_MAF=self.min_MAF, discard_perc=self.discard_perc, db_250k=db_250k, \
									array_id2row_index=snpData.array_id2row_index)
			self.computing_node(computing_parameter_obj, self.computing_node_handler)
		else:
			if getattr(self, 'output_fname', None):
				writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
			else:
				writer = None
			param_obj = PassingData(writer=writer, is_header_written=False)
			self.output_node(free_computing_nodes, param_obj, self.output_node_handler)
			del writer
		self.synchronize()	#to avoid some node early exits
	

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiLDCNVvsSNP
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
