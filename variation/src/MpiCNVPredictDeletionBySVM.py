#!/usr/bin/env mpipython
"""

Examples:
	
Description:
	2010-7-1 parallel version of CNVPredictionBySVM.py
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import getopt, csv, math, traceback
import cPickle, traceback
from Scientific import MPI
from pymodule.MPIwrapper import mpi_synchronize, MPIwrapper
from pymodule import PassingData, getListOutOfStr, figureOutDelimiter
from RunGADA import RunGADA
#import Stock_250kDB
from sets import Set
import CNVPredictionBySVM

class MpiGADA(CNVPredictionBySVM, MPIwrapper):
	__doc__ = __doc__
	option_default_dict = CNVPredictionBySVM.option_default_dict.copy()
	
	def __init__(self, **keywords):
		"""
		2009-10-28
		"""
		CNVPredictionBySVM.__init__(**keywords)
	
	def generate_params(cls, aAlpha_ls, TBackElim_ls, MinSegLen_ls, reader=None, array_id_set=None, \
					param_obj=None, debug=None, chr2start_stop_index=None):
		"""
		2010-6-5
			split the whole genome into different chromosomes
			pass data of each chromosome to a computing node along with the start_index
		2010-5-28
			return yield (parameter generator)
		2009-10-27
		"""
		param_ls = []
		for aAlpha in aAlpha_ls:
			for TBackElim in TBackElim_ls:
				for MinSegLen in MinSegLen_ls:
						param_ls.append([aAlpha, TBackElim, MinSegLen])
		sys.stderr.write("%s different (a, T, M) settings.\n"%(len(param_ls)))
		
		for row in reader:
			array_id = int(row[0])
			if array_id_set and array_id not in array_id_set:	#skip if array_id_set is specified and not in array_id_set
				continue
			ecotype_id = row[1]
			intensity_ls = row[2:]
			intensity_array = map(float, intensity_ls)
			for param in param_ls:
				for chr, start_stop_index in chr2start_stop_index.iteritems():
					start_index, stop_index = start_stop_index[:2]
					yield (param + [array_id, ecotype_id, start_index], intensity_array[start_index:stop_index+1])
	
	def computing_node_handler(self, communicator, data, param_obj):
		"""
		2009-11-23
			if self.IO_thru_file:
				call self._GADA()
			else:
				call GADA.GADA.run()
		2009-11-10
			GADA no longer uses PIPE (stdin, stdout). Specify input and output file directly.
		2009-10-28
		
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		result_ls = []
		for one_data in data:
			one_parameter, intensity_ls = one_data[:2]
			aAlpha, TBackElim, MinSegLen, array_id = one_parameter[:4]
			if self.IO_thru_file:
				tmp_input_fname = '%s_%s_A%sT%sM%s'%(param_obj.tmp_input_fname, array_id, aAlpha, TBackElim, MinSegLen)
				input_file = self.prepareGADAinput(intensity_ls, tmp_input_fname, IO_thru_file=self.IO_thru_file)
				tmp_output_fname = '%s_%s_A%sT%sM%s'%(param_obj.tmp_output_fname, array_id, aAlpha, TBackElim, MinSegLen)
				GADA_output = self._GADA(self.GADA_path, input_file, tmp_output_fname, aAlpha, \
										TBackElim, MinSegLen)
			else:
				ins = GADA.GADA(self.debug)
				GADA_output = ins.run(map(float, intensity_ls), aAlpha, TBackElim, MinSegLen)	#intensity_ls cell type is numpy.float32.
					#convert it to python float, which is recognized by GADA.GADA.run().
				self.addMedianToGADAOutput(GADA_output, intensity_ls)
				del ins
			if GADA_output:
				result_ls.append((one_parameter, GADA_output))
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return result_ls
	
	
	def output_node_handler(self, communicator, output_param_obj, data):
		"""
		2010-6-5
			pass base_start_index from one_parameter to output_GADA_output()
		2010-6-3
			ecotype_id is passed from the computing node. No need to do db fetching.
		2009-10-28
		"""
		result_ls = cPickle.loads(data)
		counter = 0
		for result in result_ls:
			one_parameter, GADA_output = result
			aAlpha, TBackElim, MinSegLen, array_id, ecotype_id, base_start_index = one_parameter[:6]
			param_combo = (aAlpha, TBackElim, MinSegLen)
			if param_combo not in output_param_obj.param_combo2writer:
				output_fname_prefix = os.path.splitext(output_param_obj.output_fname_prefix)[0]
				output_fname = '%s_A%sT%sM%s.tsv'%(output_fname_prefix, aAlpha, TBackElim, MinSegLen)
				writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
				self.output_header(writer)
				output_param_obj.param_combo2writer[param_combo] = writer
				output_param_obj.param_combo2array_id2no_of_segments[param_combo] = {}
			writer = output_param_obj.param_combo2writer[param_combo]
			no_of_segments = self.output_GADA_output(writer, GADA_output, array_id, ecotype_id, output_param_obj.chr_pos_ls, \
													output_param_obj.probe_id_ls, base_start_index=base_start_index)
			output_param_obj.param_combo2array_id2no_of_segments[param_combo][array_id] = no_of_segments
			counter += 1
		sys.stderr.write("%s GADA results were outputted.\n"%counter)
	
	def run(self):
		"""
		2009-10-28
		"""
		if self.debug:
			#for one-node testing purpose
			import pdb
			pdb.set_trace()
		
		self.communicator = MPI.world.duplicate()
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		free_computing_node_set = Set(free_computing_nodes)
		output_node_rank = self.communicator.size-1
		
		if node_rank == 0:
			db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				  				password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
			db.setup(create_tables=False)
			
			array_id2median_intensity = self.get_array_id2median_intensity(min_array_median_intensity=self.min_array_median_intensity)
			array_id2model = self.constructSVMModels(db, self.arrays_to_form_model, array_id2median_intensity,\
						minPercUnCoveredByLerContig=self.minPercUnCoveredByLerContig, cnv_method_id=self.training_cnv_method_id,\
						C=self.SVM_C, gamma=self.SVM_gamma, eps=self.SVM_eps)
			array_id2model_array_id = self.mapAnyArray2ModelArray(array_id2median_intensity, array_id2model, \
															max_median_intensity_dist=self.max_median_intensity_dist)
			
			commonData_pickle = cPickle.dumps(commonData, protocol=-1)
			sys.stderr.write("Passing data to output node %s from %s ... "%(output_node_rank, node_rank,))
			self.communicator.send(commonData_pickle, output_node_rank, 0)
			sys.stderr.write(".\n")
			
			refDataMatrix_pickle = cPickle.dumps(refDataMatrix, protocol=-1)
			for node in free_computing_nodes:	#send it to the computing_node
				sys.stderr.write("passing initial data to nodes from %s to %s ... "%(node_rank, node))
				self.communicator.send(commonData_pickle, node, 0)
				self.communicator.send(refDataMatrix_pickle, node, 0)
				sys.stderr.write(".\n")
			if len(commonData.blockDataCodedIndex_ls)==0:
				sys.stderr.write("Not a single block is formed. Exit!")
				sys.exit(0)
			del commonData, commonData_pickle, refDataMatrix, refDataMatrix_pickle
			
			
			param_ls = self.generate_params(self.aAlpha_ls, self.TBackElim_ls, self.MinSegLen_ls, array_id_set=array_id_set,\
											reader=inputWrapper.reader, debug=self.debug, \
											chr2start_stop_index=chr2start_stop_index)
		elif node_rank in free_computing_node_set:			
			pass
		else:
			# 2010-6-3 output node doesn't need db connection anymore
			db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
			db.setup(create_tables=False)
			session = db.session
			data, source, tag = self.communicator.receiveString(0, 0)
			passing_to_output_node_obj =  cPickle.loads(data)
			del data
			sys.stderr.write(".\n")
		
		self.synchronize()
		if node_rank == 0:
			#parameter_list = [param_ls, array_id2col_index, data_matrix]
			#self.input_node(parameter_list, free_computing_nodes, input_handler=self.input_handler, message_size=self.message_size)
			param_obj = PassingData(param_ls=param_ls, output_node_rank=output_node_rank, report=self.report, counter=0)
			self.inputNode(param_obj, free_computing_nodes, param_generator = param_ls, message_size=self.message_size)
		elif node_rank in free_computing_node_set:
			if self.IO_thru_file:
				tmp_input_dir = os.path.split(self.tmp_input_fname)[0]
				tmp_output_dir = os.path.split(self.tmp_output_fname)[0]
				try:
					if not os.path.isdir(tmp_input_dir):
						os.makedirs(tmp_input_dir)
					if not os.path.isdir(tmp_output_dir):
						os.makedirs(tmp_output_dir)
				except:
					 traceback.print_exc()
					 sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			computing_parameter_obj = PassingData(tmp_input_fname=self.tmp_input_fname, tmp_output_fname=self.tmp_output_fname)
			self.computing_node(computing_parameter_obj, self.computing_node_handler)
		else:
			output_dir = os.path.split(self.output_fname_prefix)[0]
			if not os.path.isdir(output_dir):
				os.makedirs(output_dir)
			param_combo2writer = {}
			param_combo2array_id2no_of_segments = {}
			output_param_obj = PassingData(param_combo2writer=param_combo2writer, output_fname_prefix=self.output_fname_prefix,\
										chr_pos_ls = passing_to_output_node_obj.chr_pos_ls,\
										probe_id_ls = passing_to_output_node_obj.probe_id_ls,\
										param_combo2array_id2no_of_segments = param_combo2array_id2no_of_segments)
			self.output_node(free_computing_nodes, output_param_obj, self.output_node_handler)
			del param_combo2writer
			self.outputSegmentStat(param_combo2array_id2no_of_segments)
		self.synchronize()	#to avoid some node early exits

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiGADA
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()