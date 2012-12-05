#!/usr/bin/env mpipython
"""

Examples:
	#run it on hpc-cmb cluster
	mpiexec ~/script/variation/src/MpiWithinArrayBlockQuantileNormalize.py 
	
	#test parallel run on desktop
	mpirun -np 5 -machinefile  /tmp/hostfile /usr/bin/mpipython ~/script/...
	
	#to test and debug (serial, rather than parallel)
	python ~/script/variation/src/MpiWithinArrayBlockQuantileNormalize.py -a 3,4 -o /tmp/array_3_4.tsv -y 1200,1612 -x 0,200 -u yh -b 
	
	#run all arrays from call method 48 on hpc-cmb
	mpiexec ~/script/variation/src/MpiWithinArrayBlockQuantileNormalize.py \
		-l 48 -o /Network/Data/250k/tmp-yh/CNV/call_48_WABQN.tsv  -u yh -r -f ~/panfs/db/raw_data
	
Description:
	2010-5-24
		This program strives to remove array-specific spatial effect.
		For each array, partition into overlapping blocks, quantile-normalize intensity from all blocks.
			Take median for probes with data from >=1 overlapping blocks.
		Output is tab-delimited format. Each row is an array. Each column is a probe.
			First row is probe id.
			2nd row is chr.
			3rd row is position.
			First column is array id.
			2nd column is ecotype id.
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import getopt, csv, math
import cPickle
from pymodule import PassingData, importNumericArray, SNPData, read_data, getListOutOfStr, figureOutDelimiter
from pymodule.SNP import NA_set
from sets import Set
from Scientific import MPI
from pymodule.MPIwrapper import MPIwrapper
import Stock_250kDB, numpy, random, cPickle
from DB_250k2Array import DB_250k2Array, probes_class
from CNVNormalize import CNVNormalize

class MpiWithinArrayBlockQuantileNormalize(MPIwrapper):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
			('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
			('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
			('dbname', 1, ): ['stock_250k', 'd', 1, '', ],\
			('db_user', 1, ): [None, 'u', 1, 'database username', ],\
			('db_passwd', 1, ):[None, 'p', 1, 'database password', ],\
			('array_id_ls', 0, ): [None, 'a', 1, 'comma or dash-separated array id list, like 61-70,81. Not specifying this means all arrays.'],\
			('call_method_id', 0, int):[0, 'l', 1, 'Restrict arrays included in this call_method. Default is no such restriction.'],\
			('array_file_directory', 0, ):['', 'f', 1, 'The results directory. Default is None. use the one given by db.'],\
			("probes_blockData_picklef", 0, ): [os.path.expanduser('~/script/variation/data/CNV/WABQN_probes_blockData_picklef'), 'e', 1, \
											'File that contains the probes & block partition info'],\
			("output_fname", 1, ): [None, 'o', 1, 'Filename to store output matrix'],\
			('blockSize', 1, int):[100, 's', 1, 'The dimension size for one block.'],\
			('jumpStep', 0, int):[None, 'j', 1, 'default is 1/4th of blockSize'],\
			('x_range', 0, ):['', 'x', 1, 'restrict probes within this x-range on the affy chip'],\
			('y_range', 0, ):['', 'y', 1, 'restrict probes within this y-range on the affy chip'],\
			('minNoOfProbesPerBlock', 0, int):[None, '', 1, 'default is 1/6th of blockSize*blockSize'],\
			('outlierBlockPercentile', 0, float):[1.0, '', 1, 'ignore this.'],\
			('message_size', 1, int):[1, 'q', 1, 'How many parameter combos one computing node should handle.'],\
			('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
			('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		if self.jumpStep is None:
			self.jumpStep = self.blockSize/4
		if self.minNoOfProbesPerBlock is None:
			self.minNoOfProbesPerBlock = self.blockSize*self.blockSize/6
		self.array_id_ls = getListOutOfStr(self.array_id_ls, data_type=int)
		self.x_range = getListOutOfStr(self.x_range, data_type=int)
		self.y_range = getListOutOfStr(self.y_range, data_type=int)
		
		self.communicator = MPI.world.duplicate()
		MPIwrapper.__init__(self, self.communicator, debug=self.debug, report=self.report)
		
	def generate_params(cls, db_250k, array_id_ls=[], call_method_id=None, array_file_directory=None):
		"""
		2010-5-24
			get the arrays
		"""
		array_info_table = Stock_250kDB.ArrayInfo.table.name
		sql_query = DB_250k2Array.generateSQLQueryToGetArrays(array_info_table, array_id_ls=array_id_ls, \
												call_method_id=call_method_id)
		rows = db_250k.metadata.bind.execute(sql_query)
		for row in rows:
			if array_file_directory:
				filename = os.path.join(array_file_directory, os.path.basename(row.filename))
			else:
				filename = row.filename
			yield (row.array_id, row.maternal_ecotype_id, filename)
	
	def getCommonData(self, db_250k, blockSize, jumpStep, x_range, y_range, minNoOfProbesPerBlock, \
						array_file_directory=None, probeType=2):
		"""
		2010-5-24
		"""
		sys.stderr.write("Getting array width ...")
		array_1 = Stock_250kDB.ArrayInfo.get(1)
		if array_file_directory:
			filename = os.path.join(array_file_directory, os.path.basename(array_1.filename))
		else:
			filename = array_1.filename
		
		returnData = DB_250k2Array.getArrayWidth(filename)
		array_width = returnData.array_width
		sys.stderr.write("Done.\n")
		
		blockBottomLeftXY_ls = self.generateBlocksWithinOneArray(array_width, blockSize=blockSize, jumpStep=jumpStep,
									x_range=x_range, y_range=y_range)
		
		probes, xy_ls, chr_pos_ls, probe_id_ls = DB_250k2Array.get_probes(db_250k.metadata.bind, \
									Stock_250kDB.Probes.table.name, snps=None, run_type=probeType, \
									x_range=x_range,\
									y_range=y_range, constructChrPos2index=False, constructXY2index=True,\
									need_xy_ls=True, need_chr_pos_ls=False, need_probe_id_ls=False)
		
		sys.stderr.write("Generating xy indices for each block ...")
		no_of_blocks = len(blockBottomLeftXY_ls)
		xy_set = set(xy_ls)
		# construct a list which stores a list of (x,y) in each block
		new_xy_set = set()	# some probes will be excluded if they are situated in blocks with too few probes
		blockDataCodedIndex_ls = []
		for k in range(no_of_blocks):
			bottomLeftXY = blockBottomLeftXY_ls[k]
			start_x, start_y = bottomLeftXY
			blockDataCodedIndex = []
			for x in xrange(start_x, start_x+blockSize):
				for y in xrange(start_y, start_y+blockSize):
					xy_pos = (x, y)
					if xy_pos in xy_set:	#make sure it's only the probes that we want
						blockDataCodedIndex.append(xy_pos)
			if len(blockDataCodedIndex)>=minNoOfProbesPerBlock:	# skip blocks that are too small
				blockDataCodedIndex_ls.append(blockDataCodedIndex)
				new_xy_set |= set(blockDataCodedIndex)
		sys.stderr.write("%s probes in %s blocks. Done.\n"%(len(new_xy_set), len(blockDataCodedIndex_ls)))
		
		sys.stderr.write("Generating new probes ...")
		new_probe_ls = [(probes.getProbeGivenXY(xy[0], xy[1]))for xy in new_xy_set]
		del probes
		new_probe_ls.sort()		# in chr,pos order
		xy2index = {}
		chr_pos_ls = []
		probe_id_ls = []
		for probe in new_probe_ls:
			xy = (probe.xpos, probe.ypos)
			xy2index[xy] = len(xy2index)
			chr_pos_ls.append((probe.chr, probe.pos))
			probe_id_ls.append(probe.probe_id)
		del new_probe_ls
		sys.stderr.write("Done.\n")
		commonData = PassingData(xy2index=xy2index, blockDataCodedIndex_ls=blockDataCodedIndex_ls, \
								array_width=array_width, probe_id_ls = probe_id_ls,\
								chr_pos_ls = chr_pos_ls)
		return commonData
	
	def prepareCommonData(self, db_250k, blockSize, jumpStep, x_range, y_range, minNoOfProbesPerBlock, \
						array_file_directory=None, probeType=2, probes_blockData_picklef=None):
		"""
		2010-5-24
		"""
		sys.stderr.write("Preparing commonData ...")
		if probes_blockData_picklef:
			if os.path.isfile(probes_blockData_picklef):	#if this file is already there, suggest to un-pickle it.
				picklef = open(probes_blockData_picklef)
				commonData = cPickle.load(picklef)
				del picklef
			else:	#if the file doesn't exist, but the filename is given, pickle commonData into it
				commonData = self.getCommonData(db_250k, blockSize, jumpStep, x_range, y_range, minNoOfProbesPerBlock, \
						array_file_directory=array_file_directory, probeType=probeType)
				picklef = open(probes_blockData_picklef, 'w')
				cPickle.dump(commonData, picklef, -1)
				picklef.close()
		else:
			commonData = self.getCommonData(db_250k, blockSize, jumpStep, x_range, y_range, minNoOfProbesPerBlock, \
						array_file_directory=array_file_directory, probeType=probeType)
		sys.stderr.write("Done.\n")
		return commonData
	
	def computing_node_handler(self, communicator, data, param_obj):
		"""
		2010-6-3
			the data received is a list of jobs to finish.
			modify it to handle >1 arrays at a time.
		2010-5-24
		
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		result_ls = []
		blockDataCodedIndex_ls = param_obj.blockDataCodedIndex_ls
		for one_data in data:
			array_id, ecotype_id, array_filename = one_data[:3]
			intensityBlockData = self.getIntensityInBlocks(array_filename, blockDataCodedIndex_ls, \
												array_width=param_obj.array_width, \
												drawHist=False, outlierBlockPercentile=1.0)
			qnorm_data_matrix = numpy.zeros([len(param_obj.xy2index), 1], numpy.float32)
			qnorm_data_matrix[:,:] = numpy.nan
			self.quantileNormalizeAllBlocks(intensityBlockData.new_blockData_ls, \
								intensityBlockData.new_blockDataCodedIndex_ls, param_obj.xy2index, \
								qnorm_data_matrix, array_index=0, returnBlockMatrix=False)
			result_ls.append((array_id, ecotype_id, qnorm_data_matrix.flatten().tolist()))
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return result_ls
	
	def writeHeader(self, writer, probe_id_ls, chr_pos_ls):
		"""
		2010-5-27
			call CNVNormalize.writeArrayXProbeHeader() to do the job.
		2010-5-24
		"""
		CNVNormalize.writeArrayXProbeHeader(writer, probe_id_ls, chr_pos_ls)
		
	def output_node_handler(self, communicator, output_param_obj, data):
		"""
		2010-5-24
		"""
		for array_id, ecotype_id, intensity_ls in cPickle.loads(data):
			writer = output_param_obj.writer
			writer.writerow([array_id, ecotype_id] + intensity_ls)
			sys.stderr.write("%s intensity data output.\n"%len(intensity_ls))
	
	@classmethod
	def extendLsByResampling(cls, data_ls, max_no_of_points):
		"""
		2010-5-19
			
		"""
		import random
		# now to make up data in the block if its number of data points is < max_no_of_points_in_one_block
		no_of_points = len(data_ls)
		if no_of_points<max_no_of_points:
			no_of_points_to_add = max_no_of_points - no_of_points
			no_of_copies = int(no_of_points_to_add/no_of_points)
			original_block_data = data_ls[:]	#[:] is to avoid reference. deep copy.
			data_ls.extend(original_block_data*no_of_copies)
			no_of_points_to_sample = no_of_points_to_add%no_of_points
			if no_of_points_to_sample>0:	# do a sampling
				data_ls.extend(random.sample(original_block_data, no_of_points_to_sample))
		return data_ls
	
	@classmethod
	def quantileNormalizeAllBlocks(cls, blockData_ls, blockDataCodedIndex_ls, xy2index, qnorm_data_matrix, array_index=None,
								returnBlockMatrix=False,):
		"""
		2010-5-19
		"""
		import numpy
		from CNVNormalize import CNVNormalize
		data_matrix = numpy.array(blockData_ls, numpy.float32).transpose()	#tranpose() is a must.
		"""
		quantile_normalize() takes each row as a probe. each column is an array/block.
		"""
		CNVNormalize.quantile_normalize(data_matrix)
		
		# aggregate data points of one probe
		codedIndex2QNormDataLs = {}
		no_of_blocks = len(blockData_ls)
		for i in range(no_of_blocks):
			blockDataCodedIndex = blockDataCodedIndex_ls[i]
			no_of_probes_in_block = len(blockDataCodedIndex)
			for j in range(no_of_probes_in_block):
				xy_pos = blockDataCodedIndex[j]
				if xy_pos not in codedIndex2QNormDataLs:
					codedIndex2QNormDataLs[xy_pos] = []
				intensity = data_matrix[j,i]
				codedIndex2QNormDataLs[xy_pos].append(intensity)
		# get the median of all qnormed intensity for one probe
		#codedIndex2Median = {}
		for xy_pos, data_ls in codedIndex2QNormDataLs.iteritems():
			row_index = xy2index.get(xy_pos)
			#codedIndex2Median[xy_pos] = numpy.median(data_ls)
			qnorm_data_matrix[row_index,array_index] = numpy.median(data_ls)
		if returnBlockMatrix:
			return data_matrix
		else:
			del data_matrix
	
	
	@classmethod
	def generateBlocksWithinOneArray(cls,array_width, blockSize=50, jumpStep=10,
									x_range=None, y_range=None):
		"""
		2010-5-17
			split out of subArrayQuantileNormalize()
		"""
		sys.stderr.write("Generating block outlines...")
		blockBottomLeftXY_ls = []
		import math
		if x_range:
			start_x = x_range[0]
			array_x_max = min(array_width, x_range[1]+1)
			array_x_span = x_range[1]-x_range[0] + 1
		else:
			start_x = 0
			array_x_max = array_width
			array_x_span = array_width
		if y_range:
			_start_y = y_range[0]
			array_y_max = min(array_width, y_range[1]+1)
			array_y_span = y_range[1]-y_range[0] + 1
		else:
			_start_y = 0
			array_y_max = array_width
			array_y_span = array_width
		no_of_blocks_in_x = int(math.ceil((array_x_span-blockSize)/float(jumpStep)))+1
		no_of_blocks_in_y = int(math.ceil((array_y_span-blockSize)/float(jumpStep)))+1
		for i in range(no_of_blocks_in_x):
			if start_x+blockSize>array_x_max:	#reduce the x to maintain the same block size
				start_x = array_x_max-blockSize
			start_y = _start_y
			for j in range(no_of_blocks_in_y):
				if start_y+blockSize>array_y_max:	#reduce the y to maintain the same block size
					start_y = array_y_max-blockSize
				blockBottomLeftXY_ls.append((start_x,start_y))
				start_y += jumpStep
			start_x += jumpStep
		sys.stderr.write("%s blocks. Done.\n"%(len(blockBottomLeftXY_ls)))
		return blockBottomLeftXY_ls
	
	
	@classmethod
	def getIntensityInBlocks(cls, array_filename, blockDataCodedIndex_ls, \
			array_width=None, drawHist=False, outlierBlockPercentile=1.0):
		"""
		2010-6-3
			import rpy in this function locally
		2010-5-24
		"""
		sys.stderr.write("Getting intensity data in blocks ...\n")
		blockData_ls = []
		import rpy
		rpy.r.library('affy')
		array = rpy.r.read_affybatch(filenames=array_filename)
		intensity_array = rpy.r.intensity(array)
		blockDataSigma_ls = []
		no_of_blocks = len(blockDataCodedIndex_ls)	# no_of_blocks has changed
		max_no_of_points_in_one_block = max([len(row) for row in blockDataCodedIndex_ls])
		for i in range(no_of_blocks):
			blockDataCodedIndex = blockDataCodedIndex_ls[i]
			blockData_ls.append([])
			for xy_pos in blockDataCodedIndex:
				xpos, ypos = xy_pos
				intensity_array_index = array_width*(array_width - xpos - 1) + ypos
				intensity = math.log10(intensity_array[intensity_array_index][0])
				blockData_ls[i].append(intensity)
			# now to make up data in the block if its number of points is < max_no_of_points_in_one_block
			blockData_ls[i] = cls.extendLsByResampling(blockData_ls[i], max_no_of_points_in_one_block)
			
			if outlierBlockPercentile<1.0 and outlierBlockPercentile>0.0:
				sigmaEstimate = numpy.std(blockData_ls[i])/numpy.sqrt(max_no_of_points_in_one_block-1)
				meanEstimate = numpy.mean(blockData_ls[i])
				blockDataSigma_ls.append(sigmaEstimate/meanEstimate)
		
		if outlierBlockPercentile<1.0 and outlierBlockPercentile>0.0:
			# only take blocks whose sigma is in bottom outlierBlockPercentile, i.e. 90%
			blockDataSigma_argsort_ls = numpy.argsort(blockDataSigma_ls)
			maxRank =  int(no_of_blocks*outlierBlockPercentile)
			new_blockData_ls = []
			new_blockDataCodedIndex_ls = []
			new_blockDataSigma_ls = []
			for i in xrange(maxRank):
				old_block_index = blockDataSigma_argsort_ls[i]
				new_blockDataSigma_ls.append(blockDataSigma_ls[old_block_index])
				new_blockData_ls.append(blockData_ls[old_block_index])
				new_blockDataCodedIndex_ls.append(blockDataCodedIndex_ls[old_block_index])
		else:
			new_blockData_ls = blockData_ls
			new_blockDataCodedIndex_ls = blockDataCodedIndex_ls
			new_blockDataSigma_ls = blockDataSigma_ls
		new_no_of_blocks = len(new_blockData_ls)
		sys.stderr.write(" %s blocks -> %s blocks.\n"%(no_of_blocks, new_no_of_blocks))
		return PassingData(new_blockData_ls=new_blockData_ls, new_blockDataCodedIndex_ls=new_blockDataCodedIndex_ls)
		
	
	def run(self):
		"""
		2009-1-22
		"""
		if self.debug:
			#for one-node testing purpose
			import pdb
			pdb.set_trace()
		
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		free_computing_node_set = Set(free_computing_nodes)
		output_node_rank = self.communicator.size-1
		
		if node_rank == 0:
			db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				  				password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
			db.setup(create_tables=False)
			param_ls = self.generate_params(db, array_id_ls=self.array_id_ls, call_method_id=self.call_method_id,\
										array_file_directory=self.array_file_directory)
			
			commonData = self.prepareCommonData(db, self.blockSize, self.jumpStep, \
						self.x_range, self.y_range, self.minNoOfProbesPerBlock, \
						array_file_directory=self.array_file_directory, probeType=2, \
						probes_blockData_picklef=self.probes_blockData_picklef)
			if self.communicator.size==1:	# single-node serial run
				blockDataCodedIndex_ls = commonData.blockDataCodedIndex_ls
				
				output_dir = os.path.split(self.output_fname)[0]
				if not os.path.isdir(output_dir):
					os.makedirs(output_dir)
				writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
				probe_id_ls = commonData.probe_id_ls
				chr_pos_ls = commonData.chr_pos_ls
				self.writeHeader(writer, probe_id_ls, chr_pos_ls)
				for array_id, ecotype_id, array_filename in param_ls:
					intensityBlockData = self.getIntensityInBlocks(array_filename, blockDataCodedIndex_ls, \
															array_width=commonData.array_width, \
															drawHist=False, outlierBlockPercentile=1.0)
					qnorm_data_matrix = numpy.zeros([len(commonData.xy2index), 1], numpy.float32)
					qnorm_data_matrix[:,:] = numpy.nan
					self.quantileNormalizeAllBlocks(intensityBlockData.new_blockData_ls, \
											intensityBlockData.new_blockDataCodedIndex_ls, commonData.xy2index, \
											qnorm_data_matrix, array_index=0, returnBlockMatrix=False)
					writer.writerow([array_id, ecotype_id] + qnorm_data_matrix.flatten().tolist())
				sys.exit(0)
			
			commonData_pickle = cPickle.dumps(commonData, protocol=-1)
			sys.stderr.write("Passing data to output node %s from %s ... "%(output_node_rank, node_rank,))
			self.communicator.send(commonData_pickle, output_node_rank, 0)
			sys.stderr.write(".\n")
			
			for node in free_computing_nodes:	#send it to the computing_node
				sys.stderr.write("passing initial data to nodes from %s to %s ... "%(node_rank, node))
				self.communicator.send(commonData_pickle, node, 0)
				sys.stderr.write(".\n")
			if len(commonData.blockDataCodedIndex_ls)==0:
				sys.stderr.write("Not a single block is formed. Exit!")
				sys.exit(0)
			del commonData, commonData_pickle
			
		elif node_rank in free_computing_node_set:
			data, source, tag = self.communicator.receiveString(0, 0)
			commonData =  cPickle.loads(data)
			if len(commonData.blockDataCodedIndex_ls)==0:
				sys.stderr.write("Not a single block is formed. Exit!")
				sys.exit(0)
			xy2index = commonData.xy2index
			blockDataCodedIndex_ls = commonData.blockDataCodedIndex_ls
			array_width = commonData.array_width
			del data, commonData
		else:
			data, source, tag = self.communicator.receiveString(0, 0)
			commonData = cPickle.loads(data)
			if len(commonData.blockDataCodedIndex_ls)==0:
				sys.stderr.write("Not a single block is formed. Exit!")
				sys.exit(0)
			probe_id_ls = commonData.probe_id_ls
			chr_pos_ls = commonData.chr_pos_ls
			del data, commonData
			
		self.synchronize()
		if node_rank == 0:
			param_obj = PassingData(param_ls=param_ls, output_node_rank=output_node_rank, report=self.report, counter=0)
			self.inputNode(param_obj, free_computing_nodes, param_generator = param_ls, message_size=self.message_size)
		elif node_rank in free_computing_node_set:
			computing_parameter_obj = PassingData(xy2index=xy2index, blockDataCodedIndex_ls=blockDataCodedIndex_ls, \
												array_width=array_width)
			self.computing_node(computing_parameter_obj, self.computing_node_handler)
		else:
			output_dir = os.path.split(self.output_fname)[0]
			if not os.path.isdir(output_dir):
				os.makedirs(output_dir)
			writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
			self.writeHeader(writer, probe_id_ls, chr_pos_ls)
			
			output_param_obj = PassingData(writer=writer)
			self.output_node(free_computing_nodes, output_param_obj, self.output_node_handler)
			del writer
		self.synchronize()	#to avoid some node early exits
		

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiWithinArrayBlockQuantileNormalize
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
