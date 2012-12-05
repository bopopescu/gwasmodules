#!/usr/bin/env mpipython
"""

Examples:
	#run it on hpc-cmb cluster 
	mpiexec ~/script/variation/src/MpiTargetVsRefBlockLTS.py ... 
	
	#test parallel run on desktop
	mpirun -np 6 -machinefile  /tmp/hostfile /usr/bin/mpipython ~/script/variation/src/MpiTargetVsRefBlockLTS.py 
	-i ~/script/variation/data/CNV/call_48_WABQN_b200_j200.tsv -o ~/script/variation/data/CNV/call_48_WABQN_b200_j100_lts.tsv 
	-u yh -p xxx -e ~/script/variation/data/CNV/LTS_probes_blockData_b200_j100_picklef -s 200 -j 100
	
	#to test and debug. use array 3 as reference. (serial, rather than parallel)
	python ~/script/variation/src/MpiTargetVsRefBlockLTS.py -i /tmp/array_3_4.tsv -a 3 -o /tmp/array_3_4_lts.tsv -y 1200,1612 -x 0,200 -u yh -b
	
	#run all arrays from call method 48 on hpc-cmb
	mpiexec ~/script/variation/src/MpiTargetVsRefBlockLTS.py \
		-i /Network/Data/250k/tmp-yh/CNV/call_48_WABQN.tsv -o /Network/Data/250k/tmp-yh/CNV/call_48_WABQN_lts.tsv -s 200 -j 100
		-u yh -f ~/panfs/db/raw_data/ -p ***
	
Description:
	2010-5-25
		This program does least-trimmed-square between reference and target to remove reference & probe-specific effect.
		It does so block-by-block in order to account for uneven-spatial distribution and the fact that linear approximation
		between reference and target is more likely to be true in local than in global.
		
		The model is y = a+bx, with a trimming fraction set at 40%. y is target intensity, x is reference.
			After fitting, y-(a+bx) is the new intensity for the target array.
		In the case of multiple referenes, tt takes median. 
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
from pymodule.algorithm import ltsFit
from sets import Set
from Scientific import MPI
from pymodule.MPIwrapper import MPIwrapper
import Stock_250kDB, numpy, random, cPickle
from DB_250k2Array import DB_250k2Array, probes_class
from CNVNormalize import CNVNormalize
from MpiWithinArrayBlockQuantileNormalize import MpiWithinArrayBlockQuantileNormalize

class MpiTargetVsRefBlockLTS(MpiWithinArrayBlockQuantileNormalize):
	__doc__ = __doc__
	option_default_dict = MpiWithinArrayBlockQuantileNormalize.option_default_dict.copy()
	option_default_dict.update({('input_fname', 1,):['', 'i', 1, 'output of MpiWithinArrayBlockQuantileNormalize.py', ]})
	option_default_dict.update({("probes_blockData_picklef", 0, ): \
							[os.path.expanduser('~/script/variation/data/CNV/LTS_probes_blockData_picklef'), 'e', 1, \
							'File that contains the probes & block partition info'],})
	option_default_dict.update({('array_id_ls', 0, ): ['1,2,43,139,145,1339', 'a', 1, 'comma or dash-separated reference array id list,'],})
	
	#option_default_dict.pop(('drivername', 1,))
	#option_default_dict.pop(('hostname', 1, ))
	#option_default_dict.pop(('dbname', 1, ))
	#option_default_dict.pop(('schema', 0, ))
	#option_default_dict.pop(('db_user', 1, ))
	#option_default_dict.pop(('db_passwd', 1, ))
	option_default_dict.pop(('call_method_id', 0, int))
	
	def __init__(self, **keywords):
		MpiWithinArrayBlockQuantileNormalize.__init__(self, **keywords)
		self.ref_array_id_set = set(self.array_id_ls)
		
	def generate_params(cls, reader=None, blockDataCodedIndex_ls=None, ref_array_id_set=None):
		"""
		2010-5-24
			read the input filename
		"""
		for row in reader:
			array_id = int(row[0])
			if array_id not in ref_array_id_set:	#skip the reference arrays
				ecotype_id = row[1]
				intensity_ls = row[2:]
				intensity_array = numpy.array(map(float, intensity_ls))
				for i in xrange(len(blockDataCodedIndex_ls)):
					blockDataCodedIndex = blockDataCodedIndex_ls[i]
					blockIntensity_ls = intensity_array.take(blockDataCodedIndex)
					yield (array_id, ecotype_id, i, blockIntensity_ls)
	
	def getCommonData(self, db_250k, blockSize, jumpStep, x_range, y_range, minNoOfProbesPerBlock, \
						array_file_directory=None, probe_id_ls=None, chr_pos_ls=None, probeType=2):
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
		
		probe_id2col_index = {}
		for i in xrange(len(probe_id_ls)):
			probe_id = probe_id_ls[i]
			probe_id2col_index[probe_id] = i
		
		probes = DB_250k2Array.get_probes(db_250k.metadata.bind, \
			Stock_250kDB.Probes.table.name, snps=None, run_type=probeType, \
			x_range=x_range,\
			y_range=y_range, constructChrPos2index=False, constructXY2index=True,\
			need_xy_ls=False, need_chr_pos_ls=False, need_probe_id_ls=False)[0]
		
		sys.stderr.write("Generating index data for each block ...")
		no_of_blocks = len(blockBottomLeftXY_ls)
		# construct a list which stores a list of col_index in each block
		included_col_index_set = set()
		blockDataCodedIndex_ls = []
		for k in range(no_of_blocks):
			bottomLeftXY = blockBottomLeftXY_ls[k]
			start_x, start_y = bottomLeftXY
			blockDataCodedIndex = []
			for x in xrange(start_x, start_x+blockSize):
				for y in xrange(start_y, start_y+blockSize):
					xy = (x, y)
					probe = probes.getProbeGivenXY(xy[0], xy[1])
					if probe and probe.probe_id in probe_id2col_index:
						col_index = probe_id2col_index.get(probe.probe_id)
						blockDataCodedIndex.append(col_index)
			if len(blockDataCodedIndex)>=minNoOfProbesPerBlock:	# skip blocks that are too small
				blockDataCodedIndex_ls.append(blockDataCodedIndex)
				included_col_index_set |= set(blockDataCodedIndex)
		sys.stderr.write("%s probes in %s blocks. Done.\n"%(len(included_col_index_set), len(blockDataCodedIndex_ls)))
		del probes
		
		commonData = PassingData(blockDataCodedIndex_ls=blockDataCodedIndex_ls, \
								array_width=array_width, probe_id_ls=probe_id_ls, chr_pos_ls=chr_pos_ls)
		return commonData
	
	def prepareCommonData(self, db_250k, blockSize, jumpStep, x_range, y_range, minNoOfProbesPerBlock, \
						array_file_directory=None, probe_id_ls=None, chr_pos_ls=None, probeType=2, \
						probes_blockData_picklef=None):
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
						array_file_directory=array_file_directory, probe_id_ls=probe_id_ls, chr_pos_ls=chr_pos_ls,\
						probeType=probeType)
				picklef = open(probes_blockData_picklef, 'w')
				cPickle.dump(commonData, picklef, -1)
				picklef.close()
		else:
			commonData = self.getCommonData(db_250k, blockSize, jumpStep, x_range, y_range, minNoOfProbesPerBlock, \
						array_file_directory=array_file_directory, probe_id_ls=probe_id_ls, chr_pos_ls=chr_pos_ls,\
						probeType=probeType)
		sys.stderr.write("Done.\n")
		return commonData
	
	def readInRefArrayData(self, input_fname, ref_array_id_set=None):
		"""
		2010-5-25
		"""
		
		sys.stderr.write("Getting data matrix for reference arrays.\n")
		reader = csv.reader(open(input_fname), delimiter=figureOutDelimiter(input_fname))
		for i in xrange(3):	# skip first 3 rows
			reader.next()
		data_matrix = []
		for row in reader:
			array_id = int(row[0])
			if array_id in ref_array_id_set:
				data_matrix.append(map(float, row[2:]))
		del reader
		data_matrix = numpy.array(data_matrix)
		sys.stderr.write("%s arrays, %s probes. Done.\n"%(data_matrix.shape[0], data_matrix.shape[1]))
		return data_matrix
	
	def ltsOneBlockAgainstAllRef(self, array_id, ecotype_id, blockIndex, blockIntensity_ls, refDataMatrix, blockDataCodedIndex_ls):
		"""
		2010-5-28
			reverse the order of reference and target. ref = a + b*target.. then new-target = a+b*target-ref
				previously. target = a + b*ref. then new-target = target - (a + b*ref)
			The order-reversal is due to Lei Li, who says the directionality matters (It might not converge in the other direction).
		2010-5-25
		"""
		no_of_refs = refDataMatrix.shape[0]
		blockDataCodedIndex = blockDataCodedIndex_ls[blockIndex]
		after_lts_data_matrix = numpy.zeros([no_of_refs, len(blockDataCodedIndex)], numpy.float32)
		after_lts_data_matrix[:,:] = numpy.nan
		
		for i in xrange(no_of_refs):
			refBlockIntensity_ls = refDataMatrix[i].take(blockDataCodedIndex)
			ltsResult = ltsFit(blockIntensity_ls, refBlockIntensity_ls)
			predictedTargetIntensity_ls = ltsResult.fit_y_ls
			after_lts_data_matrix[i,:] = predictedTargetIntensity_ls-refBlockIntensity_ls
		
		return (array_id, ecotype_id, blockIndex, numpy.median(after_lts_data_matrix, axis=0))
	
	def computing_node_handler(self, communicator, data, param_obj):
		"""
		2010-6-3
			the data received is a list of jobs to finish.
			modify it to handle >1 jobs at a time.
		2010-5-24
		
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		blockDataCodedIndex_ls = param_obj.blockDataCodedIndex_ls
		result_ls = []
		for one_data in data:
			array_id, ecotype_id, blockIndex, blockIntensity_ls = one_data[:4]
			result_ls.append(self.ltsOneBlockAgainstAllRef(array_id, ecotype_id, blockIndex, blockIntensity_ls,\
										param_obj.refDataMatrix, blockDataCodedIndex_ls))
		if self.report:
			sys.stderr.write("node %s: %s.\n"%(node_rank, repr(result_ls[:3])))
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		#sys.stderr.write("node %s: %s data points. example: %s"%(node_rank, len(result_ls[0][3]), repr(result_ls[0][3][:10])))
		return result_ls
	
	def handleComputingOutputData(self, result_ls, blockDataCodedIndex_ls, array_id2no_of_blocks_returned, \
								array_id2col_index2intensity_ls, writer):
		"""
		2010-5-24
		"""
		for array_id, ecotype_id, blockIndex, blockIntensity_ls in result_ls:
			blockDataCodedIndex = blockDataCodedIndex_ls[blockIndex]
			if array_id not in array_id2no_of_blocks_returned:
				array_id2no_of_blocks_returned[array_id] = 0
				array_id2col_index2intensity_ls[array_id] = {}
			array_id2no_of_blocks_returned[array_id] += 1
			
			for i in xrange(len(blockIntensity_ls)):
				col_index = blockDataCodedIndex[i]
				if col_index not in array_id2col_index2intensity_ls[array_id]:
					array_id2col_index2intensity_ls[array_id][col_index] = []
				array_id2col_index2intensity_ls[array_id][col_index].append(blockIntensity_ls[i])
			
			if array_id2no_of_blocks_returned[array_id] == len(blockDataCodedIndex_ls):
				sys.stderr.write("Outputting array %s intensity ...\n"%array_id)
				no_of_probes = len(array_id2col_index2intensity_ls[array_id])
				intensity_ls = ['NA']*no_of_probes
				for col_index, data_ls in array_id2col_index2intensity_ls[array_id].iteritems():
					intensity_ls[col_index] = numpy.median(data_ls)
				
				writer.writerow([array_id, ecotype_id] + intensity_ls)
				del array_id2no_of_blocks_returned[array_id]
				del array_id2col_index2intensity_ls[array_id]
				sys.stderr.write("array %s intensity output Done.\n"%array_id)
	
	def output_node_handler(self, communicator, output_param_obj, data):
		"""
		2010-5-24
		"""
		writer = output_param_obj.writer
		blockDataCodedIndex_ls = output_param_obj.blockDataCodedIndex_ls
		array_id2no_of_blocks_returned = output_param_obj.array_id2no_of_blocks_returned
		array_id2col_index2intensity_ls = output_param_obj.array_id2col_index2intensity_ls
		
		result_ls = cPickle.loads(data)
		self.handleComputingOutputData(result_ls, blockDataCodedIndex_ls, array_id2no_of_blocks_returned, \
							array_id2col_index2intensity_ls, writer)
	
	def run(self):
		"""
		2010-5-25
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
		
		# 2010-5-25 to hold final data
		array_id2no_of_blocks_returned = {}
		array_id2col_index2intensity_ls = {}
		
		if node_rank == 0:
			db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				  				password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
			db.setup(create_tables=False)
			
			reader = csv.reader(open(self.input_fname), delimiter=figureOutDelimiter(self.input_fname))
			
			probe_id_ls = reader.next()[2:]
			probe_id_ls = map(int, probe_id_ls)
			chr_ls = reader.next()[2:]
			pos_ls = reader.next()[2:]
			chr_pos_ls = zip(chr_ls, pos_ls)
			
			commonData = self.prepareCommonData(db, self.blockSize, self.jumpStep, \
						self.x_range, self.y_range, self.minNoOfProbesPerBlock, \
						array_file_directory=self.array_file_directory, probe_id_ls=probe_id_ls,\
						chr_pos_ls=chr_pos_ls, probeType=2, \
						probes_blockData_picklef=self.probes_blockData_picklef)
			param_ls = self.generate_params(reader, blockDataCodedIndex_ls=commonData.blockDataCodedIndex_ls, \
							ref_array_id_set=self.ref_array_id_set)	#must be behind prepareCommonData()
			refDataMatrix = self.readInRefArrayData(self.input_fname, ref_array_id_set=self.ref_array_id_set)
			
			if self.communicator.size==1:	# single-node serial run
				blockDataCodedIndex_ls = commonData.blockDataCodedIndex_ls
				
				output_dir = os.path.split(self.output_fname)[0]
				if not os.path.isdir(output_dir):
					os.makedirs(output_dir)
				writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
				probe_id_ls = commonData.probe_id_ls
				chr_pos_ls = commonData.chr_pos_ls
				self.writeHeader(writer, probe_id_ls, chr_pos_ls)
				for array_id, ecotype_id, blockIndex, blockIntensity_ls in param_ls:
					result_ls = self.ltsOneBlockAgainstAllRef(array_id, ecotype_id, blockIndex, blockIntensity_ls,\
												refDataMatrix, blockDataCodedIndex_ls)
					self.handleComputingOutputData(result_ls, blockDataCodedIndex_ls, array_id2no_of_blocks_returned, \
							array_id2col_index2intensity_ls, writer)
				sys.exit(0)
			
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
			
		elif node_rank in free_computing_node_set:
			data, source, tag = self.communicator.receiveString(0, 0)
			commonData =  cPickle.loads(data)
			if len(commonData.blockDataCodedIndex_ls)==0:
				sys.stderr.write("Not a single block is formed. Exit!")
				sys.exit(0)
			blockDataCodedIndex_ls = commonData.blockDataCodedIndex_ls
			del data, commonData
			
			data, source, tag,= self.communicator.receiveString(0, 0)
			refDataMatrix =  cPickle.loads(data)
			del data
		else:
			data, source, tag = self.communicator.receiveString(0, 0)
			commonData = cPickle.loads(data)
			if len(commonData.blockDataCodedIndex_ls)==0:
				sys.stderr.write("Not a single block is formed. Exit!")
				sys.exit(0)
			probe_id_ls = commonData.probe_id_ls
			chr_pos_ls = commonData.chr_pos_ls
			blockDataCodedIndex_ls = commonData.blockDataCodedIndex_ls
			del data, commonData
			
		self.synchronize()
		if node_rank == 0:
			param_obj = PassingData(param_ls=param_ls, output_node_rank=output_node_rank, report=self.report, counter=0)
			self.inputNode(param_obj, free_computing_nodes, param_generator = param_ls, message_size=self.message_size)
		elif node_rank in free_computing_node_set:
			computing_parameter_obj = PassingData(refDataMatrix=refDataMatrix, \
												blockDataCodedIndex_ls=blockDataCodedIndex_ls)
			self.computing_node(computing_parameter_obj, self.computing_node_handler, output_node_rank=output_node_rank)
		else:
			output_dir = os.path.split(self.output_fname)[0]
			if not os.path.isdir(output_dir):
				os.makedirs(output_dir)
			writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
			self.writeHeader(writer, probe_id_ls, chr_pos_ls)
			
			output_param_obj = PassingData(writer=writer, array_id2no_of_blocks_returned=array_id2no_of_blocks_returned,\
										array_id2col_index2intensity_ls=array_id2col_index2intensity_ls,\
										blockDataCodedIndex_ls=blockDataCodedIndex_ls)
			self.output_node(free_computing_nodes, output_param_obj, self.output_node_handler)
			del writer
		self.synchronize()	#to avoid some node early exits

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiTargetVsRefBlockLTS
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
