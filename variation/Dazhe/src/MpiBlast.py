#!/usr/bin/env mpipython
"""
Examples:
	# parallel run on hpc-cmb
	mpiexec ~/script/variation/Dazhe/src/MpiBlast.py -o ~/script/variation/data/CNV/ler_blast.csv -m 25
	
	# run on a single node and output the matrix to output_fname. (add -b to turn into debug mode)
	python ~/script/variation/Dazhe/src/MpiBlast.py -f /usr/bin/blastall
		-o /tmp/ler_blast.csv -p ~/script/variation/data/CNV/CNV_probelist.csv  -r -m 25
	
	# 2010-6-21 blast probes (TAIR8 coordinates) against TAIR9 sequence 
	mpiexec ~/script/variation/Dazhe/src/MpiBlast.py -o ~/script/variation/data/CNV/stock_250k_probes_vs_TAIR9.tsv
	-p ~/script/variation/data/CNV/stock_250k_probes.csv -d ~/script/variation/data/CNV/TAIR9/tair9.fas -m 24
	
Description:
	The blast script that could run in either parallel MPI or serial fashion.
	The output has these columns:
	'Chromosome', 'Position', 'Probe_ID', 'Alignment_title', "query_start", "query_end", 'Number_matches', 'Alignment_start', 'Alignment_stop'
"""
import sys, os, math
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from Scientific import MPI
import cPickle, sys, csv, traceback
import vardata
from Bio.Blast import NCBIXML, NCBIStandalone
from pymodule.MPIwrapper import mpi_synchronize, MPIwrapper
from pymodule import PassingData

class ref_genome():
	def __init__(self, datadir="/home/cmb-01/dazhemen/ref_genome/", suffix=None, no_of_chromosomes=5):
		self.datadir = datadir
		self.suffix = suffix
		self.chr = [""]*no_of_chromosomes
		self.no_of_chromosomes = no_of_chromosomes
	
	def load_chr(self):
		for i in xrange(1, self.no_of_chromosomes + 1):
			self.__load_chr__(i)
	
	def __load_chr__(self, number):
		"""
		2010-4-15
		"""
		sys.stderr.write("Loading chromosome %s ... "%number)
		if self.suffix:
			chr_fname = "chr%s.%s"%(number, self.suffix)
		else:
			chr_fname = "chr%s"%(number)
		seqf = open(os.path.join(self.datadir, chr_fname))
		seqf.readline()
		seq = ""
		for line in seqf:
			if line[0]=='>':	#2010-4-15 ignore fasta file title
				continue
			seq += line.strip()
		self.chr[number-1] = seq
		sys.stderr.write("Done.\n")
	
	def readprobe(self,chr,pos):
		"""
		2010-6-21
			fix a bug here.
				argument pos is a variable starting from 1,
				while position in self.chr[chr-1] starts from 0.
		"""
		pos = pos-1	#2010-6-21
		return self.chr[chr-1][pos-12:pos+13]

def dazheMpiBlast():
	"""
	2010-4-14
		wrap all of dazhe's old code into this function
	"""
	blast_bin_path='/home/cmb-01/yuhuang/bin/blast/bin/blastall'
	database_fname='/home/cmb-01/dazhemen/Data/Ler/Cereon_Ath_Ler.fasta'
	
	probes = vardata.vardata()
	
	comm = MPI.world.duplicate()
	ppp = 10 # probes per processor
	accs = []
	
	genome = ref_genome()
	genome.load_chr()
	probes.readfromfile("/home/cmb-01/dazhemen/CNV/CNV_probelist.csv", format=2)
	ppp = len(probes.data)/(comm.size-1) + 1
	
	print "I'm %s of %s, finished loading\n"%(comm.rank, comm.size)
	
	if comm.rank == 0: #Master node, reads the file etc
		outf = open("/home/cmb-01/dazhemen/CNV/ler_raw.csv","w")
		outf.write("Chromosome,Position,Probe_ID,Alignment_title,Number_matches\n")
		for dest in range(1,comm.size):
			data, source, tag = comm.receiveString(dest, None)
			print "I'm 0, collected data from %s\n"%source
			partial_result_list = cPickle.loads(data)
			for result_entry in partial_result_list:
				outf.write(result_entry[0].replace("_",",")+","+",".join([str(a) for a in result_entry[1:]])+"\n")
	else:
		probe_index_start = (comm.rank-1)*ppp
		result_ls = []
		if probe_index_start + ppp >= len(probes.data):
			ppp = len(probes.data) - probe_index_start
		if ppp != 0:
			tmp_blast_infname = '/home/cmb-01/dazhemen/CNV/tmp_blast/'+str(comm.rank)
			inf = open(tmp_blast_infname, 'w')
			for i in xrange(probe_index_start, probe_index_start+ppp):
				inf.write(">%s_%s_%s\n"%(probes.data[i][0],probes.data[i][1],probes.data[i][2][0])) # write the probe id
				inf.write("%s\n"%genome.readprobe(probes.data[i][0],probes.data[i][1]))
			inf.close()
			print "I'm %s, finished with generating blast file from probes %s to %s\n"%(comm.rank,probe_index_start,probe_index_start+ppp)
			result_handle, error_info = NCBIStandalone.blastall(blast_bin_path, "blastn", database_fname, tmp_blast_infname, align_view=7)
			blast_records = NCBIXML.parse(result_handle)
			print "I'm %s, finished with blasting\n"%comm.rank
			while 1:
				try:
					blast_record = blast_records.next()
				except:
					"I'm %s, finished with all records\n"%comm.rank
					break 
				no_of_hits = min(1000, len(blast_record.alignments))
				for i in range(no_of_hits):
					alignment_title = blast_record.alignments[i].title
					for hsp in blast_record.alignments[i].hsps:
						result_entry = [blast_record.query, alignment_title, 0, 0] #[query name (probe id and pos) , alignment title , number of matches, pos in contig ]
						if hsp.identities >= 15:
							result_entry[2] = hsp.identities
							result_entry[3] = hsp.sbjct_start
							result_ls.append(result_entry)   
			print "I'm %s, finished with parsing blast result\n"%comm.rank
		result_data = cPickle.dumps(result_ls)
		comm.send(result_data,0,0)

class MpiBlast(MPIwrapper):
	__doc__ = __doc__
	option_default_dict = {('database_fname', 0, ): [os.path.expanduser('~/script/variation/data/CNV/LerContig/Cereon_Ath_Ler.fasta'), 'd', 1, ''],\
						('blast_bin_path', 1, ): [os.path.expanduser('~/bin/blast/bin/blastall'), 'f', 1, 'blast binary'],\
						('probe_fname', 1, ): ['/home/cmb-01/dazhemen/CNV/CNV_probelist.csv', 'p', 1, ''],\
						('ref_genome_dir', 1, ): [os.path.expanduser('~/script/variation/data/TAIR8/'), 'u', 1, 'directory that contains ref genome sequence files. filename as chr1, chr2, ...'],\
						('ref_genome_suffix', 1, ): ['fas', 'e', 1, 'ref genome filename suffix.'],\
						('tmp_dir', 1, ): ['/tmp/', 't', 1, 'temporary directory to store blast input'],\
						('output_fname', 1, ): [None, 'o', 1, ''],\
						('block_size', 1, int):[250, 's', 1, 'number of probes each computing node will blast each time.'],\
						('min_no_of_identities', 1, int): [18, 'm', 1, 'minimum number of identities', ],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2010-4-14
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def create_init_data(self, probe_fname, ref_genome_dir):
		"""
		2010-4-14
		"""
		init_data = PassingData()
		genome = ref_genome(datadir=ref_genome_dir, suffix=self.ref_genome_suffix)
		genome.load_chr()
		probes = vardata.vardata()
		probes.readfromfile(probe_fname, format=2)
		init_data.genome = genome
		init_data.probes = probes
		return init_data
	
	def input_handler(self, parameter_list, message_size, report=0):
		"""
		2010-4-14
			prepare probe_id and probe_seq for the computing node
		"""
		if report:
			sys.stderr.write("Fetching stuff...\n")
		param_d = parameter_list[0]
		probes = param_d.probes
		genome = param_d.genome
		no_of_probes = len(probes.data)
		
		if param_d.index<no_of_probes:
			index_stop = min(no_of_probes, param_d.index+message_size)
			data_to_return = self.get_probe_id_seq_ls(genome, probes, param_d.index, index_stop)
			param_d.index = index_stop
		else:
			data_to_return = None
		if report:
			sys.stderr.write("Fetching done.\n")
		return data_to_return
	
	def get_probe_id_seq_ls(self, genome, probes, probe_index_start, probe_index_stop):
		"""
		2010-4-15
			return a list of [probe_id, probe_seq] given probe_index_start, probe_index_stop
		"""
		data_to_return = []
		for i in xrange(probe_index_start, probe_index_stop):
			probe_id = "%s_%s_%s"%(probes.data[i][0], probes.data[i][1], probes.data[i][2][0]) # write the probe id
			probe_seq = genome.readprobe(probes.data[i][0], probes.data[i][1])	# chr, pos
			data_to_return.append([probe_id, probe_seq])
		return data_to_return
	
	def blastOneBatchProbes(self, probe_id_seq_ls, blast_bin_path, database_fname, \
						tmp_blast_infname, min_no_of_identities=15, node_rank=0):
		"""
		2010-4-14
		"""
		result_ls = []
		inf = open(tmp_blast_infname, 'w')
		for probe_id, probe_seq in probe_id_seq_ls:
			inf.write(">%s\n"%probe_id) # write the probe id
			inf.write("%s\n"%probe_seq)
		inf.close()
		if self.report:
			sys.stderr.write("I'm %s, finished generating blast file for %s probes.\n"%\
						(node_rank, len(probe_id_seq_ls)))
		result_handle, error_info = NCBIStandalone.blastall(blast_bin_path, "blastn", database_fname, tmp_blast_infname, align_view=7)
		#error_info = error_info.read()	#2010-4-14 this read() causes program to hang out forever. ???
		#if error_info:
		#	sys.stderr.write("%s"%error_info)
		blast_records = NCBIXML.parse(result_handle)
		
		if self.report:
			sys.stderr.write("I'm %s, finished blasting.\n"%node_rank)
		for blast_record in blast_records:
			no_of_hits = min(1000, len(blast_record.alignments))	# top 1000 or the number of available alignments
			for i in range(no_of_hits):
				alignment_title = blast_record.alignments[i].title
				for hsp in blast_record.alignments[i].hsps:
					if hsp.identities >= min_no_of_identities:
						result_entry = [blast_record.query, alignment_title, hsp.query_start, hsp.query_end, \
									hsp.identities, hsp.sbjct_start, hsp.sbjct_end,]
						#20104-25 hsp.strand is always (None, None), hsp.frame is either (1,1) or (1, -1) when the query's end < start
							#[query name (probe id and pos) , alignment title , number of matches, pos in contig ]
						result_ls.append(result_entry)
		if self.report:
			sys.stderr.write("I'm %s, finished with %s blasts, got %s returns.\n"%\
							(node_rank, len(probe_id_seq_ls), len(result_ls)))
		return result_ls
	
	def computing_node_handler(self, communicator, data, computing_parameter_obj):
		"""
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		probe_id_seq_ls = data
		init_data = computing_parameter_obj.init_data
		tmp_blast_infname = os.path.join(computing_parameter_obj.tmp_dir, '%s_blast.fasta'%node_rank)
		result_ls = self.blastOneBatchProbes(probe_id_seq_ls, \
											computing_parameter_obj.blast_bin_path, computing_parameter_obj.database_fname, \
											tmp_blast_infname, min_no_of_identities=computing_parameter_obj.min_no_of_identities, \
											node_rank=node_rank)
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return result_ls
	
	def outputBlastReturn(self, writer, result_ls):
		for result_entry in result_ls:
			data_row = result_entry[0].split("_")
			data_row.extend(result_entry[1:])
			writer.writerow(data_row)
			#outf.write(result_entry[0].replace("_",",")+","+",".join([str(a) for a in result_entry[1:]])+"\n")
	
	def output_node_handler(self, communicator, parameter_list, data):
		"""
		"""
		outf, writer = parameter_list[:2]
		result_ls = cPickle.loads(data)
		self.outputBlastReturn(writer, result_ls)
		outf.flush()
	
	def run(self):
		"""
		2010-4-14
		"""
		self.communicator = MPI.world.duplicate()
		node_rank = self.communicator.rank
		header = ['Chromosome', 'Position', 'Probe_ID', 'Alignment_title', "query_start", "query_end", 'Number_matches', 'Alignment_start', 'Alignment_stop']
		
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		
		data_to_pickle_name_ls = []	# 2010-4-15 data to send to computing nodes. sending "genome" and "probes" takes too long.
		if node_rank == 0:
			init_data = self.create_init_data(self.probe_fname, self.ref_genome_dir)
			if self.communicator.size==1:	# 2010-4-15 only one cpu. Run in serial.
				if self.debug:
					import pdb
					pdb.set_trace()
				tmp_blast_infname = os.path.join(self.tmp_dir, '%s_blast.fasta'%0)
				probe_index_start = 0 
				outf = open(self.output_fname, 'w')
				writer = csv.writer(outf, delimiter='\t')
				writer.writerow(header)
				no_of_probes = len(init_data.probes.data)
				while probe_index_start<no_of_probes:
					probe_index_stop = min(no_of_probes, probe_index_start + self.block_size)
					probe_id_seq_ls = self.get_probe_id_seq_ls(init_data.genome, init_data.probes, probe_index_start, probe_index_stop)
					result_ls = self.blastOneBatchProbes(probe_id_seq_ls, \
												self.blast_bin_path, self.database_fname, \
												tmp_blast_infname, min_no_of_identities = self.min_no_of_identities, \
												node_rank=node_rank)
					probe_index_start = probe_index_stop
					self.outputBlastReturn(writer, result_ls)
					outf.flush()
				del writer
				sys.exit(2)
			
			for data_to_pickle_name in data_to_pickle_name_ls:
				sys.stderr.write("passing %s ... \n"%(data_to_pickle_name))
				for node in free_computing_nodes:	#send it to the computing_node
					sys.stderr.write("passing initial data to nodes from %s to %s ..."%(node_rank, node))
					data_pickle = cPickle.dumps(getattr(init_data, data_to_pickle_name), -1)
					self.communicator.send(data_pickle, node, 0)
					del data_pickle
					sys.stderr.write(" Done.\n")
				sys.stderr.write("Done.\n")
		elif node_rank in free_computing_nodes:
			init_data = PassingData()
			for data_to_pickle_name in data_to_pickle_name_ls:
				data, source, tag = self.communicator.receiveString(0, 0)
				setattr(init_data, data_to_pickle_name, cPickle.loads(data))
				del data
		else:
			pass
		
		self.synchronize()
		
		if node_rank == 0:
			param_d = PassingData()
			param_d.index = 0
			param_d.probes = init_data.probes
			param_d.genome = init_data.genome
			parameter_list = [param_d]
			self.input_node(parameter_list, free_computing_nodes, message_size=self.block_size, input_handler=self.input_handler)
		elif node_rank in free_computing_nodes:
			computing_parameter_obj = PassingData(init_data=init_data, tmp_dir=self.tmp_dir,\
												blast_bin_path = self.blast_bin_path, database_fname=self.database_fname,\
												min_no_of_identities = self.min_no_of_identities)
			self.computing_node(computing_parameter_obj, self.computing_node_handler)
		else:
			outf = open(self.output_fname, 'w')
			writer = csv.writer(outf, delimiter='\t')
			writer.writerow(header)
			parameter_list = [outf, writer,]
			self.output_node(free_computing_nodes, parameter_list, self.output_node_handler)
			del writer
			outf.close()
		
		self.synchronize()	#to avoid some node early exits
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiBlast
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
	
