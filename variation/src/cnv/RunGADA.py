#!/usr/bin/env python
"""

Examples:
	RunGADA.py -i /Network/Data/250k/tmp-yh/CNV/call_method_17_CNV_array_intensity_norm_chr4.tsv
	-o /Network/Data/250k/tmp-yh/CNV/call_method_17_CNV_array_intensity_norm_chr4_GADA_out.tsv
	
	RunGADA.py -i ./call_method_17_CNV_array_intensity_chr4_line_no_888148_1107622_norm.tsv
	-o call_method_17_CNV_array_intensity_chr4_line_no_888148_1107622_norm_w1_34.out -w1,34 -u yh
	
	# 2009-10-5 test-run GADAJRN
	chr=1;M=10; echo ~/script/variation/src/RunGADA.py -i ~/panfs/250k/CNV/call_method_17_CNV_array_intensity_norm_chr$chr.tsv
	-o ~/panfs/250k/CNV/call_17_CNV_norm_intensity_chr$chr.GADAJRN_M$M.tsv -y2 -M $M
	
Description:
	2009-10-28 program to run GADA (Pique-Regi2008, Pique-Regi2009) to infer states/amplitude of CNV probes.
	One-sample version:
		Input is output from CNVNormalize.py.
		Output is tab-delimited:
			ecotype_id, array_id, chr_pos1, chr_pos2, length, amplitude, start_probe_id, end_probe_id
	Multi-sample version:
		Matlab has to be setup beforehand (PATH, license, etc.). MATLABPATH is also setup so that GADAJRNWrap could be called.
		Input is same as one-sample version.
		Output is tab-delimited.
			ecotype_id, array_id, chr_pos1, chr_pos2, length, amplitude, start_probe_id, end_probe_id
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import csv, numpy, traceback, subprocess
from pymodule import figureOutDelimiter, getListOutOfStr
from CNVNormalize import CNVNormalize
import cStringIO
from variation.src.db import Stock_250kDB

class RunGADA(CNVNormalize):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('input_fname', 1, ): ['', 'i', 1, 'CNV intensity matrix, probe X arrays. 1st column is probe id. 2nd last col is chr. last col is pos.', ],\
							('output_fname', 1, ): ['', 'o', 1, 'self-explanatory', ],\
							('which_array_id_ls', 0, ): [None, 'w', 1, 'list of array ids indicating which arrays(s) for GADA to work on. format: 0,1-3. Work on all if not given.' ],\
							('tmp_input_fname', 1, ): ['/tmp/GADA/GADA_input', 't', 1, 'temporary file to store GADA input. non-existent folder would be created if IO_thru_file is on.'],\
							('tmp_output_fname', 1, ): ['/tmp/GADA/GADA_output', 'm', 1, 'temporary file to store GADA output. non-existent folder would be created if IO_thru_file is on.', ],\
							('aAlpha', 1, float): [0.5, 'A', 1, 'a in Gamma(a;b) the function that controls the prior for the number of breakpoints', ],\
							('TBackElim', 1, float): [4, 'T', 1, '(amp1-amp2)/stddev in GADA', ],\
							('MinSegLen', 1, int): [5, 'M', 1, 'minimum no of probes to comprise a segment in GADA', ],\
							('GADA_path', 1, ): [os.path.expanduser('~/script/variation/bin/GADA/GADA'), 'G', 1, 'path to the one-sample version GADA program'],\
							('IO_thru_file', 0, int):[0, 'I', 0, 'mode of invoking GADA. default does not involve IO. This makes GADA-calling use files.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.'],\
							('run_type', 0, int):[1, 'y', 1, 'Run type 1: one-sample version GADA (c-version), 2: multi-sample (in matlab)']}
	
	def __init__(self, **keywords):
		"""
		2008-12-05
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
		self.which_array_id_ls = getListOutOfStr(self.which_array_id_ls, data_type=int)
	
	@classmethod
	def prepareGADAinput(cls, intensity_ls, tmp_input_fname=None, IO_thru_file=False, report=False):
		"""
		2009-10-28
			return a StringIO object if IO_thru_file is false or tmp_input_fname is not
		"""
		if report:
			sys.stderr.write("Preparing input for GADA ...")
		
		if not IO_thru_file or not tmp_input_fname:
			of = cStringIO.StringIO()
		else:
			of = open(tmp_input_fname, 'w')
		no_of_rows = len(intensity_ls)
		
		for i in range(no_of_rows):
			of.write('%s\n'%intensity_ls[i])
		if report:
			sys.stderr.write("Done.\n")
		
		if not IO_thru_file or not tmp_input_fname:
			of.seek(0)	# set the current position to beginning
			return of
		else:
			of.close()
			return tmp_input_fname
	
	@classmethod
	def _GADA(cls, GADA_path, tmp_input_fname, tmp_output_fname, aAlpha=0.8, TBackElim=2, MinSegLen=2, \
			IO_thru_file=False, report=False):
		"""
		2009-11-10
			GADA no longer uses PIPE (stdin, stdout). Specify input and output file directly.
		2009-11-6
			If tmp_input_fname is a file, call Popen() in a different way.
		2009-10-28
			add argument IO_thru_file
		2009-10-26
			remove -c (segment classification into gain/loss) in invoking GADA. Segment classification is done post-GADA.
		2009-10-5
			add arguments: aAlpha=0.5, TBackElim=4, MinSegLen=5
		"""
		if report:
			sys.stderr.write("Running GADA ...")
		commandline = '%s -a %s -T %s -M %s -s -0.4 -b 0.0 -i %s -o %s'%(GADA_path, aAlpha, TBackElim, MinSegLen, tmp_input_fname, \
																		tmp_output_fname)
			# -s specifies he variance estimate for the noise. If negative it will be estimated from the provided data
			# -b Mean amplitude associated to the Neutral state. useless if -c is not specified.
			# -c Classify segments into altered state (L)oss, (N)eutral, (G)ain). If c option is not specified, only mean is returned.
		isInputFileHandler = False
		if type(tmp_input_fname) is str:
			inf = open(tmp_input_fname)
			isInputFileHandler = True
		elif type(tmp_input_fname) is cStringIO.OutputType:	# tmp_input_fname is cStringIO
			inf = tmp_input_fname
			isInputFileHandler= False
		else:
			 sys.stderr.write("GADA Input is neither file nor StringIO. skip.\n")
			 return ''
		if isInputFileHandler:
			command_handler = subprocess.Popen(commandline, shell=True, stdin=inf, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
			stdout_content, stderr_content = command_handler.communicate()
		else:
			command_handler = subprocess.Popen(commandline, shell=True, stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
			stdout_content, stderr_content = command_handler.communicate()
			# stdout_content, stderr_content = command_handler.communicate(input=inf.read())
		
		del inf
		if stderr_content:
			sys.stderr.write('stderr of %s: %s \n'%(commandline, stderr_content))
		if report:
			sys.stderr.write("Done.\n")
		"""
		if IO_thru_file:
			of = open(tmp_output_fname, 'w')
			of.write(stdout_content)
			del of
		"""
		return_content = ""
		of = open(tmp_output_fname)
		for line in of:
			return_content += line
		return return_content
		#return stdout_content
	
	
	def findOutWhichColumn(cls, col_id_ls, which_array_id_set):
		"""
		2008-12-05
		"""
		col_index_to_return_ls = []
		for i in range(len(col_id_ls)):
			col_id = int(col_id_ls[i])
			if col_id in which_array_id_set:
				col_index_to_return_ls.append(i)
		return col_index_to_return_ls
	findOutWhichColumn = classmethod(findOutWhichColumn)
	
	@classmethod
	def output_header(cls, writer, chr_pos_ls=None):
		"""
		2010-6-29
			add the median-output 
		2010-2-10
			becomes classmethod
		2009-10-26
			no more outputting the whole chr_pos_ls
		2008-12-05
		"""
		if getattr(cls, 'report', False):
			sys.stderr.write("Outputting header ...")
		header = ['ecotype_id', 'array_id', "start_probe", "end_probe", "length", "amplitude", \
				"start_probe_id", "end_probe_id", "median"]
		"""
		for chr_pos in chr_pos_ls:
			chr_pos = '_'.join(chr_pos)
			header.append(chr_pos)
		"""
		writer.writerow(header)
		if getattr(cls, 'report', False):
			sys.stderr.write("Done.\n")
	state2number = {'G': 1, 'N':0, 'L':-1}
	
	def calculateMannWhitney(self, GADA_output, intensity_ls, segment1_index, segment2_index, segment1OnTheLeft=True,\
							maxNoOfProbesInSegment1=40):
		"""
		2010-6-29
			segment2_index is the designated deletion segment.
		"""
		segment1 = GADA_output[segment1_index]
		segment2 = GADA_output[segment2_index]
		segment1_probe1_index, segment1_probe2_index, segment1_length = segment1[:3]
		if segment1_length>maxNoOfProbesInSegment1:
			segment1_length = maxNoOfProbesInSegment1
			if segment1OnTheLeft:
				segment1_probe1_index = segment1_probe2_index-maxNoOfProbesInSegment1+1
			else:
				segment1_probe2_index = segment1_probe1_index+maxNoOfProbesInSegment1-1
		segment1_intensityLs = intensity_ls[segment1_probe1_index-1:segment1_probe2_index]
		
		segment2_probe1_index, segment2_probe2_index, segment2_length = segment2[:3]
		segment2_intensityLs = intensity_ls[segment2_probe1_index-1:segment2_probe2_index]
		import rpy
		result = rpy.r.wilcox_test(segment1_intensityLs, segment2_intensityLs, alternative="g")
		rho_stat = result['statistic']['W']/(segment1_length*segment2_length)
		result['rho_stat'] = rho_stat
		return result
	
	def addMedianToGADAOutput(self, GADA_output, intensity_ls, maxDeletionMedianIntensity=-0.1):
		"""
		2010-6-29
		"""
		no_of_rows = len(GADA_output)
		for i in range(no_of_rows):
			row = GADA_output[i]
			if type(row[0])==str and (row[0].find('#')==0 or row[0].find('Start')==0):	# ignore the comments
				#skip the first line of real data. it's header "Start   Stop    Lenght  Ampl    State"
				continue
			probe1_index, probe2_index, length, amplitude = row[:4]
			probe1_index = int(probe1_index)-1
			probe2_index = int(probe2_index)-1
			median = numpy.median(intensity_ls[probe1_index:probe2_index+1])
			row.append(median)
			"""
			if median<=maxDeletionMedianIntensity:
				if i>0:
					leftWilcoxResult = self.calculateMannWhitney(GADA_output, intensity_ls, i-1, i, segment1OnTheLeft=True)
					leftWilcoxRhoStat = leftWilcoxResult['rho_stat']
					leftWilcoxPvalue = leftWilcoxResult['p.value']
				else:
					leftWilcoxRhoStat = ''
					leftWilcoxPvalue = ''
				if i<no_of_rows-1:
					rightWilcoxResult = self.calculateMannWhitney(GADA_output, intensity_ls, i+1, i, segment1OnTheLeft=False)
					rightWilcoxRhoStat = rightWilcoxResult['rho_stat']
					rightWilcoxPvalue = rightWilcoxResult['p.value']
				else:
					rightWilcoxRhoStat = ''
					rightWilcoxPvalue = ''
				row += [leftWilcoxRhoStat, leftWilcoxPvalue, rightWilcoxRhoStat, rightWilcoxPvalue]
			"""
			GADA_output[i] = row
	
	def output_GADA_output(self, writer, GADA_output, array_id, ecotype_id, chr_pos_ls=None, probe_id_ls=None,\
						base_start_index=0):
		"""
		2010-6-29
			add the median-output 
		2010-6-5
			add argument base_start_index, default=0
		2009-11-22
			if GADA_output is string, turn it into csv.reader.
			if GADA_output is a list returned from GADA.GADA.run(), nothing.
		2009-10-28
			return no_of_segments
		2009-10-26
			add chr_pos_ls, probe_id_ls to improve output
			no more output to writer_amp
		2008-12-05
		"""
		if self.report:
			sys.stderr.write("Outputting GADA output ...")
		if type(GADA_output)==str:
			GADA_output = cStringIO.StringIO(GADA_output)
			GADA_output = csv.reader(GADA_output, delimiter='\t')
		counter = 0
		for row in GADA_output:
			if type(row[0])==str and (row[0].find('#')==0 or row[0].find('Start')==0):	# ignore the comments
				#skip the first line of real data. it's header "Start   Stop    Lenght  Ampl    State"
				continue
			counter += 1
			probe1_index, probe2_index, length, amplitude = row[:4]
			probe1_index = int(probe1_index)-1 + base_start_index
			probe2_index = int(probe2_index)-1 + base_start_index
			probe1 = chr_pos_ls[probe1_index]
			probe1 = '_'.join(probe1)
			probe2 = chr_pos_ls[probe2_index]
			probe2 = '_'.join(probe2)
			
			probe1_id = probe_id_ls[probe1_index]
			probe2_id = probe_id_ls[probe2_index]
			
			new_row = [ecotype_id, array_id, probe1, probe2, length, amplitude, probe1_id, probe2_id] + row[4:]
			
			"""
			start = int(row[0])
			stop = int(row[1])
			amplitude = row[3]
			state_number = self.state2number[state]
			for i in range(start-1, stop):
				data_row.append(state_number)
				data_row_amp.append(amp)
			"""
			writer.writerow(new_row)
		#writer.writerow(data_row)
		#writer_amp.writerow(data_row_amp)
		if self.report:
			sys.stderr.write("Done.\n")
		no_of_segments = counter
		return no_of_segments
	
	def oneSampleGADA(self, ):
		"""
		import GADA	#2010-4-29 moved here to avoid import errors if this module is not compiled well
		2009-11-23
			if self.IO_thru_file:
				call self._GADA()
			else:
				call GADA.GADA.run()
		2009-9-27
			split out of run()
		"""
		data_matrix, probe_id_ls, chr_pos_ls, header = self.get_input(self.input_fname)
		array_id_ls = header[1:-2]
		array_id_ls = map(int, array_id_ls)
		if self.which_array_id_ls:
			col_index_ls = self.findOutWhichColumn(array_id_ls, set(self.which_array_id_ls))
		else:
			col_index_ls = range(data_matrix.shape[1])
		
		#output_amp_fname = '%s_amp.tsv'%os.path.splitext(self.output_fname)[0]
		#writer_amp  = csv.writer(open(output_amp_fname, 'w'), delimiter='\t')
		writer_amp = None
		#self.output_header(writer_amp, chr_pos_ls)
		writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
		self.output_header(writer, chr_pos_ls)
		
		# 2009-10-27 create directories for temp input/output files
		if self.IO_thru_file:
			tmp_input_dir = os.path.split(self.tmp_input_fname)[0]
			tmp_output_dir = os.path.split(self.tmp_output_fname)[0]
			if not os.path.isdir(tmp_input_dir):
				os.makedirs(tmp_input_dir)
			if not os.path.isdir(tmp_output_dir):
				os.makedirs(tmp_output_dir)
		
		for col_index in col_index_ls:
			array_id = array_id_ls[col_index]
			array = Stock_250kDB.ArrayInfo.get(array_id)
			if array:
				ecotype_id = array.maternal_ecotype_id
			else:
				continue
			if self.IO_thru_file:
				tmp_input_fname = '%s_%s'%(self.tmp_input_fname, array_id)
				tmp_output_fname = '%s_%s'%(self.tmp_output_fname, array_id)
				input_file = self.prepareGADAinput(data_matrix[:,col_index], tmp_input_fname, IO_thru_file=self.IO_thru_file)
				GADA_output = self._GADA(self.GADA_path, input_file, tmp_output_fname, self.aAlpha, \
										self.TBackElim, self.MinSegLen, IO_thru_file=self.IO_thru_file)
			else:
				import GADA	#2010-4-29 moved here to avoid import errors if this module is not compiled well
				ins = GADA.GADA(self.debug)
				GADA_output = ins.run(map(float, data_matrix[:,col_index]), self.aAlpha, self.TBackElim, self.MinSegLen)
				del ins
			self.output_GADA_output(writer, GADA_output, array_id, ecotype_id, chr_pos_ls, probe_id_ls)
		del writer, writer_amp
	
	def multiSampleGADA(self):
		"""
		2009-9-27
			run the matlab-version new multi-sample GADA
		"""
		data_matrix, probe_id_ls, chr_pos_ls, header = self.get_input(self.input_fname)
		array_id_ls = header[1:-2]
		array_id_ls = map(int, array_id_ls)
		if self.which_array_id_ls:
			col_index_ls = self.findOutWhichColumn(array_id_ls, set(self.which_array_id_ls))
		else:
			col_index_ls = range(data_matrix.shape[1])
		
		import time
		current_time = int(time.time()*10000)
		
		#2009-10-5 create a memmapfile
		sys.stderr.write("Converting input into matlab memmapfile ...")
		m_memmap_fname = ('%s.memmap'%os.path.splitext(self.input_fname)[0])
		if not os.path.isfile(m_memmap_fname):	# create a new one if the extant one doesn't exist
			m_memmap_fname = '/tmp/mGADA_%s_memmap.tsv'%current_time
			data_matrix = numpy.transpose(data_matrix)	#somehow, without this, matlab will read it as transposed to the original matrix
			outf = open(m_memmap_fname, 'wb')
			data_matrix.tofile(outf, format='%.8f')	# format seems to not matter.
			del outf
			del data_matrix
		sys.stderr.write("Done.\n")
		
		m_output_fname = '/tmp/mGADA_%s_output.tsv'%current_time
		
		commandline = 'matlab -nodisplay -nojvm -r "GADAJRNWrap %s %s %s %s %s %s %s"'%(m_memmap_fname,  m_output_fname, self.aAlpha, \
																					self.TBackElim,\
											self.MinSegLen, len(probe_id_ls), len(array_id_ls))
		# 2009-10-5 assume matlab is setup (PATH, license, etc.). MATLABPATH is also setup so that GADAJRNWrap could be called.
		command_handler = subprocess.Popen(commandline, shell=True, stderr=sys.stderr, stdout=sys.stdout)
		stdout_content, stderr_content = command_handler.communicate()
		if stdout_content:
			sys.stderr.write('stdout of %s: %s \n'%(commandline, stdout_content))
		if stderr_content:
			sys.stderr.write('stderr of %s: %s \n'%(commandline, stderr_content))
		if self.report:
			sys.stderr.write("Done.\n")
		
		
		sys.stderr.write("Beautifying matlab output ...")
		reader = csv.reader(open(m_output_fname, 'r'), delimiter='\t')
		writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
		self.output_header(writer)
		array_id2array = {}
		for row in reader:
			sample_index, probe1_index, probe2_index, length, amplitude = row[:5]
			array_id = array_id_ls[int(sample_index)-1]	# MATLAB index starts from 1 for list or array
			if array_id in array_id2array:
				array = array_id2array.get(array_id)
			else:
				array = Stock_250kDB.ArrayInfo.get(array_id)
				array_id2array[array_id] = array
			
			if array:
				ecotype_id = array.maternal_ecotype_id
			else:
				continue
			probe1_index = int(probe1_index)-1
			probe2_index = int(probe2_index)-1
			probe1 = chr_pos_ls[probe1_index]
			probe1 = '_'.join(probe1)
			probe2 = chr_pos_ls[probe2_index]
			probe2 = '_'.join(probe2)
			
			probe1_id = probe_id_ls[probe1_index]
			probe2_id = probe_id_ls[probe2_index]
			
			new_row = [ecotype_id, array_id, probe1, probe2, length, amplitude, probe1_id, probe2_id]
			writer.writerow(new_row)
		del writer, reader
		sys.stderr.write("Done.\n")
		
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session		
		if self.run_type==1:
			self.oneSampleGADA()
		elif self.run_type==2:
			self.multiSampleGADA()
		

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = RunGADA
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()