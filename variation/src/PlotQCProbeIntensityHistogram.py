#!/usr/bin/env python

"""
Examples:
	# draw intensity histogram for QC probes
	PlotQCProbeIntensityHistogram.py -i 390-399 -o ~/Desktop/array_intensity_histogram/QC_probe_hist -m 50000 -l

	# draw intensity histogram for SNP probes
	PlotQCProbeIntensityHistogram.py -i 498,500 -o ~/Desktop/array_intensity_histogram/SNP_probe_hist 
		-p ~/script/affy/250k_test/250kprobe_subset.txt -l

	# draw intensity histogram for all probes (QC + CNV + SNP) of Ler arrays, set xlim within 1-5
	PlotQCProbeIntensityHistogram.py -i 3,4,41,150,151 -o /tmp/all_probe_hist -m 100000 -y 3 -l -x 1,5
	
Description:
	Run type:
	1. intensity histogram of QC probes given an array id
		(run c-program ~/script/affy/sdk/calvin_files/OutputQCIntensity)
	2. intensity histogram of probes given a probe subset file
		(run c-program ~/script/affy/sdk/calvin_files/ReadatSNPtilgeno)
	3. intensity histogram of all probes in an array (QC + SNP + CNV)
		require ReadatSNPtilgeno.so (symbolic link to ~/script/affy/sdk/calvin_files/ReadatSNPtilgeno.so)
	
	Note: run_type 1 & 2 requires cdf_fname. 2 further requires probe_subset_fname is present.
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv, stat, getopt
import traceback, gc, subprocess
from pymodule import getListOutOfStr

class PlotQCProbeIntensityHistogram:
	__doc__ = __doc__
	option_default_dict = {('array_id_ls', 1, ): ['', 'i', "coma/dash-separated list of array ids. like 1,3-100"],\
			('cdf_fname', 1, ):[os.path.expanduser('~/script/variation/genotyping/250ksnp/data/atSNPtilx520433_rev2/Full/atSNPtil_geno/LibFiles/atSNPtil_geno.cdf')],\
			('figure_ofname_prefix', 1, ): ['', 'o', 1, 'Figure output filename. "%s_%s.png"%(prefix, array_id)', ],\
			('input_dir', 1, ): ['/Network/Data/250k/db/raw_data/', 'n', 1, 'Directory where cel files are stored. cel filename: "%s_raw_data.cel"%array_id.' ],\
			('no_of_bins', 1, int): [100, 'f', 1, 'Number of bins in the histogram plot'],\
			('max_intensity', 1, int,):[25000, 'm', 1, 'Maximum intensity to put into plot' ],\
			('probe_subset_fname', 0, ):['', 'p', 1, 'probe subset filename, if given, \
				will draw histogram of those probes rather than QC probes. Example: ~/script/affy/250k_test/250kprobe_subset.txt' ],\
			('run_type', 1, int):[1, 'y', 1, '1: intensity of QC probes, \
					2: intensity of probes which belong to subsets in probe_subset_fname, \
					3: intensity of all probes'],\
			('logTransform', 0, int):[0, 'l', 0, 'toggle to log10-transform probe intensity'],\
			('xlim', 0, str):['', 'x', 1, 'set the range of x-axis, like 1,5. default is automatic.'],\
			('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
			('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self,  **keywords):
		"""
		2008-07-05
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, \
														error_doc=self.__doc__, class_to_have_attr=self)
		
		self.getIntensityFuncDict = {1: self.getQCIntensityLs,
									2: self.getIntensityLsFromOutputGivenProbeSubset,
									3: self.getAllIntensityFromArray}
		
		self.xlim = getListOutOfStr(self.xlim, data_type=float)
		self.array_id_ls = getListOutOfStr(self.array_id_ls, int)
		
	def getIntensityLsFromOutputGivenProbeSubset(self, intensity_output_fname, max_intensity=None, logTransform=True):
		"""
		2010-5-6
			add argument logTransform
		2007-07-18
		"""
		sys.stderr.write("Getting intensity from %s ... "%(intensity_output_fname))
		reader = csv.reader(open(intensity_output_fname), delimiter='\t')
		reader.next()
		intensity_ls = []
		for row in reader:
			for i in range(1,len(row)):
				intensity = float(row[i])
				if max_intensity is not None and intensity>max_intensity:
					continue
				if logTransform:
					intensity = math.log10(intensity)
				intensity_ls.append(intensity)
		sys.stderr.write("Done.\n")
		return intensity_ls
	
	def getQCIntensityLs(self, qc_intensity_fname, max_intensity=None, logTransform=True):
		"""
		2010-5-6
			add argument logTransform
		2007-07-18
		"""
		reader = csv.reader(open(qc_intensity_fname), delimiter='\t')
		reader.next()
		intensity_ls = []
		for row in reader:
			intensity = float(row[3])
			if max_intensity is not None and intensity>max_intensity:
				continue
			if logTransform:
				intensity = math.log10(intensity)
			intensity_ls.append(intensity)
		return intensity_ls
	
	def getAllIntensityFromArray(self, array_cel_fname, max_intensity=None, logTransform=True):
		"""
		2010-5-6
			directly use the boost python module ReadatSNPtilgeno.
			get all probes.
		"""
		from ReadatSNPtilgeno import ReadatSNPtilgeno
		ins = ReadatSNPtilgeno(123)	# 123 is a fake ecotype id.
		old_intensity_ls = ins.returnAllProbeIntensty(array_cel_fname)
		intensity_ls = []
		for intensity in old_intensity_ls:
			if max_intensity is not None and intensity>max_intensity:
				continue
			if logTransform:
				intensity = math.log10(intensity)
			intensity_ls.append(intensity)
		return intensity_ls
	
	def _plot(self, intensity_ls, array_id, figure_ofname, no_of_bins=25, max_intensity=25000, xlim=None):
		"""
		2010-5-6
			add xlabel, ylabel
			add argument xlim, if present, set the xlim of pylab
			add max_intensity into title
		2008-07-05
		"""
		import pylab
		pylab.clf()
		pylab.xlabel('log10(intensity)')
		pylab.ylabel('count')
		pylab.title("Array: %s. max intensity: %s"%(array_id, max_intensity))
		pylab.hist(intensity_ls, no_of_bins)
		if xlim and len(xlim)==2:
			pylab.xlim(xlim)
		pylab.savefig(figure_ofname, dpi=100)
		#pylab.show()

	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		for array_id in self.array_id_ls:
			array_cel_fname = os.path.join(self.input_dir, '%s_raw_data.cel'%array_id)
			
			if self.run_type==1 or self.run_type==2:
				intensity_output_fname = '/tmp/%s_%s.tsv'%(array_id, self.run_type)
				if self.run_type==2:
					if os.path.isfile(self.probe_subset_fname):
						command = os.path.expanduser('~/script/affy/sdk/calvin_files/ReadatSNPtilgeno')
						command_ls = [command, '-i', array_cel_fname, '-d', self.cdf_fname, '-e', repr(array_id), \
											'-p', self.probe_subset_fname, '-o', intensity_output_fname ]
					else:
						sys.stderr.write("probe_subset_fname is missing for run type %s.\n"%self.run_type)
						continue
				else:
					OutputQCIntensity_command = os.path.expanduser('~/script/affy/sdk/calvin_files/OutputQCIntensity')
					command_ls = [OutputQCIntensity_command, '-i', array_cel_fname, '-d', self.cdf_fname, '-o', \
											intensity_output_fname ]
				cmd_p = subprocess.Popen(command_ls, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
				cmd_p.wait()
				cmd_out = cmd_p.stdout.read()
				print cmd_out
				cmd_err = cmd_p.stderr.read()
				sys.stderr.write(cmd_err)
			elif self.run_type==3:
				intensity_output_fname = array_cel_fname	# run_type 3 just needs the array_cel_fname
			else:
				intensity_output_fname = None
			
			if os.path.isfile(intensity_output_fname):	#2010-4-8 make sure the file exists
				intensity_ls = self.getIntensityFuncDict[self.run_type](intensity_output_fname, self.max_intensity,\
																	logTransform=self.logTransform)
				if intensity_ls:
					figure_ofname = '%s_%s.png'%(self.figure_ofname_prefix, array_id)
					self._plot(intensity_ls, array_id, figure_ofname, no_of_bins=self.no_of_bins, \
						max_intensity=self.max_intensity, \
						xlim=self.xlim)
				else:
					sys.stderr.write("intensity_ls for array %s is empty.\n"%array_id)
			else:
				sys.stderr.write("%s doesn't exist after executing %s.\n"%(intensity_output_fname, ' '.join(command_ls)))
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PlotQCProbeIntensityHistogram
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()