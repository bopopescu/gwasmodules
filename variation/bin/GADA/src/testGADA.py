#!/usr/bin/env python
"""

Examples:
	testGADA.py -i input.txt -o output_a0.5T4M5.tsv
	
	
Description:
	This is to test importing GADA.so to see if it works.
	
	Each line in the input file is just one float number for one probe.
	
	Output is tab-delimited.
		starting-probe-index, stop-probe-index, no-of-probes, amplitude
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))


from pymodule import MatrixFile
import GADA

class testGADA(object):
	__doc__ = __doc__
	option_default_dict = {('input_fname', 1, ): ['', 'i', 1, 'CNV intensity matrix, probe X arrays. 1st column is probe id. 2nd last col is chr. last col is pos.', ],\
						('output_fname', 1, ): ['', 'o', 1, 'self-explanatory', ],\
						('aAlpha', 1, float): [0.5, 'A', 1, 'a in Gamma(a;b) the function that controls the prior for the number of breakpoints', ],\
						('TBackElim', 1, int): [4, 'T', 1, '(amp1-amp2)/stddev in GADA', ],\
						('MinSegLen', 1, int): [5, 'M', 1, 'minimum no of probes to comprise a segment in GADA', ],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.'],\
							}
	
	
	def __init__(self, **keywords):
		"""
		2009-10-28
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def run(self):
		"""
		2010-6-9
		"""
		
		import sys, csv
		#input_fname = "/tmp/GADA_ATL8C17938_ATL7C28313_1_5_ATL8C17938_ATL7C28313_1_input"
		inf = MatrixFile(inputFname=self.input_fname)
		intensity_ls = []
		for row in inf:
			intensity_ls.append(float(row[0]))
		inf.close()
		sys.stderr.write("%s probes read in.\n"%len(intensity_ls))
		
		ins = GADA.GADA()
		segment_ls = ins.run(intensity_ls, self.aAlpha, self.TBackElim, self.MinSegLen)
		del ins
		writer = MatrixFile(inputFname=self.output_fname, openMode='w', delimiter='\t')
		writer.writerow(["# Parameter setting: a=%s,T=%s,MinSegLen=%s"%(self.aAlpha, self.TBackElim, self.MinSegLen)])
		header = ['Start', 'Stop', 'Length','Ampl']
		writer.writeHeader(header)
		for segment in segment_ls:
			writer.writerow(segment)
		del writer

if __name__ == '__main__':
	#test GADA
	ins = GADA.GADA(1)
	segment_ls = ins.run([1,1,1,1,0.99,0.99,1,1,0.1,0.1,0.1,0.15], 0.2, 4, 2)
	print segment_ls
	del ins

	from pymodule import ProcessOptions
	main_class = testGADA
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
