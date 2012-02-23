#!/usr/bin/env python
"""
Examples:
	%s  -o /tmp/testPrefix /tmp/overlapStat*
	
	%s 
	

Description:
	2011-9-29
		a list of input files are appended in the end
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])


bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:	   #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import csv, numpy
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter
from pymodule import yh_matplotlib


class PlotResultPeakOverlap(object):
	__doc__ = __doc__
	option_default_dict = {('outputFnamePrefix', 1, ): [None, 'o', 1, 'output filename prefix for all figures.'],\
						('title', 0, ): [None, 't', 1, 'title for the figure.'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}


	def __init__(self, inputFnameLs, **keywords):
		"""
		2011-7-11
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		self.inputFnameLs = inputFnameLs
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		no_of_result1_peaks_ls = []
		no_of_result2_peaks_ls = []
		fraction_of_result1_peaks_in_result2_ls = []
		fraction_of_result2_peaks_in_result1_ls = []
		no_of_combined_peaks_ls = []
		fraction_of_overlap_in_combined_peaks_ls = []
		
		for inputFname in self.inputFnameLs:
			reader = csv.reader(open(inputFname), delimiter=figureOutDelimiter(inputFname))
			header = reader.next()
			col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
			no_of_result1_peaks_index = col_name2index.get("no_of_result1_peaks")
			no_of_result2_peaks_index = col_name2index.get("no_of_result2_peaks")
			no_of_result1_peaks_in_result2_index = col_name2index.get("no_of_result1_peaks_in_result2")
			no_of_result2_peaks_in_result1_index = col_name2index.get("no_of_result2_peaks_in_result1")
			for row in reader:
				no_of_result1_peaks = float(row[no_of_result1_peaks_index])
				no_of_result2_peaks = float(row[no_of_result2_peaks_index])
				no_of_result1_peaks_in_result2 = float(row[no_of_result1_peaks_in_result2_index])
				no_of_result2_peaks_in_result1 = float(row[no_of_result2_peaks_in_result1_index])
				no_of_result1_peaks_ls.append(no_of_result1_peaks)
				no_of_result2_peaks_ls.append(no_of_result2_peaks)
				fraction_of_result1_peaks_in_result2_ls.append(no_of_result1_peaks_in_result2/no_of_result1_peaks)
				fraction_of_result2_peaks_in_result1_ls.append(no_of_result2_peaks_in_result1/no_of_result2_peaks)
				no_of_combined_peaks_ls.append(no_of_result1_peaks + no_of_result2_peaks)
				fraction_of_overlap_in_combined_peaks_ls.append((no_of_result1_peaks_in_result2 + no_of_result2_peaks_in_result1)/(no_of_result1_peaks + no_of_result2_peaks))
			del reader
		
		title="%s pairs"%(len(fraction_of_result1_peaks_in_result2_ls))
		if len(fraction_of_result1_peaks_in_result2_ls)>10:
			medianFraction = numpy.median(fraction_of_result1_peaks_in_result2_ls)
			title += " median %.3f"%(medianFraction)
		yh_matplotlib.drawHist(fraction_of_result1_peaks_in_result2_ls, title=title, \
						xlabel_1D="fraction of result1 peaks in result2", xticks=None, \
						outputFname="%s_hist_of_fraction_of_result1_peaks_in_result2.png"%self.outputFnamePrefix, \
						min_no_of_data_points=20, needLog=False, \
						dpi=200)
		title="%s pairs"%(len(fraction_of_result2_peaks_in_result1_ls))
		if len(fraction_of_result2_peaks_in_result1_ls)>10:
			medianFraction = numpy.median(fraction_of_result2_peaks_in_result1_ls)
			title += " median %.3f"%(medianFraction)
		yh_matplotlib.drawHist(fraction_of_result2_peaks_in_result1_ls, title=title, \
						xlabel_1D="fraction of result2 peaks in result1", xticks=None, \
						outputFname="%s_hist_of_fraction_of_result2_peaks_in_result1.png"%self.outputFnamePrefix, \
						min_no_of_data_points=20, needLog=False, \
						dpi=200)
		
		title="%s pairs"%(len(fraction_of_overlap_in_combined_peaks_ls))
		if len(fraction_of_overlap_in_combined_peaks_ls)>10:
			medianFraction = numpy.median(fraction_of_overlap_in_combined_peaks_ls)
			title += " median %.3f"%(medianFraction)
		yh_matplotlib.drawHist(fraction_of_overlap_in_combined_peaks_ls, title=title, \
						xlabel_1D="fraction of recurrent peaks in combined", xticks=None, \
						outputFname="%s_hist_of_fraction_of_recurrent_peaks_in_combined.png"%self.outputFnamePrefix, \
						min_no_of_data_points=20, needLog=False, \
						dpi=200)
		
		title="%s results"%(len(no_of_result1_peaks_ls))
		yh_matplotlib.drawScatter(no_of_result1_peaks_ls, no_of_result2_peaks_ls, \
				fig_fname="%s_no_of_peaks_result1_vs_result2.png"%self.outputFnamePrefix, \
				title=title, xlabel='No. of peaks in result1', \
				ylabel='No. of peaks in result2', dpi=300)
		
		title="%s results"%(len(no_of_result1_peaks_ls))
		yh_matplotlib.drawScatter(no_of_result1_peaks_ls, fraction_of_result1_peaks_in_result2_ls, \
				fig_fname="%s_result1_no_of_peak_vs_fraction.png"%self.outputFnamePrefix, \
				title=title, xlabel='No. of peaks in result1', \
				ylabel='Fraction found in result2', dpi=300)
		
		title="%s results"%(len(no_of_result2_peaks_ls))
		yh_matplotlib.drawScatter(no_of_result2_peaks_ls, fraction_of_result2_peaks_in_result1_ls, \
				fig_fname="%s_result2_no_of_peak_vs_fraction.png"%self.outputFnamePrefix, \
				title=title, xlabel='No. of peaks in result2', \
				ylabel='Fraction found in result1', dpi=300)
		
		title="%s pairs"%(len(fraction_of_result1_peaks_in_result2_ls))
		yh_matplotlib.drawScatter(fraction_of_result1_peaks_in_result2_ls, fraction_of_result2_peaks_in_result1_ls, \
				fig_fname="%s_1_fraction_in2_vs_2_fraction_in1.png"%self.outputFnamePrefix, \
				title=title, xlabel='result1 fraction found in result2', \
				ylabel='result2 fraction found in result1', dpi=300)
		
		title="%s pairs"%(len(no_of_combined_peaks_ls))
		yh_matplotlib.drawScatter(no_of_combined_peaks_ls, fraction_of_overlap_in_combined_peaks_ls, \
				fig_fname="%s_combined_no_of_peak_vs_fraction.png"%self.outputFnamePrefix, \
				title=title, xlabel='No. of peaks combined', \
				ylabel='Fraction recurrent', dpi=300)
		

if __name__ == '__main__':
	main_class = PlotResultPeakOverlap
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
