#!/usr/bin/env python
"""
Examples:
	#
	%s -i data/1MbBAC_as_ref/454_illu_6_sub_vs_1MbBAC.GATK.call -o /tmp/454_illu_6_sub_vs_1MbBAC.GATK.matrix -m 0.2 -n 0.4 -c1
	
	
	%s 
Description:
	Calculate pairwise distance matrix of a SNP (row) X Strain (column) matrix.
"""
import os, sys, numpy
__doc__ = __doc__%(sys.argv[0], sys.argv[0])
#2007-03-05 common codes to initiate database connection
import sys, os, math
#sys.path.insert(0, os.path.expanduser('~/lib/python'))
#sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
#sys.path.insert(0, os.path.join(os.path.expanduser('~/script/test/python')))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script/variation/src')))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

#import matplotlib; matplotlib.use("Agg")	#to avoid popup and collapse in X11-disabled environment


class CalculatePairwiseDistanceOutOfSNPXStrainMatrix(object):
	__doc__ = __doc__
	option_default_dict = {
						('inputFname', 1, ): ['', 'i', 1, 'common input file.', ],\
							('outputFname', 1, ): ['', 'o', 1, 'common output file', ],\
							('min_MAF', 1, float): [0.0, 'n', 1, 'minimum MAF for SNP filter', ],\
							('max_NA_rate', 1, float): [0.4, 'm', 1, 'maximum NA rate for SNP filter', ],\
							('convertHetero2NA', 1, int):[0, 'c', 1, 'toggle convertHetero2NA mode'],\
							('hetHalfMatchDistance', 1, float): [0.5, 'H', 1, 'distance between two half-matched genotypes. AG vs A or AG vs AC', ],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self, **keywords):
		"""
		2008-12-05
		This class is the entry to all others.
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	@classmethod
	def calculatePairwiseDistanceOutOfSNPXStrainMatrix(cls, inputFname, outputFname, snpFnameToDoFiltering=None, \
											convertHetero2NA=False,
											max_NA_rate=0.4, min_MAF=0.2, hetHalfMatchDistance=0.5):
		"""
		2011-10-20
			add argument hetHalfMatchDistance
		2011-4-7 output the pairwise distance as matrix
		2011-3-30
			add argument convertHetero2NA, max_NA_rate, min_MAF
		2011-3-29
			inputFname format:
				first row is sample id.
				first and 2nd column are the same locus id.
			genotypes are in alphabets.
		"""
		from pymodule.SNP import SNPData, transposeSNPData
		if snpFnameToDoFiltering:
			oldSNPData = SNPData(input_fname=snpFnameToDoFiltering, turn_into_array=1, ignore_2nd_column=1,\
								input_alphabet=1, turn_into_integer=1)
			snpData = transposeSNPData(oldSNPData)
			if convertHetero2NA:
				snpData = SNPData.convertHetero2NA(snpData)
			
			snpData = snpData.removeColsByNARate(snpData, max_NA_rate=max_NA_rate)
			if min_MAF>0 and min_MAF<0.5:
				snpData = snpData.removeColsByMAF(snpData, min_MAF=min_MAF)
			col_id_ls = snpData.col_id_ls
			
			
			oldSNPData = SNPData(input_fname=inputFname, turn_into_array=1, ignore_2nd_column=1,\
								input_alphabet=1, turn_into_integer=1)
			snpData = transposeSNPData(oldSNPData)
			if convertHetero2NA:
				snpData = SNPData.convertHetero2NA(snpData)
			snpData = SNPData.keepColsByColID(snpData, col_id_ls)
			
		else:
			oldSNPData = SNPData(input_fname=inputFname, turn_into_array=1, ignore_2nd_column=1,\
								input_alphabet=1, turn_into_integer=1)
			snpData = transposeSNPData(oldSNPData)
			if convertHetero2NA:
				snpData = SNPData.convertHetero2NA(snpData)
			
			snpData = snpData.removeColsByNARate(snpData, max_NA_rate=max_NA_rate)
			if min_MAF>0 and min_MAF<0.5:
				snpData = snpData.removeColsByMAF(snpData, min_MAF=min_MAF)
		
		# add outputFname to function below to output the row pairwise distance
		row_id2pairwise_dist_ls = snpData.calRowPairwiseDist(assumeBiAllelic=True, outputFname=outputFname, \
														hetHalfMatchDistance=hetHalfMatchDistance)
		cls.outputRow_id2pairwise_dist_lsInMatrix(row_id2pairwise_dist_ls, outputFname)
	
	"""
		#2011-4-7
		inputFname = os.path.expanduser('~/script/vervet/data/topConservedHG19_random_100loci_vervetSNP.tsv')
		
		inputFname = os.path.expanduser('~/mnt/hoffman2_home/script/vervet/data/1MbBAC_as_ref/454_illu_6_sub_vs_1MbBAC_minMAC3_maxMAC7_maxNoOfReadsForGenotypingError1.maxCoverage20.maxMajorAC10.maxNoOfReadsForAllSamples1000.snps')
		snpFnameToDoFiltering = os.path.expanduser('~/mnt/hoffman2_home/script/vervet/data/1MbBAC_as_ref/454_illu_6_sub_vs_1MbBAC_minMAC3_maxMAC7_maxNoOfReadsForGenotypingError1.maxCoverage20.maxMajorAC10.maxNoOfReadsForAllSamples1000_no_carribean.snps')
		snpFnameToDoFiltering = None
		convertHetero2NA = True
		max_NA_rate = 0.4
		min_MAF = 0
		outputFname ='%s.pairwiseDist.convertHetero2NA%s.minMAF%s.maxNA%s.tsv'%(os.path.splitext(inputFname)[0], \
												convertHetero2NA, min_MAF, max_NA_rate)
		VariantDiscovery.calculatePairwiseDistanceOutOfSNPXStrainMatrix(inputFname, outputFname, snpFnameToDoFiltering=snpFnameToDoFiltering,\
							convertHetero2NA=convertHetero2NA, min_MAF=min_MAF, max_NA_rate=max_NA_rate)
		sys.exit(0)
		
		#2011-7-8
		inputPrefix = os.path.expanduser("/usr/local/vervetData/ref/454/vs_MinSize200Scaffolds_by_bwasw")
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		referenceIDSet = set(['Contig0',])
		outputFname = os.path.expanduser("%s.Contig0.bam"%(inputPrefix))
		VariantDiscovery.filterAlignmentByReferenceIDs(inputFname, outputFname, referenceIDSet=referenceIDSet, readGroup='454_vs_Contig0')
		sys.exit(0)
	"""
	
	
	@classmethod
	def outputRow_id2pairwise_dist_lsInMatrix(cls, row_id2pairwise_dist_ls, outputFname):
		"""
		2011-5-14
			split out of calculatePairwiseDistanceOutOfSNPXStrainMatrix()
			it appends to outputFname, rather than overwriting.
		"""
		#2011-4-7 output the pairwise distance as matrix
		import csv, numpy
		no_of_rows = len(row_id2pairwise_dist_ls)
		row_id_ls = row_id2pairwise_dist_ls.keys()
		row_id_ls.sort()
		row_id2index = {}
		for row_id in row_id_ls:
			row_id2index[row_id] = len(row_id2index)
		
		data_matrix = numpy.zeros([no_of_rows, no_of_rows], dtype=numpy.float)
		writer = csv.writer(open(outputFname, 'a'), delimiter='\t')
		for row_id in row_id_ls:
			pairwise_dist_ls = row_id2pairwise_dist_ls.get(row_id)
			for dist in pairwise_dist_ls:
				mismatch_rate, row_id2, no_of_mismatches, no_of_non_NA_pairs = dist[:4]
				i = row_id2index.get(row_id)
				j = row_id2index.get(row_id2)
				data_matrix[i][j] = mismatch_rate
				data_matrix[j][i] = mismatch_rate
		header = [''] + row_id_ls
		writer.writerow(header)
		for i in range(no_of_rows):
			row_id = row_id_ls[i]
			data_row = [row_id] + list(data_matrix[i])
			writer.writerow(data_row)
		del writer
		
	
	def run(self,):
		if self.debug:	# 2010-4-18 enter debug mode "~/.../variation/misc.py -b"
			import pdb
			pdb.set_trace()
			debug = True
		else:
			debug =False
		
		
		snpFnameToDoFiltering = None
		"""
		convertHetero2NA = True
		max_NA_rate = 0.4
		min_MAF = 0
		outputFname ='%s.pairwiseDist.convertHetero2NA%s.minMAF%s.maxNA%s.tsv'%(os.path.splitext(self.inputFname)[0], \
												convertHetero2NA, min_MAF, max_NA_rate)
		"""
		self.calculatePairwiseDistanceOutOfSNPXStrainMatrix(self.inputFname, self.outputFname, \
							snpFnameToDoFiltering=snpFnameToDoFiltering,\
							convertHetero2NA=self.convertHetero2NA, min_MAF=self.min_MAF, max_NA_rate=self.max_NA_rate,\
							hetHalfMatchDistance=self.hetHalfMatchDistance)
		sys.exit(0)
		
	
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = CalculatePairwiseDistanceOutOfSNPXStrainMatrix
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
