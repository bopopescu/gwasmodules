#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s

Description:
	2013.05.03 a child class of BeagleGenotypeFile. used to describe Beagle likelihood file format, which looks like:
		marker alleleA alleleB 1000_709_1996093_GA_vs_524 1000_709_1996093_GA_vs_524 1000_709_1996093_GA_vs_524 1001_710_1995025_GA_vs_524 1001_710_1995025_GA_vs_524 1001_710_1995025_GA_vs_524 1002_711_2001039_GA_vs_524
		Contig791:1086 C A 1 0 0 0.9997 0.0003 0 0
		Contig791:1649 G C 1 0 0 1 0 0 0
		Contig791:4084 A C 1 0 0 1 0 0 0
		Contig791:4118 A G 1 0 0 1 0 0 0
		Contig791:4143 C A 1 0 0 1 0 0 0
		Contig791:4168 G C 1 0 0 0.9999 0.0001 0 0
		Contig791:4203 C G 1 0 0 1 0 0 0

	Example:
	
		reader = MatrixFile(inputFname='/tmp/input.txt', openMode='r')
		reader = MatrixFile('/tmp/input.txt', openMode='r')
		reader.constructColName2IndexFromHeader()
		for row in reader:
			row[reader.getColName2IndexFromHeader('KID')]
		
		inf = utils.openGzipFile(inputFname, openMode='r')
		reader = MatrixFile(inputFile=inf)
		
		#2013.2.1 writing
		writer = MatrixFile('/tmp/output.txt', openMode='w', delimiter='\t')
		writer.writeHeader(...)
		writer.writerow(row)
		writer.close()
	
	

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from pymodule.ProcessOptions import  ProcessOptions
from pymodule import utils, PassingData
from MatrixFile import MatrixFile
from BeagleGenotypeFile import BeagleGenotypeFile

parentClass = BeagleGenotypeFile
class BeagleLikelihoodFile(parentClass):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(parentClass.option_default_dict)
	option_default_dict.update({
						#('delimiter', 0, ): [' ', '', 1, 'delimiter for Beagle likelihood format is single-space'],\
						})
	def __init__(self, inputFname=None, **keywords):
		parentClass.__init__(self, inputFname=inputFname, **keywords)
	
	def constructColName2IndexFromHeader(self):
		"""
		2013.05.03
			First three column is for marker, alleleA, alleleB.
			one sample (diploid) occupies three columns.
			
		marker alleleA alleleB 1000_709_1996093_GA_vs_524 1000_709_1996093_GA_vs_524 1000_709_1996093_GA_vs_524 1001_710_1995025_GA_vs_524 1001_710_1995025_GA_vs_524 1001_710_1995025_GA_vs_524 1002_711_2001039_GA_vs_524
		Contig791:1086 C A 1 0 0 0.9997 0.0003 0 0
		Contig791:1649 G C 1 0 0 1 0 0 0
		Contig791:4084 A C 1 0 0 1 0 0 0
		"""
		self.header = self.next().genotypeLikelihoodList
		self.col_name2index = {}
		for i in xrange(len(self.header)):
			sampleID = self.header[i]
			if sampleID not in self.col_name2index:
				self.col_name2index[sampleID] = []
			self.col_name2index[sampleID].append(i)	#the index corresponds to genotypeLikelihoodList in next().
				#so it starts from 0
		return self.col_name2index
	
	def getLikelihoodListOfOneGenotypeOneSample(self, oneLocus=None, sampleID=None):
		"""
		2013.05.06
			oneLocus is output of next()
		"""
		try:
			sampleStartIndex = self.getColIndexGivenColHeader(sampleID)[0]
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
			sys.stderr.write("  sampleID %s not in beagleLikelihoodFile.\n"%(sampleID))
			raise
		tripleLikelihood = oneLocus.genotypeLikelihoodList[sampleStartIndex:sampleStartIndex+3]
		return tripleLikelihood
	
	def next(self):
		try:
			row = self.csvFile.next()
		except:
			raise StopIteration
		if not self.isRealCSV:
			row = row.strip().split()
		markerID, alleleA, alleleB = row[0:3]
		return PassingData(markerID=markerID, alleleA=alleleA, alleleB=alleleB, genotypeLikelihoodList=row[3:])

if __name__ == '__main__':
	main_class = BeagleLikelihoodFile
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()