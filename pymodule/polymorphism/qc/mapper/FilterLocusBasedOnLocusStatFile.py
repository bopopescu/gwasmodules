#!/usr/bin/env python
"""
Examples:
	%s -i folderRun/Scaffold301_splitVCF_unit1.vcf -o folderRun/Scaffold301_splitVCF_unit1.info.vcf
		--liftOverLocusMapPvalueFname ...probability.tsv
	%s 
	
	%s 
	
Description:
	2014.01.04 program that filters loci based on stats from statFname. 

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, MatrixFile, PassingData
from pymodule.yhio.VCFFile import VCFFile
from pymodule.pegasus.mapper.AbstractVCFMapper import AbstractVCFMapper

parentClass = AbstractVCFMapper
class FilterLocusBasedOnLocusStatFile(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update({
						('runType', 1, int): ['', '', 1, ' runType 1: locus missing fraction, header is (locusID, occurrence_byFixedValue).', ],\
						('statFname', 1, ): ['', '', 1, ' locus statistics file, usually one locus per line, depending on runType.', ],\
						('minValue', 0, float): [None, '', 1, 'if not None, loci whose stat below this would be removed.', ],\
						('maxValue', 0, float): [None, '', 1, 'if not None, loci whose stat above this would be removed.', ],\
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		parentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		
		self.getLocusID2StatFunctionDict = {1: self.getLocusID2MissingFraction,}
	
	def getLocusID2MissingFraction(self, inputFname=None):
		"""
		2014.01.08
			
		"""
		sys.stderr.write("Reading in the missing statistics from %s ... "%(inputFname))
		locusID2Stat = {}
		
		reader = MatrixFile(inputFname=inputFname)
		reader.constructColName2IndexFromHeader()
		locusIDIndex = reader.getColIndexGivenColHeader("locusID")
		statIndex = reader.getColIndexGivenColHeader("occurrence_byFixedValue")
		counter = 0
		for row in reader:
			locusID = row[locusIDIndex]
			chromosome, start = locusID.split('_')[:2]
			start = int(start)
			stat = float(row[statIndex])
			
			key = (chromosome, start, start)
			if key in locusID2Stat:
				if stat < locusID2Stat[key]:
					#take lowest value
					locusID2Stat[key] = stat
			else:
				locusID2Stat[key] = stat
			counter += 1
		del reader
		sys.stderr.write(" %s unique loci with missing fraction out of %s total loci.\n"%(len(locusID2Stat), counter))
		return locusID2Stat
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		
		outputDir = os.path.split(self.outputFname)[0]
		if outputDir and not os.path.isdir(outputDir):
			os.makedirs(outputDir)
		locusID2Stat = self.getLocusID2StatFunctionDict[self.runType](self.statFname)
		
		reader = VCFFile(inputFname=self.inputFname)
		writer = VCFFile(outputFname=self.outputFname, openMode='w')
		writer.metaInfoLs = reader.metaInfoLs
		writer.header = reader.header
		writer.writeMetaAndHeader()
		
		counter = 0
		real_counter = 0
		
		for vcfRecord in reader:	#assuming input VCF is sorted
			counter += 1
			key = (vcfRecord.chromosome, vcfRecord.position, vcfRecord.position)
			stat = locusID2Stat.get(key)
			if stat is None:
				continue
			
			toKeepLocus = True
			if self.minValue is not None and stat < self.minValue:
				toKeepLocus = False
			if self.maxValue is not None and stat > self.maxValue:
				toKeepLocus = False
			
			if toKeepLocus:
				real_counter += 1
				writer.writeVCFRecord(vcfRecord)
		reader.close()
		writer.close()
		if counter>0:
			fraction = real_counter/float(counter)
		else:
			fraction = -1
		sys.stderr.write("%s out of %s records, or %s, retained.\n"%(real_counter, counter, \
												fraction))

if __name__ == '__main__':
	main_class = FilterLocusBasedOnLocusStatFile
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()