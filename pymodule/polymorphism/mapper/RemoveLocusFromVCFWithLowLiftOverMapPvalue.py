#!/usr/bin/env python
"""
Examples:
	%s -i folderRun/Scaffold301_splitVCF_unit1.vcf -o folderRun/Scaffold301_splitVCF_unit1.info.vcf
		--liftOverLocusMapPvalueFname ...probability.tsv
	%s 
	
	%s 
	
Description:
	2014.01.04 program that removes loci with low mapPvalue (output of ComputeLiftOverLocusProbability.py) 
		from a lifted-over VCF (LiftOverVCFBasedOnCoordinateMap.py).

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from Bio.Seq import Seq
from pymodule import ProcessOptions, MatrixFile, PassingData
from pymodule.yhio.VCFFile import VCFFile
from pymodule.pegasus.mapper.AbstractVCFMapper import AbstractVCFMapper

parentClass = AbstractVCFMapper
class RemoveLocusFromVCFWithLowLiftOverMapPvalue(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update({
						('liftOverLocusMapPvalueFname', 1, ): ['', '', 1, ' output of ComputeLiftOverLocusProbability.py', ],\
						('minLiftOverMapPvalue', 1, float): [0.5, '', 1, 'locus with mapPvalue lower than this would be removed.', ],\
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		parentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	
	def getLocusNewID2mapPvalue(self, liftOverLocusMapPvalueFname=None):
		"""
		2014.01.04
			oldChromosome, oldStart, oldStop, oldStrand, newChromosome, newStart, newStop, mapPvalue
		"""
		sys.stderr.write("Reading in the coordinate map from %s ..."%(liftOverLocusMapPvalueFname))
		locusNewID2mapPvalue = {}
		reader = MatrixFile(inputFname=liftOverLocusMapPvalueFname)
		reader.constructColName2IndexFromHeader()
		strandIndex = reader.getColIndexGivenColHeader("oldStrand")
		newChromosomeIndex = reader.getColIndexGivenColHeader("newChromosome")
		newStartIndex = reader.getColIndexGivenColHeader("newStart")
		newStopIndex = reader.getColIndexGivenColHeader("newStop")
		mapPvalueIndex = reader.getColIndexGivenColHeader("mapPvalue")
		counter = 0
		for row in reader:
			strand = row[strandIndex]
			newChromosome = row[newChromosomeIndex]
			newStart = int(row[newStartIndex])
			newStop = int(row[newStopIndex])
			mapPvalue = float(row[mapPvalueIndex])
			
			key = (newChromosome, newStart, newStop)
			if key in locusNewID2mapPvalue:
				if mapPvalue < locusNewID2mapPvalue[key]:
					#take lowest value
					locusNewID2mapPvalue[key] = mapPvalue
			else:
				locusNewID2mapPvalue[key] = mapPvalue
			counter += 1
		del reader
		sys.stderr.write("%s unique loci with map p-value out of %s total loci.\n"%(len(locusNewID2mapPvalue), counter))
		return locusNewID2mapPvalue
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		
		outputDir = os.path.split(self.outputFname)[0]
		if outputDir and not os.path.isdir(outputDir):
			os.makedirs(outputDir)
		locusNewID2mapPvalue = self.getLocusNewID2mapPvalue(self.liftOverLocusMapPvalueFname)
		
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
			mapPvalue = locusNewID2mapPvalue.get(key)
			if mapPvalue is None:
				continue
			
			if mapPvalue >self.minLiftOverMapPvalue:
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
	main_class = RemoveLocusFromVCFWithLowLiftOverMapPvalue
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()