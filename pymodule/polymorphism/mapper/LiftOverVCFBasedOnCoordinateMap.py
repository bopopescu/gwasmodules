#!/usr/bin/env python
"""
Examples:
	%s -i folderRun/Scaffold301_splitVCF_unit1.vcf -o folderRun/Scaffold301_splitVCF_unit1.info.vcf
	
	%s 
	
	%s 
	
Description:
	2013.07.10 program that does lift over VCF variants based on output from FindSNPPositionOnNewRefFromFlanking...Output.py
		assuming input VCF is sorted

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
class LiftOverVCFBasedOnCoordinateMap(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update({
						('coordinateMapFname', 1, ): ['', '', 1, 'file that has a map between old and new coordinates. output of FindSNPPositionOnNewRefFromFlankingBlastOutput.py', ],\
						
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		parentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	
	def readInCoordinateMap(self, coordinateMapFname=None):
		"""
		2013.07.11
			querySNPID      queryStrand     queryChromosome queryStart      queryStop       queryRefBase    queryAltBase    queryAlignmentSpan
			queryAlignmentStart     queryAlignmentStop      newChr  newRefStart     newRefStop      newRefBase      targetAlignmentSpan
			targetAlignmentStart    targetAlignmentStop
		"""
		sys.stderr.write("Reading in the coordinate map from %s ..."%(coordinateMapFname))
		oldCoordinate2newCoordinateDataLs = {}
		reader = MatrixFile(inputFname=coordinateMapFname)
		reader.constructColName2IndexFromHeader()
		oldChromosomeIndex = reader.getColIndexGivenColHeader("queryChromosome")
		oldStartIndex = reader.getColIndexGivenColHeader("queryStart")
		strandIndex = reader.getColIndexGivenColHeader("queryStrand")
		oldRefBaseIndex = reader.getColIndexGivenColHeader("queryRefBase")
		oldAltBaseIndex = reader.getColIndexGivenColHeader("queryAltBase")
		
		newChromosomeIndex = reader.getColIndexGivenColHeader("newChr")
		newStartIndex = reader.getColIndexGivenColHeader("newRefStart")
		newStopIndex = reader.getColIndexGivenColHeader("newRefStop")
		newRefBaseIndex = reader.getColIndexGivenColHeader("newRefBase")
		counter = 0
		for row in reader:
			oldChromosome = row[oldChromosomeIndex]
			oldStart = int(row[oldStartIndex])
			strand = row[strandIndex]
			oldRefBase = row[oldRefBaseIndex]
			oldAltBase = row[oldAltBaseIndex]
			
			newChromosome = row[newChromosomeIndex]
			newStart = int(row[newStartIndex])
			newStop = int(row[newStopIndex])
			newRefBase = row[newRefBaseIndex]
			
			key = (oldChromosome, oldStart)
			if key not in oldCoordinate2newCoordinateDataLs:
				oldCoordinate2newCoordinateDataLs[key] = []
			oldCoordinate2newCoordinateDataLs[key].append(PassingData(strand=strand, oldRefBase=oldRefBase, \
												oldAltBase=oldAltBase, newChromosome=newChromosome, newStart=newStart,\
												newStop=newStop, newRefBase=newRefBase))
			counter += 1
		del reader
		sys.stderr.write("%s old coordinates with %s new coordinates.\n"%(len(oldCoordinate2newCoordinateDataLs), counter))
		return oldCoordinate2newCoordinateDataLs
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		
		outputDir = os.path.split(self.outputFname)[0]
		if outputDir and not os.path.isdir(outputDir):
			os.makedirs(outputDir)
		oldCoordinate2newCoordinateDataLs = self.readInCoordinateMap(self.coordinateMapFname)
		
		
		self.reader = VCFFile(inputFname=self.inputFname)
		
		self.writer = VCFFile(outputFname=self.outputFname, openMode='w')
		self.writer.metaInfoLs = self.reader.metaInfoLs
		self.writer.header = self.reader.header
		self.writer.writeMetaAndHeader()
		
		
		counter = 0
		real_counter = 0
		noOfRecordsWithMultiNewCoords = 0
		
		for vcfRecord in self.reader:	#assuming input VCF is sorted
			counter += 1
			key = (vcfRecord.chromosome, vcfRecord.position)
			newCoordinateDataLs = oldCoordinate2newCoordinateDataLs.get(key)
			if newCoordinateDataLs is None:
				continue
			if len(newCoordinateDataLs)>1:
				noOfRecordsWithMultiNewCoords += 1
				continue
			newCoordinateData = newCoordinateDataLs[0]
			vcfRecord.setChromosome(newCoordinateData.newChromosome)
			vcfRecord.setPosition(newCoordinateData.newStart)
			if newCoordinateData.strand=='-':
				newRefBase = Seq(newCoordinateData.oldRefBase).reverse_complement()
				newAltBase = Seq(newCoordinateData.oldAltBase).reverse_complement()
			else:
				newRefBase = newCoordinateData.oldRefBase
				newAltBase = newCoordinateData.oldAltBase
			
			vcfRecord.setRefAllele(newRefBase)
			vcfRecord.setAltAllele(newAltBase)
			real_counter += 1
			self.writer.writeVCFRecord(vcfRecord)
			
			
		
		self.reader.close()
		self.writer.close()
		sys.stderr.write("%s (out of %s, %s) records found new coordinates. %s records with >1 new coordinates, discarded.\n"%(real_counter, counter, \
												real_counter/float(counter), noOfRecordsWithMultiNewCoords))

if __name__ == '__main__':
	main_class = LiftOverVCFBasedOnCoordinateMap
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()