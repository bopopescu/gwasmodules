#!/usr/bin/env python
"""
Examples:
	%s -i  -o 
	
	%s 
	
	%s 
	
Description:
	2013.07.15 program that computes how many genomes/loci are left at particular switch frequency threshold if  
		chromosomes/contigs with higher switch frequency were all thrown away.
		Input is the switchPointFname output of FindSNPPositionOnNewRefFromFlankingBlastOutput.py

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, MatrixFile, PassingData
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper

parentClass = AbstractMapper
class CalculateLociAndGenomeCoveredAtEachSwitchFrequencyThreshold(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update({
				})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		parentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def readInStats(self, inputFname=None):
		"""
		2013.07.15
		"""
		sys.stderr.write("Reading stats from %s ..."%(inputFname))
		
		data_matrix = []
		
		reader = MatrixFile(inputFname)
		reader.constructColName2IndexFromHeader()
		switchFrequencyIndex = reader.getColIndexGivenColHeader("noOfSwitchPoints_by_noOfLociWithUniqueHit")
		regionSpanIndex = reader.getColIndexGivenColHeader("regionSpan")
		noOfLociIndex = reader.getColIndexGivenColHeader("#sitesInInput2")
		
		totalSpan = 0
		totalNoOfLoci = 0
		counter = 0
		for row in reader:
			counter += 1
			switchFrequency = row[switchFrequencyIndex]
			regionSpan = row[regionSpanIndex]
			noOfLoci = row[noOfLociIndex]
			if switchFrequency and regionSpan and noOfLoci:	#non-empty
				switchFrequency = float(switchFrequency)
				regionSpan = int(float(regionSpan))
				noOfLoci = int(float(noOfLoci))
				data_matrix.append([switchFrequency, regionSpan, noOfLoci])
				totalSpan += regionSpan
				totalNoOfLoci += noOfLoci
		reader.close()
		sys.stderr.write(" %s valid entries (from %s rows) with totalSpan=%s, totalNoOfLoci=%s.\n"%\
						(len(data_matrix), counter, totalSpan, totalNoOfLoci))
		return PassingData(data_matrix=data_matrix, totalSpan=totalSpan, totalNoOfLoci=totalNoOfLoci)
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		
		outputDir = os.path.split(self.outputFname)[0]
		if outputDir and not os.path.isdir(outputDir):
			os.makedirs(outputDir)
		
		
		switchPointData = self.readInStats(inputFname=self.inputFname)
		
		sys.stderr.write("Processing data ...")
		writer = MatrixFile(self.outputFname, openMode='w')
		header = ["maxSwitchFrequency", "genomeCovered", 'genomeCoveredFraction', "noOfLoci", 'noOfLociFraction']
		writer.writeHeader(header)
		
		data_matrix = switchPointData.data_matrix
		totalSpan = switchPointData.totalSpan
		totalNoOfLoci = switchPointData.totalNoOfLoci
		
		#sort it based on switchFrequency
		data_matrix.sort(reverse=True)
		maxSwitchFrequencyLs = []
		cumulativeRegionSpanLs = []
		cumulativeNoOfLociLs = []
		for i in xrange(len(data_matrix)):
			switchFrequency, regionSpan, noOfLoci = data_matrix[i]
			maxSwitchFrequencyLs.append(switchFrequency)
			if i==0:
				cumulativeRegionSpan = totalSpan-regionSpan
				
				cumulativeNoOfLoci = totalNoOfLoci - noOfLoci
			else:
				cumulativeRegionSpan = cumulativeRegionSpanLs[i-1]-regionSpan
				cumulativeNoOfLoci = cumulativeNoOfLociLs[i-1] - noOfLoci
			cumulativeRegionSpanLs.append(cumulativeRegionSpan)
			cumulativeNoOfLociLs.append(cumulativeNoOfLoci)
			writer.writerow([switchFrequency, cumulativeRegionSpan, cumulativeRegionSpan/float(totalSpan),\
							cumulativeNoOfLoci, cumulativeNoOfLoci/float(totalNoOfLoci)])
		writer.close()
		sys.stderr.write(".\n")

if __name__ == '__main__':
	main_class = CalculateLociAndGenomeCoveredAtEachSwitchFrequencyThreshold
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()