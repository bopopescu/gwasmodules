#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s --newSNPDataOutputFname ~/script/vervet/data/194SNPData/isq524CoordinateSNPData_max15Mismatch.tsv
		--originalSNPDataFname ~/script/vervet/data/194SNPData/AllSNPData.txt
		--SNPFlankSequenceFname ~/script/vervet/data/194SNPData/AllSNPFlankWithSNPMark.txt
		-i Blast/Blast194SNPFlankAgainst524_15Mismatches.2012.8.17T2334/folderBlast/blast.tsv
		-o ~/script/vervet/data/194SNPData/originalSNPID2ISQ524Coordinate_max15Mismatch.tsv
		--maxNoOfMismatches 2
		--minAlignmentSpan 10
		

Description:
	2012.10
		inputFname is FindNewRefCoordinatesGivenVCFFolderWorkflow.py's bwa-alignment output and others.
		outputFname is the map between old and new coordinates.
			originalSNPID, newChr, newRefStart, newRefStop, strand, targetAlignmentSpan.
		
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
import pysam
import numpy, re
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, SNP
from pymodule.pegasus.mapper.extractor.ExtractFlankingSequenceForVCFLoci import ExtractFlankingSequenceForVCFLoci
from pymodule import figureOutDelimiter, getColName2IndexFromHeader
from pymodule.yhio.BamFile import YHAlignedRead
from FindSNPPositionOnNewRefFromFlankingBlastOutput import FindSNPPositionOnNewRefFromFlankingBlastOutput

class FindSNPPositionOnNewRefFromFlankingBWAOutput(FindSNPPositionOnNewRefFromFlankingBlastOutput):
	__doc__ = __doc__
	option_default_dict = FindSNPPositionOnNewRefFromFlankingBlastOutput.option_default_dict.copy()
	#option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		2012.8.19
		"""
		FindSNPPositionOnNewRefFromFlankingBlastOutput.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	
	def findSNPPositionOnNewRef(self, SNPFlankSequenceFname=None, blastHitResultFname=None, bwaOutputFname=None, \
							originalSNPDataFname=None,\
							originalSNPID2NewRefCoordinateOutputFname=None, newSNPDataOutputFname=None, \
							minAlignmentSpan=10, **keywords):
		"""
		#2013.05.21 bugfix, read.qend is occasionally not available
		2012.10.14 
		2012.10.8
			argument minAlignmentSpan: the number of bases involved in the blast query-target alignment
		2012.8.19
			newSNPDataOutputFname will contain the individual X SNP matrix.
		"""
		if bwaOutputFname is None:
			bwaOutputFname = blastHitResultFname
		
		if SNPFlankSequenceFname:
			originalSNPID2positionInFlank = self.getOriginalSNPID2positionInFlank(SNPFlankSequenceFname=SNPFlankSequenceFname)
		else:
			originalSNPID2positionInFlank = None
		
		sys.stderr.write("Finding new reference coordinates for SNPs from bwa alignment output %s ... \n"%(bwaOutputFname))
		samfile = pysam.Samfile(bwaOutputFname, "rb" )
		
		counter = 0
		real_counter = 0
		queryIDSet= set()
		originalSNPID2BlastRefCoordinateLs = {}
		no_of_reads_mapped = samfile.mapped	#: not good, segmentation fault because bai file is missing
		no_of_hits_with_exception = 0
		for read in samfile:	#.fetch():
			counter += 1
			if read.is_unmapped:
				continue
			#read.mapq
			queryID = None
			try:
				yhRead = YHAlignedRead(read)
				queryID = read.qname
				queryStart = read.qstart + 1	#qstart is 0-based
				queryEnd = read.qend	#qend is 0-based but exclusive. same as 1-based but inclusive.
				
				targetChr = samfile.getrname(read.tid)
				targetStart = read.pos + 1	#pos is 0-based
				targetStop = read.aend	#aend is 0-based but exclusive. same as 1-based but inclusive.
				
				
				queryAlignmentSpan = read.qlen
				targetAlignmentSpan = read.alen
			except:	#2013.05.21 bugfix, read.qend is occasionally not available
				sys.stderr.write('Except type for query %s : %s\n'%(queryID, repr(sys.exc_info())))
				import traceback
				traceback.print_exc()
				no_of_hits_with_exception += 1
				continue
			
			queryIDSet.add(queryID)
			if read.is_reverse:
				strand = -1
				#reverse the query coordinates (pysam stores the two coordinates in ascending order regardless of strand).
				tmp = queryStart
				queryStart = queryEnd
				queryEnd = tmp
				"""
				#check whether queryStart<queryEnd or targetStart <targetStop
				if targetStart < targetStop:
					sys.stderr.write("Error: aligned to the negative strand, but targetStart (%s) < targetStop (%s).\n"%\
									(targetStart, targetStop))
					sys.exit(3)
				"""
			else:
				strand = +1
				if targetStart > targetStop:
					sys.stderr.write("Error: aligned to the negative strand, but targetStart (%s) > targetStop (%s).\n"%\
									(targetStart, targetStop))
					sys.exit(3)
			
				
			if yhRead.getNoOfIndels()==0 and queryAlignmentSpan == targetAlignmentSpan and \
					(self.minNoOfIdentities is None or yhRead.getNoOfMatches()>=self.minNoOfIdentities) and \
					(self.maxNoOfMismatches is None or yhRead.getNoOfMismatches()<=self.maxNoOfMismatches) and\
					(self.minIdentityFraction is None or yhRead.getMatchFraction()>=self.minIdentityFraction) and\
					targetAlignmentSpan>=minAlignmentSpan and queryAlignmentSpan>=minAlignmentSpan:
				if originalSNPID2positionInFlank and queryID in originalSNPID2positionInFlank:
					parseData = originalSNPID2positionInFlank.get(queryID)
					locusSpan = parseData.locusSpan	#locusSpan is the length of the locus itself, excluding the flanks
					#SNP's locusSpan is regarded as 0.
				else:
					parseData = self.parseOriginalLocusID(queryID)
					chr = parseData.chr
					start = parseData.start
					stop = parseData.stop
					if start is not None and stop is not None:
						stop = int(stop)
						start = int(start)
						locusSpan = abs(stop-start)	#=length-1
					else:
						locusSpan = None
				positionInFlank = parseData.positionInFlank
				originalRefBase = parseData.refBase
				originalAltBase = parseData.altBase
				if positionInFlank is not None and locusSpan is not None:
					positionInFlank = int(positionInFlank)
					if queryStart >queryEnd and positionInFlank<queryStart and positionInFlank>queryEnd:
						#could happen. on the opposite strand. targetStart is always bigger than targetstop
						#locus must be in the middle of queryStart and queryEnd.
						newRefStart=  targetStop - (positionInFlank-queryEnd)
						newRefStop =  targetStart + (queryStart - positionInFlank-locusSpan)
						strand = -1
					elif queryStart <queryEnd and positionInFlank>queryStart and positionInFlank<queryEnd:
						#locus must be in the middle of queryStart and queryEnd.
						newRefStart =  targetStart + (positionInFlank - queryStart)
						newRefStop =  targetStop - (queryEnd - positionInFlank-locusSpan)
						strand = + 1
					else:
						newRefStart = None
						newRefStop = None
					if newRefStart is not None and newRefStop is not None:
						if queryID not in originalSNPID2BlastRefCoordinateLs:
							originalSNPID2BlastRefCoordinateLs[queryID] = []
						newRefCoordinate = PassingData(newChr=targetChr, newRefStart=newRefStart, newRefStop=newRefStop, \
													strand=strand, targetAlignmentSpan=targetAlignmentSpan,\
													queryAlignmentSpan=queryAlignmentSpan,\
													originalRefBase=originalRefBase, originalAltBase=originalAltBase )
						originalSNPID2BlastRefCoordinateLs[queryID].append(newRefCoordinate)
						real_counter += 1
						
		sys.stderr.write(" from %s reads, no_of_hits_with_exception=%s, no_of_reads_mapped=%s, %s/%s SNPs found new-reference coordinates.\n"%\
						(counter, no_of_hits_with_exception, no_of_reads_mapped, real_counter, len(queryIDSet)))
		
		#output the mapping
		if originalSNPID2NewRefCoordinateOutputFname:
			self.outputOriginalSNPID2NewCoordinateMap(originalSNPID2BlastRefCoordinateLs=originalSNPID2BlastRefCoordinateLs, \
						originalSNPID2NewRefCoordinateOutputFname=originalSNPID2NewRefCoordinateOutputFname)
		
		if originalSNPDataFname and newSNPDataOutputFname:
			self.outputSNPDataInNewCoordinate(originalSNPDataFname=originalSNPDataFname, \
									originalSNPID2BlastRefCoordinateLs=originalSNPID2BlastRefCoordinateLs, newSNPDataOutputFname=newSNPDataOutputFname)
	

if __name__ == '__main__':
	main_class = FindSNPPositionOnNewRefFromFlankingBWAOutput
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()