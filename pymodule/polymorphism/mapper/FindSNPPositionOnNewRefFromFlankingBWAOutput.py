#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s --newSNPDataOutputFname ~/script/vervet/data/194SNPData/isq524CoordinateSNPData_max15Mismatch.tsv
		--querySNPDataFname ~/script/vervet/data/194SNPData/AllSNPData.txt
		--SNPFlankSequenceFname ~/script/vervet/data/194SNPData/AllSNPFlankWithSNPMark.txt
		-i bwa194SNPFlankAgainst524_15Mismatches.bam
		-o ~/script/vervet/data/194SNPData/originalSNPID2ISQ524Coordinate_max15Mismatch.tsv
		--maxNoOfMismatches 2
		--minAlignmentSpan 10
		

Description:
	2012.10 child of FindSNPPositionOnNewRefFromFlankingBlastOutput.
		Its main input is bwa alingment output (.bam file),
		This program is used as part of workflow inside FindNewRefCoordinatesGivenVCFFolderWorkflow.py.
		
		
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
							querySNPDataFname=None,\
							querySNPID2NewRefCoordinateOutputFname=None, newSNPDataOutputFname=None, \
							minAlignmentSpan=10, **keywords):
		"""
		2013.07.05 add queryChromosome, queryStart, queryStop, newRefBase to output
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
			querySNPID2attributes = self.getQuerySNPID2attributes(SNPFlankSequenceFname=SNPFlankSequenceFname)
		else:
			querySNPID2attributes = None
		
		sys.stderr.write("Finding new reference coordinates for SNPs from bwa alignment output %s ... \n"%(bwaOutputFname))
		samfile = pysam.Samfile(bwaOutputFname, "rb" )
		
		counter = 0
		real_counter = 0
		queryIDSet= set()
		querySNPID2NewReferenceCoordinateLs = {}
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
				queryEnd = read.qend	#qend is exclusive 0-based. same value as inclusive 1-based
				
				targetChr = samfile.getrname(read.tid)
				targetStart = read.pos + 1	#pos is 0-based leftmost coordinate
				targetStop = read.aend	#aend is 0-based, points to one past the last aligned residue. Returns None if not available.
				
				
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
				queryStrand = "-"
				#reverse the query coordinates (pysam stores the two coordinates in ascending order regardless of strand).
				#reverse them for code below
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
				queryStrand = "+"
				if targetStart > targetStop:
					sys.stderr.write("Error: aligned to the negative strand, but targetStart (%s) > targetStop (%s).\n"%\
									(targetStart, targetStop))
					sys.exit(3)
			
				
			if yhRead.getNoOfIndels()==0 and queryAlignmentSpan == targetAlignmentSpan and \
					(self.minNoOfIdentities is None or yhRead.getNoOfMatches()>=self.minNoOfIdentities) and \
					(self.maxNoOfMismatches is None or yhRead.getNoOfMismatches()<=self.maxNoOfMismatches) and\
					(self.minIdentityFraction is None or yhRead.getMatchFraction()>=self.minIdentityFraction) and\
					targetAlignmentSpan>=minAlignmentSpan and queryAlignmentSpan>=minAlignmentSpan:
				if querySNPID2attributes and queryID in querySNPID2attributes:
					parseData = querySNPID2attributes.get(queryID)
					locusSpan = parseData.locusSpan	#locusSpan is the length of the locus itself, excluding the flanks
					#SNP's locusSpan is regarded as 0.
				else:
					parseData = self.parseQueryLocusID(queryID)
					chromosome = parseData.chromosome
					start = parseData.start
					stop = parseData.stop
					if start is not None and stop is not None:
						stop = int(stop)
						start = int(start)
						locusSpan = abs(stop-start)	#=length-1
					else:
						locusSpan = None
				positionInFlank = parseData.positionInFlank
				queryRefBase = parseData.refBase
				queryAltBase = parseData.altBase
				if positionInFlank is not None and locusSpan is not None:
					positionInFlank = int(positionInFlank)
					if queryStart <queryEnd and positionInFlank>queryStart and positionInFlank<queryEnd:
						#locus must be in the middle of queryStart and queryEnd.
						newRefStart =  targetStart + (positionInFlank - queryStart)
						newRefStop =  targetStop - (queryEnd - positionInFlank-locusSpan)
						queryStrand = "+"
						#query alignment start/stop are always in ascending order, regardless of strand
						queryAlignmentStart = max(1, parseData.start - (positionInFlank-1) + (queryStart-1))
						queryAlignmentStop = queryAlignmentStart + targetAlignmentSpan-1
						
					elif queryStart >queryEnd and positionInFlank<queryStart and positionInFlank>queryEnd:
						#could happen. on the opposite strand. targetStart is always bigger than targetStop
						#locus must be in the middle of queryStart and queryEnd.
						newRefStart=  targetStop - (positionInFlank-queryEnd)
						newRefStop =  targetStart + (queryStart - positionInFlank-locusSpan)
						queryStrand = "-"
						#query alignment start/stop are always in ascending order, regardless of strand
						queryAlignmentStart = max(1, parseData.start - (positionInFlank-1) + (queryEnd-1))
						queryAlignmentStop = queryAlignmentStart + targetAlignmentSpan-1
					
					else:
						newRefStart = None
						newRefStop = None
					if newRefStart is not None and newRefStop is not None:
						if queryID not in querySNPID2NewReferenceCoordinateLs:
							querySNPID2NewReferenceCoordinateLs[queryID] = []
						#all coordinates are 1-based, inclusive.
						newRefCoordinate = PassingData(newChr=targetChr, newRefStart=newRefStart, newRefStop=newRefStop, \
												newRefBase=yhRead.read.seq[positionInFlank-1], \
												targetAlignmentSpan=targetAlignmentSpan,\
												targetAlignmentStart=targetStart,\
												targetAlignmentStop=targetStop,\
												queryStrand=queryStrand, \
												queryAlignmentSpan=queryAlignmentSpan,\
												queryAlignmentStart=queryAlignmentStart,\
												queryAlignmentStop=queryAlignmentStop,\
												queryChromosome=parseData.chromosome, \
												queryStart=parseData.start, queryStop=parseData.stop,\
												queryRefBase=queryRefBase, queryAltBase=queryAltBase )
						querySNPID2NewReferenceCoordinateLs[queryID].append(newRefCoordinate)
						real_counter += 1
						
		sys.stderr.write(" from %s reads, no_of_hits_with_exception=%s, no_of_reads_mapped=%s, %s/%s SNPs found new-reference coordinates.\n"%\
						(counter, no_of_hits_with_exception, no_of_reads_mapped, real_counter, len(queryIDSet)))
		
		if querySNPDataFname and newSNPDataOutputFname:
			self.outputSNPDataInNewCoordinate(querySNPDataFname=querySNPDataFname, \
									querySNPID2NewReferenceCoordinateLs=querySNPID2NewReferenceCoordinateLs, newSNPDataOutputFname=newSNPDataOutputFname)
		return querySNPID2NewReferenceCoordinateLs
	

if __name__ == '__main__':
	main_class = FindSNPPositionOnNewRefFromFlankingBWAOutput
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()