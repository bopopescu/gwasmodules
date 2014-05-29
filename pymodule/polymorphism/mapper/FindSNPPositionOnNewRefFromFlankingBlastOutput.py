#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s --newSNPDataOutputFname ~/script/vervet/data/194SNPData/isq524CoordinateSNPData_max15Mismatch.tsv
		--querySNPDataFname ~/script/vervet/data/194SNPData/AllSNPData.txt
		--SNPFlankSequenceFname ~/script/vervet/data/194SNPData/AllSNPFlankWithSNPMark.txt
		-i Blast/Blast194SNPFlankAgainst524_15Mismatches.2012.8.17T2334/folderBlast/blast.tsv
		-o ~/script/vervet/data/194SNPData/originalSNPID2ISQ524Coordinate_max15Mismatch.tsv
		--maxNoOfMismatches 2
		--minAlignmentSpan 10 --newSNPDataOutputFormat 2
		

Description:
	2012.8.19
		given a fasta file of flanking sequences of polymorphic loci (SNP or SVs).
		find its new positions on the blast database (new reference).
		
		inputFname is BlastWorkflow.py's output: blast SNPFlankSequenceFname (converting [A/C] to A) to a new reference.
		outputFname is the map between old and new coordinates (NOT a SNP dataset, which is --newSNPDataOutputFname).
			format 1, a tsv file with this header:
				querySNPID, queryStrand, queryChromosome, queryStart, queryStop,
				queryRefBase, queryAltBase, queryAlignmentSpan, queryAlignmentStart, queryAlignmentStop,
				newChr, newRefStart, newRefStop, \
				newRefBase, targetAlignmentSpan, targetAlignmentStart, targetAlignmentStop
			format 2, a chain file, http://genome.ucsc.edu/goldenPath/help/chain.html.
				Header Lines
					chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
				Alignment Data Lines, three required attribute values:
					size dt dq 
				NOTE: The last line of the alignment section contains only one number: 
					the ungapped alignment size of the last block.
	
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, PassingData, utils, SNP, MatrixFile
from pymodule import figureOutDelimiter, getColName2IndexFromHeader
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
from pymodule.pegasus.mapper.extractor.ExtractFlankingSequenceForVCFLoci import ExtractFlankingSequenceForVCFLoci
import numpy, re


class FindSNPPositionOnNewRefFromFlankingBlastOutput(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	#option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
							('SNPFlankSequenceFname', 0, ): ['', 'S', 1, 'in fasta format, with SNP embedded as [A/C] in sequence.\n\
	if given, use this to fetch a map between querySNPID & positionInFlank.', ],\
							('querySNPDataFname', 0, ): ['', 'l', 1, 'path to an input tsv/csv, 3-column: Sample SNP Geno.', ],\
							('newSNPDataOutputFname', 0, ): [None, '', 1, 'if given, would be a Strain X Locus Yu format polymorphism file with new coordinates.', ],\
							('newSNPDataOutputFormat', 0, int): [1, '', 1, 'format 1: SNP_ID =chr_start_stop; format 2: SNP_ID = chr_start', ],\
							('minNoOfIdentities', 0, int): [None, 'm', 1, 'minimum number of identities between a query and target', ],\
							('maxNoOfMismatches', 0, int): [None, 'a', 1, 'minimum number of mismatches between a query and target', ],\
							('minIdentityFraction', 0, float): [None, 'n', 1, 'minimum fraction of identities between a query and target', ],\
							('minAlignmentSpan', 1, int): [10, '', 1, 'minimum number of bases of the query and target that are involved in the blast alignment', ],\
							('outputFileFormat', 1, int): [1, '', 1, 'output file format. 1: ', ],\
							('chainFilename', 0, ): [None, '', 1, 'output file to store chains (UCSC chain file format)', ],\
							('switchPointFname', 0, ): ['', '', 1, 'file to store switch points information', ],\
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		2012.8.19
		"""
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		
	
	def getQuerySNPID2attributes(self, SNPFlankSequenceFname=None):
		"""
		2013.07.01 bugfix
			snpPositionInFlank is renamed to positionInFlank
		2012.10.8
			split out of findSNPPositionOnNewRef()
			the positionInFlank is 1-based.
			add one more variable, locusSpan, in the value of querySNPID2attributes.
				=0 for SNPs, =length-1 for other loci.
			
			The SNPFlankSequenceFname is in fasta format, with SNP embedded as [A/C] in sequence.
		"""
		sys.stderr.write("Deriving querySNPID2attributes from %s ...\n"%(SNPFlankSequenceFname))
		querySNPID2attributes = {}
		from Bio import SeqIO
		inf = open(SNPFlankSequenceFname, "rU")
		counter = 0
		for record in SeqIO.parse(inf, "fasta"):
			querySNPID = record.id.split()[0]	#get rid of extra comment
			parseData = self.parseQueryLocusID(querySNPID)
			chromosome = parseData.chromosome
			start = parseData.start
			stop = parseData.stop
			
			positionInFlank = record.seq.find('[') + 1

			if positionInFlank<=0:
				sys.stderr.write("Error, could not find the relative position of snp %s in the query sequence \n"%(record.id))
				raise
			else:
				refBase = record.seq[positionInFlank]
				altBase = record.seq[positionInFlank+2]
				if parseData.refBase!=refBase:
					sys.stderr.write("Reference base inferred from snp ID (%s), %s, is different from the one, %s, inferred from its flanking sequence.\n"%\
									(record.id, parseData.refBase, refBase))
					raise
				if parseData.altBase!=altBase:
					sys.stderr.write("Reference base inferred from snp ID (%s), %s, is different from the one, %s, inferred from its flanking sequence.\n"%\
									(record.id, parseData.altBase, altBase))
					raise
				querySNPID2attributes[querySNPID] = PassingData(positionInFlank=positionInFlank, locusSpan=0,\
													refBase=refBase, altBase=altBase, \
													chromosome=chromosome, start=start, stop=stop)
				#
			counter += 1
		inf.close()
		sys.stderr.write(" %s/%s SNPs with positions.\n"%(len(querySNPID2attributes), counter))
		return querySNPID2attributes
	
	def parseQueryLocusID(self, locus_id=None):
		"""
		2012.10.8
			locus_id is in the format of '%s_%s_%s_positionInFlank%s'%(chromosome, start, stop, flankingLength+1)
			output of ExtractFlankingSequenceForVCFLoci.py
		"""
		search_result = ExtractFlankingSequenceForVCFLoci.sequenceTitlePattern.search(locus_id)
		chromosome = None
		start = None
		stop = None
		refBase = None
		altBase = None
		positionInFlank = None
		if search_result:
			chromosome = search_result.group(1)
			start = int(search_result.group(2))
			stop = int(search_result.group(3))
			refBase = search_result.group(4)
			altBase = search_result.group(5)
			positionInFlank = int(search_result.group(6))
			
		return PassingData(chromosome=chromosome, start=start, stop=stop, refBase=refBase, altBase=altBase, positionInFlank=positionInFlank)
		
	
	def findSNPPositionOnNewRef(self, SNPFlankSequenceFname=None, blastHitResultFname=None, \
							querySNPDataFname=None,\
							querySNPID2NewRefCoordinateOutputFname=None, newSNPDataOutputFname=None, minAlignmentSpan=10):
		"""
		2012.10.8
			argument minAlignmentSpan: the number of bases involved in the blast query-target alignment
		2012.8.19
			newSNPDataOutputFname will contain the individual X SNP matrix.
		"""
		if SNPFlankSequenceFname:
			querySNPID2attributes = self.getQuerySNPID2attributes(SNPFlankSequenceFname=SNPFlankSequenceFname)
		else:
			querySNPID2attributes = None
		
		sys.stderr.write("Finding blast reference coordinates for SNPs from %s ... \n"%(blastHitResultFname))
		reader = csv.reader(open(blastHitResultFname), delimiter='\t')
		header =reader.next()
		col_name2index = getColName2IndexFromHeader(header)
		
		#every coordinate in blastHitResultFname is 1-based.
		"""
queryID queryStart      queryEnd        queryLength     targetChr       targetStart     targetStop      targetLength    noOfIdentities  noOfMismatches  identityPercentage
34804_309       1       417     417     Contig293       2551654 2552070 3001801 413     4       0.9904076738609112
43608_166       1       574     574     Contig269       1565599 1566170 3181654 565     9       0.9843205574912892
44412_392       2       580     580     Contig269       1776095 1776673 3181654 577     3       0.9948275862068966

		"""
		queryIDIndex = col_name2index['queryID']
		queryStartIndex = col_name2index['queryStart']
		queryEndIndex = col_name2index['queryEnd']
		
		targetChrIndex = col_name2index['targetChr']
		targetStartIndex = col_name2index['targetStart']
		targetStopIndex = col_name2index['targetStop']
		querySNPID2NewReferenceCoordinateLs = {}
		counter = 0
		real_counter = 0
		queryIDSet= set()
		for row in reader:
			queryID = row[queryIDIndex].split()[0]	##get rid of extra comment
			queryStart = int(row[queryStartIndex])
			queryEnd = int(row[queryEndIndex])
			
			targetChr = row[targetChrIndex]
			targetStart = int(row[targetStartIndex])
			targetStop = int(row[targetStopIndex])
			
			queryIDSet.add(queryID)
			
			queryAlignmentSpan = abs(queryEnd-queryStart) + 1
			targetAlignmentSpan = abs(targetStop-targetStart) + 1
			if queryAlignmentSpan == targetAlignmentSpan:
				if querySNPID2attributes and queryID in querySNPID2attributes:
					parseData = querySNPID2attributes.get(queryID)
					locusSpan = parseData.locusSpan
				else:
					parseData = self.parseQueryLocusID(queryID)
					start = parseData.start
					stop = parseData.stop
					if start is not None and stop is not None:
						stop = int(stop)
						start = int(start)
						locusSpan = abs(int(stop)-start)	#length-1
					else:
						locusSpan = None
				positionInFlank = parseData.positionInFlank
				queryRefBase = parseData.refBase
				queryAltBase = parseData.altBase
				if positionInFlank is not None and locusSpan is not None:
					positionInFlank = int(positionInFlank)
					if targetAlignmentSpan>=minAlignmentSpan and queryAlignmentSpan>=minAlignmentSpan:
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
							
							newRefCoordinate = PassingData(newChr=targetChr, newRefStart=newRefStart, newRefStop=newRefStop, \
												queryStrand=queryStrand, newRefBase="", \
												targetAlignmentSpan=targetAlignmentSpan,\
												targetAlignmentStart=targetStart,\
												targetAlignmentStop=targetStop,\
												queryAlignmentSpan=queryAlignmentSpan,\
												queryAlignmentStart=queryAlignmentStart,\
												queryAlignmentStop=queryAlignmentStop,\
												queryChromosome=parseData.chromosome, \
												queryStart=parseData.start, queryStop=parseData.stop,\
												queryRefBase=queryRefBase, queryAltBase=queryAltBase )
							querySNPID2NewReferenceCoordinateLs[queryID].append(newRefCoordinate)
							real_counter += 1
						
			counter += 1
		sys.stderr.write(" from %s blast results. %s/%s SNPs found blast-reference coordinates.\n"%\
						(counter, real_counter, len(queryIDSet)))
		
		
		if querySNPDataFname and newSNPDataOutputFname:
			self.outputSNPDataInNewCoordinate(querySNPDataFname=querySNPDataFname, \
									querySNPID2NewReferenceCoordinateLs=querySNPID2NewReferenceCoordinateLs, \
									newSNPDataOutputFname=newSNPDataOutputFname, \
									newSNPDataOutputFormat=self.newSNPDataOutputFormat)
		return querySNPID2NewReferenceCoordinateLs
	
	
	def outputQuerySNPID2NewCoordinateMap(self, querySNPID2NewReferenceCoordinateLs=None, \
						querySNPID2NewRefCoordinateOutputFname=None):
		"""
		2012.10.14
			split out of findSNPPositionOnNewRef()
		"""
		no_of_loci_with_1_newRef = 0
		no_of_loci_with_1Plus_newRef = 0
		sys.stderr.write("Outputting %s pairs in querySNPID2NewReferenceCoordinateLs to %s ..."%\
						(len(querySNPID2NewReferenceCoordinateLs), querySNPID2NewRefCoordinateOutputFname))
		writer = csv.writer(open(querySNPID2NewRefCoordinateOutputFname, 'w'), delimiter='\t')
		header = ['querySNPID', 'queryStrand', "queryChromosome", "queryStart", "queryStop", \
				'queryRefBase', 'queryAltBase',\
				'queryAlignmentSpan', 'queryAlignmentStart', 'queryAlignmentStop', \
				'newChr', 'newRefStart', 'newRefStop',  \
				"newRefBase", 'targetAlignmentSpan', 'targetAlignmentStart', 'targetAlignmentStop']
		writer.writerow(header)
		for querySNPID, newRefCoordinateLs in querySNPID2NewReferenceCoordinateLs.iteritems():
			if len(newRefCoordinateLs)==1:
				newRefCoordinate = newRefCoordinateLs[0]
				data_row = [querySNPID, newRefCoordinate.queryStrand, \
						newRefCoordinate.queryChromosome, newRefCoordinate.queryStart, newRefCoordinate.queryStop, \
						newRefCoordinate.queryRefBase, newRefCoordinate.queryAltBase, \
						newRefCoordinate.queryAlignmentSpan, newRefCoordinate.queryAlignmentStart, newRefCoordinate.queryAlignmentStop, 
						newRefCoordinate.newChr, newRefCoordinate.newRefStart, newRefCoordinate.newRefStop, \
						newRefCoordinate.newRefBase,\
						newRefCoordinate.targetAlignmentSpan, newRefCoordinate.targetAlignmentStart, newRefCoordinate.targetAlignmentStop]
				writer.writerow(data_row)
				no_of_loci_with_1_newRef += 1
			else:
				no_of_loci_with_1Plus_newRef += 1
		del writer
		sys.stderr.write("%s loci found unique new coordinates. %s loci found >1 new coordinates.\n"%\
						(no_of_loci_with_1_newRef, no_of_loci_with_1Plus_newRef))

	def outputQuerySNPID2NewCoordinateInUCSCChainFormat(self, querySNPID2NewReferenceCoordinateLs=None, \
						querySNPID2NewRefCoordinateOutputFname=None):
		"""
		2013.07.05
			if a query snp has >1 new reference coordinates, skip it.
			chain file format, http://genome.ucsc.edu/goldenPath/help/chain.html.
				Header Lines (0-based exclusive coordinates system, no actual english header though)
					chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
				Alignment Data Lines, three required attribute values:
					size dt dq
			NOTE: The last line of the alignment section contains only one number:
				the ungapped alignment size of the last block.
		"""
		no_of_loci_with_1_newRef = 0
		no_of_loci_with_1Plus_newRef = 0
		sys.stderr.write("Outputting %s pairs in querySNPID2NewReferenceCoordinateLs to %s in chain format ..."%\
						(len(querySNPID2NewReferenceCoordinateLs), querySNPID2NewRefCoordinateOutputFname))
		writer = csv.writer(open(querySNPID2NewRefCoordinateOutputFname, 'w'), delimiter='\t')
		chainID=1
		for querySNPID, newRefCoordinateLs in querySNPID2NewReferenceCoordinateLs.iteritems():
			if len(newRefCoordinateLs)==1:
				newRefCoordinate = newRefCoordinateLs[0]
				data_row = ["chain", 1, \
						newRefCoordinate.queryChromosome, newRefCoordinate.queryAlignmentSpan,
						newRefCoordinate.queryStrand,
						newRefCoordinate.queryAlignmentStart-1, newRefCoordinate.queryAlignmentStop, 
						
						newRefCoordinate.newChr, newRefCoordinate.targetAlignmentSpan, \
						"+", \
						newRefCoordinate.targetAlignmentStart-1, newRefCoordinate.targetAlignmentStop,
						chainID]
				chainID += 1
				writer.writerow(data_row)
				writer.writerow([newRefCoordinate.targetAlignmentSpan])
				writer.writerow([])
				no_of_loci_with_1_newRef += 1
			else:
				no_of_loci_with_1Plus_newRef += 1
		del writer
		sys.stderr.write("%s loci found unique new coordinates (=number of chains). %s loci found >1 new coordinates.\n"%\
						(no_of_loci_with_1_newRef, no_of_loci_with_1Plus_newRef))
	
	def outputSwitchPointInfo(self, querySNPID2NewReferenceCoordinateLs=None, outputFname=None):
		"""
		2013.07.11
			output the switch point (adjacent sites mapped to two different chromosomes) information
		"""
		
		sys.stderr.write("Converting querySNPID2NewReferenceCoordinateLs to oldCoordinateKey2newCoordinateDataLs ... ")
		oldCoordinateKey2newCoordinateDataLs = {}
		counter = 0
		for querySNPID, newRefCoordinateLs in querySNPID2NewReferenceCoordinateLs.iteritems():
			oldCoordinateKey = None
			counter += len(newRefCoordinateLs)
			for newRefCoordinate in newRefCoordinateLs:
				if oldCoordinateKey is None:
					oldCoordinateKey = (newRefCoordinate.queryChromosome, newRefCoordinate.queryStart, newRefCoordinate.queryStop)
				if oldCoordinateKey not in oldCoordinateKey2newCoordinateDataLs:
					oldCoordinateKey2newCoordinateDataLs[oldCoordinateKey] = []
				oldCoordinateKey2newCoordinateDataLs[oldCoordinateKey].append(newRefCoordinate)
		sys.stderr.write(" %s old coordinate keys with %s new coordinates.\n"%(len(oldCoordinateKey2newCoordinateDataLs),\
																		counter))
		
		sys.stderr.write("Finding switch points ...")
		counter =0
		real_counter = 0
		noOfRecordsWithMultiNewCoords = 0

		oldChromosome2SwitchData = {}
		
		oldCoordinateKeyLs = oldCoordinateKey2newCoordinateDataLs.keys()
		oldCoordinateKeyLs.sort()
		for oldCoordinateKey in oldCoordinateKeyLs:
			counter +=1
			newRefCoordinateLs = oldCoordinateKey2newCoordinateDataLs.get(oldCoordinateKey)
			
			oldChromosome = oldCoordinateKey[0]
			
			if oldChromosome not in oldChromosome2SwitchData:
				oldChromosome2SwitchData[oldChromosome] = PassingData(noOfLociWithUniqueHit=0, noOfLoci=0, \
														spanStart=oldCoordinateKey[1], \
														spanStop=oldCoordinateKey[2], noOfSwitchPoints=0,\
														previousNewChromosome=None, previousNewRefStart=None,\
														previousNewRefStop=None,\
														previousOrientationOnNewChromosome=None)
			
			switchData = oldChromosome2SwitchData[oldChromosome]
			switchData.noOfLoci += 1
			
			if len(newRefCoordinateLs)>1:
				noOfRecordsWithMultiNewCoords += 1
				continue
			
			switchData.noOfLociWithUniqueHit += 1
			newRefCoordinate = newRefCoordinateLs[0]
			
			if switchData.previousNewChromosome is not None:
				if newRefCoordinate.newChr!=switchData.previousNewChromosome:
					switchData.noOfSwitchPoints += 1
					#reset the orientation
					switchData.previousOrientationOnNewChromosome = None
					
				else:	#on the same chromosome
					currentOrientation = (newRefCoordinate.newRefStart - switchData.previousNewRefStart)>=0
					if switchData.previousOrientationOnNewChromosome is not None:
						if currentOrientation !=switchData.previousOrientationOnNewChromosome:
							switchData.noOfSwitchPoints += 1
					switchData.previousOrientationOnNewChromosome = currentOrientation
					
			#adjust the spanStop
			if newRefCoordinate.queryStop > switchData.spanStop:
				switchData.spanStop = newRefCoordinate.queryStop
					
			
			switchData.previousNewChromosome = newRefCoordinate.newChr
			switchData.previousNewRefStart = newRefCoordinate.newRefStart
			switchData.previousNewRefStop = newRefCoordinate.newRefStop
			real_counter  += 1
		if counter >0:
			fraction = real_counter/float(counter)
		else:
			fraction = -1
		sys.stderr.write("%s (out of %s, %s) records found new coordinates. %s records with >1 new coordinates, discarded.\n"%(real_counter, counter, \
																	fraction, noOfRecordsWithMultiNewCoords))
		
		
		sys.stderr.write("Outputting switch points of %s old chromosomes ..."%(len(oldChromosome2SwitchData)))
		statFile = MatrixFile(inputFname=outputFname, openMode='w', delimiter='\t')
		header = ['oldChromosome', "noOfSwitchPoints", "regionSpan", "noOfLociWithUniqueHit", "noOfSwitchesPerLocus", "noOfLoci"]
		statFile.writeHeader(header)
		noOfTotalSwitchPoints = 0
		noOfTotalLoci = 0
		for oldChromosome, switchData in oldChromosome2SwitchData.iteritems():
			if switchData.noOfLociWithUniqueHit>0:
				switchPointFraction = switchData.noOfSwitchPoints/float(switchData.noOfLociWithUniqueHit)
			else:
				switchPointFraction = -1
			data_row = [oldChromosome, switchData.noOfSwitchPoints, switchData.spanStop-switchData.spanStart+1, \
					switchData.noOfLociWithUniqueHit, switchPointFraction, len(oldCoordinateKey2newCoordinateDataLs)]
			statFile.writerow(data_row)
			noOfTotalSwitchPoints += switchData.noOfSwitchPoints
			noOfTotalLoci += switchData.noOfLociWithUniqueHit
		statFile.close()
		sys.stderr.write(' %s total switch points, %s total loci with unique hit.\n'%(noOfTotalSwitchPoints, noOfTotalLoci))
	
	def outputSNPDataInNewCoordinate(self, querySNPDataFname=None, querySNPID2NewReferenceCoordinateLs=None,\
									newSNPDataOutputFname=None, newSNPDataOutputFormat=1):
		"""
		2013.07.03 added argument newSNPDataOutputFormat
			
		2012.10.14
			split out of findSNPPositionOnNewRef()
		"""
		sys.stderr.write("Converting querySNPDataFname %s into individual X SNP format, format=%s ... "%\
						(querySNPDataFname, newSNPDataOutputFormat))
		"""
Sample  Geno    SNP
1999010 CC      cs_primer1082_247
1999068 CC      cs_primer1082_247
2000022 CT      cs_primer1082_247
2000064 CT      cs_primer1082_247
2000117 CC      cs_primer1082_247

		"""
		inf = utils.openGzipFile(querySNPDataFname)
		reader = csv.reader(inf, delimiter=figureOutDelimiter(inf))
		col_name2index = getColName2IndexFromHeader(reader.next())
		
		sampleIndex = col_name2index.get("Sample")
		genotypeIndex = col_name2index.get("Geno")
		SNPIDIndex = col_name2index.get("SNP")
		
		row_id2index = {}
		row_id_ls = []
		col_id_ls = []
		col_id2index = {}
		row_col_index2genotype = {}
		for row in reader:
			sampleID = row[sampleIndex]
			genotype = row[genotypeIndex]
			querySNPID = row[SNPIDIndex]
			if querySNPID in querySNPID2NewReferenceCoordinateLs:
				newRefCoordinateLs = querySNPID2NewReferenceCoordinateLs.get(querySNPID)
				if len(newRefCoordinateLs)==1:
					newRefCoordinate = newRefCoordinateLs[0]
					if newSNPDataOutputFormat==2:
						col_id = '%s_%s'%(newRefCoordinate.newChr, newRefCoordinate.newRefStart)
					else:
						col_id = '%s_%s_%s'%(newRefCoordinate.newChr, newRefCoordinate.newRefStart, newRefCoordinate.newRefStop)
					queryStrand = newRefCoordinate.queryStrand
					if col_id not in col_id2index:
						col_id2index[col_id] = len(col_id2index)
						col_id_ls.append(col_id)
					if sampleID not in row_id2index:
						row_id2index[sampleID] = len(row_id2index)
						row_id_ls.append(sampleID)
					if queryStrand == "-":
						genotype = SNP.reverseComplement(genotype)
					row_index = row_id2index[sampleID]
					col_index = col_id2index[col_id]
					row_col_index2genotype[(row_index, col_index)] = genotype
				else:
					continue
		data_matrix = numpy.zeros([len(row_id_ls), len(col_id2index)], dtype=numpy.int8)
		
		for row_col_index, genotype in row_col_index2genotype.iteritems():
			row_index, col_index = row_col_index[:2]
			data_matrix[row_index, col_index] = SNP.nt2number[genotype]
		sys.stderr.write("\n")
		snpData = SNP.SNPData(row_id_ls=row_id_ls, col_id_ls=col_id_ls, data_matrix=data_matrix)
		snpData.tofile(newSNPDataOutputFname)
		
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
	
		querySNPID2NewReferenceCoordinateLs = self.findSNPPositionOnNewRef(SNPFlankSequenceFname=self.SNPFlankSequenceFname, \
						blastHitResultFname=self.inputFname, querySNPDataFname=self.querySNPDataFname,\
						querySNPID2NewRefCoordinateOutputFname=self.outputFname, \
						newSNPDataOutputFname=self.newSNPDataOutputFname, minAlignmentSpan=self.minAlignmentSpan)
		
		#output the mapping
		self.outputQuerySNPID2NewCoordinateMap(querySNPID2NewReferenceCoordinateLs=querySNPID2NewReferenceCoordinateLs, \
						querySNPID2NewRefCoordinateOutputFname=self.outputFname)
		if self.chainFilename:
			self.outputQuerySNPID2NewCoordinateInUCSCChainFormat(querySNPID2NewReferenceCoordinateLs=querySNPID2NewReferenceCoordinateLs, \
						querySNPID2NewRefCoordinateOutputFname=self.chainFilename)
		if self.switchPointFname:
			self.outputSwitchPointInfo(querySNPID2NewReferenceCoordinateLs=querySNPID2NewReferenceCoordinateLs, \
									outputFname=self.switchPointFname)

if __name__ == '__main__':
	main_class = FindSNPPositionOnNewRefFromFlankingBlastOutput
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()