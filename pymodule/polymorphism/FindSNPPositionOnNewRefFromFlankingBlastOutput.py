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
	2012.8.19
		given a fasta file of flanking sequences of polymorphic loci (SNP or SVs).
		find its new positions on the blast database (new reference).
		
		inputFname is BlastWorkflow.py's output: blast SNPFlankSequenceFname (converting [A/C] to A) to a new reference.
		outputFname is the map between old and new coordinates.
			originalSNPID, newChr, newRefStart, newRefStop, strand, targetAlignmentSpan.
		
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, SNP
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
	if given, use this to fetch a map between originalSNPID & snpPositionInFlank.', ],\
							('originalSNPDataFname', 0, ): ['', 'l', 1, 'path to an input tsv/csv, 3-column: Sample SNP Geno.', ],\
							('newSNPDataOutputFname', 0, ): [None, 'x', 1, 'if given, would be a Strain X Locus Yu format polymorphism file with new coordinates.', ],\
							('minNoOfIdentities', 0, int): [None, 'm', 1, 'minimum number of identities between a query and target', ],\
							('maxNoOfMismatches', 0, int): [None, 'a', 1, 'minimum number of mismatches between a query and target', ],\
							('minIdentityFraction', 0, float): [None, 'n', 1, 'minimum fraction of identities between a query and target', ],\
							('minAlignmentSpan', 1, int): [10, '', 1, 'minimum number of bases of the query and target that are involved in the blast alignment', ],\
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		2012.8.19
		"""
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def getOriginalSNPID2positionInFlank(self, SNPFlankSequenceFname=None):
		"""
		2012.10.8
			split out of findSNPPositionOnNewRef()
			the snpPositionInFlank is 1-based.
			add one more variable, locusSpan, in the value of originalSNPID2positionInFlank.
				=0 for SNPs, =length-1 for other loci.
			
			The SNPFlankSequenceFname is in fasta format, with SNP embedded as [A/C] in sequence.
		"""
		sys.stderr.write("Deriving originalSNPID2positionInFlank from %s ...\n"%(SNPFlankSequenceFname))
		originalSNPID2positionInFlank = {}
		from Bio import SeqIO
		inf = open(SNPFlankSequenceFname, "rU")
		counter = 0
		for record in SeqIO.parse(inf, "fasta"):
			originalSNPID = record.id.split()[0]	#get rid of extra comment
			snpPositionInFlank = record.seq.find('[') + 1
			refBase = record.seq[snpPositionInFlank]
			altBase = record.seq[snpPositionInFlank+2]
			if snpPositionInFlank<=0:
				sys.stderr.write("Warning, could not find position for snp %s  in the flanking sequence \n"%(record.id))
			else:
				originalSNPID2positionInFlank[originalSNPID] = PassingData(snpPositionInFlank=snpPositionInFlank, locusSpan=0,\
																refBase=refBase, altBase=altBase)
				#
			counter += 1
		inf.close()
		sys.stderr.write(" %s/%s SNPs with positions.\n"%(len(originalSNPID2positionInFlank), counter))
		return originalSNPID2positionInFlank
	

	def parseOriginalLocusID(self, locus_id=None):
		"""
		2012.10.8
			locus_id is in the format of '%s_%s_%s_positionInFlank%s'%(chr, start, stop, flankingLength+1)
			output of ExtractFlankingSequenceForVCFLoci.py
		"""
		search_result = ExtractFlankingSequenceForVCFLoci.sequenceTitlePattern.search(locus_id)
		chr = None
		start = None
		stop = None
		refBase = None
		altBase = None
		positionInFlank = None
		if search_result:
			chr = search_result.group(1)
			start = search_result.group(2)
			stop = search_result.group(3)
			refBase = search_result.group(4)
			altBase = search_result.group(5)
			positionInFlank = search_result.group(6)
			
		return PassingData(chr=chr, start=start, stop=stop, refBase=refBase, altBase=altBase, positionInFlank=positionInFlank)
		
	
	def findSNPPositionOnNewRef(self, SNPFlankSequenceFname=None, blastHitResultFname=None, \
							originalSNPDataFname=None,\
							originalSNPID2NewRefCoordinateOutputFname=None, newSNPDataOutputFname=None, minAlignmentSpan=10):
		"""
		2012.10.8
			argument minAlignmentSpan: the number of bases involved in the blast query-target alignment
		2012.8.19
			newSNPDataOutputFname will contain the individual X SNP matrix.
		"""
		if SNPFlankSequenceFname:
			originalSNPID2positionInFlank = self.getOriginalSNPID2positionInFlank(SNPFlankSequenceFname=SNPFlankSequenceFname)
		else:
			originalSNPID2positionInFlank = None
		
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
		originalSNPID2BlastRefCoordinateLs = {}
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
				if originalSNPID2positionInFlank and queryID in originalSNPID2positionInFlank:
					parseData = originalSNPID2positionInFlank.get(queryID)
					locusSpan = parseData.locusSpan
				else:
					parseData = self.parseOriginalLocusID(queryID)
					chr = parseData.chr
					start = parseData.start
					stop = parseData.stop
					if start is not None and stop is not None:
						stop = int(stop)
						start = int(start)
						locusSpan = abs(int(stop)-start)	#length-1
					else:
						locusSpan = None
				positionInFlank = parseData.positionInFlank
				originalRefBase = parseData.refBase
				originalAltBase = parseData.altBase
				if positionInFlank is not None and locusSpan is not None:
					positionInFlank = int(positionInFlank)
					if targetAlignmentSpan>=minAlignmentSpan and queryAlignmentSpan>=minAlignmentSpan:
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
						
			counter += 1
		sys.stderr.write(" from %s blast results. %s/%s SNPs found blast-reference coordinates.\n"%\
						(counter, real_counter, len(queryIDSet)))
		
		#output the mapping
		if originalSNPID2NewRefCoordinateOutputFname:
			self.outputOriginalSNPID2NewCoordinateMap(originalSNPID2BlastRefCoordinateLs=originalSNPID2BlastRefCoordinateLs, \
						originalSNPID2NewRefCoordinateOutputFname=originalSNPID2NewRefCoordinateOutputFname)
		
		if originalSNPDataFname and newSNPDataOutputFname:
			self.outputSNPDataInNewCoordinate(originalSNPDataFname=originalSNPDataFname, \
									originalSNPID2BlastRefCoordinateLs=originalSNPID2BlastRefCoordinateLs, newSNPDataOutputFname=newSNPDataOutputFname)
	
	def outputOriginalSNPID2NewCoordinateMap(self, originalSNPID2BlastRefCoordinateLs=None, \
						originalSNPID2NewRefCoordinateOutputFname=None):
		"""
		2012.10.14
			split out of findSNPPositionOnNewRef()
		"""
		no_of_loci_with_1_newRef = 0
		no_of_loci_with_1Plus_newRef = 0
		sys.stderr.write("Outputting %s pairs in originalSNPID2BlastRefCoordinateLs to %s ..."%\
						(len(originalSNPID2BlastRefCoordinateLs), originalSNPID2NewRefCoordinateOutputFname))
		writer = csv.writer(open(originalSNPID2NewRefCoordinateOutputFname, 'w'), delimiter='\t')
		header = ['originalSNPID', 'originalRefBase', 'originalAltBase', 'newChr', 'newRefStart', 'newRefStop', 'strand', \
				'queryAlignmentSpan', 'targetAlignmentSpan']
		writer.writerow(header)
		for originalSNPID, blastRefCoordinateLs in originalSNPID2BlastRefCoordinateLs.iteritems():
			if len(blastRefCoordinateLs)==1:
				blastRefCoordinate = blastRefCoordinateLs[0]
				data_row = [originalSNPID, blastRefCoordinate.originalRefBase, blastRefCoordinate.originalAltBase, \
						blastRefCoordinate.newChr, blastRefCoordinate.newRefStart, blastRefCoordinate.newRefStop, \
						blastRefCoordinate.strand, blastRefCoordinate.queryAlignmentSpan,
						blastRefCoordinate.targetAlignmentSpan]
				writer.writerow(data_row)
				no_of_loci_with_1_newRef += 1
			else:
				no_of_loci_with_1Plus_newRef += 1
		del writer
		sys.stderr.write("%s loci found unique new coordinates. %s loci found >1 new coordinates.\n"%\
						(no_of_loci_with_1_newRef, no_of_loci_with_1Plus_newRef))
	
	def outputSNPDataInNewCoordinate(self, originalSNPDataFname=None, originalSNPID2BlastRefCoordinateLs=None,\
									newSNPDataOutputFname=None):
		"""
		2012.10.14
			split out of findSNPPositionOnNewRef()
		"""
		sys.stderr.write("Converting originalSNPDataFname %s into individual X SNP format ... "%(originalSNPDataFname))
		"""
Sample  Geno    SNP
1999010 CC      cs_primer1082_247
1999068 CC      cs_primer1082_247
2000022 CT      cs_primer1082_247
2000064 CT      cs_primer1082_247
2000117 CC      cs_primer1082_247

		"""
		inf = utils.openGzipFile(originalSNPDataFname)
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
			originalSNPID = row[SNPIDIndex]
			if originalSNPID in originalSNPID2BlastRefCoordinateLs:
				blastRefCoordinateLs = originalSNPID2BlastRefCoordinateLs.get(originalSNPID)
				if len(blastRefCoordinateLs)==1:
					blastRefCoordinate = blastRefCoordinateLs[0]
					col_id = '%s_%s_%s'%(blastRefCoordinate.newChr, blastRefCoordinate.newRefStart, blastRefCoordinate.newRefStop)
					strand = blastRefCoordinate.strand
					if col_id not in col_id2index:
						col_id2index[col_id] = len(col_id2index)
						col_id_ls.append(col_id)
					if sampleID not in row_id2index:
						row_id2index[sampleID] = len(row_id2index)
						row_id_ls.append(sampleID)
					if strand ==-1:
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
	
		self.findSNPPositionOnNewRef(SNPFlankSequenceFname=self.SNPFlankSequenceFname, \
						blastHitResultFname=self.inputFname, originalSNPDataFname=self.originalSNPDataFname,\
						originalSNPID2NewRefCoordinateOutputFname=self.outputFname, \
						newSNPDataOutputFname=self.newSNPDataOutputFname, minAlignmentSpan=self.minAlignmentSpan)
		

if __name__ == '__main__':
	main_class = FindSNPPositionOnNewRefFromFlankingBlastOutput
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()