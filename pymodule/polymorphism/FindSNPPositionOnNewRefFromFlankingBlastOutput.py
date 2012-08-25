#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -o ~/script/vervet/data/194SNPData/isq524CoordinateSNPData_max15Mismatch.tsv
		-i ~/script/vervet/data/194SNPData/AllSNPData.txt
		-x ~/script/vervet/data/194SNPData/externalSNPID2ISQ524Coordinate_max15Mismatch.tsv
		-l /Network/Data/vervet/vervetPipeline/Blast/Blast194SNPFlankAgainst524_15Mismatches.2012.8.17T2334/folderBlast/blast.tsv
		-S ~/script/vervet/data/194SNPData/AllSNPFlankWithSNPMark.txt

Description:
	2012.8.19
		given a fasta file of flanking sequences of SNPs, find its new positions on the blast database (new reference)
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, SNP
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
import numpy

class FindSNPPositionOnNewRefFromFlankingBlastOutput(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	#option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
							('SNPFlankSequenceFname', 1, ): ['', 'S', 1, 'in fasta format, with SNP embedded as [A/C]', ],\
							('blastHitMergeFname', 1, ): ['', 'l', 1, 'BlastWorkflow.py output: blast SNPFlankSequenceFname (converting [A/C] to A) to a new reference', ],\
							('externalSNPID2RefCoordinateOutputFname', 0, ): [None, 'x', 1, '', ],\
							('minNoOfIdentities', 0, int): [None, 'm', 1, 'minimum number of identities between a query and target', ],\
							('maxNoOfMismatches', 0, int): [None, 'a', 1, 'minimum number of mismatches between a query and target', ],\
							('minIdentityPercentage', 0, float): [None, 'n', 1, 'minimum percentage of identities between a query and target', ],\
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		2012.8.19
		"""
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def convert194SNPDataIntoVervetRefSNPData(self, SNPFlankSequenceFname=None, blastHitMergeFname=None, originalSNPDataFname=None,\
											externalSNPID2RefCoordinateOutputFname=None, outputFname=None):
		"""
		2012.8.19
			outputFname will contain the individual X SNP matrix.
		"""
		sys.stderr.write("Deriving externalSNPID2positionInFlank ...\n")
		externalSNPID2positionInFlank = {}
		from Bio import SeqIO
		inf = open(SNPFlankSequenceFname, "rU")
		counter = 0
		for record in SeqIO.parse(inf, "fasta"):
			externalSNPID = record.id.split()[0]	#get rid of extra comment
			snpPositionInFlank = record.seq.find('[') + 1
			if snpPositionInFlank<=0:
				sys.stderr.write("Warning, could not find position for snp %s  in the flanking sequence \n"%(record.id))
			else:
				externalSNPID2positionInFlank[externalSNPID] = snpPositionInFlank
			counter += 1
		inf.close()
		sys.stderr.write(" %s/%s SNPs with positions.\n"%(len(externalSNPID2positionInFlank), counter))
		
		sys.stderr.write("Finding blast reference coordinates for %s SNPs. \n"%(len(externalSNPID2positionInFlank)))
		import csv
		from pymodule import figureOutDelimiter, getColName2IndexFromHeader
		reader = csv.reader(open(blastHitMergeFname), delimiter='\t')
		header =reader.next()
		col_name2index = getColName2IndexFromHeader(header)
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
		externalSNPID2BlastRefCoordinate = {}
		counter = 0
		real_counter = 0
		for row in reader:
			queryID = row[queryIDIndex].split()[0]	##get rid of extra comment
			queryStart = int(row[queryStartIndex])
			queryEnd = int(row[queryEndIndex])
			
			targetChr = row[targetChrIndex]
			targetStart = int(row[targetStartIndex])
			targetStop = int(row[targetStopIndex])
			
			queryAlnSpan = abs(queryEnd-queryStart) + 1
			targetAlnSpan = abs(targetStop-targetStart) + 1
			if queryAlnSpan == targetAlnSpan:
				if queryID in externalSNPID2positionInFlank:
					positionInFlank = externalSNPID2positionInFlank.get(queryID)
					if queryStart >queryEnd:	#could happen. on the opposite strand. targetStart is always bigger than targetstop
						blastRefPosition =  targetStop - (positionInFlank-queryEnd)
						strand = -1
					else:
						blastRefPosition =  targetStart + (positionInFlank-queryStart)
						strand = + 1
					externalSNPID2BlastRefCoordinate[queryID] = (targetChr, blastRefPosition, strand)
					real_counter += 1
			counter += 1
		sys.stderr.write(" from %s blast results. %s/%s SNPs found blast-reference coordinates.\n"%\
						(counter, real_counter, len(externalSNPID2positionInFlank)))
		
		if externalSNPID2RefCoordinateOutputFname:
			sys.stderr.write("Outputting %s pairs in externalSNPID2BlastRefCoordinate to %s ..."%\
							(len(externalSNPID2BlastRefCoordinate), externalSNPID2RefCoordinateOutputFname))
			writer = csv.writer(open(externalSNPID2RefCoordinateOutputFname, 'w'), delimiter='\t')
			header = ['externalSNPID', 'targetChr', 'targetPosition']
			writer.writerow(header)
			for externalSNPID, blastRefCoordinate in externalSNPID2BlastRefCoordinate.iteritems():
				targetChr, targetPosition = blastRefCoordinate[:2]
				data_row = [externalSNPID, targetChr, targetPosition]
				writer.writerow(data_row)
			del writer
			sys.stderr.write("\n")
		
		if originalSNPDataFname and outputFname:
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
				externalSNPID = row[SNPIDIndex]
				if externalSNPID in externalSNPID2BlastRefCoordinate:
					blastRefCoordinate = externalSNPID2BlastRefCoordinate.get(externalSNPID)
					col_id = '%s_%s'%(blastRefCoordinate[0], blastRefCoordinate[1])
					strand = blastRefCoordinate[2]
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
			data_matrix = numpy.zeros([len(row_id_ls), len(col_id2index)], dtype=numpy.int8)
			
			for row_col_index, genotype in row_col_index2genotype.iteritems():
				row_index, col_index = row_col_index[:2]
				data_matrix[row_index, col_index] = SNP.nt2number[genotype]
			sys.stderr.write("\n")
			snpData = SNP.SNPData(row_id_ls=row_id_ls, col_id_ls=col_id_ls, data_matrix=data_matrix)
			snpData.tofile(outputFname)
				
			
		
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
	
		self.convert194SNPDataIntoVervetRefSNPData(SNPFlankSequenceFname=self.SNPFlankSequenceFname, \
						blastHitMergeFname=self.blastHitMergeFname, originalSNPDataFname=self.inputFname,\
						externalSNPID2RefCoordinateOutputFname=self.externalSNPID2RefCoordinateOutputFname, \
						outputFname=self.outputFname)
		

if __name__ == '__main__':
	main_class = FindSNPPositionOnNewRefFromFlankingBlastOutput
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()