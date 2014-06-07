#!/usr/bin/env python
"""
Examples:
	%s -i ~/NetworkData/vervet/db/genotype_file/method_225/94276_VCF_CAE6.sorted.vcf.gz
		--missingStatFname ./tmp/94276_VCF_CAE6.sorted_missing_stat.tsv
		--alignmentFilename ~/NetworkData/vervet/db/individual_alignment/5374_641_1986014_GA_vs_3488_by_method6_realigned0_reduced0.bam
		--alignmentMedianDepth 34
		-o ./tmp/94276_VCF_CAE6.sorted.marked.vcf
	
	%s 
	
	%s 
	
Description:
	2013.12.04 
		Two output
			one is modified VCF file.
			one is missing stat file.
				header = ["locusID", 'chromosome', 'start', 'stop', 'occurrence', 'missingReason']
		missingReason:
			1: outside range of permitted depth range
			2: low mapping quality (>10 percent  of low mapping quality reads)
			3: both

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, MatrixFile, PassingData
from pymodule.yhio.VCFFile import VCFFile
from pymodule.pegasus.mapper.AbstractVCFMapper import AbstractVCFMapper
from pymodule import SNP
import pysam
import numpy

parentClass = AbstractVCFMapper
class MarkGenotypeMissingByAlignmentQuality(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update({
					('alignmentFilename', 1, ): [None, '', 1, 'alignment file in bam format corresponding to the sample in VCF file', ],\
					('alignmentMedianDepth', 1, float): [None, '', 1, 'global median coverage of the alignment file', ],\
					('alignmentDepthFold', 0, float): [2.0, '', 1, 'genotype whose depth is outside the range of [1/foldChange, foldChange]*medianDepth', ],\
					('minMapQGoodRead', 0, float): [15, '', 1, 'the minimum mapQ for a good-quality aligned read', ],\
					('minFractionOfGoodRead', 0, float): [0.9, '', 1, 'the minimum fraction of good-quality aligned reads at one locus', ],\
					('missingStatFname', 1, ): [None, '', 1, 'file to contain the genotype missing stat output', ],\
					('sampleID', 0, ): [None, '', 1, 'name of the sample from this alignment file, default is alignment file basename', ],\
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		parentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		
		if not self.sampleID:
			self.sampleID = os.path.basename(self.alignmentFilename)
	
	def returnLocusLowMapQualityIndicator(self, alignedReadLs=None, minMapQGoodRead=2, minFractionOfGoodRead=0.9):
		"""
		2013.12.04
		"""
		totalNoOfReads = 0
		noOfGoodReads = 0.0
		medianMapQ=-10
		mapQList=[]
		for alignedRead in alignedReadLs:
			totalNoOfReads +=1
			mapQList.append(alignedRead.mapq)
			if alignedRead.mapq>=minMapQGoodRead:
				noOfGoodReads += 1
			else:
				pass
		if totalNoOfReads>0:
			fractionOfGoodRead = noOfGoodReads/(totalNoOfReads)
			medianMapQ = numpy.median(mapQList)
		else:
			fractionOfGoodRead = -1
			medianMapQ = -10
			
		if fractionOfGoodRead>=minFractionOfGoodRead:
			locusLowMapQIndicator = 0
		else:
			locusLowMapQIndicator = 2
		return PassingData(locusLowMapQIndicator=locusLowMapQIndicator, totalNoOfReads=totalNoOfReads, \
						noOfGoodReads=noOfGoodReads, fractionOfGoodRead=fractionOfGoodRead,\
						medianMapQ=medianMapQ)
	
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		outputDir = os.path.split(self.outputFname)[0]
		if outputDir and not os.path.isdir(outputDir):
			os.makedirs(outputDir)
		
		reader = VCFFile(inputFname=self.inputFname)
		
		alignmentFile = pysam.Samfile(self.alignmentFilename, "rb")
		
		writer = VCFFile(outputFname=self.outputFname, openMode='w')
		writer.metaInfoLs = reader.metaInfoLs
		writer.header = reader.header
		writer.writeMetaAndHeader()
		
		statWriter = MatrixFile(self.missingStatFname, openMode='w', delimiter='\t')
		header = ["sampleID", "locusID", 'chromosome', 'start', 'stop', 'occurrence', 'missingReason', \
				'fractionOfGoodRead', 'medianMapQ', 'totalNoOfReads']
		statWriter.writeHeader(header)
		
		counter = 0
		real_counter = 0
		minDepth = self.alignmentMedianDepth/self.alignmentDepthFold
		maxDepth = self.alignmentMedianDepth*self.alignmentDepthFold
		
		for vcfRecord in reader:
			locusID = "%s_%s"%(vcfRecord.chromosome, vcfRecord.position)
			alignedReadLs = alignmentFile.fetch(vcfRecord.chromosome, vcfRecord.position-1, vcfRecord.position+1)	#start and end in fetch() are 0-based.
			locusLowMapQData = self.returnLocusLowMapQualityIndicator(alignedReadLs=alignedReadLs,\
												minMapQGoodRead=self.minMapQGoodRead, minFractionOfGoodRead=self.minFractionOfGoodRead)
			locusLowMapQIndicator = locusLowMapQData.locusLowMapQIndicator
			depth = locusLowMapQData.totalNoOfReads
			if depth>=minDepth and depth <=maxDepth:
				locusOutOfDepthIndicator = 0 	#good
			else:
				locusOutOfDepthIndicator = 1
			
			locusLowQualityIndicator = locusOutOfDepthIndicator + locusLowMapQIndicator
			data_row = [self.sampleID, locusID, vcfRecord.chromosome, vcfRecord.position, vcfRecord.position,\
						1, locusLowQualityIndicator, locusLowMapQData.fractionOfGoodRead, \
						locusLowMapQData.medianMapQ, locusLowMapQData.totalNoOfReads]
			statWriter.writerow(data_row)
			if locusLowQualityIndicator>0:
				real_counter += 1
				#modify the VCF record
				#get sample ID column, then set its genotype missing
				vcfRecord.setGenotypeCallForOneSample(sampleID=self.sampleID, genotype="./.", convertGLToPL=True)
			#2014.1.4 output VCF record
			writer.writeVCFRecord(vcfRecord)
			counter += 1
		reader.close()
		statWriter.close()
		writer.close()
		sys.stderr.write("%s (out of %s, %s) genotypes marked missing.\n"%(real_counter, counter, \
												real_counter/float(counter)))

if __name__ == '__main__':
	main_class = MarkGenotypeMissingByAlignmentQuality
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()