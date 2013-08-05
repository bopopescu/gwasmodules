#!/usr/bin/env python
"""
Examples:
	%s 
	
	# 
	%s  -i  mendelRuntype3_s2Reduce/merged_lmendel_chr_pos.tsv  -o  mendelRuntype3_s2Reduce/medianNoOfMendelErrors_along_chromosome.tsv
		--genome_drivername=postgresql --genome_hostname=localhost --genome_dbname=vervetdb --genome_schema=genome
		--genome_db_user=yh --genome_db_passwd ...
		--chromosomeHeader Chromosome --tax_id 60711 --sequence_type_id 1 --chrOrder 1
		--samplingRate 1 --positionHeader Start --whichColumnHeader N --valueForNonPositiveYValue -1 --inputFileFormat 1

Description:
	2013.07.31 a program that calculates moving averages in windows (overlap or not) along the genome.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv, random, numpy
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.yhio.AbstractGenomeFileWalker import AbstractGenomeFileWalker
from pymodule.algorithm.RBTree import RBDict
from pymodule.yhio.CNV import CNVSegmentBinarySearchTreeKey

parentClass = AbstractGenomeFileWalker
class GenomeMovingAverageStatistics(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	#option_default_dict.update(AbstractMapper.db_option_dict.copy())
	option_default_dict.update({
					('windowSize', 0, int): [200000, '', 1, 'size of the moving window'], \
					('windowOverlapSize', 0, int): [0, '', 1, 'size of the overlap between adjacent windows'], \
					
					('run_type', 0, int): [1, '', 1, '1: median within each window; 2: mean within each window; \
	3: fraction above minimum value'], \
					('minValueForFraction', 0, float): [None, '', 1, 'the minimum value for run_type 3'],\
					
					})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		parentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
		
		#2013.07.31
		fractionFunction = lambda ls: sum([a>=self.minValueForFraction for a in ls])/float(len(ls))
		reduceType2Function = {1: numpy.median, 2: numpy.mean, 3: fractionFunction}
		self.reduceFunction = reduceType2Function.get(self.run_type, numpy.median)
	
	def setup(self, **keywords):
		"""
		2013.07.31
			construct an RBTree dictionary map between windows and their data
		"""
		parentClass.setup(self, **keywords)
		
		sys.stderr.write("Constructing segmentKey2dataLsRBDict ...")
		self.segmentKey2dataLsRBDict = RBDict()
		counter =0
		for chromosome, chromosomeSize in self.oneGenomeData.chr_id2size.iteritems():
			no_of_intervals = max(1, int(math.ceil(chromosomeSize/float(self.windowSize)))-1)
			for i in range(no_of_intervals):
				originalStartPos = i*self.windowSize + 1
				#to render adjacent intervals overlapping because trioCaller uses LD
				startPos = max(1, originalStartPos-self.windowOverlapSize)
				if i<no_of_intervals-1:
					originalStopPos = min((i+1)*self.windowSize, chromosomeSize)
				else:	#last chunk, include bp till the end
					originalStopPos = chromosomeSize
				#to render adjacent intervals overlapping because trioCaller uses LD
				stopPos = min(chromosomeSize, originalStopPos+self.windowOverlapSize)
				
				segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=chromosome, span_ls=[startPos, stopPos],\
														min_reciprocal_overlap=1.0)
				
				interval = "%s:%s-%s"%(chromosome, originalStartPos, originalStopPos)
				intervalFileBasenameSignature = '%s_%s_%s'%(chromosome, originalStartPos, originalStopPos)
				overlapInterval = "%s:%s-%s"%(chromosome, startPos, stopPos)
				overlapIntervalFileBasenameSignature = '%s_%s_%s'%(chromosome, startPos, stopPos)
				span = stopPos-startPos+1
				intervalData = PassingData(overlapInterval=overlapInterval, overlapIntervalFileBasenameSignature=overlapIntervalFileBasenameSignature,\
							interval=interval, intervalFileBasenameSignature=intervalFileBasenameSignature, \
							chr=chromosome, chromosome=chromosome, chromosomeSize=chromosomeSize,\
							originalStartPos=originalStartPos, originalStopPos=originalStopPos, \
							start=startPos, stop=stopPos, \
							overlapStart=startPos, overlapStop=stopPos, span=span, \
							dataLs=[])
				if segmentKey not in self.segmentKey2dataLsRBDict:
					self.segmentKey2dataLsRBDict[segmentKey] = []
				self.segmentKey2dataLsRBDict[segmentKey].append(intervalData)
				counter += 1
		sys.stderr.write("%s intervals in segmentKey2dataLsRBDict %s.\n"%(counter, self.segmentKey2dataLsRBDict))
		return self.segmentKey2dataLsRBDict
	
	def processRow(self, row=None, pdata=None):
		"""
		2013.07.31
		"""
		returnValue = 0
		col_name2index = getattr(pdata, 'col_name2index', None)
		y_ls = getattr(pdata, 'y_ls', None)
		if col_name2index and y_ls is not None:
			chromosomeIndex = col_name2index.get(self.chromosomeHeader, None)
			positionIndex = col_name2index.get(self.positionHeader, None)
			
			if self.whichColumnHeader:
				whichColumn = col_name2index.get(self.whichColumnHeader, None)
			elif self.whichColumn:
				whichColumn = self.whichColumn
			else:
				whichColumn = None
			if whichColumn is not None:
				yValue = row[whichColumn]
				if yValue not in self.missingDataNotation:
					yValue = self.processValue(yValue, processType=self.logY, valueForNonPositiveValue=self.valueForNonPositiveYValue)
				chromosome = row[chromosomeIndex]
				position = float(row[positionIndex])
				segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=chromosome, \
							span_ls=[position, position], \
							min_reciprocal_overlap=0.0000001, )
				node_ls = []
				self.segmentKey2dataLsRBDict.findNodes(segmentKey, node_ls=node_ls)
				for node in node_ls:
					for intervalData in node.value:
						intervalData.dataLs.append(yValue)
				returnValue = 1
		return returnValue
	
	def processHeader(self, header=None, pdata=None, rowDefinition=None):
		"""
		2013.07.31
			override this to output custom header
		"""
		header = ["chromosome", "start", "end", "noOfEntries", "score"]
		self._writeHeader(header=header, pdata=pdata, rowDefinition=rowDefinition)
		
	
	def reduce(self, **keywords):
		"""
		2012.10.15
			run after all files have been walked through
		"""
		for node in self.segmentKey2dataLsRBDict:
			for oneData in node.value:
				if len(oneData.dataLs)>0:
					reduce_value = self.reduceFunction(oneData.dataLs)
					data_row=  [node.key.chromosome, node.key.start, node.key.stop, len(oneData.dataLs), reduce_value]
					self.writer.writerow(data_row)
		
	

if __name__ == '__main__':
	main_class = GenomeMovingAverageStatistics
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()