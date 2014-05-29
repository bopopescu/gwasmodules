#!/usr/bin/env python
"""
Examples:
	%s -i -o 
	
	%s 
	
	%s 
	
Description:
	2013.11.24 Input is a coordinate map file for loci (before and after LiftOver).
		This program computes probability that each locus is correctly mapped to the target coordinate system.
		
		Output:
		1. main output, probability for each SNP being correctly mapped.
			oldChromosome, oldStart, oldStop, newChromosome, newStart, newStop, mapPvalue
		2. side output, delta of old and new intervals.
			oldChromosome, oldStart, oldStop, newChromosome, newStart, newStop, intervalDelta.

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, MatrixFile, PassingData
from pymodule.yhio.AbstractMatrixFileWalker import AbstractMatrixFileWalker
import numpy
from scipy.stats import norm

parentClass = AbstractMatrixFileWalker
class ComputeLiftOverLocusProbability(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update({
						('locusIntervalDeltaOutputFname', 1, ): ['', '', 1, 'file that would contain delta of intervals from old and new coordinate system. \
	Used to check if normal distribution on each chromosome. Output format: oldChromosome, oldStart, oldStop, newChromosome, newStart, newStop, intervalDelta.', ],\
						('startPosition', 0, int):[None, '', 1, 'probability for loci whose start positions are bigger than this argument would be computed.\
	Model parameters are estimated using all input data. This argument is used to avoid edge/boundary effect.'],\
						('stopPosition', 0, int):[None, '', 1, 'probability for loci whose stop positions are less than this argument would be computed.\
	Model parameters are estimated using all input data. This argument is used to avoid edge/boundary effect.'],\
						
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		parentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)

	def setup(self, **keywords):
		"""
		noOfTotalIntervals = 0
		noOfCrossChromosomeIntervals = 0
		
		targetChromosome 2 mapData
			intervalDeltaList	=> median
			orientation  (queryStrand)
				0=forward
				1=backward
			mean	=> using 80% of data (sort the delta list, then take 10% to 90% of the list)
			stddev	=> if stddev is zero, use 1.
		
		locusKey (oldChromosome, oldStart, oldStop) 2 mapData
			targetCoordinate (newChromosome, newStart, newStop).
			leftIntervalDelta: None = boundary
			rightIntervalDelta: None = boundary, 10E10 = cross chromosome
			
			probability: max( P(SNP_i_left_interval), P(SNP_i_right_interval)).
				P(interval):
					If one interval is on the same chromosome,  P(target-chromosome)*P(interval delta size)
					If not, P(chromosome-cross event). 
			
		Not implemented: for a whole genome input (rather than a window),
			an RBTree of windows should be used to counter regional effect.
		
		2013.11.24
			run before anything is run
		"""
		AbstractMatrixFileWalker.setup(self, **keywords)
		
		self.noOfTotalIntervals = 0.0
		self.noOfCrossChromosomeIntervals = 0.0 #make it float for division
		
		self.targetChromosome2mapData = {}
		self.locusKey2mapData = {}
		self.previousLocusData = None
		
		#write header for the main output
		header = ['oldChromosome', 'oldStart', 'oldStop', 'oldStrand', 'newChromosome', 'newStart', 'newStop', 'mapPvalue']
		self.writer.writerow(header)
		self.invariantPData.headerOutputted = True	#avoid double header output
		
		#open the other writer and write header
		self.sideOutput = MatrixFile(self.locusIntervalDeltaOutputFname, openMode='w', delimiter='\t')
		header = ['oldChromosome', 'oldStart', 'oldStop', 'oldStrand', 'newChromosome', 'newStart', 'newStop', 'intervalDelta']
		self.sideOutput.writeHeader(header)
		

	def processRow(self, row=None, pdata=None):
		"""
		2012.10.7
		"""
		returnValue = 1
		self.col_name2index = getattr(pdata, 'col_name2index', None)
		queryStrandIndex = self.col_name2index.get("queryStrand")
		
		queryChromosomeIndex = self.col_name2index.get("queryChromosome")
		queryStartIndex = self.col_name2index.get("queryStart")
		queryStopIndex = self.col_name2index.get("queryStop")
		
		newChrIndex = self.col_name2index.get("newChr")
		newRefStartIndex = self.col_name2index.get("newRefStart")
		newRefStopIndex = self.col_name2index.get("newRefStop")
		
		queryStrand = row[queryStrandIndex]
		queryChromosome = row[queryChromosomeIndex]
		queryStart = int(row[queryStartIndex])
		queryStop = int(row[queryStopIndex])
		
		newChr = row[newChrIndex]
		newRefStart = int(row[newRefStartIndex])
		newRefStop = int(row[newRefStopIndex])
		
		#create current locus data
		locusKey = (queryChromosome, queryStart, queryStop)
		currentLocusData = PassingData(locusKey=locusKey, queryStrand=queryStrand, queryChromosome=queryChromosome,\
					queryStart=queryStart, queryStop=queryStop, \
					newChr=newChr, newRefStart=newRefStart, newRefStop=newRefStop)
		
		#insert entry into locusKey2mapData
		self.locusKey2mapData[locusKey] = PassingData(locusData = currentLocusData, leftIntervalDelta=None,\
										rightIntervalDelta=None, mapProbability=None)
		if self.previousLocusData is not None:
			#calculate interval delta
			if self.previousLocusData.newChr != currentLocusData.newChr:
				intervalDelta = 10E10
				self.noOfCrossChromosomeIntervals += 1
			else:
				querySpan = currentLocusData.queryStart - currentLocusData.queryStop
				targetSpan = currentLocusData.newRefStart - currentLocusData.newRefStop
				if queryStrand=='+':
					intervalDelta = targetSpan - querySpan
				else:
					intervalDelta = targetSpan + querySpan
				# insert it into self.targetChromosome2mapData
				if currentLocusData.newChr not in self.targetChromosome2mapData:
					self.targetChromosome2mapData[currentLocusData.newChr] = PassingData(intervalDeltaList=[],\
																			orientation=queryStrand,\
																			mean=None,\
																			stddev=None,\
																			probability=None)
				self.targetChromosome2mapData[currentLocusData.newChr].intervalDeltaList.append(intervalDelta)
			
			#output to the side
			self.sideOutput.writerow([currentLocusData.queryChromosome,\
					currentLocusData.queryStart, currentLocusData.queryStop, currentLocusData.queryStrand, \
					currentLocusData.newChr, currentLocusData.newRefStart, currentLocusData.newRefStop, intervalDelta])
			
			#assign it as right interval delta of previous locus
			self.locusKey2mapData[self.previousLocusData.locusKey].rightIntervalDelta = intervalDelta
			
			# assign it as left interval delta of current locus.
			self.locusKey2mapData[locusKey].leftIntervalDelta = intervalDelta
			
			self.noOfTotalIntervals +=1
			
		self.previousLocusData = currentLocusData
		return returnValue
	
	def calculateLocusMapProbabilityGivenIntervalDelta(self, intervalDelta=None, targetChromosomeMapData=None, crossChromosomeProbability=None):
		"""
		2013.11.25
		"""
		mapProbability = 1
		if intervalDelta is not None:
			if intervalDelta==10E10:
				mapProbability *= crossChromosomeProbability
			else:
				lessThanGivenValueProb = norm.cdf(intervalDelta, loc=targetChromosomeMapData.mean, scale=targetChromosomeMapData.stddev)
				if intervalDelta>targetChromosomeMapData.mean:	#two-sided p-value
					deltaProb = 2*(1- lessThanGivenValueProb)
				else:
					deltaProb = 2*lessThanGivenValueProb
				mapProbability *= targetChromosomeMapData.probability * deltaProb
				
		return mapProbability
	
	def reduce(self, **keywords):
		"""
		2012.10.15
			run after all files have been walked through
		"""
		counter = 0
		real_counter = 0
		
		locusKeyList = self.locusKey2mapData.keys()
		locusKeyList.sort()
		
		sys.stderr.write("%s target chromosomes, %s cross-chromosome intervals, %s total intervals .\n "%\
						(len(self.targetChromosome2mapData), self.noOfCrossChromosomeIntervals, self.noOfTotalIntervals))
		
		if self.noOfTotalIntervals>0:
			sys.stderr.write("Running estimates for each target chromosome ... ")
			#estimates for each chromosome
			self.crossChromosomeProbability = float(self.noOfCrossChromosomeIntervals)/self.noOfTotalIntervals
			for targetChromosome in self.targetChromosome2mapData:
				mapData = self.targetChromosome2mapData.get(targetChromosome)
				#overall probability for an interval to be on this chromosome
				if len(mapData.intervalDeltaList)==0:	#just one crossing event
					mapData.probability = 1/float(self.noOfTotalIntervals)
				else:
					mapData.probability = len(mapData.intervalDeltaList)/float(self.noOfTotalIntervals)
				#estimate mean and stddev
				mapData.intervalDeltaList.sort()
				startIndex = max(0, int(len(mapData.intervalDeltaList)*0.1))
				stopIndex = max(int(len(mapData.intervalDeltaList)*0.9)+1, 1)
				if startIndex>=stopIndex:
					stopIndex = startIndex + 1
				robustDataList = mapData.intervalDeltaList[startIndex:stopIndex]
				
				stddev = 1
				if len(robustDataList)>0:
					mapData.mean = numpy.mean(robustDataList)
					if len(robustDataList)>1:
						stddev = numpy.std(robustDataList)
				else:
					mapData.mean = 0
				if stddev == 0:
					stddev = 1
				mapData.stddev = stddev
			sys.stderr.write(".\n")
			
			#output
			sys.stderr.write("Output %s SNPs with map p-value ..."%(len(locusKeyList)))
			for locusKey in locusKeyList:
				counter += 1
				locusMapData = self.locusKey2mapData.get(locusKey)
				locusData = locusMapData.locusData
				if locusMapData.leftIntervalDelta!=None:
					leftProbability = self.calculateLocusMapProbabilityGivenIntervalDelta(intervalDelta=locusMapData.leftIntervalDelta, \
														targetChromosomeMapData=self.targetChromosome2mapData.get(locusData.newChr),\
														crossChromosomeProbability=self.crossChromosomeProbability)
				else:
					leftProbability = 0
				if locusMapData.rightIntervalDelta!=None:
					rightProbability = self.calculateLocusMapProbabilityGivenIntervalDelta(intervalDelta=locusMapData.rightIntervalDelta, \
												targetChromosomeMapData=self.targetChromosome2mapData.get(locusData.newChr),\
												crossChromosomeProbability=self.crossChromosomeProbability)
				else:
					rightProbability = 0
				mapProbability = max(leftProbability, rightProbability)
				data_row = [locusData.queryChromosome,\
					locusData.queryStart, locusData.queryStop, locusData.queryStrand, \
					locusData.newChr, locusData.newRefStart, locusData.newRefStop, mapProbability]
				self.writer.writerow(data_row)
				real_counter += 1
			sys.stderr.write("\n")
		else:	#single SNP (give a low probability)
			sys.stderr.write("Zero intervals, output %s SNPs with 0.001 map p-value ..."%(len(locusKeyList)))
			for locusKey in locusKeyList:
				counter += 1
				locusMapData = self.locusKey2mapData.get(locusKey)
				locusData = locusMapData.locusData
				mapProbability = 0.001
				data_row = [locusData.queryChromosome,\
					locusData.queryStart, locusData.queryStop, locusData.queryStrand, \
					locusData.newChr, locusData.newRefStart, locusData.newRefStop, mapProbability]
				self.writer.writerow(data_row)
				real_counter += 1
			sys.stderr.write("\n")
		
		if counter >0:
			fraction = float(real_counter)/float(counter)
		else:
			fraction = -1
		sys.stderr.write("%s/%s (%.3f) outputted.\n"%(real_counter, counter, fraction))
		
		self.sideOutput.close()
		#close the self.invariantPData.writer
		AbstractMatrixFileWalker.reduce(self, **keywords)
	
if __name__ == '__main__':
	main_class = ComputeLiftOverLocusProbability
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()