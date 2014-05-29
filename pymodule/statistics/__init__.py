#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i /tmp/Contig315_StKitts_vs_Nevis.tsv --xColumnHeader=StKitts --whichColumnHeader=Nevis
		-s 1.0 -o /tmp/Contig315_StKitts_vs_Nevis.2D.png
	

Description:
	2012.10.12
		program to estimate how many outliers off the y=x axis.
		1. hard cutoff. abs(y-x)<=minDelta
		2. model y-x as a normal distribution, estimate its mean/variance
			then add them up as chi-squared statistic.
	If "-i ..." is given, it is regarded as one of the input files (plus the ones in trailing arguments). 
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:	   #64bit
#	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
#	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
#else:   #32bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import numpy, random
from pymodule import getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter
from pymodule import yh_matplotlib

def calculateChiSqStatOfDeltaVector(dataVector=None, mean=None, std=None):
	"""
	2012.10.14
		adapted from vervet/src/pedigree/DetectWrongLabelByCompKinshipVsIBD.DetectWrongLabelByCompKinshipVsIBD.calculateChiSqStatOfDeltaVector()
	2012.8.22
		fisher's method in combining pvalue of z-score (normalized abs(kinship-ibd)).
	"""
	import rpy
	medianAbsDelta = None
	noOfNonMissing = None
	chiSqStat = 0
	noOfNonMissing = 0
	chiSqMinusLogPvalue = None
	for i in xrange(len(dataVector)):
		delta = dataVector[i]
		if delta is not None:
			Z = abs((delta-mean)/std)	#2012.8.23 take absolute value, since it's always P(X>a), negative z-score gets wrong portion.
			logPvalue = rpy.r.pnorm(Z, lower_tail = rpy.r.FALSE, log=rpy.r.TRUE)	#the latter log is natural log.
			#should use two-tail, rather than one-tail
			logPvalue += math.log(2)
			
			#pvalue = utils.getZScorePvalue(Z)	#numerical underflow, pvalue=0
			chiSqStat += -2*logPvalue	#should use natural log
			noOfNonMissing += 1
	if noOfNonMissing>0:
		chiSqMinusLogPvalue = -rpy.r.pchisq(chiSqStat, 2*noOfNonMissing, lower_tail = rpy.r.FALSE, log=rpy.r.TRUE)
		#chiSqPvalue = stats.chisqprob(chiSqStat, df=2*noOfNonMissing)
		#equivalent to = 1- stats.chi2.cdf(chiSqStat, df) = Prob (x<chiSqStat)
	return PassingData(chiSqStat=chiSqStat, noOfNonMissing=noOfNonMissing, chiSqMinusLogPvalue=chiSqMinusLogPvalue)

def estimateMeanStdFromData(dataVector=None, excludeTopFraction=0.2):
	"""
	2012.10.14
		adapted from vervet/src/pedigree/DetectWrongLabelByCompKinshipVsIBD.DetectWrongLabelByCompKinshipVsIBD.estimateAbsDeltaMeanStd()
	2012.8.22
	"""
	sys.stderr.write("Estimating mean&std using the middle %.1f%% of data (n=%s) ..."%\
					((1-excludeTopFraction)*100, len(dataVector)))
	noOfRows = len(dataVector)
	import numpy
	# 2012.8.22 draw some histogram to check what data looks like
#		if len(dataVector)>10:
#			outputFname = '%s_kinship_ibd_hist.png'%(self.outputFnamePrefix)
#			yh_matplotlib.drawHist(dataVector, title='', \
#							xlabel_1D="kinship-ibd", xticks=None, \
#							outputFname=outputFname, min_no_of_data_points=10, \
#							needLog=True, \
#							dpi=200, min_no_of_bins=25)
	#dataVector = map(abs, dataVector)	#2012.8.23 no abs
	dataVector.sort()
	startIndex = min(0, int(len(dataVector)*(excludeTopFraction/2))-1)
	stopIndex = int(len(dataVector)*(1-excludeTopFraction/2))
	dataVector = dataVector[startIndex:stopIndex]
	
	data_mean = numpy.mean(dataVector)
	data_std = numpy.std(dataVector)
	
	sys.stderr.write(" mean=%.3f, std=%.3f.\n"%(data_mean, data_std))
	return PassingData(mean=data_mean, std=data_std)


def getZScorePvalue(zscore=None, twoSided=False):
	"""
	2012.10.15
		was in pymodule/utils.py
	2012.8.22
		becasue this import wouldn't work. hard to remember:
	
			>>> import scipy
			>>> scipy.stats.norm.sf
			Traceback (most recent call last):
			  File "<stdin>", line 1, in <module>
			AttributeError: 'module' object has no attribute 'stats'
		
		zscore could also be a vector (list)
	"""
	import scipy.stats as stats
	pvalue = stats.norm.sf(zscore)
	if twoSided:
		pvalue  = pvalue* 2
	return pvalue

def reOrderListOfListByFirstMember(listOfList=None):
	"""
	2012.10.29
		this function reorders every list in listOfList in ascending order using values in the 1st list of listOfList.
		And every member of listOfList will be turned into numpy.float array.
	2012.10.25
	"""
	firstList = listOfList[0]
	x_ar = numpy.array(firstList, numpy.float)
	#sort x_ar and y_ar must be in the order of x_ar
	indexOfOrderList = numpy.argsort(x_ar)
	returnListOfList = []
	for ls in listOfList:
		ar = numpy.array(ls, numpy.float)
		ar = ar[indexOfOrderList]
		returnListOfList.append(ar)
	return PassingData(listOfList=returnListOfList)

def splineFit(x_ls=None, y_ls=None, no_of_steps=100, needReorderData=True):
	"""
	2012.10.25
	http://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html
	class scipy.interpolate.UnivariateSpline(x, y, w=None, bbox=[None, None], k=3, s=None)[source]
		One-dimensional smoothing spline fit to a given set of data points.
		
		Fits a spline y=s(x) of degree k to the provided x, y data. s specifies the number of knots by specifying a smoothing condition.
		
		Parameters :	
		x : array_like
			1-D array of independent input data. Must be increasing.
		y : array_like
			1-D array of dependent input data, of the same length as x.
		w : array_like, optional
			Weights for spline fitting. Must be positive. If None (default), weights are all equal.
			bbox : array_like, optional
			2-sequence specifying the boundary of the approximation interval. If None (default), bbox=[x[0], x[-1]].
		k : int, optional
			Degree of the smoothing spline. Must be <= 5.
		s : float or None, optional
			Positive smoothing factor used to choose the number of knots. Number of knots will be increased until the smoothing condition is satisfied:
			sum((w[i]*(y[i]-s(x[i])))**2,axis=0) <= s
			If None (default), s=len(w) which should be a good value if 1/w[i] is an estimate of the standard deviation of y[i].
			If 0, spline will interpolate through all data points.

	"""
	import scipy.interpolate
	if needReorderData:
		sortData = reOrderListOfListByFirstMember(listOfList=[x_ls, y_ls])
		x_ls = sortData.listOfList[0]
		y_ls = sortData.listOfList[1]
	
	sp = scipy.interpolate.UnivariateSpline(x_ls, y_ls,k=5)
	step = (x_ls[-1]-x_ls[0])/float(no_of_steps)
	n_x_ls = numpy.arange(x_ls[0], x_ls[-1], step)
	n_y_ls = map(sp, n_x_ls)
	return PassingData(x_ls=n_x_ls, y_ls=n_y_ls)

def movingAverage(listOfList=None, no_of_steps=100, needReorderData=True, reduceType=1, minNoOfTotal=10,\
				minValueForFraction=0.8):
	"""
	2012.10.29
		a function that calculates moving average of every list in listOfList.
		stepping is applied to the first member in listOfList.
		all member lists are of the same length.
		reduceType
			1: median
			2: mean
			3: fraction that is >=minValueForFraction, only applied to lists after the first one, numpy.median for X
	2012.10.26
	"""
	if needReorderData:
		sortData = reOrderListOfListByFirstMember(listOfList=listOfList)
		listOfList = sortData.listOfList
	firstList = sortData.listOfList[0]
	step = (firstList[-1]-firstList[0])/float(no_of_steps)
	stepIndex2Data = {}
	noOfLists = len(listOfList)
	for j in xrange(len(firstList)):
		x = firstList[j]
		stepIndex = int(x/step)	#figure which bracket/bag all the data from this column should fall into.
		if stepIndex not in stepIndex2Data:
			stepIndex2Data[stepIndex] = PassingData(listOfList=[])
			for i in xrange(noOfLists):
				stepIndex2Data[stepIndex].listOfList.append([])
		#y = y_ls[j]
		#stepIndex2Data[stepIndex].listOfList[0].append(x)
		for i in xrange(noOfLists):
			valueAtThatList = listOfList[i][j]
			stepIndex2Data[stepIndex].listOfList[i].append(valueAtThatList)
	
	stepIndexList = stepIndex2Data.keys()
	stepIndexList.sort()
	
	fractionFunction = lambda ls: sum([a>=minValueForFraction for a in ls])/float(len(ls))
	reduceType2Function = {1: numpy.median, 2: numpy.mean, 3: fractionFunction}
	reduceFunction = reduceType2Function.get(reduceType, numpy.median)
	
	n_x_ls = []
	returnListOfList = []
	for i in xrange(noOfLists):
		returnListOfList.append([])
	
	for stepIndex in stepIndexList:
		data = stepIndex2Data.get(stepIndex)
		subListOfList = data.listOfList
		if len(subListOfList[0])<minNoOfTotal:
			continue
		
		for i in xrange(noOfLists):
			if i==0 and reduceType==3:
				_reduceFunction = numpy.median
			else:
				_reduceFunction = reduceFunction
			
			reduce_value = _reduceFunction(subListOfList[i])
			returnListOfList[i].append(reduce_value)
	return PassingData(listOfList=returnListOfList)


def calculate7NumberSummaryForOneList(ls, returnObj=None):
	"""
	2012.10.29 moved from variation/src/common.py
	2009-11-25
		calculate a 7-number summary stats for a given list
	"""
	from scipy import stats	# for scoreatpercentile/percentileatscore to get quartiles
	if returnObj is None:
		returnObj = PassingData()
	returnObj.minimum = numpy.min(ls)
	returnObj.first_decile = stats.scoreatpercentile(ls, 10)
	returnObj.lower_quartile = stats.scoreatpercentile(ls, 25)
	returnObj.median = stats.scoreatpercentile(ls, 50)
	returnObj.upper_quartile = stats.scoreatpercentile(ls, 75)
	returnObj.last_decile = stats.scoreatpercentile(ls, 90)
	returnObj.maximum = numpy.max(ls)
	return returnObj

class NumberContainer(object):
	"""
	#2013.05.24 a structure to scale/normalize any range/series of numbers
		to make them [0,1]
	"""
	def __init__(self, minValue=None, maxValue=None):
		self.minValue = minValue
		self.maxValue = maxValue
		self.valueList = []
	
	def addOneValue(self, value=None):
		if value is not None:
			self.valueList.append(value)
			if self.minValue is None or self.minValue > value:
				self.minValue = value
			if self.maxValue is None or self.maxValue < value:
				self.maxValue = value
	
	def normalizeValue(self, value=None):
		"""
		2013.05.24
		"""
		if self.minValue is not None and self.maxValue is not None:
			return (value-self.minValue)/float(self.maxValue-self.minValue)
		else:
			raise

class DiscreteProbabilityMassContainer(object):
	"""
	Examples:
		probabilityMassContainer = DiscreteProbabilityMassContainer(object2proabilityMassDict=self.originalIndividualID2representativeData)
		sampledIndividualID = probabilityMassContainer.sampleObject()
	
	2013.05.26
		function to do sampling
	"""
	def __init__(self, object2proabilityMassDict=None):
		"""
		2013.05.26
		"""
		from pymodule.algorithm.RBTree import RBDict
		self.rbDict = RBDict()
		self.totalProbabilityMass = 1	#default
		if object2proabilityMassDict is not None:
			self._constructFromDiscreteProbabilityMassDict(dc=object2proabilityMassDict)
	
	def _constructFromDiscreteProbabilityMassDict(self, dc=None):
		"""
		2013.05.28
			dc is a structure with object name as key, and object probability mass (normalized or not) as value. i.e.
				{"1978001":0.5, "1980001":1.5}
				
			argument probabilityNormalized: whether the sum of all values in dc adds up to 1.
		"""
		from pymodule.yhio.CNV import CNVSegmentBinarySearchTreeKey
		startProbMass = 0.0
		for discreteVariable, probabilityMass in dc.iteritems():
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome="1", span_ls=[startProbMass, startProbMass+probabilityMass], \
												min_reciprocal_overlap=0.001, isDataDiscrete=False)
						#min_reciprocal_overlap=1: must be complete overlap in order for two objects occupying same key
			self.rbDict[segmentKey] = discreteVariable
			startProbMass += probabilityMass
		self.totalProbabilityMass = startProbMass
		sys.stderr.write("%s\n"%(repr(self.rbDict)))
	
	def sampleObject(self):
		"""
		"""
		from pymodule.yhio.CNV import CNVSegmentBinarySearchTreeKey
		u = random.random()*self.totalProbabilityMass
		key = CNVSegmentBinarySearchTreeKey(chromosome="1", span_ls=[u], \
												min_reciprocal_overlap=0.0000001)
		#randint.(0,noOfTotalRows-1)
		
		node = self.rbDict.findNode(key)
		if node:
			return node.value
		else:
			return None
	
		