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

from pymodule import getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter
from pymodule import yh_matplotlib
import numpy, random

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