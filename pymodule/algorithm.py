#!/usr/bin/env python
"""
2009-11-2
	module for various algorithms
"""

import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from ProcessOptions import  ProcessOptions
import copy

def listSubsets(element_ls, subset_size=None):
	"""
	Example:
		element_ls = [1,6,7]
		listSubsets(element_ls, subset_size=1)
	
		listSubsets(element_ls, subset_size=2)
	
		listSubsets(element_ls)
	
	2009-11-2
		list all possible subsets of element_ls. If subset_size is given, only produce subsets of that size.
			Otherwise, all of them.
		It's backtracking algorithm at play.
	"""
	sys.stderr.write("Listing all subsets of a list of %s elements ... "%(len(element_ls)))
	no_of_elements = len(element_ls)
	candidate_ls_ls = [[0,1]]	# 0 or 1 in k-th candidate_ls_ls signifies whether element_ls[k] would be included or not.
	#for i in range(no_of_elements):
	#	candidate_ls_ls.append([0, 1])
	
	one_solution=[]
	solution_ls=[]
	k = 0
	while k>=0:
		while len(candidate_ls_ls[k])>0:
			next_element = candidate_ls_ls[k].pop()
			if len(one_solution)==no_of_elements:
				one_solution[k] = next_element
			else:
				one_solution.append(next_element)
			k += 1
			if k==no_of_elements:
				if subset_size is None or sum(one_solution)==subset_size:
					one_subset = []
					for i in range(no_of_elements):
						if one_solution[i]==1:
							one_subset.append(element_ls[i])
					solution_ls.append(one_subset)		# python data by reference. have to use copy to sever the tie.
				break
			if len(candidate_ls_ls)<=k:	# still growing
				candidate_ls_ls.append([0,1])
			else:	# fully grown, now just replace the candidate list
				candidate_ls_ls[k] = [0, 1]
		k -= 1
	sys.stderr.write("Done.\n")
	return solution_ls

def simpleCandidateSetGenerator(k, element_ls=range(2), max_solution_size=5):
	"""
	2010-4-21
		candidates will be from element_ls as far as k is less than max_solution_size.
		used as argument candidateSetGenerator by enumerateCombinations().
	"""
	if k<max_solution_size:
		return element_ls[:]	#[:] is a must. otherwise, range(2) is used by reference.
	else:
		return []


def enumerateCombinations(candidateSetGenerator=None):
	"""
	Examples:
		# 2010-4-21 directly use simpleCandidateSetGenerator
		enumerateCombinations(simpleCandidateSetGenerator)
		
		# 2010-4-21 modify the max_solution_size to be 4
		f = lambda x: simpleCandidateSetGenerator(x, element_ls=range(2), max_solution_size=4)
		enumerateCombinations(f)
	
	2010-4-21
		This function enumerates combinations given candidate sets at each position.
		The backtrack algorithm will stop when the candidateSetGenerator returns nothing.
		candidateSetGenerator accepts argument k, index of the element in one solution.
	"""
	sys.stderr.write("Enumerating combinations of candidates ...")
	candidate_set_ls = [candidateSetGenerator(0)]	# first candidate
	
	one_solution=[]
	solution_ls = []
	k = 0
	while k>=0:
		while len(candidate_set_ls[k])>0:
			next_element = candidate_set_ls[k].pop()
			if len(one_solution)>k:
				one_solution[k] = next_element
			else:
				one_solution.append(next_element)
			k += 1
			candidate_set = candidateSetGenerator(k)
			if not candidate_set:	# no more candidate set, this is a solution.
				solution_ls.append(one_solution[:])	# [:] is a must. otherwise, Every single item in solution_ls will be same.
				# because each one is a reference to one_solution. 
			if len(candidate_set_ls)<=k:	# still growing
				candidate_set_ls.append(candidate_set)
			else:	# fully grown, now just replace the candidate list
				candidate_set_ls[k] = candidate_set
		k -= 1
	sys.stderr.write("Done.\n")
	return solution_ls

def ltsFit(x_ls, y_ls, fractionUsed=0.6, startX=1, stopX=5):
	"""
	2010-6-1
		solve the computing node hang-up (I/O stuck) issue by adding these:
			import ROOT
			try:	# 2010-5-31 old version (5.18.0) doesn't have IgnoreCommandLineOptions.
				ROOT.PyConfig.IgnoreCommandLineOptions = True	#otherwise
				# Warning in <TApplication::GetOptions>: file <output file by -o > has size 0, skipping
			except:
				pass
			try:	# 2010-5-31  disable .StartGuiThread
				ROOT.PyConfig.StartGuiThread = 0
			except:
				pass
	2010-5-30
		return chiSquare as well
	2010-5-21
		use ROOT to do least trimmed square (LTS) fitting:
			fit the y=a+bx with trimming fraction = 1-fractionUsed.
	
	Example:
	
	import numpy
	x_ls = numpy.array(range(100), numpy.float)
	y_ls = x_ls/2.
	for i in range(len(y_ls)):
		import random
		new_y = random.random()-0.5
		y_ls[i] += new_y
	
	# mess up some portion of y
	for i in range(5):
		import random
		new_y = random.random()
		new_y_index = random.sample(range(100),1)
		y_ls[new_y_index[0]] = new_y
	import numpy
	x_ls = numpy.array([ 2.64884758,  3.51235008,  2.83090925,  3.41229248,  3.01451969,\
    2.49899888,  3.69988108,  2.74896216,  3.05307841,  3.75705409,\
    3.08653784,  3.10703993,  3.61071348,  3.21285319,  2.91460752,\
    3.53737831,  3.06333303,  3.35391617,  3.43568516,  3.34429312,\
    3.31576061,  2.8007164 ,  2.73639655,  3.14690256,  3.10174704,\
    2.80888581,  2.72754121,  2.90064001,  3.19270658,  3.50596333,\
    2.61804676,  3.18127131,  3.27542663,  3.09586573], dtype=numpy.float32)	# numpy.float32 is not supported by ROOT
	y_ls = numpy.array([ 2.52827311,  3.27265358,  2.36172366,  2.95760489,  2.50920248,\
    2.3443923 ,  3.23502254,  2.35410833,  2.50582743,  2.48501062,\
    2.82510138,  2.70799541,  2.43136382,  2.76342535,  2.45178652,\
    3.08224201,  2.26481771,  2.7387805 ,  3.23274207,  2.82769203,\
    2.25042009,  2.56702638,  2.4082365 ,  2.44793224,  2.65127802,\
    2.57460976,  2.43136382,  2.39005065,  2.70027065,  3.04452848,\
    2.28555727,  2.71933126,  2.6468935 ,  2.54157925], dtype=numpy.float32)
    
	fit_y_ls = ltsFit(x_ls, y_ls)
	
	import pylab
	pylab.plot(x_ls, y_ls, '.')
	pylab.plot(x_ls, fit_y_ls, '.')
	pylab.legend(['raw data','fitted'])
	pylab.show()
	sys.exit(0)
	
	"""
	import ROOT
	try:	# 2010-5-31 old version (5.18.0) doesn't have IgnoreCommandLineOptions.
		ROOT.PyConfig.IgnoreCommandLineOptions = True	#otherwise
		# Warning in <TApplication::GetOptions>: file <output file by -o > has size 0, skipping
	except:
		pass
	try:	# 2010-5-31  disable .StartGuiThread
		ROOT.PyConfig.StartGuiThread = 0
	except:
		pass
	
	#ROOT.gROOT.Reset()	# 2010-5-31 dont' know what this is  for.
	ROOT.gROOT.SetBatch(True)	#to avoid interative mode (drawing canvas and etc.)	
	from ROOT import TFormula, TF1, TGraph
	import numpy
	lm = TF1('lm', 'pol1', startX, stopX)	#[0]+[1]*x is essentially same as pol1 but option rob in Fit() only works with pol1.
	#ROOT is very dtype-sensitive. numpy.float32 won't work.
	if hasattr(x_ls, 'dtype') and x_ls.dtype==numpy.float:
		pass
	else:
		sys.stderr.write('converting x_ls')
		x_ls = numpy.array(x_ls, dtype=numpy.float)
		sys.stderr.write(".\n")
	if hasattr(y_ls, 'dtype') and y_ls.dtype==numpy.float:
		pass
	else:
		sys.stderr.write('converting y_ls')
		y_ls = numpy.array(y_ls, dtype=numpy.float)
		sys.stderr.write(".\n")
	gr = TGraph(len(x_ls), x_ls, y_ls)
	gr.Fit(lm, "+rob=%s"%fractionUsed)
	fit = gr.GetFunction('lm')
	chiSquare = fit.GetChisquare();
	fit_y_ls = []
	for x in x_ls:
		fit_y_ls.append(fit.Eval(x))
	from utils import PassingData
	return PassingData(fit_y_ls=fit_y_ls, chiSquare=chiSquare)

if __name__ == '__main__':
	#import pdb
	#pdb.set_trace()
	
	element_ls = [1,6,7]
	print listSubsets(element_ls, subset_size=1)
	
	print listSubsets(element_ls, subset_size=2)
	
	print listSubsets(element_ls)