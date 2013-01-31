#!/usr/bin/env python
"""
2013.1.27 coordinates are all 1-based and inclusive (??), start=1,stop=100 means [1,100]?? 
2009-10-31
	module for CNV (copy-number-variation)-related functions & classes
"""

import os, sys, math
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import fileinput
import numpy
from pymodule.ProcessOptions import ProcessOptions
from pymodule.yhio.SNP import GenomeWideResult, DataObject
from pymodule.utils import getColName2IndexFromHeader, dict_map, importNumericArray, figureOutDelimiter, PassingData

def get_overlap_ratio(span1_ls, span2_ls):
	"""
	2013.1.27 add total_span, overlapFraction, and return PassingData (data structure)
	2012.5.17 bugfix. overlap_length can't be negative (hard-set to zero).
	2010-8-18
		figure out the overlap coordinates
	2010-8-2
		add +1 in "max(0, segment_stop_pos - qc_start+1)" to recognize even 1-bp overlapping
		bug fix: span1_ls and/or span2_ls have same value on start & stop.
	2010-7-30
		doc: if two spans have no overlap, it returns 0,0
	2010-2-11
		swap overlapFraction1 and overlapFraction2 so that they match the span1_ls and span2_ls
	2009-12-13
		calculate the two overlap ratios for two segments
	"""
	segment_start_pos, segment_stop_pos = span1_ls
	qc_start, qc_stop = span2_ls
	#2010-8-18 figure out the overlap coordinates
	overlap_start_pos = max(qc_start, segment_start_pos)
	total_start_pos = min(qc_start, segment_start_pos)
	
	overlap_stop_pos = min(qc_stop, segment_stop_pos)
	total_stop_pos = max(qc_stop, segment_stop_pos)
	
	total_span = abs(total_stop_pos-total_start_pos) + 1	#2013.1.27 add +1 ?
	
	if overlap_start_pos>overlap_stop_pos:
		overlap_start_pos = None
		overlap_stop_pos = None
	overlap_length = max(0, segment_stop_pos - qc_start+1) - max(0, segment_stop_pos - qc_stop) - \
		max(0, segment_start_pos - qc_start)
		# accomodates 6 scenarios
	overlap_length = float(overlap_length)
	if segment_stop_pos == segment_start_pos:	#2010-8-2 bugfix
		if overlap_length>0:
			overlapFraction1 = 1.0
		else:
			overlap_length = 0	#2012.5.17 bugfix. can't be negative
			overlapFraction1 = 0.0
	else:
		overlapFraction1 = overlap_length/(segment_stop_pos-segment_start_pos+1)	#2013.1.27 add +1?
	if qc_stop == qc_start:
		if overlap_length>0:
			overlapFraction2 = 1.0
		else:
			overlap_length = 0#2012.5.17 bugfix. can't be negative
			overlapFraction2 = 0.0
	else:
		overlapFraction2 = overlap_length/(qc_stop-qc_start+1)	#2013.1.27 add +1 ?
	overlapFraction = float(overlap_length)/total_span	#2013.1.27
	return PassingData(overlapFraction1=overlapFraction1, overlapFraction2=overlapFraction2, \
					overlap_length=overlap_length, \
					overlapFraction=overlapFraction,\
					overlap_start_pos=overlap_start_pos, overlap_stop_pos=overlap_stop_pos, \
					total_start_pos=total_start_pos, total_stop_pos=total_stop_pos,\
					total_span = total_span)

def is_reciprocal_overlap(span1_ls, span2_ls, min_reciprocal_overlap=0.6):
	"""
	2009-12-12
		return True if both overlap ratios are above the min_reciprocal_overlap
	"""
	overlapData = get_overlap_ratio(span1_ls, span2_ls)
	overlapFraction1 = overlapData.overlapFraction1
	overlapFraction2 = overlapData.overlapFraction2
	if overlapFraction1>=min_reciprocal_overlap and overlapFraction2>=min_reciprocal_overlap:
		return True
	else:
		return False

class CNVSegmentBinarySearchTreeKey(object):
	"""
	2010-8-2
		add argument keywords to __init__()
	2009-12-12
		a key designed to represent a CNV segment in the node of a binary search tree (BinarySearchTree.py) or RBTree (RBTree.py),
		
		It has custom comparison function based on the is_reciprocal_overlap() function.
	"""
	def __init__(self, chromosome=None, span_ls=None, min_reciprocal_overlap=0.6, **keywords):
		self.chromosome = chromosome
		
		"""
		if len(span_ls)==2 and span_ls[0]==span_ls[1]:	#2010-8-2 if start and stop are same position, reduce span_ls.
			span_ls = [span_ls[0]]
		"""
		
		self.span_ls = span_ls
		self.min_reciprocal_overlap = min_reciprocal_overlap
		
		self.start = self.span_ls[0]
		if len(self.span_ls)>1:
			self.stop = self.span_ls[1]
			if self.start is not None and self.stop is not None:	#2012.12.31 bugfix
				self.span = abs(self.stop-self.start)
			else:
				self.span = 0
		else:
			self.stop = self.start
			self.span = 0
		if self.stop<self.start:
			sys.stderr.write("Warning: not supposed to happen. stop %s is smaller than start %s.\n"%(self.stop, self.start))
		for key, value in keywords.iteritems():	#2010-8-2
			setattr(self, key, value)
	
	def __lt__(self, other):
		"""
		2009-12-12
			whether self is less than other, in chromosomal order
		"""
		if self.chromosome==other.chromosome:
			if len(self.span_ls)==1:
				if len(other.span_ls)==1:
					return self.span_ls[0]<other.span_ls[0]
				elif len(other.span_ls)>1:
					return self.span_ls[0]<other.span_ls[0]
				else:
					return None
			elif len(self.span_ls)>1:
				if len(other.span_ls)==1:
					return self.span_ls[0]<other.span_ls[0]
				elif len(other.span_ls)>1:
					overlapData = get_overlap_ratio(self.span_ls, other.span_ls)
					overlapFraction1 = overlapData.overlapFraction1
					overlapFraction2 = overlapData.overlapFraction2
					min_overlap = min(overlapFraction1, overlapFraction2)
					if min_overlap>=self.min_reciprocal_overlap:	#should be equal,
						return False
					#if overlapFraction1>0 and overlapFraction1<=1 and  overlapFraction2>0 and overlapFraction2<=1:	# there's overlap between two segments 
					#	return False
					else:
						return self.span_ls[0]<other.span_ls[0]	# whether the stop of this segment is ahead of the start of other
				else:
					return None
		else:
			return self.chromosome<other.chromosome
	
	def __le__(self, other):
		"""
		2009-12-12
			whether self is less than or equal to other, in chromosomal order
		"""
		if self.chromosome==other.chromosome:
			if len(self.span_ls)==1:	# self is a point
				if len(other.span_ls)==1:
					return self.span_ls[0]<=other.span_ls[0]
				elif len(other.span_ls)>1:
					return self.span_ls[0]<=other.span_ls[0]
				else:
					return None
			elif len(self.span_ls)>1:
				if len(other.span_ls)==1:
					return self.span_ls[0]<=other.span_ls[0]
				elif len(other.span_ls)>1:
					is_overlap = is_reciprocal_overlap(self.span_ls, other.span_ls, min_reciprocal_overlap=self.min_reciprocal_overlap)
					if is_overlap:
						return True
					else:
						return self.span_ls[0]<=other.span_ls[0]	# whether the stop of this segment is ahead of the start of other
				else:
					return None
		else:
			return self.chromosome<other.chromosome
		
		
	def __eq__(self, other):
		"""
		2010-2-11
			fix a bug when len(self.span_ls)>1 and len(other.span_ls)==1. it was completely wrong before.
		2009-12-12
		"""
		if self.chromosome==other.chromosome:
			if len(self.span_ls)==1:
				if len(other.span_ls)==1:
					return self.span_ls[0]==other.span_ls[0]
				elif len(other.span_ls)>1:
					return self.span_ls[0]>=other.span_ls[0] and self.span_ls[0]<=other.span_ls[1]	# equal if self is within the "other" segment
				else:
					return None
			elif len(self.span_ls)>1:
				if len(other.span_ls)==1:	# self is a segment. other is a point position.
					return self.span_ls[0]<=other.span_ls[0] and self.span_ls[1]>=other.span_ls[0]	# if self includes the point position, yes it's equal
				elif len(other.span_ls)>1:
					# need to calculate min_reciprocal_overlap
					return is_reciprocal_overlap(self.span_ls, other.span_ls, min_reciprocal_overlap=self.min_reciprocal_overlap)
					#return self.span_ls[1]<other.span_ls[0]	# whether the stop of this segment is ahead of the start of other
				else:
					return None
		else:
			return False
				
	def __ne__(self, other):
		"""
		2009-12-12
			whether self is not equal to other
		"""
		return not self.__eq__(other)
		
	def __ge__(self, other):
		"""
		2009-12-12
		"""
		if self.chromosome==other.chromosome and len(self.span_ls)>1 and len(other.span_ls)>1:
			overlapData = get_overlap_ratio(self.span_ls, other.span_ls)
			overlapFraction1 = overlapData.overlapFraction1
			overlapFraction2 = overlapData.overlapFraction2
			if overlapFraction1>=min_reciprocal_overlap and overlapFraction2>=min_reciprocal_overlap:
				return True
			#elif overlapFraction1>0 and overlapFraction1<=1 and  overlapFraction2>0 and overlapFraction2<=1:	# there's overlap between two segments 
			#	return False
			else:
				return not self.__lt__(other)
		else:
			return not self.__lt__(other)
		
	def __gt__(self, other):
		"""
		2009-12-12
		"""
		if self.chromosome==other.chromosome and len(self.span_ls)>1 and len(other.span_ls)>1:
			overlapData = get_overlap_ratio(self.span_ls, other.span_ls)
			overlapFraction1 = overlapData.overlapFraction1
			overlapFraction2 = overlapData.overlapFraction2
			min_overlap = min(overlapFraction1, overlapFraction2)
			if min_overlap>=self.min_reciprocal_overlap:	#should be equal,
				return False
			#if overlapFraction1>0 and overlapFraction1<=1 and  overlapFraction2>0 and overlapFraction2<=1:	# there's overlap between two segments 
			#	return False
			else:
				return not self.__le__(other)
		else:
			return not self.__le__(other)

	def __str__(self):
		"""
		2009-12-13
		"""
		return "chromosome: %s, span_ls: %s"%(self.chromosome, repr(self.span_ls))
	
	def getKey(self):
		"""
		2012.11.20 return a tuple that is dictionary-able.
		"""
		return (self.chromosome, self.start, self.stop)


class SegmentTreeNodeKey(CNVSegmentBinarySearchTreeKey):
	"""
	2010-1-28
		tree node key which counts both is_reciprocal_overlap()=True and "other" being wholly embedded in "self" as equal
		
		similar purpose as the function leftWithinRightAlsoEqualCmp().
		
		One strange thing (???) about who is "self", who is "other" in the __eq__() below.
			In the situation you have a new CNVSegmentBinarySearchTreeKey key (call it "segmentKey"),
			and want to test if this "segmentKey" is in the tree or not as in "if segmentKey in tree: ...".
			 
			If the tree is comprised of CNVSegmentBinarySearchTreeKey instances, self=segmentKey, other=nodes in the tree.
			If the tree is comprised of SegmentTreeNodeKey instances, self=nodes in the tree, other=segmentKey.
				That's why here used the condition that "other"'s overlap ratio (overlapFraction2) ==1..
			
		
	"""
	def __init__(self, chromosome=None, span_ls=None, min_reciprocal_overlap=0.6):
		CNVSegmentBinarySearchTreeKey.__init__(self, chromosome=chromosome, span_ls=span_ls, min_reciprocal_overlap=min_reciprocal_overlap)
	
	def __eq__(self, other):
		"""
		2010-1-28
		"""
		return rightWithinLeftAlsoEqual(self, other)


class CNVCompareBySmallOverlapRatio(object):
	"""
	2011-3-16
		renamed from CNVCompare to CNVCompareByMinRecipOverlapRatio
	2010-8-16
		a class to compare two CNVSegmentBinarySearchTreeKey according to any arbitrary min_reciprocal_overlap
	"""
	def __init__(self, min_reciprocal_overlap = 0.4):
		"""
		2010-8-16
		"""
		self.min_reciprocal_overlap = min_reciprocal_overlap
	
	def cmp(self, key1, key2):
		if self.eq(key1, key2, min_reciprocal_overlap=self.min_reciprocal_overlap):
			return 0
		elif self.lt(key1, key2, min_reciprocal_overlap=self.min_reciprocal_overlap):
			return -1
		else:
			return +1
	
	@classmethod
	def eq(cls, key1, key2, min_reciprocal_overlap=0.4):
		if key1.chromosome==key2.chromosome:
			if len(key1.span_ls)==1:
				if len(key2.span_ls)==1:
					return key1.span_ls[0]==key2.span_ls[0]
				elif len(key2.span_ls)>1:
					return key1.span_ls[0]>=key2.span_ls[0] and key1.span_ls[0]<=key2.span_ls[1]	# equal if self is within the "other" segment
				else:
					return None
			elif len(key1.span_ls)>1:
				if len(key2.span_ls)==1:	# self is a segment. other is a point position.
					return key1.span_ls[0]<=key2.span_ls[0] and key1.span_ls[1]>=key2.span_ls[0]	# if self includes the point position, yes it's equal
				elif len(key2.span_ls)>1:
					# need to calculate min_reciprocal_overlap
					return is_reciprocal_overlap(key1.span_ls, key2.span_ls, min_reciprocal_overlap=min_reciprocal_overlap)
				else:
					return None		#2010-8-16 shouldn't be here
		else:
			return False
		
		return True
	
	@classmethod
	def lt(cls, key1, key2, min_reciprocal_overlap=0.4):
		if key1.chromosome==key2.chromosome:
			if len(key1.span_ls)==1:
				if len(key2.span_ls)==1:
					return key1.span_ls[0]<key2.span_ls[0]
				elif len(key2.span_ls)>1:
					return key1.span_ls[0]<key2.span_ls[0]
				else:
					return None
			elif len(key1.span_ls)>1:
				if len(key2.span_ls)==1:
					return key1.span_ls[0]<key2.span_ls[0]
				elif len(key2.span_ls)>1:
					overlapData = get_overlap_ratio(key1.span_ls, key2.span_ls)
					overlapFraction1 = overlapData.overlapFraction1
					overlapFraction2 = overlapData.overlapFraction2
					recip_overlap = min(overlapFraction1, overlapFraction2)
					if recip_overlap>=min_reciprocal_overlap:	#should be equal,
						return False
					
					#if overlapFraction1>0 and overlapFraction1<=1 and  overlapFraction2>0 and overlapFraction2<=1:	# there's overlap between two segments 
					#	return False
					else:
						return key1.span_ls[0]<key2.span_ls[0]	# whether the start of key1 is ahead of the start of key2
				else:
					return None	#shouldn't be here
		else:
			return key1.chromosome<key2.chromosome
	
	@classmethod
	def gt(cls, key1, key2, min_reciprocal_overlap=0.4):
		"""
		2010-8-16
			not used
		"""
		pass

#2011-3-24
CNVCompare = CNVCompareBySmallOverlapRatio


class CNVCompareByOverlapLen(object):
	"""
	2011-3-16
		a class to compare two CNVSegmentBinarySearchTreeKey according to any arbitrary min_overlap_len
	"""
	def __init__(self, min_overlap_len = 1):
		"""
		2011-3-16
		"""
		self.min_overlap_len = min_overlap_len
	
	def cmp(self, key1, key2):
		if self.eq(key1, key2, min_overlap_len=self.min_overlap_len):
			return 0
		elif self.lt(key1, key2, min_overlap_len=self.min_overlap_len):
			return -1
		else:
			return +1
	
	@classmethod
	def eq(cls, key1, key2, min_overlap_len=1):
		if key1.chromosome==key2.chromosome:
			if len(key1.span_ls)==1:
				if len(key2.span_ls)==1:
					return key1.span_ls[0]==key2.span_ls[0]
				elif len(key2.span_ls)>1:
					return key1.span_ls[0]>=key2.span_ls[0] and key1.span_ls[0]<=key2.span_ls[1]	# equal if self is within the "other" segment
				else:
					return None
			elif len(key1.span_ls)>1:
				if len(key2.span_ls)==1:	# self is a segment. other is a point position.
					return key1.span_ls[0]<=key2.span_ls[0] and key1.span_ls[1]>=key2.span_ls[0]	# if self includes the point position, yes it's equal
				elif len(key2.span_ls)>1:
					# need to calculate min_reciprocal_overlap
					overlapData = get_overlap_ratio(key1.span_ls, key2.span_ls)
					overlapFraction1 = overlapData.overlapFraction1
					overlapFraction2 = overlapData.overlapFraction2
					overlap_length = overlapData.overlap_length
					if overlap_length>=min_overlap_len:	#should be equal,
						return True
					else:
						return False
				else:
					return None		#2010-8-16 shouldn't be here
		else:
			return False
		
		return True
	
	@classmethod
	def lt(cls, key1, key2, min_overlap_len=0.4):
		"""
		2011-3-16
		"""
		if key1.chromosome==key2.chromosome:
			if len(key1.span_ls)==1:
				if len(key2.span_ls)==1:
					return key1.span_ls[0]<key2.span_ls[0]
				elif len(key2.span_ls)>1:
					return key1.span_ls[0]<key2.span_ls[0]
				else:
					return None
			elif len(key1.span_ls)>1:
				if len(key2.span_ls)==1:
					return key1.span_ls[0]<key2.span_ls[0]
				elif len(key2.span_ls)>1:
					overlapData = get_overlap_ratio(key1.span_ls, key2.span_ls)
					overlapFraction1 = overlapData.overlapFraction1
					overlapFraction2 = overlapData.overlapFraction2
					overlap_length = overlapData.overlap_length
					if overlap_len>=min_overlap_len:	#should be equal,
						return False
					else:
						return key1.span_ls[0]<key2.span_ls[0]	# whether the start of key1 is ahead of the start of key2
				else:
					return None	#shouldn't be here
		else:
			return key1.chromosome<key2.chromosome
	
	@classmethod
	def gt(cls, key1, key2, min_reciprocal_overlap=0.4):
		"""
		2011-3-16
			not used
		"""
		pass

class CNVCompareByBigOverlapRatio(object):
	"""
	2011-3-16
		similar to CNVCompareBySmallOverlapRatio, but use max(overlapFraction1, overlapFraction2) to compare aginst min_overlap_len, \
			rather than min(overlapFraction1, overlapFraction2)
	2010-8-16
		a class to compare two CNVSegmentBinarySearchTreeKey according to any arbitrary min_reciprocal_overlap
	"""
	def __init__(self, min_reciprocal_overlap = 0.4):
		"""
		2010-8-16
		"""
		self.min_reciprocal_overlap = min_reciprocal_overlap
	
	def cmp(self, key1, key2):
		if self.eq(key1, key2, min_reciprocal_overlap=self.min_reciprocal_overlap):
			return 0
		elif self.lt(key1, key2, min_reciprocal_overlap=self.min_reciprocal_overlap):
			return -1
		else:
			return +1
	
	@classmethod
	def eq(cls, key1, key2, min_reciprocal_overlap=0.4):
		if key1.chromosome==key2.chromosome:
			if len(key1.span_ls)==1:
				if len(key2.span_ls)==1:
					return key1.span_ls[0]==key2.span_ls[0]
				elif len(key2.span_ls)>1:
					return key1.span_ls[0]>=key2.span_ls[0] and key1.span_ls[0]<=key2.span_ls[1]	# equal if self is within the "other" segment
				else:
					return None
			elif len(key1.span_ls)>1:
				if len(key2.span_ls)==1:	# self is a segment. other is a point position.
					return key1.span_ls[0]<=key2.span_ls[0] and key1.span_ls[1]>=key2.span_ls[0]	# if self includes the point position, yes it's equal
				elif len(key2.span_ls)>1:
					overlapData = get_overlap_ratio(key1.span_ls, key2.span_ls)
					overlapFraction1 = overlapData.overlapFraction1
					overlapFraction2 = overlapData.overlapFraction2
					overlap_length = overlapData.overlap_length
					recip_overlap = max(overlapFraction1, overlapFraction2)
					if recip_overlap>=min_reciprocal_overlap:	#should be equal,
						return True
					else:
						return False
				else:
					return None		#2010-8-16 shouldn't be here
		else:
			return False
		
		return True
	
	@classmethod
	def lt(cls, key1, key2, min_reciprocal_overlap=0.4):
		if key1.chromosome==key2.chromosome:
			if len(key1.span_ls)==1:
				if len(key2.span_ls)==1:
					return key1.span_ls[0]<key2.span_ls[0]
				elif len(key2.span_ls)>1:
					return key1.span_ls[0]<key2.span_ls[0]
				else:
					return None
			elif len(key1.span_ls)>1:
				if len(key2.span_ls)==1:
					return key1.span_ls[0]<key2.span_ls[0]
				elif len(key2.span_ls)>1:
					overlapData = get_overlap_ratio(key1.span_ls, key2.span_ls)
					overlapFraction1 = overlapData.overlapFraction1
					overlapFraction2 = overlapData.overlapFraction2
					recip_overlap = max(overlapFraction1, overlapFraction2)
					if recip_overlap>=min_reciprocal_overlap:	#should be equal,
						return False
					else:
						return key1.span_ls[0]<key2.span_ls[0]	# whether the start of key1 is ahead of the start of key2
				else:
					return None	#shouldn't be here
		else:
			return key1.chromosome<key2.chromosome
	
	@classmethod
	def gt(cls, key1, key2, min_reciprocal_overlap=0.4):
		"""
		2010-8-16
			not used
		"""
		pass

def leftWithinRightAlsoEqual(key1, key2):
	"""
	2010-1-28
		Besides CNVSegmentBinarySearchTreeKey.__eq__(), if key1 is embedded in key2, it's also regarded as equal.
		this function is solely used in leftWithinRightAlsoEqualCmp().
	"""
	equalResult = key1.__eq__(key2)
	if equalResult:
		return equalResult
	else:
		if len(key1.span_ls)==2 and len(key2.span_ls)==2:
			overlapData = get_overlap_ratio(key1.span_ls, key2.span_ls)
			overlapFraction1 = overlapData.overlapFraction1
			overlapFraction2 = overlapData.overlapFraction2
			if overlapFraction1==1.:	# 2010-1-28 added the overlapFraction2==1.
				return True
			else:
				return equalResult
		else:
			return equalResult

def rightWithinLeftAlsoEqual(key1, key2):
	"""
	2010-1-28
		Besides CNVSegmentBinarySearchTreeKey.__eq__(), if key2 is embedded in key1, it's also regarded as equal.
		this function is solely used in leftWithinRightAlsoEqualCmp().
	"""
	equalResult = key1.__eq__(key2)
	if equalResult:
		return equalResult
	else:
		if len(key1.span_ls)==2 and len(key2.span_ls)==2:
			overlapData = get_overlap_ratio(key1.span_ls, key2.span_ls)
			overlapFraction1 = overlapData.overlapFraction1
			overlapFraction2 = overlapData.overlapFraction2
			overlap_length = overlapData.overlap_length
					
			if overlapFraction2==1.:	# 2010-1-28 added the overlapFraction2==1.
				return True
			else:
				return equalResult
		else:
			return equalResult

def leftWithinRightAlsoEqualCmp(key1, key2):
	"""
	2010-1-28
		a cmp function for RBDict to compare CNVSegmentBinarySearchTreeKey or SegmentTreeNodeKey.
		
		If key1 is embedded in key2, it's also regarded as equal.
		
		similar purpose as the class SegmentTreeNodeKey()
		
		add this to RBDict as in "tree = RBDict(cmpfn=leftWithinRightAlsoEqualCmp)".
			key1 is the new entry or query. key2 is the one that is already in RBDict. 
	"""
	if leftWithinRightAlsoEqual(key1, key2):
		return 0
	elif key1 > key2:
		return 1
	else:  #key1 < key2
		return -1

def rightWithinLeftAlsoEqualCmp(key1, key2):
	"""
	2010-1-28
		a cmp function for RBDict to compare CNVSegmentBinarySearchTreeKey or SegmentTreeNodeKey.
		If key2 is embedded in key1, it's also regarded as equal.
		
		similar purpose as the class SegmentTreeNodeKey().
		
		add this to RBDict as in "tree = RBDict(cmpfn=rightWithinLeftAlsoEqualCmp)".
			key1 is the new entry or query. key2 is the one that is already in RBDict.
		 
	
	"""
	if rightWithinLeftAlsoEqual(key1, key2):
		return 0
	elif key1 > key2:
		return 1
	else:  #key1 < key2
		return -1

def getCNVDataFromFileInGWA(input_fname_ls, array_id, max_amp=-0.33, min_amp=-0.33, min_size=50, min_no_of_probes=None, \
						report=False, chr=None, start=None, stop=None, filter_gada_by_cutoff=True):
	"""
	2010-6-5
		data_obj.genome_wide_result_name = gwr_name
		add ecotyp_id into data_obj.comment
	2010-6-3
		add argument chr, start, stop to restrict fragments within those
	2009-10-31
		get deletion (below max_amp) or duplication (above min_amp) from files (output by RunGADA.py)
	"""
	sys.stderr.write("Getting CNV calls for array %s, min_size %s, min_no_of_probes %s from %s ..."%\
					(array_id, min_size, min_no_of_probes, repr(input_fname_ls)))
	
	gwr_name = "(a-id %s)"%(array_id)
	gwr = GenomeWideResult(name=gwr_name)
	gwr.data_obj_ls = []	#list and dictionary are crazy references.
	gwr.data_obj_id2index = {}
	genome_wide_result_id = id(gwr)
	
	amp_ls = []
	array_id2array = {}
	counter = 0
	real_counter = 0
	no_of_segments = 0
	input_handler = fileinput.input(input_fname_ls)
	header = input_handler.readline().strip().split('\t')
	col_name2index = getColName2IndexFromHeader(header)
	ecotype_id = None
	for line in input_handler:
		if line.find("array_id")!=-1:
			continue
		line = line.strip()
		row = line.split('\t')
		cnv_array_id = int(row[col_name2index['array_id']])
		cnv_ecotype_id = int(row[col_name2index.get('ecotype_id', col_name2index['array_id'])])
		counter += 1
		if cnv_array_id==array_id:
			no_of_segments += 1
			if ecotype_id is None:
				ecotype_id = cnv_ecotype_id
			start_probe = row[col_name2index['start_probe']].split('_')	# split chr_pos
			start_probe = map(int, start_probe)
			start_probe_id = row[col_name2index.get('start_probe_id', col_name2index['start_probe'])]
			
			stop_probe = row[col_name2index['end_probe']].split('_')
			stop_probe = map(int, stop_probe)
			if start_probe[0]!=stop_probe[0]:	#fragments on different chromosome. ignore.
				continue
			end_probe_id = row[col_name2index.get('end_probe_id', col_name2index['end_probe'])]
			
			no_of_probes = int(row[col_name2index['length']])
			if min_no_of_probes is not None and no_of_probes<min_no_of_probes:
				continue
			amplitude = float(row[col_name2index['amplitude']])
			segment_chromosome = start_probe[0]
			segment_start_pos = start_probe[1]-12
			segment_stop_pos = stop_probe[1]+12
			segment_length = abs(segment_stop_pos-segment_start_pos)
			if chr is not None and segment_chromosome!=chr:
				continue
			if start is not None and (segment_start_pos<start or segment_stop_pos<start):
				continue
			if stop is not None and (segment_start_pos>stop or segment_stop_pos>stop):
				continue
			if min_size is not None and segment_length<min_size:
				continue
			if amplitude<=max_amp:
				color = 'r'
			elif amplitude>=min_amp:
				color = 'b'
			else:
				color = 'g'
			if amplitude<=max_amp or amplitude>=min_amp or filter_gada_by_cutoff == False:
				real_counter += 1
				data_obj = DataObject(chromosome=segment_chromosome, position=segment_start_pos, stop_position=segment_stop_pos, \
									value=amplitude, color=color)
				data_obj.comment = 'ecotype-id %s, start probe-id %s, end probe-id %s, no of probes %s, size=%s'%\
							(ecotype_id, start_probe_id, end_probe_id, no_of_probes, segment_length)
				data_obj.genome_wide_result_name = gwr_name
				data_obj.genome_wide_result_id = genome_wide_result_id
				gwr.add_one_data_obj(data_obj)
				
		if report and counter%10000==0:
			sys.stderr.write('%s%s\t%s\t%s'%('\x08'*80, counter, no_of_segments, real_counter))
	sys.stderr.write("\n")
	
	if gwr.max_value<3:	# insertion at y=3
		gwr.max_value=3
	if gwr.min_value>-1:	# deletion at y = -1
		gwr.min_value = -1
	gwr.name = '%s '%ecotype_id +  gwr.name
	setattr(gwr, 'ecotype_id', ecotype_id)
	sys.stderr.write(" %s segments. Done.\n"%(len(gwr.data_obj_ls)))
	return gwr

def mergeTwoSegments(segment, next_segment):
	"""
	2010-7-28
		copied from CNV.mergeTwoCNVQCSegments() in variation/src/misc.py
		
	2010-7-19
		called by getCNVQCDataFromDB()
		
		segment = (row.chromosome, row.start, row.stop, \
				row.size_affected, row.no_of_probes_covered, row.copy_number, row.id, row.ecotype_id)
	"""
	new_segment_stop = max(segment[2], next_segment[2])
	
	# to get the no_of_probes_covered for the merged segment
	if segment[4] is not None and next_segment[4] is not None:
		overlapping_len = float(segment[2]-next_segment[1])	# it could be negative, which would end up increasing the no_of_probes_covered
		import numpy
		# estimate the number of probes in the overlapping region for each segment, take the average
		no_of_probes_in_the_overlapping_of_s1 = int(round(segment[4]*overlapping_len/segment[3]))
		no_of_probes_in_the_overlapping_of_s2 = int(round(next_segment[4]*overlapping_len/next_segment[3]))
		no_of_probes_in_the_overlapping = numpy.mean([no_of_probes_in_the_overlapping_of_s1, no_of_probes_in_the_overlapping_of_s2])
		no_of_probes_covered = segment[4] + next_segment[4] - no_of_probes_in_the_overlapping
	else:
		no_of_probes_covered = None
	
	new_segment = [segment[0], segment[1], new_segment_stop, new_segment_stop-segment[1]+1, no_of_probes_covered] + \
				segment[5:]
	return new_segment

def mergeOverlappingAmongSegments(segment_ls):
	"""
	2010-7-28
		called by turnSegmentGWRIntoRBDict()
	"""
	segment_ls.sort()	#make sure in chromosmal order
	new_segment_ls = []
	for i in xrange(len(segment_ls)):
		segment = segment_ls[i]
		if len(new_segment_ls)>0:
			previous_segment = new_segment_ls[-1]
		else:
			previous_segment = None
		
		if previous_segment is not None and previous_segment[0]==segment[0] and previous_segment[2]>=segment[1]-1:
			previous_segment = mergeTwoSegments(previous_segment, segment)
			new_segment_ls[-1] = previous_segment
		else:
			new_segment_ls.append(segment)
	sys.stderr.write("%s segments merged into %s non-overlapping segments.\n"%(len(segment_ls), len(new_segment_ls)))
	return new_segment_ls

def turnSegmentGWRIntoRBDict(gwr, extend_dist=20000, min_reciprocal_overlap=0.6, report=True):
	"""
	2010-7-28
		add functionality to merge overlapping segments.
		assuming data_obj_ls in gwr is in chromosomal order.
	2010-3-17
		extend_dist is used to enlarge the segments in each data_obj of gwr,
	"""
	sys.stderr.write("Turning a segment-gwr (start-stop style) into an RBDict ...\n")

	segment_ls = []
	for data_obj in gwr.data_obj_ls:
		start = max(data_obj.position-extend_dist, 0)
		stop = data_obj.stop_position+extend_dist
		segment = [data_obj.chromosome, start, stop, stop-start+1, None]	#the 4th bit (=None) is no_of_probes_covered.
		segment_ls.append(segment)
	
	segment_ls = mergeOverlappingAmongSegments(segment_ls)
	
	from RBTree import RBDict	# 2010-1-26 RBDict is more efficiency than binary_tree.
	rbDict = RBDict(cmpfn=leftWithinRightAlsoEqualCmp)

	for segment in segment_ls:
		chromosome, start, stop = segment[:3]
		segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=chromosome, span_ls=[start, stop], \
												min_reciprocal_overlap=min_reciprocal_overlap)
		rbDict[segmentKey] = segment
	if report:
		print "\tDepth of rbDict: %d" % (rbDict.depth())
		print "\tOptimum Depth: %f (%d) (%f%% depth efficiency)" % (rbDict.optimumdepth(), math.ceil(rbDict.optimumdepth()),
															  math.ceil(rbDict.optimumdepth()) / rbDict.depth())
	sys.stderr.write("%s objects converted.\n"%len(rbDict))
	return rbDict

def getProbeIntensityDataFromProbeXArray(input_fname, data_type=numpy.float32, \
						nonIntensityColumnLabelSet = set(['probes_id', 'chromosome', 'position']),\
						array_id=None):
	"""
	2010-7-28
		bug-fix:
			make sure _array_id and array_id are of the same type.
	2010-6-28
		add argument array_id to get intensity only for that array.
	2010-5-23
		add argument nonIntensityColumnLabelSet to find out which columns harbor intensity data.
		Non-intensity columns (probes_id, chromosome, position) are not bound to any particular column.
	2010-3-18
		copied from CNVNormalize.get_input()
	2009-10-28
		switch the default data_type to numpy.float32 to save memory on 64bit machines
	2009-9-28
		add argument data_type to specify data type of data_matrix.
		default is numpy.float (numpy.float could be float32, float64, float128 depending on the architecture).
			numpy.double is also fine.
	2009-5-18
		become classmethod
	"""
	sys.stderr.write("Getting Probe X Array intensity data from %s ..."%input_fname)
	import csv, subprocess
	reader = csv.reader(open(input_fname), delimiter=figureOutDelimiter(input_fname))
	commandline = 'wc -l %s'%input_fname
	command_handler = subprocess.Popen(commandline, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	stdout_content, stderr_content = command_handler.communicate()
	if stderr_content:
		sys.stderr.write('stderr of %s: %s \n'%(commandline, stderr_content))
	no_of_rows = int(stdout_content.split()[0])-1
	
	header = reader.next()
	intensity_col_index_ls = []
	type_of_given_array_id = type(array_id)	#2010-7-28 
	for i in xrange(len(header)):
		label = header[i]
		if label not in nonIntensityColumnLabelSet:	# it's array_id if it's not in nonIntensityColumnLabelSet.
			_array_id = type_of_given_array_id(label)	#2010-7-28 cast to whatever type array_id is of.
			if array_id is not None and _array_id!=array_id:	#2010-6-28	skip this row
				continue
			intensity_col_index_ls.append(i)
	no_of_cols = len(intensity_col_index_ls)
	col_name2index = getColName2IndexFromHeader(header)
	data_matrix = numpy.zeros([no_of_rows, no_of_cols], data_type)
	probe_id_ls = []
	chr_pos_ls = []
	i=0
	for row in reader:
		probe_id = row[col_name2index["probes_id"]]
		probe_id_ls.append(probe_id)
		chromosome = row[col_name2index["chromosome"]]
		pos = row[col_name2index["position"]]
		chr_pos_ls.append((chromosome, pos))
		for j in intensity_col_index_ls:
			try:
				data_matrix[i][j-1] = float(row[j])
			except:
				data_matrix[i][j-1] = numpy.nan
		i += 1
	del reader
	sys.stderr.write("Done.\n")
	return data_matrix, probe_id_ls, chr_pos_ls, header

def getProbeIntensityDataFromArrayXProbe(input_fname, data_type=numpy.float32, \
						nonIntensityColumnLabelSet = None, array_id=None, no_of_arrays=None):
	"""
	2010-7-28
		add argument no_of_arrays to skip the "wc -l" command
		bug-fix:
			make sure _array_id and array_id are of the same type.
	2010-6-28
		add argument array_id to get intensity only for that array.
	2010-5-25
		similar to getProbeIntensityDataFromProbeXArray, but the input is in Probe X Array. 
	"""
	sys.stderr.write("Getting Array X Probe intensity data from %s for array %s ...\n"%(input_fname, array_id))
	import csv, subprocess
	reader = csv.reader(open(input_fname), delimiter=figureOutDelimiter(input_fname))
	if array_id is not None:	#2010-7-28
		no_of_arrays = 1
	elif no_of_arrays is not None:
		no_of_arrays = no_of_arrays
	else:
		commandline = 'wc -l %s'%input_fname
		command_handler = subprocess.Popen(commandline, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
		stdout_content, stderr_content = command_handler.communicate()
		if stderr_content:
			sys.stderr.write('stderr of %s: %s \n'%(commandline, stderr_content))
		no_of_arrays = int(stdout_content.split()[0])-3
	
	header = reader.next()
	probe_id_ls = header[2:]
	chr_ls = reader.next()[2:]
	pos_ls = reader.next()[2:]
	chr_pos_ls = zip(chr_ls, pos_ls)
	
	no_of_probes = len(header)-2
	data_matrix = numpy.zeros([no_of_probes, no_of_arrays], data_type)
	array_id_ls = []
	i=0
	type_of_given_array_id = type(array_id)#2010-7-28
	line_number = 0
	for row in reader:
		line_number += 1
		sys.stderr.write('%s%s'%("\x08"*80, line_number))
		
		if array_id is not None and len(array_id_ls)>0:	#2010-8-6 already got the array, jump out.
			break
		
		_array_id = type_of_given_array_id(row[0])	#2010-7-28 cast to whatever type array_id is of.
		if array_id is not None and _array_id!=array_id:	#2010-6-28	skip this row
			continue
		ecotype_id = row[1]
		
		try:
			data_intensity_ls = map(float, row[2:])
			data_matrix[:,i] = data_intensity_ls
			array_id_ls.append(_array_id)
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
			sys.stderr.write("array %s, ecotype %s, row %s"%(_array_id, ecotype_id, i))
			
		i += 1
	header = ['probes_id'] + array_id_ls + ['chromosome', 'position']
	del reader
	sys.stderr.write("  Done.\n")
	return data_matrix, probe_id_ls, chr_pos_ls, header

def getProbeIntensityData(input_fname, data_type=numpy.float32, \
						nonIntensityColumnLabelSet = set(['probes_id', 'chromosome', 'position']),\
						array_id = None):
	"""
	2010-6-28
		add argument array_id to get intensity only for that array.
	2010-5-30
		use header[0] to guess the orientation, and call getProbeIntensityDataFromArrayXProbe() or getProbeIntensityDataFromProbeXArray()
	2010-5-23
		add argument nonIntensityColumnLabelSet to find out which columns harbor intensity data.
		Non-intensity columns (probes_id, chromosome, position) are not bound to any particular column.
	2010-3-18
		copied from CNVNormalize.get_input()
	2009-10-28
		switch the default data_type to numpy.float32 to save memory on 64bit machines
	2009-9-28
		add argument data_type to specify data type of data_matrix.
		default is numpy.float (numpy.float could be float32, float64, float128 depending on the architecture).
			numpy.double is also fine.
	2009-5-18
		become classmethod
	"""
	sys.stderr.write("Testing which type of intensity data ... ")
	import csv, subprocess
	reader = csv.reader(open(input_fname), delimiter=figureOutDelimiter(input_fname))
	header = reader.next()
	del reader
	if header[0]=='':
		# nothing there. it's Array X Probe format
		return getProbeIntensityDataFromArrayXProbe(input_fname, data_type=data_type, array_id=array_id)
	else:
		return getProbeIntensityDataFromProbeXArray(input_fname, data_type=data_type, \
						nonIntensityColumnLabelSet = nonIntensityColumnLabelSet, array_id=array_id)

class TilingProbeIntensityData(object):
	"""
	2010-7-28
		a data structure which makes accessing the intensity data of any array from the big file easy and memory-efficient.
		
	"""
	def __init__(self, input_fname = None, min_reciprocal_overlap=None):
		self.input_fname = input_fname
		self.min_reciprocal_overlap = min_reciprocal_overlap
		self.array_id2index = {}
		self.intensity_matrix = []
		self.array_id2gwr = {}
		self.probe_id_ls = None
		self.chr_pos_ls = None
	
	def getIntensityForOneArrayInGWRGivenRBDict(self, array_id, rbDict=None, additionalTitle=None):
		"""
		2010-7-28
		"""
		if array_id not in self.array_id2index:	#first check if it's fetched already.		
			data_matrix, probe_id_ls, chr_pos_ls, header = getProbeIntensityData(self.input_fname, array_id=array_id)
			
			array_id_ls = header[1:-2]
			array_id_ls = map(int, array_id_ls)
			if len(array_id_ls)==0:	#nothing found for this array in the file
				sys.stderr.write("No data for array %s in %s.\n"%(array_id, self.input_fname))
				return None
			
			if self.probe_id_ls is None:
				self.probe_id_ls = probe_id_ls
			if self.chr_pos_ls is None:
				self.chr_pos_ls = chr_pos_ls
			
			self.array_id2index[array_id] = len(self.array_id2index)
			self.intensity_matrix.append(data_matrix[:,0])		# data_matrix is in probe X array format.
			
		from SNP import GenomeWideResult, DataObject
		gwr_name = 'array %s probe intensity'%array_id
		if additionalTitle:
			gwr_name += ' for %s'%(additionalTitle)
		
		gwr = GenomeWideResult(name=gwr_name)
		# 2010-3-18 custom
		gwr.array_id = array_id
		#gwr.ecotype_id = array.maternal_ecotype_id
		#gwr.nativename = ecotype_nativename
		
		genome_wide_result_id = id(gwr)
		
		array_index = self.array_id2index.get(array_id)
		
		no_of_rows = len(self.probe_id_ls)
		for i in range(no_of_rows):
			chr_pos = self.chr_pos_ls[i]
			chr, pos = map(int, chr_pos)
			cnvSegmentKey = CNVSegmentBinarySearchTreeKey(chromosome=chr, span_ls=[pos],\
														min_reciprocal_overlap=self.min_reciprocal_overlap)
			if rbDict is None or cnvSegmentKey in rbDict:
				probeIntensity = self.intensity_matrix[array_index][i]
				data_obj = DataObject(chromosome=chr, position=pos, value=probeIntensity)
				data_obj.comment = ''
				data_obj.genome_wide_result_name = gwr_name
				data_obj.genome_wide_result_id = genome_wide_result_id
				gwr.add_one_data_obj(data_obj)
		sys.stderr.write(" %s probes. Done.\n"%(len(gwr.data_obj_ls)))
		return gwr

def fetchIntensityInGWAWithinRBDictGivenArrayIDFromTilingIntensity(tilingIntensityData, array_id, rbDict, gwr_name=None,\
																min_reciprocal_overlap=0.6):
	"""
	2010-3-18
		tilingIntensityData is of type SNPData.
		
	"""
	sys.stderr.write("Getting intensity data within the chosen segments for array %s ..."%array_id)
	col_index = tilingIntensityData.col_id2col_index.get(array_id)
	if col_index is None:
		sys.stderr.write("Error: No tiling intensity for array %s.\n"%array_id)
		return None
	
	from SNP import GenomeWideResult, DataObject
	
	gwr = GenomeWideResult(name=gwr_name)
	# 2010-3-18 custom
	gwr.array_id = array_id
	#gwr.ecotype_id = array.maternal_ecotype_id
	#gwr.nativename = ecotype_nativename
	
	genome_wide_result_id = id(gwr)
		
	no_of_rows = len(tilingIntensityData.row_id_ls)
	for i in range(no_of_rows):
		chr_pos = tilingIntensityData.row_id_ls[i]
		chr, pos = map(int, chr_pos)
		cnvSegmentKey = CNVSegmentBinarySearchTreeKey(chromosome=chr, span_ls=[pos],\
													min_reciprocal_overlap=min_reciprocal_overlap)
		if cnvSegmentKey in rbDict:
			probeIntensity = tilingIntensityData.data_matrix[i][col_index]
			data_obj = DataObject(chromosome=chr, position=pos, value=probeIntensity)
			data_obj.comment = ''
			data_obj.genome_wide_result_name = gwr_name
			data_obj.genome_wide_result_id = genome_wide_result_id
			gwr.add_one_data_obj(data_obj)
	sys.stderr.write(" %s probes. Done.\n"%(len(gwr.data_obj_ls)))
	return gwr


class ArrayXProbeFileWrapper(object):
	"""
	2010-6-3
		a wrapper around CNV input file, which is in array X probe format.
		It will read 3 lines of headers first and then let users to access file through self.reader.
	"""
	def __init__(self, input_fname):
		"""
		"""
		self.input_fname = input_fname
		import csv
		from pymodule import figureOutDelimiter
		self.reader = csv.reader(open(self.input_fname), delimiter=figureOutDelimiter(self.input_fname))
		self.probe_id_ls, self.chr_pos_ls = self.getHeader(self.reader)
	
	@classmethod
	def getHeader(cls, reader):
		sys.stderr.write("Getting probe id, chr, pos info ...")
		probe_id_ls = reader.next()[2:]
		probe_id_ls = map(int, probe_id_ls)
		chr_ls = reader.next()[2:]
		pos_ls = reader.next()[2:]
		chr_pos_ls = zip(chr_ls, pos_ls)
		sys.stderr.write(".\n")
		return probe_id_ls, chr_pos_ls

def get_chr2start_stop_index(chr_pos_ls):
	"""
	2010-6-5
		get the boundary of all chromosomes based on chr_pos_ls
	"""
	sys.stderr.write("Getting chr2start_stop_index ...")
	no_of_probes = len(chr_pos_ls)
	chr2start_stop_index = {}
	old_chr = None
	for i in xrange(no_of_probes):
		chr = chr_pos_ls[i][0]
		if chr not in chr2start_stop_index:
			chr2start_stop_index[chr] = [i]
		if old_chr is not None and chr!=old_chr:
			chr2start_stop_index[old_chr].append(i-1)
		
		old_chr = chr
	chr2start_stop_index[old_chr].append(i)	#append the the stop index to the last chromosome
	sys.stderr.write("Done.\n")
	return chr2start_stop_index

def readQuanLongPECoverageIntoGWR(input_fname, additionalTitle=None, windowSize=100, chr=None, start=None, stop=None):
	"""
	2010-11-23
		This function reads 100bp coverage data from files given by Quan and turns it into instance of GenomeWideResult.
			example file: "variation/data/CNV/QuanLongPE/coverage4Yu/algustrum.8230.Chr1.coverage"
			format is two-column, tab-delimited. First column is position/100. 2nd column is coverage.:
				1       0.0
				2       0.0
				3       2.84
		
		the coverage data is split into different chromosomes. so argument chr is not used here.
	"""
	from CNV import CNVSegmentBinarySearchTreeKey
	from SNP import GenomeWideResult, DataObject
	import os, sys, re, csv
	sys.stderr.write("Reading Quan's coverage data from %s ... "%(input_fname))
	chrPattern = re.compile(r'.*Chr(\d+).coverage')
	chrPatternSearchResult = chrPattern.search(input_fname)
	if chrPatternSearchResult:
		chr = int(chrPatternSearchResult.group(1))
	else:
		sys.stderr.write("Can't parse chromosome out of %s. aborted.\n"%input_fname)
	reader = csv.reader(open(input_fname), delimiter='\t')
	gwr_name = 'log10(coverage) %s'%os.path.basename(input_fname)
	if additionalTitle:
		gwr_name += " "+ additionalTitle
	
	coverage_gwr = GenomeWideResult(name=gwr_name)
	
	genome_wide_result_id = id(coverage_gwr)
	
	for row in reader:
		position, coverage = row[:2]
		position = int(position)*windowSize-windowSize/2	#coverage is caculated in 100-base window.
		if start is not None and position<start:
			continue
		if stop is not None and position>stop:
			continue
		coverage = float(coverage)
		if coverage<=0.0:
			coverage = -3
		else:
			coverage = math.log10(coverage)
		data_obj = DataObject(chromosome=chr, position=position, value=coverage)
		data_obj.comment = ''
		data_obj.genome_wide_result_name = gwr_name
		data_obj.genome_wide_result_id = genome_wide_result_id
		coverage_gwr.add_one_data_obj(data_obj)
	sys.stderr.write("Done.\n")
	return coverage_gwr


def findCorrespondenceBetweenTwoCNVRBDict(rbDict1=None, rbDict2=None):
	"""
	2013.1.27
	"""
	dc1Length = len(rbDict1)
	dc2Length = len(rbDict2)
	sys.stderr.write("Finding correspondence between two CNV RB dictionaries, %s and %s nodes in input1, input2 respectively ... "%\
					(dc1Length, dc2Length))
	nodePairList = []	#each cell is a tuple of (node1, node2, overlap)
	compareIns = CNVCompare(min_reciprocal_overlap=0.0000001)	#to detect any overlap
	setOfMatchedNodesInRBDict2 = set()
	for input1Node in rbDict1:
		targetNodeLs = []
		rbDict2.findNodes(input1Node.key, node_ls=targetNodeLs, compareIns=compareIns)
		if targetNodeLs:
			for input2Node in targetNodeLs:
				overlapData = get_overlap_ratio(input1Node.key.span_ls, input2Node.key.span_ls)
				overlapFraction = overlapData.overlapFraction
				nodePairList.append((input1Node, input2Node, overlapFraction))
				setOfMatchedNodesInRBDict2.add(input2Node)
		else:
			nodePairList.append((input1Node, None, None))
	sys.stderr.write(" %s input2 nodes matched. "%(len(setOfMatchedNodesInRBDict2)))
	
	no_of_unmatched_input2_nodes_matched_to_input1 = 0
	for input2Node in rbDict2:
		if input2Node not in setOfMatchedNodesInRBDict2:
			targetNodeLs = []
			rbDict1.findNodes(input2Node.key, node_ls=targetNodeLs, compareIns=compareIns)
			if targetNodeLs:	#this should not happen
				sys.stderr.write("Warning: after first scan node from rbDict2, key=%s, value=%s, still has hits from rbDict1.\n"%\
								(str(input2Node.key), str(input2Node.value[0])))
				no_of_unmatched_input2_nodes_matched_to_input1 += 1
				
				for input1Node in targetNodeLs:
					overlapData = get_overlap_ratio(input1Node.key.span_ls, input2Node.key.span_ls)
					overlapFraction = overlapData.overlapFraction
					nodePairList.append((input1Node, input2Node, overlapFraction))
			else:
				nodePairList.append((None, input2Node, None))
	sys.stderr.write(" %s un-matched input2 nodes have matches from input1. found %s corresponding pairs.\n"%\
					(no_of_unmatched_input2_nodes_matched_to_input1, len(nodePairList)))
	return nodePairList


# test program if this file is run
if __name__ == "__main__":
	import os, sys, math
	
	if len(sys.argv)>=2 and sys.argv[1]=='-b':	# 2010-6-17 enter debug mode "CNV.py -b"
		import pdb
		pdb.set_trace()
		debug = True
	else:
		debug = False
	
	cnv_ls = [[1, (2323,2600)], [2,(50000,)], [3,(43214,78788)], [5,(150,500)], [5,(500,950)], [5, (43241, 43242)]]
	no_of_cnvs = len(cnv_ls)
	min_reciprocal_overlap = 0.6
	
	#from BinarySearchTree import binary_tree
	#tree = binary_tree()
	from RBTree import RBDict	#2010-1-26 binary_tree and RBDict are swappable. but RBDict is more efficient (balanced).
	tree = RBDict(cmpfn=leftWithinRightAlsoEqualCmp)	# 2010-1-28 use the custom cmpfn if you want the case that left within right is regarded as equal as well.  
	tree = RBDict(cmpfn=rightWithinLeftAlsoEqualCmp)
	
	for cnv in cnv_ls:
		segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=cnv[0], span_ls=cnv[1], min_reciprocal_overlap=min_reciprocal_overlap)
		tree[segmentKey] = cnv
	
	print "Binary Tree Test\n"
	print "Node Count: %d" % len(tree)
	print "Depth: %d" % tree.depth()
	print "Optimum Depth: %f (%d) (%f%% depth efficiency)" % (tree.optimumdepth(), math.ceil(tree.optimumdepth()),
															  math.ceil(tree.optimumdepth()) / tree.depth())
	
	print "Efficiency: %f%% (total possible used: %d, total wasted: %d)" % (tree.efficiency() * 100,
																			  len(tree) / tree.efficiency(),
																			  (len(tree) / tree.efficiency()) - len(tree))
	"""
	print "Min: %s" % repr(tree.min())
	print "Max: %s" % repr(tree.max())
	
	print "List of Layers:\n\t" + repr(tree.listlayers()) + "\n"
	print "\"Recursive\" List:\n\t" + repr(tree.listrecursive()) + "\n"
	print "List of Keys:\n\t" + repr(tree.listkeys()) + "\n"
	print "List of Data:\n\t" + repr(tree.listdata()) + "\n"
	print "List of Nodes:\n\t" + repr(tree.listnodes()) + "\n"
	print "Dictionary:\n\t" + repr(tree.dict()) + "\n"
	print "Formatted Tree:\n" + tree.formattree() + "\n"
	print "Formatted Tree (Root in Middle):\n" + tree.formattreemiddle() + "\n"
	"""
	test_cnv_ls = [	[2,(50000,)], [3,(43214,43219)], [3,(43214,78788)], [5, (43242,)], [5,(144,566)], [5, (144, 1000)], [5, (50000, 70000)]]
	for test_cnv in test_cnv_ls:
		segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=test_cnv[0], span_ls=test_cnv[1], \
												min_reciprocal_overlap=min_reciprocal_overlap)
		print "segmentKey", segmentKey
		if segmentKey in tree:
			targetSegment = tree.get(segmentKey)
			#if targetSegment:
			print "\tIn tree with target", targetSegment
			node_ls = []
			tree.findNodes(segmentKey, node_ls)
			print "\tfindNodes()", [node.value for node in node_ls]
		else:
			print "\tNot in tree"
	
	
	"""
	for i in range(no_of_cnvs):
		for j in range(i+1, no_of_cnvs):
			cnv1 = cnv_ls[i]
			cnv2 = cnv_ls[j]
			segmentKey1 = CNVSegmentBinarySearchTreeKey(chromosome=cnv1[0], span_ls=cnv1[1], \
					min_reciprocal_overlap=min_reciprocal_overlap)
			segmentKey2 = CNVSegmentBinarySearchTreeKey(chromosome=cnv2[0], span_ls=cnv2[1], \
					min_reciprocal_overlap=min_reciprocal_overlap)
			print segmentKey1, "vs", segmentKey2
			print ">", segmentKey1>segmentKey2
			print ">=", segmentKey1>=segmentKey2
			print "<", segmentKey1<segmentKey2
			print "<=", segmentKey1<=segmentKey2
			print "==", segmentKey1==segmentKey2
	"""