#!/usr/bin/env python
"""
2012.1.3
	module for graph algorithms/data structures
"""

import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from pymodule.utils import PassingData
import copy
import numpy
from networkx import Graph, DiGraph
import networkx as nx

class GraphWrapper(Graph):
	def __init__(self, data=None, **keywords):
		"""
		2013.1.3
			wrapper around networkx.Graph
		"""
		Graph.__init__(self, data=None, **keywords)
		
	def recursiveTrimOutOfSetLeafNodes(self, nodeIdSet=None, recursiveDepth=0, total_no_of_nodes_removed=0):
		"""
		2013.1.2
			recursively remove all nodes with no out-degree.
			for directed graph. nodes with out_degree=0 are leaf nodes.
		"""
		if recursiveDepth==0:
			sys.stderr.write("Trimming leaf nodes that are not in given set (%s elements) ... \n"%\
						(len(nodeIdSet)))
		else:
			sys.stderr.write("%sdepth=%s"%('\x08'*80, recursiveDepth))
		no_of_nodes_removed = 0
		for n in self.nodes():
			if n not in nodeIdSet and self.out_degree(n)==0:
				self.remove_node(n)
				no_of_nodes_removed += 1
		total_no_of_nodes_removed += no_of_nodes_removed
		if no_of_nodes_removed>0:
			self.recursiveTrimOutOfSetLeafNodes(nodeIdSet=nodeIdSet, recursiveDepth=recursiveDepth+1,\
											total_no_of_nodes_removed=total_no_of_nodes_removed)
		else:
			sys.stderr.write("\n %s nodes removed."%total_no_of_nodes_removed)
		if recursiveDepth==0:
			sys.stderr.write("\n")
	
	def recursiveRemoveUniDegreeNodes(self):
		"""
		2013.1.2
			recursively remove all nodes whose degree (in + out) <= 1
		"""
		self.recursiveRemoveNodesByDegree(maxDegree=1, recursiveDepth=0)
	
	def recursiveRemoveNodesByDegree(self, maxDegree=0, recursiveDepth=0, total_no_of_nodes_removed=0, report=False):
		"""
		2013.1.2
			recursively remove all nodes whose degree (in + out) is below or equal to maxDegree 
		"""
		if report:
			if recursiveDepth==0:
				sys.stderr.write("Recursively remove nodes whose degree is <=%s ...\n"%maxDegree)
			else:
				sys.stderr.write("%sdepth=%s"%('\x08'*80, recursiveDepth))
		no_of_nodes_removed = 0
		for n in self.nodes():
			if self.degree(n)<=maxDegree:
				self.remove_node(n)
				no_of_nodes_removed += 1
		total_no_of_nodes_removed += no_of_nodes_removed
		if no_of_nodes_removed>0:
			self.recursiveRemoveNodesByDegree(maxDegree=maxDegree, recursiveDepth=recursiveDepth+1,\
											total_no_of_nodes_removed=total_no_of_nodes_removed)
		elif report:	#final round
				sys.stderr.write("\n %s nodes removed."%total_no_of_nodes_removed)
		if recursiveDepth==0 and report:
			sys.stderr.write("\n")
	
	def findAllLeafNodes(self):
		"""
		2013.1.3
		"""
		sys.stderr.write("Finding all leaf nodes ...")
		nodeList = []
		for n in self.nodes():
			if self.out_degree(n)==0:
				nodeList.append(n)
		sys.stderr.write(" %s leaves \n"%(len(nodeList)))
		return nodeList
	
	def calculateNodeHierarchyLevel(self):
		"""
		2013.1.3
			the hierarchy level for each node is the shortest path to all leaf nodes.
		"""
		
		leafNodes = set(self.findAllLeafNodes())
		sys.stderr.write("Calculating node hierarchy level ...")
		node2HierarchyLevel = {}
		#leaf nodes' hierarchy level=0
		
		for source in self.nodes():
			if source in leafNodes:
				node2HierarchyLevel[source] = 0
			else:
				level = None
				for target in leafNodes:
					try:
						l = nx.astar_path_length(self, source, target)
						if level is None:
							level = l
						elif l<level:
							level = l
					except:
						#no path between source and target
						pass
				node2HierarchyLevel[source] = level
		sys.stderr.write("%s nodes with hierarchy level.\n"%(len(node2HierarchyLevel)))
		return node2HierarchyLevel
	
	def findAllFounders(self):
		"""
		2013.3.5
			opposite of findAllLeafNodes()
		"""
		sys.stderr.write("Finding all founders ...")
		nodeList = []
		for n in self.nodes():
			if self.in_degree(n)==0:
				nodeList.append(n)
		sys.stderr.write(" %s founders \n"%(len(nodeList)))
		return nodeList
	
	def orderMembersByDistanceToFounders(self, founderSet=None):
		"""
		2013.3.5
			similar to calculateNodeHierarchyLevel
				but calculates shortest path to all founders (not leaf nodes)
					and the distance of the longest (not shortest) shortest path is the level.
				level of founder = 0
				level of founders' kids = 1
			This is good to simulate genotype of pedigree from founders to everyone.
		"""
		if founderSet is None:
			founderSet = set(self.findAllFounders())
		
		sys.stderr.write("Ordering members of the pedigree based on distance to founders ...")
		founderDistance2NodeList = {}
		#leaf nodes' hierarchy level=0
		counter = 0
		for source in self.nodes():
			if source in founderSet:
				founderDistance2NodeList[source] = 0
			else:
				level = None
				for target in founderSet:
					try:
						l = nx.astar_path_length(self, source, target)
						if level is None:
							level = l
						elif l>level:
							level = l
					except:
						#no path between source and target
						pass
				founderDistance2NodeList[source] = level
			counter += 1
		sys.stderr.write("%s different hierarchy level.\n"%(len(founderDistance2NodeList)))
		return founderDistance2NodeList
	
	
class DiGraphWrapper(DiGraph, GraphWrapper):
	def __init__(self, data=None, **keywords):
		"""
		2013.3.5, put GraphWrapper, behind DiGraph in inheritance
		2013.1.3
		"""
		DiGraph.__init__(self, data=None, **keywords)
		

if __name__ == '__main__':
	import pdb
	pdb.set_trace()
	
