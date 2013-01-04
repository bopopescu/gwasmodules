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
		
	def trimOutOfSetLeafNodes(self, nodeIdSet=None):
		"""
		2013.1.2
			recursively remove all nodes with no out-degree.
			for directed graph. nodes with out_degree=0 are leaf nodes.
		"""
		no_of_nodes_removed = 0
		for n in self.nodes():
			if n not in nodeIdSet and self.out_degree(n)==0:
				self.remove_node(n)
				no_of_nodes_removed += 1
		if no_of_nodes_removed>0:
			self.trimOutOfSetLeafNodes(nodeIdSet=nodeIdSet)
		
	def removeUniDegreeNodes(self):
		"""
		2013.1.2
			recursively remove all nodes whose degree (in + out) <= 1
		"""
		self.removeNodesByDegree(maxDegree=1)
	
	def removeNodesByDegree(self, maxDegree=0):
		"""
		2013.1.2
			recursively remove all nodes whose degree (in + out) is below or equal to maxDegree 
		"""
		no_of_nodes_removed = 0
		for n in self.nodes():
			if self.degree(n)<=maxDegree:
				self.remove_node(n)
				no_of_nodes_removed += 1
		if no_of_nodes_removed>0:
			self.removeNodesByDegree(maxDegree=maxDegree)
	
	def findAllLeafNodes(self):
		"""
		2013.1.3
		"""
		nodeList = []
		for n in self.nodes():
			if self.out_degree(n)==0:
				nodeList.append(n)
		return nodeList
	
	def calculateNodeHierarchyLevel(self):
		"""
		2013.1.3
			the hierarchy level for each node is the shortest path to all leaf nodes.
		"""
		leafNodes = set(self.findAllLeafNodes())
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
		return node2HierarchyLevel
	
class DiGraphWrapper(GraphWrapper, DiGraph):
	def __init__(self, data=None, **keywords):
		"""
		2013.1.3
		"""
		DiGraph.__init__(self, data=None, **keywords)
		

if __name__ == '__main__':
	import pdb
	pdb.set_trace()
	
