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
	
	
	def _recursiveCalculateNodeDistanceToLeaf(self, leafNodes=None, node2distanceToLeaf=None):
		"""
		2013.10.16
		"""
		newLeafNodeSet = set()
		for leafNode in leafNodes:
			if leafNode not in node2distanceToLeaf:
				node2distanceToLeaf[leafNode] = 0
			for parent in self.predecessors(leafNode):
				node2distanceToLeaf[parent] = node2distanceToLeaf[leafNode] + 1
				newLeafNodeSet.add(parent)
		if newLeafNodeSet:
			self._recursiveCalculateNodeDistanceToLeaf(newLeafNodeSet, node2distanceToLeaf)
		return node2distanceToLeaf
	
	def calculateNodeHierarchyLevel(self):
		"""
		2013.1.3
			the hierarchy level for each node is the shortest path to all leaf nodes.
		"""
		
		leafNodes = set(self.findAllLeafNodes())
		sys.stderr.write("Calculating node hierarchy level ...")
		node2distanceToLeaf = {}
		#leaf nodes' hierarchy level=0
		self._recursiveCalculateNodeDistanceToLeaf(leafNodes, node2distanceToLeaf)
		"""
		for source in self.nodes():
			if source in leafNodes:
				node2distanceToLeaf[source] = 0
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
				node2distanceToLeaf[source] = level
		"""
		sys.stderr.write("%s nodes with hierarchy level.\n"%(len(node2distanceToLeaf)))
		return node2distanceToLeaf
	
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
		nodeID2distanceToFounder = {}
		#leaf nodes' hierarchy level=0
		counter = 0
		for source in self.nodes():
			if source in founderSet:
				nodeID2distanceToFounder[source] = 0
			else:
				level = None
				#find the longest path to the founder
				for founder in founderSet:
					try:
						l = nx.astar_path_length(self, founder, source)
						if level is None:
							level = l
						elif l>level:
							level = l
					except:
						#no path between source and founder
						pass
				nodeID2distanceToFounder[source] = level
			counter += 1
		sys.stderr.write("%s different hierarchy level.\n"%(len(nodeID2distanceToFounder)))
		return nodeID2distanceToFounder
	
	
class DiGraphWrapper(DiGraph, GraphWrapper):
	def __init__(self, data=None, **keywords):
		"""
		2013.3.5, put GraphWrapper, behind DiGraph in inheritance
		2013.1.3
		"""
		DiGraph.__init__(self, data=None, **keywords)
		self._undirectedG = None
	
	@property
	def undirectedGraph(self):
		"""
		2013.10.01
			assuming graph does not change anymore.
			
		"""
		if self._undirectedG is None:
			self._undirectedG = self.to_undirected()
		return self._undirectedG
	
	def getShortestPathInUndirectedVersion(self, nodeID1=None, nodeID2=None):
		"""
		2013.10.01
			get undirected version of graph, get shortest path between two nodes,
			then add edge direction to each edge
			
		"""
		shortestPath = nx.astar_path(self.undirectedGraph, source=nodeID1, target=nodeID2)
		pathWithDirection = []
		for i in xrange(len(shortestPath)-1):
			if self.has_edge(shortestPath[i], shortestPath[i+1]):
				edgeDirection = +1
			elif self.has_edge(shortestPath[i+1], shortestPath[i]):
				edgeDirection = -1
			else:
				sys.stderr.write("ERROR: Edge %s -> %s exists in undirected graph, but it (or its reverse) does not exist in directed version.\n"%
								(shortestPath[i], shortestPath[i+1]))
				raise
			pathWithDirection.append((shortestPath[i], shortestPath[i+1], edgeDirection))
		return pathWithDirection

if __name__ == '__main__':
	import pdb
	pdb.set_trace()
	
