#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s

Description:
	2013.05.03 a child class of MatrixFile. used to describe pedigree file in plink format (delimiter could be ' ' or tab), which looks like:
		The header line does not exist. Six columns are: FamilyID IndividualID PaternalID MaternalID Sex(1=male; 2=female; other=unknown) Phenotype.
		
		F0_1 990_699_1984014_GA_vs_524copy1 0 0 2 1
		F0_1 1513_641_1986014_GA_vs_524copy1 0 0 1 1
		F0_1 984_693_1996027_GA_vs_524copy1 1513_641_1986014_GA_vs_524copy1 990_699_1984014_GA_vs_524copy1 1 1
		F1_1 1582_1672_1993040_GA_vs_524copy1 0 0 1 1
		F1_1 1917_2966_1992045_GA_vs_524copy1 0 0 2 1
		F1_1 1931_2980_2000040_GA_vs_524copy1 1582_1672_1993040_GA_vs_524copy1 1917_2966_1992045_GA_vs_524copy1 1 1
	
	Example:
	
		plinkPedigreeFile = PlinkPedigreeFile(inputFname=self.pedigreeFname)
		familyStructureData = plinkPedigreeFile.getFamilyStructurePlinkWay()
		
		plinkPedigreeFile = PlinkPedigreeFile(inputFname=self.pedigreeFname)
		for nodeID in plinkPedigreeFile.pedigreeGraph:
			...
		plinkPedigreeFile.close()
	
	

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
import networkx as nx
from pymodule.ProcessOptions import  ProcessOptions
from pymodule import utils, PassingData
from pymodule.algorithm import graph
from MatrixFile import MatrixFile

class PlinkPedigreeFile(MatrixFile):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(MatrixFile.option_default_dict)
	option_default_dict.update({
						('dummyIndividualNamePrefix', 0, ): ['dummy', '', 1, 'the prefix to name a dummy parent (TrioCaller format). The suffix is its order among all dummies.'],\
							
						})
	def __init__(self, inputFname=None, **keywords):
		MatrixFile.__init__(self, inputFname=inputFname, **keywords)
		
		self.familyID2MemberList= {}
		self.familySize2SampleIDList = {}
		
		#2013.07.19
		self._pedigreeGraph = None
		self._childNodeSet = None
	
	@property
	def pedigreeGraph(self):
		"""
		2013.07.19
		"""
		if self._pedigreeGraph is None:
			self.pedigreeGraph = None
		return self._pedigreeGraph
	
	@pedigreeGraph.setter
	def pedigreeGraph(self, argument_ls=[]):
		"""
		2013.07.19
		"""
		graphData = self.getPedigreeGraph()
		self._pedigreeGraph = graphData.DG
		self._childNodeSet = graphData.childNodeSet
	
	def getPedigreeGraph(self,):
		"""
		2013.05.03
		"""
		sys.stderr.write("Getting pedigree graph from this file %s ..."%(self.inputFname))
		
		self._resetInput()
		#self.constructColName2IndexFromHeader()	#there is no header.
		DG = graph.DiGraphWrapper()
		childNodeSet = set()
		counter = 0
		for row in self:
			DG.add_node(row.individualID)	#in case this guy has no parents, then won't be added via add_edge()
			childNodeSet.add(row.individualID)
			if row.paternalID!='0' and row.paternalID.find(self.dummyIndividualNamePrefix)!=0:
				DG.add_edge(row.paternalID, row.individualID)
			if row.maternalID!='0' and row.maternalID.find(self.dummyIndividualNamePrefix)!=0:
				DG.add_edge(row.maternalID, row.individualID)
			counter += 1
		sys.stderr.write("%s children, %s nodes. %s edges. %s connected components.\n"%(\
							len(childNodeSet), DG.number_of_nodes(), DG.number_of_edges(), \
							nx.number_connected_components(DG.to_undirected())))
		return PassingData(DG=DG, childNodeSet=childNodeSet)
	
	def getFamilyStructure(self):
		"""
		2013.07.19
		"""
		sys.stderr.write("Finding unique pairs (singletons or groups) of parents ...\n ")
		noOfParents2FamilyData = {}
		for nodeID in self.pedigreeGraph:
			parents = self.pedigreeGraph.predecessors(nodeID)
			noOfParents = len(parents)
			if noOfParents not in noOfParents2FamilyData:
				noOfParents2FamilyData[noOfParents] = PassingData(parentTupleSet=set(), parentIDSet=set(), childIDSet=set())
			parents.sort()
			noOfParents2FamilyData[noOfParents].parentTupleSet.add(tuple(parents))
			noOfParents2FamilyData[noOfParents].parentIDSet = noOfParents2FamilyData[noOfParents].parentIDSet.union(set(parents))
			noOfParents2FamilyData[noOfParents].childIDSet.add(nodeID)
		
		noOfNuclearFamilies = noOfParents2FamilyData.get(2, 0)
		
		self._reportFamilyStructure(noOfParents2FamilyData)
		return PassingData(noOfParents2FamilyData=noOfParents2FamilyData)
	
	def getFamilyStructurePlinkWay(self, ):
		"""
		2013.07.19
			Plink ignores individuals that do not have independent entries (only show up as parents of others).
			
			Plink Mendel (error detection) works only on families where both parents have their own independent entries in the file.
				This function calculates the number of them (=len(noOfParents2FamilyData[2].parentTupleSet) ).
			
			#founders = number of lines (entries/individuals) where both parents are 0 (do not think in graph way).
			#non-founders = number of lines (entries/individuals), one or both parents are NOT 0.
				non-founders were classified into two
					#non-founders with 2 parents in N1 nuclear families = both their parents are included in the file with independent entries (AKA genotyped).
						= len(noOfParents2FamilyData[2].childIDSet)
						N1 = len(noOfParents2FamilyData[2].parentTupleSet)
						
					#non-founders without 2 parents in N2 nuclear families = the rest
						= len(noOfParents2FamilyData[1].childIDSet) + len(noOfParents2FamilyData[0].childIDSet)
						N2 = len(noOfParents2FamilyData[1].parentTupleSet) + len(noOfParents2FamilyData[0].parentTupleSet)
			
		"""
		sys.stderr.write("Getting number of unique parent-set that both parents ")
		pGraph = self.pedigreeGraph
		self._resetInput()
		individualIDSet = set()
		founderIndividualIDSet = set()
		nonFounderIndividualIDSet = set()
		for row in self:
			if row.paternalID=='0' and row.maternalID=='0':
				founderIndividualIDSet.add(row.individualID)
			else:
				nonFounderIndividualIDSet.add(row.individualID)
			
			individualIDSet.add(row.individualID)
		self._resetInput()
		
		noOfParents2FamilyData = {}
		for nodeID in self.pedigreeGraph:
			if nodeID in individualIDSet and nodeID in nonFounderIndividualIDSet:	#must have an independent entry
					#and exclude founders
				parents = self.pedigreeGraph.predecessors(nodeID)
				parents.sort()
				#calculate no of parents in the plink way, both parents must have independent entries
				noOfParents = 0
				for parentID in parents:
					if parentID in individualIDSet:
						noOfParents += 1
				if noOfParents not in noOfParents2FamilyData:
					noOfParents2FamilyData[noOfParents] = PassingData(parentTupleSet=set(), parentIDSet=set(), childIDSet=set())
				
				parentTupleList = []
				for parentID in parents:
					if parentID!='0':
						parentTupleList.append(parentID)
						noOfParents2FamilyData[noOfParents].parentIDSet.add(parentID)
				
				noOfParents2FamilyData[noOfParents].parentTupleSet.add(tuple(parentTupleList))
				noOfParents2FamilyData[noOfParents].childIDSet.add(nodeID)
		self._reportFamilyStructure(noOfParents2FamilyData)
		sys.stderr.write("\t %s founders, %s non-founders.\n"%(len(founderIndividualIDSet), len(nonFounderIndividualIDSet)))
		return PassingData(noOfParents2FamilyData=noOfParents2FamilyData, founderIndividualIDSet=founderIndividualIDSet,\
						nonFounderIndividualIDSet=nonFounderIndividualIDSet)
	
	def _reportFamilyStructure(self, noOfParents2FamilyData=None):
		"""
		2013.07.19
		"""
		sys.stderr.write("\t%s\t%s\t%s\t%s\n"%("parentSetSize", "noOfFamilies", "noOfParents", "noOfKids"))
		for noOfParents, familyData in noOfParents2FamilyData.iteritems():
			parentIDSet = familyData.parentIDSet
			childIDSet = familyData.childIDSet
			sys.stderr.write("\t%s\t%s\t%s\t%s\n"%(noOfParents, len(familyData.parentTupleSet), len(parentIDSet), len(childIDSet)))
		
	
	def next(self):
		try:
			row = self.csvFile.next()
		except:
			raise StopIteration
		if not self.isRealCSV:
			row = row.strip().split()
		familyID, individualID, paternalID, maternalID, sex, phenotype = row
		return PassingData(familyID=familyID, individualID=individualID, paternalID=paternalID, \
						maternalID=maternalID, sex=sex, phenotype=phenotype)

if __name__ == '__main__':
	main_class = PlinkPedigreeFile
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()