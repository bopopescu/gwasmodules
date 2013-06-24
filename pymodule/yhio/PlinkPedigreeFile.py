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
	
		reader = MatrixFile(inputFname='/tmp/input.txt', openMode='r')
		reader = MatrixFile('/tmp/input.txt', openMode='r')
		reader.constructColName2IndexFromHeader()
		for row in reader:
			row[reader.getColName2IndexFromHeader('KID')]
		
		inf = utils.openGzipFile(inputFname, openMode='r')
		reader = MatrixFile(inputFile=inf)
		
		#2013.2.1 writing
		writer = MatrixFile('/tmp/output.txt', openMode='w', delimiter='\t')
		writer.writeHeader(...)
		writer.writerow(row)
		writer.close()
	
	

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
	
	def getPedigreeGraph(self,):
		"""
		2013.05.03
		"""
		sys.stderr.write("Getting pedigree graph from this file %s ..."%(self.inputFname))
		
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