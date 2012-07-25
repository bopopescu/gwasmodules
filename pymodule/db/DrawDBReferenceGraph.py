#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i ~/script/variation/src/Stock_250kDB.py -O /tmp/Stock_250kDB

Description:
	2012.7.9
		This program constructs a database reference graph out of the source code (Elixir db declarative code),
			each node representing one table, each directed edge representing foreign key from one table to another;
			and plots the reference graph in a hierarchical manner (layout=graphviz's dot).
		Each component of the directed graph has its own png output.
		Comments will NOT be skipped in parsing. So watch out.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
import networkx as nx

class DrawDBReferenceGraph(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.pop(('outputFname', 0, ))
	
	def __init__(self,  **keywords):
		"""
		"""
		AbstractMapper.__init__(self, **keywords)
	
	def parseAndGenerateReferenceGraph(self, inputFname=None):
		"""
		2012.7.9
		"""
		sys.stderr.write("Parsing and generating reference graph of %s ..."%(inputFname))
		import networkx as nx
		import re
		tableNamePattern = re.compile(r'^class +(?P<tableName>\w+)\(.*Entity.*\)')
		many2OnePattern = re.compile(r"""ManyToOne\(['"]((%s.)|)(?P<targetTableName>\w+)['"]""")
		many2ManyPattern = re.compile(r"""ManyToMany\(['"]((%s.)|)(?P<targetTableName>\w+)['"].*tablename=["'](?P<middleTableName>\w+)['"],""")
		DG=nx.DiGraph()
		inf = open(inputFname, 'r')
		currentTableName = None
		for line in inf:
			if tableNamePattern.search(line):
				currentTableName = tableNamePattern.search(line).group('tableName')
				DG.add_node(currentTableName)
			elif many2OnePattern.search(line):
				targetTableName = many2OnePattern.search(line).group('targetTableName')
				if not currentTableName or not targetTableName:
					sys.stderr.write("something wrong: currentTableName (%s) is None or empty while targetTableName (%s) is found.\n"%
									(currentTableName, targetTableName))
					sys.exit(3)
				DG.add_edge(currentTableName, targetTableName)
				
			elif many2ManyPattern.search(line):
				targetTableName = many2ManyPattern.search(line).group('targetTableName')
				middleTableName = many2ManyPattern.search(line).group('middleTableName')
				
				if not currentTableName or not targetTableName or not middleTableName:
					sys.stderr.write("something wrong: currentTableName (%s), targetTableName (%s) or middleTableName (%s) is None or empty.\n"%
									(currentTableName, targetTableName, middleTableName))
					sys.exit(3)
				DG.add_edge(middleTableName, targetTableName)
				DG.add_edge(middleTableName, currentTableName)
				
		
		sys.stderr.write("%s nodes. %s edges. %s connected components.\n"%(DG.number_of_nodes(), DG.number_of_edges(), \
															nx.number_connected_components(DG.to_undirected())))
		return DG
	
	def drawTableGraph(self, DG=None, outputFnamePrefix=None):
		"""
		"""
		sys.stderr.write("Drawing the reference graph of all tables in the graph ...\n")
		
		import pylab, numpy
		import matplotlib as mpl
		counter = 0
		for unDirectedSubG in nx.connected_component_subgraphs(DG.to_undirected()):
			#connected_component_subgraphs() can't work on directed graph.
			#nx.draw_circular(DG,with_labels=False, alpha=0.5)
			pylab.clf()
			#get rid of all cushion around the figure
			axe_pvalue = pylab.axes([0, 0, 1.0, 1.0], frameon=False)	#left gap, bottom gap, width, height
			pylab.figure(axe_pvalue.figure.number)
			
			#pylab.axis("off")
			#pylab.figure(figsize=(100, 60))
			layout = 'dot'
			"""
			pos = nx.graphviz_layout(DG, prog=layout)
			nx.draw_networkx_edges(DG, pos, alpha=0.9, width=0.8)
			
			nx.draw_networkx_nodes(DG, pos, alpha=0.9, width=0, linewidths=0.5, cmap=mpl.cm.jet, vmin=0, vmax=1.0)
			"""
			subG = DG.subgraph(unDirectedSubG.nodes())
			#check toplogical sortable
			sys.stderr.write("is component %s a DAG? %s.\n"%(counter, nx.is_directed_acyclic_graph(subG)))
		
			nx.draw_graphviz(subG, prog=layout, with_labels=True,node_size=10, font_size=1, width=0.2, linewidths=0.2, alpha=0.5)
			pylab.savefig('%s_graphviz_%s_component_%s.png'%(outputFnamePrefix, layout, counter), dpi=700)
			counter += 1
		sys.stderr.write(".\n")
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		DG = self.parseAndGenerateReferenceGraph(inputFname=self.inputFname)
		self.drawTableGraph(DG, outputFnamePrefix = self.outputFnamePrefix)

if __name__ == '__main__':
	main_class = DrawDBReferenceGraph
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()