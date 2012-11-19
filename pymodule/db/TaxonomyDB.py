#!/usr/bin/env python
"""
Examples:
	#setup database in postgresql
	TaxonomyDB.py -u crocea -k taxonomy
	
	# 2012.6.6 setup taxonomy schema in vervetdb.
	TaxonomyDB.py -u yh -k taxonomy -d vervetdb -v postgresql
	
	#setup database in mysql
	%s -v mysql -u yh -z papaya -d taxonomy -k ""
	
	# 2012.6.6 setup & import NCBI tax dump into db
	%s -u yh -k taxonomy -d vervetdb -v postgresql -i /usr/local/research_data/NCBI/taxonomy
	
Description:
	2012.6.6 (funcitonal now)
	This is a wrapper for the taxonomy database (taxdump in ftp://ftp.ncbi.nih.gov/pub/taxonomy/), built on top of elixir.
	It also imports the taxdump into db if given a folder of un-tarred dump files (-i).
"""
import sys, os
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

from sqlalchemy.engine.url import URL
from elixir import Unicode, DateTime, String, Integer, UnicodeText, Text
from elixir import Entity, Field, using_options, using_table_options
from elixir import OneToMany, ManyToOne, ManyToMany, OneToOne
from elixir import setup_all, session, metadata, entities
from elixir.options import using_table_options_handler	#using_table_options() can only work inside Entity-inherited class.
from datetime import datetime
from sqlalchemy.schema import ThreadLocalMetaData, MetaData
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy import UniqueConstraint, create_engine
from sqlalchemy import and_, or_, not_


from db import ElixirDB

__session__ = scoped_session(sessionmaker(autoflush=False, autocommit=True))
#__metadata__ = ThreadLocalMetaData() #2008-11-04 not good for pylon

__metadata__ = MetaData()

class Citation(Entity):
	"""
	2012.6.6
	"""
	cit_id = Field(Integer)	#- the unique id of citation
	cit_key = Field(Text)	# -- citation key
	pubmed_id = Field(Integer)	#  -- unique id in PubMed database (0 if not in PubMed)
	medline_id = Field(Integer)		#  -- unique id in MedLine database (0 if not in MedLine)
	url = Field(Text)	# -- URL associated with citation
	text = Field(Text)	#-- any text (usually article name and authors).
		#-- The following characters are escaped in this text by a backslash:
		#-- newline (appear as "\n"),
		#-- tab character ("\t"),
		#-- double quotes ('\"'),
		#-- backslash character ("\\").
	tax_list = ManyToMany("Node", tablename='citation2node', local_colname='citation_id')
		# -- list of node ids separated by a single space
	created_by = Field(String(256))
	updated_by = Field(String(256))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='citation')
	using_table_options(mysql_engine='InnoDB')

class Delnode(Entity):
	"""
	2012.6.6
	"""
	tax_id = Field(Integer, unique=True)	# -- taxonomy database division id
	created_by = Field(String(256))
	updated_by = Field(String(256))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='delnode')
	using_table_options(mysql_engine='InnoDB')

class Division(Entity):
	"""
	2012.6.6
	"""
	id = Field(Integer, primary_key=True)	# -- taxonomy database division id
	code = Field(Text)	# -- GenBank division code (three characters)
		# e.g. BCT, PLN, VRT, MAM, PRI...
	name = Field(Text)
	comments = Field(Text)
	created_by = Field(String(256))
	updated_by = Field(String(256))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='division')
	using_table_options(mysql_engine='InnoDB')

class Gencode(Entity):
	"""
	"""
	id = Field(Integer, primary_key=True)	# -- GenBank genetic code id
	abbreviation = Field(Text)	# -- genetic code name abbreviation
	name = Field(Text)	# -- genetic code name
	code = Field(Text)	#-- translation table for this genetic code
	starts = Field(Text)	#	-- start codons for this genetic code
	comment = Field(Text)
	created_by = Field(String(256))
	updated_by = Field(String(256))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='gencode')
	using_table_options(mysql_engine='InnoDB')

class Merged(Entity):
	"""
	2012.6.6
	"""
	old_tax_id = Field(Integer)		#-- id of nodes which has been merged
	new_tax = ManyToOne('%s.Node'%__name__, colname='new_tax_id', ondelete='CASCADE', onupdate='CASCADE')
		#-- id of nodes which is result of merging
	created_by = Field(String(256))
	updated_by = Field(String(256))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='merged')
	using_table_options(mysql_engine='InnoDB')


class Name(Entity):
	"""
	2012.6.6
		taxonomy name
	"""
	tax = ManyToOne('%s.Node'%__name__, colname='tax_id', ondelete='CASCADE', onupdate='CASCADE')
	name_txt = Field(Text)
	unique_name = Field(Text)
	name_class = Field(Text)
	created_by = Field(String(256))
	updated_by = Field(String(256))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='name')
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('tax_id','name_txt','name_class'))
	

class Node(Entity):
	"""
	2012.6.6
		taxonomy node and its parent node
	"""
	tax_id = Field(Integer, unique=True, primary_key=True)	#-- node id in GenBank taxonomy database
	parent_tax_id = Field(Integer)	#2012.6.7 some parent taxonomy node might not exist at all.
	#parent_tax = ManyToOne('%s.Node'%__name__, colname='parent_tax_id', ondelete='CASCADE', onupdate='CASCADE')
	rank = Field(Text)	#-- rank of this node (superkingdom, kingdom, ...)
	embl_code = Field(Text)	#-- locus-name prefix; not unique
	division = ManyToOne('%s.Division'%__name__, colname='division_id', ondelete='CASCADE', onupdate='CASCADE')
		#see division.dmp file
	inherited_div_flag = Field(Integer)	#  (1 or 0) -- 1 if node inherits division from parent
	genetic_code = ManyToOne('%s.Gencode'%__name__, colname='genetic_code_id', ondelete='CASCADE', onupdate='CASCADE')
		# -- see gencode.dmp file
	inherited_GC_flag = Field(Integer)	#(1 or 0) -- 1 if node inherits genetic code from parent
	mitochondrial_genetic_code = ManyToOne('%s.Gencode'%__name__, colname='mitochondrial_genetic_code_id', \
										ondelete='CASCADE', onupdate='CASCADE')
		# -- see gencode.dmp file
	inherited_MGC_flag = Field(Integer) #(1 or 0)	#  -- 1 if node inherits mitochondrial gencode from parent
	GenBank_hidden_flag = Field(Integer)	# (1 or 0) -- 1 if name is suppressed in GenBank entry lineage
	hidden_subtree_root_flag = Field(Integer)		#(1 or 0) -- 1 if this subtree has no sequence data yet
	comments = Field(Text)
	citation_list = ManyToMany("Citation", tablename='citation2node', local_colname='tax_id')
	created_by = Field(String(256))
	updated_by = Field(String(256))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='node')
	using_table_options(mysql_engine='InnoDB')

class TaxonomyDB(ElixirDB):
	__doc__ = __doc__
	option_default_dict = ElixirDB.option_default_dict.copy()
	option_default_dict[('drivername', 1,)][0] = 'postgresql'
	option_default_dict[('database', 1,)][0] = 'taxonomy'
	option_default_dict.update({
			('inputFolder', 0, ):[None, 'i', 1, 'where taxdump.tar.gz was un-tarred. it should contain citations.dmp, division.dmp, \n\
	gencode.dmp, names.dmp, delnodes.dmp, merged.dmp, nodes.dmp. \n\
	If this argument is provided, this program will import them into a taxonomy database.\n'],\
			('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
			})	
	def __init__(self, **keywords):
		"""
		2008-10-08
			simplified further by moving db-common lines to ElixirDB
		2008-07-09
		"""
		from ProcessOptions import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		self.setup_engine(metadata=__metadata__, session=__session__, entities=entities)
		self.ncbiTaxDumpFileDelimiter = '\t|\t'
		
		self._scientific_name2tax_id = None
		self._tax_id2scientific_name = None
	
	def splitLine(self, line=None):
		line = line.strip()
		if line[-2:]=='\t|':	#the last field has this ending
			line = line[:-2]
		row = line.split(self.ncbiTaxDumpFileDelimiter)
		return row
	
	def importMergedTaxNodes(self, inputFile=None):
		"""
		2012.6.6
		"""
		sys.stderr.write('Importing merged taxonomy nodes ...')
		counter = 0
		for line in inputFile:
			row = self.splitLine(line)
			old_tax_id = int(row[0])
			new_tax_id = int(row[1])
			db_entry = Merged(old_tax_id=old_tax_id, new_tax_id=new_tax_id)
			self.session.add(db_entry)
			counter += 1
		self.session.flush()
		sys.stderr.write("%s merged nodes added into db.\n"%(counter))
	
	def importTaxNodes(self, inputFile=None, divisionID2division=None, gencodeID2gencode=None):
		"""
		2012.6.6
			create a structure of parentTaxonomy2Child first to insure parent taxonomy is in db first
		"""
		sys.stderr.write('Importing taxonomy nodes ...\n')
		counter = 0
		import networkx as nx
		DG = nx.DiGraph()	#use this to store the relationship among them
		DG.add_weighted_edges_from([(1,2,0.5), (3,1,0.75)])
		row_ls = []
		taxID2rowIndex = {}
		for line in inputFile:
			row = self.splitLine(line)
			
			tax_id = int(row[0])
			parent_tax_id = int(row[1])
			taxID2rowIndex[tax_id] = len(row_ls)
			row_ls.append(row)
			if parent_tax_id==tax_id:	#the first taxonomy node has itself as the parent, which makes the graph not a DAG.
				sys.stderr.write("Warning: taxonomy ID %s's parent is itself. Not added in the directed graph.\n"%(tax_id))
			else:
				DG.add_edge(parent_tax_id, tax_id)
		G = DG.to_undirected()
		sys.stderr.write(" %s nodes, %s edges, %s connected components.\n"%(DG.number_of_nodes(), DG.number_of_edges(),\
															nx.number_connected_components(G)))

		#first check if this is DAG , it should be.
		if not nx.is_directed_acyclic_graph(DG):
			sys.stderr.write("Error: tax node graph is not a directed acyclic graph. it should be.\n")
			#sys.exit(3)
		# then save the nodes in the topological sort order (parents get saved before children).
		taxIDInTopoSortOrder = nx.topological_sort(DG)
		for tax_id in taxIDInTopoSortOrder:
			if tax_id not in taxID2rowIndex:
				sys.stderr.write("Warning: taxonomy ID %s has no real data in the input file.\n"%(tax_id))
				continue
			rowIndex = taxID2rowIndex.get(tax_id)
			row = row_ls[rowIndex]
			tax_id = int(row[0])
			parent_tax_id = int(row[1])
			rank = row[2]
			embl_code = row[3]
			division_id = int(row[4])
			
			division = divisionID2division[division_id]
			
			inherited_div_flag = int(row[5])
			genetic_code_id = int(row[6])
			genetic_code = gencodeID2gencode[genetic_code_id]
			
			inherited_GC_flag = int(row[7])
			mitochondrial_genetic_code_id = int(row[8])
			mitochondrial_genetic_code = gencodeID2gencode[mitochondrial_genetic_code_id]
			
			inherited_MGC_flag = int(row[9])
			GenBank_hidden_flag = int(row[10])
			hidden_subtree_root_flag = int(row[11])
			comments = row[12].decode('latin1')
			
			node = Node(tax_id=tax_id, parent_tax_id=parent_tax_id, rank=rank, embl_code=embl_code, inherited_div_flag=inherited_div_flag,\
					inherited_GC_flag=inherited_GC_flag, inherited_MGC_flag=inherited_MGC_flag, GenBank_hidden_flag=GenBank_hidden_flag,\
					hidden_subtree_root_flag=hidden_subtree_root_flag, comments=comments)
			node.division = division
			node.genetic_code = genetic_code
			node.mitochondrial_genetic_code = mitochondrial_genetic_code
			self.session.add(node)
			counter += 1
			if counter%5000==0:
				sys.stderr.write("%s%s"%('\x08'*40, counter))
				self.session.flush()
		sys.stderr.write("%s%s\n"%('\x08'*40, counter))
		self.session.flush()
		sys.stderr.write("%s nodes added into db.\n"%(counter))
	
	def importDeletedNodes(self, inputFile=None):
		"""
		2012.6.6
		"""
		sys.stderr.write('Importing deleted nodes ...\n')
		counter = 0
		for line in inputFile:
			row = self.splitLine(line)
			tax_id = int(row[0])
			db_entry = Delnode(tax_id=tax_id)
			self.session.add(db_entry)
			self.session.flush()
			counter += 1
			if counter%5000==0:
				sys.stderr.write("%s%s"%('\x08'*40, counter))
		sys.stderr.write("%s%s\n"%('\x08'*40, counter))
		sys.stderr.write("%s added into db.\n"%(counter))
	
	def importTaxNames(self, inputFile=None):
		"""
		2012.6.6
		"""
		sys.stderr.write('Importing taxonomy names ...\n')
		counter = 0
		for line in inputFile:
			row = self.splitLine(line)
			tax_id = int(row[0])
			name_txt = row[1]
			unique_name = row[2]
			name_class = row[3]
			db_entry = Name(tax_id=tax_id, name_txt=name_txt, unique_name=unique_name, name_class=name_class)
			self.session.add(db_entry)
			counter += 1
			if counter%5000==0:
				sys.stderr.write("%s%s"%('\x08'*40, counter))
				self.session.flush()
		sys.stderr.write("%s%s\n"%('\x08'*40, counter))
		self.session.flush()
		sys.stderr.write("%s added into db.\n"%(counter))
	
	def importGeneticCode(self, inputFile=None):
		"""
		2012.6.6
		"""
		sys.stderr.write('Importing genetic codes ...')
		counter = 0
		gencodeID2gencode = {}
		for line in inputFile:
			row = self.splitLine(line)
			gencode_id = int(row[0])
			abbr = row[1]
			name = row[2]
			code = row[3]
			starts = row[4]
			gencode = Gencode(id=gencode_id, abbreviation=abbr, name=name, code=code, starts=starts)
			self.session.add(gencode)
			self.session.flush()
			gencodeID2gencode[gencode.id] = gencode
			counter += 1
		sys.stderr.write("%s added into db.\n"%(counter))
		return gencodeID2gencode
	
	
	def importDivision(self, inputFile=None):
		"""
		2012.6.6
		"""
		sys.stderr.write('Importing divisions ...')
		counter = 0
		divisionID2division = {}
		for line in inputFile:
			row = self.splitLine(line)
			division = Division(id=int(row[0]), code=row[1], name=row[2], comments=row[3])
			self.session.add(division)
			self.session.flush()
			divisionID2division[division.id] = division
			counter += 1
		sys.stderr.write("%s added into db.\n"%(counter))
		return divisionID2division
	
	def importCitation(self, inputFile=None):
		"""
		2012.6.6
		"""
		sys.stderr.write('Importing citations ...\n')
		counter = 0
		for line in inputFile:
			row = self.splitLine(line)
			cit_id = int(row[0])
			cit_key = row[1].decode('latin1')
			pubmed_id = int(row[2])
			medline_id = int(row[3])
			url = row[4].decode('latin1')
			if len(row)>5:
				text = row[5].decode('latin1')	#the input file contains utf-8 characters beyond 127 ('ascii' codec limit). 
				#.encode('ascii','replace')
			else:
				text = None
			if len(row)>6:
				taxIDList = row[6].split()
				taxIDList = map(int, taxIDList)
			else:
				taxIDList = []
			db_entry = Citation(cit_id=cit_id, cit_key=cit_key, pubmed_id=pubmed_id, medline_id=medline_id, url=url,\
							text=text)
			for taxID in taxIDList:
				tax = Node.get(taxID)
				if tax is None:
					sys.stderr.write("Error: taxonomy ID %s not in table node.\n"%(taxID))
					sys.exit(2)
				db_entry.tax_list.append(tax)
			self.session.add(db_entry)
			counter += 1
			if counter%5000==0:
				sys.stderr.write("%s%s"%('\x08'*40, counter))
				self.session.flush()
		sys.stderr.write("%s%s\n"%('\x08'*40, counter))
		self.session.flush()
		sys.stderr.write("%s added into db.\n"%(counter))
	
	
	def importNCBITaxDump(self, inputFolder=None, ):
		"""
		2012.6.6
			
		"""
		sys.stderr.write("Importing NCBI tax dump files from %s ...\n"%(inputFolder))
		citationsFile = open(os.path.join(inputFolder, 'citations.dmp'), 'r')
		divisionFile = open(os.path.join(inputFolder, 'division.dmp'), 'r')
		gencodeFile = open(os.path.join(inputFolder, 'gencode.dmp'), 'r')
		namesFile = open(os.path.join(inputFolder, 'names.dmp'), 'r')
		delnodesFile = open(os.path.join(inputFolder, 'delnodes.dmp'), 'r')
		mergedFile = open(os.path.join(inputFolder, 'merged.dmp'), 'r')
		nodesFile = open(os.path.join(inputFolder, 'nodes.dmp'), 'r')
		openedFileList = [citationsFile, divisionFile, gencodeFile, namesFile, delnodesFile, mergedFile, nodesFile]
		
		gencodeID2gencode = self.importGeneticCode(gencodeFile)
		divisionID2division = self.importDivision(divisionFile)
		self.importTaxNodes(nodesFile, gencodeID2gencode=gencodeID2gencode, divisionID2division=divisionID2division)
		
		self.importTaxNames(namesFile)
		self.importMergedTaxNodes(mergedFile)
		self.importDeletedNodes(delnodesFile)
		self.importCitation(citationsFile)
		
		for openedFile in openedFileList:
			openedFile.close()
		
		sys.stderr.write("\n")
	
	def run(self):
		"""
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
	
		if self.inputFolder:	#only pick the file if the output file is not empty.
			self.session.begin()
			self.importNCBITaxDump(self.inputFolder)
			if self.commit:
				self.session.flush()
				self.session.commit()
			else:
				self.session.rollback()
	
	@property
	def scientific_name2tax_id(self):
		"""
		2012.8.28
			copied from pymodule/utils.py
		2012.6.6
			update it to get table names from TaxonomyDB
		"""
		if self._scientific_name2tax_id is None:
			self.fillUpMapBetweenTaxIDAndScientificName()
		return self._scientific_name2tax_id
	
	def returnTaxIDGivenScientificName(self, scientific_name=None):
		"""
		2012.8.28
			copied from pymodule/utils.py
		"""
		return self.scientific_name2tax_id.get(scientific_name)
	
	def returnTaxIDGivenSentence(self, sentence=None):
		"""
		2012.8.28
			copied from pymodule/utils.py
		2008-07-29
		"""
		tax_id_to_return = None
		for scientific_name, tax_id in self.scientific_name2tax_id.iteritems():
			if sentence.find(scientific_name)>=0:
				tax_id_to_return = tax_id
				break
		return tax_id_to_return
	
	def fillUpMapBetweenTaxIDAndScientificName(self):
		"""
		2012.8.29
			remove " and o.rank='species'" from the query condition
		2012.8.29
		"""
		self._scientific_name2tax_id = {}
		self._tax_id2scientific_name = {}
		curs = self.metadata.bind
		rows = curs.execute("SELECT n.name_txt, n.tax_id FROM %s.%s n, %s.%s o where n.name_class='scientific name' \
				and n.tax_id=o.tax_id"%(self.schema, Name.table.name, \
															self.schema, Node.table.name))
		#rows = curs.fetchall()
		for row in rows:
			scientific_name = row.name_txt
			tax_id = row.tax_id
			self._tax_id2scientific_name[tax_id] = scientific_name
			self._scientific_name2tax_id[scientific_name] = tax_id
		sys.stderr.write("%s entries in _tax_id2scientific_name. %s entries in _scientific_name2tax_id.\n"%\
						(len(self._tax_id2scientific_name), len(self._scientific_name2tax_id)))

	@property
	def tax_id2scientific_name(self):
		"""
		2012.8.28
		"""
		if self._tax_id2scientific_name is None:
			self.fillUpMapBetweenTaxIDAndScientificName()
		return self._tax_id2scientific_name
	
	def returnScientificNameGivenTaxID(self, tax_id=None):
		"""
		2012.8.28
			there's a cache involved in the 1st query.
		"""
		tax_id2scientific_name = self.tax_id2scientific_name
		return tax_id2scientific_name.get(tax_id)
	
	def getScientificNameGivenTaxID(self, tax_id=None):
		"""
		2012.8.28
			no cache, jsut db query
		"""
		curs = self.metadata.bind
		rows = curs.execute("SELECT n.name_txt, n.tax_id FROM %s.%s n, %s.%s o where n.name_class='scientific name' \
				and n.tax_id=o.tax_id and n.tax_id=%s"%(self.schema, Name.table.name, \
																	self.schema, Node.table.name, tax_id))
		row = rows.fetchone()
		if row:
			scientific_name = row.name_txt
			return scientific_name
		else:
			return None
	
	def getTaxIDGivenScientificName(self, scientific_name=None):
		"""
		2012.8.29
			no cache, just db query
		"""
		curs = self.metadata.bind
		rows = curs.execute("SELECT n.name_txt, n.tax_id FROM %s.%s n, %s.%s o where n.name_class='scientific name' \
				and n.tax_id=o.tax_id and n.name_txt='%s'"%(self.schema, Name.table.name,\
																		self.schema, Node.table.name, scientific_name))
		row = rows.fetchone()
		if row:
			return row.tax_id
		else:
			return None
	
if __name__ == '__main__':
	import sys, os, math
	bit_number = math.log(sys.maxint)/math.log(2)
	if bit_number>40:       #64bit
		sys.path.insert(0, os.path.expanduser('~/lib64/python'))
		sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
	else:   #32bit
		sys.path.insert(0, os.path.expanduser('~/lib/python'))
		sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

	from pymodule import ProcessOptions
	main_class = TaxonomyDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.setup()
	instance.run()
