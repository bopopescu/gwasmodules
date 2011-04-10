#!/usr/bin/env python
"""

Examples:
	#postgresql
	TAIRGeneXML2GenomeDB.py -c -i *.xgs
	
	TAIRGeneXML2GenomeDB.py -z localhost -k genome -c -i TAIR9_genome_release/Tair9_XML/ch1.xml -r

	#mysql
	TAIRGeneXML2GenomeDB.py -v mysql -u yh -z papaya -d genome_tair -k "" -c -i TAIR9_genome_release/Tair9_XML/ch1.xml -r
	
	TAIRGeneXML2GenomeDB.py -v mysql -u yh -z banyan -d genome_tair -r  ~/script/variation/data/TAIR9/TAIR9_genome_release/Tair9_XML/ch[2-5].*
	
Description:
	Put TAIR xml genome release files into db.
	This program checks whether an entry is in db or not only for AnnotAssembly, Gene, EntrezgeneMapping, EntrezgeneType,
		Gene2go, GeneCommentaryType, but not for GeneCommentary and GeneSegment.
	If program stops half way and needs to be re-run, all prior inserted data has to be purged MANUALLY.
	
	xml files are downloaded from ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR9_genome_release/Tair9_XML/
"""

import sys, getopt, os
sys.path += [os.path.join(os.path.expanduser('~/script/annot/bin'))]
from codense.common import db_connect
import xml.etree.cElementTree as ElementTree
from GenomeDB import GenomeDatabase, Gene, EntrezgeneMapping, SequenceType, EntrezgeneType, \
	GeneSegment, GeneCommentaryType, GeneCommentary, AnnotAssembly, Gene2go
from datetime import datetime
from pymodule import PassingData
import traceback
import datetime
#import argparse	#in python 2.7
from optparse import OptionParser

def usernameOptionCallBack(option, opt_str, value, parser):
	#2010-8-13
	# useless. called only when -u is supplied, which misses the whole point.
	option_name = option.dest
	if not value:
		value = raw_input("%s: "%opt_str)
	setattr(parser.values, option_name, value)

def passwordOptionCallBack(option, opt_str, value, parser):
	option_name = option.dest
	if not value:
		value = getpass.getpass("%s: "%opt_str)
	setattr(parser.values, option_name, value)


class TAIRGeneXML2GenomeDB:
	__doc__ = __doc__
	
	usage = "usage: %prog [options] inputfiles \n\n-h/--help for help."
	parser = OptionParser(usage=usage, description='', epilog=__doc__)
	#parser.add_option('inputfiles', type=int, nargs='+',
	#				help='input files in XML. can detect gzipped file by .gz in filename.')
	parser.add_option('-v', '--drivername', default="postgres",
					help='which type of database? mysql or postgres [%default]')
	parser.add_option('-z', '--hostname', dest='hostname', default="localhost",
					help='hostname of the db server [%default]')
	parser.add_option('-d', '--dbname', dest='dbname', default="graphdb",
					help='database name [%default]')
	parser.add_option('-k', '--schema', dest='schema', default="",
					help='database schema (only for postgres) [%default]')
	parser.add_option('-u', '--db_user', dest='db_user', 
					help='database username [%default]')
	parser.add_option('-p', '--db_passwd', 
					help='database password [%default]')
	parser.add_option('-x', '--tax_id', default=3702, type=int,
					help='taxonomy ID for the organism [%default]')
	parser.add_option('-c', '--commit', default=False, action='store_true',
					help='commit db transaction. The program commits regardless of this argument.[%default]')
	parser.add_option('-b', '--debug', default=False, action='store_true',
					help='toggle debug mode [%default]')
	parser.add_option('-r', '--report', default=False, action='store_true',
					help='toggle report, more verbose stdout/stderr. [%default]')
	
	option_default_dict = {('drivername', 1,):['postgres', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['graphdb', 'd', 1, 'database name', ],\
							('schema', 0, ): ['', 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('tax_id', 1, int): [3702, '', 1, 'taxonomy ID', ],\
							('inputfiles', 1, ):[None, 'i', 1, 'comma-separated input filenames, or a list of files'],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self,  argv):
		"""
		2008-07-29
			use option_default_dict
		"""
		(options, args) = self.parser.parse_args(argv)
		for option in self.parser.option_list:
			option_name = option.dest
			if option_name:
				setattr(self, option_name, getattr(options, option_name))
		self.inputfiles = args[1:]
		if len(self.inputfiles) <1:
			self.parser.error("Please supply input files.\n")
		
		#from pymodule import ProcessOptions
		#self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		#if type(self.inputfiles)==str:
		#	self.inputfiles = self.inputfiles.split(',')
	
	def getGeneCommentaryType(self, session, type_name=None, type_id=None):
		"""
		2010-8-12
		"""
		if type_id is not None:
			gene_commentary_type = GeneCommentaryType.get(type_id)
		else:
			gene_commentary_type = GeneCommentaryType.query.filter_by(type=type_name).first()
		if not gene_commentary_type:
			if self.debug:
				sys.stderr.write("\t Gene-commentary_type=%s not in db yet.\n"%type_name)
			gene_commentary_type = GeneCommentaryType(type=type_name)
			if type_id is not None:
				gene_commentary_type.id = type_id
			session.add(gene_commentary_type)
			session.flush()
		return gene_commentary_type
	
	def addGeneCommentary(self, element, entrezgene_mapping, session, param_obj=None):
		"""
		2010-8-12
		"""
		if not hasattr(param_obj, 'no_of_total'):
			setattr(param_obj, 'no_of_total', 0)
		if not hasattr(param_obj, 'no_of_into_db'):
			setattr(param_obj, 'no_of_into_db', 0)
		
		model_start = int(element.findtext('COORDSET/END5'))	#2010-8-18 this usually corresponds to the CDS start-stop, not EXONs
		model_stop = int(element.findtext('COORDSET/END3'))
		pub_locus = element.findtext('PUB_LOCUS')
		
		
		CDS_box_ls = []
		EXON_box_ls = []
		UTR_box_ls = []
		for exon_element in element.findall('EXON'):
			#gene segment for this exon
			start = int(exon_element.findtext('COORDSET/END5'))
			stop = int(exon_element.findtext('COORDSET/END3'))
			EXON = (min(start, stop), max(start, stop))
			EXON_box_ls.append(EXON)
			
			# if CDS exist, get all boxes, construct another gene commentary, and each CDS is a gene segment
			for cds_element in exon_element.findall("CDS"):
				start = int(cds_element.findtext('COORDSET/END5'))
				stop = int(cds_element.findtext('COORDSET/END3'))
				CDS = (min(start, stop), max(start, stop))	#always smaller as start. strand information in EntrezgeneMapping
				CDS_box_ls.append(CDS)
			utr_elements = exon_element.find("UTRS")
			if utr_elements:
				for utr_element in utr_elements:
					start = int(utr_element.findtext('COORDSET/END5'))
					stop = int(utr_element.findtext('COORDSET/END3'))
					UTR = (min(start, stop), max(start, stop))	#always smaller as start. strand information in EntrezgeneMapping
					UTR_box_ls.append(UTR)
		
		EXON_box_ls.sort()
		if len(EXON_box_ls)>0:	#2010-8-18 replace the CDS start-stop with EXON start-stop.
			model_start = EXON_box_ls[0][0]
			model_stop = EXON_box_ls[-1][1]
		gene_commentary = GeneCommentary(label=pub_locus, \
						text=entrezgene_mapping.gene.description, start=model_start, stop=model_stop)
		# start is always smaller than stop. Strand information is in table EntrezgeneMapping
		gene_commentary.gene = entrezgene_mapping
		
		commentary_type_name = entrezgene_mapping.entrezgene_type.type
		
		gene_commentary_type = self.getGeneCommentaryType(session, type_name=commentary_type_name)
		
		gene_commentary.gene_commentary_type = gene_commentary_type
		session.add(gene_commentary)
		
		# 2010-8-18 get all the introns
		no_of_exons = len(EXON_box_ls)
		for i in xrange(no_of_exons):
			EXON = EXON_box_ls[i]
			start, stop = EXON[:2]
			gene_segment = GeneSegment(start=min(start, stop), stop=max(start, stop))
			gene_segment.gene_commentary = gene_commentary
			gene_commentary_type = self.getGeneCommentaryType(session, type_name='exon')
			gene_segment.gene_commentary_type = gene_commentary_type
			session.add(gene_segment)
			param_obj.no_of_into_db += 1
			
			if i>0:
				previous_EXON = EXON_box_ls[i-1]
				current_EXON = EXON_box_ls[i]
				intron_start = previous_EXON[1]+1
				intron_stop = current_EXON[0]-1
				
				gene_segment = GeneSegment(start=intron_start, stop=intron_stop)
				gene_segment.gene_commentary = gene_commentary
				gene_commentary_type = self.getGeneCommentaryType(session, type_name='intron')
				gene_segment.gene_commentary_type = gene_commentary_type
				session.add(gene_segment)
				param_obj.no_of_into_db += 1
		
		
		CDS_box_ls.sort()
		if CDS_box_ls:
			whole_CDS_start = CDS_box_ls[0][0]
			whole_CDS_stop = CDS_box_ls[-1][1]
			protein_commentary = GeneCommentary(label=pub_locus, \
						start=whole_CDS_start, stop=whole_CDS_stop)
			# start is always smaller than stop. Strand information is in table EntrezgeneMapping
			protein_commentary.gene = entrezgene_mapping
			protein_commentary.gene_commentary = gene_commentary
			
			gene_commentary_type = self.getGeneCommentaryType(session, type_name="peptide")
			
			protein_commentary.gene_commentary_type = gene_commentary_type
			session.add(protein_commentary)
			param_obj.no_of_into_db += 1
			
			for CDS in CDS_box_ls:
				gene_segment = GeneSegment(start=CDS[0], stop=CDS[1])
				gene_segment.gene_commentary = protein_commentary
				gene_commentary_type = self.getGeneCommentaryType(session, type_name='CDS')
				gene_segment.gene_commentary_type = gene_commentary_type
				session.add(gene_segment)
				param_obj.no_of_into_db += 1
			
			if UTR_box_ls:
				#2010-8-18 add all UTR boxes to this protein_commentary
				for UTR in UTR_box_ls:
					if entrezgene_mapping.strand =='+1':
						if UTR[1]<whole_CDS_start:
							type_name = '5UTR'
						elif UTR[0]>whole_CDS_stop:
							type_name = '3UTR'
						else:
							pass	#shouldn't be here
					else:
						if UTR[1]<whole_CDS_start:
							type_name = '3UTR'
						elif UTR[0]>whole_CDS_stop:
							type_name = '5UTR'
						else:
							pass	#shouldn't be here
					gene_segment = GeneSegment(start=UTR[0], stop=UTR[1])
					gene_segment.gene_commentary = protein_commentary
					gene_commentary_type = self.getGeneCommentaryType(session, type_name=type_name)
					gene_segment.gene_commentary_type = gene_commentary_type
					session.add(gene_segment)
					param_obj.no_of_into_db += 1
	
	def addGene(self, session, annot_assembly, element, type_of_gene=None, param_obj=None):
		"""
		2010-8-11
		"""
		if not hasattr(param_obj, 'no_of_total'):
			setattr(param_obj, 'no_of_total', 0)
		if not hasattr(param_obj, 'no_of_into_db'):
			setattr(param_obj, 'no_of_into_db', 0)
		
		gene_synonym = element.find("GENE_SYNONYM")
		feat_name = element.findtext("FEAT_NAME")
		com_name = element.findtext("COM_NAME")
		start = int(element.findtext('COORDSET/END5'))
		stop = int(element.findtext('COORDSET/END3'))
		
		date_string = element.findtext("DATE")
		try:	#2010-8-12
			modification_date = datetime.datetime.strptime(date_string, '%b %d %Y %I:%M%p')
		except:
			modification_date = None
		
		if start<=stop:
			strand = "+1"
		else:
			# reverse start and stop
			strand = "-1"
			tmp = start
			start = stop
			stop = tmp
		
		gene_info = element.find("GENE_INFO")
		
		locustag = gene_info.findtext("PUB_LOCUS")
		gene_symbol = gene_info.findtext("GENE_SYM")
		if not gene_symbol:
			gene_symbol = locustag
		com_name = gene_info.findtext("COM_NAME")
		comment = gene_info.findtext("COMMENT")
		pub_comment = gene_info.findtext("PUB_COMMENT")
		if comment:
			pub_comment += "\n %s"%comment
		if type_of_gene=='protein-coding':
			is_pseudogene = gene_info.findtext("IS_PSEUDOGENE")
			if is_pseudogene=='1':
				type_of_gene = 'pseudo'
		else:
			type_of_gene = element.tag	#non-protein-coding genes have their tags as gene types (like rRNA,)
		
		gene = Gene.query.filter_by(tax_id=annot_assembly.tax_id).filter_by(chromosome=annot_assembly.chromosome).\
			filter_by(locustag=locustag).filter_by(strand=strand).filter_by(start=start).filter_by(stop=stop).\
			filter_by(type_of_gene=type_of_gene).first()
		if not gene:
			gene = Gene(gene_symbol=gene_symbol, tax_id=annot_assembly.tax_id, locustag=locustag, \
				chromosome=annot_assembly.chromosome,\
				description = pub_comment, full_name_from_nomenclature_authority=com_name, type_of_gene=type_of_gene, \
				modification_date=modification_date, strand=strand, start=start, stop=stop)
			session.add(gene)
		else:
			param_obj.no_of_genes_already_in_db += 1
		
		entrezgene_type = EntrezgeneType.query.filter_by(type=type_of_gene).first()	#query the db to see if it exists or not
		if not entrezgene_type:
			entrezgene_type = EntrezgeneType(type=type_of_gene)
			session.add(entrezgene_type)
			session.flush()
			param_obj.no_of_into_db += 1
		
		gene_ontology  = gene_info.find("GENE_ONTOLOGY")
		if gene_ontology:
			for go_id_elem in gene_ontology.findall('GO_ID'):
				go_id = go_id_elem.get("ASSIGNMENT")
				evidence = go_id_elem.find("GO_EVIDENCE/EV_CODE").get("CODE")
				evidence_source = go_id_elem.findtext("GO_EVIDENCE/EVIDENCE")
				go_term = go_id_elem.findtext("GO_TERM")
				category = go_id_elem.findtext("GO_TYPE").split()[1]	#component in "cellular component" is enough
				if gene.gene_id:
					gene2go = Gene2go.query.filter_by(gene_id=gene.gene_id).filter_by(go_id=go_id).filter_by(evidence=evidence).\
						filter_by(category=category).filter_by(tax_id=annot_assembly.tax_id).first()
				else:
					gene2go = None
				if not gene2go:
					gene2go = Gene2go(tax_id=annot_assembly.tax_id, go_id=go_id, evidence=evidence, category=category, \
									go_description=go_term, pubmed_ids=evidence_source)
					gene2go.gene = gene
					session.add(gene2go)
					session.flush()
					param_obj.no_of_into_db += 1
				else:
					gene2go.pubmed_ids += '; %s'%evidence_source
					param_obj.no_of_gene2go_already_in_db += 1
		
		if gene.gene_id:
			entrezgene_mapping = EntrezgeneMapping.query.filter_by(gene_id=gene.gene_id).first()
		else:
			entrezgene_mapping = None
		if not entrezgene_mapping:
			entrezgene_mapping = EntrezgeneMapping(tax_id=annot_assembly.tax_id, chromosome=annot_assembly.chromosome,\
											start=start, stop=stop, strand=strand,)
			entrezgene_mapping.genomic_annot_assembly = annot_assembly
			entrezgene_mapping.gene = gene
			entrezgene_mapping.entrezgene_type = entrezgene_type
			
			session.add(entrezgene_mapping)
			param_obj.no_of_into_db += 1
		else:
			param_obj.no_of_entrezgene_mappings_already_in_db += 1
		
		for model_elem in element.findall("MODEL"):
			self.addGeneCommentary(model_elem, entrezgene_mapping, session, param_obj=param_obj)
		
		param_obj.no_of_total += 1
		if param_obj.no_of_total%1000==0:
			try:
				session.flush()
				#session.expunge_all()
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
				import pdb
				pdb.set_trace()
			
			report = getattr(param_obj, 'report', False)
			if report:
				sys.stderr.write('%s%s into db out of %s. %s gene(s) already_in_db.\n\
						%s entrezgene_mapping(s) already_in_db. %s gene2go(s) already in db.\n'%\
						('\x08'*100, param_obj.no_of_into_db, param_obj.no_of_total, param_obj.no_of_genes_already_in_db,\
						param_obj.no_of_entrezgene_mappings_already_in_db, param_obj.no_of_gene2go_already_in_db))
		
	def parse_xml_file(self, session, filename, tax_id=None, param_obj=None):
		"""
		2010-8-11
			
		"""
		counter = 0
		real_counter = 0
		if filename[-2:]=='gz':	#2010-8-12
			import gzip
			fHandler = gzip.open(filename)
		else:
			fHandler = open(filename)
		
		for event, elem in ElementTree.iterparse(fHandler):
			if elem.tag == 'ASSEMBLY':
				assembly_elem = elem
				chromosome = assembly_elem.get("CHROMOSOME")
				current_date = assembly_elem.get("CURRENT_DATE")
				
				start = int(assembly_elem.findtext('COORDSET/END5'))
				stop = int(assembly_elem.findtext('COORDSET/END3'))
				organism = assembly_elem.findtext("HEADER/ORGANISM")
				lineage = assembly_elem.findtext("HEADER/LINEAGE")
				date_last_touched = assembly_elem.findtext("HEADER/SEQ_LAST_TOUCHED/DATE")
				clone_name = assembly_elem.findtext("HEADER/CLONE_NAME")
				comment = 'CURRENT_DATE %s, SEQ_LAST_TOUCHED %s, lineage %s, %s'%(current_date, date_last_touched, lineage, organism)
				
				annot_assembly = AnnotAssembly.query.filter_by(acc_ver=clone_name).first()
				if not annot_assembly:
					annot_assembly = AnnotAssembly(acc_ver=clone_name, chromosome=chromosome, start=start, stop=stop, tax_id=tax_id, \
												sequence_type_id=1, comment=comment)
					session.add(annot_assembly)
					session.flush()
				
				protein_elements = elem.find('GENE_LIST/PROTEIN_CODING')
				RNA_elements = elem.find('GENE_LIST/RNA_GENES')
				TE_elements = elem.find('GENE_LIST/TRANSPOSABLE_ELEMENT_GENES')
				
				for protein_element in protein_elements:
					self.addGene(session, annot_assembly, protein_element, type_of_gene='protein-coding',\
								param_obj=param_obj)
					real_counter += 1
				for element in RNA_elements:
					self.addGene(session, annot_assembly, element, param_obj=param_obj)
					real_counter += 1
				for element in TE_elements:
					self.addGene(session, annot_assembly, element, param_obj=param_obj)
					real_counter += 1
				
				#repeat_elements = elem.find("REPEAT_LIST")
				#misc_elements = elem.find("MISC_INFO")
				
				#release memory
				elem.clear()
				counter += 1
			if self.report and counter%2000==0 and counter>0:
				sys.stderr.write("%s\t%s/%s"%('\x08'*20, counter, real_counter))
		sys.stderr.write('%s%s into db out of %s. %s gene(s) already_in_db.\n \
					%s entrezgene_mapping(s) already_in_db. %s gene2go(s) already in db.\n'%\
					('\x08'*100, param_obj.no_of_into_db, param_obj.no_of_total, param_obj.no_of_genes_already_in_db,\
					param_obj.no_of_entrezgene_mappings_already_in_db, param_obj.no_of_gene2go_already_in_db))
	
	def run(self):
		"""
		11-13-05 
			--db_connect()
			--parse_entrezgene_xml_file()
				--is_gi_valid_in_annot_assembly_table()
				--find_info_dict()
					--return_location_list()
				--submit_to_entrezgene_mapping_table()
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(self.inputfiles))
		db = GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)	#2010-6-22
		session = db.session
		param_obj = PassingData(session=db.session, no_of_genes_already_in_db=0, no_of_entrezgene_mappings_already_in_db=0,\
					no_of_total=0, no_of_into_db=0, report=self.report, no_of_commentaries_already_in_db=0,\
					no_of_gene_segments_already_in_db=0, no_of_gene2go_already_in_db=0)
		for f in self.inputfiles:
			sys.stderr.write("%d/%d:\t%s\n"%(self.inputfiles.index(f)+1,len(self.inputfiles),f))
			self.parse_xml_file(session, f, tax_id=self.tax_id, param_obj=param_obj)
		
		session.flush()
		if self.commit:
			session.commit()
		else:
			session.rollback()

if __name__ == '__main__':
	#import pdb
	#pdb.set_trace()
	#from pymodule import ProcessOptions
	main_class = TAIRGeneXML2GenomeDB
	instance = TAIRGeneXML2GenomeDB(sys.argv)
	#po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	#instance = main_class(**po.long_option2value)
	instance.run()
