#!/usr/bin/env python
"""

Examples:
	FindCNVContext.py -u yh -c

Description:
	program to find the context (nearby genes) of a CNV. It fills results into db table cnv_context.
	
	Suck all the TU(transcription unit)s into a RBDict allowing any degree of overlap between boxes,
		even entirely identical (assuming the UTR/exon/intron structures of two identical TUs will be different).
		For each TU, construct a RBDict with its UTR/exon/intron structure. This RBDict serves as the value in the RBDict.
	For any given CNV, find all elements in RBDict that have any overlap with the CNV.
		for each element, find its UTR/exon/intron RBDict, and find out whether CNV overlaps with which exon, which intron, which UTR.


"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, getopt
import warnings, traceback
from pymodule import ProcessOptions, PassingData, GenomeDB
import Stock_250kDB
from pymodule.db import formReadmeObj
from pymodule.CNV import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
from pymodule.RBTree import RBDict

class FindCNVContext(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'stock_250k database name', ],\
							('genome_dbname', 1, ): ['genome_tair', 'g', 1, 'genome database name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('cnv_method_id', 1, int): [None, 'm', 1, 'construct contexts for CNVs from this cnv_method_id'],\
							('output_fname', 0, ): [None, 'o', 1, 'if given, QC results will be outputed into it.'],\
							('max_distance', 0, int): [20000, 'x', 1, "maximum distance allowed between a CNV and a gene"],\
							('tax_id', 0, int): [3702, '', 1, 'Taxonomy ID to get gene position and coordinates.'],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-08-19
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def createGenomeRBDict(self, genome_db, tax_id=3702, max_distance=20000, debug=False):
		"""
		2010-8-17
		"""
		sys.stderr.write("Creating a RBDict for all genes from organism %s ... \n"%tax_id)
		genomeRBDict = RBDict()
		query = GenomeDB.EntrezgeneMapping.query.filter_by(tax_id=tax_id)
		counter = 0
		real_counter = 0
		for row in query:
			try:	# convert to integer except when "C" or "M"/mitochondria is encountered.
				chromosome = int(row.chromosome)	#integer chromosomes should be converted as CNV.chromosome is integer.
			except:
				chromosome = row.chromosome
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=chromosome, \
							span_ls=[max(1, row.start - max_distance), row.stop + max_distance], \
							min_reciprocal_overlap=1, strand=row.strand, gene_id=row.gene_id,\
							gene_start=row.start, gene_stop=row.stop)	#2010-8-17 any overlap is tolerated.
			genomeRBDict[segmentKey] = []
			counter += 1
			for gene_commentary in row.gene_commentaries:
				if not gene_commentary.gene_commentary_id:
					# ignore gene_commentary that are derived from other gene_commentaries. 
					# they'll be handled within the parental gene_commentary.
					geneRBDict = RBDict()
					geneRBDict.gene_commentary_id = gene_commentary.id
					#gene_commentary.construct_annotated_box()
					box_ls = gene_commentary.constructAnnotatedBox()
					#box_ls=gene_commentary.box_ls
					no_of_boxes = len(box_ls)
					
					numberPorter = PassingData(cds_number = 0,\
											intron_number = 0,\
											utr_number = 0,\
											exon_number = 0)
					for i in xrange(no_of_boxes):
						if row.strand == "-1":	#reverse
							box = box_ls[-i-1]
						else:
							box = box_ls[i]
						start, stop, box_type, is_translated, gene_segment_id = box[:5]
						numberVariableName = None
						if box_type=='3UTR' or box_type=='5UTR':
							numberPorter.utr_number += 1
							numberVariableName = 'utr_number'
						elif box_type=='CDS':
							numberPorter.cds_number += 1
							numberVariableName = 'cds_number'
						elif box_type=='intron':
							numberPorter.intron_number += 1
							numberVariableName = 'intron_number'
						elif box_type=='exon':
							numberPorter.exon_number += 1
							numberVariableName = 'exon_number'
						genePartKey = CNVSegmentBinarySearchTreeKey(chromosome=chromosome, span_ls=[start, stop], \
														min_reciprocal_overlap=1, label=box_type, cds_number=None,\
														intron_number=None, utr_number=None,\
														exon_number=None, gene_segment_id=gene_segment_id)
									#2010-8-17 any overlap is tolerated.
						if numberVariableName is not None:	#set the specific number
							setattr(genePartKey, numberVariableName, getattr(numberPorter, numberVariableName, None))
						geneRBDict[genePartKey] = None
						real_counter += 1
					genomeRBDict[segmentKey].append(geneRBDict)
			if counter%1000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*100, counter, real_counter))
				if debug:
					break
		sys.stderr.write("%s%s\t%s\n"%('\x08'*100, counter, real_counter))
		sys.stderr.write("%s Done.\n"%(str(genomeRBDict)))
		return genomeRBDict
	
	def findCNVcontext(self, db_250k, genomeRBDict, cnv_method_id=None, compareIns=None, max_distance=50000, debug=0,
					param_obj=None):
		"""
		2010-8-18
		"""
		sys.stderr.write("Finding CNV context ... \n")
		session = db_250k.session
		TableClass = Stock_250kDB.CNV
		query = TableClass.query.filter_by(cnv_method_id=cnv_method_id)
		for row in query:
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=row.chromosome, \
							span_ls=[row.start, row.stop], \
							min_reciprocal_overlap=0.0000001, )	#min_reciprocal_overlap doesn't matter here.
				# it's decided by compareIns.
			node_ls = []
			genomeRBDict.findNodes(segmentKey, node_ls=node_ls, compareIns=compareIns)
			for node in node_ls:
				geneSegKey = node.key
				geneRBDictLs = node.value
				# geneSegKey.span_ls expands 20kb upstream or downstream of the gene.
				overlap1, overlap2, overlap_length, overlap_start_pos, overlap_stop_pos = get_overlap_ratio(segmentKey.span_ls, \
														[geneSegKey.gene_start, geneSegKey.gene_stop])[:5]
				if overlap_length>0:	#use fraction of length as coordinates.
					gene_length = geneSegKey.gene_stop - geneSegKey.gene_start + 1
					try:
						if geneSegKey.strand == '+1':
							term5_disp_pos = abs(overlap_start_pos - geneSegKey.gene_start)/float(gene_length)
							term3_disp_pos = abs(overlap_stop_pos - geneSegKey.gene_start + 1)/float(gene_length)
						else:
							term5_disp_pos = abs(geneSegKey.gene_stop - overlap_stop_pos)/float(gene_length)
							term3_disp_pos = abs(geneSegKey.gene_stop - overlap_start_pos + 1)/float(gene_length)
					except:
						import pdb
						pdb.set_trace()
				else:	#no overlap at all
					term3_disp_pos = None
					if geneSegKey.strand == '+1':
						if row.stop<=geneSegKey.gene_start:	#upstream
							term5_disp_pos = row.stop - geneSegKey.gene_start
						elif row.start>=geneSegKey.gene_stop:	# downstream
							term5_disp_pos = row.start - geneSegKey.gene_stop
					else:
						if row.stop<=geneSegKey.gene_start:	#downstream
							term5_disp_pos = geneSegKey.gene_start - row.stop
						elif row.start>=geneSegKey.gene_stop:	# upstream
							term5_disp_pos = geneSegKey.gene_stop - row.start
				cnv_context = Stock_250kDB.CNVContext.query.filter_by(cnv_id=row.id).filter_by(gene_id=geneSegKey.gene_id).first()
				
				if cnv_context:
					param_obj.no_of_cnv_contexts_already_in_db += 1
				else:
					cnv_context = Stock_250kDB.CNVContext(cnv_id=row.id, gene_id = geneSegKey.gene_id, \
													gene_strand=geneSegKey.strand, term5_disp_pos=term5_disp_pos, \
													term3_disp_pos=term3_disp_pos,\
													overlap_length=overlap_length, \
													overlap_fraction_in_cnv=overlap1, overlap_fraction_in_gene=overlap2)
					session.add(cnv_context)
					param_obj.no_of_into_db += 1
				param_obj.no_of_total_contexts += 1
				
				for geneRBDict in geneRBDictLs:
					gene_box_node_ls = []
					geneRBDict.findNodes(segmentKey, node_ls=gene_box_node_ls, compareIns=compareIns)
					for gene_box_node in gene_box_node_ls:
						gene_box_key = gene_box_node.key
						overlap1, overlap2, overlap_length, overlap_start_pos, overlap_stop_pos = get_overlap_ratio(segmentKey.span_ls, \
												gene_box_key.span_ls)[:5]
						cnv_annotation = Stock_250kDB.CNVAnnotation.query.filter_by(cnv_id=row.id).filter_by(cnv_context_id=cnv_context.id).\
							filter_by(gene_commentary_id= geneRBDict.gene_commentary_id).\
							filter_by(gene_segment_id= gene_box_key.gene_segment_id).first()
						if cnv_annotation:
							param_obj.no_of_cnv_annotations_already_in_db += 1
						else:
							cnv_annotation = Stock_250kDB.CNVAnnotation(cnv_id=row.id, gene_commentary_id = geneRBDict.gene_commentary_id, \
												gene_segment_id=gene_box_key.gene_segment_id, label=gene_box_key.label, \
												utr_number = gene_box_key.utr_number, cds_number = gene_box_key.cds_number, \
												intron_number = gene_box_key.intron_number, exon_number = gene_box_key.exon_number,\
												overlap_length=overlap_length, \
												overlap_fraction_in_cnv=overlap1, overlap_fraction_in_gene=overlap2)
							cnv_annotation.cnv_context = cnv_context
							session.add(cnv_annotation)
							param_obj.no_of_into_db += 1
						param_obj.no_of_total_annotations += 1
				
				if param_obj.no_of_into_db>2000:
					session.flush()
					param_obj.no_of_into_db = 0
					sys.stderr.write("\t %s/%s CNVContext(s) & %s/%s CNVAnnotation(s) already in db.\n"%(\
								param_obj.no_of_cnv_contexts_already_in_db, param_obj.no_of_total_contexts, \
								param_obj.no_of_cnv_annotations_already_in_db, param_obj.no_of_total_annotations))
		session.flush()
		session.expunge_all()
		sys.stderr.write("\t %s/%s CNVContext(s) & %s/%s CNVAnnotation(s) already in db.\n"%(\
								param_obj.no_of_cnv_contexts_already_in_db, param_obj.no_of_total_contexts, \
								param_obj.no_of_cnv_annotations_already_in_db, param_obj.no_of_total_annotations))
		
	def run(self):
		"""
		2010-8-15
			create a RBDict of gene-forms [chr, start, stop] with min_overlap_ratio=1.
				value is a sub-RBDict of the gene structure (UTR, non-UTR-exon, intron)
			
			given any CNV, use RBDict.findNodes() to find all gene-forms.
				WATCH: use an alternative comparison function.
			
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		genome_db = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.genome_dbname, )
		genome_db.setup(create_tables=False)
		
		genomeRBDict = self.createGenomeRBDict(genome_db, tax_id=self.tax_id, max_distance=self.max_distance, \
											debug=self.debug)
		
		
		db_250k = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user, password=self.db_passwd, \
									hostname=self.hostname, database=self.dbname)
		db_250k.setup(create_tables=False)
		
		param_obj = PassingData(no_of_total_annotations=0, session=db_250k.session, \
					cnv_method_id=self.cnv_method_id, no_of_total_contexts=0, no_of_into_db=0, report=self.report,\
					no_of_cnv_contexts_already_in_db=0, no_of_cnv_annotations_already_in_db=0)
		compareIns = CNVCompare(min_reciprocal_overlap=0.0000001)	#any overlap is an overlap
		self.findCNVcontext(db_250k, genomeRBDict, cnv_method_id=self.cnv_method_id, compareIns=compareIns, \
					max_distance=self.max_distance, debug=self.debug, param_obj=param_obj)
		db_250k.session.flush()
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = FindCNVContext
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
