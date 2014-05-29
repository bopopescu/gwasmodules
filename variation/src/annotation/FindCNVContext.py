#!/usr/bin/env python
"""

Examples:
	# Find contexts for cnv method 20 in table CNV
	FindCNVContext.py -z banyan --db_user yh --genome_schema genome -m 20
	
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
from pymodule.db import formReadmeObj
from pymodule import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
from pymodule import RBDict
from variation.src import AbstractVariationMapper
from variation.src import Stock_250kDB

class FindCNVContext(AbstractVariationMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVariationMapper.option_default_dict.copy()
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update(AbstractVariationMapper.genome_db_option_dict.copy())
	
	option_default_dict.update({
							('cnv_method_id', 1, int): [None, 'm', 1, 'construct contexts for CNVs from this cnv_method_id'],\
							('genomeRBDictPickleFname', 1, ): ['', '', 1, 'The file to contain pickled genomeRBDict.'],\
							('max_distance', 0, int): [20000, 'x', 1, "maximum distance allowed between a CNV and a gene"],\
							('tax_id', 0, int): [3702, '', 1, 'Taxonomy ID to get gene position and coordinates.'],\
							})
	
	def __init__(self,  **keywords):
		"""
		2008-08-19
		"""
		AbstractVariationMapper.__init__(self, inputFnameLs=None, **keywords)
	
	
	def findCNVcontext(self, db_250k, genomeRBDict, cnv_method_id=None, compareIns=None, max_distance=50000, debug=0,
					param_obj=None):
		"""
		2011-3-25
			cast row.chromosome (from db) into str type.
		2010-10-3
			bug fixed: (chr, start, stop) is not unique. There are genes with the same coordinates.
		2010-8-18
		"""
		sys.stderr.write("Finding CNV context ... \n")
		session = db_250k.session
		TableClass = Stock_250kDB.CNV
		query = TableClass.query.filter_by(cnv_method_id=cnv_method_id)
		for row in query:
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=str(row.chromosome), \
							span_ls=[row.start, row.stop], \
							min_reciprocal_overlap=0.0000001, )	#min_reciprocal_overlap doesn't matter here.
				# it's decided by compareIns.
			node_ls = []
			genomeRBDict.findNodes(segmentKey, node_ls=node_ls, compareIns=compareIns)
			for node in node_ls:
				geneSegKey = node.key
				for oneGeneData in node.value:
					# geneSegKey.span_ls expands 20kb upstream or downstream of the gene.
					overlapData = get_overlap_ratio(segmentKey.span_ls, \
															[oneGeneData.gene_start, oneGeneData.gene_stop])
					overlapFraction1 = overlapData.overlapFraction1
					overlapFraction2 = overlapData.overlapFraction2
					overlap_length = overlapData.overlap_length
					overlap_start_pos = overlapData.overlap_start_pos
					overlap_stop_pos = overlapData.overlap_stop_pos 
					if overlap_length>0:	#use fraction of length as coordinates.
						gene_length = oneGeneData.gene_stop - oneGeneData.gene_start + 1
						try:
							if oneGeneData.strand == '+1':
								term5_disp_pos = abs(overlap_start_pos - oneGeneData.gene_start)/float(gene_length)
								term3_disp_pos = abs(overlap_stop_pos - oneGeneData.gene_start + 1)/float(gene_length)
							else:
								term5_disp_pos = abs(oneGeneData.gene_stop - overlap_stop_pos)/float(gene_length)
								term3_disp_pos = abs(oneGeneData.gene_stop - overlap_start_pos + 1)/float(gene_length)
						except:
							import pdb
							pdb.set_trace()
					else:	#no overlap at all
						term3_disp_pos = None
						if oneGeneData.strand == '+1':
							if row.stop<=oneGeneData.gene_start:	#upstream
								term5_disp_pos = row.stop - oneGeneData.gene_start
							elif row.start>=oneGeneData.gene_stop:	# downstream
								term5_disp_pos = row.start - oneGeneData.gene_stop
						else:
							if row.stop<=oneGeneData.gene_start:	#downstream
								term5_disp_pos = oneGeneData.gene_start - row.stop
							elif row.start>=oneGeneData.gene_stop:	# upstream
								term5_disp_pos = oneGeneData.gene_stop - row.start
					cnv_context = Stock_250kDB.CNVContext.query.filter_by(cnv_id=row.id).filter_by(gene_id=oneGeneData.gene_id).first()
					
					if cnv_context:
						param_obj.no_of_cnv_contexts_already_in_db += 1
					else:
						cnv_context = Stock_250kDB.CNVContext(cnv_id=row.id, gene_id = oneGeneData.gene_id, \
														gene_strand=oneGeneData.strand, term5_disp_pos=term5_disp_pos, \
														term3_disp_pos=term3_disp_pos,\
														overlap_length=overlap_length, \
														overlap_fraction_in_cnv=overlapFraction1, overlap_fraction_in_gene=overlapFraction2)
						session.add(cnv_context)
						param_obj.no_of_into_db += 1
					param_obj.no_of_total_contexts += 1
					
					for geneCommentaryRBDict in oneGeneData.geneCommentaryRBDictLs:
						gene_box_node_ls = []
						geneCommentaryRBDict.findNodes(segmentKey, node_ls=gene_box_node_ls, compareIns=compareIns)
						for gene_box_node in gene_box_node_ls:
							gene_box_key = gene_box_node.key
							overlapData = get_overlap_ratio(segmentKey.span_ls, gene_box_key.span_ls)
							overlapFraction1 = overlapData.overlapFraction1
							overlapFraction2 = overlapData.overlapFraction2
							overlap_length = overlapData.overlap_length
							overlap_start_pos = overlapData.overlap_start_pos
							overlap_stop_pos = overlapData.overlap_stop_pos
							cnv_annotation = Stock_250kDB.CNVAnnotation.query.filter_by(cnv_id=row.id).filter_by(cnv_context_id=cnv_context.id).\
								filter_by(gene_commentary_id= geneCommentaryRBDict.gene_commentary_id).\
								filter_by(gene_segment_id= gene_box_key.gene_segment_id).first()
							if cnv_annotation:
								param_obj.no_of_cnv_annotations_already_in_db += 1
							else:
								cnv_annotation = Stock_250kDB.CNVAnnotation(cnv_id=row.id, \
												gene_commentary_id = geneCommentaryRBDict.gene_commentary_id, \
												gene_segment_id=gene_box_key.gene_segment_id, label=gene_box_key.label, \
												utr_number = gene_box_key.utr_number, cds_number = gene_box_key.cds_number, \
												intron_number = gene_box_key.intron_number, exon_number = gene_box_key.exon_number,\
												overlap_length=overlap_length, \
												overlap_fraction_in_cnv=overlapFraction1, overlap_fraction_in_gene=overlapFraction2)
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
		
		genome_db = GenomeDB.GenomeDatabase(drivername=self.genome_drivername, db_user=self.genome_db_user,
						db_passwd=self.genome_db_passwd, hostname=self.genome_hostname, dbname=self.genome_dbname, \
						schema=self.genome_schema)
		genome_db.setup(create_tables=False)
		
		"""
		genome_db = GenomeDB.GenomeDatabase(drivername=self.drivername, db_user=self.db_user,
						db_passwd=self.db_passwd, hostname=self.hostname,\
						database=self.genome_dbname, schema=self.genome_schema)
		genome_db.setup(create_tables=False)
		"""
		genomeRBDict = genome_db.dealWithGenomeRBDict(self.genomeRBDictPickleFname, tax_id=self.tax_id, \
													max_distance=self.max_distance, debug=self.debug)
		#genome_db.createGenomeRBDict(tax_id=self.tax_id, max_distance=self.max_distance, \
		#									debug=self.debug)
		#2011-3-10 temporary: exit early , only for the genomeRBDictPickleFname
		#sys.exit(3)
		
		db_250k = self.db_250k
		
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
