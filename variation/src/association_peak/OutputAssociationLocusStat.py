#!/usr/bin/env python
"""
Examples:
	%s  -i 1-35,39-48,57-82,158-159,161-179,182-186,272-283,314-351,362-380,418-589
		-l 80,75 -u yh -z banyan -s 2,3,1,3 -c
	
	%s  -i ... -l 32,80,57,75 -u yh -z banyan -s 2,3,1,3 -c

Description:
	2012.9.28
		program that outputs statistics of association result peaks:
			associationPeakID-chr-start-stop, cumulativeStart,cumulativeStop, connectivity,
			recurrence-in-CNV, recurrence-in-SNP,  median p-value in CNV, median p-value in SNP,
			recurrence-in-CNV-EMMA, recurrence-in-SNP-EMMA,  median p-value in CNV-EMMA, median p-value in SNP-EMMA
	
	Restrict result peaks to to phenotype_id_ls, analysis 
	1. get all association peaks => result_peak_ls => result (filter: phenotype, analysis, call method)
		
	
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
import numpy
import networkx as nx
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
#from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
from pymodule import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
from pymodule import RBDict
from pymodule import yh_matplotlib, GenomeDB, utils
from variation.src.mapper.AbstractVariationMapper import AbstractVariationMapper
from variation.src import Stock_250kDB

from AssociationPeak2AssociationLocus import AssociationPeak2AssociationLocus

class OutputAssociationLocusStat(AssociationPeak2AssociationLocus):
	__doc__ = __doc__
	option_default_dict = AssociationPeak2AssociationLocus.option_default_dict.copy()
	option_default_dict.update(AssociationPeak2AssociationLocus.db_option_dict)
	option_default_dict.update({
							('call_method_id_ls', 1, ): ['', 'l', 1, 'loci from these peaks will be outputted', ],\
							('analysis_method_id_ls', 1, ): ['1,32', 'a', 1, 'loci from these peaks will be outputted', ],\
							('phenotype_id_ls', 1, ): ['', 'i', 1, 'loci from these peaks will be outputted', ],\
							('min_overlap_ratio', 1, float): ['0.05', 'm', 1, 'minimum overlap ratio, overlap length/total' ],\
							('result_peak_type_id_ls', 1, ): ['', 's', 1, 'list of result peak type IDs, in this order:\
								snp-KW-gwas,snp-EMMA-gwas,deletion-KW-gwas,deletion-EMMA-gwas', ],\
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		AssociationPeak2AssociationLocus.__init__(self, inputFnameLs, **keywords)
		
		self.call_method_id_ls = getListOutOfStr(self.call_method_id_ls, data_type=int)
		self.analysis_method_id_ls = getListOutOfStr(self.analysis_method_id_ls, data_type=int)
		self.phenotype_id_ls = getListOutOfStr(self.phenotype_id_ls, data_type=int)
		self.phenotype_id_ls.sort()
		self.result_peak_type_id_ls = getListOutOfStr(self.result_peak_type_id_ls, data_type=int)
		
		self.phenotype_id_set = set(self.phenotype_id_ls)
		self.call_method_id_set = set(self.call_method_id_ls)
		self.analysis_method_id_set = set(self.analysis_method_id_ls)
	
	connectDB = AbstractVariationMapper.connectDB
	def run(self):
		"""
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		#without commenting out db_vervet connection code. schema "genome" wont' be default path.
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		#chrOrder=1 is to order chromosomes alphabetically
		oneGenomeData = db_genome.getOneGenomeData(tax_id=3702, chr_gap=0, chrOrder=1, sequence_type_id=1)
		chr_id2cumu_start = oneGenomeData.chr_id2cumu_start
		
		db_250k = self.db_250k
		
		db_250k.session.begin()
		
		writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		header=['associationLocusID', 'chromosome|string', 'cumulativeStart', 'cumulativeStop', 'connectivity']
		
		i = 0
		block_size = 5000
		TableClass = Stock_250kDB.AssociationLocus
		rows = TableClass.query.offset(i).limit(block_size)
		locus_type_id2index = {1:0, 2:1}	#to order the data in the output, 1 is SNP. 2 is CNV 
		self.analysis_method_id_ls.sort()
		analysis_method_id2index = {}	#to order the data in the output
		for i in xrange(len(self.analysis_method_id_ls)):
			analysis_method_id = self.analysis_method_id_ls[i]
			analysis_method_id2index[analysis_method_id] = len(analysis_method_id2index)
		
		locus_type_id_analysis_method_id2orderIndex = {}
		dataStartColIndex = len(header)
		for locus_type_id in [1,2]:
			for analysis_method_id in self.analysis_method_id_ls:
				header.append('frequency-locusType%s-%s'%(locus_type_id, analysis_method_id))
				header.append('medianPvalue-locusType%s-%s'%(locus_type_id, analysis_method_id))
				key = (locus_type_id, analysis_method_id)
				locus_type_id_analysis_method_id2orderIndex[key] = len(locus_type_id_analysis_method_id2orderIndex)
				
		writer.writerow(header)
		
		noOfLocusType = 2
		while rows.count()!=0:
			for row in rows:
				i += 1
				locusTypeID2analysisMethod2Data = {}
				for result_peak in row.result_peak_ls:
					call_method_id = result_peak.result.call_method_id
					phenotype_method_id = result_peak.result.phenotype_method_id
					analysis_method_id = result_peak.result.analysis_method_id
					locus_type_id = result_peak.result.call_method.locus_type_id
					if call_method_id in self.call_method_id_set and phenotype_method_id in self.phenotype_id_set \
						and analysis_method_id in self.analysis_method_id_set:
						
						if locus_type_id not in locusTypeID2analysisMethod2Data:
							locusTypeID2analysisMethod2Data[locus_type_id] = {}
						
						
						analysisMethod2Data = locusTypeID2analysisMethod2Data[locus_type_id]
						if analysis_method_id not in analysisMethod2Data:
							analysisMethod2Data[analysis_method_id] = PassingData(peak_score_ls =[], phenotype_id_set=set())
						analysisMethod2Data[analysis_method_id].peak_score_ls.append(result_peak.peak_score)
						analysisMethod2Data[analysis_method_id].phenotype_id_set.add(phenotype_method_id)
				
				if locusTypeID2analysisMethod2Data:	#not empty
					associationLocusID = 'a%s_chr%s_%s_%s'%(row.id, row.chromosome, row.start, row.stop)
					cumuStart = chr_id2cumu_start.get(row.chromosome)
					data_row = [associationLocusID, row.chromosome, cumuStart+row.start, cumuStart+row.stop, row.connectivity]
					for j in xrange(dataStartColIndex, len(header)):	#to fill it up
						data_row.append(0)
					for locus_type_id, analysisMethod2Data in locusTypeID2analysisMethod2Data.iteritems():
						for analysis_method_id, pdata in analysisMethod2Data.iteritems():
							key = (locus_type_id, analysis_method_id)
							orderIndex = locus_type_id_analysis_method_id2orderIndex[key]
							data_row[dataStartColIndex+orderIndex*2] = len(pdata.phenotype_id_set)
							data_row[dataStartColIndex+orderIndex*2+1] = numpy.median(pdata.peak_score_ls)
					writer.writerow(data_row)
			rows = TableClass.query.offset(i).limit(block_size)
		del writer


if __name__ == '__main__':
	main_class = OutputAssociationLocusStat
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()		