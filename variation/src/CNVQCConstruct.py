#!/usr/bin/env python
"""

Examples:
	# Clark et al. Science 2007's data
	CNVQCConstruct.py -i ~/script/variation/data/CNV/PERL_deletions_in_PRPs.txt -s Clark2007a -t deletion -m PCR -u yh -n 1 -c
	
	# Korbinian Schneeberger & Stephan Ossowski Paired-end
	CNVQCConstruct.py -i ~/script/variation/data/CNV/PERL_deletions_in_PRPs.txt -s SchneebergerOssowski2009
		-t deletion -m PairedEndSolexa -u yh -y 2 -c
	
	# Bob schmitz data
	CNVQCConstruct.py -i ~/script/variation/data/CNV/cvi_chr1_indels -s BobSchmitz2009 -t DEL -m PairedEndSolexa -o Cvi-0
		-u yh -y 3 -c
	
	# Ler-contig derived deletions
	CNVQCConstruct.py -i /tmp/Ler-deletions.tsv -s LerContigCoveredDeltion -t DEL -m PCR -o Ler-1 -u yh -y 4 -c -n 1
	
	# Ler-contig derived deletions/duplications. Deletions could also be uncovered sequences, rather than real deletions.
	# "-t DEL" is useless.
	CNVQCConstruct.py -i /tmp/Ler-copy-number.tsv -s LerContigCNV -t DEL -m PCR -o Ler-1 -u yh -y 4 -c -n1
	
	# 2010-6-20 Ler-contig derived normal segments by running nucmer and CNV.LerContig.discoverLerDeletionAndSpanFromMummerCoordsOutput()
	# "-t DEL" is useless.
	CNVQCConstruct.py -i data/CNV/mummer/TAIR8_LerContig_maxmatch_maxgap40_mincluster60_min_identity_perc0.9...ref_covered.tsv
		-s LerContigNormalMaxMatch -t DEL -m PCR -o Ler-1 -u yh -y 4 -c -z banyan
	
	CNVQCConstruct.py -i data/CNV/mummer/TAIR8_LerContig_mum_maxgap40_mincluster60_min_identity_perc0.9_..._ref_covered.tsv
		-s LerContigNormalMum -t DEL -m PCR -o Ler-1 -u yh -y 4 -c -z banyan
	
	
	# Ler-contig span coverage over Col ref genome (the corresponding Col start-stop positions for a Ler contig)
	# it's discovered by running CNV.LerContig.discoverLerContigSpanOverCol() from variation/src/misc.py
	# CNV type (-t) is useless.
	CNVQCConstruct.py -i ~/script/variation/data/CNV/Ler-span-over-Col-mdr0.4-mld50000.tsv -s LerContigCoverage
		-t unknown -m PCR -o Ler-1 -u yh -y 5 -c -e 1
	
	# put MpiBlast.py output into table SequenceFragment2Probe
	CNVQCConstruct.py -i ~/script/variation/data/CNV/ler_blast_2.tsv -s LerContigCoverage
		-t unknown -m PCR -o Ler-1 -u yh -y 6 -c -n1
	
	# 2010-6-20 Ler-contig span coverage over Col ref genome (the corresponding Col start-stop positions for a Ler contig)
	# it's discovered by nucmer + show-coords of the mummer package. First 3 lines of input file will be skipped.
	# CNV type (-t) is useless.
	CNVQCConstruct.py -i ~/script/variation/data/CNV/mummer/TAIR8_LerContig_maxmatch_maxgap40_mincluster60.coords
		-s LerContigCoverage
		-t unknown -m PCR -o Ler-1 -u yh -y 7 -c -e 2 -z banyan
	
	# 2010-6-25 BreakDancerMax's output from Quan Long's PE-derived deletions
	~/script/variation/src/CNVQCConstruct.py -i algustrum.8230.bkd -s QuanLong75bpPairEnd -t unknown -m PairedEndSolexaAndBreakDancer 
	-o algustrum.8230 -u yh -y 8 -c -z banyan
	
	# 2010-6-25 data derived (by CNV.QuanLongPairedEndSolexaData.discoverDeletionsFromCoverageData() ) 
	# from Quan Long's PE-derived coverage data
	# format: original_id, ecotype_id, chromosome, start, stop, size_affected, copy_number, mean_coverage
	~/script/variation/src/CNVQCConstruct.py -i coverage4Yu_deletions.tsv -s QuanLong75bpPairEnd -t DEL -m PairedEndSolexaCoverage 
		-o randomID -u yh -y 9 -c -z banyan
	
	#2010-7-31 Lyrata-sequence-derived normal segments by running nucmer and 
	#	CNV.LerContig.discoverLerDeletionAndSpanFromMummerCoordsOutput()
	#	"-t DEL" is useless.
	CNVQCConstruct.py -i data/CNV/mummer/TAIR9_Lyrata_maxmatch_maxgap40_mincluster60_min_identity_perc0.9_ref_covered.tsv
		-s TAIR9LyrataNormalMaxMatch -t DEL -m PCR -o Lyrata -u yh -y 4 -c -z banyan
	
	#2010-8-2 put output of CNV.Lyrata.getLyrataNormalANDDeletionFromQuanFileInTAIR9() into db
	# It's basically the portion of A. thaliana that could be matched to A lyrata.
	CNVQCConstruct.py -i ~/script/variation/data/lyrata/at_ancestor_t9.normaSegment.tsv
		-s TAIR9LyrataNormalByTinaAndQuan -t DEL -m PCR -o Lyrata -u yh -y 4 -c -z banyan
	
Description:
	2009-10-28 Fill up relevant db tables with CNV QC data from different sources.
	
	Argument -t "cnv_type" is useless for run_type 3 & 4. it depends on whether size_affected is positive or negative or copy_number=0 or >1.
	
	run_type 4: Ler-contig (or other raw sequence) derived CNVs. First blast CNV probes against contig. 
		Then Run CNV.LerContig.discoverLerDeletionDuplication() 
			or CNV.LerContig.discoverLerDeletionAndSpanFromMummerCoordsOutput() from variation/src/misc.py
		No need to run CNV.fillDeletionInCNVQCProbeCalls() to set no_of_probes_covered. It's in the input file.
	
	run_type 7: The data is sent to table SequenceFragmentRefPos, same as run_type 5. Differentiated through column version.
		First 4 lines of input are skipped. The 4th line is header.
		no_of_lines_to_skip is set to 4 internally regardless of the argument given by the user.
	
	run_type 8:
		First 2 lines of input are skipped. The 2nd line is header.
		no_of_lines_to_skip is set to 2 internally regardless of the argument given by the user.
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import csv, numpy, traceback, subprocess
from pymodule import figureOutDelimiter, getListOutOfStr, getColName2IndexFromHeader
import Stock_250kDB

class CNVQCConstruct(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('input_fname', 1, ): ['', 'i', 1, 'input file.', ],\
							('data_source', 1, ): ['', 's', 1, 'a short name for data source'],\
							('cnv_type', 1, ): ['', 't', 1, 'which type of CNV? INS, DEL, etc. in table cnv_type or will be inserted if new.'],\
							('cnv_method', 0, ): ['', 'm', 1, 'method used to derive CNV. If not provided, it is ignored.', ],\
							('no_of_lines_to_skip', 1, int): [1, 'n', 1, 'Number of lines to skip in the beginning of the input file. Header is the first of such lines.', ],\
							('original_id', 0, ): ['', 'o', 1, 'only used for run_type 3 & 4. cuz data is from one accession and not specifying which one.', ],\
							('version', 1, int): [1, 'e', 1, 'corresponds to SequenceFragmentRefPos.version. only for run_type 5 & 7'],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.'],\
							('run_type', 0, int):[1, 'y', 1, 'Run type 1: Clark et al. 2007, 2: Korbinian Schneeberger & Stephan Ossowski Paired-end, \
								3: Bob Schmitz, \
								4: Ler-contig derived CNVs, 5: Ler-contig span over Col genome (table SequenceFragmentRefPos),\
								6: MpiBlast.py output into SequenceFragment2Probe, \
								7: Ler-contig span info derived from nucmer coords output (table SequenceFragmentRefPos),\
								8: Quan Long BreakDancer output (table CNVQCCall),\
								9: data derived from Quan Long coverage data (table CNVQCCall)']}

	def __init__(self, **keywords):
		"""
		2009-10-28
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
	
	cnv_type2cnv_type_obj = {}
	chr_name_tax_id2chr_obj = {}
	original_id_data_source2acc_obj = {}

	@classmethod
	def getDBObj(cls, session, table_class, short_name):
		"""
		2009-10-28
		"""
		sys.stderr.write("Getting/Creating a %s object ..."%(table_class.table.name))
		db_obj = table_class.query.filter_by(short_name=short_name).first()
		if not db_obj:
			db_obj = table_class(short_name=short_name)
			session.add(db_obj)
			session.flush()
		sys.stderr.write("Done.\n")
		return db_obj
	
	def getCNVTypeObj(self, session, short_name):
		"""
		2010-6-24
			a short cut
		"""
		if short_name in self.cnv_type2cnv_type_obj:
			cnv_type_obj = self.cnv_type2cnv_type_obj.get(short_name)
		else:
			cnv_type_obj = self.getDBObj(session, Stock_250kDB.CNVType, short_name)
			self.cnv_type2cnv_type_obj[short_name] = cnv_type_obj
		return self.cnv_type2cnv_type_obj[short_name]
	
	@classmethod
	def getChromosomeObjFromDB(cls, session, name, tax_id=3702):
		"""
		2010-6-24
		"""
		table_class = Stock_250kDB.Chromosome
		sys.stderr.write("Getting/Creating a %s object ..."%(table_class.table.name))
		db_obj = table_class.query.filter_by(name=name).filter_by(tax_id=tax_id).first()
		if not db_obj:
			db_obj = table_class(name=name, tax_id=tax_id)
			session.add(db_obj)
			session.flush()
		sys.stderr.write("Done.\n")
		return db_obj
	
	def getChromosomeObj(self, session, name, tax_id=3702):
		"""
		2010-6-24
		"""
		chr_name_tax_id = (name, tax_id)
		if chr_name_tax_id in self.chr_name_tax_id2chr_obj:
			chr_obj = self.chr_name_tax_id2chr_obj.get(chr_name_tax_id)
		else:
			chr_obj = self.getChromosomeObjFromDB(session, name, tax_id)
			self.chr_name_tax_id2chr_obj[chr_name_tax_id] = chr_obj
		return self.chr_name_tax_id2chr_obj[chr_name_tax_id]
	
	@classmethod
	def getCNVQCAccessionObjFromDB(cls, session, original_id, data_source_obj, ecotype_id=None):
		"""
		2010-7-22
			if ecotype_id is not None, use it part of query to get CNVQCAccession object.
		2010-6-24
			renamed to getCNVQCAccessionObjFromDB.
			getCNVQCAccessionObj() becomes a function calling this one.
			add optional argument ecotype_id.
		2009-10-28
		"""
		if ecotype_id is not None:
			db_obj = Stock_250kDB.CNVQCAccession.query.filter_by(ecotype_id=ecotype_id).\
				filter_by(data_source_id=data_source_obj.id).first()
		else:
			db_obj = Stock_250kDB.CNVQCAccession.query.filter_by(original_id=original_id).\
				filter_by(data_source_id=data_source_obj.id).first()
		if not db_obj:
			db_obj = Stock_250kDB.CNVQCAccession(original_id=original_id, ecotype_id=ecotype_id)
			db_obj.data_source = data_source_obj
			session.add(db_obj)
			session.flush()
		return db_obj
	
	def getCNVQCAccessionObj(self, session, original_id, data_source_obj, ecotype_id=None):
		"""
		2010-6-24
		"""
		acc_key = (original_id, data_source_obj.short_name)
		if acc_key in self.original_id_data_source2acc_obj:
			acc_obj = self.original_id_data_source2acc_obj.get(acc_key)
		else:
			acc_obj = self.getCNVQCAccessionObjFromDB(session, original_id, data_source_obj, ecotype_id=ecotype_id)
			self.original_id_data_source2acc_obj[acc_key] = acc_obj
		return self.original_id_data_source2acc_obj[acc_key]
	
	@classmethod
	def getSequenceFragmentObj(cls, session, short_name, accession_obj):
		"""
		2010-1-28
		"""
		db_obj = Stock_250kDB.SequenceFragment.query.filter_by(short_name=short_name).filter_by(accession_id=accession_obj.id).first()
		if not db_obj:
			db_obj = Stock_250kDB.SequenceFragment(short_name=short_name)
			db_obj.accession = accession_obj
			session.add(db_obj)
			session.flush()
		return db_obj
	
	def generateCNVQCCallObjFromClark2007(self, session, row, data_source_obj, cnv_type_obj, cnv_method_obj=None, original_id=None):
		"""
		2009-10-29
			run_type 1 (clark et al 2007) format:
				Accession
				Chromosome
				Start deletion position
				End deletion position
				Deletion length (bp)
				Forward primer for validation
				Reverse primer for validation
		"""
		original_id, chromosome, start, stop, size_affected = row[:5]
		score = None
		chromosome = int(chromosome)
		start = int(start)
		stop = int(stop)
		size_affected = int(size_affected)
		
		acc_obj = self.getCNVQCAccessionObj(session, original_id, data_source_obj)
		cnv_qc_call = Stock_250kDB.CNVQCCall(chromosome=chromosome, start=start, stop=stop, \
											size_affected=size_affected, score=score)
		cnv_qc_call.accession = acc_obj
		cnv_qc_call.cnv_type = cnv_type_obj
		cnv_qc_call.cnv_method = cnv_method_obj
		return cnv_qc_call
	
	def generateCNVQCCallObjFromSchneebergerOssowski(self, session, row, data_source_obj, cnv_type_obj, cnv_method_obj=None, original_id=None):
		"""
		2009-10-29
			run_type 2 (Schneeberger & Ossowski) format:
				SampleID
				Number of read pairs supporting deletion
				chromosome
				begin (it is not the beginning of the deletion but the beginning of a region in which the deletion is located)
				end (dito)
				length (as supported by the variation of the distances of the alignments)
				length (as supported by the reference)
				number of positions without core-alignments between begin and end (only looking at non-repetitive)
				number of positions without core-alignments between begin and end
				p-value (likeliness that the given insert size distribution generated the stretched mate pairs)
		"""
		original_id, no_of_read_pairs, chromosome, start, stop, length_supported_by_var_of_dist_of_align, size_affected = row[:7]
		score = float(row[-1])	# Data from Korbinian and Stephan has pvalues.
		chromosome = int(chromosome)
		start = int(start)
		stop = int(stop)
		size_affected = int(size_affected)
		acc_obj = self.getCNVQCAccessionObj(session, original_id, data_source_obj)
		comment = '#read pairs: %s, length by dist var: %s, #positions with no non-repeat-reads: %s, #positions with no reads: %s'%\
				(no_of_read_pairs, length_supported_by_var_of_dist_of_align, row[7], row[8])
		cnv_qc_call = Stock_250kDB.CNVQCCall(chromosome=chromosome, start=start, stop=stop, \
											size_affected=size_affected, score=score, comment=comment)
		cnv_qc_call.accession = acc_obj
		cnv_qc_call.cnv_type = cnv_type_obj
		cnv_qc_call.cnv_method = cnv_method_obj
		return cnv_qc_call
	
	def generateCNVQCCallObjFromBobSchmitzData(self, session, row, data_source_obj, cnv_type_obj, cnv_method_obj=None, original_id=None):
		"""
		2009-10-29
			cnv_type_obj is useless. it depends on whether size_affected is positive or negative.
			
		"""
		chromosome, start_stop, contig_id, contig_start_stop, size_affected = row[:5]
		chromosome = int(chromosome)
		start, stop = start_stop.split(' - ')
		start = int(start)
		stop = int(stop)
		size_affected = int(size_affected)
		if size_affected<0:
			cnv_type = "DEL"
		elif size_affected>0:
			cnv_type = "INS"
		else:
			cnv_type = "NORMAL"
		if cnv_type in self.cnv_type2cnv_type_obj:
			cnv_type_obj = self.cnv_type2cnv_type_obj.get(cnv_type)
		else:
			cnv_type_obj = self.getDBObj(session, Stock_250kDB.CNVType, cnv_type)
			self.cnv_type2cnv_type_obj[cnv_type] = cnv_type_obj
		
		comment = '%s %s'%(contig_id, contig_start_stop)
		acc_obj = self.getCNVQCAccessionObj(session, original_id, data_source_obj)
		cnv_qc_call = Stock_250kDB.CNVQCCall(chromosome=chromosome, start=start, stop=stop, \
											size_affected=abs(size_affected), comment=comment)
		cnv_qc_call.accession = acc_obj
		cnv_qc_call.cnv_type = cnv_type_obj
		cnv_qc_call.cnv_method = cnv_method_obj
		return cnv_qc_call
	
	def generateCNVQCCallObjFromLerContigDerivedCNVs(self, session, row, data_source_obj, cnv_type_obj, cnv_method_obj=None, \
													original_id=None, col_name2index=None):
		"""
		2010-8-2 no probe id information. it's not based on CNV probes.
		2009-12-08
			each row is output of
			
				1. CNV.LerContig.discoverLerDeletionAndSpanFromMummerCoordsOutput() of misc.py
				2. CNV.LerContig.discoverLerDeletionDuplication() from misc.py
			
			cnv_type_obj is useless. it depends on whether copy_number is 0 or >1.
			
		"""
		start_probe_id, start_chr_pos, stop_probe_id, stop_chr_pos, no_of_probes, length, copy_number = row[:7]
		if 'contig' in col_name2index:
			contig_id = row[col_name2index.get('contig')]
		else:
			contig_id = None
		
		if 'size_difference' in col_name2index:
			size_difference = int(row[col_name2index.get('size_difference')])
		else:
			size_difference = None
		
		start_chr_pos = start_chr_pos.split('_')
		start_chr_pos = map(int, start_chr_pos)
		stop_chr_pos = stop_chr_pos.split('_')
		stop_chr_pos = map(int, stop_chr_pos)
		if no_of_probes and no_of_probes!='NA':	#2010-6-17
			no_of_probes = int(no_of_probes)
		else:
			no_of_probes = None
		copy_number = float(copy_number)
		
		chromosome = start_chr_pos[0]
		if start_probe_id =='NA' or not start_probe_id:	#2010-8-2 no probe id information. it's not based on CNV probes.
			start = start_chr_pos[1]
			stop = stop_chr_pos[1]
		else:
			start = start_chr_pos[1] - 12
			stop = stop_chr_pos[1] + 12
		size_affected = size_difference
		if copy_number==0:
			cnv_type = "DEL"
		elif copy_number>1:
			cnv_type = "DUP"
		elif copy_number==1:
			cnv_type = "NORMAl"
		else:
			cnv_type = 'unknown'
		cnv_type_obj = self.getCNVTypeObj(session, cnv_type)
		
		acc_obj = self.getCNVQCAccessionObj(session, original_id, data_source_obj)
		if acc_obj.id and cnv_type_obj.id and cnv_method_obj.id:
			query = Stock_250kDB.CNVQCCall.query.filter_by(accession_id=acc_obj.id).filter_by(chromosome=chromosome).\
				filter_by(start=start).filter_by(stop=stop).filter_by(cnv_type_id=cnv_type_obj.id).\
				filter_by(cnv_method_id=cnv_method_obj.id)
			if query.count()==0:
				new_obj = True
			else:
				new_obj = False
		else:
			new_obj = True
		
		if new_obj:
			cnv_qc_call = Stock_250kDB.CNVQCCall(chromosome=chromosome, start=start, stop=stop, \
											size_affected=abs(size_affected), no_of_probes_covered=no_of_probes,\
											copy_number=copy_number, comment=contig_id)
			cnv_qc_call.accession = acc_obj
			cnv_qc_call.cnv_type = cnv_type_obj
			cnv_qc_call.cnv_method = cnv_method_obj
		else:
			cnv_qc_call = query.first()
			cnv_qc_call.comment += ', %s'%contig_id
		return cnv_qc_call
	
	def generateCNVQCCallObjFromQuanLongBreakDancerOutput(self, session, row, data_source_obj, cnv_type_obj, cnv_method_obj=None,\
													original_id=None, col_name2index=None):
		"""
		2010-6-24
			It's based on 20X pair-end short-sequence data. 75bp per fragment.
			
			Beginning of file looks like:
			
			#vastervik.1.9058.sort.bam	mean:359.380	std:46.520	uppercutoff:544.850	lowercutoff:172.780	readlen:110.000	library:NA	reflen:119665076	...
			#Chr1	Pos1	Orientation1	Chr2	Pos2	Orientation2	Type	Size	Score	num_Reads	num_Reads_lib	Allele_frequency	Version	Run_Param
			Chr1	582	0+10-	Chr1	1942	0+10-	INV	1233	99	10	NA|10	0.24	BreakDancerMax-0.0.1r81	
		"""
		chr = self.getChromosomeObj(session, row[col_name2index['#Chr1']], tax_id=3702)
		start = int(row[col_name2index['Pos1']])
		stopChromosome = self.getChromosomeObj(session, row[col_name2index['Chr2']], tax_id=3702)
		stop = int(row[col_name2index['Pos2']])
		
		if chr.id==stopChromosome.id and start>stop:	# reverse the coordinates
			tmp = stop
			stop = start
			start = tmp
			sys.stderr.write("start, stop reversal corrected. %s.\n"%('\t'.join(row)))
		
		Orientation1 = row[col_name2index.get('Orientation1')]
		Orientation2 = row[col_name2index.get('Orientation2')]
		num_Reads_lib = row[col_name2index.get('num_Reads_lib')]
		Allele_frequency = row[col_name2index.get('Allele_frequency')]
		versionRunProgram = row[col_name2index.get('Version')]

		size_affected = int(row[col_name2index.get('Size')])
		
		score = float(row[col_name2index['Score']])
		num_Reads = int(row[col_name2index['num_Reads']])
		cnv_type_obj = self.getCNVTypeObj(session, row[col_name2index['Type']])
		
		comment = "num_Reads %s, Orientation1 %s, Orientation2 %s, num_Reads_lib %s, Allele_frequency %s, versionRunProgram %s"%\
				(num_Reads, Orientation1, Orientation2, num_Reads_lib, Allele_frequency, versionRunProgram,)
		
		"""
		original_id is the prefix of the filename from Quan Long. It's comprised of ecotype-name and ecotype-id. looks like:
			rev_1.8369 or algustrum.8230
		
		"""
		ecotype_name, ecotype_id = original_id.split('.')[:2]
		ecotype_id = int(ecotype_id)
		acc_obj = self.getCNVQCAccessionObj(session, ecotype_name, data_source_obj, ecotype_id=ecotype_id)
		
		if acc_obj.id and cnv_type_obj.id and cnv_method_obj.id and chr.id and stopChromosome.id:
			query = Stock_250kDB.CNVQCCall.query.filter_by(accession_id=acc_obj.id).filter_by(chromosome=chr.id).\
				filter_by(start=start).filter_by(stop=stop).filter_by(cnv_type_id=cnv_type_obj.id).\
				filter_by(stop_chromosome=stopChromosome.id).filter_by(cnv_method_id=cnv_method_obj.id)
			if query.count()==0:
				new_obj = True
			else:
				new_obj = False
		else:
			new_obj = True
		
		cnv_qc_call = None
		if new_obj:
			cnv_qc_call = Stock_250kDB.CNVQCCall(start=start, stop=stop, \
											size_affected=size_affected, \
											comment=comment, score=score)
			cnv_qc_call.chr = chr
			cnv_qc_call.stopChromosome = stopChromosome
			cnv_qc_call.accession = acc_obj
			cnv_qc_call.cnv_type = cnv_type_obj
			cnv_qc_call.cnv_method = cnv_method_obj
		return cnv_qc_call
	
	def generateCNVQCCallObjFromQuanLongCoverageDerived(self, session, row, data_source_obj, cnv_type_obj=None, \
												cnv_method_obj=None, original_id=None, col_name2index=None):
		"""
		2010-7-22
			the input is coverage data post-processed by CNV.QuanLongPairedEndSolexaData.discoverDeletionOrNonDeletionFromCoverageData()
		"""
		original_id, ecotype_id, chromosome, start, stop, size_affected, copy_number, mean_coverage = row[:8]
		ecotype_id = int(ecotype_id)
		start = int(start)
		stop = int(stop)
		size_affected = int(size_affected)
		copy_number = int(copy_number)
		mean_coverage = float(mean_coverage)
		score = mean_coverage
		
		chr = self.getChromosomeObj(session, chromosome, tax_id=3702)
		
		acc_obj = self.getCNVQCAccessionObj(session, original_id, data_source_obj, ecotype_id=ecotype_id)
		if acc_obj.id and cnv_type_obj.id and cnv_method_obj.id and chr.id:
			query = Stock_250kDB.CNVQCCall.query.filter_by(accession_id=acc_obj.id).filter_by(chromosome=chr.id).\
				filter_by(start=start).filter_by(stop=stop).filter_by(cnv_type_id=cnv_type_obj.id).\
				filter_by(cnv_method_id=cnv_method_obj.id)
			if query.count()==0:
				new_obj = True
			else:
				new_obj = False
		else:
			new_obj = True
		
		cnv_qc_call = None
		if new_obj:
			cnv_qc_call = Stock_250kDB.CNVQCCall(start=start, stop=stop, \
										size_affected=size_affected, score=score, copy_number=copy_number)
			cnv_qc_call.chr = chr
			cnv_qc_call.accession = acc_obj
			cnv_qc_call.cnv_type = cnv_type_obj
			cnv_qc_call.cnv_method = cnv_method_obj
		return cnv_qc_call
	
	def generateSequenceFragmentRefPosObjFromLerContigSpansOverCol(self, session, row, data_source_obj, cnv_type_obj=None, \
													cnv_method_obj=None, \
													original_id=None, col_name2index=None, version=1):
		"""
		2010-6-24
			renamed from generateCNVQCCallObjFromLerContigSpansOverCol() to 
				generateSequenceFragmentRefPosObjFromLerContigSpansOverCol()
		2010-6-14
			add version=1 to SequenceFragmentRefPos
		2010-1-27
			each row is output of CNV.discoverLerContigSpanOverCol() from misc.py.
			
		"""
		start_probe_id, start_chr_pos, stop_probe_id, stop_chr_pos, no_of_probes, length, \
				copy_number, contig_id, contig_start, contig_stop, size_difference = row[:11]
		
		start_probe_id = int(start_probe_id)
		start_chr_pos = start_chr_pos.split('_')
		start_chr_pos = map(int, start_chr_pos)
		stop_probe_id = int(stop_probe_id)
		stop_chr_pos = stop_chr_pos.split('_')
		stop_chr_pos = map(int, stop_chr_pos)
		
		no_of_probes = int(no_of_probes)
		copy_number = copy_number
		contig_start = int(contig_start)
		contig_stop = int(contig_stop)
		size_difference = int(size_difference)
		
		chromosome = start_chr_pos[0]
		start = start_chr_pos[1]
		stop = stop_chr_pos[1]
		
		acc_obj = self.getCNVQCAccessionObj(session, original_id, data_source_obj)
		sequence_fragment_obj = self.getSequenceFragmentObj(session, contig_id, acc_obj)
		
		#2010-1-28 check db whether it's there or not
		# it needs to be checked because sometimes, more than one contig has same Col spanning in terms of probe start&stop
		if sequence_fragment_obj.id:
			query = Stock_250kDB.SequenceFragmentRefPos.query.filter_by(start_probe_id=start_probe_id).\
				filter_by(stop_probe_id=stop_probe_id).\
				filter_by(sequence_fragment_id=sequence_fragment_obj.id).\
				filter_by(fragment_start=contig_start).filter_by(fragment_stop=contig_stop).filter_by(version=version)
				# 2010-6-14 add version
			if query.count()==0:
				new_obj = True
			else:
				new_obj = False
		else:
			new_obj = True
		
		if new_obj:
			sequence_fragment_ref_pos = Stock_250kDB.SequenceFragmentRefPos(chromosome=chromosome, start=start, stop=stop, \
										size_difference=size_difference,\
										start_probe_id=start_probe_id, stop_probe_id=stop_probe_id, \
										no_of_probes_covered=no_of_probes,\
										fragment_start=contig_start, fragment_stop=contig_stop, version=version)
			# 2010-6-14 add version
			sequence_fragment_ref_pos.sequence_fragment = sequence_fragment_obj
		else:
			sequence_fragment_ref_pos = None
		return sequence_fragment_ref_pos

	def generateSequenceFragmentRefPosObjFromNucmerLerContigSpansOverCol(self, session, row, data_source_obj, cnv_type_obj=None, \
													cnv_method_obj=None, \
													original_id=None, col_name2index=None, version=2, comment=''):
		"""
		2010-6-24
			renamed from generateCNVQCCallObjFromNucmerLerContigSpansOverCol() to
				generateSequenceFragmentRefPosObjFromNucmerLerContigSpansOverCol()
		2010-6-14
			similar to generateSequenceFragmentRefPosObjFromLerContigSpansOverCol() but using results from nucmer
		"""
		ref_start_col_index = col_name2index.get('[S1]')
		ref_stop_col_index = col_name2index.get('[E1]')
		query_start_col_index = col_name2index.get('[S2]')
		query_stop_col_index = col_name2index.get('[E2]')
		
		ref_len_col_index = col_name2index.get('[LEN 1]')
		query_len_col_index = col_name2index.get('[LEN 2]')
		
		identity_col_index = col_name2index.get('[% IDY]')
		ref_query_tag_col_index = col_name2index.get('[TAGS]')
		
		start = int(row[ref_start_col_index])
		stop = int(row[ref_stop_col_index])
		contig_start = int(row[query_start_col_index])
		contig_stop = int(row[query_stop_col_index])
		
		ref_len = int(row[ref_len_col_index])
		query_len = int(row[query_len_col_index])
		
		identity_perc = row[identity_col_index]	# it looks like 99.87 = percentage X 100. 
		ref_query_tag_ls = row[ref_query_tag_col_index:]
		chr = ref_query_tag_ls[0]
		chr = int(chr[3])
		contig_id = ref_query_tag_ls[1]
		
		comment = '%s %%IDY %s'%(comment, identity_perc)
		
		size_difference = ref_len-query_len
		
		
		acc_obj = self.getCNVQCAccessionObj(session, original_id, data_source_obj)
		sequence_fragment_obj = self.getSequenceFragmentObj(session, contig_id, acc_obj)
		
		#2010-1-28 check db whether it's there or not
		# it needs to be checked because sometimes, more than one contig has same Col spanning in terms of probe start&stop
		if sequence_fragment_obj.id:
			query = Stock_250kDB.SequenceFragmentRefPos.query.filter_by(chromosome=chr).filter_by(start=start).\
				filter_by(stop=stop).filter_by(sequence_fragment_id=sequence_fragment_obj.id).\
				filter_by(fragment_start=contig_start).filter_by(fragment_stop=contig_stop).filter_by(version=version)
			if query.count()==0:
				new_obj = True
			else:
				new_obj = False
		else:
			new_obj = True
		
		if new_obj:
			sequence_fragment_ref_pos = Stock_250kDB.SequenceFragmentRefPos(chromosome=chr, start=start, stop=stop, \
										size_difference=size_difference,\
										fragment_start=contig_start, fragment_stop=contig_stop, version=version)
			sequence_fragment_ref_pos.sequence_fragment = sequence_fragment_obj
			sequence_fragment_ref_pos.comment = comment
		else:
			sequence_fragment_ref_pos = None
		return sequence_fragment_ref_pos

	def generateSequenceFragment2ProbeObj(self, session, row, data_source_obj, cnv_type_obj=None, \
													cnv_method_obj=None, \
													original_id=None, col_name2index=None):
		"""
		2010-4-15
			each row is output of MpiBlast.py
			
		"""
		contig_id = row[col_name2index['Alignment_title']]
		contig_id = contig_id.split()[1]	# 2010-1-28, take "ATL8C9990" in the case of "ATL8C9990 ATL7C121_1"
		probe_id = int(row[col_name2index['Probe_ID']])
		chr = int(row[col_name2index['Chromosome']])
		pos = int(row[col_name2index['Position']])
		query_start = int(row[col_name2index['query_start']])
		query_end = int(row[col_name2index['query_end']])
		no_of_identities = int(row[col_name2index['Number_matches']])
		fragment_start = int(row[col_name2index['Alignment_start']])
		fragment_stop = int(row[col_name2index['Alignment_stop']])
		
		start = pos-12 + (query_start-1)
		stop = pos + 12 - (25- query_end)
		
		acc_obj = self.getCNVQCAccessionObj(session, original_id, data_source_obj)
		sequence_fragment_obj = self.getSequenceFragmentObj(session, contig_id, acc_obj)
		
		#2010-1-28 check db whether it's there or not
		# it needs to be checked because sometimes, more than one contig has same Col spanning in terms of probe start&stop
		if sequence_fragment_obj.id:
			query = Stock_250kDB.SequenceFragment2Probe.query.filter_by(probe_id=probe_id).\
				filter_by(sequence_fragment_id=sequence_fragment_obj.id).\
				filter_by(fragment_start=fragment_start).filter_by(fragment_stop=fragment_stop)
			if query.count()==0:
				new_obj = True
			else:
				new_obj = False
		else:
			new_obj = True
		
		if new_obj:
			sequence_fragment2probe = Stock_250kDB.SequenceFragment2Probe(chromosome=chr, start=start, stop=stop, \
										probe_id=probe_id, no_of_identities=no_of_identities,\
										fragment_start=fragment_start, fragment_stop=fragment_stop)
			sequence_fragment2probe.sequence_fragment = sequence_fragment_obj
		else:
			sequence_fragment2probe = None
		return sequence_fragment2probe

	def putQCIntoDB(self, session, input_fname, no_of_lines_to_skip, data_source_obj, cnv_type_obj, cnv_method_obj=None, \
				run_type=1, original_id=None, version=1):
		"""
		2009-10-28
		"""
		sys.stderr.write("Putting QC data into database ... \n")
		reader = csv.reader(open(input_fname), delimiter=figureOutDelimiter(input_fname))
		input_file_basename = os.path.basename(input_fname)
		
		if run_type==7:
			# 2010-6-14 need to skip first 4 lines (3 comment-lines + 1 header)  for nucmer coords file
			no_of_lines_to_skip = 4
		elif run_type==8:	# skip 2 lines (1 comment-line + 1 header) for breakdancer output from Quan Long 
			no_of_lines_to_skip = 2
		
		col_name2index = None
		
		for i in range(no_of_lines_to_skip):
			header = reader.next()
		
		# the last line to be skipped will be the header
		if col_name2index is None:
			col_name2index = getColName2IndexFromHeader(header)
		
		counter = 0
		for row in reader:
			if run_type==1:
				db_obj = self.generateCNVQCCallObjFromClark2007(session, row, data_source_obj, cnv_type_obj, cnv_method_obj)
			elif run_type==2:
				db_obj = self.generateCNVQCCallObjFromSchneebergerOssowski(session, row, data_source_obj, cnv_type_obj, cnv_method_obj)
			elif run_type==3:
				db_obj = self.generateCNVQCCallObjFromBobSchmitzData(session, row, data_source_obj, cnv_type_obj,\
																		cnv_method_obj, original_id=original_id)
			elif run_type==4:
				db_obj = self.generateCNVQCCallObjFromLerContigDerivedCNVs(session, row, data_source_obj, cnv_type_obj, \
															cnv_method_obj=cnv_method_obj, \
															original_id=original_id, col_name2index=col_name2index)
			
			elif run_type==5:
				db_obj = self.generateSequenceFragmentRefPosObjFromLerContigSpansOverCol(session, row, data_source_obj, cnv_type_obj, \
															cnv_method_obj=cnv_method_obj, \
															original_id=original_id, col_name2index=col_name2index, version=version)
			elif run_type==6:
				db_obj = self.generateSequenceFragment2ProbeObj(session, row, data_source_obj, cnv_type_obj, \
														cnv_method_obj=cnv_method_obj, \
														original_id=original_id, col_name2index=col_name2index)
			elif run_type==7:
				db_obj = self.generateSequenceFragmentRefPosObjFromNucmerLerContigSpansOverCol(session, row, data_source_obj, \
													cnv_type_obj, \
													cnv_method_obj=cnv_method_obj, \
													original_id=original_id, col_name2index=col_name2index,\
													version=version, comment=input_file_basename)
			elif run_type==8:
				db_obj = self.generateCNVQCCallObjFromQuanLongBreakDancerOutput(session, row, data_source_obj, \
													cnv_type_obj, cnv_method_obj=cnv_method_obj,\
													original_id=original_id, col_name2index=col_name2index)
			elif run_type==9:
				db_obj = self.generateCNVQCCallObjFromQuanLongCoverageDerived(session, row, data_source_obj, \
												cnv_type_obj=cnv_type_obj, \
												cnv_method_obj=cnv_method_obj, 
												col_name2index=col_name2index)
			else:
				sys.stderr.write("Run type %s not supported.\n"%run_type)
			if db_obj:
				session.add(db_obj)
				session.flush()
			counter += 1
			if counter%5000==0:
				sys.stderr.write("%s%s"%('\x08'*40, counter))
		session.flush()
		sys.stderr.write("%s records. Done.\n"%counter)
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		
		if self.run_type==3 and  not self.original_id:
			sys.stderr.write("Need original_id for run_type %s.\n"%self.run_type)
			sys.exit(2)
		
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				  			password=self.db_passwd, hostname=self.hostname, database=self.dbname, 
				   			schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		#session.begin()
		
		data_source_obj = self.getDBObj(session, Stock_250kDB.DataSource, self.data_source)
		cnv_type_obj = self.getDBObj(session, Stock_250kDB.CNVType, self.cnv_type)
		if self.cnv_method:
			cnv_method_obj = self.getDBObj(session, Stock_250kDB.CNVMethod, self.cnv_method)
		else:
			cnv_method_obj = None
		self.putQCIntoDB(session, self.input_fname, self.no_of_lines_to_skip, data_source_obj, cnv_type_obj, \
						cnv_method_obj=cnv_method_obj, run_type=self.run_type, original_id=self.original_id,\
						version=self.version)
		
		if self.commit:
			session.commit()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = CNVQCConstruct
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()