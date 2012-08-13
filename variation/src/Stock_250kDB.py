#!/usr/bin/env python
"""
Examples:
	#setup database in postgresql
	Stock_250kDB.py -v postgres -u crocea -z localhost -d graphdb -k public
	
	#setup database in mysql
	Stock_250kDB.py -u yh -z papaya.usc.edu
	
Description:
	2008-07-09
	This is a wrapper for the stock_250k database, build on top of elixir. supposed to supercede the table definitions in mysql.sql.
"""
import sys, os, math
from sqlalchemy.types import LargeBinary

bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from sqlalchemy.engine.url import URL
from elixir import Unicode, DateTime, String, Integer, UnicodeText, Text, Boolean, Float, Binary, Enum
from elixir import Entity, Field, using_options, using_table_options
from elixir import OneToMany, ManyToOne, ManyToMany
from elixir import setup_all, session, metadata, entities
from elixir.options import using_table_options_handler	#using_table_options() can only work inside Entity-inherited class.
from sqlalchemy import UniqueConstraint, create_engine
from sqlalchemy.schema import ThreadLocalMetaData, MetaData
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy import and_, or_, not_

from datetime import datetime

from pymodule import SNPData, PassingData
from pymodule.db import ElixirDB, TableClass
import os
import hashlib 


__session__ = scoped_session(sessionmaker(autoflush=False, autocommit=True))
#__metadata__ = ThreadLocalMetaData()	#2008-10 not good for pylon

__metadata__ = MetaData()

class README(Entity, TableClass):
	#2008-08-07
	title = Field(String(2000))
	description = Field(String(60000))
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='readme', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class Phenotype(Entity):
	"""
	2009-8-12
		add created_by, updated_by, date_created, date_updated
		add unique constraint ('ecotype_id', 'method_id', 'replicate')
	2009-6-2
		point phenotype_method to the correct class name: PhenotypeMethod, rather than Phenotype
	"""
	ecotype_id = Field(Integer, nullable=False)
	value = Field(Float)
	replicate = Field(Integer)
	phenotype_method = ManyToOne('%s.PhenotypeMethod'%__name__, colname='method_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='phenotype', metadata=__metadata__, session=__session__)
	using_table_options(UniqueConstraint('ecotype_id', 'method_id', 'replicate'))
	using_table_options(mysql_engine='InnoDB')

class PhenotypeAvg(Entity):
	"""
	2009-8-18
		rename stdev to stddev
	2009-7-30
		add created_by, updated_by, date_created, date_updated
		add unique constraint (ecotype_id, method_id)
	2009-6-2
		point phenotype_method to the correct class name: PhenotypeMethod, rather than Phenotype
	"""
	ecotype_id = Field(Integer, nullable=False)
	value = Field(Float)
	stddev = Field(Float)
	sample_size = Field(Integer)
	ready_for_publication = Field(Integer, default=0)
	phenotype_method = ManyToOne('%s.PhenotypeMethod'%__name__, colname='method_id', ondelete='CASCADE', onupdate='CASCADE')
	readme = ManyToOne("%s.README"%__name__, colname='readme_id', ondelete='CASCADE', onupdate='CASCADE')
	transformed_value = Field(Float)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='phenotype_avg', metadata=__metadata__, session=__session__)
	using_table_options(UniqueConstraint('ecotype_id', 'method_id'))
	using_table_options(mysql_engine='InnoDB')

class BiologyCategory(Entity):
	"""
	2012.3.22
		add column preferred_gene_list
	#2008-08-21
	"""
	short_name = Field(String(256), unique=True)
	description = Field(String(8192))
	gene_list_type_ls = OneToMany("%s.GeneListType"%__name__)
	phenotype_method_ls = OneToMany("%s.PhenotypeMethod"%__name__)
	preferred_list_type_id = Field(Integer)	#2012.7.9 not a foreign key anymore as it causes cyclic dependency (check GeneListType)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='biology_category', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	
	def returnGeneListIDList(self):
		"""
		2012.3.23
		
		"""
		if self.preferred_list_type_id:
			return [self.preferred_list_type_id]
		else:
			list_type_id_ls = []
			for list_type in self.gene_list_type_ls:
				list_type_id_ls.append(list_type.id)
			return list_type_id_ls

class GeneListType(Entity):
	short_name = Field(String(256))
	biology_category = ManyToOne("%s.BiologyCategory"%__name__, colname='biology_category_id', ondelete='CASCADE', onupdate='CASCADE')
	type = ManyToOne("%s.GeneListSuperType"%__name__, colname='super_type_id', ondelete='CASCADE', onupdate='CASCADE')
	original_filename = Field(String(760), unique=True)	#for unique constraint in mysql, max key length is 767 bytes
	description = Field(String(8192))
	gene_list = OneToMany('%s.GeneList'%__name__)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_gene_list_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('short_name', 'super_type_id'))

class GeneListSuperType(Entity):
	#2008-08-22
	short_name = Field(String(256), unique=True)
	description = Field(String(8192))
	list_type_ls = OneToMany(GeneListType)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_list_super_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class GeneList(Entity):
	"""
	2012.3.26
		add ncbi_gene_id
	"""
	gene_id = Field(Integer)
	ncbi_gene_id = Field(Integer)	#2012.3.26
	list_type = ManyToOne('%s.GeneListType'%__name__, colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE', inverse='gene_list')
	original_name = Field(String(128))
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_gene_list', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('gene_id', 'list_type_id'))

class GenomeWideResultMethod(Entity):
	"""
	2010-10-10
	"""
	short_name = Field(String(256), unique=True)
	description = Field(String(8192))
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='genome_wide_result_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class GenomeMarker(Entity):
	"""
	2011-4-19
		change type of chromosome from int to string
	2010-10-10
	"""
	chromosome = Field(String(256))	#2011-4-19
	start = Field(Integer)
	stop = Field(Integer)
	description = Field(String(512), deferred=True)
	object = Field(LargeBinary(1000000), deferred=True)	#a python dictionary to store other attributes
	created_by = Field(String(128), deferred=True)
	updated_by = Field(String(128), deferred=True)
	date_created = Field(DateTime, default=datetime.now, deferred=True)
	date_updated = Field(DateTime, deferred=True)
	using_options(tablename='genome_marker', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('chromosome', 'start', 'stop'))
	

class GenomeWideResult(Entity):
	"""
	2010-10-10
		a general table to store all kinds of genomic data
	"""
	method = ManyToOne('%s.GenomeWideResultMethod'%__name__, colname='method_id', ondelete='CASCADE', onupdate='CASCADE')
	marker = ManyToOne('%s.GenomeMarker'%__name__, colname='marker_id', ondelete='CASCADE', onupdate='CASCADE')
	value = Field(Float)
	rank = Field(Integer)
	object = Field(LargeBinary(13217728), deferred=True)	#a python dictionary to store other attributes
	created_by = Field(String(128), deferred=True)
	updated_by = Field(String(128), deferred=True)
	date_created = Field(DateTime, default=datetime.now, deferred=True)
	date_updated = Field(DateTime, deferred=True)
	using_options(tablename='genome_wide_result', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('method_id', 'marker_id'))
	

class Snps(Entity):
	"""
	2011-4-19
		change type of chromosome from int to string
	2011-4-15
		add column cnv, for linking up
	2010-10-19
		add column include_after_qc, Field(Integer):
			indicating whether this SNP is still included in final dataset after QC.
	2010-6-17
		add argument tair8_chromosome, tair8_position
	"""
	name = Field(String(200), nullable = False, deferred=True)
	chromosome = Field(String(256))
	position = Field(Integer)
	end_position = Field(Integer)
	allele1 = Field(String(1))
	allele2 = Field(String(2))
	tair8_chromosome = Field(String(256))
	tair8_position = Field(Integer)
	include_after_qc = Field(Integer, default=0)
	locus_type = ManyToOne('%s.LocusType'%__name__, colname='locus_type_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200), deferred=True)
	updated_by = Field(String(200), deferred=True)
	date_created = Field(DateTime, default=datetime.now, deferred=True)
	date_updated = Field(DateTime, deferred=True)
	using_options(tablename='snps', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('name', 'chromosome', 'position', 'end_position', 'locus_type_id'))

class LocusType(Entity):
	"""
	2012.3.7
		used to separate different sets of loci in table Snps
	"""
	short_name = Field(String(256), unique=True)
	description = Field(String(8192))
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='locus_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class SnpsContext(Entity):
	snp = ManyToOne('%s.Snps'%__name__, colname='snps_id', ondelete='CASCADE', onupdate='CASCADE')
	disp_pos = Field(Integer)
	gene_id = Field(Integer)
	gene_strand = Field(String(1))
	left_or_right = Field(String(200), deferred=True)
	disp_pos_comment = Field(String(2000), deferred=True)
	created_by = Field(String(200), deferred=True)
	updated_by = Field(String(200), deferred=True)
	date_created = Field(DateTime, default=datetime.now, deferred=True)
	date_updated = Field(DateTime, deferred=True)
	using_options(tablename='snps_context', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('snps_id', 'gene_id'))

class CallMethod(Entity):
	"""
	2012.3.9
		add column locus_type
	2012.2.28
		add column no_of_accessions, no_of_loci
	2010-2-22
		add accession_set_id
		make column parent_id a foreign key to CallMethod
		add a unique constraint (parent_id, accession_set_id)
	2009-6-17
		add column parent_id
	2009-3-12
		add column filename
	2009-2-22
		add field unique_ecotype
	2009-2-6
		add eigen_value_ls
	2009-1-7
		add various columns to denote QC parameters
	"""
	short_name = Field(String(20))
	parent_call_method = ManyToOne('%s.CallMethod'%__name__, colname='parent_id', ondelete='CASCADE', onupdate='CASCADE')
	accession_set = ManyToOne('%s.AccessionSet'%__name__, colname='accession_set_id', ondelete='CASCADE', onupdate='CASCADE')
	locus_type = ManyToOne('%s.LocusType'%__name__, colname='locus_type_id', ondelete='CASCADE', onupdate='CASCADE')
	min_oligo_call_prob = Field(Float)
	max_array_mismatch_rate = Field(Float)
	max_array_NA_rate = Field(Float)
	max_snp_mismatch_rate = Field(Float)
	max_snp_NA_rate = Field(Float)
	npute_window_size = Field(Integer)
	no_of_accessions = Field(Integer)	#2012.2.28
	no_of_loci = Field(Integer)	#2012.2.28
	imputed = Field(Boolean)
	unique_ecotype = Field(Boolean)
	avg_array_mismatch_rate = Field(Float)
	filename = Field(String(2024))
	method_description = Field(String(8000))
	data_description = Field(String(8000))
	comment = Field(String(8000))
	call_info_ls = OneToMany("CallInfo")
	eigen_value_ls = OneToMany("%s.CallMethodEigenValue"%__name__, order_by='which_eigen')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='call_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('parent_id', 'accession_set_id', 'locus_type_id'))
	
class PhenotypeMethod(Entity):
	"""
	2009-9-1
		add min_value, stddev to facilitate transformations that utilize these two values
		 
	2009-5-11
		add no_of_accessions, growth_condition, phenotype_scoring, citations
	"""
	short_name = Field(String(40), unique=True)
	only_first_96 = Field(Boolean, default=0)
	biology_category = ManyToOne("BiologyCategory", colname='biology_category_id', ondelete='CASCADE', onupdate='CASCADE')
	no_of_accessions = Field(Integer);
	growth_condition = Field(String(8000))
	phenotype_scoring = Field(String(8000))
	citations = Field(String(8000))
	method_description = Field(String(8000))
	data_description = Field(String(8000))
	min_value = Field(Float)
	stddev = Field(Float)
	comment = Field(String(8000))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	data_type = Field(String(200))
	transformation_description = Field(String(8000))
	owner = ManyToOne("Users",colname='owner_id')
	access = Field(Enum(('public','restricted'),default='restricted'))
	groups = ManyToMany('Groups',tablename='acl_phenotype_method_groups',remote_colname='groups_id')
	users = ManyToMany('Users',tablename='acl_phenotype_method_users',remote_colname='users_id')

	using_options(tablename='phenotype_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	
	@classmethod
	def getPhenotypesForACL(cls,biology_category_id = None,user = None):
		#
		query = PhenotypeMethod.query
		if biology_category_id is not None:
			query = query.filter(PhenotypeMethod.biology_category_id==biology_category_id)

		clause = and_(PhenotypeMethod.access == 'public')
		if user is not None:
			if user.isAdmin == 'Y':
				return query
			clause = or_(clause,PhenotypeMethod.owner == user, PhenotypeMethod.users.any(Users.id == user.id),PhenotypeMethod.groups.any(Groups.id.in_([group.id for group in user.groups])))
		query = query.filter(clause)
		return query
	
	def checkACL(self,user):
		if self.access == 'public':
			return True
		if user is None:
			return False
		if user.isAdmin == 'Y':
			return True
		if self.owner == user: 
			return True
		if user in self.users:
			return True
		if [group in self.groups for group in user.groups]: 
			return True
		return False
	
	def getProperShortName(self):
		"""
		2012.3.11
			replace <i>, </i> with emptry string
		"""
		tmp_str = self.short_name.replace("<i>", "")
		tmp_str = tmp_str.replace("</i>", "")
		tmp_str = tmp_str.replace('<', '')
		tmp_str = tmp_str.replace('>', '')
		tmp_str = tmp_str.replace('/', '')
		return tmp_str
		
		
class AnalysisMethod(Entity):
	"""
	2012.6.5
		add column association_test_type
	2009-5-14
		add min_MAF
	2008-09-16
		add smaller_score_more_significant
	"""
	short_name = Field(String(120))
	method_description = Field(String(8000))
	min_maf = Field(Float)
	smaller_score_more_significant = Field(Integer)
	association_test_type = Field(Integer)	#2012.6.5 corresponds to Association.py's test_type
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='analysis_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	
class ResultsMethodType(Entity):
	short_name = Field(String(30), unique=True)
	method_description = Field(String(8000))
	data_description = Field(String(8000))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='results_method_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class AccessionSet(Entity):
	"""
	2010-1-15
		the table which records the set of accessions, like first 96, 192, 384, ...
	"""
	short_name = Field(String(60), unique=True)
	description = Field(String(8000))
	comment = Field(String(8000))
	ecotype_ls = OneToMany('%s.AccessionSet2Ecotype'%__name__)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='accession_set', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class AccessionSet2Ecotype(Entity):
	"""
	2010-1-15
		records which ecotype goes with which accession set
	"""
	accession_set = ManyToOne('%s.AccessionSet'%__name__, colname='accession_set_id', ondelete='CASCADE', onupdate='CASCADE')
	ecotype_id = Field(Integer)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='accession_set2ecotype', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('accession_set_id', 'ecotype_id'))

class ResultsMethod(Entity):
	"""
	2012.3.7
		add column locus_type
	2011-2-6
		add column cnv_method to hold CNV association results
		change unique constraint to include cnv_method_id
	2010-9-21
		add remove_outliers, pseudo_heritability, transformation_parameters
	2010-1-22
		add unique constraint
	"""
	short_name = Field(String(60), unique=True)
	filename = Field(String(1000), unique=True)
	original_filename = Field(Text)
	method_description = Field(String(8000))
	data_description = Field(String(8000))
	comment = Field(String(8000))
	phenotype_method = ManyToOne('%s.PhenotypeMethod'%__name__, colname='phenotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	call_method = ManyToOne('%s.CallMethod'%__name__, colname='call_method_id', ondelete='CASCADE', onupdate='CASCADE')
	results_method_type = ManyToOne('%s.ResultsMethodType'%__name__, colname='results_method_type_id', ondelete='CASCADE', onupdate='CASCADE')
	analysis_method = ManyToOne('%s.AnalysisMethod'%__name__, colname='analysis_method_id', ondelete='CASCADE', onupdate='CASCADE')
	transformation_method = ManyToOne('%s.TransformationMethod'%__name__, colname='transformation_method_id', ondelete='CASCADE', onupdate='CASCADE')
	cnv_method = ManyToOne("%s.CNVMethod"%__name__, colname='cnv_method_id', ondelete='CASCADE', onupdate='CASCADE')
	locus_type = ManyToOne('%s.LocusType'%__name__, colname='locus_type_id', ondelete='CASCADE', onupdate='CASCADE')
	rm_json = OneToMany('%s.ResultsMethodJson'%__name__)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	remove_outliers = Field(Integer, default=0)
	pseudo_heritability = Field(Float)
	transformation_parameters = Field(String(11))
	using_options(tablename='results_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('call_method_id', 'phenotype_method_id', \
										'results_method_type_id', 'analysis_method_id', 'cnv_method_id'))
	
	def getDateStampedFilename(self):
		"""
		2012.3.21
			xxx.tsv => xxx.2012_3_21.tsv
		"""
		from datetime import datetime
		lastModDatetime = datetime.fromtimestamp(os.stat(self.filename).st_mtime)
		prefix, suffix = os.path.splitext(self.filename)
		newFilename = '%s.%s_%s_%s%s'%(prefix, lastModDatetime.year, lastModDatetime.month,\
									lastModDatetime.day, suffix)
		return newFilename
class ResultsMethodJson(Entity):
	"""
	2009-5-4
		table to store json structures sent over to clients
	"""
	result = ManyToOne('ResultsMethod', colname='results_id', ondelete='CASCADE', onupdate='CASCADE')
	no_of_top_snps = Field(Integer)
	min_MAF = Field(Float)
	json_data = Field(LargeBinary(134217728), deferred=True)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='results_method_json', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('results_id', 'no_of_top_snps', 'min_MAF'))

class ResultLandscape(Entity):
	"""
	2011-3-28
		table to store the landscape of association result
	"""
	result = ManyToOne('ResultsMethod', colname='result_id', ondelete='CASCADE', onupdate='CASCADE')
	start_locus = ManyToOne('Snps', colname='start_locus_id', ondelete='CASCADE', onupdate='CASCADE')
	stop_locus = ManyToOne('Snps', colname='stop_locus_id', ondelete='CASCADE', onupdate='CASCADE')
	no_of_loci = Field(Integer)	#number of loci in between start_locus and stop_locus
	neighbor_distance = Field(Integer)
	comment = Field(Text)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='result_landscape', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('result_id', 'start_locus_id', 'stop_locus_id', 'neighbor_distance'))

class ResultPeak(Entity):
	"""
	2012.6.22 add column association_locus
	2011-4-19
		table to store the peaks of association result
	"""
	result = ManyToOne('%s.ResultsMethod'%__name__, colname='result_id', ondelete='CASCADE', onupdate='CASCADE')
	result_peak_type = ManyToOne('%s.ResultPeakType'%__name__, colname='result_peak_type_id', ondelete='CASCADE', onupdate='CASCADE')
	chromosome = Field(String(256))
	start = Field(Integer)
	stop = Field(Integer)
	start_locus = ManyToOne('%s.Snps'%__name__, colname='start_locus_id', ondelete='CASCADE', onupdate='CASCADE')
	stop_locus = ManyToOne('%s.Snps'%__name__, colname='stop_locus_id', ondelete='CASCADE', onupdate='CASCADE')
	no_of_loci = Field(Integer)	#number of loci in between stop_locus of start_bridge and start_locus of stop_bridge,
	# including the two as well.
	peak_locus = ManyToOne('%s.Snps'%__name__, colname='peak_locus_id', ondelete='CASCADE', onupdate='CASCADE')
	peak_score = Field(Float)
	association_locus = ManyToOne('%s.AssociationLocus'%__name__, colname='association_locus_id', ondelete='CASCADE', onupdate='CASCADE')
	comment = Field(Text)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='result_peak', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('result_id', 'result_peak_type_id', 'chromosome', 'start', 'stop'))


class AssociationLocus(Entity):
	#2012.6.22
	chromosome = Field(String(256))
	start = Field(Integer)
	stop = Field(Integer)
	threshold = Field(Float)
	connectivity = Field(Float)
	no_of_peaks = Field(Integer)
	comment = Field(Text)
	result_peak_ls = OneToMany("%s.ResultPeak"%(__name__))
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='association_locus', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('chromosome', 'start', 'stop', 'threshold'))


class AssociationLocusStat(Entity):
	#2012.6.22
	association_locus = ManyToOne('%s.AssociationLocus'%__name__, colname='association_locus_id', ondelete='CASCADE', onupdate='CASCADE')
	no_of_phenotypes = Field(Integer)
	biology_category = ManyToOne("%s.BiologyCategory"%__name__, colname='biology_category_id', ondelete='CASCADE', onupdate='CASCADE')
	call_method = ManyToOne('%s.CallMethod'%__name__, colname='call_method_id', ondelete='CASCADE', onupdate='CASCADE')
	analysis_method = ManyToOne('%s.AnalysisMethod'%__name__, colname='analysis_method_id', ondelete='CASCADE', onupdate='CASCADE')
	result_peak_type = ManyToOne('%s.ResultPeakType'%__name__, colname='result_peak_type_id', ondelete='CASCADE', onupdate='CASCADE')
	#result_peak_type_id = Field(Integer)
	comment = Field(Text)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='association_locus_stat', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('association_locus_id', 'biology_category_id', 'call_method_id', 'analysis_method_id', \
										'result_peak_type_id'))



class ResultPeakType(Entity):
	"""
	2012.2.15
		add field, min_MAF
	2011-4-19
		type for ResultPeak
	"""
	short_name = Field(String(30), unique=True)
	description = Field(Text)
	min_MAF = Field(Float)
	min_score = Field(Float)
	neighbor_distance = Field(Integer)
	max_neighbor_distance = Field(Integer)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='result_peak_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('min_MAF', 'min_score', 'neighbor_distance', 'max_neighbor_distance'))

class ResultPeakOverlap(Entity):
	"""
	2012.2.22
		
	"""
	result1 = ManyToOne('ResultsMethod', colname='result1_id', ondelete='CASCADE', onupdate='CASCADE')
	result2 = ManyToOne('ResultsMethod', colname='result2_id', ondelete='CASCADE', onupdate='CASCADE')
	
	result1_peak_type = ManyToOne('ResultPeakType', colname='result1_peak_type_id', ondelete='CASCADE', onupdate='CASCADE')
	result2_peak_type = ManyToOne('ResultPeakType', colname='result2_peak_type_id', ondelete='CASCADE', onupdate='CASCADE')
	peak_padding = Field(Integer)	#this is the extension of peak width in checking peak overlap to accommodate LD
	
	no_of_result1_peaks = Field(Integer)
	no_of_result2_peaks = Field(Integer)
	no_of_result1_peaks_in_result2 = Field(Integer)
	no_of_result2_peaks_in_result1 = Field(Integer)
	
	comment = Field(Text)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='result_peak_overlap', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('result1_id', 'result2_id', 'result1_peak_type_id', 'result2_peak_type_id', 'peak_padding'))

class Results(Entity):
	"""
	2009-2-22
		change column name of "results_method" from results_method_id to results_id
	2009-2-6
		a table to store data derived from filenames stored in results_method
		add more fields
	"""
	snp = ManyToOne('Snps', colname='snps_id', ondelete='CASCADE', onupdate='CASCADE')
	result = ManyToOne('ResultsMethod', colname='results_id', ondelete='CASCADE', onupdate='CASCADE')
	score = Field(Float)
	rank = Field(Integer)
	beta = Field(Float)
	maf = Field(Float)
	mac = Field(Integer)
	genotype_var_perc = Field(Float)
	correlation = Field(Float)
	odds_ratio = Field(Float)
	statistic = Field(Float)
	object = Field(LargeBinary(134217728), deferred=True)	#a python dictionary to store other attributes
	using_options(tablename='results', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('results_id', 'snps_id'))
	
	@classmethod
	def getResultsForACL(cls,snp,call_method_id,user = None):
		query = Results.query.filter_by(snps_id=snp.id).join(ResultsMethod).join(PhenotypeMethod).filter(ResultsMethod.call_method_id == call_method_id)
		clause = and_(PhenotypeMethod.access == 'public')
		if user is not None:
			if user.isAdmin == 'Y':
				return query
			clause = or_(clause,PhenotypeMethod.owner == user, PhenotypeMethod.users.any(Users.id == user.id),PhenotypeMethod.groups.any(Groups.id.in_([group.id for group in user.groups])))
		query = query.filter(clause)
		return query

class CandidateGeneRankSumTestResult(Entity):
	"""
	2008-10-09
		add test_type = Field(Integer) to conform to CandidateGeneRankSumTestResultMethod
		rename results_by_gene to result
		rename results_by_gene_id to results_id
	2008-09-16
		table linked to results_by_gene, not results_method
	2008-07-17
	"""
	result = ManyToOne('ResultsByGene', colname='results_id', ondelete='CASCADE', onupdate='CASCADE')
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE')
	statistic = Field(Float)
	pvalue = Field(Float)
	min_distance = Field(Integer)
	get_closest = Field(Integer)
	min_MAF = Field(Float)
	max_pvalue_per_gene = Field(Integer)
	candidate_sample_size = Field(Integer)
	non_candidate_sample_size = Field(Integer)
	test_type = Field(Integer)
	comment = Field(Text)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_gene_rank_sum_test_rbg', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('results_method_id', 'list_type_id'))

class CandidateGeneRankSumTestResultMethod(Entity):
	"""
	2009-9-16
		link test_type to CandidateGeneTopSNPTestRMType
	2008-10-09
		similar CandidateGeneRankSumTestResult linked to results_method

	"""
	result = ManyToOne('ResultsMethod', colname='results_id', ondelete='CASCADE', onupdate='CASCADE')
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE')
	statistic = Field(Float)
	pvalue = Field(Float)
	min_distance = Field(Integer)
	get_closest = Field(Integer)
	min_MAF = Field(Float)
	max_pvalue_per_gene = Field(Integer)
	candidate_sample_size = Field(Integer)
	non_candidate_sample_size = Field(Integer)
	#test_type = Field(Integer)
	test_type = ManyToOne('CandidateGeneTopSNPTestRMType', colname='test_type_id', ondelete='CASCADE', onupdate='CASCADE')
	comment = Field(Text)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_gene_rank_sum_test_rm', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('results_id', 'list_type_id', 'min_distance', 'get_closest', 'min_MAF', 'test_type_id'))

class CandidateGeneRankSumTestResultMethodType(Entity):
	"""
	2009-9-16
		reference for the test_type_id field of CandidateGeneRankSumTestResultMethod
	"""
	gw_looping_type = Field(Integer)	#whether chromosomes shall be shuffled or not before SNPs are looped
	allow_two_sample_overlapping = Field(Integer)
	results_type = Field(Integer)
	test_type = ManyToOne('AnalysisMethod', colname='test_type_id', ondelete='CASCADE', onupdate='CASCADE')
	null_distribution_type = ManyToOne('NullDistributionType', colname='null_distribution_type_id', ondelete='CASCADE', onupdate='CASCADE')
	how_to_handle_rank = Field(String(200))
	comment = Field(Text)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_gene_rst_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('gw_looping_type', \
										'allow_two_sample_overlapping', 'results_type', 'test_type_id',\
										'null_distribution_type_id', 'how_to_handle_rank'))


class ResultsByGene(Entity):
	"""
	2008-11-14
		more sensible unique constraints
	2008-09-15
		modify it to contain gene-based results(files) deriving from snp-based results in table ResultsMethod
	2008-08-27
		add readme_id
		modify unique constraint to include readme_id
	2008-07-19
	"""
	short_name = Field(String(60))
	filename = Field(String(767), unique=True)
	min_distance = Field(Integer)
	get_closest = Field(Integer)
	min_MAF = Field(Float)
	comment = Field(Text)
	results_method = ManyToOne('ResultsMethod', colname='results_method_id', ondelete='CASCADE', onupdate='CASCADE')
	rank_test = OneToMany("CandidateGeneRankSumTestResult")
	top_snp_test = OneToMany("CandidateGeneTopSNPTest")
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='results_by_gene', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('results_method_id', 'min_distance', 'get_closest', 'min_MAF'))

class SnpsQC(Entity):
	"""
	2008-08-07
	"""
	snp = ManyToOne('Snps', colname='snps_id', ondelete='CASCADE', onupdate='CASCADE')
	min_probability = Field(Float)
	max_call_info_mismatch_rate = Field(Float)
	snps_name = Field(String(200), deferred=True)
	tg_snps_name = Field(String(200), deferred=True)
	NA_rate = Field(Float)
	no_of_NAs = Field(Integer)
	no_of_totals = Field(Integer)
	relative_NA_rate = Field(Float)
	relative_no_of_NAs = Field(Integer)
	relative_no_of_totals = Field(Integer)
	mismatch_rate = Field(Float)
	no_of_mismatches = Field(Integer)
	no_of_non_NA_pairs = Field(Integer)
	qc_method = ManyToOne("%s.QCMethod"%__name__, colname='qc_method_id', ondelete='CASCADE', onupdate='CASCADE')
	call_method = ManyToOne("CallMethod", colname='call_method_id', ondelete='CASCADE', onupdate='CASCADE')
	readme = ManyToOne("%s.README"%__name__, colname='readme_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200), deferred=True)
	updated_by = Field(String(200), deferred=True)
	date_created = Field(DateTime, default=datetime.now, deferred=True)
	date_updated = Field(DateTime, deferred=True)
	using_options(tablename='snps_qc', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class QCMethod(Entity):
	"""
	2009-5-20
		add column snps_table, ignore_het
	"""
	short_name = Field(String(30), unique=True)
	data1_type = Field(String(30), nullable=False)
	data2_type = Field(String(30), nullable=False)
	snps_table = Field(String(1000))	#db table holds information about all SNPs of data2_type
	ignore_het = Field(Integer)	# whether this QC ignores heterozygous in both data1_type and data2_type
	method_description = Field(String(8000))
	data_description = Field(String(8000))
	comment = Field(String(8000))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='qc_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class ContaminantType(Entity):
	"""
	2008-08-07
	"""
	short_name = Field(String(100), unique=True)
	description = Field(Text)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='contaminant_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class ArrayInfo(Entity):
	"""
	2009-11-25
		add array_quartile_ls
	2009-3-25
		add median_intensity
	"""
	name = Field(String(40))
	filename = Field(String(1000))
	original_filename = Field(String(1000))
	description = Field(String(2000))
	ecotype_id = Field(Integer)
	maternal_ecotype_id = Field(Integer)
	paternal_ecotype_id = Field(Integer)
	contaminant_type = ManyToOne("%s.ContaminantType"%__name__, colname='contaminant_type_id', ondelete='CASCADE', onupdate='CASCADE')
	md5sum = Field(String(100))
	experimenter = Field(String(200))
	median_intensity = Field(Float)
	samples = Field(String(20))
	dna_amount = Field(String(20))
	S260_280 = Field(Float)
	total_vol = Field(String(20))
	hyb_vol = Field(String(20))
	seed_source = Field(String(100))
	method_name = Field(String(250))
	readme = ManyToOne("%s.README"%__name__, colname='readme_id', ondelete='CASCADE', onupdate='CASCADE')
	array_quartile_ls = OneToMany("%s.ArrayQuartile"%__name__)
	created_by = Field(String(200), deferred=True)
	updated_by = Field(String(200), deferred=True)
	date_created = Field(DateTime, default=datetime.now, deferred=True)
	date_updated = Field(DateTime, deferred=True)
	using_options(tablename='array_info', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class ArrayQuartile(Entity):
	"""
	2010-3-1
		add array_quartile_outlier_ls as a link to ArrayQuartileOutlier
	2009-11-24
		add dlrspread (DLRSpread) = IQR(dLR) / 4*erfinv(0.5), where dLR (derivative log ratio) is an array of
		differences between log ratios of adjacent probes, erfinv is the Inverse Error
		Function and IQR is Inter Quartile Range. Calculated by DB_250k2Array.py.
	2009-11-20
		a table storing quartiles of arrays in chosen segments
		DB_250k2Array.py fills this table.
	"""
	array = ManyToOne("%s.ArrayInfo"%__name__, colname='array_id', ondelete='CASCADE', onupdate='CASCADE')
	start_probe = ManyToOne("%s.Probes"%__name__, colname='start_probe_id', ondelete='CASCADE', onupdate='CASCADE')
	stop_probe = ManyToOne("%s.Probes"%__name__, colname='stop_probe_id', ondelete='CASCADE', onupdate='CASCADE')
	no_of_probes = Field(Integer)
	minimum = Field(Float)
	first_decile = Field(Float)
	lower_quartile = Field(Float)
	median = Field(Float)
	upper_quartile = Field(Float)
	last_decile = Field(Float)
	maximum = Field(Float)
	dlrspread = Field(Float)	# 2009-11-24 calculated using tiling probes only
	
	dlr_minimum = Field(Float)	# 2009-11-25	7-number summary for dlr (derivative log ratio)
	dlr_first_decile = Field(Float)
	dlr_lower_quartile = Field(Float)
	dlr_median = Field(Float)
	dlr_upper_quartile = Field(Float)
	dlr_last_decile = Field(Float)
	dlr_maximum = Field(Float)
	
	array_quartile_outlier_ls = OneToMany("%s.ArrayQuartileOutlier"%__name__)
	
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='array_quartile', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('array_id', 'start_probe_id', 'stop_probe_id'))

class ArrayQuartileOutlier(Entity):
	"""
	2009-11-20
		probes whose intensity is either below lower_quartile-1.5IQR or above upper_quartile+1.5IQR
	"""
	array_quartile = ManyToOne("%s.ArrayQuartile"%__name__, colname='array_quartile_id', ondelete='CASCADE', onupdate='CASCADE')
	probe = ManyToOne("%s.Probes"%__name__, colname='probe_id', ondelete='CASCADE', onupdate='CASCADE')
	intensity = Field(Float)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='array_quartile_outlier', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('array_quartile_id', 'probe_id'))

class CallInfo(Entity):
	"""
	2009-2-2
		add pc_value_ls
	"""
	filename = Field(String(1000))
	description = Field(String(2000))
	array = ManyToOne("%s.ArrayInfo"%__name__, colname='array_id', ondelete='CASCADE', onupdate='CASCADE')
	NA_rate = Field(Float)
	readme = ManyToOne("%s.README"%__name__, colname='readme_id', ondelete='CASCADE', onupdate='CASCADE')
	call_method = ManyToOne('%s.CallMethod'%__name__, colname='method_id', ondelete='CASCADE', onupdate='CASCADE')
	call_qc_ls = OneToMany("%s.CallQC"%__name__)
	pc_value_ls = OneToMany("%s.CallInfoPCValue"%__name__, order_by='which_pc')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='call_info', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class CallQC(Entity):
	call_info = ManyToOne("%s.CallInfo"%__name__, colname='call_info_id', ondelete='CASCADE', onupdate='CASCADE')
	min_probability = Field(Float)
	ecotype_id = Field(Integer)
	duplicate = Field(Integer)
	tg_ecotype_id = Field(Integer)
	tg_duplicate = Field(Integer)
	NA_rate = Field(Float)
	no_of_NAs = Field(Integer)
	no_of_totals = Field(Integer)
	mismatch_rate = Field(Float)
	no_of_mismatches = Field(Integer)
	no_of_non_NA_pairs = Field(Integer)
	call_method = ManyToOne('CallMethod', colname='call_method_id', ondelete='CASCADE', onupdate='CASCADE')
	readme = ManyToOne("%s.README"%__name__, colname='readme_id', ondelete='CASCADE', onupdate='CASCADE')
	qc_method = ManyToOne("%s.QCMethod"%__name__, colname='qc_method_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='call_qc', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class Probes(Entity, TableClass):
	"""
	2011-4-19
		change type of chromosome from int to string
	2010-6-17
		add argument tair8_chromosome, tair8_position
	2010-5-23
		add field Tair9Copy, filled in by this function from misc.py
			tair9_blast_result_fname = os.path.expanduser("~/script/variation/data/CNV/cnv_probe_blast_against_tair9.tsv")
			CNV.markCNVProbeTAIR9CopyNumber(db_250k, tair9_blast_result_fname)
	"""
	snp = ManyToOne('Snps', colname='snps_id', ondelete='CASCADE', onupdate='CASCADE')
	seq = Field(String(25), deferred=True)
	chromosome = Field(String(256))
	position = Field(Integer)
	allele = Field(String(1))
	strand = Field(String(20))
	xpos = Field(Integer)
	ypos = Field(Integer)
	direction = Field(String(20), deferred=True)
	gene = Field(String(50), deferred=True)
	RNA = Field(String(50), deferred=True)
	tu = Field(String(50), deferred=True)
	flank = Field(String(50), deferred=True)
	expressedClones = Field(Float)
	totalClones = Field(Integer)
	multiTranscript = Field(String(50), deferred=True)
	LerDel = Field(String(50), deferred=True)
	LerCopy = Field(Integer)
	LerSNPdelL = Field(Integer)
	LerSNPdelR = Field(Integer)
	LerSNPpos = Field(Integer)
	Tair9Copy = Field(Integer)	# 2010-5-23 the copy number of a probe in TAIR9 version Col-0
	promoter = Field(Boolean)
	utr5 = Field(Boolean)
	utr3 = Field(Boolean)
	intron = Field(Boolean)
	intergenic = Field(Boolean)
	downstream = Field(Boolean)
	cda = Field(Boolean)
	tair8_chromosome = Field(String(256))
	tair8_position = Field(Integer)
	created_by = Field(String(200), deferred=True)
	updated_by = Field(String(200), deferred=True)
	date_created = Field(DateTime, default=datetime.now, deferred=True)
	date_updated = Field(DateTime, deferred=True)
	using_options(tablename='probes', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class CandidateGeneTopSNPTest(Entity):
	"""
	2008-10-15
		add candidate_sample_size, non_candidate_sample_size, starting_rank, test_type to be compatible with CandidateGeneTopSNPTestRM
	2008-09-16
		table linked to results_by_gene, not results_method
	2008-08-20
	"""
	result = ManyToOne('ResultsByGene', colname='results_id', ondelete='CASCADE', onupdate='CASCADE')
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE')
	pvalue = Field(Float)
	min_distance = Field(Integer)
	get_closest = Field(Integer)
	min_MAF = Field(Float)
	candidate_sample_size = Field(Integer)
	non_candidate_sample_size = Field(Integer)
	no_of_top_candidate_genes = Field(Integer)
	no_of_top_genes = Field(Integer)
	no_of_top_snps = Field(Integer)
	starting_rank = Field(Integer)
	test_type = Field(Integer)
	comment = Field(Text)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_gene_top_snp_test', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class CandidateGeneTopSNPTestRMType(Entity):
	"""
	2008-10-26
		move min_distance to CandidateGeneTopSNPTestRM
	2008-10-22
		hierarchical for CandidateGeneTopSNPTestRM
	"""
	get_closest = Field(Integer)
	min_MAF = Field(Float)
	allow_two_sample_overlapping = Field(Integer)
	results_type = Field(Integer)
	test_type = ManyToOne('AnalysisMethod', colname='test_type_id', ondelete='CASCADE', onupdate='CASCADE')
	null_distribution_type = ManyToOne('NullDistributionType', colname='null_distribution_type_id', ondelete='CASCADE', onupdate='CASCADE')	#2008-10-30 useless
	comment = Field(Text)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_gene_top_snp_test_rm_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('get_closest', 'min_MAF', \
										'allow_two_sample_overlapping', 'results_type', 'test_type_id',\
										'null_distribution_type_id'))
	
class CandidateGeneTopSNPTestRM(Entity):
	"""
	2008-10-31
		add pvalue_gw_looping, pvalue_random_gene_list (defunct, check doc at the side of their definitions).
		null_distribution_type_id in CandidateGeneTopSNPTestRMType doesn't mean much.
		NULL data under different simulations stored in TopSNPTestRMNullData.
	2008-10-26
		add min_distance, no_of_tests_passed, no_of_tests
	2008-10-23
		restructure it for the new varieties of tests and link to CandidateGeneTopSNPTestRMType
	2008-10-15
		similar to CandidateGeneTopSNPTest but stores test results from ResultsMethod
	"""
	result = ManyToOne('ResultsMethod', colname='results_id', ondelete='CASCADE', onupdate='CASCADE')
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE')
	pvalue = Field(Float)
	pvalue_gw_looping = Field(Float)	# 2008-10-30	defunct. CandidateGeneTopSNPTestRMNullData stores NULL data under different NULL distribution.
										# all pvalues are merged into one entry.
	pvalue_random_gene_list = Field(Float)	# 2008-10-30 defunct.
	candidate_sample_size = Field(Integer)
	non_candidate_sample_size = Field(Integer)
	candidate_gw_size = Field(Integer)
	non_candidate_gw_size = Field(Integer)
	no_of_top_candidate_genes = Field(Integer)
	no_of_top_genes = Field(Integer)
	no_of_top_snps = Field(Integer)
	min_distance = Field(Integer)
	starting_rank = Field(Integer)
	no_of_tests_passed = Field(Integer)
	no_of_tests = Field(Integer)
	max_score = Field(Float)	#2008-10-22 keep record of max/min score in among these snps
	min_score = Field(Float)
	type = ManyToOne('CandidateGeneTopSNPTestRMType', colname='type_id', ondelete='CASCADE', onupdate='CASCADE')
	null_data_ls = OneToMany('TopSNPTestRMNullData')
	comment = Field(Text)
	using_options(tablename='candidate_gene_top_snp_test_rm', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('results_id', 'list_type_id', 'no_of_top_snps', \
										'min_distance', 'starting_rank', 'type_id'))

class CandidateGeneTopSNPTestRG(Entity):
	"""
	2008-11-11
		a table similar to CandidateGeneTopSNPTestRM. but candidate_sample_size, non_candidate_sample_size and etc refer to distinct genes instead of not SNPs.
		
	"""
	result = ManyToOne('ResultsMethod', colname='results_id', ondelete='CASCADE', onupdate='CASCADE')
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE')
	pvalue = Field(Float)
	pvalue_gw_looping = Field(Float)
	pvalue_random_gene_list = Field(Float)
	candidate_sample_size = Field(Integer)
	non_candidate_sample_size = Field(Integer)
	candidate_gw_size = Field(Integer)
	non_candidate_gw_size = Field(Integer)
	no_of_top_candidate_genes = Field(Integer)
	no_of_top_genes = Field(Integer)
	no_of_top_snps = Field(Integer)
	min_distance = Field(Integer)
	starting_rank = Field(Integer)
	no_of_tests_passed = Field(Integer)
	no_of_tests = Field(Integer)
	max_score = Field(Float)
	min_score = Field(Float)
	type = ManyToOne('CandidateGeneTopSNPTestRMType', colname='type_id', ondelete='CASCADE', onupdate='CASCADE')
	comment = Field(Text)
	using_options(tablename='candidate_gene_top_snp_test_rg', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('results_id', 'list_type_id', 'no_of_top_snps', \
										'min_distance', 'starting_rank', 'type_id'))



class CandidateVsNonRatioPlot(Entity):
	"""
	2008-11-04
		table to store candidate-ratio/non-candidate-ratio curve plots in database
	"""
	type = ManyToOne('CandidateGeneTopSNPTestRMType', colname='type_id', ondelete='CASCADE', onupdate='CASCADE')
	result = ManyToOne('ResultsMethod', colname='results_id', ondelete='CASCADE', onupdate='CASCADE')
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE')
	png_thumbnail = Field(LargeBinary(64217728), deferred=True)
	png_data = Field(LargeBinary(134217728), deferred=True)
	svg_data = Field(LargeBinary(length=134217728), deferred=True)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_vs_non_ratio_plot', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('type_id', 'results_id', 'list_type_id'))

class TopSNPTestRMNullData(Entity):
	"""
	2008-10-31
		table to store NULL data of CandidateGeneTopSNPTestRM from different simulations under different NULL distributions
	"""
	observed = ManyToOne('CandidateGeneTopSNPTestRM', colname='observed_id', ondelete='CASCADE', onupdate='CASCADE')
	candidate_sample_size = Field(Integer)
	candidate_gw_size = Field(Integer)
	run_no = Field(Integer)
	null_distribution_type = ManyToOne('NullDistributionType', colname='null_distribution_type_id', ondelete='CASCADE', onupdate='CASCADE')
	using_options(tablename='top_snp_test_rm_null_data', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('observed_id', 'run_no', 'null_distribution_type_id'))
	
class TransformationMethod(Entity):
	"""
	2010-08-18
	"""
	name = Field(String(30))
	description = Field(Text)
	formular = Field(String(100))
	function = Field(String(20))
	using_options(tablename='transformation_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	
class SNPRegionPlotType(Entity):
	"""
	2008-10-06
	"""
	short_name = Field(String(256), unique=True)
	description = Field(String(8192))
	snp_region_plot_ls = OneToMany('SNPRegionPlot')
	#ManyToMany('SNPRegionPlot', tablename='snp_region_plot2type')
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='snp_region_plot_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class SNPRegionPlot(Entity):
	"""
	2011-4-19
		change type of chromosome from int to string
	2009-5-1
		add column call_method_id
	2008-10-21
		rename img_data to png_data and add svg_data
	2008-10-16
		fix an error in the definitions of img_data and original_filename. they were swapped
	2008-10-06
		table to store binary SNP region plots
	"""
	chromosome = Field(String(256))
	start = Field(Integer)
	stop = Field(Integer)
	png_data = Field(LargeBinary(134217728), deferred=True)
	svg_data = Field(LargeBinary(length=134217728), deferred=True)
	center_snp_position = Field(Integer)
	original_filename = Field(Text)	#no bigger than 128M
	call_method = ManyToOne('CallMethod', colname='call_method_id', ondelete='CASCADE', onupdate='CASCADE')
	phenotype_method = ManyToOne('PhenotypeMethod', colname='phenotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	plot_type = ManyToOne('SNPRegionPlotType', colname='plot_type_id', ondelete='CASCADE', onupdate='CASCADE')
	"""
	plot_types = ManyToMany('SNPRegionPlotType', colname='plot_type_id', ondelete='CASCADE', onupdate='CASCADE', \
						tablename='snp_region_plot2type')	#local_side, remote_side
	"""
	plot2gene_ls = OneToMany("SNPRegionPlotToGene")
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='snp_region_plot', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('chromosome', 'start', 'stop', 'phenotype_method_id', 'plot_type_id'))

class SNPRegionPlotToGene(Entity):
	"""
	2008-10-06
	"""
	snp_region_plot = ManyToOne('SNPRegionPlot', colname='plot_id', ondelete='CASCADE', onupdate='CASCADE')
	gene_id = Field(Integer)
	using_options(tablename='snp_region_plot2gene', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('plot_id', 'gene_id'))

class LD(Entity):
	"""
	2008-10-15
		table to store Linkage Disequilibrium data (output of MpiLD.py)
	"""
	chr1 = Field(Integer)
	pos1 = Field(Integer)
	chr2 = Field(Integer)
	pos2 = Field(Integer)
	d = Field(Float)
	d_prime = Field(Float)
	r2 = Field(Float)
	snp1_maf = Field(Float)
	snp2_maf = Field(Float)
	no_of_pairs = Field(Integer)
	call_method_id = Field(Integer)
	#call_method = ManyToOne('CallMethod', colname='call_method_id', ondelete='CASCADE', onupdate='CASCADE')
	using_options(tablename='ld', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class NullDistributionType(Entity):
	"""
	2008-10-21
		a table to store definition for different types null distributions
	"""
	short_name = Field(String(256), unique=True)
	description = Field(String(8192))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='null_distribution_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class ScoreRankHistogramType(Entity):
	"""
	2008-10-30
		connect it to "ResultsGene" via ManyToMany()
	2008-10-21
		add null_distribution_type_id
	2008-10-15
		type of score/rank histograms in ScoreRankHistogram
	"""
	call_method = ManyToOne('CallMethod', colname='call_method_id', ondelete='CASCADE', onupdate='CASCADE')
	min_distance = Field(Integer)
	get_closest = Field(Integer)
	min_MAF = Field(Float)
	allow_two_sample_overlapping = Field(Integer)
	results_type = Field(Integer)
	null_distribution_type = ManyToOne('NullDistributionType', colname='null_distribution_type_id', ondelete='CASCADE', onupdate='CASCADE')
	score_rank_hist_ls = OneToMany('ScoreRankHistogram')
	results_gene_ls = ManyToMany("ResultsGene", tablename='results_gene2type', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='score_rank_histogram_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('call_method_id', 'min_distance', 'get_closest', 'min_MAF', \
										'allow_two_sample_overlapping', 'results_type', 'null_distribution_type_id'))

class ScoreRankHistogram(Entity):
	"""
	2008-10-15
		table to store score/rank histogram divided by candidate gene list. output of CheckCandidateGeneRank.py
	"""
	phenotype_method = ManyToOne('PhenotypeMethod', colname='phenotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE')
	hist_type = ManyToOne('ScoreRankHistogramType', colname='hist_type_id', ondelete='CASCADE', onupdate='CASCADE')
	score_hist = Field(LargeBinary(length=134217728), deferred=True)
	score_hist_svg = Field(LargeBinary(length=134217728), deferred=True)
	rank_hist = Field(LargeBinary(length=134217728), deferred=True)
	rank_hist_svg = Field(LargeBinary(length=134217728), deferred=True)
	original_filename = Field(Text)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='score_rank_histogram', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('phenotype_method_id', 'list_type_id', 'hist_type_id'))

class MAFVsScorePlot(Entity):
	"""
	2008-10-17
		table to store pvalue vs MAF plots
	"""
	results_method = ManyToOne('ResultsMethod', colname='results_method_id', ondelete='CASCADE', onupdate='CASCADE')
	png_data = Field(LargeBinary(134217728), deferred=True)
	svg_data = Field(LargeBinary(length=134217728), deferred=True)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='maf_vs_score_plot', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('results_method_id'))

class ResultsGene(Entity):
	"""
	2009-4-7
		add rank into the unique combo. rank could be different for the same SNP under different MAF scenarios.
	2008-10-30
		remove order_by='rank' in types = ManyToMany('ScoreRankHistogramType', ...)
	2008-10-28
		store genes associated with SNPs in results_method under certain association rule
	"""
	snp = ManyToOne('Snps', colname='snps_id', ondelete='CASCADE', onupdate='CASCADE')
	disp_pos = Field(Integer)
	gene_id = Field(Integer)
	score = Field(Float)
	rank = Field(Float)
	result = ManyToOne('ResultsMethod', colname='results_id', ondelete='CASCADE', onupdate='CASCADE')
	types = ManyToMany('ScoreRankHistogramType', \
						tablename='results_gene2type')	#local_side, remote_side
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='results_gene', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('snps_id', 'gene_id', 'results_id', 'rank'))


class FTPathway(Entity):
	"""
	2008-11-14
		table storing flowering time pathways
	"""
	short_name = Field(String(256), unique=True)
	description = Field(String(8192))
	gene_ls = OneToMany('FTGene')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='ft_pathway', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')


class FTPathwayRelationship(Entity):
	"""
	2008-11-14
		table storing flowering time pathway relationships
	"""
	pathway1 = ManyToOne('FTPathway', colname='pathway1_id', ondelete='CASCADE', onupdate='CASCADE')
	pathway2 = ManyToOne('FTPathway', colname='pathway2_id', ondelete='CASCADE', onupdate='CASCADE')
	relationship_type_id = Field(Integer)
	arrow_start_point_loc = Field(String(512))
	arrow_end_point_loc = Field(String(512))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='ft_pathway_relationship', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('pathway1_id', 'pathway2_id', 'relationship_type_id'))


class FTGene(Entity):
	"""
	2008-11-14
		table storing flowering time pathway genes
	"""
	gene_id = Field(Integer)
	pathway = ManyToOne('FTPathway', colname='pathway_id', ondelete='CASCADE', onupdate='CASCADE')
	paper = Field(String(512))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='ft_gene', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('gene_id', 'pathway_id'))


class SNPAnnotation(Entity):
	"""
	2009-01-05
		table to store annotation of SNPs, like synonymous, non-synonymous, ...
		information finer than SnpsContext
	"""
	snp = ManyToOne('Snps', colname='snps_id', ondelete='CASCADE', onupdate='CASCADE')
	gene_id = Field(Integer)
	gene_commentary_id = Field(Integer)
	snp_annotation_type = ManyToOne('SNPAnnotationType', colname='snp_annotation_type_id', ondelete='CASCADE', onupdate='CASCADE')
	which_exon_or_intron = Field(Integer)
	pos_within_codon = Field(Integer)
	comment = Field(String(512))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='snp_annotation', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('snps_id', 'gene_commentary_id', 'snp_annotation_type_id'))

class SNPAnnotationType(Entity):
	"""
	2009-01-05
		table to store types of SNP annotation, like synonymous, non-synonymous, ...
	"""
	short_name = Field(String(256), unique=True)
	description = Field(String(8192))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='snp_annotation_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')


class CallInfoPCValue(Entity):
	"""
	2009-2-2
		table to store principle component values after PCA is applied on a certain call method dataset
	"""
	call_info = ManyToOne('CallInfo', colname='call_info_id', ondelete='CASCADE', onupdate='CASCADE')
	which_pc = Field(Integer)
	pc_value = Field(Float)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='call_info_pc_value', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('call_info_id', 'which_pc'))


class CallMethodEigenValue(Entity):
	"""
	2009-2-2
		table to store eigen values after PCA is applied on a certain call method dataset
	"""
	call_method = ManyToOne('CallMethod', colname='call_method_id', ondelete='CASCADE', onupdate='CASCADE')
	which_eigen = Field(Integer)
	eigen_value = Field(Float)
	variance_perc = Field(Float)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='call_method_eigen_value', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('call_method_id', 'which_eigen'))

class PhenotypeHistPlot(Entity):
	"""
	2009-2-6
		phenotype histogram
	"""
	phenotype_method = ManyToOne('PhenotypeMethod', colname='phenotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	hist_thumb = Field(LargeBinary(64217728), deferred=True)
	hist_plot = Field(LargeBinary(134217728), deferred=True)
	hist_plot_fname = Field(String(1000))
	hist_log_thumb = Field(LargeBinary(64217728), deferred=True)
	hist_log_plot = Field(LargeBinary(134217728), deferred=True)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='phenotype_hist_plot', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('phenotype_method_id'))

class CallPhenotypeQQPlot(Entity):
	"""
	2009-2-6
		qq plot
	"""
	call_method = ManyToOne('CallMethod', colname='call_method_id', ondelete='CASCADE', onupdate='CASCADE')
	phenotype_method = ManyToOne('PhenotypeMethod', colname='phenotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	qq_thumb = Field(LargeBinary(64217728), deferred=True)
	qq_plot = Field(LargeBinary(134217728), deferred=True)
	qq_plot_fname = Field(String(1000))
	qq_log_thumb = Field(LargeBinary(64217728), deferred=True)
	qq_log_plot = Field(LargeBinary(134217728), deferred=True)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='call_phenotype_qq_plot', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('call_method_id', 'phenotype_method_id'))

class CmpEnrichmentOfTwoAnalysisMethods(Entity):
	"""
	2009-10-2
		add column pvalue, for the enrichment_ratio under certain null distribution based on type_id
	2009-4-16 to be filled up by PlotCmpTwoAnalysisMethods.py
	"""
	result1 = ManyToOne('ResultsMethod', colname='results_id1', ondelete='CASCADE', onupdate='CASCADE')
	result2 = ManyToOne('ResultsMethod', colname='results_id2', ondelete='CASCADE', onupdate='CASCADE')
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE')
	type = ManyToOne('ScoreRankHistogramType', colname='type_id', ondelete='CASCADE', onupdate='CASCADE')
	r1_min_score = Field(Float)
	r1_max_score = Field(Float)
	r2_min_score = Field(Float)
	r2_max_score = Field(Float)
	threshold_type = ManyToOne('CmpEnrichmentOfTwoAnalysisMethodsType', colname='threshold_type_id', ondelete='CASCADE', onupdate='CASCADE')
	candidate_sample_size = Field(Integer)
	non_candidate_sample_size = Field(Integer)
	candidate_gw_size = Field(Integer)
	non_candidate_gw_size = Field(Integer)
	enrichment_ratio = Field(Float)
	pvalue = Field(Float)	# 2009-10-2 pvalue for the enrichment_ratio under certain null distribution based on type_id
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='cmp_enrichment_of_two_analysis_methods', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('results_id1', 'results_id2', 'list_type_id', 'type_id', 'r1_min_score', \
										'r1_max_score', 'r2_min_score', 'r2_max_score', 'threshold_type_id'))

class CmpEnrichmentOfTwoAnalysisMethodsType(Entity):
	"""
	2009-5-4 table to store threshold types of CmpEnrichmentOfTwoAnalysisMethods
	"""
	r1_threshold = Field(Float)
	r2_threshold = Field(Float)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='cmp_enrichment_of_two_analysis_methods_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('r1_threshold', 'r2_threshold'))


class CNVQCAccession(Entity):
	"""
	2009-10-26
		to keep track of where the accessions are from
	"""
	original_id = Field(String(200))
	alternative_id = Field(String(2000))
	ecotype_id = Field(Integer)
	data_source = ManyToOne('DataSource', colname='data_source_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='cnv_qc_accession', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('original_id', 'data_source_id'))

class DataSource(Entity):
	"""
	2009-10-26
		records of where the data is from
	"""
	short_name = Field(String(200), unique=True)
	url = Field(String(2000))
	address = Field(String(2000))
	contact = Field(String(2000))
	email = Field(String(2000))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='data_source', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class CNV(Entity):
	"""
	2011-4-19
		change type of chromosome from int to string
	2010-7-28
		table storing all types of copy number variation
	"""
	chromosome = Field(String(256))
	start = Field(Integer)
	stop = Field(Integer)
	start_probe = ManyToOne("%s.Probes"%__name__, colname='start_probe_id', ondelete='CASCADE', onupdate='CASCADE')
	stop_probe = ManyToOne("%s.Probes"%__name__, colname='stop_probe_id', ondelete='CASCADE', onupdate='CASCADE')
	no_of_probes_covered = Field(Integer)
	size_affected = Field(Integer)
	score = Field(Float)
	score_std = Field(Float)
	frequency = Field(Float)
	fractionNotCoveredByLyrata = Field(Float)
	comment = Field(String(8124), deferred=True)
	cnv_type =  ManyToOne('%s.CNVType'%__name__, colname='cnv_type_id', ondelete='CASCADE', onupdate='CASCADE')
	cnv_method = ManyToOne('%s.CNVMethod'%__name__, colname='cnv_method_id', ondelete='CASCADE', onupdate='CASCADE')
	cnv_call_ls = ManyToMany("CNVCall", tablename='cnv2cnv_call', ondelete='CASCADE', onupdate='CASCADE')	#2010-7-31
	cnv_array_call_ls = OneToMany("CNVArrayCall")
	created_by = Field(String(200), deferred=True)
	updated_by = Field(String(200), deferred=True)
	date_created = Field(DateTime, default=datetime.now, deferred=True)
	date_updated = Field(DateTime, deferred=True)
	using_options(tablename='cnv', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('chromosome', 'start', 'stop', 'cnv_type_id', 'cnv_method_id'))

class CNVArrayCall(Entity):
	"""
	2010-7-28
		
	"""
	array = ManyToOne("%s.ArrayInfo"%__name__, colname='array_id', ondelete='CASCADE', onupdate='CASCADE')
	cnv = ManyToOne("%s.CNV"%__name__, colname='cnv_id', ondelete='CASCADE', onupdate='CASCADE')
	cnv_method = ManyToOne('%s.CNVMethod'%__name__, colname='cnv_method_id', ondelete='CASCADE', onupdate='CASCADE')
	score = Field(Float)
	fractionDeletedInPECoverageData = Field(Float)	#2010-8-1
	comment = Field(String(8124), deferred=True)
	created_by = Field(String(200), deferred=True)
	updated_by = Field(String(200), deferred=True)
	date_created = Field(DateTime, default=datetime.now, deferred=True)
	date_updated = Field(DateTime, deferred=True)
	using_options(tablename='cnv_array_call', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('array_id', 'cnv_id', 'cnv_method_id'))

class CNVMethod(Entity):
	"""
	2012.2.28
		add column no_of_accessions, no_of_loci
	2011-2-17
		add column filename
	2009-10-26
		which type of method is used to infer CNVs
	"""
	short_name = Field(String(200), unique=True)
	no_of_accessions = Field(Integer)	#2012.2.28
	no_of_loci = Field(Integer)	#2012.2.28
	description = Field(String(8000))
	filename = Field(Text)
	comment = Field(String(8000))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='cnv_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class CNVType(Entity):
	"""
	2009-10-26
		type of CNV, insertion, deletion, duplication
	"""
	short_name = Field(String(200), unique=True)
	description = Field(String(8000))
	comment = Field(String(8000))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='cnv_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class CNVQCCall(Entity):
	"""
	2010-6-24
		add chromosome_stop, and change the unique constraint
	2009-12-8
		add copy_number to distinguish CNVs that have multiple copies
	2009-10-26
		the CNV data collected from other sources to verify our calls
	"""
	accession = ManyToOne('CNVQCAccession', colname='accession_id', ondelete='CASCADE', onupdate='CASCADE')
	chr = ManyToOne('Chromosome', colname='chromosome', ondelete='CASCADE', onupdate='CASCADE')	#chr must be different from field name.
	start = Field(Integer)
	stopChromosome = ManyToOne('Chromosome', colname='stop_chromosome', ondelete='CASCADE', onupdate='CASCADE')
	stop = Field(Integer)
	size_affected = Field(Integer)
	score = Field(Float)
	no_of_probes_covered = Field(Integer)
	copy_number = Field(Integer)
	sequence = Field(LargeBinary(length=134217728), deferred=True)
	comment = Field(String(8124))
	cnv_type = ManyToOne('CNVType', colname='cnv_type_id', ondelete='CASCADE', onupdate='CASCADE')
	cnv_method = ManyToOne('CNVMethod', colname='cnv_method_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='cnv_qc_call', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('accession_id', 'chromosome', 'start', 'stop_chromosome', 'stop',\
										'cnv_type_id', 'cnv_method_id'))

class CNVQCProbeCall(Entity):
	"""
	2011-4-19
		change type of chromosome from int to string
	2009-10-26
		the probe-based CNV QC data, from
			1. check the probes against the data from CNVQCCall
			2. blast the probes against some other sequence data, 2010 or Ler contigs
	"""
	accession = ManyToOne('CNVQCAccession', colname='accession_id', ondelete='CASCADE', onupdate='CASCADE')
	probe = ManyToOne("Probes", colname='probe_id', ondelete='CASCADE', onupdate='CASCADE')
	chromosome = Field(String(256))
	position = Field(Integer)
	size_affected = Field(Integer)
	target_position = Field(Integer)
	score = Field(Float)
	comment = Field(String(8124))
	cnv_qc_call = ManyToOne("CNVQCCall", colname='cnv_qc_call_id', ondelete='CASCADE', onupdate='CASCADE')
	cnv_type = ManyToOne('CNVType', colname='cnv_type_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='cnv_qc_probe_call', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('accession_id', 'probe_id', 'target_position', 'cnv_type_id', 'cnv_qc_call_id'))

class CNVCall(Entity):
	"""
	2011-4-19
		change type of chromosome from int to string
	2010-7-1
		add column probability
	2010-6-29
		add column median_intensity
	2010-3-16
		call_info_id becomes array_id
		add column no_of_probes_covered
	2009-10-26
		the CNV from the tiling part of the 250k affy array
	"""
	array = ManyToOne("%s.ArrayInfo"%__name__, colname='array_id', ondelete='CASCADE', onupdate='CASCADE')
	chromosome = Field(String(256))
	start = Field(Integer)
	stop = Field(Integer)
	start_probe = ManyToOne("%s.Probes"%__name__, colname='start_probe_id', ondelete='CASCADE', onupdate='CASCADE')
	stop_probe = ManyToOne("%s.Probes"%__name__, colname='stop_probe_id', ondelete='CASCADE', onupdate='CASCADE')
	no_of_probes_covered = Field(Integer)
	size_affected = Field(Integer)
	amplitude = Field(Float)
	median_intensity = Field(Float)
	score = Field(Float)
	probability = Field(Float)	# 2010-7-1
	percUnCoveredByLerContig = Field(Float)	#2010-6-27
	fractionDeletedInPECoverageData = Field(Float)	#2010-7-23
	comment = Field(String(8124))
	cnv_type =  ManyToOne('%s.CNVType'%__name__, colname='cnv_type_id', ondelete='CASCADE', onupdate='CASCADE')
	cnv_method = ManyToOne('%s.CNVMethod'%__name__, colname='cnv_method_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='cnv_call', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('array_id', 'start_probe_id', 'stop_probe_id', 'cnv_type_id', 'cnv_method_id'))

class CNVContext(Entity):
	"""
	2010-8-18
		update on a few columns
	2010-7-31
		context for CNVs
	"""
	cnv = ManyToOne('CNV', colname='cnv_id', ondelete='CASCADE', onupdate='CASCADE')
	term5_disp_pos = Field(Float)	#[0,1) is for within gene fraction, <=-1 is upstream. >=1 is downstream
	term3_disp_pos = Field(Float)	#(0,1] is for within gene fraction, null if no overlap/outside.
	gene_id = Field(Integer)
	gene_strand = Field(String(1))
	overlap_length = Field(Integer)
	overlap_fraction_in_gene = Field(Float)
	overlap_fraction_in_cnv = Field(Float)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='cnv_context', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('cnv_id', 'gene_id'))

class CNVAnnotation(Entity):
	"""
	2010-8-18
		store details of CNVContext, which exon, which intron. which CDS, which UTR, etc.
	"""
	cnv = ManyToOne('CNV', colname='cnv_id', ondelete='CASCADE', onupdate='CASCADE')
	cnv_context = ManyToOne('CNVContext', colname='cnv_context_id', ondelete='CASCADE', onupdate='CASCADE')
	gene_commentary_id = Field(Integer)
	gene_segment_id = Field(Integer)
	label = Field(Text)
	utr_number = Field(Integer)
	cds_number = Field(Integer)
	intron_number = Field(Integer)
	exon_number = Field(Integer)
	overlap_length = Field(Integer)
	overlap_fraction_in_gene = Field(Float)
	overlap_fraction_in_cnv = Field(Float)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='cnv_annotation', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('cnv_id', 'cnv_context_id', 'gene_commentary_id', 'gene_segment_id'))


class SequenceFragment(Entity):
	"""
	2010-1-27
		table that keeps track of sequence fragments such as Ler contigs
	"""
	accession = ManyToOne('CNVQCAccession', colname='accession_id', ondelete='CASCADE', onupdate='CASCADE')
	short_name = Field(String(256))
	size = Field(Integer)
	description = Field(String(8124))
	comment = Field(String(8124))
	sequence = Field(LargeBinary(length=134217728), deferred=True)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='sequence_fragment', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('accession_id', 'short_name'))


class SequenceFragmentRefPos(Entity):
	"""
	2011-4-19
		change type of chromosome from int to string
	2010-12-6
		change field copy_number to cnv_type
	2010-6-14
		change the unqiue constraint
		add column version to differentiate between different version of SequenceFragmentRefPos data
	2010-1-27
		the position of sequence fragments mapped onto the reference genome (which is Col)
	"""
	sequence_fragment = ManyToOne('%s.SequenceFragment'%__name__, colname='sequence_fragment_id', ondelete='CASCADE', onupdate='CASCADE')
	chromosome = Field(String(256))
	start = Field(Integer)
	stop = Field(Integer)
	size_difference = Field(Integer)
	start_probe = ManyToOne("%s.Probes"%__name__, colname='start_probe_id', ondelete='CASCADE', onupdate='CASCADE')
	stop_probe = ManyToOne("%s.Probes"%__name__, colname='stop_probe_id', ondelete='CASCADE', onupdate='CASCADE')
	no_of_probes_covered = Field(Integer)
	fragment_start = Field(Integer)
	fragment_stop = Field(Integer)
	cnv_type = ManyToOne('CNVType', colname='cnv_type_id', ondelete='CASCADE', onupdate='CASCADE')
	version = Field(Integer)
	comment = Field(String(8124))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='sequence_fragment_ref_pos', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('sequence_fragment_id', 'start', 'stop', 'fragment_start', 'fragment_stop', 'version'))


class SequenceFragment2Probe(Entity):
	"""
	2011-4-19
		change type of chromosome from int to string
	2010-4-15
		table recording the matching probes of sequence fragments
	"""
	sequence_fragment = ManyToOne('%s.SequenceFragment'%__name__, colname='sequence_fragment_id', ondelete='CASCADE', onupdate='CASCADE')
	fragment_start = Field(Integer)
	fragment_stop = Field(Integer)
	probe = ManyToOne("%s.Probes"%__name__, colname='probe_id', ondelete='CASCADE', onupdate='CASCADE')
	chromosome = Field(String(256))
	start = Field(Integer)
	stop = Field(Integer)
	no_of_identities = Field(Integer)
	comment = Field(String(8124))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='sequence_fragment2probe', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('sequence_fragment_id', 'probe_id', 'fragment_start', 'fragment_stop'))



class AssociationOverlappingType(Entity):
	"""
	2009-11-30
		add column no_of_methods = number of association methods involved in this overlapping type
	2009-11-2
		
	"""
	short_name = Field(String(200), unique=True)
	no_of_methods = Field(Integer)
	description = Field(String(8000))
	comment = Field(String(8000))
	analysis_method_ls = ManyToMany("AnalysisMethod", tablename='overlapping_type2method', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='association_overlapping_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class AssociationOverlappingStat(Entity):
	"""
	2009-11-2
	"""
	call_method = ManyToOne('CallMethod', colname='call_method_id', ondelete='CASCADE', onupdate='CASCADE')
	phenotype_method = ManyToOne('PhenotypeMethod', colname='phenotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	no_of_top_snps = Field(Integer)
	no_of_overlapping_snps = Field(Integer)
	overlapping_type = ManyToOne('AssociationOverlappingType', colname='overlapping_type_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='association_overlapping_stat', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('call_method_id', 'phenotype_method_id', 'overlapping_type_id', 'no_of_top_snps'))

"""
#2008-10-29 'results_gene2type' automatically generated by ManyToMany is sub-optimal because it uses MyIAM engine and doesn't allow foreign key.
#two ways around this
#1. set "default-storage-engine=InnoDB" in mysqld section in my.cnf
#OR
#2. create the table manually with same columns + innodb engine + foreign keys 
class ResultsGene2Type(Entity):
	results_gene = ManyToOne('ResultsGene', colname='results_gene_id', ondelete='CASCADE', onupdate='CASCADE')
	score_rank_histogram_type = ManyToOne('ScoreRankHistogramType', colname='score_rank_histogram_type_id', ondelete='CASCADE', onupdate='CASCADE')
	using_options(tablename='results_gene2type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
"""

"""
class TopSNPTestRMPlotType(Entity):
	#2008-10-20
	#	type for TopSNPTestRMPlot
	no_of_snps = Field(Integer)
	test_type = Field(Integer)
	results_type = Field(Integer)
	min_distance = Field(Integer)
	get_closest = Field(Integer)
	min_MAF = Field(Float)
	
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='top_snp_test_rm_plot_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('no_of_snps', 'test_type', 'results_type', \
										'min_distance', 'get_closest', 'min_MAF'))

class TopSNPTestRMPlot(Entity):
	#2008-10-20
	#	store plots based on data from CandidateGeneTopSNPTestRM
	
	png_data = Field(Binary(134217728), deferred=True)
	svg_data = Field(Binary(length=134217728), deferred=True)
	results_method = ManyToOne('ResultsMethod', colname='results_method_id', ondelete='CASCADE', onupdate='CASCADE')
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE')
	plot_type = ManyToOne('TopSNPTestRMPlotType', colname='plot_type_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='top_snp_test_rm_plot', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('results_method_id', 'list_type_id', 'plot_type_id'))
"""

class WebHelp(Entity):
	"""
	2009-4-26 table to store help document of the web interface (accessible on the web thru a controller)
	"""
	short_name = Field(String(256), unique=True)
	title = Field(String(1024))
	content = Field(Text(8192*20))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='web_help', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class QCCrossMatch(Entity):
	"""
	2009-6-9
		table recording the cross-matching results, used to identify mis-labeling, contaminant, etc.
	"""
	ecotype_id = Field(Integer)
	call_info = ManyToOne("CallInfo", colname='call_info_id', ondelete='CASCADE', onupdate='CASCADE')
	vs_ecotype_id = Field(Integer)
	mismatch_rate = Field(Float)
	no_of_mismatches = Field(Integer)
	no_of_non_NA_pairs = Field(Integer)
	qc_method = ManyToOne("%s.QCMethod"%__name__, colname='qc_method_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='qc_cross_match', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('ecotype_id', 'call_info_id', 'vs_ecotype_id', 'qc_method_id'))
	
	
class Groups(Entity):
	id = Field(Integer,primary_key=True)
	name = Field(Unicode(30),required=True)
	users = ManyToMany('Users',tablename='acl_groupmembers',local_colname='groups_id')
	phenotype_methods = ManyToMany("PhenotypeMethod",tablename='acl_phenotype_method_groups')
	using_table_options(mysql_engine='InnoDB',useexisting=True)
	using_options(tablename='acl_groups',metadata=__metadata__, session=__session__)
	def __repr__(self):
		return (u'<Group: name=%s>' % self.name).encode('utf-8')


class Users(Entity):
	id = Field(Integer,primary_key=True)
	email = Field(String(100), required=True)
	_password = Field(String(40), required = True,colname='password',synonym='password')
	realname = Field(Unicode(50), required = True)
	organisation = Field(Unicode(100))
	isAdmin = Field(Enum(('Y','N')),default='N',required=True )
	created = Field(DateTime, default=datetime.now)
	groups = ManyToMany('Groups',tablename='acl_groupmembers',local_colname='users_id')
	phenotype_methods = ManyToMany("PhenotypeMethod",tablename='acl_phenotype_method_users')
	using_options(tablename='acl_users',metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	
	def validate_password(self, password):
		"""
        Check the password against existing credentials.
        
        :param password: the password that was provided by the user to
            try and authenticate. This is the clear text version that we will
            need to match against the hashed one in the database.
        :type password: unicode object.
        :return: Whether the password is valid.
        :rtype: bool
        
        """
		hashed_pass = hashlib.sha1()
		hashed_pass.update(self._password[:8] + password)
		return self._password[8:] == hashed_pass.hexdigest()

	def _set_password(self, password):
		"""encrypts password on the fly using the encryption
        algo defined in the configuration
        """
		self._password = self._encrypt_password(password)

	def _get_password(self):
		"""returns password
        """
		return self._password

	password = property(_get_password,_set_password)

	@classmethod
	def _encrypt_password(cls,  password):
		"""Hash the given password with the specified algorithm. Valid values
        for algorithm are 'md5' and 'sha1'. All other algorithm values will
        be essentially a no-op."""
		hashed_password = password

		if isinstance(password, unicode):
			password_8bit = password.encode('UTF-8')

		else:
			password_8bit = password


		salt = hashlib.sha1()
		salt.update(os.urandom(60))
		salt_text = salt.hexdigest()
		hash = hashlib.sha1()
		hash.update(salt_text[:8] + password_8bit)
		hashed_password = salt_text[:8] + hash.hexdigest()
		print '*'*20,  salt_text[:8], " ", hashed_password[:8]




		if not isinstance(hashed_password, unicode):
			hashed_password = hashed_password.decode('UTF-8')

		return hashed_password

class Chromosome(Entity):
	"""
	2011-4-25
		add column rank (integer) to order chromosomes
	2010-6-24
		a dedicated table to translate between chromosome names and numbers
	"""
	name = Field(String(255))
	description = Field(String(6000))
	tax_id = Field(Integer)
	rank = Field(Integer)	#2011-4-19, within a genome, rank=1 is the first chromosome and so on.
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_table_options(UniqueConstraint('name', 'tax_id'))
	using_options(tablename='chromosome', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class Stock_250kDB(ElixirDB):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ):['localhost', 'z', 1, 'hostname of the db server', ],\
							('database', 1, ):['stock_250k', 'd', 1, '',],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('username', 1, ):[None, 'u', 1, 'database username',],\
							('password', 1, ):[None, 'p', 1, 'database password', ],\
							('port', 0, ):[None, 'o', 1, 'database port number'],\
							('pool_recycle', 0, int):[3600, 'l', 1, 'the length of time to keep connections open before recycling them.'],\
							('sql_echo', 0, bool):[False, 's', 0, 'display SQL Statements'],\
							('echo_pool', 0, bool):[False, 'e', 0, 'if True, the connection pool will log all checkouts/checkins to the logging stream, which defaults to sys.stdout.'],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-10-06
			add option 'pool_recycle' to recycle connection. MySQL typically close connections after 8 hours.
			__metadata__.bind = create_engine(self._url, pool_recycle=self.pool_recycle)
		2008-07-09
		"""
		from pymodule import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
												class_to_have_attr=self)
		
		if self.echo_pool:	#2010-9-19 passing echo_pool to create_engine() causes error. all pool log disappeared.
			#2010-9-19 Set up a specific logger with our desired output level
			import logging
			#import logging.handlers
			logging.basicConfig()
			#LOG_FILENAME = '/tmp/sqlalchemy_pool_log.out'
			my_logger = logging.getLogger('sqlalchemy.pool')
	
			# Add the log message handler to the logger
			#handler = logging.handlers.RotatingFileHandler(LOG_FILENAME, maxBytes=1000000, backupCount=5)
			# create formatter
			#formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
			# add formatter to handler
			#handler.setFormatter(formatter)
			#my_logger.addHandler(handler)
			my_logger.setLevel(logging.DEBUG)
		self.setup_engine(metadata=__metadata__, session=__session__, entities=entities)
		
		
		self._cnv_id2chr_pos = {}	#2011-2-24
		self._cnv_method_id = None	#2011-2-24
		self._chr_pos2snp_id = {}	#2011-2-28
		self._snp_id2chr_pos = {}	#2011-2-28
	
	def setup(self, create_tables=True):
		"""
		2008-09-07
			expose option create_tables, default=True. assign it to False if no new table is to be created.
		"""
		setup_all(create_tables=create_tables)	#create_tables=True causes setup_all to call elixir.create_all(), which in turn calls metadata.create_all()
		#2008-08-26 setup_all() would setup other databases as well if they also appear in the program. Seperate this to be envoked after initialization
		# to ensure the metadata of other databases is setup properly.
	
	@property
	def data_dir(self, ):
		"""
		2012.3.23
			(learnt from VervetDB)
			get the master directory in which all files attached to this db are stored.
		"""
		dataDirEntry = README.query.filter_by(title='data_dir').first()
		if not dataDirEntry or not dataDirEntry.description or not os.path.isdir(dataDirEntry.description):
			# todo: need to test dataDirEntry.description is writable to the user
			sys.stderr.write("data_dir not available in db or not accessible on the harddisk. Raise exception.\n")
			raise
			return None
		else:
			return dataDirEntry.description
	
	def reScalePathByNewDataDir(self, filePath="", newDataDir=""):
		"""
		2012.3.23
			in case that the whole /Network/Data/250k/db is stored in a different place (=newDataDir)
				how to rescale the filePath ( stored in the database tables) to reflect its new path.
		"""
		if not newDataDir:	#newDataDir is nothing.
			return filePath
		else:
			oldDataDir = self.data_dir()
			if filePath.find(oldDataDir)==0:
				relativePath = filePath[len(oldDataDir)-1:]
				newPath = os.path.join(newDataDir, relativePath)
				return newPath
			else:
				sys.stderr.write("Warning: %s doesn't include old data dir %s. Return Nothing.\n"%(filePath, oldDataDir))
				return None
	
	#@classmethod
	def getGWA(self, call_method_id, phenotype_method_id, analysis_method_id, results_directory=None, min_MAF=0.1, \
			construct_chr_pos2index=False, pdata=None):
		"""
		2012.3.23
			not a classmethod anymore
		2010-3-14
			become a classmethod
		2009-3-12
			convenient function to get genome-wide-association result from database
		"""
		rm = ResultsMethod.query.filter_by(call_method_id=call_method_id).filter_by(phenotype_method_id=phenotype_method_id).\
					filter_by(analysis_method_id=analysis_method_id).first()
		#from GeneListRankTest import GeneListRankTest
		return self.getResultMethodContent(rm.id, results_directory=results_directory, min_MAF=min_MAF, 
																construct_chr_pos2index=construct_chr_pos2index,\
																pdata=pdata)
	
	@classmethod
	def getSNPMatrix(cls, call_method_id, ignore_2nd_column=True):
		"""
		2010-3-14
			become a classmethod
		2009-3-12
			given a call_method_id, fetch a whole SNP dataset
		"""
		cm = CallMethod.get(call_method_id)
		snpData = SNPData(input_fname=cm.filename, turn_into_array=1, ignore_2nd_column=ignore_2nd_column)	#use 1st column (ecotype id) as main ID
		return snpData
	
	def getSNP(self, chromosome=None, start=None, stop=None, locus_type_id=None):
		"""
		2012.3.20
			add locus_type_id
		2011-4-22
			get Snps object based on chromosome, start, stop or else create it if not in db
		"""
		query = Snps.query.filter_by(chromosome=chromosome).filter_by(position=start)
		if stop:
			query = query.filter(Snps.end_position==stop)
		db_entry = query.first()
		no_of_entries = query.count()
		if no_of_entries>1:
			sys.stderr.write("Warning: %s db entries for (chr=%s, start=%s, stop=%s).\n"%(no_of_entries, chromosome, start, stop))
		if not db_entry:
			name = '%s_%s'%(chromosome, start)
			if stop:
				name += '_%s'%(stop)
			db_entry = Snps(name=name, chromosome=chromosome, position=start, \
							end_position=stop, locus_type_id=locus_type_id)
			self.session.add(db_entry)
		return db_entry
	
	def dealWithSNPInfo(self, snpInfoPickleFname, locus_type_id=1):
		"""
		2012.3.19 a wrapper around getSNPInfo() to pickle the data structure
		"""
		sys.stderr.write("Dealing with SNPInfo ...")
		from pymodule.SNP import SNPInfo
		from pymodule.CNV import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
		from pymodule.RBTree import RBDict
		import cPickle
		if snpInfoPickleFname:
			if os.path.isfile(snpInfoPickleFname):	#if this file is already there, suggest to un-pickle it.
				picklef = open(snpInfoPickleFname)
				snpInfo = cPickle.load(picklef)
				del picklef
			else:	#if the file doesn't exist, but the filename is given, pickle data into it
				snpInfo = self.getSNPInfo(locus_type_id=locus_type_id)
				#pickle the object
				picklef = open(snpInfoPickleFname, 'w')
				cPickle.dump(snpInfo, picklef, -1)
				picklef.close()
		else:
			snpInfo = self.getSNPInfo(locus_type_id=locus_type_id)
		sys.stderr.write(" %s snps, %s distinct loci in locusRBDict. %s loci have the same genomic coordinates.\n"%\
						(len(snpInfo.snps_id2index), len(snpInfo.locusRBDict), snpInfo.no_of_same_position_loci))
		return snpInfo
	
	def getSNPInfo(self, locus_type_id=1):
		"""
		2012.3.19
			moved from variation.src.DrawSNPRegion.DrawSNPRegion
		2012.3.7
			add end_position into chr_pos
			todo :
				#. have to somehow separate snp-gwas loci and cnv-gwas loci , because they overlap
				#. generate RB dictionary to store all loci for fast retrieval based on (chr, start, stop)
		2009-2-18
			replace PassingData with class SNPInfo from pymodule.SNP to wrap up all data
		2009-1-22
			order by chromosome, position
			no chr_pos2snps_id, can get snps_id from chr_pos2index => data_ls
		2009-1-5
			add chr_pos2snps_id
		2008-09-24
			in order
		"""
		sys.stderr.write("Getting info of locus_type_id=%s loci in chromosomal order ...  "%(locus_type_id))
		chr_pos_ls = []
		data_ls = []
		chr_pos2index = {}
		snps_id2index = {}
		i = 0
		block_size = 50000
		
		from pymodule.SNP import SNPInfo
		from pymodule.CNV import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
		from pymodule.RBTree import RBDict
		locusRBDict = RBDict()
		
		where_ls = ['chromosome is not null', 'position is not null']
		if locus_type_id:
			where_ls.append('locus_type_id=%s'%(locus_type_id))
		
		rows = self.metadata.bind.execute("select id, chromosome, position, end_position, allele1, allele2, locus_type_id from %s \
				where  %s order by chromosome, position, end_position"%(Snps.table.name, ' and '.join(where_ls)))
		#.query.offset(i).limit(block_size)
		#while rows.count()!=0:
		prev_outputted_i = ""
		no_of_same_position_loci = 0
		for row in rows:
			obj = PassingData()
			for item in row.items():
				setattr(obj, item[0], item[1])
			if obj.end_position is None:
				obj.end_position = obj.position	#this operation couldn't happen on row (a RowProxy object). 
				#*** AttributeError: 'RowProxy' object has no attribute 'end_position'
			chr_pos = (obj.chromosome, obj.position, obj.end_position)
			chr_pos_ls.append(chr_pos)
			chr_pos2index[chr_pos] = len(chr_pos2index)
			snps_id2index[obj.id] = len(snps_id2index)
			data_ls.append(obj)
			
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=obj.chromosome, \
							span_ls=[obj.position, obj.end_position], \
							min_reciprocal_overlap=1,)
							#2010-8-17 overlapping keys are regarded as separate instances as long as they are not identical.
			if segmentKey not in locusRBDict:
				locusRBDict[segmentKey] = []
			else:
				no_of_same_position_loci += 1
			locusRBDict[segmentKey].append(snps_id2index[obj.id])	#value is a list of the index
			i += 1
		#	if self.debug and i>40000:
		#		break
		#	rows = Stock_250kDB.Snps.query.offset(i).limit(block_size)
			if i%10000==0:
				sys.stderr.write("%s%s"%('\x08'*len(str(prev_outputted_i)), i))
				prev_outputted_i = i
		sys.stderr.write("%s%s"%('\x08'*len(str(prev_outputted_i)), i))
		snpInfo = SNPInfo(locus_type_id=locus_type_id)
		#snpInfo.chr_pos_ls = chr_pos_ls
		snpInfo.chr_pos2index = chr_pos2index
		snpInfo.snps_id2index = snps_id2index
		snpInfo.data_ls = data_ls
		snpInfo.locusRBDict = locusRBDict
		snpInfo.no_of_same_position_loci = no_of_same_position_loci
		sys.stderr.write(" %s snps, %s distinct loci in locusRBDict. %s loci have the same genomic coordinates.\n"%\
						(len(snpInfo.snps_id2index), len(snpInfo.locusRBDict), snpInfo.no_of_same_position_loci))
		return snpInfo
	
	def getNonCNVSnpsID2CNVSnpsID(self, cnv_locus_type_id=2):
		"""
		2012-3.20
			create a dictionary between locus_type_id=1 and locus_type_id=2 on the condition value of (chromosome, postion) pairs.
		"""
		sys.stderr.write("Getting a map between non-CNV Snps.id (locus_type_id!=%s) and cnv-like Snps.id (locus_type_id=%s) ..."%\
						(cnv_locus_type_id, cnv_locus_type_id))
		dc= {}
		query = self.metadata.bind.execute("select s1.id as sid, s2.id as tid from (select id,chromosome,position from %s where locus_type_id!=%s) as s1, \
			(select id, chromosome, position from %s where locus_type_id=%s) as s2 \
			where s1.chromosome=s2.chromosome and s1.position=s2.position"%(Snps.table.name, cnv_locus_type_id, Snps.table.name, cnv_locus_type_id))
		for row in query:
			dc[row.sid] = row.tid
		
		sys.stderr.write("%s entries.\n"%(len(dc)))
		return dc
	
	def getSNPsID2SpecificTypeSNPsID(self, locus_type_id=1, priorTAIRVersion=False):
		"""
		2012.4.24
			Some SNP db ID were mixed up with recombination locus db ID / CNV locus db ID or vice versa.
			
			An example is some old SNP datasets in which tair8_chromosome, tair8_position are shared by one SNP and one recombination locus.
			In the eventual SNP dataset (locus was marked with DB ID), instead of SNP db ID, some SNPs have recombination locus DB ID.
			
			This function provides a map to translate these wrong DB IDs back to the correct ones
				based on tair8_chromosome, tair8_position from table snps.
			This function is similar to getNonCNVSnpsID2CNVSnpsID() and is more flexible.
		"""
		if priorTAIRVersion==True:
			chromosome_table_fieldname = 'tair8_chromosome'
			position_table_fieldname = 'tair8_position'
			endPositionFieldName = 'end_position'
		else:
			chromosome_table_fieldname = 'chromosome'
			position_table_fieldname = 'position'
			endPositionFieldName = 'end_position'
		sys.stderr.write("Getting a map from Snps.id to Snps.id with locus_type_id=%s based on matching (%s, %s)  ..."%\
						(locus_type_id, chromosome_table_fieldname, position_table_fieldname))
		dc= {}
		query = self.metadata.bind.execute("select s1.id as sid, s2.id as tid from (select id, %s, %s from %s) as s1, \
			(select id, %s, %s from %s where locus_type_id=%s) as s2 \
			where s1.%s=s2.%s and s1.%s=s2.%s"%(chromosome_table_fieldname, position_table_fieldname, Snps.table.name, \
			chromosome_table_fieldname, position_table_fieldname, Snps.table.name, locus_type_id, \
			chromosome_table_fieldname, chromosome_table_fieldname, position_table_fieldname, position_table_fieldname))
		counter = 0
		real_counter = 0
		for row in query:
			counter += 1
			if row.sid!=row.tid:
				real_counter += 1
			
			dc[row.sid] = row.tid
		
		sys.stderr.write("%s out of %s entries were linked to a differing Snps.id.\n"%(real_counter, len(dc)))
		return dc
	
	@classmethod
	def getCNVQCInGWA(cls, accession_id=None, cnv_type_id=None, min_size=None, min_no_of_probes=None,\
					chr=None, start=None, stop=None, cnv_method_id=None):
		"""
		2010-7-22
			add argument cnv_method_id
		2010-6-18
			replace (TableClass.size_affected>=min_size) with (TableClass.stop-TableClass.start+1)>=min_size
		2010-4-28
			add argument chr, start, stop
		2010-3-14
			become a classmethod
		2009-10-30
			get CNV QC calls from db in GWA format
		"""
		sys.stderr.write("Getting CNVQC calls for accession %s, method_id %s, cnv_type %s, min_size %s, min_no_of_probes %s ..."%\
						(accession_id, cnv_method_id, cnv_type_id, min_size, min_no_of_probes))
		TableClass = CNVQCCall
		query = TableClass.query.filter_by(accession_id=accession_id)
		if min_size is not None:
			query = query.filter((TableClass.stop-TableClass.start+1)>=min_size)
		if cnv_type_id is not None:
			query = query.filter_by(cnv_type_id=cnv_type_id)
		if min_no_of_probes is not None:
			query = query.filter(TableClass.no_of_probes_covered>=min_no_of_probes)
		if cnv_method_id is not None:
			query = query.filter_by(cnv_method_id=cnv_method_id)
		
		# 2010-4-28
		query = cls.limitQueryByChrPosition(query, TableClass, chr, start, stop)
			
		from pymodule import GenomeWideResult, DataObject
		cnv_qc_accession = CNVQCAccession.get(accession_id)
		
		gwr_name = "CNVQCCall %s (%s) type %s, method %s"%(cnv_qc_accession.original_id, cnv_qc_accession.id, \
										cnv_type_id, cnv_method_id)
		if cnv_qc_accession.ecotype_id:
			gwr_name += " e-id: %s"%(cnv_qc_accession.ecotype_id)
		gwr = GenomeWideResult(name=gwr_name)
		gwr.data_obj_ls = []	#list and dictionary are crazy references.
		gwr.data_obj_id2index = {}
		# 2010-3-18 custom
		gwr.array_id = None
		gwr.cnv_qc_accession_id = cnv_qc_accession.id
		gwr.ecotype_id = cnv_qc_accession.ecotype_id
		gwr.nativename = cnv_qc_accession.original_id
		
		genome_wide_result_id = id(gwr)
		
		for row in query:
			data_obj = DataObject(chromosome=row.chromosome, position=row.start, stop_position=row.stop, value=row.cnv_type_id*2-3)
			data_obj.comment = 'cnv_type: %s, cnv_method %s, size_affected %s, no_of_probes_covered %s, copy_number %s, score %s.\n'%\
								(row.cnv_type_id, row.cnv_method_id, row.size_affected, row.no_of_probes_covered, row.copy_number,\
								row.score)
			if row.comment:
				data_obj.comment += '\n' + row.comment
			data_obj.genome_wide_result_name = gwr_name
			data_obj.genome_wide_result_id = genome_wide_result_id
			gwr.add_one_data_obj(data_obj)
		if gwr.max_value<3:	# insertion at y=3
			gwr.max_value=3
		if gwr.min_value>-1:	# deletion at y = -1
			gwr.min_value = -1
		sys.stderr.write(" %s segments. Done.\n"%(len(gwr.data_obj_ls)))
		return gwr
	
	def getCNVInGWA(self, cnv_method_id=None,cnv_type_id=None, chr=None, start=None, \
						stop=None, min_size=None, min_no_of_probes=None):
		"""
		2010-8-5
		"""
		sys.stderr.write("Getting CNV calls for cnv_type %s, min_size %s, min_no_of_probes %s method_id %s ..."%\
						(cnv_type_id, min_size, min_no_of_probes, cnv_method_id))
		TableClass = CNV
		query = TableClass.query
		if min_size is not None:
			query = query.filter((TableClass.stop-TableClass.start+1)>=min_size)
		if cnv_type_id is not None:
			query = query.filter_by(cnv_type_id=cnv_type_id)
		if min_no_of_probes is not None:
			query = query.filter(TableClass.no_of_probes_covered>=min_no_of_probes)
		if cnv_method_id is not None:
			query = query.filter_by(cnv_method_id=cnv_method_id)
		
		# 2010-4-28
		query = cls.limitQueryByChrPosition(query, TableClass, chr, start, stop)
			
		from pymodule import GenomeWideResult, DataObject
		
		gwr_name = "CNV type %s, method %s"%(cnv_type_id, cnv_method_id)
		
		gwr = GenomeWideResult(name=gwr_name)
		gwr.data_obj_ls = []	#list and dictionary are crazy references.
		gwr.data_obj_id2index = {}
		# 2010-3-18 custom
		gwr.array_id = None
		
		genome_wide_result_id = id(gwr)
		
		for row in query:
			data_obj = DataObject(chromosome=row.chromosome, position=row.start, stop_position=row.stop, value=row.cnv_type_id*2-3)
			data_obj.comment = 'cnv_type: %s, cnv_method %s, size_affected %s, no_of_probes_covered %s, score %s.\n'%\
								(row.cnv_type_id, row.cnv_method_id, row.size_affected, row.no_of_probes_covered,\
								row.score)
			if row.comment:
				data_obj.comment += '\n' + row.comment
			data_obj.genome_wide_result_name = gwr_name
			data_obj.genome_wide_result_id = genome_wide_result_id
			gwr.add_one_data_obj(data_obj)
		if gwr.max_value<3:	# insertion at y=3
			gwr.max_value=3
		if gwr.min_value>-1:	# deletion at y = -1
			gwr.min_value = -1
		sys.stderr.write(" %s segments. Done.\n"%(len(gwr.data_obj_ls)))
		return gwr
	
	def getCNVArrayCallInGWA(self, array_id=None, cnv_method_id=None,cnv_type_id=None, \
							chr=None, start=None, stop=None, \
							min_size=None, min_no_of_probes=None):
		"""
		2010-8-5
		"""
		sys.stderr.write("Getting CNVArrayCall for array %s, method_id %s, cnv_type %s, min_size %s, min_no_of_probes %s ..."%\
						(array_id, cnv_method_id, cnv_type_id, min_size, min_no_of_probes))
		TableClass = CNVArrayCall
		query = TableClass.query.filter_by(array_id=array_id)
		if min_size is not None:
			query = query.filter(TableClass.cnv.has((CNV.stop-CNV.start+1)>=min_size))
		if cnv_type_id is not None:
			query = query.filter(TableClass.cnv.has(cnv_type_id=cnv_type_id))
		if min_no_of_probes is not None:
			query = query.filter(TableClass.cnv.has(CNV.no_of_probes_covered>=min_no_of_probes))
		if cnv_method_id is not None:
			query = query.filter_by(cnv_method_id=cnv_method_id)
		if chr:
			query = query.filter(TableClass.cnv.has(chromosome=chr))
		
		if start and stop:
			query = query.filter(or_(and_(TableClass.cnv.has(CNV.start>=start), TableClass.cnv.has(CNV.start<=stop)), \
									and_(TableClass.cnv.has(CNV.stop>=start), TableClass.cnv.has(CNV.stop<=stop))))
		
		from pymodule import GenomeWideResult, DataObject
		array = ArrayInfo.get(array_id)
		
		try:
			rows = CNVArrayCall.table.metadata.bind.execute("select * from view_array where array_id=%s"%array_id)
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()[0]))
			traceback.print_exc()
			sys.stderr.write("Error in running select * from view_array where array_id=%s via CNVArrayCall.table.metadata.bind.execute\n"%\
							array_id)
			rows = self.metadata.bind.execute("select * from view_array where array_id=%s"%array_id)
		ecotype_nativename = None
		for row in rows:
			ecotype_nativename = row.maternal_nativename
			break
		
		gwr_name = "CNVArrayCall type %s of %s (a-id %s, e-id %s) by method %s"%(cnv_type_id, ecotype_nativename, \
												array.id, array.maternal_ecotype_id, cnv_method_id)
		gwr = GenomeWideResult(name=gwr_name)
		gwr.data_obj_ls = []	#list and dictionary are crazy references.
		gwr.data_obj_id2index = {}
		# 2010-3-18 custom
		gwr.array_id = array.id
		gwr.ecotype_id = array.maternal_ecotype_id
		gwr.nativename = ecotype_nativename
		
		genome_wide_result_id = id(gwr)
		
		for row in query:
			data_obj = DataObject(chromosome=row.cnv.chromosome, position=row.cnv.start, stop_position=row.cnv.stop, value=row.cnv.cnv_type_id*2-3)
			data_obj.comment = 'cnv_type: %s, cnv_method %s, size_affected %s, no_of_probes_covered %s, score %s.\n'%\
								(row.cnv.cnv_type_id, row.cnv_method_id, row.cnv.size_affected, row.cnv.no_of_probes_covered, \
								row.score)
			if row.comment:
				data_obj.comment += '\n' + row.comment
			data_obj.genome_wide_result_name = gwr_name
			data_obj.genome_wide_result_id = genome_wide_result_id
			gwr.add_one_data_obj(data_obj)
		if gwr.max_value<3:	# insertion at y=3
			gwr.max_value=3
		if gwr.min_value>-1:	# deletion at y = -1
			gwr.min_value = -1
		sys.stderr.write(" %s segments. Done.\n"%(len(gwr.data_obj_ls)))
		return gwr
		
	
	@classmethod
	def getSequenceRefPosInGWA(cls, accession_id=None, min_size=None, min_no_of_probes=None,\
							chr=None, start=None, stop=None, version=1):
		"""
		2010-6-18
			change default of min_size to None
		2010-6-17
			add argument version
		2010-4-28
			add argument chr, start, stop
		2010-3-13
			get ref genome position data of sequence fragments from one accession from db in GWA format
		"""
		sys.stderr.write("Getting positions of sequences from accession %s on ref genome, min_size %s, min_no_of_probes %s ..."%\
						(accession_id, min_size, min_no_of_probes))
		TableClass = SequenceFragmentRefPos
		query = TableClass.query.filter(SequenceFragmentRefPos.sequence_fragment.has(accession_id=accession_id)).\
			filter_by(version=version)
		
		if min_size is not None:
			query = query.filter((TableClass.stop-TableClass.start+1)>=min_size)
		if min_no_of_probes is not None:
			query = query.filter(TableClass.no_of_probes_covered>=min_no_of_probes)
		
		# 2010-4-28
		query = cls.limitQueryByChrPosition(query, TableClass, chr, start, stop)
		
		from pymodule import GenomeWideResult, DataObject
		cnv_qc_accession = CNVQCAccession.get(accession_id)
		
		gwr_name = "RefCoverage %s (%s)"%(cnv_qc_accession.original_id, cnv_qc_accession.id)
		if cnv_qc_accession.ecotype_id:
			gwr_name += " e-id: %s"%(cnv_qc_accession.ecotype_id)
		gwr = GenomeWideResult(name=gwr_name)
		gwr.data_obj_ls = []	#list and dictionary are crazy references.
		gwr.data_obj_id2index = {}
		# 2010-3-18 custom
		gwr.array_id = None
		gwr.cnv_qc_accession_id = cnv_qc_accession.id
		gwr.ecotype_id = cnv_qc_accession.ecotype_id
		gwr.nativename = cnv_qc_accession.original_id
		
		genome_wide_result_id = id(gwr)
		
		for row in query:
			data_obj = DataObject(chromosome=row.chromosome, position=row.start, stop_position=row.stop, value=1,\
								target_seq_id=row.sequence_fragment_id, target_start=row.fragment_start, \
								target_stop=row.fragment_stop)
			data_obj.comment = 'genome span %s, size_difference %s, fragment id %s (%s, size=%s), fragment start %s, fragment stop %s, fragment span %s.\n'%\
										(row.stop-row.start, row.size_difference, row.sequence_fragment_id, \
										row.sequence_fragment.short_name, row.sequence_fragment.size, \
										row.fragment_start, row.fragment_stop, abs(row.fragment_stop-row.fragment_start))
			if row.comment:
				data_obj.comment += ' ' + row.comment
			data_obj.genome_wide_result_name = gwr_name
			data_obj.genome_wide_result_id = genome_wide_result_id
			gwr.add_one_data_obj(data_obj)
		
		# 2010-3-11 set a proper range on the y-axis.
		gwr.max_value= 2
		gwr.min_value = -1
		sys.stderr.write(" %s segments. Done.\n"%(len(gwr.data_obj_ls)))
		return gwr
	
	@classmethod
	def getOneSequenceRefPosInGWA(cls, sequence_fragment_id=None, min_size=None, min_no_of_probes=None, \
								chr=None, start=None, stop=None, version=1):
		"""
		2010-11-15
			add argument chr, start, stop to restrict the segments to that space
			report the version argument as well
		2010-6-18
			change default of min_size to None
		2010-6-14
			add version
		2010-4-15
			get ref genome position data of only one sequence fragment from db in GWA format
				fragment positions are treated as y-axis values.
				called by GenomeBrowser.py
		"""
		sys.stderr.write("Getting positions of sequences from sequence fragment %s on ref genome, min_size %s, \
						min_no_of_probes %s version %s ..."%\
						(sequence_fragment_id, min_size, min_no_of_probes, version))
		TableClass = SequenceFragmentRefPos
		query = TableClass.query.filter(TableClass.sequence_fragment.has(id=sequence_fragment_id)).filter_by(version=version)
		query = cls.limitQueryByChrPosition(query, TableClass, chr, start, stop)
		if min_size is not None:
			query = query.filter((TableClass.stop-TableClass.start+1)>=min_size)
		if min_no_of_probes is not None:
			query = query.filter(TableClass.no_of_probes_covered>=min_no_of_probes)
		
		from pymodule import GenomeWideResult, DataObject
		sequence_fragment = SequenceFragment.get(sequence_fragment_id)
		cnv_qc_accession = sequence_fragment.accession
		gwr_name = "Fragment %s (%s) from %s"%(sequence_fragment.short_name, sequence_fragment.id, cnv_qc_accession.original_id)
		if cnv_qc_accession.ecotype_id:
			gwr_name += " e-id: %s"%(cnv_qc_accession.ecotype_id)
		gwr = GenomeWideResult(name=gwr_name)
		gwr.data_obj_ls = []	#list and dictionary are crazy references.
		gwr.data_obj_id2index = {}
		# 2010-3-18 custom
		gwr.array_id = None
		gwr.cnv_qc_accession_id = cnv_qc_accession.id
		gwr.ecotype_id = cnv_qc_accession.ecotype_id
		gwr.nativename = cnv_qc_accession.original_id
		gwr.sequence_fragment_id = sequence_fragment_id
		
		genome_wide_result_id = id(gwr)
		
		for row in query:
			data_obj = DataObject(chromosome=row.chromosome, position=row.start, stop_position=row.stop, value=1,\
								target_seq_id=row.sequence_fragment_id, target_start=row.fragment_start, \
								target_stop=row.fragment_stop)
			data_obj.comment = 'genome span %s, size_difference %s, fragment id %s (%s, size=%s), fragment start %s, fragment stop %s, fragment span %s.'%\
										(row.stop-row.start, row.size_difference, row.sequence_fragment_id, \
										row.sequence_fragment.short_name, row.sequence_fragment.size, \
										row.fragment_start, row.fragment_stop, abs(row.fragment_stop-row.fragment_start))
			if row.comment:
				data_obj.comment += ' ' + row.comment
			data_obj.genome_wide_result_name = gwr_name
			data_obj.genome_wide_result_id = genome_wide_result_id
			gwr.add_one_data_obj(data_obj)
			
			min_value = min(row.fragment_start, row.fragment_stop)
			max_value = max(row.fragment_start, row.fragment_stop)
			if gwr.min_value is None or min_value<gwr.min_value:
				gwr.min_value = min_value
			if gwr.max_value is None or max_value>gwr.max_value:
				gwr.max_value = max_value
		sys.stderr.write(" %s segments. Done.\n"%(len(gwr.data_obj_ls)))
		return gwr
	
	def getSequenceRefPosInMultiGWA(self, accession_id=None, min_size=None, min_no_of_probes=None,\
							chr=None, start=None, stop=None, version=1):
		"""
		2010-11-15
			bugfix: pass argument version to getOneSequenceRefPosInGWA()
		2010-6-14
			add version
		2010-4-28
			for each sequence fragment from accession_id that meets the requirement
				get ref genome position data of only one sequence fragment from db in GWA format
				-- getOneSequenceRefPosInGWA()
			return a list of gwrs
				
			fragment positions are treated as y-axis values.
			called by GenomeBrowser.py
		"""
		
		sql_string = "select distinct t.sequence_fragment_id from %s t, %s s where s.id=t.sequence_fragment_id and t.version=%s"%\
					(SequenceFragmentRefPos.table.name, SequenceFragment.table.name, version)
		
		condition_ls = []
		if accession_id:
			condition_ls.append("s.accession_id=%s"%accession_id)
		if chr and start and stop:
			condition_ls.append("t.chromosome=%s"%chr)
			condition_ls.append("((t.start>=%s and t.start<=%s) or (t.stop>=%s and t.stop<=%s))"%\
							(start, stop, start, stop))
		if min_size is not None:
			condition_ls.append("(t.stop-t.start+1)>=%s"%min_size)
		if min_no_of_probes is not None:
			condition_ls.append("t.no_of_probes_covered>=%s"%min_no_of_probes)
		
		condition = ' and '.join(condition_ls)
		sql_string += ' and %s'%condition
		print sql_string
		rows = self.metadata.bind.execute(sql_string)
		count = 0
		
		sequence_fragment_id_set = set()
		genome_wide_result_ls = []
		for row in rows:
			gwr = self.getOneSequenceRefPosInGWA(row.sequence_fragment_id, min_size=min_size, \
										min_no_of_probes=min_no_of_probes, chr=chr, start=start, stop=stop, \
										version=version)
			genome_wide_result_ls.append(gwr)
		return genome_wide_result_ls
	
	@classmethod
	def getOneSequence2ProbeInGWA(cls, sequence_fragment_id=None, min_size=None, min_no_of_probes=None, \
								min_no_of_identities=25):
		"""
		2010-4-18
			get all the matching probes of only one sequence fragment from db in GWA format
				min_size and min_no_of_probes are not used.
				fragment positions are treated as y-axis values.
				called by GenomeBrowser.py
		"""
		sys.stderr.write("Getting matching probes for sequence fragment %s on ref genome, min_no_of_identities%s ..."%\
						(sequence_fragment_id, min_no_of_identities))
		TableClass = SequenceFragment2Probe
		query = TableClass.query.filter(TableClass.sequence_fragment.has(id=sequence_fragment_id))
		
		if min_no_of_probes is not None:
			query = query.filter(TableClass.no_of_identities>=min_no_of_identities)
		
		from pymodule import GenomeWideResult, DataObject
		sequence_fragment = SequenceFragment.get(sequence_fragment_id)
		cnv_qc_accession = sequence_fragment.accession
		gwr_name = "Fragment %s (%s) from %s"%(sequence_fragment.short_name, sequence_fragment.id, cnv_qc_accession.original_id)
		if cnv_qc_accession.ecotype_id:
			gwr_name += " e-id: %s"%(cnv_qc_accession.ecotype_id)
		gwr = GenomeWideResult(name=gwr_name)
		gwr.data_obj_ls = []	#list and dictionary are crazy references.
		gwr.data_obj_id2index = {}
		# 2010-3-18 custom
		gwr.array_id = None
		gwr.cnv_qc_accession_id = cnv_qc_accession.id
		gwr.ecotype_id = cnv_qc_accession.ecotype_id
		gwr.nativename = cnv_qc_accession.original_id
		gwr.sequence_fragment_id = sequence_fragment_id
		
		genome_wide_result_id = id(gwr)
		
		for row in query:
			data_obj = DataObject(chromosome=row.chromosome, position=row.start, stop_position=row.stop, value=1,\
								target_seq_id=row.sequence_fragment_id, target_start=row.fragment_start, \
								target_stop=row.fragment_stop)
			data_obj.comment = 'probe id %s, no_of_identities %s, fragment id %s (%s, size=%s), fragment start %s, fragment stop %s, fragment span %s.'%\
										(row.probe_id, row.no_of_identities, row.sequence_fragment_id, \
										row.sequence_fragment.short_name, row.sequence_fragment.size, \
										row.fragment_start, row.fragment_stop, abs(row.fragment_stop-row.fragment_start))
			if row.comment:
				data_obj.comment += ' ' + row.comment
			data_obj.genome_wide_result_name = gwr_name
			data_obj.genome_wide_result_id = genome_wide_result_id
			gwr.add_one_data_obj(data_obj)
			
			min_value = min(row.fragment_start, row.fragment_stop)
			max_value = max(row.fragment_start, row.fragment_stop)
			if gwr.min_value is None or min_value<gwr.min_value:
				gwr.min_value = min_value
			if gwr.max_value is None or max_value>gwr.max_value:
				gwr.max_value = max_value
		sys.stderr.write(" %s segments. Done.\n"%(len(gwr.data_obj_ls)))
		return gwr
	
	def getSequence2ProbeInMultiGWA(self, accession_id=None, min_size=None, min_no_of_probes=None, min_no_of_identities=25,\
							chr=None, start=None, stop=None):
		"""
		2010-4-28
			for each sequence fragment from accession_id that meets the requirement
				get all the matching probes of one sequence fragment from db in GWA format
				-- getOneSequence2ProbeInGWA()
			return a list of gwrs
			
			min_size and min_no_of_probes are not used.
			fragment positions are treated as y-axis values.
			called by GenomeBrowser.py
		"""
		
		sql_string = "select distinct t.sequence_fragment_id from %s t, %s s where s.id=t.sequence_fragment_id "%\
					(SequenceFragment2Probe.table.name, SequenceFragment.table.name)
		
		condition_ls = []
		if accession_id:
			condition_ls.append("s.accession_id=%s"%accession_id)
		if chr and start and stop:
			condition_ls.append("t.chromosome=%s"%chr)
			condition_ls.append("((t.start>=%s and t.start<=%s) or (t.stop>=%s and t.stop<=%s))"%\
							(start, stop, start, stop))
		condition = ' and '.join(condition_ls)
		sql_string += ' and %s'%condition
		print sql_string
		rows = self.metadata.bind.execute(sql_string)
		count = 0
		
		sequence_fragment_id_set = set()
		genome_wide_result_ls = []
		for row in rows:
			gwr = self.getOneSequence2ProbeInGWA(row.sequence_fragment_id, min_size=min_size, \
												min_no_of_probes=min_no_of_probes, \
												min_no_of_identities=min_no_of_identities)
			genome_wide_result_ls.append(gwr)
		return genome_wide_result_ls
	
	@classmethod
	def getResultsAndTestResultsClass(cls, top_snp_test_type_id=None, results_type=None):
		"""
		2010-2-24
			this function helps other functions to figure out which ResultsClass and TestResultClass 
				if top_snp_test_type_id is not None, get the corresponding CandidateGeneTopSNPTestRMType, then get the results_type
				use results_type to figure out which ResultsClass and TestResultClass to use
		"""
		if top_snp_test_type_id is not None:
			test_type = CandidateGeneTopSNPTestRMType.get(top_snp_test_type_id)
			results_type = test_type.results_type
		if results_type==3:
			ResultsClass = ResultsMethod
			TestResultClass = CandidateGeneTopSNPTestRG
		elif results_type==2:
			ResultsClass = ResultsByGene
			TestResultClass = CandidateGeneTopSNPTest
		elif results_type==1:
			ResultsClass = ResultsMethod
			TestResultClass = CandidateGeneTopSNPTestRM
		else:	# default
			ResultsClass = None
			TestResultClass = None
		return (ResultsClass, TestResultClass)
	
	def getCNVCallInGWA(self, array_id=None, cnv_type_id=None, cnv_method_id=None, min_size=None, min_no_of_probes=None,\
					chr=None, start=None, stop=None):
		"""
		2010-7-28
			fix a bug. argument cnv_method_id wasn't used before.
		2010-6-18
			replace (TableClass.size_affected>=min_size) with (TableClass.stop-TableClass.start+1)>=min_size
		2010-4-28
			add argument chr, start, stop
		2010-3-17
			get CNV calls from db in GWA format
		"""
		sys.stderr.write("Getting CNV calls for array_id %s, cnv_type %s, cnv method %s, min_size %s, min_no_of_probes %s chr %s start %s stop %s ..."%\
						(array_id, cnv_type_id, cnv_method_id, min_size, min_no_of_probes, chr, start, stop ))
		TableClass = CNVCall
		query = TableClass.query.filter_by(array_id=array_id)
		if min_size is not None:
			query = query.filter((TableClass.stop-TableClass.start+1)>=min_size)
		if cnv_type_id is not None:
			query = query.filter_by(cnv_type_id=cnv_type_id)
		if min_no_of_probes is not None:
			query = query.filter(TableClass.no_of_probes_covered>=min_no_of_probes)
		if cnv_method_id is not None:
			query = query.filter_by(cnv_method_id=cnv_method_id)
		
		# 2010-4-28
		query = self.limitQueryByChrPosition(query, TableClass, chr, start, stop)
			
		
		from pymodule import GenomeWideResult, DataObject
		array = ArrayInfo.get(array_id)
		
		try:
			rows = CNVCall.table.metadata.bind.execute("select * from view_array where array_id=%s"%array_id)
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()[0]))
			traceback.print_exc()
			sys.stderr.write("Error in running select * from view_array where array_id=%s via CNVCall.table.metadata.bind.execute\n"%\
							array_id)
			rows = self.metadata.bind.execute("select * from view_array where array_id=%s"%array_id)
		ecotype_nativename = None
		for row in rows:
			ecotype_nativename = row.maternal_nativename
			break
		
		gwr_name = "CNV type %s of %s (a-id %s, e-id %s) by method %s"%(cnv_type_id, ecotype_nativename, \
												array.id, array.maternal_ecotype_id, cnv_method_id)
		gwr = GenomeWideResult(name=gwr_name)
		gwr.data_obj_ls = []	#list and dictionary are crazy references.
		gwr.data_obj_id2index = {}
		# 2010-3-18 custom
		gwr.array_id = array.id
		gwr.ecotype_id = array.maternal_ecotype_id
		gwr.nativename = ecotype_nativename
		
		genome_wide_result_id = id(gwr)
		
		for row in query:
			data_obj = DataObject(chromosome=row.chromosome, position=row.start, stop_position=row.stop, \
								value=row.amplitude)
			data_obj.comment = 'cnv_type: %s, cnv_method: %s, size_affected %s, no_of_probes_covered %s, amplitude %s.\n'%\
								(row.cnv_type.short_name, row.cnv_method.short_name, row.size_affected, row.no_of_probes_covered, row.amplitude)
			data_obj.comment += 'start_probe_id %s, stop_probe_id %s.\n'%\
								(row.start_probe_id, row.stop_probe_id)
			if row.comment:
				data_obj.comment += '\n'+row.comment
			data_obj.genome_wide_result_name = gwr_name
			data_obj.genome_wide_result_id = genome_wide_result_id
			gwr.add_one_data_obj(data_obj)
		if gwr.max_value<3:		# insertion at y=3
			gwr.max_value=3
		if gwr.min_value>-1:	# deletion at y = -1
			gwr.min_value = -1
		sys.stderr.write(" %s segments. Done.\n"%(len(gwr.data_obj_ls)))
		return gwr
	
	def getSequenceFragmentRefPosFromDBInRBDict(self, data_source_id=None, ecotype_id=None, \
								min_QC_segment_size=None, min_no_of_probes=None, min_reciprocal_overlap=0.6, version=1):
		"""
		2010-6-14
			add argument version
		2010-3-17 copied from CNV.getLerContigSpanDataFromDB() of misc.py
		2010-1-28
			This data provides information regarding which parts of Col genome are actually covered by the Ler contigs,
			which is key to get an accurate FPR for the tiling-array deletions because you'll know for sure this deletion
			is a false positive only when the deletion falls into the Ler contig coverage.
			
		"""
		sys.stderr.write("Getting SequenceFragmentRefPos data ... \n")		
		sql_string = "select a.ecotype_id, p.chromosome, p.start, p.stop, p.size_difference, p.no_of_probes_covered,\
					p.sequence_fragment_id, p.id \
					from %s p, %s f, %s a where f.accession_id=a.id and p.sequence_fragment_id=f.id and p.version=%s"%\
					(SequenceFragmentRefPos.table.name, SequenceFragment.table.name, \
					CNVQCAccession.table.name, version)
		if data_source_id is not None:
			sql_string += " and a.data_source_id=%s"%data_source_id
		if ecotype_id is not None:
			sql_string += " and a.ecotype_id=%s"%ecotype_id
		if min_no_of_probes is not None:
			sql_string += " and p.no_of_probes_covered>=%s"%min_no_of_probes
		if min_QC_segment_size is not None:
			sql_string += " and p.stop-p.start>=%s"%min_QC_segment_size
		
		sql_string += " order by RAND()"	# 2010-1-26 random ordering to optimize the binary_tree, not needed for RBDict.
		
		rows = self.metadata.bind.execute(sql_string)
		count = 0
		ecotype_id2span_data = {}
		from pymodule.RBTree import RBDict	# 2010-1-26 RBDict is more efficient than binary_tree.
		from pymodule.CNV import CNVSegmentBinarySearchTreeKey, leftWithinRightAlsoEqualCmp
		for row in rows:
			if row.ecotype_id not in ecotype_id2span_data:
				ecotype_id2span_data[row.ecotype_id] = RBDict(cmpfn=leftWithinRightAlsoEqualCmp)	# 2010-1-28 left within right is also equal.
			
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=row.chromosome, span_ls=[row.start, row.stop], \
													min_reciprocal_overlap=min_reciprocal_overlap)
			ecotype_id2span_data[row.ecotype_id][segmentKey] = (row.chromosome, row.start, row.stop, row.size_difference, \
															row.no_of_probes_covered, row.sequence_fragment_id, row.id)
			
			count += 1
		
		import math
		for ecotype_id, tree in ecotype_id2span_data.iteritems():
			print "\tDepth of Ecotype %s's tree: %d holding %s items." % (ecotype_id, tree.depth(), len(tree))
			print "\tOptimum Depth: %f (%d) (%f%% depth efficiency)" % (tree.optimumdepth(), math.ceil(tree.optimumdepth()),
															  math.ceil(tree.optimumdepth()) / tree.depth())
		
		sys.stderr.write("\t%s span data points for %s ecotypes. Done.\n"%(count, len(ecotype_id2span_data)))
		return ecotype_id2span_data
	
	def getArrayIDLsGivenEcotypeID(self, ecotype_id):
		"""
		2010-3-18
			type of ecotype_id could be str or integer
		"""
		rows = self.metadata.bind.execute("select * from stock_250k.view_array where maternal_ecotype_id=%s and maternal_ecotype_id=paternal_ecotype_id"%ecotype_id)
		array_id_ls = []
		for row in rows:
			array_id_ls.append(row.array_id)
		return array_id_ls
	
	def getAccessionIDInfoGivenSQLQuery(self, sql_string):
		"""
		2010-4-27
			for auto completion in GenomeBrowser.py
		"""
		rows = self.metadata.bind.execute(sql_string)
		count = 0
		
		list_id_ls = []
		list_id2index = {}
		list_label_ls = []
		for row in rows:
			list_id_ls.append(row.id)
			data_source = DataSource.get(row.data_source_id)
			label = '%s_%s_%s_%s'%(row.id, row.original_id, row.ecotype_id, data_source.short_name)
			list_label_ls.append(label)
			list_id2index[row.id] = len(list_id2index)
		from pymodule import PassingData
		return PassingData(list_id_ls=list_id_ls, list_id2index=list_id2index, list_label_ls=list_label_ls)
	
	def getAccessionIDInfoWithDataInSequenceFragmentRefPos(self,):
		"""
		2010-4-27
			for auto completion in GenomeBrowser.py
		"""
		sql_string = "select distinct a.* \
				from %s p, %s f, %s a where f.accession_id=a.id and p.sequence_fragment_id=f.id"%\
				(SequenceFragmentRefPos.table.name, SequenceFragment.table.name, \
				CNVQCAccession.table.name)
		return self.getAccessionIDInfoGivenSQLQuery(sql_string)
	
	def getAccessionIDInfoWithDataInSequenceFragment2Probe(self,):
		"""
		2010-4-27
			for auto completion in GenomeBrowser.py
		"""
		sql_string = "select distinct a.* \
				from %s p, %s f, %s a where f.accession_id=a.id and p.sequence_fragment_id=f.id"%\
				(SequenceFragment2Probe.table.name, SequenceFragment.table.name, \
				CNVQCAccession.table.name)
		return self.getAccessionIDInfoGivenSQLQuery(sql_string)
	
	def getAccessionIDInfoWithDataInCNVQCCall(self,):
		"""
		2010-4-27
			for auto completion in GenomeBrowser.py
		"""
		sql_string = "select distinct a.* \
				from %s q, %s a where q.accession_id=a.id "%\
				(CNVQCCall.table.name, CNVQCAccession.table.name)
		return self.getAccessionIDInfoGivenSQLQuery(sql_string)
	
	def getArrayIDInfoWithDataInGivenTable(self, TableClass=CNVCall):
		"""
		2010-8-5
			change name from getArrayIDInfoWithDataInCNVCall to getArrayIDInfoWithDataInGivenTable
			add argument TableClass
		2010-4-27
			for auto completion in GenomeBrowser.py
		"""
		sql_string = "select distinct v.* from %s c, view_array v where c.array_id=v.array_id "%\
					(TableClass.table.name)
		rows = self.metadata.bind.execute(sql_string)
		count = 0
		
		list_id_ls = []
		list_id2index = {}
		list_label_ls = []
		for row in rows:
			list_id_ls.append(row.array_id)
			label = '%s_%s_%s'%(row.array_id, row.maternal_nativename, row.maternal_ecotype_id)
			list_label_ls.append(label)
			list_id2index[row.array_id] = len(list_id2index)
		from pymodule import PassingData
		return PassingData(list_id_ls=list_id_ls, list_id2index=list_id2index, list_label_ls=list_label_ls)
	
	@classmethod
	def limitQueryByChrPosition(cls, query, TableClass, chr=None, start=None, stop=None):
		"""
		2010-7-28
			chr and (start&stop) conditions get separated.
			One could just specify chr to get stuff on that particular chromosome.
		# 2010-4-28
			common function invoked by other member functions.
		"""
		if chr:
			query = query.filter(TableClass.chromosome==chr)
		
		if start and stop:
			query = query.filter(or_(and_(TableClass.start>=start, TableClass.start<=stop), \
									and_(TableClass.stop>=start, TableClass.stop<=stop)))
		return query
	
	def getSeqFragmentRefPosVersionInfo(self):
		"""
		2010-6-14
			for GenomeBrowser.py to autocomplete comboboxentry_seq_ref_pos_version
		"""
		sql_string = "select distinct version from %s srp"%\
					(SequenceFragmentRefPos.table.name)
		rows = self.metadata.bind.execute(sql_string)
		count = 0
		
		list_id_ls = []
		list_id2index = {}
		list_label_ls = []
		for row in rows:
			SFRP_object = SequenceFragmentRefPos.query.filter_by(version=row.version).first()
			list_id_ls.append(SFRP_object.version)
			label = '%s %s'%(row.version, SFRP_object.comment)
			list_label_ls.append(label)
			list_id2index[row.version] = len(list_id2index)
		from pymodule import PassingData
		return PassingData(list_id_ls=list_id_ls, list_id2index=list_id2index, list_label_ls=list_label_ls)
	
	def getIDShortNameInfoGivenSQLQuery(self, sql_string):
		"""
		2010-7-23
		"""
		rows = self.metadata.bind.execute(sql_string)
		count = 0
		
		list_id_ls = []
		list_id2index = {}
		list_label_ls = []
		for row in rows:
			list_id_ls.append(row.id)
			label = '%s %s'%(row.id, row.short_name)
			list_label_ls.append(label)
			list_id2index[row.id] = len(list_id2index)
		from pymodule import PassingData
		return PassingData(list_id_ls=list_id_ls, list_id2index=list_id2index, list_label_ls=list_label_ls)
	
	def getCNVMethodOrTypeInfoDataInCNVCallOrQC(self, TableClass=CNVMethod, CNVTableClass=CNVCall, extra_condition=None):
		"""
		2010-10-19
			add argument extra_condition
		2010-7-23
			set TableClass = CNVType to get type info.
			set CNVTableClass=CNVQCCall to get info for the QC table
		"""
		if TableClass==CNVType:
			field_name = 'cnv_type_id'
		elif TableClass==CNVMethod:
			field_name = 'cnv_method_id'
		else:
			field_name = 'cnv_method_id'
		where_sql_ls = [ "m.id=s.%s"%(field_name)]
		if extra_condition is not None:
			where_sql_ls.append(extra_condition)
		sql_string = "select distinct m.* from %s m, %s s where %s order by m.id"%\
					(TableClass.table.name, CNVTableClass.table.name, " and ".join(where_sql_ls))
		return self.getIDShortNameInfoGivenSQLQuery(sql_string)
	
	def getSNPChrPos2IDLs(self, priorTAIRVersion=False, locus_type_id=None):
		"""
		2012.3.9
			add argument locus_type_id
		2012.3.7
			key is changed from (chr, pos) to (chr, pos, end_position)
		2011-5-6
			returns both maps
		"""
		dict_to_return = {}
		snp_id2chr_pos = {}
		chr_pos_2snp_id = {}
		if priorTAIRVersion==True:
			chromosome_table_fieldname = 'tair8_chromosome'
			position_table_fieldname = 'tair8_position'
			endPositionFieldName = 'end_position'
		else:
			chromosome_table_fieldname = 'chromosome'
			position_table_fieldname = 'position'
			endPositionFieldName = 'end_position'
		
		where_sentence = 'where %s is not null and %s is not null'%(chromosome_table_fieldname, position_table_fieldname)
		if locus_type_id is not None:
			where_sentence += ' and locus_type_id=%s'%(locus_type_id)
		rows = self.metadata.bind.execute("select id, %s as chromosome, %s as position, %s from %s %s "%\
								(chromosome_table_fieldname, position_table_fieldname, endPositionFieldName, Snps.table.name, \
								where_sentence, ))
		for row in rows:
			if row.end_position is None:	#2012.3.7
				end_position = row.position
			else:
				end_position = row.end_position
			key = (row.chromosome, row.position, end_position)
			chr_pos_2snp_id[key] = row.id
			snp_id2chr_pos[row.id] = key
		dict_to_return['snp_id2chr_pos'] = snp_id2chr_pos
		dict_to_return['chr_pos_2snp_id'] = chr_pos_2snp_id
		return dict_to_return
	
	
	def getSNPChrPos2ID(self, keyType=1, priorTAIRVersion=False, locus_type_id=None):
		"""
		2012.3.9
			add argument locus_type_id
		2011-5-4
			bugfix: when priorTAIRVersion=True, alias the field names to chromosome and position
				(instead of tair8_chromosome, tair8_position)
		2011-2-15
			add argument priorTAIRVersion
				if =true:
					(chr,pos) would be from (tair8_chromosome, tair8_position)
				else:
					(chr,pos) would be from (chromosome, position)
		2011-2-15	
			report how many entries were fetched
		2011-1-24
			add status reporting
		2010-10-13
			return a dictionary which translates (chr, pos) to snps.id or vice versa.
			used in Calls2DB_250k.py
			
			keyType
				1: (chr,pos) is key. id is value
				2: id is key, (chr, pos) is value.
		"""
		sys.stderr.write("Getting a map between Snps.id and (chr,pos), keyType %s, priorTAIRVersion=%s, locus_type_id=%s ..."%\
						(keyType, priorTAIRVersion, locus_type_id))
		data = self.getSNPChrPos2IDLs(priorTAIRVersion, locus_type_id=locus_type_id)
		
		if keyType==1:
			dict_to_return = data.get('chr_pos_2snp_id')
		else:
			dict_to_return = data.get('snp_id2chr_pos')
		sys.stderr.write("%s entries.\n"%(len(dict_to_return)))
		return dict_to_return
	
	def getSNPID2ChrPos(self, priorTAIRVersion=False, locus_type_id=None):
		"""
		2012.3.9
			add argument locus_type_id
		2011-10-14 add argument priorTAIRVersion
		2010-10-13
			return a dictionary which translates snps.id to (chr, pos).
			used in DB_250k2data.py
			
		"""
		return self.getSNPChrPos2ID(keyType=2, priorTAIRVersion=priorTAIRVersion, locus_type_id=locus_type_id)
	
	
	def get_db_id_given_chr_pos2db_id(self, snp_id):
		"""
		2011-2-27
			snp_id could be chr_pos or just Snps.id
			
			called by Calls2DB_250k.py and others.
		"""
		snp_id = snp_id.split('_')[:2]
		chr_pos = tuple(map(int, snp_id))
		if len(chr_pos)==1 or (len(chr_pos)>=2 and chr_pos[1]==0):
			db_id = chr_pos[0]
		else:
			db_id = self.chr_pos2snp_id.get(chr_pos)
		return db_id
	
	def get_chr_pos_given_db_id2chr_pos(self, snp_id=None):
		"""
		2011-2-27
			snp_id could be chr_pos or just Snps.id
			
			called by QC_250k.py and Calls2DB_250k.py
		"""
		snp_id = snp_id.split('_')[:2]
		chr_pos = tuple(map(int, snp_id))
		if len(chr_pos)==1 or (len(chr_pos)>=2 and chr_pos[1]==0):
			db_id = chr_pos[0]
			chr_pos = self.snp_id2chr_pos.get(db_id)
		if chr_pos:
			snp_id = '%s_%s'%(chr_pos[0], chr_pos[1])
			return snp_id
		else:
			return None
	
	@property
	def snp_id2chr_pos(self,):
		"""
		2011-2-24
			convenient function
		"""
		if not self._snp_id2chr_pos:
			self.snp_id2chr_pos = False
		return self._snp_id2chr_pos
	
	@snp_id2chr_pos.setter
	def snp_id2chr_pos(self, argument_ls=[False, None]):
		"""
		2012.3.9
			replace argument priorTAIRVersion with argument_ls
			usage:
				db.snp_id2chr_pos = True	#if you need to set priorTAIRVersion to True and leave locus_type_id to default (None).
				db.snp_id2chr_pos = (True, None)	#if you need to set priorTAIRVersion to True and locus_type_id=None
				chr_pos = db.snp_id2chr_pos.get(1)
		2011-2-28
		"""
		argData = self.handleSNPChrPos2IDArguments(argument_ls)
		self._snp_id2chr_pos = self.getSNPChrPos2ID(keyType=2, priorTAIRVersion=argData.priorTAIRVersion,\
										locus_type_id=argData.locus_type_id)
	
	@property
	def chr_pos2snp_id(self,):
		"""
		2011-2-24
			convenient function
		"""
		if not self._chr_pos2snp_id:	#set it if it's empty
			self.chr_pos2snp_id = False
		return self._chr_pos2snp_id
	
	
	@chr_pos2snp_id.setter
	def chr_pos2snp_id(self, argument_ls=[False, None]):
		"""
		2012.3.9
			replace argument priorTAIRVersion with argument_ls
			usage:
				db.chr_pos2snp_id = True	#if you need to set priorTAIRVersion to True and leave locus_type_id to default (None).
				db.chr_pos2snp_id = (True,1)	# if you need to set priorTAIRVersion to True and locus_type_id=1
				snp_id = db.chr_pos2snp_id.get((1,657))
		2011-2-28
		"""
		argData = self.handleSNPChrPos2IDArguments(argument_ls)
		self._chr_pos2snp_id = self.getSNPChrPos2ID(keyType=1, priorTAIRVersion=argData.priorTAIRVersion,\
												locus_type_id=argData.locus_type_id)
	
	def handleSNPChrPos2IDArguments(self, argument_ls=[False, None]):
		"""
		2012.3.9
			called by chr_pos2snp_id or snp_id2chr_pos:
				db.snp_id2chr_pos = True	#if you need to set priorTAIRVersion to True and leave locus_type_id to default (None).
				db.snp_id2chr_pos = (True, None)	#if you need to set priorTAIRVersion to True and locus_type_id=None
				chr_pos = db.snp_id2chr_pos.get(1)
				
				db.chr_pos2snp_id = True	#if you need to set priorTAIRVersion to True and leave locus_type_id to default (None).
				db.chr_pos2snp_id = (True,1)	# if you need to set priorTAIRVersion to True and locus_type_id=1
				snp_id = db.chr_pos2snp_id.get((1,657))
		"""
		priorTAIRVersion = False
		locus_type_id = None
		if type(argument_ls)==list or type(argument_ls)==tuple:
			priorTAIRVersion = argument_ls[0]
			if len(argument_ls)>1:
				locus_type_id = argument_ls[1]
		else:	#2012.3.9 single argument (old way).
			priorTAIRVersion = argument_ls
		return PassingData(priorTAIRVersion=priorTAIRVersion, locus_type_id=locus_type_id)
	
	@property
	def cnv_id2chr_pos(self):
		"""
		2011-2-24
			get self._cnv_id2chr_pos
		"""
		if not self._cnv_id2chr_pos:
			sys.stderr.write("Warning: cnv_id2chr_pos has not been initialized. Run 'cnvMethodID=20;self.cnv_id2chr_pos=cnvMethodID;'.\n")
		return self._cnv_id2chr_pos
	
	@cnv_id2chr_pos.setter
	def cnv_id2chr_pos(self, cnv_method_id=None):
		"""
		2011-2-24
			setter of cnv_id2chr_pos
			
			usage:
				if db._cnv_method_id!=20:
					db.cnv_id2chr_pos = 20
				chr_pos = db.cnv_id2chr_pos.get(222244)
		"""
		sys.stderr.write("Getting a map between CNV.id and (chr,pos), cnv method %s ..."%(cnv_method_id))
		self._cnv_id2chr_pos = {}
		self._cnv_method_id = cnv_method_id
		rows = self.metadata.bind.execute("select id, chromosome, start, stop from %s where cnv_method_id=%s"\
										%(CNV.table.name, cnv_method_id))
		for row in rows:
			key = row.id
			value = (row.chromosome, row.start, row.stop)
			self._cnv_id2chr_pos[key] = value
		sys.stderr.write("%s entries. Done.\n"%(len(self._cnv_id2chr_pos)))
	
	#@classmethod
	def getResultMethodContent(self, results_method_id, results_directory=None, min_MAF=0.1, construct_chr_pos2index=False, \
						pdata=None, min_value_cutoff=None,):
		"""
		2012.3.23
			results_directory is equivalent to /Network/Data/250k/db/, not /Network/Data/250k/db/results/type_1/ (before).
			also use self.reScalePathByNewDataDir() to get the updated path.
		2011-3-10
			moved from GeneListRankTest
		2010-3-8
			add argument min_value_cutoff and pass it to getGenomeWideResultFromFile()
		2010-2-2
			add the db object "rm" to genome_wide_result to make other db info accessible
		2008-10-23
			if pdata doesn't have construct_chr_pos2index defined. otherwise, pdata overrides the option.
		2008-10-22
			before deciding do_log10_transformation based on analysis_method, try to get it from pdata
		2008-10-21
			add pdata to conceal the passing of chr_pos2index to getGenomeWideResultFromFile()
		2008-10-15
			cache the genome_wide_result under cls.genome_wide_result, if cls.genome_wide_result.results_id==rm.id, directly return that.
		2008-09-24
			add option construct_chr_pos2index
		2008-09-16
			if result_fname is not a file, return None
		2008-09-15
			use field smaller_score_more_significant from ResultsMethod to set do_log10_transformation
		2008-08-13
			split from getGeneID2MostSignificantHit()
		"""
		#genome_wide_result = getattr(cls, 'genome_wide_result', None)
		do_log10_transformation = getattr(pdata, 'do_log10_transformation', None)
		rm = ResultsMethod.get(results_method_id)
		from pymodule import getGenomeWideResultFromFile
		# 2011-3-21 no more caching
		#if genome_wide_result is not None and genome_wide_result.results_id==rm.id:
		#	return genome_wide_result
		
		if rm.analysis_method_id==13: #Huh -Bjarni
			sys.stderr.write("Analysis method id=%s is not supported.\n"%rm.analysis_method_id)
			return None
		"""
		if results_directory:	#given a directory where all results are.
			result_fname = os.path.join(results_directory, os.path.basename(rm.filename))
		else:
			result_fname = rm.filename
		"""
		result_fname = self.reScalePathByNewDataDir(filePath=rm.filename, newDataDir=results_directory)
		
		if do_log10_transformation is None:
			#based on the analysis method id, whether do -log() or not. it'll affect the later step of taking maximum pvalue out of SNPs associated with one gene
			if hasattr(rm, 'analysis_method'):
					if rm.analysis_method.smaller_score_more_significant==1:
						do_log10_transformation = True
					else:
						do_log10_transformation = False
			else:
				return None
		if pdata is None:
			pdata = PassingData()
		
		#2011-3-21 assign min_MAF to pdata only when pdata.min_MAF is None or doesn't exist.
		min_MAF_request = getattr(pdata, 'min_MAF', None)
		if min_MAF_request is None:
			pdata.min_MAF = min_MAF
		if not hasattr(pdata, 'construct_chr_pos2index'):	#if pdata doesn't have construct_chr_pos2index defined. otherwise, pdata overrides the option.
			pdata.construct_chr_pos2index = construct_chr_pos2index
		if os.path.isfile(result_fname):
			genome_wide_result = getGenomeWideResultFromFile(result_fname, do_log10_transformation=do_log10_transformation, \
													pdata=pdata, min_value_cutoff=min_value_cutoff)
			genome_wide_result.results_id = rm.id
			genome_wide_result.rm = rm	# 2010-2-2 add the db object "rm" to genome_wide_result to make other db info accessible
		else:
			sys.stderr.write("Skip. %s doesn't exist.\n"%result_fname)
			genome_wide_result = None
		self.genome_wide_result = genome_wide_result
		return genome_wide_result
	
	def getResultLs(self, call_method_id=None, analysis_method_id_ls=[], phenotype_method_id_ls=[], \
				call_method_id_ls=[], cnv_method_id=None):
		"""
		2012.6.28
			handle the query it differently when datasetConditionLS has only one entry
		2011-10-16
			add argument call_method_id_ls & cnv_method_id
		2011-10-12
			bugfix in filter by call_method_id
		2011-5-9
			given constraints, find all association results from ResultsMethod
		"""
		query = ResultsMethod.query
		datasetConditionLS = []
		if call_method_id:
			datasetConditionLS.append(ResultsMethod.call_method_id==call_method_id)
		if call_method_id_ls:
			datasetConditionLS.append(ResultsMethod.call_method_id.in_(call_method_id_ls))
		if cnv_method_id:
			datasetConditionLS.append(ResultsMethod.cnv_method_id==cnv_method_id)
		
		if datasetConditionLS:
			if len(datasetConditionLS)>1:
				query = query.filter(or_(*datasetConditionLS))
			elif len(datasetConditionLS)==1:	#only one
				query = query.filter(datasetConditionLS[0])
		
		if analysis_method_id_ls:
			query = query.filter(ResultsMethod.analysis_method_id.in_(analysis_method_id_ls))
		if phenotype_method_id_ls:
			query = query.filter(ResultsMethod.phenotype_method_id.in_(phenotype_method_id_ls))
		return query
	
	def getGeneList(self, list_type_id):
		"""
		2011-3-16
			moved from GeneListRankTest.py
		"""
		sys.stderr.write("Getting gene_list %s ... "%list_type_id)
		rows = GeneList.query.filter_by(list_type_id=list_type_id)
		candidate_gene_list = []
		for row in rows:
			candidate_gene_list.append(row.gene_id)
		sys.stderr.write("%s genes. Done.\n"%(len(candidate_gene_list)))
		return candidate_gene_list
	
	list_type_id2candidate_gene_list_info = {}
	def dealWithCandidateGeneList(self, list_type_id, return_set=False):
		"""
		2011-3-16
			moved from GeneListRankTest.py
		2008-10-30
			turn this into a classmethod
		2008-10-23
			add option return_set
		2008-10-15
			to make caching candidate gene list possible
		"""
		if list_type_id not in self.list_type_id2candidate_gene_list_info:	#internal cache
			candidate_gene_list = self.getGeneList(list_type_id)
			self.list_type_id2candidate_gene_list_info[list_type_id] = PassingData(candidate_gene_list=candidate_gene_list, \
																	candidate_gene_set=set(candidate_gene_list))
		if return_set:
			return self.list_type_id2candidate_gene_list_info[list_type_id].candidate_gene_set
		else:
			return self.list_type_id2candidate_gene_list_info[list_type_id].candidate_gene_list
	
	def getPhenotypeMethodLsGivenBiologyCategoryID(self, biology_category_id=None, access=None):
		"""
		2012.2.22
			add argument access. If None, no filter based on this.
				1=public
				2=restricted
		2011-10-17
		"""
		sys.stderr.write("Getting list of phenotype method with biology category id=%s, access=%s ..."%(biology_category_id, access))
		query = PhenotypeMethod.query
		if biology_category_id is not None:
			query = query.filter_by(biology_category_id=biology_category_id)
		if access is not None:
			query = query.filter_by(access=access)
		phenotype_method_ls = []
		for row in query:
			phenotype_method_ls.append(row)
		
		sys.stderr.write("%s phenotypes.\n"%(len(phenotype_method_ls)))
		return phenotype_method_ls
	
	def getResultPeakList(self, result_id_ls=[], result_peak_type_id=None):
		"""
		2012.3.10
			return a list of result_peak objects + unique biology categories
			
		"""
		sys.stderr.write("Filtering a list of %s results by checking ResultPeak (result_peak_type_id=%s)..."%(len(result_id_ls), \
																							result_peak_type_id))
		TableClass = ResultPeak
		query = TableClass.query.filter(TableClass.result_id.in_(result_id_ls))
		if result_peak_type_id is not None:
			query = query.filter_by(result_peak_type_id=result_peak_type_id)
		return query
	
	def constructRBDictFromResultPeak(self, result_id, result_peak_type_id, peakPadding=10000):
		"""
		2012.6.24
			add attributes to result_peakRBDict 
		2012.3.19
			moved from variation.src.TwoGWASPeakOverlap
		2011-10-16
		"""
		sys.stderr.write("Constructing RBDict for peaks from result %s, (peak type %s) ..."%(result_id, result_peak_type_id))
		from pymodule.CNV import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
		from pymodule.RBTree import RBDict
		result_peakRBDict = RBDict()
		result_peakRBDict.result_id = result_id	#2012.6.22
		result_peakRBDict.result_peak_type_id = result_peak_type_id	#2012.6.22
		result_peakRBDict.peakPadding = peakPadding
		query = ResultPeak.query.filter_by(result_id=result_id).filter_by(result_peak_type_id=result_peak_type_id)
		counter = 0
		real_counter = 0
		for row in query:
			counter += 1
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=row.chromosome, \
							span_ls=[max(1, row.start - peakPadding), row.stop + peakPadding], \
							min_reciprocal_overlap=1, result_peak_id=row.id)
							#2010-8-17 overlapping keys are regarded as separate instances as long as they are not identical.
			if segmentKey not in result_peakRBDict:
				result_peakRBDict[segmentKey] = []
			result_peakRBDict[segmentKey].append(row)
		sys.stderr.write("%s peaks. Done.\n"%counter)
		return result_peakRBDict
	
	def filterResultIDLsBasedOnResultPeak(self, result_id_ls, result_peak_type_id=None):
		"""
		2012.3.22
			moved from PairwiseGWASPeakOverlapPipeline.py
		2011-10-17
			return a list of ResultsMethod.id that have entries in ResultPeak with result_peak_type_id
		"""
		sys.stderr.write("Filtering a list of %s results by checking ResultPeak (result_peak_type_id=%s)..."%(len(result_id_ls), \
																							result_peak_type_id))
		new_result_id_ls = []
		for result_id in result_id_ls:
			firstPeak = ResultPeak.query.filter_by(result_id=result_id).filter_by(result_peak_type_id=result_peak_type_id).first()
			if firstPeak:
				new_result_id_ls.append(result_id)
		sys.stderr.write(" %s entries left.\n"%(len(new_result_id_ls)))
		return new_result_id_ls
	
	def getAssociationLocus(self, chromosome=None, start=None, stop=None, no_of_peaks=None, connectivity=None,\
						threshold=None, result_peak_ls=None):
		"""
		2012.6.28
		"""
		association_locus = AssociationLocus.query.filter_by(chromosome=chromosome).filter_by(start=start).filter_by(stop=stop).\
			filter_by(threshold=threshold).first()
		
		if not association_locus:
			association_locus = AssociationLocus(chromosome=chromosome, start=start, stop=stop, threshold=threshold,\
							no_of_peaks=no_of_peaks, connectivity=connectivity)
			association_locus.result_peak_ls = result_peak_ls
			self.session.add(association_locus)
			#self.session.flush()
		return association_locus
	
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = Stock_250kDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.setup()
	import pdb
	pdb.set_trace()
	
	#2010-8-6 complex query
	TableClass = CNVArrayCall
	array_id = 612
	query = TableClass.query.filter_by(array_id=array_id)
	min_size = 100
	query = query.filter(TableClass.cnv.has((CNV.stop-CNV.start+1)>=min_size))
	cnv_type_id = 1
	query = query.filter(TableClass.cnv.has(cnv_type_id=cnv_type_id))
	min_no_of_probes = 5
	query = query.filter(TableClass.cnv.has(CNV.no_of_probes_covered>=min_no_of_probes))
	cnv_method_id = 22
	query = query.filter_by(cnv_method_id=cnv_method_id)
	chr = 1
	query = query.filter(TableClass.cnv.has(chromosome=chr))
	start =25000
	stop = 2500000
	query = query.filter(or_(and_(TableClass.cnv.has(CNV.start>=start), TableClass.cnv.has(CNV.start<=stop)), \
								and_(TableClass.cnv.has(CNV.stop>=start), TableClass.cnv.has(CNV.stop<=stop))))
	counter = 0
	for row in query:
		print row.id, row.array_id, row.cnv_id
		counter += 1
		if counter>=10:
			break
	
	import sqlalchemy
	s = sqlalchemy.sql.select([AssociationOverlappingStat.table.c.call_method_id.distinct()])
	connection = instance.metadata.bind
	#result = connection.execute(s).fetchmany(3)
	result = connection.execute(s)
	for r in result:
		print r
		print dir(r)
		print r.call_method_id
		
	rows = GeneListType.query.all()
	for row in rows[:10]:
		print row.gene_list[0].list_type_id
