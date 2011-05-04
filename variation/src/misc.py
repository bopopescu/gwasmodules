#!/usr/bin/env python
"""
Examples:
	# enter debug mode
	misc.py -b
	
	# change the db connection setting 
	misc.py -z localhost -u yh
	
Description:
	It's a file with all sorts of classes that haven't gone on to be independent.
	
	All those classes are manually chosen to be run in Main.run().
"""
import os, sys, numpy

class DisplaySNPMatrix(object):
	"""
	2007-03-05
		show an image visualizing SNP data
	2007-06-05
		set aspect='auto' in imshow(), the default (pylab.image.rcParams['image.aspect'])='equal', which is bad
	"""
	def display_snp_matrix(input_fname, output_fname=None, need_sort=0, need_savefig=0, xlabel='', ylabel=''):
		import csv, Numeric, pylab
		reader = csv.reader(open(input_fname), delimiter='\t')
		header = reader.next()
		data_matrix = []
		for row in reader:
			data_row = row[2:]
			data_row = map(int, data_row)
			data_matrix.append(data_row)
		del reader
		data_matrix.reverse()	#2007-03-06 reverse() due to the imshow()'s y axis starting from bottom
		if need_sort:
			data_matrix.sort()
		data_matrix = Numeric.array(data_matrix)
		
		pylab.clf()
		pylab.imshow(data_matrix, aspect='auto', interpolation='nearest')	#2007-06-05
		pylab.colorbar()
		if xlabel:
			pylab.xticks([data_matrix.shape[1]/2], [xlabel])
		if ylabel:
			pylab.yticks([data_matrix.shape[0]/2], [ylabel])
		if need_savefig:
			pylab.savefig('%s.eps'%output_fname, dpi=300)
			pylab.savefig('%s.svg'%output_fname, dpi=300)
			pylab.savefig('%s.png'%output_fname, dpi=300)
		pylab.show()
	
	def make_snp_matrix_legend(value_ls, label_ls, output_fname=None):
		"""
		2007-10-25
			to pair with display_snp_matrix()
		"""
		import numpy, pylab
		label_ls_copy = label_ls[:]
		data_matrix = numpy.zeros([len(value_ls), 1], numpy.int)
		for i in range(len(value_ls)):
			data_matrix[i,0] = value_ls[i]
		label_ls_copy.reverse()	#pylab put the label on starting from the bottom of the picture
		pylab.clf()
		pylab.imshow(data_matrix, aspect='auto', interpolation='nearest')
		pylab.yticks(range(len(value_ls)), label_ls_copy, fontsize=60, verticalalignment='bottom', horizontalalignment='right')
		pylab.xticks([],[])
		if output_fname:
			pylab.savefig('%s.eps'%output_fname, dpi=300)
			pylab.savefig('%s.svg'%output_fname, dpi=300)
			pylab.savefig('%s.png'%output_fname, dpi=300)
		pylab.show()
	
	"""
	2007-10-25
	value_ls = range(3)
	label_ls = ['NA', 'A', 'C']
	make_snp_matrix_legend(value_ls, label_ls)
	
	"""


class DB149(object):
	
	"""
	2007-04-09
	"""
	def draw_SNP_gap_histogram(curs, snp_locus_table, output_fname, need_savefig=0):
		SNP_gap_ls = []
		prev_chromosome = None
		prev_position = None
		curs.execute("select chromosome, position from %s order by chromosome, position"%snp_locus_table)
		rows = curs.fetchall()
		for row in rows:
			chromosome, position = row
			if prev_chromosome==None or prev_position==None:
				prev_chromosome = chromosome
				prev_position = position
			elif chromosome==prev_chromosome:
				SNP_gap_ls.append(position-prev_position)
			prev_chromosome = chromosome
			prev_position = position
		import pylab
		pylab.clf()
		pylab.hist(SNP_gap_ls, 20)
		pylab.title("hist of SNP gap")
		if need_savefig:
			pylab.savefig('%s.eps'%output_fname, dpi=300)
			pylab.savefig('%s.svg'%output_fname, dpi=300)
			pylab.savefig('%s.png'%output_fname, dpi=300)
		pylab.show()
		return SNP_gap_ls
	
	@classmethod
	def blast_snp_segment(cls, curs, snp_locus_table, output_fname, database_fname, flanking_seq_length=12, \
						max_no_of_hits_to_be_outputted=3, blast_bin_path=os.path.expanduser('~/bin/blast/bin/blastall'), \
						annot_assembly_table='sequence.annot_assembly', \
						raw_sequence_table='sequence.raw_sequence', tmp_blast_infname='/tmp/blast_input'):
		"""
		2007-04-30
		"""
		import sys
		from Bio.Blast import NCBIXML, NCBIStandalone
		sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
		from codense.common import get_sequence_segment
		curs.execute("select s.id, s.acc, s.chromosome, s.position, a.gi from %s s, %s a\
			where a.tax_id=s.tax_id and a.chromosome=s.chromosome"%(snp_locus_table, annot_assembly_table))
		rows = curs.fetchall()
		counter = 0
		outf = open(output_fname, 'w')
		for row in rows:
			snp_locus_id, acc, chromosome, position, genomic_gi = row
			genomic_start = position-flanking_seq_length
			genomic_stop = position + flanking_seq_length
			seq = get_sequence_segment(curs, genomic_gi, genomic_start, genomic_stop, annot_assembly_table, raw_sequence_table)
			if seq:
				for candidate_snp in ['A', 'C', 'G', 'T']:
					seq = seq[:flanking_seq_length] + candidate_snp + seq[flanking_seq_length+1:]	#str doesn't support item assignment
					inf = open(tmp_blast_infname, 'w')
					inf.write('>%s(snp_locus_id=%s) candidate_snp=%s\n'%(acc, snp_locus_id, candidate_snp))
					inf.write('%s\n'%seq)
					del inf
					result_handle, error_info = NCBIStandalone.blastall(blast_bin_path, "blastn", database_fname, tmp_blast_infname, align_view=7)	#align_view=7 toggles the xml output
					blast_parser = NCBIXML.BlastParser()	#the parser has to be re-instanced every time otherwise the blast records parsed in the previous round stay in the parser
					blast_records = blast_parser.parse(result_handle)
					del blast_parser
					no_of_hits = min(max_no_of_hits_to_be_outputted, len(blast_records.alignments))
					outf.write("%s (id=%s) candidate_snp=%s\n"%(acc, snp_locus_id, candidate_snp))
					outf.write("query sequence = %s\n"%seq)
					for i in range(no_of_hits):
						outf.write("\t%s\n"%(blast_records.alignments[i].title))
						for hsp in blast_records.alignments[i].hsps:
							outf.write("\t\tquery_start=%s, sbjct_start=%s, frame=%s, identities=%s/%s E=%s\n"%(hsp.query_start, hsp.sbjct_start, hsp.frame, hsp.identities, 2*flanking_seq_length+1, hsp.expect))
							outf.write("\t\t%s\n"%hsp.query)
							outf.write("\t\t%s\n"%hsp.match)
							outf.write("\t\t%s\n"%hsp.sbjct)
							outf.write("\n")
						outf.write("\n")
					outf.write("\n")
			sys.stderr.write("%s%s"%('\x08'*10, counter))
			counter += 1
		del outf
	
	"""
	my_blast_db = os.path.expanduser("~/bin/blast/db/Arabidopsis_thaliana.main_genome.fasta")
	DB149.blast_snp_segment(curs, 'snp_locus', './blast_149snps_vs_thaliana_len_25.txt', my_blast_db)
	
	
	my_blast_db = os.path.expanduser("~/bin/blast/db/Arabidopsis_lyrata.main_genome.scaffolds.fasta")
	DB149.blast_snp_segment(curs, 'snp_locus', './blast_149snps_vs_lyrata_len_51.txt', my_blast_db, flanking_seq_length=25)
	"""
	
	"""
	2007-04-30
	"""
	def fill_snp_locus_table_with_25mer_thaliana_call(curs, snp_locus_table, annot_assembly_table='sequence.annot_assembly', \
		raw_sequence_table='sequence.raw_sequence', need_commit=0):
		import sys
		sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
		from codense.common import get_sequence_segment
		curs.execute("select s.id, s.acc, s.chromosome, s.position, a.gi from %s s, %s a\
			where a.tax_id=s.tax_id and a.chromosome=s.chromosome"%(snp_locus_table, annot_assembly_table))
		rows = curs.fetchall()
		counter = 0
		for row in rows:
			snp_locus_id, acc, chromosome, position, genomic_gi = row
			genomic_start = position - 12
			genomic_stop = position + 12
			seq = get_sequence_segment(curs, genomic_gi, genomic_start, genomic_stop, annot_assembly_table, raw_sequence_table)
			if seq:
				thaliana_call = seq[12]
				curs.execute("update %s set thaliana_call='%s', flanking_25mer='%s' where id=%s"%(snp_locus_table, thaliana_call, seq, snp_locus_id))
			sys.stderr.write("%s%s"%('\x08'*10, counter))
			counter += 1
		if need_commit:
			curs.execute("end")
	
	"""
	fill_snp_locus_table_with_25mer_thaliana_call(curs, 'dbsnp.snp_locus', need_commit=1)
	"""
	
	
	"""
	2007-03-05
		add the position info given by Yan Li from Borevitz Lab
	"""
	def fill_snp_locus_table_with_position_info(curs, input_fname, output_table, need_commit):
		import csv
		reader = csv.reader(open(input_fname))
		snp_locus_acc_list = reader.next()[1:]
		snp_locus_position_list = reader.next()[1:]
		for i in range(len(snp_locus_acc_list)):
			snp_locus_acc = snp_locus_acc_list[i]
			snp_locus_acc = snp_locus_acc.replace('.', ' ')
			snp_locus_position = snp_locus_position_list[i]
			chromosome, position = snp_locus_position.split('-')
			curs.execute("update %s set chromosome='%s', position=%s where acc='%s'"%(output_table, chromosome, position, snp_locus_acc))
		if need_commit:
			curs.execute("end")
	
	
	"""
	fill_snp_locus_table_with_position_info(curs, './script/variation/data/snp_position.csv', 'dbsnp.snp_locus', need_commit=1)
	"""
	
	"""
	2007-09-21
		check duplicate calls in the database
	"""
	def check_inconsistent_duplicate_calls(curs, ecotype_table, calls_table, debug=0):
		"""
		2007-09-22
			buggy, don't use it. use check_inconsistent_duplicate_calls2()
		"""
		sys.stderr.write("Checking inconsistent duplicate calls ...\n")
		from common import nt2number
		no_of_strain_snp_pairs = 0
		no_of_inconsistent = 0
		inconsistency = 0
		old_call_number = 0
		old_strain_snp_pair = None
		curs.execute("select e.nativename, c.snpid, c.call1, c.callhet from %s e, %s c where e.id=c.ecotypeid order by nativename, snpid"%(ecotype_table, calls_table))
		rows = curs.fetchall()
		counter = 0
		from sets import Set
		strain_snp_pair_set = Set()
		inconsistent_dup_strain_snp_pair_set = Set()
		if debug:
			import pdb
			pdb.set_trace()
		for row in rows:
			nativename, snpid, call, callhet = row
			nativename = nativename.upper()	#bug here. same name appears in >1 forms differing in cases
			call = call.upper()
			if callhet:
				callhet.upper()
				call = call+callhet
			call_number = nt2number[call]
			strain_snp_pair = (nativename, snpid)
			if strain_snp_pair != old_strain_snp_pair:
				no_of_strain_snp_pairs += 1
				no_of_inconsistent += inconsistency
				if inconsistency == 1:
					inconsistent_dup_strain_snp_pair_set.add(old_strain_snp_pair)
				old_strain_snp_pair = strain_snp_pair
				inconsistency = 0
			else:
				if inconsistency == 0 and old_call_number!=0 and call_number!=0 and old_call_number!=call_number:	#if inconsistency==1, don't change it. comparison only between same pairs
					inconsistency = 1
			old_call_number = call_number
			counter += 1
			strain_snp_pair_set.add(strain_snp_pair)
		#don't miss out the last one
		no_of_inconsistent += inconsistency
		old_strain_snp_pair = strain_snp_pair
		sys.stderr.write('%s%s'%('\x08'*20, counter))
		sys.stderr.write("\nDone.\n")
		print "%s strain snp pairs"%(len(strain_snp_pair_set))
		print 'inconsistent ratio: %s/%s=%s'%(no_of_inconsistent, no_of_strain_snp_pairs, float(no_of_inconsistent)/no_of_strain_snp_pairs)
		return strain_snp_pair_set, inconsistent_dup_strain_snp_pair_set
	
	def check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=1, debug=0):
		"""
		use_ecotypeid_as_strainname controls whether it's nativename or ecotypeid in strain_snp_pair
		"""
		sys.stderr.write("Checking inconsistent duplicate calls ...\n")
		if strainname_type==1:
			print 'distinct on (nativename,snpid) pair:'
		elif strainname_type == 2:
			print 'distinct on (ecotypeid,snpid) pair:'
		elif strainname_type == 3:
			print 'distinct on (nativename, stockparent ,snpid) pair:'
		else:
			sys.stderr.write("unsupported strainname_type %s\n"%strainname_type)
			return None, None
		from common import nt2number
		no_of_strain_snp_pairs = 0
		no_of_inconsistent = 0
		old_strain_snp_pair = None
		duplicated_times = 0
		if strainname_type == 1:
			curs.execute("select e.nativename, c.snpid, c.call1, c.callhet from %s e, %s c where e.id=c.ecotypeid order by nativename, snpid"%(ecotype_table, calls_table))
		elif strainname_type == 2:
			curs.execute("select e.id, c.snpid, c.call1, c.callhet from %s e, %s c where e.id=c.ecotypeid order by nativename, snpid"%(ecotype_table, calls_table))
		elif strainname_type == 3:
			curs.execute("select e.nativename, e.stockparent, c.snpid, c.call1, c.callhet from %s e, %s c where e.id=c.ecotypeid order by nativename, stockparent, snpid"%(ecotype_table, calls_table))
		rows = curs.fetchall()
		counter = 0
		from sets import Set
		strain_snp_pair_set = Set()
		inconsistent_dup_strain_snp_pair_set = Set()
		duplicate_call_type2counter = {}
		if debug:
			import pdb
			pdb.set_trace()
		for row in rows:
			if strainname_type==3:
				nativename, stockparent, snpid, call, callhet = row
			elif strainname_type==1:
				nativename, snpid, call, callhet = row
			elif strainname_type==2:
				ecotypeid, snpid, call, callhet = row
			if strainname_type!=2:
				nativename = nativename.upper()	#bug here. same name appears in >1 forms differing in cases
			call = call.upper()
			if callhet:
				callhet.upper()
				call = call+callhet
			call_number = nt2number[call]
			if strainname_type==1:
				strain_snp_pair = (nativename, snpid)
			elif strainname_type==2:
				strain_snp_pair = (ecotypeid, snpid)
			elif strainname_type==3:
				strain_snp_pair = (nativename, stockparent, snpid)
			if old_strain_snp_pair==None:	#first time
				old_strain_snp_pair = strain_snp_pair
				duplicated_times += 1
			elif strain_snp_pair != old_strain_snp_pair:
				if duplicated_times>1:
					no_of_strain_snp_pairs += 1
					if len(duplicate_call_type2counter)>1:
						no_of_inconsistent += 1
						inconsistent_dup_strain_snp_pair_set.add(old_strain_snp_pair)
				old_strain_snp_pair = strain_snp_pair
				duplicated_times = 1	#the current position is counted as 1
				duplicate_call_type2counter = {}
			else:
				duplicated_times += 1
			if call_number!=0:	#have to be put last
				if call_number not in duplicate_call_type2counter:
					duplicate_call_type2counter[call_number] = 0
				duplicate_call_type2counter[call_number] += 1
			counter += 1
			strain_snp_pair_set.add(strain_snp_pair)
		#don't miss out the last one
		if duplicated_times>1:
			no_of_strain_snp_pairs += 1
			if len(duplicate_call_type2counter)>1:
				no_of_inconsistent += 1
		old_strain_snp_pair = strain_snp_pair
		print "%s strain snp pairs"%counter
		print "%s distinct strain snp pairs"%(len(strain_snp_pair_set))
		print 'inconsistent ratio: %s/%s=%s'%(no_of_inconsistent, no_of_strain_snp_pairs, float(no_of_inconsistent)/no_of_strain_snp_pairs)
		sys.stderr.write("Done.\n")
		return strain_snp_pair_set, inconsistent_dup_strain_snp_pair_set
	
	"""
		#ecotype_table = 'ecotype'
		#calls_table = 'calls'
		#strain_snp_pair_set1, inconsistent_dup_strain_snp_pair_set1 = DB149.check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table)
		#strain_snp_pair_set2, inconsistent_dup_strain_snp_pair_set2 = DB149.check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=2)
		#strain_snp_pair_set3, inconsistent_dup_strain_snp_pair_set3 = DB149.check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=3)
		
	"""
	
	def get_nativename_snp_distinct_set(curs, ecotype_table, calls_table):
		"""
		2007-09-22 find out the bug caused by the case-insensitive problem of sql.
		"""
		curs.execute("select distinct e.nativename, c.snpid from %s e, %s c where e.id=c.ecotypeid" %(ecotype_table, calls_table))
		from sets import Set
		strain_snp_pair_set = Set()
		rows = curs.fetchall()
		for row in rows:
			nativename, snpid = row
			strain_snp_pair = (nativename, snpid)
			strain_snp_pair_set.add(strain_snp_pair)
		return strain_snp_pair_set
	"""
	ecotype_table = 'ecotype'
	calls_table = 'calls'
	strain_snp_pair_set, inconsistent_dup_strain_snp_pair_set = check_inconsistent_duplicate_calls(curs, ecotype_table, calls_table)
	strain_snp_pair_set1, inconsistent_dup_strain_snp_pair_set1 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table)
	strain_snp_pair_set2, inconsistent_dup_strain_snp_pair_set2 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=2)
	strain_snp_pair_set3, inconsistent_dup_strain_snp_pair_set3 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=3)
	
	nativename_snp_distinct_set = get_nativename_snp_distinct_set(curs, ecotype_table, calls_table)
	
	strain_snp_pair_set - nativename_snp_distinct_set gives 1937 different pairs. turns out that sql's "distinct e.nativename, c.snpid" ignores the cases , like Kz-9 and KZ-9 are same. however python doesn't do that.
	"""

class Data149_PopulationStructure(object):
	#2007-06-17 function to calculate the great circle distance of a sphere given latitude and longitude
	# http://en.wikipedia.org/wiki/Great-circle_distance
	def cal_great_circle_distance(lat1, lon1, lat2, lon2, earth_radius=6372.795):
		import math
		lat1_rad = lat1*math.pi/180
		lon1_rad = lon1*math.pi/180
		lat2_rad = lat2*math.pi/180
		lon2_rad = lon2*math.pi/180
		long_diff = abs(lon1_rad-lon2_rad)
		sin_lat1 = math.sin(lat1_rad)
		cos_lat1 = math.cos(lat1_rad)
		sin_lat2 = math.sin(lat2_rad)
		cos_lat2 = math.cos(lat2_rad)
		spheric_angular_diff = math.atan2(math.sqrt(math.pow(cos_lat2*math.sin(long_diff),2) + math.pow(cos_lat1*sin_lat2-sin_lat1*cos_lat2*math.cos(long_diff),2)), sin_lat1*sin_lat2+cos_lat1*cos_lat2*math.cos(long_diff))
		return earth_radius*spheric_angular_diff
	
	#2007-06-17 draw location of all strains onto a map
	def test_draw_strain_location(pic_area=[-180,-90,180,90]):
		import pylab
		from matplotlib.toolkits.basemap import Basemap
		pylab.clf()
		fig = pylab.figure()
		m = Basemap(llcrnrlon=pic_area[0],llcrnrlat=pic_area[1],urcrnrlon=pic_area[2],urcrnrlat=pic_area[3],\
	            resolution='l',projection='mill')
		#m.drawcoastlines()
		m.drawparallels(pylab.arange(-9,90,30), labels=[1,1,1,1])
		m.drawmeridians(pylab.arange(-180,180,60), labels=[1,1,1,1])
		m.fillcontinents()
		m.drawcountries()
		m.drawstates()
		
		import MySQLdb
		db = MySQLdb.connect(db="stock",host='natural.uchicago.edu', user='iamhere', passwd='iamhereatusc')
		c = db.cursor()
		c.execute("""select latitude, longitude from ecotype where latitude is not null and longitude is not null""")
		lat_lon_ls = c.fetchall()
		euc_coord1_ls = []
		euc_coord2_ls = []
		for lat, lon in lat_lon_ls:
			euc_coord1, euc_coord2 = m(lon, lat)	#longitude first, latitude 2nd
			euc_coord1_ls.append(euc_coord1)
			euc_coord2_ls.append(euc_coord2)
		m.scatter(euc_coord1_ls, euc_coord2_ls, 25, marker='o', zorder=10)
		pylab.title("worldwide distribution of %s strains"%(len(lat_lon_ls)))
		pylab.show()
	
	#2007-06-17 calculate pairwise distance among strains by calling cal_great_circle_distance()
	def cal_pairwise_distance_of_strains():
		import MySQLdb, pylab
		db = MySQLdb.connect(db="stock",host='natural.uchicago.edu', user='iamhere', passwd='iamhereatusc')
		c = db.cursor()
		c.execute("""select latitude, longitude from ecotype where latitude is not null and longitude is not null""")
		lat_lon_ls = c.fetchall()
		distance_ls = []
		no_of_sites = len(lat_lon_ls)
		for i in range(no_of_sites):
			for j in range(i+1, no_of_sites):
				distance_ls.append(cal_great_circle_distance(lat_lon_ls[i][0], lat_lon_ls[i][1], lat_lon_ls[j][0], lat_lon_ls[j][1]))
		pylab.clf()
		pylab.hist(distance_ls, 20)
		pylab.title("histogram of distances between %s strains(km)"%no_of_sites)
		pylab.show()
	
	"""
	#2007-06-18, calll cal_great_circle_distance()
	import MySQLdb
	conn = MySQLdb.connect(db="stock",host='natural.uchicago.edu', user='iamhere', passwd='iamhereatusc')
	cursor = conn.cursor()
	cursor.execute("select latitude, longitude from ecotype where latitude is not null and longitude is not null")
	lat_lon_ls = cursor.fetchall()
	"""
	
	
	"""
	2007-09-13
	following functions to investigate inter/intra-population identity
	"""
	def construct_identity_pair_ls(input_fname):
		"""
		2007-09-13
			input_fname is a data frame file with snps coded in integers. 1st row is the ecotype id in table ecotype
		"""
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
		no_of_strains = len(strain_acc_list)
		no_of_snps = len(header)-2
		identity_pair_ls = []
		for i in range(no_of_strains):
			for j in range(i+1, no_of_strains):
				no_of_same_cols = 0
				for k in range(no_of_snps):
					if data_matrix[i][k] == data_matrix[j][k] or data_matrix[i][k]==0 or data_matrix[j][k]==0:
						no_of_same_cols += 1
				if no_of_same_cols == no_of_snps:
					identity_pair_ls.append([strain_acc_list[i], strain_acc_list[j]])
		return identity_pair_ls
	
	
	def construct_graph_out_of_identity_pair(identity_pair_ls):
		"""
		2007-09-18
			construct a graph based on identity pairs in identity_pair_ls
			identity_pair_ls's entries will be converted to integers
		"""
		sys.stderr.write("Constructing graph out of identity_pair_ls ...")
		import networkx as nx
		g = nx.XGraph()
		for identity_pair in identity_pair_ls:
			identity_pair = map(int, identity_pair)
			g.add_edge(identity_pair[0], identity_pair[1], 1)
		sys.stderr.write("Done.\n")
		return g
	
	def get_singleton_strain_id_ls(g, input_fname):
		"""
		2007-10-01
			node in g is in ecotype id format
			so strain_acc (1st column) in input_fname should be integer (ecotype id)
		"""
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
		sys.stderr.write("Getting singleton_strain_id_ls ...")
		from sets import Set
		graph_node_set = Set(g.nodes())
		singleton_strain_id_ls = []
		for strain_acc in strain_acc_list:
			strain_acc = int(strain_acc)
			if strain_acc not in graph_node_set:
				singleton_strain_id_ls.append(strain_acc)
		sys.stderr.write("Done.\n")
		return singleton_strain_id_ls
	
	def construct_site_graph_out_of_strain_graph(g, ecotypeid2pos, node_weight_by_no_of_nodes=0):
		"""
		2007-09-20
			#reduce the drawing by omitting identical edges (linking same GPS positions), construct a graph whose node is (lat,lon) and then plot this graph.
		2007-10-01
			add node_weight_by_no_of_nodes
		"""
		sys.stderr.write("Constructing site graph out of strain graph...")
		site2pos = {}
		site2weight = {}
		import networkx as nx
		site_g = nx.XGraph()
		for e in g.edges():
			pos1 = (round(ecotypeid2pos[e[0]][0],2), round(ecotypeid2pos[e[0]][1],2))
			pos2 = (round(ecotypeid2pos[e[1]][0],2), round(ecotypeid2pos[e[1]][1],2))
			if not node_weight_by_no_of_nodes:	#2007-10-01
				if pos1 not in site2weight:
					site2weight[pos1] = 0
				site2weight[pos1] += 1
				if pos2 not in site2weight:
					site2weight[pos2] = 0
				site2weight[pos2] += 1
			if pos1==pos2:
				site_g.add_node(pos1)
			else:
				if not site_g.has_edge(pos1, pos2):
					site_g.add_edge(pos1, pos2, 0)
				site_g.adj[pos1][pos2] += 1
				site_g.adj[pos2][pos1] += 1
		if node_weight_by_no_of_nodes:	#2007-10-01
			for n in g.nodes():
				pos = (round(ecotypeid2pos[n][0],2), round(ecotypeid2pos[n][1],2))
				if pos not in site2weight:
					site2weight[pos] = 0
				site2weight[pos] += 1
		for n in site_g:
			site2pos[n] = n
		sys.stderr.write("Done.\n")
		return site_g, site2weight, site2pos
	
	def collapseStrainsWithSamePos(ecotypeid_ls, ecotypeid2pos):
		"""
		2007-10-09
		"""
		sys.stderr.write("Collapsing Strains with same GPS ...")
		site2pos = {}
		site2weight = {}
		for ecotypeid in ecotypeid_ls:
			pos = (round(ecotypeid2pos[ecotypeid][0],2), round(ecotypeid2pos[ecotypeid][1],2))
			if pos not in site2weight:
				site2weight[pos] = 0
				site2pos[pos] = pos
			site2weight[pos] += 1
		sys.stderr.write("Done.\n")
		return site2weight, site2pos
	
	
	def draw_strains_on_map(ecotypeid_ls, ecotypeid2pos, pic_title,  pic_area=[-130,10,140,70], output_fname_prefix=None, label_ls=None, need_collapse_strains_with_same_pos=0):
		"""
		2007-10-09
		2007-11-08
			add label_ls, need_collapse_strains_with_same_pos
		"""
		import os, sys
		sys.stderr.write("Drawing strains on a map ...\n")
		import pylab, math
		from matplotlib.toolkits.basemap import Basemap
		pylab.clf()
		fig = pylab.figure()
		fig.add_axes([0.05,0.05,0.9,0.9])	#[left, bottom, width, height]
		m = Basemap(llcrnrlon=pic_area[0],llcrnrlat=pic_area[1],urcrnrlon=pic_area[2],urcrnrlat=pic_area[3],\
		resolution='l',projection='mill')
		
		
		sys.stderr.write("\tDrawing nodes ...")
		euc_coord1_ls = []
		euc_coord2_ls = []
		diameter_ls = []
		ax=pylab.gca()
		if need_collapse_strains_with_same_pos:
			site2weight, site2pos = collapseStrainsWithSamePos(ecotypeid_ls, ecotypeid2pos)
			for n in site2pos:
				lat, lon = site2pos[n]
				euc_coord1, euc_coord2 = m(lon, lat)	#longitude first, latitude 2nd
				euc_coord1_ls.append(euc_coord1)
				euc_coord2_ls.append(euc_coord2)
				diameter_ls.append(8*math.sqrt(site2weight[n]))
		else:
			for i in range(len(ecotypeid_ls)):
				ecotypeid = ecotypeid_ls[i]
				lat, lon = ecotypeid2pos[ecotypeid]
				euc_coord1, euc_coord2 = m(lon, lat)
				euc_coord1_ls.append(euc_coord1)
				euc_coord2_ls.append(euc_coord2)
				if label_ls:
					ax.text(euc_coord1, euc_coord2, str(label_ls[i]), size=2, alpha=0.5, horizontalalignment='center', verticalalignment='center', zorder=12)
				diameter_ls.append(3)
		import numpy
		#diameter_ls = numpy.array(diameter_ls)
		m.scatter(euc_coord1_ls, euc_coord2_ls, diameter_ls, marker='o', color='r', alpha=0.4, zorder=10, faceted=False)
		sys.stderr.write("Done.\n")
		
		#m.drawcoastlines()
		m.drawparallels(pylab.arange(-90,90,30), labels=[1,1,0,1])
		m.drawmeridians(pylab.arange(-180,180,30), labels=[1,1,0,1])
		m.fillcontinents()
		m.drawcountries()
		m.drawstates()
		
		pylab.title(pic_title)
		if output_fname_prefix:
			pylab.savefig('%s.eps'%output_fname_prefix, dpi=600)
			#pylab.savefig('%s.svg'%output_fname_prefix, dpi=600)
			pylab.savefig('%s.png'%output_fname_prefix, dpi=600)
		sys.stderr.write("Done.\n")
	
	def display_matrix_of_component(input_fname, ecotypeid_ls, ecotypeid2pos, output_fname=None, need_sort=0, need_savefig=0):
		"""
		2007-09-20
			display the data from that component
			
		"""
		import csv, Numeric, pylab
		cc_ecotypeid_pos = []
		for ecotypeid in ecotypeid_ls:
			cc_ecotypeid_pos.append(ecotypeid2pos[ecotypeid])
		import numpy
		argsort_index = numpy.argsort(cc_ecotypeid_pos, 0)	#watch it's two dimensional
		ecotypeid2row_index = {}
		ytick_label_ls = []
		cc_size = len(ecotypeid_ls)
		for i in range(cc_size):
			ecotypeid_ls_index = argsort_index[i][1]	#sort based on longitude
			ecotypeid = ecotypeid_ls[ecotypeid_ls_index]
			ecotypeid2row_index[ecotypeid] = i
			ytick_label_ls.append('%s (%.2f, %.2f)'%(ecotypeid, ecotypeid2pos[ecotypeid][0], ecotypeid2pos[ecotypeid][1]))
		reader = csv.reader(open(input_fname), delimiter='\t')
		header = reader.next()
		data_matrix = [0]*cc_size
		for row in reader:
			ecotypeid = int(row[0])
			if ecotypeid in ecotypeid2row_index:
				data_row = row[2:]
				data_row = map(int, data_row)
				data_matrix[ecotypeid2row_index[ecotypeid]] = data_row
		del reader
		data_matrix.reverse()	#2007-03-06 reverse() due to the imshow()'s y axis starting from bottom
		if need_sort:
			data_matrix.sort()
		data_matrix = Numeric.array(data_matrix)
		
		pylab.clf()
		pylab.imshow(data_matrix, aspect='auto', interpolation='nearest')	#2007-06-05
		pylab.colorbar()
		pylab.yticks(range(cc_size), ytick_label_ls)
		if need_savefig:
			pylab.savefig('%s.eps'%output_fname, dpi=300)
			pylab.savefig('%s.svg'%output_fname, dpi=300)
			pylab.savefig('%s.png'%output_fname, dpi=300)
		pylab.show()
	
	def get_ecotypeid2pos(curs, ecotype_table, with_data_affiliated=0):
		"""
		2007-09-18
		2007-10-09
			add with_data_affiliated
		"""
		sys.stderr.write("Getting ecotypeid2pos from %s ..."%ecotype_table)
		ecotypeid2pos = {}
		if with_data_affiliated:
			curs.execute("select distinct e.id, e.latitude, e.longitude from %s e, calls c where e.latitude is not null and e.longitude is not null and e.id=c.ecotypeid"%ecotype_table)
		else:
			curs.execute("select id, latitude, longitude from %s where latitude is not null and longitude is not null"%ecotype_table)
		rows = curs.fetchall()
		for row in rows:
			ecotypeid, latitude, longitude = row
			ecotypeid2pos[ecotypeid] = (latitude, longitude)
		sys.stderr.write("Done.\n")
		return ecotypeid2pos
	
	def get_popid2pos_size_and_ecotypeid2pop_id(curs, popid2ecotypeid_table, selected=0):
		"""
		2007-09-13
		2007-09-24
			doesn't have to be selected
		"""
		sys.stderr.write("Getting popid2pos_size & ecotypeid2pop_id ...")
		popid2pos_size = {}
		ecotypeid2pop_id = {}
		curs.execute("select popid, latitude, longitude, ecotypeid from %s where selected=%s"%(popid2ecotypeid_table, selected))
		rows = curs.fetchall()
		for row in rows:
			popid, latitude, longitude, ecotypeid = row
			if popid not in popid2pos_size:
				popid2pos_size[popid] = [None, 0]
			if popid2pos_size[popid][0]==None:
				popid2pos_size[popid][0] = (latitude, longitude)
			popid2pos_size[popid][1] += 1
			ecotypeid2pop_id[ecotypeid] = popid
		sys.stderr.write("Done.\n")
		return popid2pos_size, ecotypeid2pop_id
	
	def construct_pop_identity_graph_from_strain_identity_ls(curs, popid2ecotypeid_table, identity_pair_ls):
		"""
		2007-09-18
		"""
		import os, sys
		popid2pos_size, ecotypeid2pop_id = get_popid2pos_size_and_ecotypeid2pop_id(curs, popid2ecotypeid_table)
		
		no_of_inter_pop_identities = 0
		no_of_valid_identities = 0	#valid means both strains belong to some population
		import networkx as nx
		g = nx.XGraph()
		popid2intra_pop_connections = {}
		for e in identity_pair_ls:
			e = map(int, e)
			popid1 = ecotypeid2pop_id.get(e[0])
			popid2 = ecotypeid2pop_id.get(e[1])
			if popid1 and popid2:	#both are classified into populations
				no_of_valid_identities += 1
				if popid1 not in popid2intra_pop_connections:
					popid2intra_pop_connections[popid1] = 0
					g.add_node(popid1)
				if popid2 not in popid2intra_pop_connections:
					popid2intra_pop_connections[popid2] = 0
					g.add_node(popid2)
				if popid1!=popid2:
					if not g.has_edge(popid1, popid2):
						g.add_edge(popid1, popid2, 0)
					g.adj[popid1][popid2] += 1	#increase the edge weight
					g.adj[popid2][popid1] += 1
					no_of_inter_pop_identities += 1
				else:
					popid2intra_pop_connections[popid1] += 1
		print '%s/%s=%s inter population identities'%(no_of_inter_pop_identities, no_of_valid_identities, float(no_of_inter_pop_identities)/no_of_valid_identities)
		return g, popid2intra_pop_connections, popid2pos_size, ecotypeid2pop_id
	
	def get_pic_area(pos_ls):
		"""
		2007-10-02
			get a good pic_area
		"""
		min_lat = pos_ls[0][0]
		max_lat = pos_ls[0][0]
		min_lon = pos_ls[0][1]
		max_lon = pos_ls[0][1]
		for i in range(1, len(pos_ls)):
			pos = pos_ls[i]
			if pos[0]<min_lat:
				min_lat = pos[0]
			if pos[0]>max_lat:
				max_lat = pos[0]
			if pos[1]<min_lon:
				min_lon = pos[1]
			if pos[1]>max_lon:
				max_lon = pos[1]
		pic_area = [-180, -90, 180, 90]	#this is the boundary
		if min_lon-5>pic_area[0]:
			pic_area[0] = min_lon-5
		if min_lat-5>pic_area[1]:
			pic_area[1] = min_lat-5
		if max_lon+5<pic_area[2]:
			pic_area[2] = max_lon+5
		if max_lat+5<pic_area[3]:
			pic_area[3] = max_lat+5
		return pic_area
	
	def draw_graph_on_map(g, node2weight, node2pos, pic_title,  pic_area=[-130,10,140,70], output_fname_prefix=None, need_draw_edge=0):
		"""
		2007-09-13
			identity_pair_ls is a list of pairs of strains (ecotype id as in table ecotype)
		2007-10-08
			correct a bug in 4*diameter_ls, diameter_ls has to be converted to array first.
			sqrt the node weight, 8 times the original weight
		"""
		import os, sys
		sys.stderr.write("Drawing graph on a map ...\n")
		import pylab, math
		from matplotlib.toolkits.basemap import Basemap
		pylab.clf()
		fig = pylab.figure()
		fig.add_axes([0.05,0.05,0.9,0.9])	#[left, bottom, width, height]
		m = Basemap(llcrnrlon=pic_area[0],llcrnrlat=pic_area[1],urcrnrlon=pic_area[2],urcrnrlat=pic_area[3],\
		resolution='l',projection='mill')
		
		sys.stderr.write("\tDrawing nodes ...")
		euc_coord1_ls = []
		euc_coord2_ls = []
		diameter_ls = []
		for n in g.nodes():
			lat, lon = node2pos[n]
			euc_coord1, euc_coord2 = m(lon, lat)	#longitude first, latitude 2nd
			euc_coord1_ls.append(euc_coord1)
			euc_coord2_ls.append(euc_coord2)
			diameter_ls.append(math.sqrt(node2weight[n]))
		import numpy
		diameter_ls = numpy.array(diameter_ls)
		m.scatter(euc_coord1_ls, euc_coord2_ls, 8*diameter_ls, marker='o', color='r', alpha=0.4, zorder=12, faceted=False)
		sys.stderr.write("Done.\n")
		
		if need_draw_edge:
			sys.stderr.write("\tDrawing edges ...")
			ax=pylab.gca()
			for popid1, popid2, no_of_connections in g.edges():
				lat1, lon1 = node2pos[popid1]
				lat2, lon2 = node2pos[popid2]
				x1, y1 = m(lon1, lat1)
				x2, y2 = m(lon2, lat2)
				ax.plot([x1,x2],[y1,y2], 'g', linewidth=math.log(no_of_connections+1)/2, alpha=0.2, zorder=10)
			sys.stderr.write("Done.\n")
		
		#m.drawcoastlines()
		m.drawparallels(pylab.arange(-90,90,30), labels=[1,1,0,1])
		m.drawmeridians(pylab.arange(-180,180,30), labels=[1,1,0,1])
		m.fillcontinents()
		m.drawcountries()
		m.drawstates()
		
		pylab.title(pic_title)
		if output_fname_prefix:
			pylab.savefig('%s_site_network.eps'%output_fname_prefix, dpi=600)
			pylab.savefig('%s_site_network.svg'%output_fname_prefix, dpi=600)
			pylab.savefig('%s_site_network.png'%output_fname_prefix, dpi=600)
		sys.stderr.write("Done.\n")
	
	def check_cross_ocean_components(g, cc, ecotypeid2pos, longi_watershed=-30):
		from sets  import Set
		cross_ocean_cc2lon_pair_set = {}
		for i in range(len(cc)):
			gs = g.subgraph(cc[i])
			for e in gs.edges():
				l1 = ecotypeid2pos[e[0]][1]
				l2 = ecotypeid2pos[e[1]][1]
				if (l1<longi_watershed and l2>longi_watershed) or (l1>longi_watershed and l2<longi_watershed):
					if i not in cross_ocean_cc2lon_pair_set:
						cross_ocean_cc2lon_pair_set[i] = Set()
					cross_ocean_cc2lon_pair_set[i].add((l1, l2))
		return cross_ocean_cc2lon_pair_set
	
	def get_clique_ls_indexed_by_cc_number(g):
		"""
		2007-10-01 each component is partitioned into cliques
			the resulting clique_ls is put in the order of connected components
			like clique_ls_indexed_by_cc_number[0] is clique_ls from connected component 0
		"""
		sys.stderr.write("Getting clique_ls_indexed_by_cc_number ...\n")
		import networkx as nx
		from pymodule.PartitionGraphIntoCliques import PartitionGraphIntoCliques
		PartitionGraphIntoCliques_ins = PartitionGraphIntoCliques(algorithm_type=2)
		g_cc_ls = nx.connected_components(g)
		clique_ls_indexed_by_cc_number = []
		for i in range(len(g_cc_ls)):
			sys.stderr.write("%s\t%s"%('\x08'*20, i))
			cc_g = g.subgraph(g_cc_ls[i])
			PartitionGraphIntoCliques_ins.partition(cc_g.copy())
			clique_ls_indexed_by_cc_number.append(PartitionGraphIntoCliques_ins.clique_ls)
		sys.stderr.write("Done.\n")
		return g_cc_ls, clique_ls_indexed_by_cc_number
	
	def get_haplotype_size2number(clique_ls_indexed_by_cc_number):
		"""
		2007-10-02
		"""
		sys.stderr.write("Getting haplotype_size2number ...")
		haplotype_size2number = {}
		for i in range(len(clique_ls_indexed_by_cc_number)):
			for clique in clique_ls_indexed_by_cc_number[i]:
				haplotype_size = len(clique)
				if haplotype_size not in haplotype_size2number:
					haplotype_size2number[haplotype_size] = 0
				haplotype_size2number[haplotype_size] += 1
		sys.stderr.write("Done.\n")
		return haplotype_size2number
	
	def get_ordered_clique_ls(clique_ls_indexed_by_cc_number):
		"""
		2007-10-02
			return cliques with size in descending order.
		"""
		sys.stderr.write("Getting ordered_clique_ls ...")
		clique_ls = []
		ordered_clique_ls = []
		import numpy
		for i in range(len(clique_ls_indexed_by_cc_number)):
			for clique in clique_ls_indexed_by_cc_number[i]:
				clique_ls.append(clique)
		clique_ls_argsort = list(numpy.argsort(map(len, clique_ls)))
		clique_ls_argsort.reverse()	#from big to small
		for index in clique_ls_argsort:
			ordered_clique_ls.append(clique_ls[index])
		sys.stderr.write("Done.\n")
		return ordered_clique_ls
	
	def get_strains_within_longitude_span(ecotypeid_ls, ecotypeid2pos, longitude_span=[-12, 50]):
		"""
		2007-10-04
			check the map roughly
			[-12,50] is europe
			[-130,-60] is north america
			[60,90] is central asia
			[130, 150] is japan
			
		"""
		no_of_strains = 0
		strain_ls = []
		for ecotypeid in ecotypeid_ls:
			pos = ecotypeid2pos[ecotypeid]
			if pos[1]>=longitude_span[0] and pos[1]<=longitude_span[1]:
				no_of_strains += 1
				strain_ls.append(ecotypeid)
		print "#strains with %s is %s"%(longitude_span, no_of_strains)
		return strain_ls
	
	def draw_histogram_of_haplotype_size(haplotype_size2number, output_fname_prefix=None):
		"""
		2007-10-02
		"""
		import pylab
		pylab.clf()
		haplotype_size_ls = []
		count_ls = []
		for haplotype_size, count in haplotype_size2number.iteritems():
			haplotype_size_ls.append(haplotype_size)
			count_ls.append(count)
		pylab.title("haplotype_size bar chart")
		pylab.xlabel("haplotype_size")
		pylab.ylabel("count")
		pylab.bar(haplotype_size_ls, count_ls)
		if output_fname_prefix:
			pylab.savefig('%s.eps'%output_fname_prefix, dpi=300)
			pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
			pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		pylab.show()
	
	def calGenoDistAndGeoDist(data_matrix_fname, ecotypeid2pos, longitude_span=[-180, 180], dist_range=[0.0,1.0]):
		"""
		2008-01-30
			add dist_range to restrict the output pairs within that genotype distance
		2008-01-29
			return strain_pair_ls as well
			add a general cutoff
		2007-10-03
			plot to see relationship between genotype distance and geographic distance
			-cal_great_circle_distance()
		2007-10-04
			add longitude_span
		"""
		sys.stderr.write("Calculating geno, geo distance ...\n")
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(data_matrix_fname)
		import Numeric
		strain_pair_ls = []
		geno_dist_ls = []
		geo_dist_ls = []
		no_of_NA_pairs = 0
		no_of_strains = len(data_matrix)
		no_of_snps = len(data_matrix[0])
		no_of_pairs = 0
		no_of_included_pairs = 0
		for i in range(no_of_strains):
			for j in range(i+1, no_of_strains):
				pos1 = ecotypeid2pos[int(strain_acc_list[i])]
				pos2 = ecotypeid2pos[int(strain_acc_list[j])]
				if pos1[1]>=longitude_span[0] and pos1[1]<=longitude_span[1] and pos2[1]>=longitude_span[0] and pos2[1]<=longitude_span[1]:
					no_of_pairs += 1
					no_of_valid_pairs = 0.0
					no_of_matching_pairs = 0.0
					for k in range(no_of_snps):
						if data_matrix[i][k]!=0 and data_matrix[j][k]!=0:
							no_of_valid_pairs += 1
							if data_matrix[i][k] == data_matrix[j][k]:
								no_of_matching_pairs += 1
					if no_of_valid_pairs!=0:
						distance = 1 - no_of_matching_pairs/no_of_valid_pairs
						if distance>=dist_range[0] and distance<=dist_range[1]:
							geno_dist_ls.append(distance)
							geo_dist = cal_great_circle_distance(pos1[0], pos1[1], pos2[0], pos2[1])
							geo_dist_ls.append(geo_dist)
							strain_pair_ls.append([int(strain_acc_list[i]), int(strain_acc_list[j])])
							no_of_included_pairs += 1
					else:
						no_of_NA_pairs += 1
		print "out of %s pairs, %s are NA, %s pairs' distance is within range(%s)."%(no_of_pairs, no_of_NA_pairs, no_of_included_pairs, dist_range)
		return strain_pair_ls, geno_dist_ls, geo_dist_ls
	
	def calGenoDistAndGeoDistBetweenTwoAreas(data_matrix_fname, ecotypeid2pos, longitude_span1=[-180, 180], longitude_span2=[-180, 180]):
		"""
		2007-10-05
			modified from calGenoDistAndGeoDist with two longitude spans
		"""
		sys.stderr.write("Calculating geno, geo distance ...\n")
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(data_matrix_fname)
		import Numeric
		geno_dist_ls = []
		geo_dist_ls = []
		no_of_NA_pairs = 0
		no_of_strains = len(data_matrix)
		no_of_snps = len(data_matrix[0])
		no_of_pairs = 0
		for i in range(no_of_strains):
			for j in range(i+1, no_of_strains):
				pos1 = ecotypeid2pos[int(strain_acc_list[i])]
				pos2 = ecotypeid2pos[int(strain_acc_list[j])]
				if (pos1[1]>=longitude_span1[0] and pos1[1]<=longitude_span1[1] and pos2[1]>=longitude_span2[0] and pos2[1]<=longitude_span2[1]) or (pos2[1]>=longitude_span1[0] and pos2[1]<=longitude_span1[1] and pos1[1]>=longitude_span2[0] and pos1[1]<=longitude_span2[1]):
					no_of_pairs += 1
					no_of_valid_pairs = 0.0
					no_of_matching_pairs = 0.0
					for k in range(no_of_snps):
						if data_matrix[i][k]!=0 and data_matrix[j][k]!=0:
							no_of_valid_pairs += 1
							if data_matrix[i][k] == data_matrix[j][k]:
								no_of_matching_pairs += 1
					if no_of_valid_pairs!=0:
						distance = 1 - no_of_matching_pairs/no_of_valid_pairs
						geno_dist_ls.append(distance)
						geo_dist = cal_great_circle_distance(pos1[0], pos1[1], pos2[0], pos2[1])
						geo_dist_ls.append(geo_dist)
					else:
						no_of_NA_pairs += 1
		print "out of %s pairs, %s are NA"%(no_of_pairs, no_of_NA_pairs)
		return geno_dist_ls, geo_dist_ls
	
	def plotHist(ls, title=None, output_fname_prefix=None, no_of_bins=40):
		import pylab
		pylab.clf()
		pylab.hist(ls, no_of_bins)
		if title:
			pylab.title(title)
		if output_fname_prefix:
			pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
			pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		pylab.show()
	
	def plotXY(xls, yls, title=None, xlabel=None, ylabel=None, output_fname_prefix=None):
		import pylab
		pylab.clf()
		
		pylab.title(title)
		pylab.xlabel(xlabel)
		pylab.ylabel(ylabel)
		pylab.plot(xls, yls, '.')
		if output_fname_prefix:
			pylab.savefig('%s.eps'%output_fname_prefix, dpi=300)
			pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
			pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		pylab.show()
	
	def outputGraphInCliquerFormat(graph, output_fname):
		"""
		2007-10-05
			to test if cliquer suffers the same problem
		"""
		of = open(output_fname, 'w')
		of.write("p edge %s %s\n"%(graph.number_of_nodes(), graph.number_of_edges()))
		node_name2index = {}
		for n in graph.nodes():
			node_name2index[n] = len(node_name2index)+1
		for e in graph.edges():
			of.write("e %s %s\n"%(node_name2index[e[0]], node_name2index[e[1]]))
		del of
	
	def get_ecotypeid_ls_given_popid(curs, popid2ecotypeid_table, popid):
		ecotypeid_ls = []
		curs.execute("select ecotypeid from %s where popid=%s and selected=1"%(popid2ecotypeid_table, popid))
		rows = curs.fetchall()
		for row in rows:
			ecotypeid_ls.append(row[0])
		return ecotypeid_ls
	
	"""
	2008-01-25
		following functions to investigate how identity spreads against distance in different regions.
	"""
	def find_span_index_out_of_span_ls(pos, longitude_span_ls):
		"""
		2008-01-25
			find the index of the longitude_span in which pos drops in
		"""
		span_index = -1
		for i in range(len(longitude_span_ls)):
			if pos[1]>=longitude_span_ls[i][0] and pos[1]<=longitude_span_ls[i][1]:
				span_index = i
		return span_index
	
	def get_span_pair_count_ls(input_fname, ecotypeid2pos, span_ls):
		"""
		2008-01-30
			calculate the number of possible pairs for a specific span type (either single span or two crossing-spans)
			for normalization purpose in group_identity_pair_according_to_region_distance()
		"""
		sys.stderr.write("Getting span_pair_count_ls ...")
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
		strain_acc_list = map(int, strain_acc_list)
		#1st check how many strains in each span
		span_count_ls = [0]*len(span_ls)
		for ecotypeid in strain_acc_list:
			pos = ecotypeid2pos[ecotypeid]
			for i in range(len(span_ls)):
				if pos[1]>=span_ls[i][0] and pos[1]<=span_ls[i][1]:
					span_count_ls[i] += 1
		#2nd calculate the number of pairs in each span
		span_pair_count_ls = []
		for i in range(len(span_count_ls)):
			span_pair_count_ls.append(span_count_ls[i]*span_count_ls[i]/2.0)
		#3rd calculate the number of pairs crossing two spans
		#check all combinations of span_ls
		for i in range(len(span_ls)):
			for j in range(i+1, len(span_ls)):
				span_pair_count_ls.append(span_count_ls[i]*span_count_ls[j])
		sys.stderr.write("Done.\n")
		return span_pair_count_ls
	
	def group_identity_pair_according_to_region_distance(strain_pair_ls, geno_dist_ls, geo_dist_ls, ecotypeid2pos, span_ls, span_label_ls, geo_distance_window_size=100, max_identity_geno_dist=0.0, span_pair_count_ls=[]):
		"""
		2008-01-30
			add span_pair_count_ls to normalize the counts
		2008-01-29
			strain_pair_ls, geno_dist_ls, geo_dist_ls replace identity_pair_ls
			use max_identity_geno_dist to determine whether to call it identity or not
		2008-01-25
			to group the identity_pair_ls into different regions and distance windows.
			geo_distance_window_size is in unit km
		"""
		sys.stderr.write("Grouping identity pairs according to region and geographic distance ...\n")
		import math
		list_of_window_no2count = []
		#1. get the row labels (label for the span type of each identity_pair)
		#2. link the span_index_pair to the index of list_of_window_no2count
		#3. initiate list_of_window_no2count
		row_label_ls =[] #region designation
		span_index_pair2list_index = {}
		for i in range(len(span_ls)):
			row_label_ls.append(span_label_ls[i])
			span_index_pair = (i,i)
			span_index_pair2list_index[span_index_pair] = i
			list_of_window_no2count.append({})
		for i in range(len(span_ls)):
			for j in range(i+1, len(span_ls)):
				row_label_ls.append('%s - %s'%(span_label_ls[i], span_label_ls[j]))
				span_index_pair2list_index[(i,j)] = len(row_label_ls)-1 #both order have the same index
				span_index_pair2list_index[(j,i)] = len(row_label_ls)-1
				list_of_window_no2count.append({})
	
		#record all the counts
		max_window_no = -1 #used to determine the final maxtrix shape
		for i in range(len(strain_pair_ls)):
			strain_pair = strain_pair_ls[i]
			geno_dist = geno_dist_ls[i]
			geo_dist = geo_dist_ls[i]
			if geno_dist<=max_identity_geno_dist:
				pos1 = ecotypeid2pos[strain_pair[0]]
				pos2 = ecotypeid2pos[strain_pair[1]]
				span_index1 = find_span_index_out_of_span_ls(pos1, span_ls)
				span_index2 = find_span_index_out_of_span_ls(pos2, span_ls)
				if span_index1>=0 and span_index2>=0:	#both data are in span_ls
					window_no = math.floor(geo_dist/geo_distance_window_size) #take the floor as the bin number
					if window_no>max_window_no:
						max_window_no = window_no
					list_index = span_index_pair2list_index[(span_index1, span_index2)]
					window_no2count = list_of_window_no2count[list_index]
					if window_no not in window_no2count:
						window_no2count[window_no] = 0
					window_no2count[window_no] += 1
		#transform list_of_window_no2count to matrix
		import numpy
		if max_window_no>=0: #-1 means nothing recorded
			matrix_of_counts_by_region_x_geo_dist = numpy.zeros([len(list_of_window_no2count), max_window_no+1], numpy.float) #+1 because window_no starts from 0
			for i in range(len(list_of_window_no2count)):
				for j in range(max_window_no+1):
					if j in list_of_window_no2count[i]:
						if span_pair_count_ls:
							matrix_of_counts_by_region_x_geo_dist[i][j] = list_of_window_no2count[i][j]/float(span_pair_count_ls[i])	#2008-01-30 normalize
						else:
							matrix_of_counts_by_region_x_geo_dist[i][j] = list_of_window_no2count[i][j]
		else:
			matrix_of_counts_by_region_x_geo_dist = None
		#generate the labels for the columns of matrix_of_counts_by_region_x_geo_dist
		col_label_ls = [] #geo distance window designation
		col_index_ls = [] #the X-axis value
		for i in range(max_window_no+1):
			col_index = (i+0.5)*geo_distance_window_size #right in the middle of the window/bin
			col_index_ls.append(col_index)
			col_label_ls.append('%skm'%int(col_index)) #km is unit
		sys.stderr.write("Done.\n")
		return matrix_of_counts_by_region_x_geo_dist, row_label_ls, col_label_ls, col_index_ls
	
	def draw_table_as_bar_chart(matrix, row_label_ls, col_label_ls, col_index_ls=None, output_fname_prefix=None):
		"""
		2008-01-30
			bars from different rows no longer sit on one another. get them on the same ground and one after another.
		2008-01-25
			originated from http://matplotlib.sourceforge.net/screenshots/table_demo.py
		"""
		sys.stderr.write("Drawing Table as bar char ...")
		import matplotlib
		import pylab
		from pymodule.colours import get_colours
		pylab.clf()
		pylab.axes([0.1, 0.3, 0.8, 0.6])   # leave room below the axes for the table
		# Get some pastel shades for the colours
		colours = get_colours(len(row_label_ls))
		colours = list(colours)
		#colours.reverse()
		rows = len(matrix) #matrix.shape[0] will also do, but in case matrix is a 2D list
		ind = pylab.arange(len(col_label_ls))  # the x locations for the groups
		cellText = []
		width = 0.2    # the width of the bars
		#width = (col_index_ls[1]-col_index_ls[0])/2     # the width of the bars is half the window size
		yoff = pylab.array([0.0] * len(col_label_ls)) # the bottom values for stacked bar chart
		bar_ins_ls = []
		step = 0.32
		for row in xrange(rows):
		    bar_ins = pylab.bar(ind+step*row, matrix[row], width, bottom=yoff, color=colours[row])
		    bar_ins_ls.append(bar_ins)
		    #yoff = yoff + matrix[row]
		    cellText.append(matrix[row])
		# Add a table at the bottom of the axes
		colours.reverse()
		cellText.reverse()
		row_label_ls_reverse = row_label_ls[:]
		row_label_ls_reverse.reverse()
		the_table = pylab.table(cellText=cellText, rowLabels=row_label_ls_reverse, rowColours=colours, colLabels=col_label_ls)
		pylab.ylabel("Counts")
		pylab.xlabel("Geographic Distance")
		#vals = arange(0, 2500, 500)
		#pylab.yticks(vals*1000, ['%d' % val for val in vals])
		#pylab.xticks([])
		pylab.title('Counts of identity pairs by region X geo distance')
		
		pylab.legend([bar_ins[0] for bar_ins in bar_ins_ls], row_label_ls, shadow=True)
		if output_fname_prefix:
			pylab.savefig('%s_identity_bar_char_by_region_x_geo_dist.svg'%output_fname_prefix, dpi=600)
			pylab.savefig('%s_identity_bar_char_by_region_x_geo_dist.png'%output_fname_prefix, dpi=600)
		pylab.show()
		sys.stderr.write("Done.\n")
	
	"""
	from misc import *
	
	input_fname = 'script/variation/stock20071008/data_d110_c0_5.tsv'
	identity_pair_ls = construct_identity_pair_ls(input_fname)
	#store it in a file as it takes a long time to generate
	import cPickle
	of = open('%s_identity_pair_ls'%os.path.splitext(input_fname)[0], 'w')
	cPickle.dump(identity_pair_ls, of)
	del of
	
	import cPickle
	inf = open('%s_identity_pair_ls'%os.path.splitext(input_fname)[0], 'r')
	identity_pair_ls = cPickle.load(inf)
	del inf
	
	#codes to group strains into populations
	popid2ecotypeid_table = 'popid2ecotypeid_25'
	pop_iden_g, popid2intra_pop_connections, popid2pos_size, ecotypeid2pop_id = construct_pop_identity_graph_from_strain_identity_ls(curs, popid2ecotypeid_table, identity_pair_ls)
	
	popid2pos = {}
	for popid, pos_size in popid2pos_size.iteritems():
		popid2pos[popid] = pos_size[0]
	
	import math
	popid2weight = {}
	for popid, intra_pop_connections in popid2intra_pop_connections.iteritems():
		popid2weight[popid] = 8*math.log(intra_pop_connections+1)
	
	draw_graph_on_map(pop_iden_g, popid2weight, popid2pos, 'inter population identity map '+popid2ecotypeid_table, output_fname_prefix='identity_map1')
	
	##graph of strains
	strain_iden_g = construct_graph_out_of_identity_pair(identity_pair_ls)
	ecotypeid2pos = get_ecotypeid2pos(curs, 'ecotype')
	ecotypeid2weight = {}
	for n in strain_iden_g:
		ecotypeid2weight[n] = 1
	
	import networkx as nx
	strain_iden_g_cc = nx.connected_components(strain_iden_g)
	
	#2007-10-09 draw all strains on map
	output_fname_prefix = '%s'%os.path.splitext(input_fname)[0]
	strain_map_fname_prefix = '%s_strain_map'%output_fname_prefix
	draw_strains_on_map(ecotypeid2pos.keys(), ecotypeid2pos, 'map of strains',  get_pic_area(ecotypeid2pos.values()), strain_map_fname_prefix)
	
	ecotypeid2pos_with_data = get_ecotypeid2pos(curs, 'ecotype', 1)
	draw_strains_on_map(ecotypeid2pos_with_data.keys(), ecotypeid2pos_with_data, 'map of strains with snp data',  get_pic_area(ecotypeid2pos_with_data.values()), '%s_with_data'%strain_map_fname_prefix)
	
	#draw cc of strain_iden_g
	for i in range(len(strain_iden_g_cc)):
		gs = strain_iden_g.subgraph(strain_iden_g_cc[i])
		site_g, site2weight, site2pos = construct_site_graph_out_of_strain_graph(gs, ecotypeid2pos)
		for site, weight in site2weight.iteritems():
			site2weight[site] =  8*math.log(weight+1)
		
		display_matrix_of_component(input_fname, strain_iden_g_cc[i], ecotypeid2pos, output_fname='ecotype_identity_map_cc%s'%i, need_sort=0, need_savefig=1)
		draw_graph_on_map(site_g, site2weight, site2pos, 'ecotype identity map cc%s'%i, pic_area=[-130,34,35,65], output_fname_prefix='ecotype_identity_cc%s'%i)
	
	
	#2007-10-02 use clique to represent unique haplotype
	strain_iden_g_cc, clique_ls_indexed_by_cc_number = get_clique_ls_indexed_by_cc_number(strain_iden_g)
	haplotype_size2number = get_haplotype_size2number(clique_ls_indexed_by_cc_number)
	singleton_strain_id_ls = get_singleton_strain_id_ls(strain_iden_g, input_fname)
	haplotype_size2number[1] += len(singleton_strain_id_ls)	#need to make up the singletons lost in graph
	hist_output_fname_prefix = '%s_haplotype_size_bar'%os.path.splitext(input_fname)[0]
	draw_histogram_of_haplotype_size(haplotype_size2number, hist_output_fname_prefix)
	ordered_clique_ls = get_ordered_clique_ls(clique_ls_indexed_by_cc_number)
	for i in range(len(ordered_clique_ls)):
		gs = strain_iden_g.subgraph(ordered_clique_ls[i])
		site_g, site2weight, site2pos = construct_site_graph_out_of_strain_graph(gs, ecotypeid2pos, node_weight_by_no_of_nodes=1)
		display_matrix_of_component(input_fname, ordered_clique_ls[i], ecotypeid2pos, output_fname='haplotype_%s_size%s'%(i,len(ordered_clique_ls[i])), need_sort=0, need_savefig=1)
		draw_graph_on_map(site_g, site2weight, site2pos, 'haplotype %s size %s'%(i,len(ordered_clique_ls[i])), pic_area=get_pic_area(site2pos.values()), output_fname_prefix='haplotype_%s_size%s_map'%(i,len(ordered_clique_ls[i])))
	
	europe_lon_span = [-12,50]
	norame_lon_span = [-130,-60]
	centralasia_lon_span = [60,90]
	japan_lon_span = [130, 150]
	norame_strain_ls = get_strains_within_longitude_span(ordered_clique_ls[0], ecotypeid2pos, norame_lon_span)
	japan_strain_ls = get_strains_within_longitude_span(ordered_clique_ls[0], ecotypeid2pos, japan_lon_span)
	
	
	cross_atlantic_cc2lon_pair_set = check_cross_ocean_components(strain_iden_g, ordered_clique_ls, ecotypeid2pos, longi_watershed=-30)
	cross_eurasia_cc2lon_pair_set = check_cross_ocean_components(strain_iden_g, ordered_clique_ls, ecotypeid2pos, longi_watershed=70)
	
	
	#2007-10-03 to see whether geographic distance correlates with genotype distance
	import numpy
	strain_pair_ls, geno_dist_ls, geo_dist_ls =  calGenoDistAndGeoDist(input_fname, ecotypeid2pos)
	
	#store them into cPickle
	import cPickle
	of = open('%s_strain_pair_geno_dist_geo_dist_ls'%os.path.splitext(input_fname)[0], 'w')
	cPickle.dump([strain_pair_ls, geno_dist_ls, geo_dist_ls], of)
	del of
	
	geo_dist_ls = numpy.array(geo_dist_ls)
	geno_dist_ls = numpy.array(geno_dist_ls)
	
	output_fname_prefix = '%s'%os.path.splitext(input_fname)[0]
	geo_output_fname_prefix = '%s_geo_distance_hist'%output_fname_prefix
	plotHist(geo_dist_ls, title="Histogram of non-NA geographic distances", output_fname_prefix=geo_output_fname_prefix)
	
	import random
	index_ls = range(len(geo_dist_ls))
	index_selected_ls = random.sample(index_ls, 10000)
	geo_dist_selected_ls = geo_dist_ls[index_selected_ls]
	geno_dist_selected_ls = geno_dist_ls[index_selected_ls]
	
	geno_vs_geo_output_fname_prefix = '%s_geno_vs_geo_dist'%output_fname_prefix
	plotXY(geo_dist_selected_ls, geno_dist_selected_ls, title='genotype vs geographic distance', xlabel="geographic dist", ylabel='genotype dist', output_fname_prefix=geno_vs_geo_output_fname_prefix)
	plot_LD(geo_dist_selected_ls, geno_dist_selected_ls,  title='genotype vs geographic distance', xlabel="geographic dist", ylabel='genotype dist', output_fname_prefix=geno_vs_geo_output_fname_prefix)
	
	#2007-10-04 divide data into continents to see whether geographic distance correlates with genotype distance
	europe_lon_span = [-12,50]
	norame_lon_span = [-130,-60]
	centralasia_lon_span = [60,90]
	japan_lon_span = [130, 150]
	def sample_geno_geo_correlation(geno_dist_ls, geo_dist_ls, output_fname_prefix):
		import numpy
		geo_dist_ls = numpy.array(geo_dist_ls)
		geno_dist_ls = numpy.array(geno_dist_ls)
		import random
		index_ls = range(len(geo_dist_ls))
		index_selected_ls = random.sample(index_ls, 10000)
		geo_dist_selected_ls = geo_dist_ls[index_selected_ls]
		geno_dist_selected_ls = geno_dist_ls[index_selected_ls]
		plot_LD(geo_dist_selected_ls, geno_dist_selected_ls,  title='genotype vs geographic distance', xlabel="geographic dist", ylabel='genotype dist', output_fname_prefix=output_fname_prefix)
	
	eur_strain_pair_ls, eur_geno_dist_ls, eur_geo_dist_ls =  calGenoDistAndGeoDist(input_fname, ecotypeid2pos, europe_lon_span)
	eur_geno_vs_geo_output_fname_prefix = '%s_geno_vs_geo_dist_eur'%output_fname_prefix
	sample_geno_geo_correlation(eur_geno_dist_ls, eur_geo_dist_ls, eur_geno_vs_geo_output_fname_prefix)
	
	noramer_strain_pair_ls, noramer_geno_dist_ls, noramer_geo_dist_ls =  calGenoDistAndGeoDist(input_fname, ecotypeid2pos, norame_lon_span)
	noramer_geno_vs_geo_output_fname_prefix = '%s_geno_vs_geo_dist_noramer'%output_fname_prefix
	sample_geno_geo_correlation(noramer_geno_dist_ls, noramer_geo_dist_ls, noramer_geno_vs_geo_output_fname_prefix)
	
	eur_noramer_strain_pair_ls, eur_noramer_geno_dist_ls, eur_noramer_geo_dist_ls =  calGenoDistAndGeoDistBetweenTwoAreas(input_fname, ecotypeid2pos, europe_lon_span, norame_lon_span)
	eur_noramer_geno_vs_geo_output_fname_prefix = '%s_geno_vs_geo_dist_eur_noramer'%output_fname_prefix
	sample_geno_geo_correlation(eur_noramer_geno_dist_ls, eur_noramer_geo_dist_ls, eur_noramer_geno_vs_geo_output_fname_prefix)
	
	#2007-10-01 parition cc into cliques
	from PartitionGraphIntoCliques import PartitionGraphIntoCliques
	PartitionGraphIntoCliques_ins = PartitionGraphIntoCliques(0)
	for for i in range(len(strain_iden_g_cc)):
		g_cc = strain_iden_g.subgraph(strain_iden_g_cc[i])
		PartitionGraphIntoCliques_ins.partition(g_cc.copy())
		map(len, PartitionGraphIntoCliques_ins.clique_ls)
		for j in range(len(PartitionGraphIntoCliques_ins.clique_ls)):
			gs = g_cc.subgraph(PartitionGraphIntoCliques_ins.clique_ls[j])
			site_g, site2weight, site2pos = construct_site_graph_out_of_strain_graph(gs, ecotypeid2pos)
			for site, weight in site2weight.iteritems():
				site2weight[site] =  8*math.log(weight+1)
			
			display_matrix_of_component(input_fname, PartitionGraphIntoCliques_ins.clique_ls[j], ecotypeid2pos, output_fname='ecotype_identity_map_cc%sc%s'%(i,j), need_sort=0, need_savefig=1)
			draw_graph_on_map(site_g, site2weight, site2pos, 'ecotype identity map cc%s clique%s'%(i,j), pic_area=[-130,34,120,65], output_fname_prefix='ecotype_identity_cc%sc%s'%(i,j))
	
	#2008-01-25 to investigate how identity spreads against distance in different regions.
	europe_lon_span = [-12,50]
	norame_lon_span = [-130,-60]
	centralasia_lon_span = [60,90]
	japan_lon_span = [130, 150]
	
	span_ls = [europe_lon_span, norame_lon_span]
	span_label_ls = ['europe', 'nor america']
	
	dist_range = [0.0, 0.4]
	strain_pair_ls, geno_dist_ls, geo_dist_ls =  calGenoDistAndGeoDist(input_fname, ecotypeid2pos, dist_range=dist_range)
	
	def float2str(input):
		input = repr(input)
		input = input.replace('.', '_')
		return input
	
	import cPickle
	output_fname = '%s_strain_pair_geno_dist_%s_%s_geo_dist_ls'%(os.path.splitext(input_fname)[0], dist_range[0], dist_range[1])
	output_fname = output_fname.replace('.', '_')
	of = open(output_fname, 'w')
	cPickle.dump([strain_pair_ls, geno_dist_ls, geo_dist_ls], of)
	del of
	
	matrix_of_counts_by_region_x_geo_dist, row_label_ls, col_label_ls, col_index_ls = group_identity_pair_according_to_region_distance(strain_pair_ls, geno_dist_ls, geo_dist_ls, ecotypeid2pos, span_ls, span_label_ls, geo_distance_window_size=100, max_identity_geno_dist=0.0)
	draw_table_as_bar_chart(matrix_of_counts_by_region_x_geo_dist, row_label_ls, col_label_ls, col_index_ls)
	"""
	
	"""
	2007-09-24
		magnus's idea to reduce fake identities
	"""
	def get_data_matrix_sorted_by_NA_perc(input_fname):
		"""
		2007-09-24
			input_fname is a data frame file with snps coded in integers. 1st row is the ecotype id in table ecotype
		"""
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
		no_of_strains = len(strain_acc_list)
		no_of_snps = len(header)-2
		data_matrix_sorted_by_NA_perc = []
		for i in range(no_of_strains):
			no_of_NAs = 0
			for j in range(no_of_snps):
				if data_matrix[i][j] == 0:
					no_of_NAs += 1
			data_matrix_sorted_by_NA_perc.append([float(no_of_NAs)/no_of_snps, data_matrix[i], strain_acc_list[i]])
		data_matrix_sorted_by_NA_perc.sort()
		return data_matrix_sorted_by_NA_perc
	
	def get_strain_id2index_from_data_matrix_sorted_by_NA_perc(data_matrix_sorted_by_NA_perc):
		sys.stderr.write("Getting strain_id2index from data_matrix_sorted_by_NA_perc ...")
		strain_id2index = {}
		no_of_strains = len(data_matrix_sorted_by_NA_perc)
		no_of_snps = len(data_matrix_sorted_by_NA_perc[0][1])
		for i in range(no_of_strains):
			ecotypeid = data_matrix_sorted_by_NA_perc[i][2]
			strain_id2index[ecotypeid] = i
		sys.stderr.write("Done.\n")
		return strain_id2index
	
	
	def get_unique_haplotype_set(data_matrix_sorted_by_NA_perc):
		sys.stderr.write("Getting unique haplotype set ...")
		no_of_strains = len(data_matrix_sorted_by_NA_perc)
		no_of_snps = len(data_matrix_sorted_by_NA_perc[0][1])
		from sets import Set
		unique_haplotype_set = Set()
		for i in range(no_of_strains):
			no_of_prev_identity_strains = 0
			for j in range(i):
				no_of_same_cols = 0
				for k in range(no_of_snps):
					call1 = data_matrix_sorted_by_NA_perc[i][1][k]
					call2 = data_matrix_sorted_by_NA_perc[j][1][k]
					if call1 == call2 or call1==0 or call2==0:
						no_of_same_cols += 1
				if no_of_same_cols == no_of_snps:
					no_of_prev_identity_strains += 1
			if no_of_prev_identity_strains==0:
				unique_haplotype_set.add(data_matrix_sorted_by_NA_perc[i][2])
		sys.stderr.write("Done.\n")
		return unique_haplotype_set
	
	
	def assign_remainder_strains_to_unique_haplotype_set(unique_haplotype_set, data_matrix_sorted_by_NA_perc, strain_id2index):
		sys.stderr.write("Assigning remainder strains to unique_haplotype_set ...")
		identity_pair_ls = []
		no_of_strains = len(data_matrix_sorted_by_NA_perc)
		no_of_snps = len(data_matrix_sorted_by_NA_perc[0][1])
		for i in range(no_of_strains):
			ecotypeid = data_matrix_sorted_by_NA_perc[i][2]
			if ecotypeid not in unique_haplotype_set:	#skip strains in unique_haplotype_set
				prev_identity_strains = []
				for seedid in unique_haplotype_set:
					seed_data_row = data_matrix_sorted_by_NA_perc[strain_id2index[seedid]][1]
					no_of_same_cols = 0
					for k in range(no_of_snps):
						call1 = data_matrix_sorted_by_NA_perc[i][1][k]
						call2 = seed_data_row[k]
						if call1 == call2 or call1==0 or call2==0:
							no_of_same_cols += 1
					if no_of_same_cols == no_of_snps:
						prev_identity_strains.append(seedid)
				if len(prev_identity_strains)==1:	#uniquely assigned
					identity_pair_ls.append([ecotypeid, prev_identity_strains[0]])
		sys.stderr.write("Done.\n")
		return identity_pair_ls
	
	def construct_identity_pair_ls_given_strain_id_list(strain_id_list, data_matrix_sorted_by_NA_perc, strain_id2index):
		sys.stderr.write("Constructing identity_pair_ls given strain_id_list ...")
		identity_pair_ls = []
		no_of_strains = len(strain_id_list)
		no_of_snps = len(data_matrix_sorted_by_NA_perc[0][1])
		for i in range(no_of_strains):
			ecotypeid1 = repr(strain_id_list[i])	#convert int to char
			index_of_ecotypeid1 = strain_id2index[ecotypeid1]
			for j in range(i+1, no_of_strains):
				ecotypeid2 = repr(strain_id_list[j])
				index_of_ecotypeid2 = strain_id2index[ecotypeid2]
				no_of_same_cols = 0
				for k in range(no_of_snps):
					call1 = data_matrix_sorted_by_NA_perc[index_of_ecotypeid1][1][k]
					call2 = data_matrix_sorted_by_NA_perc[index_of_ecotypeid2][1][k]
					if call1 == call2 or call1==0 or call2==0:
						no_of_same_cols += 1
				if no_of_same_cols == no_of_snps:
					identity_pair_ls.append([ecotypeid1, ecotypeid2])
		sys.stderr.write("Done.\n")
		return identity_pair_ls
	
	
	"""
	input_fname = 'script/variation/stock20070919/data_d110_c0_5.tsv'
	data_matrix_sorted_by_NA_perc = get_data_matrix_sorted_by_NA_perc(input_fname)
	
	unique_haplotype_set = get_unique_haplotype_set(data_matrix_sorted_by_NA_perc)
	
	strain_id2index = get_strain_id2index_from_data_matrix_sorted_by_NA_perc(data_matrix_sorted_by_NA_perc)
	
	identity_pair_ls2 = assign_remainder_strains_to_unique_haplotype_set(unique_haplotype_set, data_matrix_sorted_by_NA_perc, strain_id2index)
	
	import cPickle
	of = open('%s_identity_pair_ls2'%os.path.splitext(input_fname)[0], 'w')
	cPickle.dump(identity_pair_ls2, of)
	del of
	popid2ecotypeid_table = 'popid2ecotypeid_25'
	pop_iden_g2, popid2intra_pop_connections, popid2pos_size, ecotypeid2pop_id = construct_pop_identity_graph_from_strain_identity_ls(curs, popid2ecotypeid_table, identity_pair_ls2)
	
	popid2pos = {}
	for popid, pos_size in popid2pos_size.iteritems():
		popid2pos[popid] = pos_size[0]
	
	import math
	popid2weight = {}
	for popid, intra_pop_connections in popid2intra_pop_connections.iteritems():
		popid2weight[popid] = 8*math.log(intra_pop_connections+1)
	
	draw_graph_on_map(pop_iden_g2, popid2weight, popid2pos, 'inter population identity map '+popid2ecotypeid_table, output_fname_prefix='identity_map2')
	
	
	strain_iden_g2 = construct_graph_out_of_identity_pair(identity_pair_ls2)
	
	ecotypeid2pos = get_ecotypeid2pos(curs, 'ecotype')
	ecotypeid2weight = {}
	for n in strain_iden_g:
		ecotypeid2weight[n] = 1
	
	import networkx as nx
	strain_iden_g_cc2 = nx.connected_components(strain_iden_g2)
	
	i=11
	
	
	for i in range(len(strain_iden_g_cc2)):
		gs = strain_iden_g2.subgraph(strain_iden_g_cc2[i])
		site_g, site2weight, site2pos = construct_site_graph_out_of_strain_graph(gs, ecotypeid2pos)
		for site, weight in site2weight.iteritems():
			site2weight[site] =  8*math.log(weight+1)
		
		display_matrix_of_component(input_fname, strain_iden_g_cc2[i], ecotypeid2pos, output_fname='ecotype_identity_tcc%s'%i, need_sort=0, need_savefig=1)
		draw_graph_on_map(site_g, site2weight, site2pos, 'ecotype identity map cc%s'%i, output_fname_prefix='ecotype_identity_tcc%s'%i)
	
	
	#2007-09-25 complete the identity graph
	identity_pair_ls2_complete = construct_identity_pair_ls_given_strain_id_list(strain_iden_g2.nodes(), data_matrix_sorted_by_NA_perc, strain_id2index)
	
	pop_iden_g2_complete, popid2intra_pop_connections, popid2pos_size, ecotypeid2pop_id = construct_pop_identity_graph_from_strain_identity_ls(curs, popid2ecotypeid_table, identity_pair_ls2_complete)
	
	popid2pos = {}
	for popid, pos_size in popid2pos_size.iteritems():
		popid2pos[popid] = pos_size[0]
	
	import math
	popid2weight = {}
	for popid, intra_pop_connections in popid2intra_pop_connections.iteritems():
		popid2weight[popid] = 8*math.log(intra_pop_connections+1)
	
	draw_graph_on_map(pop_iden_g2_complete, popid2weight, popid2pos, 'inter population identity map '+popid2ecotypeid_table, output_fname_prefix='identity_map2_complete')
	
	strain_iden_g2_complete = construct_graph_out_of_identity_pair(identity_pair_ls2_complete)
	strain_iden_g_cc2_complete = nx.connected_components(strain_iden_g2_complete)
	
	for i in range(len(strain_iden_g_cc2_complete)):
		gs = strain_iden_g2_complete.subgraph(strain_iden_g_cc2_complete[i])
		site_g, site2weight, site2pos = construct_site_graph_out_of_strain_graph(gs, ecotypeid2pos)
		for site, weight in site2weight.iteritems():
			site2weight[site] =  8*math.log(weight+1)
		
		display_matrix_of_component(input_fname, strain_iden_g_cc2_complete[i], ecotypeid2pos, output_fname='ecotype_identity_cc2_complete%s'%i, need_sort=0, need_savefig=1)
		draw_graph_on_map(site_g, site2weight, site2pos, 'ecotype identity map cc%s'%i, output_fname_prefix='ecotype_identity_cc2_complete%s'%i)
	
	
	cross_ocean_cc_set2_atlantic = check_cross_ocean_components(strain_iden_g2_complete, strain_iden_g_cc2_complete, ecotypeid2pos)
	
	cross_ocean_cc_set2_eurasia = check_cross_ocean_components(strain_iden_g2_complete, strain_iden_g_cc2_complete, ecotypeid2pos, 60)
	"""



class LD:
	"""
	2007-09-26
		calculate LD (r^2)
	2008-01-15
		put everything into class
	"""
	def __init__(self):
		pass
	
	def group_ordered_snps_into_chr_snp_2layer_ls(self, curs, snp_acc_list, snp_locus_table='snps'):
		"""
		2007-09-26
			assume snps are already in chromosome, position order, just need to find out
			where to stop the chromosome
		"""
		sys.stderr.write("Grouping ordered snps into chr_snp_2layer_ls ...")
		old_chromosome = -1
		chr_snp_2layer_ls = []
		snp_position_ls = []	#no chromosome, just position.
		chr_ls = []
		for i in range(len(snp_acc_list)):
			curs.execute("select chromosome, position from %s where snpid='%s'"%(snp_locus_table, snp_acc_list[i]))
			rows = curs.fetchall()
			chromosome, position = rows[0]
			snp_position_ls.append(position)
			if old_chromosome == -1:	#1st encounter
				old_chromosome = chromosome
			elif chromosome!=old_chromosome:
				chr_snp_2layer_ls.append(chr_ls)
				chr_ls = []
				old_chromosome = chromosome
			chr_ls.append(i)
		chr_snp_2layer_ls.append(chr_ls)
		sys.stderr.write("Done.\n")
		return chr_snp_2layer_ls, snp_position_ls
	
	def fill_in_snp_allele2index(self, diploid_allele, allele2index):
		from common import number2nt, nt2number
		if diploid_allele>4:
			nt = number2nt[diploid_allele]
			allele1 = nt2number[nt[0]]
			allele2 = nt2number[nt[1]]
		else:
			allele1 = allele2 = diploid_allele
		if allele1 not in allele2index:
			allele2index[allele1] = len(allele2index)
		if allele2 not in allele2index:
			allele2index[allele2] = len(allele2index)
		return allele1, allele2
	
	def calculate_LD(self, input_fname, curs, snp_locus_table='snps', debug=0, check_bit_ls=[1,0,0,0,0,0,0,0,0,0], chr_ls=[]):
		"""
		exclude pairs with one or two NAs
		exclude pairs both of who are heterozygous calls (can't figure out the phase)
		(only one of pairs is heterozygous is all right)
		2008-01-23 add chr_ls
		"""
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
		if debug:
			import pdb
			pdb.set_trace()
		import random
		no_of_strains = len(strain_acc_list)
		snp_acc_list = header[2:]
		chr_snp_2layer_ls, snp_position_ls = self.group_ordered_snps_into_chr_snp_2layer_ls(curs, snp_acc_list, snp_locus_table)
		no_of_chrs = len(chr_snp_2layer_ls)
		if len(chr_ls)==0:	#2008-01-23 if there's nothing in chr_ls
			chr_ls = range(no_of_chrs)
		r_square_ls = []
		distance_ls = []
		allele_freq_ls = []
		snp_pair_ls = []
		D_ls = []
		D_prime_ls = []
		import Numeric
		counter = 0
		period = len(check_bit_ls)
		for c in chr_ls:
			no_of_snps = len(chr_snp_2layer_ls[c])
			for i in range(no_of_snps):
				if len(check_bit_ls)>1:	#2008-01-23 no shuffling if there's only 1 bit
					random.shuffle(check_bit_ls)
				for j in range(i+1, no_of_snps):
					snp1_index = chr_snp_2layer_ls[c][i]
					snp2_index = chr_snp_2layer_ls[c][j]
					distance = abs(snp_position_ls[snp1_index]-snp_position_ls[snp2_index])
					if distance>10000:	#2008-01-15 skip SNPs two far
						continue
					counter += 1
					if check_bit_ls[counter%period]==0:	#2008-01-15 skip some portion of them randomly
						continue
					counter_matrix = Numeric.zeros([2,2])
					snp1_allele2index = {}
					snp2_allele2index = {}
					for k in range(no_of_strains):
						snp1_allele = data_matrix[k][snp1_index]
						snp2_allele = data_matrix[k][snp2_index]
						if snp1_allele!=0 and snp2_allele!=0 and not (snp1_allele>4 and snp2_allele>4):
							snp1_allele1, snp1_allele2 = self.fill_in_snp_allele2index(snp1_allele, snp1_allele2index)
							snp2_allele1, snp2_allele2 = self.fill_in_snp_allele2index(snp2_allele, snp2_allele2index)
							counter_matrix[snp1_allele2index[snp1_allele1],snp2_allele2index[snp2_allele1]] += 1
							counter_matrix[snp1_allele2index[snp1_allele2],snp2_allele2index[snp2_allele2]] += 1
					PA = sum(counter_matrix[0,:])
					Pa = sum(counter_matrix[1,:])
					PB = sum(counter_matrix[:,0])
					Pb = sum(counter_matrix[:,1])
					total_num = float(PA+Pa)
					try:
						PA = PA/total_num
						Pa = Pa/total_num
						PB = PB/total_num
						Pb = Pb/total_num
						PAB = counter_matrix[0,0]/total_num
						D = PAB-PA*PB
						PAPB = PA*PB
						PAPb = PA*Pb
						PaPB = Pa*PB
						PaPb = Pa*Pb
						Dmin = max(-PAPB, -PaPb)
						Dmax = min(PAPb, PaPB)
						if D<0:
							D_prime = D/Dmin
						else:
							D_prime = D/Dmax
						r2 = D*D/(PA*Pa*PB*Pb)
					except:	#2008-01-23 exceptions.ZeroDivisionError, Dmin or Dmax could be 0 if one of(-PAPB, -PaPb)  is >0 or <0
						sys.stderr.write('Unknown except, ignore: %s\n'%repr(sys.exc_info()[0]))
						continue
					allele_freq = (min(PA, Pa),min(PB, Pb))
					D_ls.append(D)
					D_prime_ls.append(D_prime)
					r_square_ls.append(r2)
					distance_ls.append(distance)
					allele_freq_ls.append(allele_freq)
					snp_pair_ls.append((snp_acc_list[i], snp_acc_list[j]))
		return D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls
	
	
	def plot_LD(self, x_ls, y_ls, title, xlabel, ylabel, max_dist=0, output_fname_prefix=None):
		import pylab
		import numpy
		import scipy.interpolate
		#x_ls = numpy.array(x_ls)
		#y_ls = numpy.array(y_ls)
		x_argsort_ls = numpy.argsort(x_ls)
		new_x_ls = []
		new_y_ls = []
		for i in range(len(x_argsort_ls)):
			if max_dist>0:
				if x_ls[x_argsort_ls[i]]<=max_dist:
					new_x_ls.append(x_ls[x_argsort_ls[i]])
					new_y_ls.append(y_ls[x_argsort_ls[i]])
			else:
				new_x_ls.append(x_ls[x_argsort_ls[i]])
				new_y_ls.append(y_ls[x_argsort_ls[i]])
		
		sp = scipy.interpolate.UnivariateSpline(new_x_ls,new_y_ls)
		step = (new_x_ls[-1]-new_x_ls[0])/100
		n_x_ls = numpy.arange(new_x_ls[0], new_x_ls[-1], step)
		n_y_ls = map(sp, n_x_ls)
		pylab.clf()
		pylab.title(title)
		pylab.xlabel(xlabel)
		pylab.ylabel(ylabel)
		pylab.plot(new_x_ls, new_y_ls, '.')
		pylab.plot(n_x_ls, n_y_ls)
		if output_fname_prefix:
			if max_dist:
				output_fname_prefix = '%s_%s'%(output_fname_prefix, max_dist)
			pylab.savefig('%s.eps'%output_fname_prefix, dpi=450)
			pylab.savefig('%s.svg'%output_fname_prefix, dpi=450)
			pylab.savefig('%s.png'%output_fname_prefix, dpi=450)
	
	def get_snp_acc2LD_data(self, LD_ls, distance_ls, snp_pair_ls):
		"""
		2008-01-24
			the previus model doesn't give good results.
			curve and real data are quite far apart. it's easily disrupted by points with low LD. axis distance=0 (the separating line) sometimes sits in the middle, like 5kb.
			so try Alex Platt's suggestion:
				r^2 = ax^b => log(r^2) = log(a) + b*log(x)

		2008-01-23
			process data for linear model fitting
			r^2 = 1/(a+bx) (page 483 of Hartl2007) => 1/r^2 = a+bx
		"""
		sys.stderr.write("Getting snp_acc2LD_data ...")
		snp_acc2LD_data = {}
		import math
		for i in range(len(snp_pair_ls)):
			snp_pair = snp_pair_ls[i]
			for snp_acc in snp_pair:
				if snp_acc not in snp_acc2LD_data:
					snp_acc2LD_data[snp_acc] = [[], []]
				if LD_ls[i]>0 and distance_ls[i]>0:
					snp_acc2LD_data[snp_acc][0].append(math.log(LD_ls[i]))
					snp_acc2LD_data[snp_acc][1].append([1, math.log(distance_ls[i])])
		sys.stderr.write("Done.\n")
		return snp_acc2LD_data
	
	def LD_linear_fitting(self, snp_acc2LD_data):
		"""
		2008-01-23
			actual linear fitting
		"""
		sys.stderr.write("Linear fitting for LD ...")
		from annot.bin.module_cc.linear_model import linear_model
		linear_model_ins = linear_model()
		snp_acc2LD_linear_fitting = {}
		for snp_acc, LD_data in snp_acc2LD_data.iteritems():
			LD_ls, X_ls = LD_data
			if len(LD_ls)>5:
				linear_model_ins.prepare_data(LD_ls, X_ls)
				linear_model_ins.run()
				coeff_list = linear_model_ins.coefficients()
				chisq_tuple = linear_model_ins.chisq_return()
				snp_acc2LD_linear_fitting[snp_acc] = [coeff_list, chisq_tuple]
				linear_model_ins.cleanup()
		del linear_model_ins
		sys.stderr.write("Done.\n")
		return snp_acc2LD_linear_fitting
	
	def check_one_snp_acc_LD_decay(self, snp_acc2LD_data, snp_acc2LD_linear_fitting, snp_acc):
		"""
		2008-01-24
			follow the model change in get_snp_acc2LD_data()
		2008-01-23
			unfold the data processed in snp_acc2LD_data
			to see whether the theoretical curve (blue line) fits real data (red dots).
		"""
		import math
		a = snp_acc2LD_linear_fitting[snp_acc][0][0]
		b = snp_acc2LD_linear_fitting[snp_acc][0][1]
		theoretical_curve_func = lambda x: math.exp(a)*math.pow(x,b)
		x_ls = range(500,10000,100)
		y_ls = map(theoretical_curve_func, x_ls)
		import pylab
		pylab.clf()
		pylab.plot(x_ls, y_ls)
		
		inverse_func = lambda x: math.exp(x)
		r2_ls = map(inverse_func, snp_acc2LD_data[snp_acc][0])
		distance_ls = [math.exp(row[1]) for row in snp_acc2LD_data[snp_acc][1]]
		pylab.plot(distance_ls, r2_ls, 'r.')
		pylab.show()
	
	"""
from variation.src.misc import LD
input_fname = os.path.expanduser('~/script/variation/stock20070919/data_d110_c0_5.tsv')
snp_locus_table = 'snps'
LD_ins = LD()
D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls = LD_ins.calculate_LD(input_fname, curs, snp_locus_table)

title = 'LD decay'
xlabel = 'Distance'
ylabel = r'$r^2$'
output_fname_prefix = '%s_LD_r2'%os.path.splitext(input_fname)[0]
LD_ins.plot_LD(distance_ls, r_square_ls, title, xlabel, ylabel, 0, output_fname_prefix)
LD_ins.plot_LD(distance_ls, r_square_ls, title, xlabel, ylabel, 500000, output_fname_prefix)
LD_ins.plot_LD(distance_ls, r_square_ls, title, xlabel, ylabel, 1000000, output_fname_prefix)
LD_ins.plot_LD(distance_ls, r_square_ls, title, xlabel, ylabel, 5000000, output_fname_prefix)

ylabel = 'D'
output_fname_prefix = '%s_LD_D'%os.path.splitext(input_fname)[0]
LD_ins.plot_LD(distance_ls, D_ls, title, xlabel, ylabel, 0, output_fname_prefix)
LD_ins.plot_LD(distance_ls, D_ls, title, xlabel, ylabel, 500000, output_fname_prefix)
LD_ins.plot_LD(distance_ls, D_ls, title, xlabel, ylabel, 1000000, output_fname_prefix)
LD_ins.plot_LD(distance_ls, D_ls, title, xlabel, ylabel, 5000000, output_fname_prefix)

ylabel = "D'"
output_fname_prefix = '%s_LD_D_prime'%os.path.splitext(input_fname)[0]
LD_ins.plot_LD(distance_ls, D_prime_ls, title, xlabel, ylabel, 0, output_fname_prefix)
LD_ins.plot_LD(distance_ls, D_prime_ls, title, xlabel, ylabel, 500000, output_fname_prefix)
LD_ins.plot_LD(distance_ls, D_prime_ls, title, xlabel, ylabel, 1000000, output_fname_prefix)
LD_ins.plot_LD(distance_ls, D_prime_ls, title, xlabel, ylabel, 5000000, output_fname_prefix)


#2008-01-15 250k data
input_fname = os.path.expanduser('~/script/variation/genotyping/250ksnp/data/data_250k.tsv')
snp_locus_table = 'snps_250k'
LD_ins = LD()
#2008-01-23 chromosome 1 and no random skipping
D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls = LD_ins.calculate_LD(input_fname, curs, snp_locus_table, check_bit_ls=[1], chr_ls=[1])
D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls = LD_ins.calculate_LD(input_fname, curs, snp_locus_table, check_bit_ls=[1]+range(100), chr_ls=[2])

pickle_fname = os.path.expanduser('~/.pickle/LD_ins')
of = open(pickle_fname, 'w')
cPickle.dump([D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls], of)
del of

#2008-01-23 load the computed data by cPickle
import cPickle
inf = open(pickle_fname, 'r')
D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls = cPickle.load(inf)
del inf

snp_acc2LD_data = LD_ins.get_snp_acc2LD_data(r_square_ls, distance_ls, snp_pair_ls)
snp_acc2LD_linear_fitting = LD_ins.LD_linear_fitting(snp_acc2LD_data)

def check_one_snp_acc_LD_decay(snp_acc2LD_data, snp_acc2LD_linear_fitting, snp_acc):
	import math
	a = snp_acc2LD_linear_fitting[snp_acc][0][0]
	b = snp_acc2LD_linear_fitting[snp_acc][0][1]
	theoretical_curve_func = lambda x: math.exp(a)*math.pow(x,b)
	x_ls = range(500,10000,100)
	y_ls = map(theoretical_curve_func, x_ls)
	import pylab
	pylab.clf()
	pylab.plot(x_ls, y_ls)
	
	inverse_func = lambda x: math.exp(x)
	r2_ls = map(inverse_func, snp_acc2LD_data[snp_acc][0])
	distance_ls = [math.exp(row[1]) for row in snp_acc2LD_data[snp_acc][1]]
	pylab.plot(distance_ls, r2_ls, 'r.')
	pylab.show()


LD_ins.check_one_snp_acc_LD_decay(snp_acc2LD_data, snp_acc2LD_linear_fitting,'1_10219_T_A')
LD_ins.check_one_snp_acc_LD_decay(snp_acc2LD_data, snp_acc2LD_linear_fitting,'1_3102_A_G')

for snp_acc in snp_acc2LD_linear_fitting:
	LD_ins.check_one_snp_acc_LD_decay(snp_acc2LD_data, snp_acc2LD_linear_fitting, snp_acc)
	print snp_acc2LD_linear_fitting[snp_acc]
	raw_input(":")
	"""
	
	def plotLDHist(input_fname, output_fname_prefix, which_LD_statistic=1, min_r2=0.95, discard_perc=0.999, no_of_bins=50, \
				   pairType=0, normed=False, debug=False):
		"""
		2009-4-5
			add argument discard_perc which controls the percentage of pairs discarded randomly to reduce computational load
			add argument no_of_bins
			pairType=0, all
			pairType=1, intra-chromosomal
			pairType=2, inter-chromosomal
		02/20/09
			histogram of LD under cutoff by min_r2
		"""
		import csv, random
		from pymodule import getColName2IndexFromHeader, PassingData
		from variation.src.DrawSNPRegion import DrawSNPRegion, LD_statistic
		reader = csv.reader(open(input_fname), delimiter='\t')
		col_name2index = getColName2IndexFromHeader(reader.next())
		xylim_data = PassingData(xlim = [-1,-1], ylim=[-1,-1])
		LD_ls = []
		counter = 0
		real_counter = 0
		for row in reader:
			counter += 1
			snp1 = row[col_name2index['snp1']].split('_')
			snp1 = tuple(map(int, snp1))
			snp2 = row[col_name2index['snp2']].split('_')
			snp2 = tuple(map(int, snp2))
			if pairType==1 and snp1[0]!=snp2[0]:	#intra-chromosomal, but the two snps are from different chromosomes
				continue
			if pairType==2 and snp1[0]==snp2[0]:	#inter-chromosomal, but the two snps are from the same chromosome
				continue
			if discard_perc==0:
				u = 1
			else:	#sample a uniform unless discard_perc is not 0
				u = random.random()
			if u<discard_perc:
				continue
			
			allele1_freq = float(row[col_name2index['allele1_freq']])
			allele2_freq = float(row[col_name2index['allele2_freq']])
			#if allele1_freq>=min_MAF and allele2_freq>=min_MAF:	#meet the minimum minor-allele-frequency
			LD_stat = float(row[col_name2index[LD_statistic.get_name(which_LD_statistic)]])
			LD_stat = abs(LD_stat)
			if LD_stat>=min_r2:
				LD_ls.append(LD_stat)
				real_counter+=1
			if debug and counter==50000:
				break
			if real_counter%5000==0:
				sys.stderr.write("%s%s\t\t%s"%('\x08'*40, real_counter, counter))
		del reader
		
		import pylab
		pylab.clf()
		pylab.title('Histogram of %s'%LD_statistic.get_name(which_LD_statistic))
		pylab.ylabel('frequency')
		pylab.hist(LD_ls, no_of_bins, alpha=0.5, normed=normed)
		pylab.savefig('%s_%s_min_r2_%s_discard_perc_%s_hist.png'%\
					  (output_fname_prefix, LD_statistic.get_name(which_LD_statistic), min_r2, discard_perc), dpi=300)
	
	"""
input_fname = os.path.expanduser('~/panfs/250k/dataset/call_method_29_LD_D_prime_m0.8.tsv')
output_fname_prefix = os.path.expanduser('%s'%os.path.splitext(input_fname)[0])
plotLDHist(input_fname, output_fname_prefix, min_r2=0.95)

input_fname = os.path.expanduser('~/panfs/250k/dataset/call_method_29_LD_D_prime_m0_p.997.tsv')
output_fname_prefix = os.path.expanduser('%s_intra_chr'%os.path.splitext(input_fname)[0])
plotLDHist(input_fname, output_fname_prefix, which_LD_statistic=1, min_r2=0, discard_perc=0.994, no_of_bins=200, pairType=1, normed=True)
output_fname_prefix = os.path.expanduser('%s_inter_chr'%os.path.splitext(input_fname)[0])
plotLDHist(input_fname, output_fname_prefix, which_LD_statistic=1, min_r2=0, discard_perc=0.994, no_of_bins=200, pairType=2, normed=True)

	"""
	
	def filterLD(input_fname, output_fname, which_LD_statistic=1, min_r2=0.95, debug=False):
		"""
		2009-3-6
			select entries from input_fname that are above min_r2 into output_fname
		"""
		import csv
		from pymodule import getColName2IndexFromHeader
		from variation.src.DrawSNPRegion import DrawSNPRegion, LD_statistic
		reader = csv.reader(open(input_fname), delimiter='\t')
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		
		col_name2index = getColName2IndexFromHeader(reader.next())
		counter = 0
		real_counter= 0
		LD_ls = []
		counter = 0
		real_counter = 0
		for row in reader:
			snp1 = row[col_name2index['snp1']].split('_')
			snp1 = tuple(map(int, snp1))
			snp2 = row[col_name2index['snp2']].split('_')
			snp2 = tuple(map(int, snp2))
			allele1_freq = float(row[col_name2index['allele1_freq']])
			allele2_freq = float(row[col_name2index['allele2_freq']])
			#if allele1_freq>=min_MAF and allele2_freq>=min_MAF:	#meet the minimum minor-allele-frequency
			LD_stat = float(row[col_name2index[LD_statistic.get_name(which_LD_statistic)]])
			LD_stat = abs(LD_stat)
			if LD_stat>=min_r2:
				#LD_ls.append(LD_stat)
				writer.writerow(row)
				real_counter+=1
			counter += 1
			if debug and counter==50000:
				break
			if real_counter%2000==0:
				sys.stderr.write("%s%s\t\t%s"%('\x08'*40, real_counter, counter))
		del reader, writer
		
	"""
input_fname = os.path.expanduser('~/panfs/250k/call_method_29_LD_D_prime_m0.8.tsv')
min_r2=0.85
output_fname = '%s_min_r2_%s.tsv'%(os.path.expanduser('%s'%os.path.splitext(input_fname)[0]), min_r2)
filterLD(input_fname, output_fname, min_r2=min_r2)
	
	"""
	
	def getSNPPairDistWithLDAboveMin(input_fname, which_LD_statistic=1, min_r2=0.95, debug=False):
		"""
		2009-3-6
			
		"""
		import csv
		from pymodule import getColName2IndexFromHeader
		from variation.src.DrawSNPRegion import DrawSNPRegion, LD_statistic
		reader = csv.reader(open(input_fname), delimiter='\t')
		
		col_name2index = getColName2IndexFromHeader(reader.next())
		counter = 0
		real_counter= 0
		LD_ls = []
		counter = 0
		real_counter = 0
		inter_chr_pair_ls = []
		snp_pair_dist_ls = []
		for row in reader:
			snp1 = row[col_name2index['snp1']].split('_')
			snp1 = tuple(map(int, snp1))
			snp2 = row[col_name2index['snp2']].split('_')
			snp2 = tuple(map(int, snp2))
			allele1_freq = float(row[col_name2index['allele1_freq']])
			allele2_freq = float(row[col_name2index['allele2_freq']])
			#if allele1_freq>=min_MAF and allele2_freq>=min_MAF:	#meet the minimum minor-allele-frequency
			LD_stat = float(row[col_name2index[LD_statistic.get_name(which_LD_statistic)]])
			LD_stat = abs(LD_stat)
			if LD_stat>=min_r2:
				#LD_ls.append(LD_stat)
				real_counter+=1
				if snp1[0]!=snp2[0]:
					inter_chr_pair_ls.append((snp1, snp2))
				else:
					snp_pair_dist = abs(snp1[1]-snp2[1])
					snp_pair_dist_ls.append(snp_pair_dist)
			
			counter += 1
			if debug and counter==50000:
				break
			if real_counter%2000==0:
				sys.stderr.write("%s%s\t\t%s"%('\x08'*40, real_counter, counter))
		del reader
		return inter_chr_pair_ls, snp_pair_dist_ls
	
	"""
input_fname = os.path.expanduser('~/panfs/250k/call_method_29_LD_D_prime_m0.8_min_r2_0.85.tsv')
min_r2 = 0.95
inter_chr_pair_ls, snp_pair_dist_ls = getSNPPairDistWithLDAboveMin(input_fname, min_r2=min_r2)
	"""
	class Counter(object):
		"""
		2009-4-8
			a class for calSNPPairStats()
		"""
		def __init__(self):
			self.distance_ls = []
			self.lowest_freq_ls = []
			self.no_of_intrachr_pairs = 0
			self.no_of_interchr_pairs = 0

	def calSNPPairStats(input_fname, output_fname_prefix, discard_perc=0.999, no_of_bins=50, \
					normed=False, debug=False):
		"""
		2009-4-8
			given the data which records the frequency of each allelic combination from two SNPs
				count how many SNP pairs with 4 combos, how many with 3, how many with 2.
					and among each category, how many are inter-chromosomal, how many are intra.
						for the intra pairs, draw histogram of distance
					also among each category, draw histogram of the frequency of the least-frequent combo
		"""
		import csv, random
		reader = csv.reader(open(input_fname), delimiter='\t')
		counter_dc_idx = 9
		from pymodule import getColName2IndexFromHeader, PassingData
		from variation.src.DrawSNPRegion import DrawSNPRegion, LD_statistic
		reader = csv.reader(open(input_fname), delimiter='\t')
		col_name2index = getColName2IndexFromHeader(reader.next())
		counter = 0
		real_counter= 0
		
		no_of_combos2counter = {}
		for row in reader:
			counter_dict = eval(row[counter_dc_idx])
			no_of_combos = len(counter_dict)
			if no_of_combos not in no_of_combos2counter:
				no_of_combos2counter[no_of_combos] = Counter()
			
			counter += 1
			snp1 = row[col_name2index['snp1_id']].split('_')
			snp1 = tuple(map(int, snp1))
			snp2 = row[col_name2index['snp2_id']].split('_')
			snp2 = tuple(map(int, snp2))
			
			if discard_perc==0:
				u = 1
			else:	#sample a uniform unless discard_perc is not 0
				u = random.random()
			
			if snp1[0]!=snp2[0]:	#intra-chromosomal, but the two snps are from different chromosomes
				no_of_combos2counter[no_of_combos].no_of_interchr_pairs += 1
				
			if snp1[0]==snp2[0]:	#inter-chromosomal, but the two snps are from the same chromosome
				no_of_combos2counter[no_of_combos].no_of_intrachr_pairs += 1
				if u>=discard_perc:
					snp_pair_dist = abs(snp1[1]-snp2[1])
					no_of_combos2counter[no_of_combos].distance_ls.append(snp_pair_dist)
			
			if u>=discard_perc:
				lowest_no_between_two_snps = min(counter_dict.values())
				no_of_total_accessions = sum(counter_dict.values())
				lowest_freq = lowest_no_between_two_snps/float(no_of_total_accessions)
				no_of_combos2counter[no_of_combos].lowest_freq_ls.append(lowest_freq)
				real_counter+=1
			if debug and counter==50000:
				break
			if real_counter%5000==0:
				sys.stderr.write("%s%s\t\t%s"%('\x08'*40, real_counter, counter))
		del reader
		sys.stderr.write("%s%s\t\t%s\n"%('\x08'*40, real_counter, counter))
		
		import pylab
		for no_of_combos, counter_obj in no_of_combos2counter.iteritems():
			no_of_total_pairs = counter_obj.no_of_interchr_pairs + counter_obj.no_of_intrachr_pairs
			perc = counter_obj.no_of_interchr_pairs/float(no_of_total_pairs)
			print '%s/%s=%.5f inter-chromosomal SNP pairs for pairs with %s combos'%(counter_obj.no_of_interchr_pairs, no_of_total_pairs, perc, no_of_combos)
			pylab.clf()
			pylab.title('Histogram of Intra-chromosal Dist For SNPPairs with %s combos'%no_of_combos)
			pylab.ylabel('frequency')
			pylab.hist(counter_obj.distance_ls, no_of_bins, alpha=0.5, normed=normed)
			pylab.savefig('%s_%s_combos_discard_perc_%s_dist_hist.png'%\
					  (output_fname_prefix, no_of_combos, discard_perc), dpi=300)
			pylab.clf()
			pylab.title('Histogram of Frequency of Least-freq Combo For SNPPairs with %s combos'%no_of_combos)
			pylab.ylabel('frequency')
			pylab.hist(counter_obj.lowest_freq_ls, no_of_bins, alpha=0.5, normed=normed)
			pylab.savefig('%s_%s_combos_discard_perc_%s_lowest_freq_hist.png'%\
					  (output_fname_prefix, no_of_combos, discard_perc), dpi=300)
"""
input_fname = os.path.expanduser('~/panfs/250k/InterSNPCount/SNPpair_7_FT_22C.tsv')
output_fname_prefix = os.path.expanduser('%s'%os.path.splitext(input_fname)[0])
calSNPPairStats(input_fname, output_fname_prefix, discard_perc=0.4, no_of_bins=200, normed=True, debug=True)

"""
			

class Trio(object):
	"""
	2007-03-18
		the input_fname is output of MpiTrioAncestryInference.py
		for each strain, it outputs no_of_pairs, no_of_distinct_parents, avg_no_of_jumps
	"""
	def TrioAncestrySummary(input_fname, output_fname):
		import sys, os, csv
		from sets import Set
		sys.stderr.write("Reading trio ancestry information ...")
		reader = csv.reader(open(input_fname), delimiter='\t')
		strain_index2data = {}	#value is list of [no_of_pairs, Set_of_parents, no_of_jump_list]
		for row in reader:
			i,j,k = row[:3]
			no_of_jumps = int(row[-1])
			i = int(i)
			j = int(j)
			k = int(k)
			if k not in strain_index2data:
				strain_index2data[k] = [0, Set(), []]
			strain_index2data[k][0] += 1
			strain_index2data[k][1].add(i)
			strain_index2data[k][1].add(j)
			strain_index2data[k][2].append(no_of_jumps)
		del reader
		sys.stderr.write('done\n')
		
		sys.stderr.write("Outputting summary data ...")
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		for strain_index, data in strain_index2data.iteritems():
			no_of_distinct_parents = len(data[1])
			avg_no_of_jumps = sum(data[2])/float(data[0])
			writer.writerow([strain_index, data[0], no_of_distinct_parents, avg_no_of_jumps])
		del writer
		sys.stderr.write('done\n')
	
	"""
	TrioAncestrySummary('./script/variation/data/justin_data_filtered.trio_ances', './script/variation/data/justin_data_filtered.trio_ances.summary')
	"""
	
	"""
	2007-03-18
	"""
	def draw_histogram_of_data_from_specified_column(input_fname, column=-1, no_of_bins=20):
		import pylab, csv
		data_ls = []
		reader = csv.reader(open(input_fname), delimiter='\t')
		for row in reader:
			data_ls.append(float(row[column]))
		pylab.clf()
		pylab.hist(data_ls, no_of_bins)
		pylab.show()
	
	"""
	draw_histogram_of_data_from_specified_column('./script/variation/data/justin_data_filtered.trio_ances.summary', 1)
	"""
	
	"""
	2007-03-19
	"""
	def get_one_ancestry_line(trio_ances_fname, triplet):
		import csv, sys
		sys.stderr.write("Getting one ancestry line...")
		reader = csv.reader(open(trio_ances_fname), delimiter='\t')
		ancestry_ls = None
		for row in reader:
			this_triplet = row[:3]
			this_triplet = map(int, this_triplet)
			if this_triplet[0] == triplet[0] and this_triplet[1]==triplet[1] and this_triplet[2] == triplet[2]:
				ancestry_ls = row[3:-2]	#-2 to discard the last chromosome separator and no_of_jumps
				break
		sys.stderr.write("Done.\n")
		del reader
		return ancestry_ls
	
	def DrawTrioAncestry(data_matrix_fname, trio_ances_fname, triplet=[]):
		"""
		2007-04-17 double the layers to accomodate the heterozygotes
		"""
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(data_matrix_fname)
		if not triplet:
			print 'Error: triplet is not specified. '
			return
		ancestry_ls = get_one_ancestry_line(trio_ances_fname, triplet)
		import Numeric
		sub_data_matrix = Numeric.zeros([12, len(ancestry_ls)])
		from common import nt2number, number2nt
		no_of_chrs = 0
		for i in range(len(ancestry_ls)):
			if ancestry_ls[i]=='||':
				no_of_chrs += 1
			else:
				#first fill up the last 6 encoding lines
				sub_data_matrix[6,i] = 1
				sub_data_matrix[7,i] = 1
				sub_data_matrix[8,i] = 2
				sub_data_matrix[9,i] = 2
				ancestry_bit = int(ancestry_ls[i])
				if ancestry_bit == 2:
					sub_data_matrix[10,i] = 1
					sub_data_matrix[11,i] = 2
				else:
					sub_data_matrix[10,i] = sub_data_matrix[11,i] = ancestry_bit+1	#the original coding is 0 or 1, but 0 is reserved for NA
				for j in range(len(triplet)):	#second deal with the 1st 3 lines with real data
					SNP_call = data_matrix[triplet[j]][i-no_of_chrs]
					if SNP_call>4:
						SNP_call_1 = nt2number[number2nt[SNP_call][0]]
						sub_data_matrix[2*j,i] = SNP_call_1
						SNP_call_2 = nt2number[number2nt[SNP_call][1]]
						sub_data_matrix[2*j+1,i] = SNP_call_2
					else:
						sub_data_matrix[2*j,i] = sub_data_matrix[2*j+1,i] = SNP_call
		import pylab
		pylab.clf()
		pylab.imshow(sub_data_matrix, interpolation='nearest')
		pylab.colorbar()
		ytick_label_ls = []
		for i in range(2):	#repeat the labels
			for j in range(len(triplet)):
				ytick_label_ls.append(strain_acc_list[triplet[j]])	#2007-04-17 double it
				ytick_label_ls.append(strain_acc_list[triplet[j]])
		ytick_loc_ls = range(12)
		ytick_loc_ls.reverse()
		pylab.yticks(ytick_loc_ls, ytick_label_ls)
		pylab.show()
	
	"""
	2007-10-16
	"""
	def get_snp_index2pos(snp_acc_list, curs, snps_table='stock20071008.snps'):
		sys.stderr.write("Getting snp_index2pos ...")
		snp_index2pos = {}
		for i in range(len(snp_acc_list)):
			curs.execute("select chromosome, position from %s where snpid='%s'"%(snps_table, snp_acc_list[i]))
			rows = curs.fetchall()
			chromosome, position = rows[0]
			snp_index2pos[i] = (chromosome, position)
		sys.stderr.write("Done.\n")
		return snp_index2pos
	
	def get_shared_block_ls(row1, row2, snp_index2pos):
		shared_block_ls = []
		shared_block = []
		for i in range(len(row1)):
			if row1[i] == row2[i]:	#same allele
				if i==0 or snp_index2pos[i-1][0]==snp_index2pos[i][0]:	#1st snp or previous snp and this one are on same chromosome
					shared_block.append(i)
				elif shared_block:	#same allele, but previous snp is on different chromosome
					shared_block_ls.append(shared_block)
					shared_block = []	#clear it up
					shared_block.append(i)
			elif shared_block:	#different allele but shared_block is not empty
				shared_block_ls.append(shared_block)	#store the shared_block
				shared_block = []	#clear it up
		if row1[-1]==row2[-1]:	#last snp is same. then the shared_block is not put in yet.
			shared_block_ls.append(shared_block)
		return shared_block_ls
	
	def get_all_pairs_with_max_shared_block_length_ls(strain_acc_list, snp_index2pos, data_matrix):
		sys.stderr.write("Getting all_pairs_with_max_shared_block_length_ls ...")
		all_pairs_with_max_shared_block_length_ls = []
		no_of_strains = len(strain_acc_list)
		for i in range(no_of_strains):
			for j in range(i+1, no_of_strains):
				shared_block_ls = get_shared_block_ls(data_matrix[i], data_matrix[j], snp_index2pos)
				shared_block_length_ls = map(len, shared_block_ls)
				shared_block_length_ls.sort()
				all_pairs_with_max_shared_block_length_ls.append([shared_block_length_ls[-1], (i,j)])	#take the maximum shared block length
		all_pairs_with_max_shared_block_length_ls.sort()
		sys.stderr.write("Done.\n")
		return all_pairs_with_max_shared_block_length_ls
	
	def DrawTwoRowAndSharedBlock(data_matrix, strain_index_pair, strain_acc_list, snp_index2pos):
		"""
		2007-10-16
			used to check get_shared_block_ls() is correct
		"""
		import Numeric
		from common import nt2number, number2nt
		data_row1 = []
		data_row2 = []
		shared_block_row = []
		no_of_chrs = 0
		row1 = data_matrix[strain_index_pair[0]]
		row2 = data_matrix[strain_index_pair[1]]
		for i in range(len(row1)):
			#for shared_block_row: -1 is chromosome separator, 1 is same, 0 is different.
			if i!=0 and snp_index2pos[i-1][0]!=snp_index2pos[i][0]:	#a different chromosome
				shared_block_row.append(-1)
				data_row1.append(-1)	#-1 is chromosome separator
				data_row2.append(-1)
			if row1[i]==row2[i]:
				shared_block_row.append(1)
			else:
				shared_block_row.append(0)
			data_row1.append(row1[i])
			data_row2.append(row2[i])
		sub_data_matrix = Numeric.array([data_row1, data_row2, shared_block_row])
		import pylab
		pylab.clf()
		pylab.imshow(sub_data_matrix, interpolation='nearest')
		pylab.colorbar()
		ytick_label_ls = []
		for strain_index in strain_index_pair:
			ytick_label_ls.append(strain_acc_list[strain_index])
		ytick_label_ls.append('block')
		pylab.yticks([2,1,0], ytick_label_ls)	#it's reversed.
		pylab.show()
	
	
	def DrawSharedBlock_ls(shared_block_ls, snp_index2pos, chr_id2cumu_size):
		"""
		2008-02-01
			add comments
		2007-10-16
		"""
		import pylab
		pylab.clf()
		#draw the chromosome separator as green circles
		for chr_id, cumu_size in chr_id2cumu_size.iteritems():
			pylab.plot([cumu_size], [1], 'go')
		#draw the snp first as red sticks
		for snp_index,pos in snp_index2pos.iteritems():
			chr_id, chr_pos = pos
			cumu_chr_pos = chr_id2cumu_size[chr_id-1]+chr_pos
			pylab.plot([cumu_chr_pos], [1], 'r|')
		#draw the blocks as lines crossing those red sticks
		for shared_block in shared_block_ls:
			if len(shared_block)>1:
				starting_snp_index = shared_block[0]
				ending_snp_index = shared_block[-1]
				
				chr_id = snp_index2pos[starting_snp_index][0]
				starting_chr_pos = snp_index2pos[starting_snp_index][1]
				ending_chr_pos = snp_index2pos[ending_snp_index][1]
				#add the cumulative chromosome size
				cumu_starting_chr_pos = chr_id2cumu_size[chr_id-1]+starting_chr_pos
				cumu_ending_chr_pos = chr_id2cumu_size[chr_id-1]+ending_chr_pos
				#draw the block
				pylab.plot([cumu_starting_chr_pos, cumu_ending_chr_pos], [1,1], c='b')
		#chromosome separator on the x axis
		xtick_label_ls = []
		xtick_loc_ls = []
		chr_id_ls = chr_id2cumu_size.keys()
		chr_id_ls.sort()
		for i in range(1,len(chr_id_ls)):
			chr_id = chr_id_ls[i]
			xtick_loc_ls.append(chr_id2cumu_size[chr_id])
			xtick_label_ls.append('%s:%s'%(chr_id, chr_id2size[chr_id]))
		pylab.xticks(xtick_loc_ls, xtick_label_ls)	#it's reversed.
		pylab.show()

	"""
	#2007-10-16 check shared block between pairwise strains
	
	input_fname = 'script/variation/stock20071008/data_d110_c0_5_d001.tsv'
	from FilterStrainSNPMatrix import FilterStrainSNPMatrix
	FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
	header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
	snp_acc_list = header[2:]
	snp_index2pos = get_snp_index2pos(snp_acc_list, curs, snps_table='stock20071008.snps')
	all_pairs_with_max_shared_block_length_ls = get_all_pairs_with_max_shared_block_length_ls(strain_acc_list, snp_index2pos, data_matrix)
	
	#2007-10-16 test to draw the diagrams
	from common import get_chr_id2size, get_chr_id2cumu_size
	chr_id2size = get_chr_id2size(curs)
	data_ls = get_chr_id2cumu_size(chr_id2size)
	chr_id2cumu_size, chr_gap = data_ls[:2]
	shared_block_ls = get_shared_block_ls(data_matrix[1], data_matrix[2], snp_index2pos)
	DrawSharedBlock_ls(shared_block_ls, snp_index2pos, chr_id2cumu_size)
	
	"""
	"""
	2007-03-21
		The stat is d(i,j)-|d(i,k)-d(j,k)| for k is a child of i and j
		the bigger this value is, the more divergent between two parents and child is in equal distance to two parents
	"""
	def cal_trio_stat(trio_ances_fname, distance_matrix, output_fname, pic_fname_prefix, need_savefig=0, report=0):
		import csv
		import Numeric
		no_of_strains = distance_matrix.shape[0]
		triplet2trio_stat = {}
		trio_stat_ls = []
		reader = csv.reader(open(trio_ances_fname), delimiter='\t')
		counter = 0
		for row in reader:
			triplet = row[:3]
			triplet = map(int, triplet)
			i,j,k = triplet
			trio_stat = distance_matrix[i,j] - abs(distance_matrix[i,k] - distance_matrix[j,k])
			triplet2trio_stat[tuple(triplet)] = trio_stat
			trio_stat_ls.append(trio_stat)
			counter += 1
			if report and counter%3000 == 0:
				sys.stderr.write("%s%s"%('\x08'*20, counter))
		if report:
			sys.stderr.write("%s%s\n"%('\x08'*20, counter))
		del reader
		sys.stderr.write("outputting trio_stat...")
		writer = csv.writer(open(output_fname,'w'), delimiter='\t')
		for triplet, trio_stat in triplet2trio_stat.iteritems():
			writer.writerow(list(triplet)+[trio_stat])
		del writer
		sys.stderr.write("Done\n")
		import pylab
		pylab.clf()
		pylab.hist(trio_stat_ls, 20)
		pylab.title("hist of d(i,j)-|d(i,k)-d(j,k)|")
		if need_savefig:
			pylab.savefig('%s_trio_stat_his.eps'%pic_fname_prefix, dpi=300)
			pylab.savefig('%s_trio_stat_hist.svg'%pic_fname_prefix, dpi=300)
			pylab.savefig('%s_trio_stat_hist.png'%pic_fname_prefix, dpi=300)
		pylab.show()
		return triplet2trio_stat
	
	"""
	triplet2trio_stat = cal_trio_stat('./script/variation/data/justin_data_filtered.trio_ances', distance_matrix, './script/variation/data/justin_data_filtered.trio_ances.trio_stat', './script/variation/data/justin_data_filtered.trio_ances', need_savefig=1, report=1)
	"""

class Data2010(object):
	"""
	2007-11-08
		get 192 strains for suzi and draw them on the map
	"""
	def get_ecotypeid_ls_nativename_ls_of_192_strains(curs):
		ecotypeid_ls = []
		nativename_ls = []
		curs.execute("select b.ecotypeid, e.nativename from ecotype e, batch_ecotype b where b.ecotypeid=e.id and b.batchid=2")
		rows = curs.fetchall()
		for row in rows:
			ecotypeid, nativename = row
			ecotypeid_ls.append(ecotypeid)
			nativename_ls.append(nativename)
		return ecotypeid_ls, nativename_ls
	
	
	"""
	ecotypeid_ls, nativename_ls = get_ecotypeid_ls_nativename_ls_of_192_strains(curs)
	ecotypeid2pos = get_ecotypeid2pos(curs, 'ecotype')
	draw_strains_on_map(ecotypeid_ls, ecotypeid2pos, pic_title='192 strains',  pic_area=[-130,10,140,70], output_fname_prefix='/tmp/suzi_192', label_ls=nativename_ls, need_collapse_strains_with_same_pos=0)
	
	"""
	
	@classmethod
	def getFRIAlignment(cls, output_fname, alignment_id=1843):
		"""
		2008-07-31
			output alignment in a format for DrawSNPMatrix.py to draw SNP matrix plot
		"""
		from variation.src.AtDB import AtDB, Sequence, Alignment
		db = AtDB(hostname='localhost')
		rows = Sequence.query.filter_by(alignment=alignment_id).order_by(Sequence.accession).all()
		"""
		#2008-07-31 output in fasta format
		outf = open(output_fname, 'w')
		is_target_alignment_outputted = 0
		for row in rows:
			if not is_target_alignment_outputted:
				outf.write('>ref\n')
				outf.write('%s\n'%row.alignment_obj.target)
				is_target_alignment_outputted = 1
			outf.write('>%s %s\n'%(row.accession, row.accession_obj.name))
			outf.write('%s\n'%(row.bases))
		del outf
		"""
		#2008-08-01 output in a yh SNP matrix format for DrawSNPMatrix.py
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		is_target_alignment_outputted = 0
		for row in rows:
			if not is_target_alignment_outputted:
				header_row = ['name', 1]
				one_row = ['ref', 1]
				for i in range(len(row.alignment_obj.target)):
					base = row.alignment_obj.target[i]
					one_row.append(base)
					header_row.append(i+1)
				writer.writerow(header_row)
				writer.writerow(one_row)
				is_target_alignment_outputted = 1
			one_row = ['%s %s'%(row.accession, row.accession_obj.name), 1]
			for base in row.bases:
				one_row.append(base)
			writer.writerow(one_row)
		del writer
	
	"""
	getFRIAlignment('/tmp/alignment_1843.fasta')
	"""

class JBDataGWA(object):
	"""
	2010-4-20
		split out of GWA(), functions written to analyze seasoning flowering time data from Justin Borevitz lab, collected by Postdoc Yan Li.
	"""
	@classmethod
	def removeRowsNotInTargetSNPData(cls, individualID_ls, lineID_ls):
		"""
		2010-3-22
			multiple individuals in individualID_ls share the same ID (lineID).
			The purpose of this function is to create a matrix mapping each individual to a line, 
				which would be used to create a new kinship matrix.
		"""
		sys.stderr.write("Creating individual-to-line incidence matrix ...")
		no_of_replicates = len(individualID_ls)
		no_of_lines = len(lineID_ls)
		import numpy
		Z = numpy.zeros([no_of_replicates, no_of_lines])	# incidence matrix
		no_of_replicates_with_line = 0
		for i in range(no_of_replicates):
			individualID = individualID_ls[i]
			for j in range(no_of_lines):
				if individualID==lineID_ls[j]:
					Z[i,j] = 1
					no_of_replicates_with_line += 1
					break
		sys.stderr.write("%s out of %s without any corresponding line found. Done.\n"%(no_of_replicates-no_of_replicates_with_line, no_of_replicates))
		return Z
	
	@classmethod
	def encodeEnvrionmentData(cls, data_in_str, encode_dict={}):
		"""
		2010-3-22
		"""
		encoded_data_ls = []
		for data_point in data_in_str:
			code = encode_dict.get(data_point)
			if code is None:
				sys.stderr.write("Warning: code is None in encoding.\n")
			encoded_data_ls.append(code)
		return encoded_data_ls
	
	@classmethod
	def extractCovariateDataFromJBData(cls, JBData, snp_id_to_be_included_ls=[], backupData=None,\
									env_variable_coding_ls = []):
		"""
		called by checkEpistasisInJBLabData()
		
		2010-4-20
			use env_variable_coding_ls to aggregate all encoding together and shrink the code.
				it's a list of (env_variable, new_variable_name, encode_dict), i.e. ('planting', 'planting', {'spring':0, 'summer':1})
		2010-3-31
			add backupData which is SNP-only. If one snp_id is not found in backupData, search JBData.
		2010-3-22
			extracting SNP data and environmental data from JBData
		"""
		sys.stderr.write("Extracting covariate data ...")
		import numpy
		from pymodule.SNP import SNPData
		
		snp_id_to_be_included_set = set(snp_id_to_be_included_ls)	# to remove replicates
		unique_snp_id_to_be_included_ls = list(snp_id_to_be_included_set)	# a unique list of snp ids
		unique_snp_id_to_be_included_ls.sort()	# in chromosomal order. hope so.
		
		new_col_id_ls = []
		transposed_data_matrix = []
		for col_id in unique_snp_id_to_be_included_ls:
			if backupData is not None:
				col_index = backupData.col_id2col_index.get(col_id)
			else:
				col_index = None
			if col_index is not None:
				new_col_id_ls.append(col_id)
				transposed_data_matrix.append(map(int, backupData.data_matrix[:,col_index]))
			else:
				col_index = JBData.col_id2col_index.get(col_id)
				if col_index is not None:
					new_col_id_ls.append(col_id)
					transposed_data_matrix.append(map(int, JBData.data_matrix[:,col_index]))
		
		#lat_col_index = JBData.col_id2col_index['lat']
		#new_col_id_ls.append('lat')
		#transposed_data_matrix.append(map(float, JBData.data_matrix[:,lat_col_index]))
		
		for env_variable, new_variable_name, encode_dict in env_variable_coding_ls:
			variable_col_index = JBData.col_id2col_index[env_variable]
			binary_ls =  cls.encodeEnvrionmentData(JBData.data_matrix[:,variable_col_index], encode_dict=encode_dict)
			new_col_id_ls.append(new_variable_name)
			transposed_data_matrix.append(binary_ls)
		
		transposed_data_matrix = numpy.array(transposed_data_matrix)
		data_matrix = numpy.transpose(transposed_data_matrix)
		
		covariateData = SNPData(row_id_ls=JBData.row_id_ls, col_id_ls=new_col_id_ls, data_matrix=data_matrix)
		sys.stderr.write("Done.\n")
		return covariateData
	
	@classmethod
	def generateBaseFormula(cls, snp_id_to_be_included_ls=[], environment_variate_name_ls = ['planting', 'loc', 'TopOrNot'],\
						interaction_snp_id_in_base_formula_ls = None, includeGXE=True, includeEXE=False, includeE=True,\
						includeG=True):
		"""
		2010-4-18
			argument includeE: whether to inlcude environmental variables or not
			includeGXE: whether to inlcude GXE interactions or not
			includeEXE: whether to inlcude EXE interactions or not
			includeG: whether to include SNPs or not
		2010-3-23
			a base formula involving snp, snpXE, EXE, EXEXE.
		"""
		base_formula = []
		
		for snp_id in snp_id_to_be_included_ls:
			if includeG:
				base_formula.append((snp_id,))	# append((snp_id)) has a bug. it's equal as append(snp_id)
			if includeGXE:
				# snp X environment  (only planting + loc)
				for environment_variate_name in environment_variate_name_ls[:2]:
					if environment_variate_name:	# 2010-3-26 not empty
						base_formula.append((snp_id, environment_variate_name))
		
		if includeE:
			# add covariates for all environments
			for environment_variate_name in environment_variate_name_ls:
				if environment_variate_name:	# 2010-3-26 not empty
					base_formula.append((environment_variate_name,))
		
		if includeEXE:

			# pair interaction among the first 3 environments
			no_of_environments = 3	# len(environment_variate_name_ls)
			for i in range(no_of_environments):
				for j in range(i+1, no_of_environments):
					if environment_variate_name_ls[i] and environment_variate_name_ls[j]:	# not empty
						base_formula.append((environment_variate_name_ls[i], environment_variate_name_ls[j]))
			
			# triple interaction among first 3 environments
			no_of_environments = 3
			for i in range(no_of_environments):
				for j in range(i+1, no_of_environments):
					for k in range(j+1, no_of_environments):
						if environment_variate_name_ls[i] and environment_variate_name_ls[j] and environment_variate_name_ls[k]:	# not empty
							base_formula.append((environment_variate_name_ls[i], environment_variate_name_ls[j],\
												environment_variate_name_ls[k]))
		
		# 2010-3-26
		if interaction_snp_id_in_base_formula_ls:
			no_of_interactions = len(interaction_snp_id_in_base_formula_ls)
			for i in range(no_of_interactions):
				snp1 = interaction_snp_id_in_base_formula_ls[i][0]
				snp2 = interaction_snp_id_in_base_formula_ls[i][1]
				base_formula.append((snp1, snp2))
				if includeGXE:
					# include 3-way interactions (SNP epistasis X environment) 
					for new_covariate_name in environment_variate_name_ls[:2]:	#only planting & loc
						if new_covariate_name:		# not empty
							base_formula.append((snp1, snp2, new_covariate_name))
				
		return base_formula
	
	@classmethod
	def generateVariateMatrix(cls, covariateData, base_formula):
		"""
		2010-3-23
			generate a covariate matrix based on base_formula
			
		"""
		sys.stderr.write("Generating covariate matrix ...")
		import numpy
		no_of_rows, no_of_cols = covariateData.data_matrix.shape
		genotype_matrix = numpy.zeros([no_of_rows, len(base_formula)], dtype=numpy.int)
		for k in range(len(base_formula)):
			covariate_ls = base_formula[k]
			if type(covariate_ls)==str:	# it's just one variate
				no_of_covariates = 1
				first_variate_name = covariate_ls
			else:
				no_of_covariates = len(covariate_ls)
				first_variate_name = covariate_ls[0]
			first_variate_index = covariateData.col_id2col_index[first_variate_name]
			covariate = covariateData.data_matrix[:, first_variate_index]
			covariate = covariate.reshape([no_of_rows, 1])
			for i in range(1, no_of_covariates):	# multiply all the other covariates
				new_covariate_name = covariate_ls[i]
				new_covariate_index = covariateData.col_id2col_index[new_covariate_name]
				new_covariate = covariateData.data_matrix[:, new_covariate_index]
				new_covariate = new_covariate.reshape([no_of_rows, 1])
				covariate = new_covariate*covariate
			genotype_matrix[:, k] = covariate[:,0]
		
		sys.stderr.write("Done.\n")
		return genotype_matrix
	
	@classmethod
	def _drawEffectPvalueBarChart(cls, axe_effect, axe_pvalue, covariate_name2index, variate_name_order_ls, coeff_list, \
						beta_interval_delta_ls, coeff_p_value_list, covariate_name=None, width=0.2, single_variate2x_pos={},\
						color='r', x_pos_offset=0, drawn_covariate_name_set=set(), starting_left=0):
		"""
		2010-4-6
			add drawn_covariate_name_set to know which have been plotted, which haven't
		2010-4-1
			called by drawEffectPvalueBarChart()
		"""
		sys.stderr.write("Drawing one set of barcharts for single variate with %s ..."%covariate_name)
		import math
		from pymodule import yh_matplotlib
		xticks_ls = []
		left_ls = []
		height_ls = []
		pvalue_height_ls = []
		interval_ls = []
		left = starting_left
		for single_variate in variate_name_order_ls:
			if covariate_name is not None:
				variate = tuple(list(single_variate) + [covariate_name])
			else:
				variate = single_variate
			drawn_covariate_name_set.add(variate)
			variate_index = covariate_name2index.get(variate)
			if variate_index is not None and variate_index<len(coeff_list):
				left += 1
				if single_variate in single_variate2x_pos:
					x_pos = single_variate2x_pos[single_variate]
				else:
					x_pos = left
					single_variate2x_pos[single_variate] = x_pos
				
				xticks_ls.append(variate)
				left_ls.append(x_pos+x_pos_offset)	# 2010-4-6 difference between left and x_pos is the x_pos_offset
				height_ls.append(coeff_list[variate_index])
				if len(beta_interval_delta_ls)>variate_index:	# 2010-4-23 beta_interval_delta_ls could be nothing
					beta_interval_delta = beta_interval_delta_ls[variate_index]
				else:
					beta_interval_delta = 0
				interval_ls.append(abs(beta_interval_delta))	# 2010-4-23 beta_interval_delta has to be positive
				pvalue = coeff_p_value_list[variate_index]
				if pvalue>0:
					log_pvalue = -math.log10(pvalue)
				elif pvalue==0:
					log_pvalue = 30
				else:
					log_pvalue = pvalue - 1
				pvalue_height_ls.append(log_pvalue)
		if len(left_ls)>0 and len(height_ls)>0:
			bar_effect = axe_effect.bar(left_ls, height_ls, width, color=color, yerr=interval_ls, linewidth=0)
			bar_pvalue = axe_pvalue.bar(left_ls, pvalue_height_ls, width, color=color, linewidth=0)
			yh_matplotlib.autoLabelBarChartWithHeight(axe_effect, bar_effect)
			yh_matplotlib.autoLabelBarChartWithHeight(axe_pvalue, bar_pvalue)
		else:
			bar_effect = None
			bar_pvalue = None
		from pymodule import PassingData
		return_data = PassingData(bar_effect=bar_effect, bar_pvalue=bar_pvalue, left_ls=left_ls, xticks_ls=xticks_ls)
		sys.stderr.write("Done.\n")
		return return_data
		
	@classmethod
	def drawEffectPvalueBarChart(cls, association_result, output_fname_prefix, covariate_name_tuple_ls=[],\
								variate_name_order_ls = [('1_3978063',), ('4_264496',), ('4_268809',), ('4_269962',),
														('5_3188328',), ('4_264496', '5_3188328'), ('4_268809', '5_3188328'), \
														('4_269962', '5_3188328'), ('4_429928',),\
														
														('5_25376551',), ('4_1356197',), ('4_199214',), ('4_387727',),\
														('4_158958',), ('5_18620282',), ('4_286905',), ('2_8516520',), ('3_9340928',),\
														
														('1_24345319',), ('1_24345319', '5_3188328'), \
														
														('planting',), ('loc',), ('TopOrNot',), ('MiddleOrNot',)],\
								snp_id2gene_name ={'1_3978063':'SFT1', '1_24345319':'FT', '2_8516520': 'SFT2',\
												'3_9340928': 'SFT3', '4_268809': 'Fler',\
												'4_269962': 'Fcol', \
												'4_199214': 'CRP', '4_264496':'FRI', '4_158958': 'SFT4.1',\
												'4_286905': 'SFT4.2', '4_387727': 'SFT4.3', '4_429928': 'SFT4.4', \
												'4_1356197': 'SFT4.5', '5_3188328': 'FLC', '5_18620282':'DOG1', \
												'5_25376551': 'SFT5', 'planting': 'S', 'loc':'L', 'TopOrNot':'T',\
												'MiddleOrNot': 'M', 'intercept':'int'},\
								drawIntercept=False):
		"""
		2010-4-1
			draw bar chart for effect size (beta), pvalue
			
			the order of covariates in covariate_name_tuple_ls have to match the order in pdata.coeff_list
		"""
		sys.stderr.write("Drawing barcharts for effect size and pvalue  ...\n")
		covariate_name2index = {}
		for i in range(len(covariate_name_tuple_ls)):
			covariate_name = covariate_name_tuple_ls[i]
			covariate_name = tuple(covariate_name)
			covariate_name2index[covariate_name] = i
		
		pdata = association_result
		snp_index = pdata.snp_index
		test_variate_name = snp_index
		S_square = getattr(pdata, 'S_square', None)
		beta_interval_delta_ls = getattr(pdata, 'beta_interval_delta_ls', None)
		if beta_interval_delta_ls is None:
			beta_interval_delta_ls = []
		coeff_p_value_list = getattr(pdata, 'coeff_p_value_list', None)
		if coeff_p_value_list is None:
			coeff_p_value_list = []
		coeff_list = getattr(pdata, 'coeff_list', None)
		if coeff_list is None:
			coeff_list = []
		
		import math
		import pylab
		pylab.clf()
		
		axe_effect = pylab.axes([0.05, 0.05, 0.9, 0.4], frameon=False)
		axe_effect.grid(True, alpha=0.3)
		axe_effect.set_title('effect size')
		axe_pvalue = pylab.axes([0.05, 0.55, 0.9, 0.4], frameon=False, sharex=axe_effect)
		axe_pvalue.grid(True, alpha=0.3)
		axe_pvalue.set_title('-log(pvalue)')
		
		width = 0.2
		legend_ls = []
		bar_effect_ls = []
		bar_pvalue_ls = []
		single_variate2x_pos = {}
		drawn_covariate_name_set = set()
		
		return_data = cls._drawEffectPvalueBarChart(axe_effect, axe_pvalue, covariate_name2index, variate_name_order_ls, \
							coeff_list, beta_interval_delta_ls, coeff_p_value_list, covariate_name=None, width=width, \
							single_variate2x_pos=single_variate2x_pos, color='r', \
							drawn_covariate_name_set=drawn_covariate_name_set)
		bar_effect = return_data.bar_effect
		bar_pvalue = return_data.bar_pvalue
		left_ls = return_data.left_ls
		xticks_ls = return_data.xticks_ls
		bar_effect_ls.append(bar_effect[0])
		bar_pvalue_ls.append(bar_pvalue[0])
		legend_ls.append('SNP/Environment')
		
		
		return_data = cls._drawEffectPvalueBarChart(axe_effect, axe_pvalue, covariate_name2index, variate_name_order_ls, \
							coeff_list, beta_interval_delta_ls, coeff_p_value_list, covariate_name='planting', width=width, \
							single_variate2x_pos=single_variate2x_pos, color='g', x_pos_offset=width, \
							drawn_covariate_name_set=drawn_covariate_name_set)
		bar_effect = return_data.bar_effect
		bar_pvalue = return_data.bar_pvalue
		if bar_effect and bar_pvalue:
			bar_effect_ls.append(bar_effect[0])
			bar_pvalue_ls.append(bar_pvalue[0])
			legend_ls.append('SNP/Environment X season')
		
		return_data = cls._drawEffectPvalueBarChart(axe_effect, axe_pvalue, covariate_name2index, variate_name_order_ls, \
							coeff_list, beta_interval_delta_ls, coeff_p_value_list, covariate_name='loc', width=width, \
							single_variate2x_pos=single_variate2x_pos, color='b', x_pos_offset=2*width, \
							drawn_covariate_name_set=drawn_covariate_name_set)
		bar_effect = return_data.bar_effect
		bar_pvalue = return_data.bar_pvalue
		if bar_effect and bar_pvalue:
			bar_effect_ls.append(bar_effect[0])
			bar_pvalue_ls.append(bar_pvalue[0])
			legend_ls.append('SNP/Environment X location')
		
		#axe_effect.legend(bar_effect_ls, legend_ls)
		axe_pvalue.legend(bar_pvalue_ls, legend_ls)
		
		
		if not drawIntercept:
			drawn_covariate_name_set.add(('intercept',))	# fake it as drawn already
		
		# draw all the remaining variates in covariate_name2index (triple environmental interaction and others if exist)
		undrawn_variate_name_set = set(covariate_name2index.keys()) - drawn_covariate_name_set
		undrawn_variate_name_ls = list(undrawn_variate_name_set)
		undrawn_variate_name_ls.sort()
		#triple_env_variate_tuple = tuple(environment_variate_name_ls[:3])
		starting_left = left_ls[-1]	# the previous rightmost left serves as the new starting point
		#x_pos = left_ls[-1] + 1
		#single_variate2x_pos[triple_env_variate_tuple] = x_pos	# manually set the x_pos to 1 + the last position
		return_data = cls._drawEffectPvalueBarChart(axe_effect, axe_pvalue, covariate_name2index, undrawn_variate_name_ls, \
							coeff_list, beta_interval_delta_ls, coeff_p_value_list, covariate_name=None, width=width, \
							single_variate2x_pos=single_variate2x_pos, color='k', x_pos_offset=0, \
							drawn_covariate_name_set=drawn_covariate_name_set, starting_left=starting_left)
		if return_data.bar_effect:
			left_ls.extend(return_data.left_ls)
			xticks_ls.extend(return_data.xticks_ls)
		
		# improve the xticks
		official_xticks_ls = []
		for variate_name_tuple in xticks_ls:
			official_name_tuple = []
			for variate_name in variate_name_tuple:
				official_name = snp_id2gene_name.get(variate_name)
				if official_name is None:
					official_name = variate_name
				official_name_tuple.append(official_name)
			final_official_name = ':'.join(official_name_tuple)
			official_xticks_ls.append(final_official_name)
		#pylab.xticks(left_ls, official_xticks_ls)
		
		axe_effect.set_xticks(left_ls)
		axe_effect.set_xticklabels(official_xticks_ls, fontsize='xx-small')
		axe_pvalue.set_xticklabels(official_xticks_ls, fontsize='xx-small')
		
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		sys.stderr.write("Done.\n")
	
	@classmethod
	def getPhenotypeLsOutOfJBData(cls, JBData, phenotypeColName='DTF', logPhenotype=False, data_type=float):
		"""
		2010-4-16
			call JBData.getColDataGivenColName() instead
		2010-4-5
			get the phenotype list out of JBData
		"""
		import math
		phenotype_ls = JBData.getColDataGivenColName(phenotypeColName, convert_type=float)
		if logPhenotype:
			phenotype_ls = map(math.log10, phenotype_ls)
		return phenotype_ls
	
	@classmethod
	def writeHeaderJBData(cls, output_fname, base_formula, GXE_environment_variate_name_ls=[], special_interaction_snp_id_ls=None, \
						run_genome_scan=True, extraVariateNameLs = ['S_square', 'var_perc', ]):
		"""
		2010-4-18
			add argument extraVariateNameLs
		2010-4-5
			called by checkEpistasisInJBLabData()
			This function outputs a header above the tsv output-table by checkEpistasisInJBLabData().
			Association.output_multi_variate_lm_results() would output the actual data.
		"""
		import csv
		writer = csv.writer(open(output_fname,'w'), delimiter='\t')
		base_formula_name_ls = ['-'.join(covariate) for covariate in base_formula]
		if run_genome_scan:
			if special_interaction_snp_id_ls:
				test_variate_name = '-'.join(special_interaction_snp_id_ls)
			else:
				test_variate_name = 'test_variate' 
			header = [test_variate_name, '%s-pvalue'%test_variate_name, '%s-beta'%test_variate_name, \
					'%s-interval_delta'%test_variate_name]
			# 2010-3-24 this is for the interaction between the test_variate & environments, if exists
			# add them in the end
			for new_covariate_name in GXE_environment_variate_name_ls:
				if new_covariate_name:
					base_formula_name_ls.append('%s-%s'%(test_variate_name, new_covariate_name))
		else:
			header = []
		for extra_variate_name in extraVariateNameLs:	#2010-4-18
			header.append(extra_variate_name)
		header.extend(['intercept-pvalue', 'intercept-beta', 'intercept-interval_delta'])
		
		for variate_name in base_formula_name_ls:
			header.append('%s-pvalue'%variate_name)
			header.append('%s-beta'%variate_name)
			header.append('%s-interval_delta'%variate_name)
		writer.writerow(header)
		return writer
	
	@classmethod
	def drawScanAroundFRIJBLabData(cls, db, input_fname, output_fname_prefix, need_svg=False, ylim_type=3,
								snp_to_highlight={'4_268809':'r', '4_269962':'b'}):
		"""
		2010-4-7
			input_fname is output by JBDataGWA.checkEpistasisInJBLabData()
			add argument ylim_type:
				1: ylim = ax.get_ylim(); ax.set_ylim([0, ylim[1]])
				2: ax.set_ylim([min_y, max_y])
			draw a pvalue vs chromosomal position plot, like those GWAS plots.
			Multi-chromosome functionality is taken into account upon writing but not tested.
			In reality, it's used to draw association pvalues around FRI region.
			Each cofactor would also have a genomic-scan plot, drawn in a different axe. 
		"""
		import csv
		from pymodule.utils import figureOutDelimiter
		
		sys.stderr.write("Reading data from %s ..."%(input_fname))
		if not os.path.isfile(input_fname):
			sys.stderr.write("Error. File doesn't exist.\n")
			return
		reader = csv.reader(open(input_fname,'r'), delimiter=figureOutDelimiter(input_fname))
		
		import re, math
		variable_name_pvalue_pattern = re.compile(r'(\w+)-pvalue')	# to find the variable name 
		header = reader.next()
		# find the column index which contains the pvalue for the variable
		variable_name2col_index = {}
		variable_name2chr2pvalue_ls = {}	# data structure to store the pvalues
		variable_name_in_order_ls = []	# the order the plots will be in
		for i in range(len(header)):
			header_name = header[i]
			re_result = variable_name_pvalue_pattern.search(header_name)
			if re_result:
				variable_name = re_result.group(1)
				variable_name_in_order_ls.append(variable_name)
				variable_name2col_index[variable_name] = i
				variable_name2chr2pvalue_ls[variable_name] = {}		
		
		chr2pos_ls = {}
		snp_to_highlight2chr_and_data_index = {}	# find out where the highlight SNPs are
		for row in reader:
			chr_pos = row[0]
			snp_id = chr_pos
			chr_pos = chr_pos.split('_')
			chr_pos = map(int, chr_pos)
			chr, pos = chr_pos[:2]
			if chr not in chr2pos_ls:
				chr2pos_ls[chr] = []
			chr2pos_ls[chr].append(pos)
			for variable_name, col_index in variable_name2col_index.iteritems():
				if col_index>=len(row):
					continue
				if chr not in variable_name2chr2pvalue_ls[variable_name]:
					variable_name2chr2pvalue_ls[variable_name][chr] = []
				pvalue = float(row[col_index])
				if pvalue>0:
					logPvalue = -math.log10(pvalue)
				else:
					logPvalue = -1
				variable_name2chr2pvalue_ls[variable_name][chr].append(logPvalue)
				
			if snp_id in snp_to_highlight:
				snp_to_highlight2chr_and_data_index[snp_id] = [chr, len(chr2pos_ls[chr])-1]
		del reader
		sys.stderr.write("Done.\n")
		
		if len(chr2pos_ls)>1:	#multiple chromosomes
			from variation.src.common import get_chr_id2size, get_chr_id2cumu_size
			chr_id_int2size = get_chr_id2size(db.metadata.bind)
			chr_id2cumu_size, chr_gap, chr_id_ls = get_chr_id2cumu_size(chr_id_int2size, chr_gap=0)
		
		# remove variables who don't have any real pvalue data
		variable_name_with_data_in_order_ls = []
		for variable_name in variable_name_in_order_ls:
			chr2pvalue_ls = variable_name2chr2pvalue_ls.get(variable_name)
			if chr2pvalue_ls:
				variable_name_with_data_in_order_ls.append(variable_name)
		variable_name_in_order_ls = variable_name_with_data_in_order_ls
		
		import pylab
		pylab.clf()
		
		fig = pylab.figure(figsize=(8,10))
		
		variable_name2axe = {}
		first_axe = None
		inter_axe_gap = 0.05
		axe_height = (1.0-inter_axe_gap)/len(variable_name_in_order_ls)-inter_axe_gap
		import numpy
		for i in range(len(variable_name_in_order_ls)):
			variable_name = variable_name_in_order_ls[i]
			axe = pylab.axes([0.1, inter_axe_gap*(i+1)+axe_height*i, 0.85, axe_height], frameon=False, sharex=first_axe)
			if first_axe is None:
				first_axe = axe
			axe.set_title(variable_name)
			axe.grid(True, alpha=0.3)
			
			chr_ls = chr2pos_ls.keys()
			chr_ls.sort()
			max_y = None
			min_y = None
			chr2pvalue_ls = variable_name2chr2pvalue_ls[variable_name]
			for chr in chr_ls:
				pos_ls = chr2pos_ls[chr]
				x_ls = numpy.array(pos_ls)
				if len(chr_ls)>1:	# add the chromosome size only if there are multiple chromosomes
					x_ls += chr_id2cumu_size[chr]-chr_id_int2size[chr]
				pvalue_ls = chr2pvalue_ls[chr]
				
				if pos_ls and pvalue_ls:
					if max_y is None:
						max_y = max(pvalue_ls)
					else:
						max_y = max(max_y, max(pvalue_ls))
					if min_y is None:
						min_y = min(pvalue_ls)
					else:
						min_y = min(min_y, min(pvalue_ls))
					axe.plot(x_ls, pvalue_ls, '.', color='c', markeredgewidth=0, markersize=5, alpha=0.8)
					
					for snp_id, chr_and_data_index in snp_to_highlight2chr_and_data_index.iteritems():
						chr_for_this_snp, data_index = chr_and_data_index[:2]
						if chr_for_this_snp == chr:
							snp_x_pos = x_ls[data_index]
							snp_pvalue = pvalue_ls[data_index]
							axe.plot([snp_x_pos], [snp_pvalue], '.', color=snp_to_highlight[snp_id], \
									markeredgewidth=0, markersize=5, alpha=0.8)
			if ylim_type==1:
				ylim = axe.get_ylim()
				axe.set_ylim([0, ylim[1]])
			elif ylim_type==2:
				axe.set_ylim([min_y, max_y])
			axe.set_ylabel("-log(P-value)")
			
			variable_name2axe[variable_name] = axe
		
		#separate each chromosome
		#for chr in chr_ls[:-1]:
		#	print chr
		#	ax.axvline(chr_id2cumu_size[chr], linestyle='--', color='k', linewidth=0.8)
		

		
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		if need_svg:
			pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
	"""
	#2010-4-7
	
	# 2010-4-7 draw genome scan around FRI plot for all data together
	start_snp_id = '4_1'
	stop_snp_id = '4_1700000'
	
	for logPhenotype in [True, False]:
		for snp_id_to_be_included_ls in [[], ['4_268809'], ['4_269962'], ['4_268809', '4_269962']]:
			cholesky_inverse_fname = None
			included_snp_ids = '_'.join(snp_id_to_be_included_ls)
			filename_extra_ls = [included_snp_ids, start_snp_id, stop_snp_id]
			if logPhenotype:
				filename_extra_ls.append('logPhenotype')
			input_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/DaysToFlower16replicates_around_FRI_cofactor_%s.tsv'%
											('_'.join(filename_extra_ls)))
			output_fname_prefix = os.path.splitext(input_fname)[0]
			JBDataGWA.drawScanAroundFRIJBLabData(db_250k, input_fname, output_fname_prefix)
	sys.exit(0)
	
	# 2010-4-7 draw genome scan around FRI plot for each environment
	loc_value_ls = ['spain', 'sweden']
	planting_value_ls = ['spring', 'summer']
	start_snp_id = '4_1'
	stop_snp_id = '4_1700000'
	
	for logPhenotype in [True, False]:
		for planting_value in planting_value_ls:
			for loc_value in loc_value_ls:
				for snp_id_to_be_included_ls in [[], ['4_268809'], ['4_269962'], ['4_268809', '4_269962']]:
					cholesky_inverse_fname = None
					included_snp_ids = '_'.join(snp_id_to_be_included_ls)
					filename_extra_ls = [included_snp_ids, start_snp_id, stop_snp_id]
					if logPhenotype:
						filename_extra_ls.append('logPhenotype')
					input_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/DaysToFlower16replicates_around_FRI_%s_%s_cofactor_%s.tsv'%
													(loc_value, planting_value, '_'.join(filename_extra_ls)))
					output_fname_prefix = os.path.splitext(input_fname)[0]
					JBDataGWA.drawScanAroundFRIJBLabData(db_250k, input_fname, output_fname_prefix)
	sys.exit(0)
	
	"""

	@classmethod
	def getCholeskyInverseData(cls, L_inverse_fname, kinship_output_fname=None, genotype_fname_to_generate_kinship=None,\
								JBData=None, logPhenotype=False, preEMMAXCofactor_ls=None, env_variable_coding_ls=None,\
								needKinship=False, vg=None, ve=None):
		"""
		2010-8-22
			return JBData as well in case it's changed in getExpandedKinship
		2010-4-23
			add argument needKinship. whether need a kinship matrix which matches the JBData.
				Usually it's more than just loading the matrix out of kinship_output_fname.
		2010-4-16 
			split out of checkEpistasisInJBLabData(). and called by it.
			If L_inverse_fname exits, get L_inverse_Data out of it right away.
			Otherwise, generate kinship, estimate variance, and get choleskty inverse, store the L_inverse in the file.
		"""
		sys.stderr.write("Getting inverse of cholesky decomposition of kinship ...")
		from pymodule.SNP import SNPData
		from Association import Association
		from pymodule import PassingData
		import numpy, os
		Z = None
		if os.path.isfile(L_inverse_fname):
			L_inverse_Data = SNPData(input_fname=L_inverse_fname, ignore_2nd_column=1, matrix_data_type=float, turn_into_array=1)
			L_inverse = L_inverse_Data.data_matrix
			if needKinship:
				newKinshipData, JBData, Z = Association.getExpandedKinship(kinship_output_fname, genotype_fname_to_generate_kinship, JBData)
			else:
				newKinshipData = None
		else:
			newKinshipData, JBData, Z = Association.getExpandedKinship(kinship_output_fname, genotype_fname_to_generate_kinship, JBData)
			if vg is None or ve is None:
				phenotype_ls = cls.getPhenotypeLsOutOfJBData(JBData, logPhenotype=logPhenotype,)
				non_NA_genotype_ls = []
				if preEMMAXCofactor_ls:
					JBDataEnvEncoded = cls.extractCovariateDataFromJBData(JBData, env_variable_coding_ls=env_variable_coding_ls)
					for i in range(len(preEMMAXCofactor_ls)):
						factorName = preEMMAXCofactor_ls[i]
						non_NA_genotype_ls.append(JBDataEnvEncoded.getColDataGivenColName(factorName))
				non_NA_genotype_ls = numpy.array(non_NA_genotype_ls)
				preEMMAX_data = Association.preEMMAX(phenotype_ls, newKinshipData.data_matrix, non_NA_genotype_ls=non_NA_genotype_ls, \
											debug=True, Z=None)
				variance_matrix = preEMMAX_data.variance_matrix
				non_NA_phenotype_ar = preEMMAX_data.non_NA_phenotype_ar
				L_inverse = preEMMAX_data.L_inverse
				
				"""
				#2010-9-7 plot the likelihood  and d(likelihood) curve
				one_emma_rs = preEMMAX_data.one_emma_rs
				logdelta = one_emma_rs.logdelta
				LL = one_emma_rs.LL
				dLL = one_emma_rs.dLL
				import pylab
				pylab.clf()
				pylab.plot(logdelta, LL)
				pylab.show()
				import pdb
				pdb.set_trace()
				pylab.clf()
				pylab.plot(logdelta, dLL)
				pylab.show()
				"""
			else:
				no_of_rows, no_of_cols = newKinshipData.data_matrix.shape
				variance_matrix = vg*newKinshipData.data_matrix + ve*numpy.identity(no_of_rows)
				L = numpy.linalg.cholesky(variance_matrix)
				L_inverse = numpy.linalg.inv(L)
		
			L_inverse_Data = SNPData(row_id_ls=JBData.row_id_ls, col_id_ls=JBData.row_id_ls, data_matrix=L_inverse)
			L_inverse_Data.tofile(L_inverse_fname)
		sys.stderr.write("Done.\n")
		return PassingData(L_inverse_Data=L_inverse_Data , kinship_matrix = newKinshipData.data_matrix, JBData = JBData, Z=Z)
	
	@classmethod
	def outputPhenotypePrediction(cls, output_fname, covariateData=None, phenotype_ls=[], \
								genotype_matrix=None, coeff_list=[], env_variable_reverse_code_ls=[], logPhenotype=False):
		"""
		2010-4-20
			called by checkEpistasisInJBLabData to output predicted, observed phenotypes and their difference
			
			covariateData is SNPData structure holding data from all the singleton covariates (no interactions),
				from which genotype_matrix is derived.
			phenotype_ls is the phenotype list (log-transformed or not depends on the argument logPhenotype passed from checkEpistasisInJBLabData)
			genotype_matrix is a matrix of all covariates excluding the intercept.
			coeff_list is a list of coefficients including the intercept.
			
			env_variable_reverse_code_ls is used to construct the col_id_ls, labels for all the columns.
		"""
		sys.stderr.write("Outputting predicted, observed phenotypes ...")
		import numpy, math
		# genotype_matrix doesn't include the all-1 column reserved for the intercept
		non_intercept_effect_array = numpy.dot(genotype_matrix, coeff_list[1:])	# get the effect size for every row, except the intercept
		effect_size_array = non_intercept_effect_array + coeff_list[0]	# add intercept to every row's effect
		if logPhenotype:	# if original phenotype is logged, then power it back
			effect_size_array = [math.pow(10, effect_size) for effect_size in effect_size_array]
			phenotype_ls = [math.pow(10, phenotype_value) for phenotype_value in phenotype_ls]
		
		no_of_environments = len(env_variable_reverse_code_ls)
		col_id_ls = []	# the names for the outputted columns (phenotype data)
		
		# enumerate combinations of different environments
		from pymodule.algorithm import simpleCandidateSetGenerator, enumerateCombinations
		f = lambda x: simpleCandidateSetGenerator(x, element_ls=range(2), max_solution_size=no_of_environments)
		environment_combination_ls = enumerateCombinations(f)
		env_combination2starting_col_index = {}
		for i in range(len(environment_combination_ls)):
			env_combination = environment_combination_ls[i]
			start_col_index = i*3
			env_combination2starting_col_index[tuple(env_combination)] = start_col_index
			single_env_name_ls = [env_variable_reverse_code_ls[k][2][env_combination[k]] for k in range(len(env_combination))]
			env_combination_name = ''.join(single_env_name_ls)
			col_id_ls.append('%s_%s'%(start_col_index+1, env_combination_name))	# add a number in the front for a fake phenotype id
			col_id_ls.append('%s_%s'%(start_col_index+2, env_combination_name+'Pred'))
			col_id_ls.append('%s_%s'%(start_col_index+3, env_combination_name+'Dif'))
		
		# index of environment variables in covariateData 
		env_index_in_covariateData = []
		for env_variable, new_variable_name, reverse_code_dict in env_variable_reverse_code_ls:
			col_index = covariateData.col_id2col_index[new_variable_name]
			env_index_in_covariateData.append(col_index)
		
		data_matrix = []
		row_id_ls = covariateData.row_id_ls
		no_of_rows = len(row_id_ls)
		new_row_id_ls = []
		row_id2row_index = {}
		for i in range(no_of_rows):
			row_id = row_id_ls[i]
			if row_id not in row_id2row_index:
				new_row_id_ls.append(row_id)
				row_id2row_index[row_id] = len(row_id2row_index)
				data_matrix.append(['NA']*len(col_id_ls))
			row_index = row_id2row_index.get(row_id)
			predictedPhenotype = effect_size_array[i]
			observedPhenotype = phenotype_ls[i]
			diffPhenotype = predictedPhenotype - observedPhenotype
			env_setting_for_this_row = [covariateData.data_matrix[i][col_index] for col_index in env_index_in_covariateData]
			start_col_index = env_combination2starting_col_index.get(tuple(env_setting_for_this_row))
			if start_col_index is not None:
				data_matrix[row_index][start_col_index] = observedPhenotype
				data_matrix[row_index][start_col_index+1] = predictedPhenotype
				data_matrix[row_index][start_col_index+2] = diffPhenotype
		
		from pymodule import SNPData
		phenotypeData = SNPData(row_id_ls=new_row_id_ls, col_id_ls=col_id_ls, data_matrix=data_matrix)
		phenotypeData.tofile(output_fname)
		sys.stderr.write("Done.\n")
	
	@classmethod
	def checkEpistasisInJBLabData(cls, phenotype_genotype_fname, genotype_fname_to_generate_kinship, output_fname, \
								vg=None, ve=None, kinship_fname=None, cholesky_inverse_fname=None, includeInteraction=False, \
								snp_id_to_be_included_ls=['1_24345319', '1_3978063', '2_8516520', '3_9340928', \
								'4_1356197',  '4_158958', '4_199214', '4_264496', '4_286905', '4_387727', \
								'4_429928',  '5_18620282', '5_25376551', '5_3188328', '4_268809', '4_269962'],\
								report = True, debug=False, with_base_formula=True,\
								interaction_snp_id_in_base_formula_ls=[],\
								special_interaction_snp_id_ls=[],\
								planting_value=None, loc_value=None, run_genome_scan=0, useYanSNPData=False, logPhenotype=False,\
								start_snp_id=None, stop_snp_id=None, drawIntercept=False, includeGXE=False, includeEXE=False,\
								addEnvAsCofactorInPreEMMAX=False, includeE=True, includeG=True, run_type=4):
		"""
		2010-4-23
			add argument run_type (same as Association.linear_model but default is 4.).
				1: pure_linear_model
				2: emma
				3: pure_linear_model via R
				4: generalized least square with specified variance matrix
		2010-4-15
		argument run_genome_scan
			1/True: genome-wide one-SNP scan with covariates from snp_id_to_be_included_ls
			0/False: just one model
			2: epistasis scan 
		argument addEnvAsCofactorInPreEMMAX:
			whether to add environmental variables (season, location, TopOrNot, MiddleOrNot) if available as cofactors in estimating vg & ve & L_inverse_Data
		start_snp_id & stop_snp_id are used to specify the range for genome scan.
			
		
		snp_id_to_be_included_ls=['X1_24345319', 'X1_3978063', 'X2_8516520', 'X3_9340928', \
						'X4_1356197',  'X4_158958', 'X4_199214', 'X4_264496', 'X4_286905', 'X4_387727', \
						'X4_429928',  'X5_18620282', 'X5_25376551', 'X5_3188328'],\
		
		2010-3-22
			do pairwise detection for all 
			
			argument planting_value & loc_value are used (if specified) to throw out non-matching rows (ecotypes)
		"""
		from pymodule.SNP import SNPData
		from Association import Association
		import numpy, math
		from pymodule.utils import addExtraToFilenamePrefix
		
		# construct a new output file name
		if start_snp_id:
			output_fname = addExtraToFilenamePrefix(output_fname, start_snp_id)
		if stop_snp_id:
			output_fname = addExtraToFilenamePrefix(output_fname, stop_snp_id)
		if logPhenotype:
			output_fname = addExtraToFilenamePrefix(output_fname, "logPhenotype")
		if os.path.isfile(output_fname):
			sys.stderr.write("Error: File %s already exits.\n"%output_fname)
			return None
		
		JBData = SNPData(input_fname=phenotype_genotype_fname, turn_into_array=1, ignore_2nd_column=1, \
							data_starting_col=2, turn_into_integer=False)
		
		L_inverse_fname_additional_parts = []
		if planting_value:
			L_inverse_fname_additional_parts.append(planting_value)
			JBData = SNPData.keepRowsWhoseOneColMatchValue(JBData, 'planting', planting_value)
			environment_variate_name_ls = ['']	# excluding planting variate if it's already filtered.
		else:
			environment_variate_name_ls = ['planting']
		if loc_value:
			L_inverse_fname_additional_parts.append(loc_value)
			JBData = SNPData.keepRowsWhoseOneColMatchValue(JBData, 'loc', loc_value)
			environment_variate_name_ls.append('')	# excluding loc variate if it's already filtered.
		else:
			environment_variate_name_ls.append('loc')
		environment_variate_name_ls.append('TopOrNot')
		environment_variate_name_ls.append('MiddleOrNot')
		#environment_variate_name_ls.append('lat')
		
		# 2010-9-2 temporary to exclude environmental variables altogether
		#environment_variate_name_ls = []
		
		if addEnvAsCofactorInPreEMMAX:
			preEMMAXCofactor_ls = []
			for env_name in environment_variate_name_ls:
				if env_name:
					preEMMAXCofactor_ls.append(env_name)
		else:
			preEMMAXCofactor_ls = None
		
		if kinship_fname:
			kinship_output_fname = kinship_fname
		else:
			kinship_output_fname = os.path.splitext(genotype_fname_to_generate_kinship)[0]+'_kinship.tsv'
		
		
		Z = None	# to check whether Z would be computed or not
		
		phenotype_genotype_base_fname = os.path.splitext(os.path.basename(phenotype_genotype_fname))[0]
		
		if cholesky_inverse_fname:
			L_inverse_fname = cholesky_inverse_fname
		else:
			if logPhenotype:
				L_inverse_fname_additional_parts.append('logPhenotype')
			L_inverse_fname = "%s_%s_L_inverse.tsv"%(os.path.splitext(kinship_output_fname)[0], \
													phenotype_genotype_base_fname)
			if preEMMAXCofactor_ls:
				L_inverse_fname_additional_parts.append("preEMMAXcofactor")
				L_inverse_fname_additional_parts.extend(preEMMAXCofactor_ls)
			if L_inverse_fname_additional_parts:
				L_inverse_fname = addExtraToFilenamePrefix(L_inverse_fname, '_'.join(L_inverse_fname_additional_parts))
		
		# 2010-4-20 a list denoting how to encode those environmental variables in JBData
		env_variable_coding_ls = [('loc', 'loc', {'spain':0, 'sweden':1}),\
								('planting', 'planting', {'spring':0, 'summer':1}),\
								('shelfHeight', 'TopOrNot', {'bottom':0, 'bottom-middle':0, 'top':1, 'top-middle':1}),\
								('shelfHeight', 'MiddleOrNot', {'bottom':0, 'bottom-middle':1, 'top':0, 'top-middle':1})]
		
		# 2010-9-2 temporary , remove env variables altogether
		#env_variable_coding_ls = []
		
		# 2010-4-21 the reverse of env_variable_coding_ls used in outputPhenotypePrediction
		env_variable_reverse_code_ls = []
		for env_variable, new_variable_name, encode_dict in env_variable_coding_ls:
			if new_variable_name=='TopOrNot':
				reverse_code_dict = {0:'Bot', 1:'Top'}
			elif new_variable_name=='MiddleOrNot':
				reverse_code_dict = {0:'', 1:'Mid'}
			else:	# planting and loc is just reverse of encode_dict
				reverse_code_dict = {}
				for key, value in encode_dict.iteritems():
					key = key[0].upper() + key[1:]	# upper-case the first letter
					reverse_code_dict[value] = key[:3]	#First 3 letters
			env_variable_reverse_code_ls.append((env_variable, new_variable_name, reverse_code_dict))
		
		if run_type==2:
			needKinship = True		# EMMA needs kinshipMatrix
		else:
			needKinship = False
		returnData = cls.getCholeskyInverseData(L_inverse_fname, kinship_output_fname, \
											genotype_fname_to_generate_kinship, \
											JBData=JBData, logPhenotype=logPhenotype, \
											preEMMAXCofactor_ls=preEMMAXCofactor_ls, \
											env_variable_coding_ls=env_variable_coding_ls,
											needKinship=needKinship,vg=vg, ve=ve)
		L_inverse_Data = returnData.L_inverse_Data
		kinship_matrix = returnData.kinship_matrix
		JBData = returnData.JBData	#2010-8-22
		Z = returnData.Z # 2010-9-1
		if useYanSNPData:
			selectSNPDataWithReplicates = None	# 2010-4-2 use SNP data from yan's file
			JBData = JBData.removeRowsNotInTargetSNPData(L_inverse_Data)
		elif not snp_id_to_be_included_ls:	#if snp_id_to_be_included_ls is empty, selectSNPData would be empty.:
			selectSNPDataWithReplicates = None
			JBData = JBData.removeRowsNotInTargetSNPData(L_inverse_Data)
		else:
			# 2010-3-31 create a snpData with only selected columns
			selectSNPData = SNPData(input_fname=genotype_fname_to_generate_kinship, turn_into_array=1, ignore_2nd_column=1, \
							col_id_key_set=set(snp_id_to_be_included_ls))
			selectSNPData, allele_index2allele_ls = selectSNPData.convert2Binary(row_id_as_major_allele="6909")	# Col-0 as the 
			# expand the selectSNPData to include replicates
			if Z is None:
				JBData = JBData.removeRowsNotInTargetSNPData(selectSNPData)
				Z = Association.createIndividualToLineIncidenceMatrix(JBData.row_id_ls, selectSNPData.row_id_ls)
			snp_matrix_with_replicates = numpy.dot(Z, selectSNPData.data_matrix)
			selectSNPDataWithReplicates = SNPData(row_id_ls=JBData.row_id_ls, col_id_ls=selectSNPData.col_id_ls, \
												data_matrix=snp_matrix_with_replicates)
			
			genotype_fname_to_generate_kinship_base_fname = os.path.splitext(os.path.basename(genotype_fname_to_generate_kinship))[0]
			selectSNPDataWithReplicates_genotype_fname = '%s_genotype_from_%s.tsv'%(os.path.splitext(output_fname)[0], \
																				genotype_fname_to_generate_kinship_base_fname)
			if not os.path.isfile(selectSNPDataWithReplicates_genotype_fname):
				selectSNPDataWithReplicates.tofile(selectSNPDataWithReplicates_genotype_fname)
		
		L_inverse = L_inverse_Data.data_matrix
		phenotype_ls = cls.getPhenotypeLsOutOfJBData(JBData, logPhenotype=logPhenotype, )
		
		phenotype_variance = numpy.var(phenotype_ls)	# 2010-4-18 to calculate the variance explained without cholesky transformation
		non_NA_phenotype_ar = numpy.array(phenotype_ls)
		non_NA_phenotype_ar = numpy.dot(L_inverse, non_NA_phenotype_ar)	# numpy.dot and numpy.inner has subtle difference.
		covariateData = cls.extractCovariateDataFromJBData(JBData, snp_id_to_be_included_ls, \
														backupData = selectSNPDataWithReplicates,\
														env_variable_coding_ls=env_variable_coding_ls)
		
		#if run_genome_scan:
		#	base_formula = []
		#	environment_variate_name_ls = ['']*len(environment_variate_name_ls)
		#else:
		base_formula = cls.generateBaseFormula(snp_id_to_be_included_ls, environment_variate_name_ls, \
											interaction_snp_id_in_base_formula_ls=interaction_snp_id_in_base_formula_ls,\
											includeGXE=includeGXE, includeEXE=includeEXE, includeE=includeE,\
											includeG=includeG)
		
		base_variate_matrix = cls.generateVariateMatrix(covariateData, base_formula)
		
		no_of_rows, no_of_covariates = covariateData.data_matrix.shape
		results = []
		counter = 0
		real_counter = 0
		
		no_of_snps = len(snp_id_to_be_included_ls)
		
		if includeGXE:
			GXE_environment_variate_name_ls = environment_variate_name_ls
		else:
			GXE_environment_variate_name_ls = []
		extraVariateNameLs = ['S_square', 'var_perc','var_perc_real', 'heritability',]	#2010-4-18 explicitly set the extra columns to be outputted
		writer = cls.writeHeaderJBData(output_fname, base_formula, \
							GXE_environment_variate_name_ls=GXE_environment_variate_name_ls,\
							special_interaction_snp_id_ls=special_interaction_snp_id_ls,\
							run_genome_scan=run_genome_scan, extraVariateNameLs=extraVariateNameLs)
		if run_genome_scan==1 or run_genome_scan==True:
			snpData = SNPData(input_fname=genotype_fname_to_generate_kinship, turn_into_array=1, ignore_2nd_column=1)
			newSnpData, allele_index2allele_ls = snpData.convert2Binary(row_id_as_major_allele="6909")
			snpData = newSnpData
			Z = Association.createIndividualToLineIncidenceMatrix(JBData.row_id_ls, snpData.row_id_ls)
			results = []
			no_of_cols = snpData.data_matrix.shape[1]
			if start_snp_id and stop_snp_id:
				start_snp = start_snp_id.split('_')
				start_snp = map(int, start_snp)
				stop_snp = stop_snp_id.split('_')
				stop_snp = map(int, stop_snp)
				from PlotGroupOfSNPs import PlotGroupOfSNPs
				snp_region_tup = [start_snp[0], start_snp[1], stop_snp[0], stop_snp[1],]
				top_snp_data = PlotGroupOfSNPs.getTopSNPData(snp_region_tup=snp_region_tup, chr_pos_ls=snpData.col_id_ls, \
															sortChrPosLs=False)
				snp_id_ls = top_snp_data.snp_id_ls
			else:
				snp_id_ls = snpData.col_id_ls
			for snp_id in snp_id_ls:
				j = snpData.col_id2col_index[snp_id]
				genotype_ls = snpData.data_matrix[:,j]
				genotype_matrix = genotype_ls.reshape([len(genotype_ls), 1])
				genotype_matrix = numpy.dot(Z, genotype_matrix)
				if base_formula:
					genotype_matrix = numpy.hstack([base_variate_matrix] + [genotype_matrix])
				snp_id = snpData.col_id_ls[j]
				
				pdata = Association.linear_model(genotype_matrix, non_NA_phenotype_ar, min_data_point=3, \
												snp_index=snp_id, \
									kinship_matrix=kinship_matrix, eig_L=None, run_type=run_type, \
									counting_and_NA_checking=True,\
									variance_matrix=None, lower_triangular_cholesky_inverse=L_inverse)
				
				# counting_and_NA_checking=True to get MAF and MAC. doesn't matter if it's False.
				if pdata is not None:
					results.append(pdata)
					Association.output_multi_variate_lm_results([pdata], writer, run_genome_scan=run_genome_scan,\
															extraVariateNameLs=extraVariateNameLs)
					real_counter += 1
				counter += 1
				if counter%2000==0:
					sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
			sys.stderr.write("%s\t%s\t%s\n"%('\x08'*40, counter, real_counter))
			return
		elif run_genome_scan==0 or run_genome_scan==False:
			genotype_matrix = base_variate_matrix	#numpy.hstack(covariate_ls)	
			pdata = Association.linear_model(genotype_matrix, non_NA_phenotype_ar, min_data_point=3, \
								snp_index=None, \
								kinship_matrix=kinship_matrix, eig_L=None, run_type=run_type, counting_and_NA_checking=True,\
								variance_matrix=None, lower_triangular_cholesky_inverse=L_inverse)
			if pdata is not None:
				geno_effect_var = numpy.var(numpy.dot(genotype_matrix, pdata.coeff_list[1:]))	# genotype_matrix doesn't include intercept.
				pdata.var_perc_real = geno_effect_var/phenotype_variance	#2010-4-18
				Association.output_multi_variate_lm_results([pdata], writer, run_genome_scan=run_genome_scan,\
														extraVariateNameLs=extraVariateNameLs)
				try:
					output_fname_prefix = os.path.splitext(output_fname)[0]
					cls.drawEffectPvalueBarChart(pdata, output_fname_prefix, \
											covariate_name_tuple_ls=[('intercept',)] + base_formula, drawIntercept=drawIntercept)
				except:
					sys.stderr.write("Error encountered during drawEffectPvalueBarChart().\n")
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
				# 2010-4-21
				output_fname = '%s_phenotype_predicted_logged.tsv'%output_fname_prefix
				if not os.path.isfile(output_fname):
					cls.outputPhenotypePrediction(output_fname, covariateData, phenotype_ls, \
									genotype_matrix, coeff_list=pdata.coeff_list, \
									env_variable_reverse_code_ls=env_variable_reverse_code_ls, logPhenotype=logPhenotype)
					# 2010-4-21 you could set logPhenotype above to False to output logged phenotypes
		elif run_genome_scan==2:	#2010-4-18 epistasis scan looking for significant interactions given the model
			if special_interaction_snp_id_ls:	#2010-3-26
				interaction_snp_id_ls = special_interaction_snp_id_ls
			else:
				interaction_snp_id_ls = snp_id_to_be_included_ls
			no_of_snps = len(interaction_snp_id_ls)
			output_fname_prefix = os.path.splitext(output_fname)[0]
			
			for i in range(no_of_snps):
				for j in range(i+1, no_of_snps):
					snp1 = interaction_snp_id_ls[i]
					snp2 = interaction_snp_id_ls[j]
					i = covariateData.col_id2col_index[snp1]
					j = covariateData.col_id2col_index[snp2]
					
					covariate_ls = []
					covariate1 = covariateData.data_matrix[:,i]
					covariate1 = covariate1.reshape([no_of_rows, 1])
					#covariate_ls.append(covariate1)
					
					covariate2 = covariateData.data_matrix[:,j]
					covariate2 = covariate2.reshape([no_of_rows, 1])
					#covariate_ls.append(covariate2)
					
					#if includeInteraction:
					interaction_covariate = covariate1*covariate2
					
					additional_variates_in_formula = []
					"""
					# 2010-3-24 include 3-way interactions (SNP epistasis X environment)
					for new_covariate_name in environment_variate_name_ls[:2]:	#only planting & loc
						if new_covariate_name:	# 2010-3-26 not empty
							new_covariate_index = covariateData.col_id2col_index[new_covariate_name]
							new_covariate = covariateData.data_matrix[:, new_covariate_index]
							new_covariate = new_covariate.reshape([no_of_rows, 1])
							covariate = new_covariate*interaction_covariate
							covariate_ls.append(covariate)
							additional_variates_in_formula.append((snp1, snp2, new_covariate_name))
					"""
					covariate_ls.append(interaction_covariate)
					additional_variates_in_formula.append((snp1, snp2))
					
					genotype_matrix = numpy.hstack([base_variate_matrix] + covariate_ls)	#numpy.hstack(covariate_ls)	
					pdata = Association.linear_model(genotype_matrix, non_NA_phenotype_ar, min_data_point=3, \
										snp_index='%s-%s'%(snp1, snp2), \
										kinship_matrix=kinship_matrix, eig_L=None, run_type=run_type, \
										counting_and_NA_checking=True, \
										variance_matrix=None, lower_triangular_cholesky_inverse=L_inverse)
					if pdata is not None:
						geno_effect_var = numpy.var(numpy.dot(genotype_matrix, pdata.coeff_list[1:]))	# genotype_matrix doesn't include intercept.
						pdata.var_perc_real = geno_effect_var/phenotype_variance	#2010-4-18
						results.append(pdata)
						
						Association.output_multi_variate_lm_results([pdata], writer, extraVariateNameLs=extraVariateNameLs)
						_output_fname_prefix = '%s_%s_%s'%(output_fname_prefix, snp1, snp2)
						cls.drawEffectPvalueBarChart(pdata, _output_fname_prefix, \
													covariate_name_tuple_ls=[('intercept',)] + base_formula + additional_variates_in_formula,\
													drawIntercept=drawIntercept)
						real_counter += 1
					
					counter += 1
					if report and counter%100==0:
						sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
				if debug and real_counter>0:
					break
			if report:
				sys.stderr.write("%s\t%s\t%s\n"%('\x08'*40, counter, real_counter))
		del writer
		sys.stderr.write("Done.\n")
		return
		#sys.exit(0)
	
	"""
		# 2010-3-22
		### test
		#phenotype_genotype_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF.tsv')
		#genotype_fname_to_generate_kinship = '/Network/Data/250k/tmp-yh/250k_data/call_method_49_core482_test.tsv'
		#output_fname = os.path.expanduser('/tmp/d4DTF_pairwise_epistasis_using_core482_test.tsv')
		#vg=1451.803; ve=331.4423;
		#JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, output_fname, \
		#							vg=vg, ve=ve, debug=True, includeInteraction=True)
		#sys.exit(0)
		
		###### generate kinship, estimate vg, ve, etc. using core482 from call 49
		phenotype_genotype_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF.tsv')
		genotype_fname_to_generate_kinship = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482.tsv')
		#kinship_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/K.tsv')
		kinship_fname = None
		output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF_core482_new_K_new_vg_ve_pairwise_cofactor.tsv')
		#vg=1451.803; ve=331.4423; 	# 2010-3-24 estimates from ~/script/variation/data/JBLabSeasonFlowering/K.tsv
		vg=None
		ve=None
		vg=1534.10171646	# 2010-3-24 estimates from ~/mnt/panfs/250k/dataset/call_method_49_core482.tsv
		ve=499.58683376
		#JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, output_fname, vg=vg, ve=ve, kinship_fname=kinship_fname)
		
		output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF_core482_new_K_new_vg_ve_full_model_plus_snp_environment_epistasis.tsv')
		JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, output_fname, vg=vg, ve=ve, \
				includeInteraction=True, kinship_fname=kinship_fname)
		sys.exit(0)
		
		
		###### use existing core469 kinship and old vg and ve estimates (from JBLabSeasonFlowering.R)
		kinship_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/K.tsv')
		output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF_core469_K_old_vg_ve_pairwise_cofactor.tsv')
		vg=1451.803; ve=331.4423;
		#vg=None
		#ve=None
		JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, output_fname, vg=vg, ve=ve, kinship_fname=kinship_fname)
		
		#output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF_core469_K_old_vg_ve_pairwise_epistasis.tsv')
		#JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, output_fname, vg=vg, ve=ve, \
		#		includeInteraction=True, kinship_fname=kinship_fname)
		
		
		# 2010-3-26 final full model with FRI-FLC, and FT-FLC
		cholesky_inverse_fname = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482_kinship_7224_replicates_L_inverse.tsv')
		output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF_core482_new_K_new_vg_ve_full_model_plus_FT_FLC_and_FRI_FLC_epistasis.tsv')
		interaction_snp_id_in_base_formula_ls = ['X4_264496', 'X5_3188328']	# FRI-FLC
		special_interaction_snp_id_ls = ['X1_24345319','X5_3188328']	# 2010-3-26 FT-FLC
		JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, output_fname, vg=vg, ve=ve, \
				includeInteraction=True, kinship_fname=kinship_fname, cholesky_inverse_fname=cholesky_inverse_fname,\
				interaction_snp_id_in_base_formula_ls = interaction_snp_id_in_base_formula_ls, \
				special_interaction_snp_id_ls=special_interaction_snp_id_ls)
		
		# 2010-4-1 final full model with FRI-FLC, and FT-FLC
		cholesky_inverse_fname = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482_kinship_7224_replicates_L_inverse.tsv')
		cholesky_inverse_fname = None
		vg=None
		ve=None
		output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/DaysToFlower16replicates_full_model_plus_FT_FLC_and_FRI_FLC_epistasis_with_snp_from_core482_with_FRI_del_impute.tsv')
		interaction_snp_id_in_base_formula_ls = ['4_264496', '5_3188328']	# FRI-FLC
		special_interaction_snp_id_ls = ['1_24345319','5_3188328']	# 2010-3-26 FT-FLC
		snp_id_to_be_included_ls=['1_24345319', '1_3978063', '2_8516520', '3_9340928', \
									'4_1356197',  '4_158958', '4_199214', '4_264496', '4_286905', '4_387727', \
									'4_429928',  '5_18620282', '5_25376551', '5_3188328']
		JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, output_fname, vg=vg, ve=ve,\
				snp_id_to_be_included_ls=snp_id_to_be_included_ls,\
				includeInteraction=True, kinship_fname=kinship_fname, cholesky_inverse_fname=cholesky_inverse_fname,\
				interaction_snp_id_in_base_formula_ls = interaction_snp_id_in_base_formula_ls, \
				special_interaction_snp_id_ls=special_interaction_snp_id_ls)
		
		
		# 2010-4-7 full model without GXE , GXG
		phenotype_genotype_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/DaysToFlower16replicates_tg_ecotypeid.tsv')
		genotype_fname_to_generate_kinship = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482_with_FRI_del_chr_order_one_time_impute_yu_format.tsv')
		kinship_fname = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482_kinship.tsv')	#two FRI allele shall not make a difference.
		cholesky_inverse_fname = None
		vg=None
		ve=None
		output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/DaysToFlower16replicates_full_model_without_GXE_GXG_core482_with_FRI_del_impute.tsv')
		interaction_snp_id_in_base_formula_ls = []
		snp_id_to_be_included_ls=['1_24345319', '1_3978063', '2_8516520', '3_9340928', \
								'4_1356197',  '4_158958', '4_199214', '4_264496', '4_286905', '4_387727', \
								'4_429928',  '5_18620282', '5_25376551', '5_3188328']
		JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, output_fname, vg=vg, ve=ve,\
				snp_id_to_be_included_ls=snp_id_to_be_included_ls,\
				includeInteraction=True, kinship_fname=kinship_fname, cholesky_inverse_fname=cholesky_inverse_fname,\
				interaction_snp_id_in_base_formula_ls = interaction_snp_id_in_base_formula_ls, \
				useYanSNPData=False, logPhenotype=True, run_genome_scan=False, drawIntercept=True)
		sys.exit(0)
		
		#2010-3-26 run genome scan phenotype ~ SNP
		cholesky_inverse_fname = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482_kinship_7224_replicates_L_inverse.tsv')
		output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF_core482_new_K_new_vg_ve_gwa.tsv')	
		JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, \
									output_fname, vg=vg, ve=ve, \
									kinship_fname=kinship_fname, \
									cholesky_inverse_fname = cholesky_inverse_fname,\
									run_genome_scan=True)
		sys.exit(0)
		
		#2010-4-6 run genome scan phenotype ~ SNP
		cholesky_inverse_fname = None
		output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/DaysToFlower16replicates_core482_gwa.tsv')
		JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, \
									output_fname, vg=vg, ve=ve, \
									kinship_fname=kinship_fname, \
									cholesky_inverse_fname = cholesky_inverse_fname,\
									run_genome_scan=True)
		sys.exit(0)
		
		
		# 2010-3-27 full model with old K with new vg and ve estimates
		vg=None
		ve=None
		vg=1451.80295945; ve=331.442279994; 	# 2010-3-24 estimates from ~/script/variation/data/JBLabSeasonFlowering/K.tsv
		kinship_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/K.tsv')
		output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF_core469_old_K_new_vg_ve_full_model_plus_FT_FLC_and_FRI_FLC_epistasis.tsv')
		interaction_snp_id_in_base_formula_ls = ['X4_264496', 'X5_3188328']	# FRI-FLC
		special_interaction_snp_id_ls = ['X1_24345319','X5_3188328']	# 2010-3-26 FT-FLC
		cholesky_inverse_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/K_new_vg_ve_L_inverse.tsv')
		JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, output_fname, vg=vg, ve=ve, \
				includeInteraction=True, kinship_fname=kinship_fname, cholesky_inverse_fname=cholesky_inverse_fname, \
				interaction_snp_id_in_base_formula_ls = interaction_snp_id_in_base_formula_ls, \
				special_interaction_snp_id_ls=special_interaction_snp_id_ls)
		sys.exit(0)
		
		##### 2010-3-26 add FRI-FLC interaction to 4 block-only full model
		interaction_snp_id_in_base_formula_ls = ['4_264496', '5_3188328']	# FRI-FLC
		special_interaction_snp_id_ls = ['1_24345319','5_3188328']	# 2010-3-26 FT-FLC
		vg=None
		ve=None
		loc_value_ls = ['spain', 'sweden']
		planting_value_ls = ['spring', 'summer']
		for planting_value in planting_value_ls:
			for loc_value in loc_value_ls:
				cholesky_inverse_fname = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482_kinship_%s_%s_L_inverse_1.tsv'%(loc_value, planting_value))
				output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF_core482_new_K_new_vg_ve_FT_FLC_and_FRI_FLC_epistasis_%s_%s.tsv'%(loc_value, planting_value))	
				JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, \
											output_fname, vg=vg, ve=ve, \
											includeInteraction=True, kinship_fname=kinship_fname, \
											cholesky_inverse_fname = cholesky_inverse_fname,\
											interaction_snp_id_in_base_formula_ls = interaction_snp_id_in_base_formula_ls, \
											special_interaction_snp_id_ls=special_interaction_snp_id_ls,\
											planting_value=planting_value, loc_value=loc_value)
		
		sys.exit(0)
		
		
		#### 2010-4-7 verify using Yan's newsst SNP data & kinship
		phenotype_genotype_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/data/a2.tsv')
		
		genotype_fname_to_generate_kinship = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482_with_FRI_del_chr_order_one_time_impute_yu_format.tsv')
		kinship_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/data/K1b.tsv')
		output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/results/a2DTF_K1b_full_model_without_GXE_GXG_YanData.tsv')
		vg=None
		ve=None
		cholesky_inverse_fname = None
		interaction_snp_id_in_base_formula_ls = []
		snp_id_to_be_included_ls=['1_24345319', '1_3978063', '2_8516520', '3_9340928', \
								'4_1356197',  '4_158958', '4_199214', '4_264496', '4_286905', '4_387727', \
								'4_429928',  '5_18620282', '5_25376551', '5_3188328']
		
		JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, output_fname, vg=vg, ve=ve,\
				snp_id_to_be_included_ls=snp_id_to_be_included_ls,\
				includeInteraction=True, kinship_fname=kinship_fname, cholesky_inverse_fname=cholesky_inverse_fname,\
				interaction_snp_id_in_base_formula_ls = interaction_snp_id_in_base_formula_ls, \
				useYanSNPData=True, logPhenotype=False, run_genome_scan=False, drawIntercept=True)
		sys.exit(0)
		
		#############################
		#####2010-4-28 final code
		# 2010-4-6
		phenotype_genotype_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF.tsv')
		#phenotype_genotype_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF_tg_ecotypeid.tsv')
		phenotype_genotype_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/DaysToFlower16replicates_tg_ecotypeid.tsv')
		
		genotype_fname_to_generate_kinship = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482.tsv')
		genotype_fname_to_generate_kinship = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482_with_FRI_del_impute.tsv')
		genotype_fname_to_generate_kinship = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482_with_FRI_del_chr_order_one_time_impute_yu_format.tsv')
		#kinship_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/K.tsv')
		kinship_fname = None
		kinship_fname = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482_kinship.tsv')	#two FRI allele shall not make a difference.
		output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF_core482_new_K_new_vg_ve_pairwise_cofactor.tsv')
		#vg=1451.803; ve=331.4423; 	# 2010-3-24 estimates from ~/script/variation/data/JBLabSeasonFlowering/K.tsv
		vg=None
		ve=None
		vg=1534.10171646	# 2010-3-24 estimates from ~/mnt/panfs/250k/dataset/call_method_49_core482.tsv
		ve=499.58683376
		#JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, output_fname, vg=vg, ve=ve, kinship_fname=kinship_fname)
		
		
		# 2010-4-1 final full model with FRI-FLC, and FT-FLC
		cholesky_inverse_fname = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482_kinship_7224_replicates_L_inverse.tsv')
		cholesky_inverse_fname = None
		vg=None
		ve=None
		#output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF_core482_intersection_yanSNPData_new_K_new_vg_ve_full_model_plus_FT_FLC_and_FRI_FLC_epistasis.tsv')
		#output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF_core482_new_K_new_vg_ve_full_model_plus_FT_FLC_and_FRI_FLC_epistasis.tsv')
		#output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF_full_model_without_GXE_YanSNP.tsv')
		interaction_snp_id_in_base_formula_ls = [('4_264496', '5_3188328'), ('1_24345319','5_3188328')]	# FRI-FLC, FT-FLC
		interaction_snp_id_in_base_formula_ls = [('1_24345319','5_3188328')]	# FT-FLC
		interaction_snp_id_in_base_formula_ls = []
		#interaction_snp_id_in_base_formula_ls = ['4_264496', '5_3188328']	# FRI-FLC
		#special_interaction_snp_id_ls = ['1_24345319','5_3188328']	# 2010-3-26 FT-FLC
		snp_id_to_be_included_ls=['1_24345319', '1_3978063', '2_8516520', '3_9340928', \
								'4_1356197',  '4_158958', '4_199214', '4_264496', '4_286905', '4_387727', \
								'4_429928',  '5_18620282', '5_25376551', '5_3188328']
		
		snp_id_to_be_included_ls=['1_24345319', '1_3978063', '2_8516520', '3_9340928', \
								'4_1356197',  '4_158958', '4_199214', '4_264496',  '4_268809', '4_269962', '4_286905', '4_387727', \
								'4_429928',  '5_18620282', '5_25376551', '5_3188328']
		snp_id_to_be_included_ls=['1_24345319', '1_3978063', '2_8516520', '3_9340928', \
								'4_1356197',  '4_268809', '4_269962', '5_18620282', '5_25376551', '5_3188328']	
		
		
		snp_id_to_be_included_ls=['1_24345319', '1_3978063', '2_8516520', '3_9340928', \
								'4_1356197',  '4_158958', '4_268809', '4_269962', '4_387727', \
								'5_18620282', '5_25376551', '5_3188328']
		
		#2010-4-23 use original EMMA, rather EMMAX
		run_type = 2
		output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/DaysToFlower16replicates_full_model_9_run_type_%s.tsv'%run_type)
		JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, output_fname, vg=vg, ve=ve,\
				snp_id_to_be_included_ls=snp_id_to_be_included_ls,\
				includeInteraction=True, kinship_fname=kinship_fname, cholesky_inverse_fname=cholesky_inverse_fname,\
				interaction_snp_id_in_base_formula_ls = interaction_snp_id_in_base_formula_ls, \
				useYanSNPData=False, logPhenotype=True, run_genome_scan=False, drawIntercept=True, includeGXE=True, \
				includeEXE=True, includeE=True, includeG=True, run_type=run_type)
		sys.exit(0)
		
		
		# 2010-4-16 include E cofactors in preEMMAX (estimating vg and ve)
		output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/DaysToFlower16replicates_full_model_9_preEMMAX_w_Ecofactors.tsv')	
		JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, output_fname, vg=vg, ve=ve,\
				snp_id_to_be_included_ls=snp_id_to_be_included_ls,\
				includeInteraction=True, kinship_fname=kinship_fname, cholesky_inverse_fname=cholesky_inverse_fname,\
				interaction_snp_id_in_base_formula_ls = interaction_snp_id_in_base_formula_ls, \
				useYanSNPData=False, logPhenotype=True, run_genome_scan=False, drawIntercept=True, includeGXE=False, \
				includeEXE=False, addEnvAsCofactorInPreEMMAX=True)
		output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/DaysToFlower16replicates_full_model_9_preEMMAX_w_Ecofactors_GXE.tsv')	
		JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, output_fname, vg=vg, ve=ve,\
				snp_id_to_be_included_ls=snp_id_to_be_included_ls,\
				includeInteraction=True, kinship_fname=kinship_fname, cholesky_inverse_fname=cholesky_inverse_fname,\
				interaction_snp_id_in_base_formula_ls = interaction_snp_id_in_base_formula_ls, \
				useYanSNPData=False, logPhenotype=True, run_genome_scan=False, drawIntercept=True, includeGXE=True, \
				includeEXE=True, addEnvAsCofactorInPreEMMAX=True)
		sys.exit(0)
		
		#### 2010-4-15 epistasis scan: run_genome_scan=2
		output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/DaysToFlower16replicates_full_model_9_epistasis_scan.tsv')	
		JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, output_fname, vg=vg, ve=ve,\
				snp_id_to_be_included_ls=snp_id_to_be_included_ls,\
				kinship_fname=kinship_fname, cholesky_inverse_fname=cholesky_inverse_fname,\
				interaction_snp_id_in_base_formula_ls = interaction_snp_id_in_base_formula_ls, \
				useYanSNPData=False, logPhenotype=True, run_genome_scan=2, drawIntercept=True, includeGXE=False, includeEXE=False)
		sys.exit(0)
		
		
		special_interaction_snp_id_ls = []
		loc_value_ls = ['spain', 'sweden']
		planting_value_ls = ['spring', 'summer']
		logPhenotype = True
		common_prefix = 'DaysToFlower16replicates_full_model_9'
		for planting_value in planting_value_ls:
			# take data from the same planting_value
			cholesky_inverse_fname = None
			output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/%s_%s.tsv'%(common_prefix, planting_value))
			JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, \
										output_fname, vg=vg, ve=ve, \
										snp_id_to_be_included_ls=snp_id_to_be_included_ls,\
										includeInteraction=True, kinship_fname=kinship_fname, \
										cholesky_inverse_fname = cholesky_inverse_fname,\
										interaction_snp_id_in_base_formula_ls = interaction_snp_id_in_base_formula_ls, \
										special_interaction_snp_id_ls=special_interaction_snp_id_ls,\
										planting_value=planting_value, loc_value=None,\
										logPhenotype=logPhenotype, run_genome_scan=False, drawIntercept=True, \
										includeGXE=False, includeEXE=True)
			for loc_value in loc_value_ls:
				# take data from the same loc_value
				cholesky_inverse_fname = None
				output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/%s_%s.tsv'%\
												(common_prefix, loc_value,))
				JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, \
											output_fname, vg=vg, ve=ve, \
											snp_id_to_be_included_ls=snp_id_to_be_included_ls,\
											includeInteraction=True, kinship_fname=kinship_fname, \
											cholesky_inverse_fname = cholesky_inverse_fname,\
											interaction_snp_id_in_base_formula_ls = interaction_snp_id_in_base_formula_ls, \
											special_interaction_snp_id_ls=special_interaction_snp_id_ls,\
											planting_value=None, loc_value=loc_value,\
											logPhenotype=logPhenotype, run_genome_scan=False, drawIntercept=True, \
											includeGXE=False, includeEXE=True)
				
				#cholesky_inverse_fname = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482_kinship_%s_%s_L_inverse_1.tsv'%(loc_value, planting_value))
				cholesky_inverse_fname = None
				output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/%s_%s_%s.tsv'%\
												(common_prefix, loc_value, planting_value))
				JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, \
											output_fname, vg=vg, ve=ve, \
											snp_id_to_be_included_ls=snp_id_to_be_included_ls,\
											includeInteraction=True, kinship_fname=kinship_fname, \
											cholesky_inverse_fname = cholesky_inverse_fname,\
											interaction_snp_id_in_base_formula_ls = interaction_snp_id_in_base_formula_ls, \
											special_interaction_snp_id_ls=special_interaction_snp_id_ls,\
											planting_value=planting_value, loc_value=loc_value,\
											logPhenotype=logPhenotype, run_genome_scan=False, drawIntercept=True, \
											includeGXE=False, includeEXE=True)
		
		sys.exit(0)
		
		
		# 2010-4-6 run genome scan around FRI for all environments together, log & non-log phenotype
		start_snp_id = '4_1'
		stop_snp_id = '4_1700000'
		for logPhenotype in [True, False]:
			for snp_id_to_be_included_ls in [['4_268809'], ['4_269962'], ['4_268809', '4_269962'], []]:
				cholesky_inverse_fname = None
				included_snp_ids = '_'.join(snp_id_to_be_included_ls)
				output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/DaysToFlower16replicates_around_FRI_cofactor_%s.tsv'%
												(included_snp_ids))
				JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, \
											output_fname, vg=vg, ve=ve, \
											snp_id_to_be_included_ls=snp_id_to_be_included_ls,\
											kinship_fname=kinship_fname,\
											cholesky_inverse_fname = cholesky_inverse_fname,\
											start_snp_id = start_snp_id, stop_snp_id = stop_snp_id,\
											logPhenotype=logPhenotype, run_genome_scan=True)
		sys.exit(0)
		
		# 2010-4-6 run genome scan around FRI for 4 environments, log & non-log phenotype
		loc_value_ls = ['spain', 'sweden']
		planting_value_ls = ['spring', 'summer']
		start_snp_id = '4_1'
		stop_snp_id = '4_1700000'
		logPhenotype = True
		for planting_value in planting_value_ls:
			for loc_value in loc_value_ls:
				for snp_id_to_be_included_ls in [[], ['4_268809'], ['4_269962'], ['4_268809', '4_269962']]:
					cholesky_inverse_fname = None
					included_snp_ids = '_'.join(snp_id_to_be_included_ls)
					output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/DaysToFlower16replicates_around_FRI_%s_%s_cofactor_%s.tsv'%
													(loc_value, planting_value, included_snp_ids))
					JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, \
												output_fname, vg=vg, ve=ve, \
												snp_id_to_be_included_ls=snp_id_to_be_included_ls,\
												kinship_fname=kinship_fname,\
												cholesky_inverse_fname = cholesky_inverse_fname,\
												start_snp_id = start_snp_id, stop_snp_id = stop_snp_id,\
												planting_value=planting_value, loc_value=loc_value, \
												logPhenotype=logPhenotype, run_genome_scan=True)
		sys.exit(0)
		
		
		
		##### 2010-3-26 add FRI-FLC interaction to 4 block-only full model
		snp_id_to_be_included_ls=['1_24345319', '1_3978063', '2_8516520', '3_9340928', \
									'4_1356197',  '4_158958', '4_199214', '4_264496', '4_286905', '4_387727', \
									'4_429928',  '5_18620282', '5_25376551', '5_3188328']
		interaction_snp_id_in_base_formula_ls = [('4_264496', '5_3188328'), ('1_24345319','5_3188328')]	# FRI-FLC, FT-FLC
		interaction_snp_id_in_base_formula_ls = []	#2010-4-7
		special_interaction_snp_id_ls = []	# 2010-3-26 
		vg=None
		ve=None
		loc_value_ls = ['spain', 'sweden']
		planting_value_ls = ['spring', 'summer']
		for planting_value in planting_value_ls:
			for loc_value in loc_value_ls:
				for logPhenotype in [True, False]:
					#cholesky_inverse_fname = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482_kinship_%s_%s_L_inverse_1.tsv'%(loc_value, planting_value))
					cholesky_inverse_fname = None
					output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/DaysToFlower16replicates_full_model_without_GXE_GXG_%s_%s.tsv'%\
													(loc_value, planting_value))
					JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, \
												output_fname, vg=vg, ve=ve, \
												snp_id_to_be_included_ls=snp_id_to_be_included_ls,\
												includeInteraction=True, kinship_fname=kinship_fname, \
												cholesky_inverse_fname = cholesky_inverse_fname,\
												interaction_snp_id_in_base_formula_ls = interaction_snp_id_in_base_formula_ls, \
												special_interaction_snp_id_ls=special_interaction_snp_id_ls,\
												planting_value=planting_value, loc_value=loc_value,\
												logPhenotype=logPhenotype, run_genome_scan=False, drawIntercept=True)
		
		sys.exit(0)
		
		#################### 2010-9-6 for PNAS revisions
		phenotype_genotype_fname = os.path.expanduser('~/Downloads/BLUP_381_spSpring.csv')
		phenotype_genotype_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering20100820/DaysToFlower16replicates_tg_ecotypeid.tsv')
		
		genotype_fname_to_generate_kinship = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering20100820/call_method_49_core482_with_FRI_del_chr_order_one_time_impute_yu_format.tsv')
		kinship_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering20100820/call_method_49_core482_kinship.tsv')	
		vg=None
		ve=None
		interaction_snp_id_in_base_formula_ls = []
		snp_id_to_be_included_ls=['1_24345319', '1_3978063', '2_8516520', '3_9340928', \
								'4_1356197',  '4_158958', '4_268809', '4_269962', '4_387727', \
								'5_18620282', '5_25376551', '5_3188328']
		
		special_interaction_snp_id_ls = []	# 2010-3-26 
		vg=None
		ve=None
		
		
		#snp_id_to_be_included_ls = []
		loc_value_ls = ['spain', 'sweden']
		planting_value_ls = ['spring', 'summer']
		for planting_value in planting_value_ls:
			for loc_value in loc_value_ls:
				for logPhenotype in [True]:
					#cholesky_inverse_fname = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482_kinship_%s_%s_L_inverse_1.tsv'%(loc_value, planting_value))
					cholesky_inverse_fname = None
					output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering20100820/DaysToFlower16replicates_noSNP_%s_%s.tsv'%\
													(loc_value, planting_value))
					JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, \
												output_fname, vg=vg, ve=ve, \
												snp_id_to_be_included_ls=snp_id_to_be_included_ls,\
												includeInteraction=True, kinship_fname=kinship_fname, \
												cholesky_inverse_fname = cholesky_inverse_fname,\
												interaction_snp_id_in_base_formula_ls = interaction_snp_id_in_base_formula_ls, \
												special_interaction_snp_id_ls=special_interaction_snp_id_ls,\
												planting_value=planting_value, loc_value=loc_value,\
												logPhenotype=logPhenotype, run_genome_scan=False, drawIntercept=True, run_type=4)
		
		sys.exit(0)
		
		# 2010-9-2 get snp_id from some file 
		snp_id_to_be_included_ls = []
		import csv
		reader = csv.reader(open(os.path.expanduser("~/Downloads/69SNPs.csv")))
		reader.next()
		for row in reader:
			snp_id_to_be_included_ls.append(row[0])
		
		for logPhenotype in [True]:
			#cholesky_inverse_fname = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482_kinship_%s_%s_L_inverse_1.tsv'%(loc_value, planting_value))
			cholesky_inverse_fname = None
			output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering20100820/DaysToFlower16replicates_69SNP.tsv')
			JBDataGWA.checkEpistasisInJBLabData(phenotype_genotype_fname, genotype_fname_to_generate_kinship, \
										output_fname, vg=vg, ve=ve, \
										snp_id_to_be_included_ls = snp_id_to_be_included_ls,\
										includeInteraction = True, kinship_fname = kinship_fname, \
										cholesky_inverse_fname = cholesky_inverse_fname,\
										interaction_snp_id_in_base_formula_ls = interaction_snp_id_in_base_formula_ls, \
										special_interaction_snp_id_ls = special_interaction_snp_id_ls,\
										planting_value = None, loc_value = None,\
										logPhenotype=logPhenotype, run_genome_scan=False, drawIntercept=True)
		sys.exit(0)
		
		
	"""
	
	
	@classmethod
	def splitJBd4DTFIntoPhenotype(cls, JBd4DTF_fname, output_fname):
		"""
		2010-3-30
			split JBlab's one-column DTF data into multi-column phenotype file (loc X planting X shelfHeight, 2 X 2 X 4)
				so that ~/script/variation/src/PutPhenotypeIntoDB.py can put it into db.
		"""
		sys.stderr.write("Splitting the phenotype from JBd4DTF_fname into multiple (=16) phenotypes (multi-column one output file) ... ")
		from pymodule.SNP import SNPData
		JBData = SNPData(input_fname=JBd4DTF_fname, turn_into_array=1, ignore_2nd_column=1, \
						    data_starting_col=2, turn_into_integer=False)
		ecotype_id2row_index = {}
		unique_ecotype_id_ls = []
		for row_id in JBData.row_id_ls:
			if row_id not in ecotype_id2row_index:
				unique_ecotype_id_ls.append(row_id)
				ecotype_id2row_index[row_id] = len(ecotype_id2row_index)
		
		phenotype_name2col_index = {}
		phenotype_name_ls = []
		phenotype_col_index = JBData.col_id2col_index['DTF']
		planting_col_index = JBData.col_id2col_index['planting']
		loc_col_index = JBData.col_id2col_index['loc']
		shelfHeight_col_index = JBData.col_id2col_index['shelfHeight']
		
		no_of_rows = len(unique_ecotype_id_ls)
		import numpy
		phenotype_matrix = numpy.copy(JBData.data_matrix[:no_of_rows, :1])	# dimension is not right, but to get the correct dtype.
		phenotype_matrix.resize([no_of_rows, 16])	# 16 columns due to 16 environmental combinations
		# use numpy.copy to preserve the old data type
		phenotype_matrix[:,:] = 'NA'	# default is 'NA'
		phenotype_ls = JBData.data_matrix[:, phenotype_col_index]
		for i in range(len(JBData.row_id_ls)):
			ecotype_id = JBData.row_id_ls[i]
			row_index = ecotype_id2row_index[ecotype_id]
			phenotype_name = '%s-%s-%s'%(JBData.data_matrix[i, loc_col_index], JBData.data_matrix[i, planting_col_index], JBData.data_matrix[i, shelfHeight_col_index])
			if phenotype_name not in phenotype_name2col_index:
				phenotype_name2col_index[phenotype_name] = len(phenotype_name2col_index)
				phenotype_name_ls.append(phenotype_name)
			col_index = phenotype_name2col_index[phenotype_name]
			phenotype_matrix[row_index, col_index] = phenotype_ls[i]
		
		phenotypeData = SNPData(row_id_ls=unique_ecotype_id_ls, col_id_ls=phenotype_name_ls, data_matrix=phenotype_matrix)
		phenotypeData.tofile(output_fname)
		
	"""
	# 2010-3-31
	JBd4DTF_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF.tsv')
	JBd4DTF_fname = os.path.expanduser('/tmp/DaysToFlower16replicates.csv')	# from google docs
	output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF_multi_phenotype_split.tsv')
	JBDataGWA.splitJBd4DTFIntoPhenotype(JBd4DTF_fname, output_fname)
	
	"""
	

class GWA(object):
	"""
	2009-4-29
		check the distance between top genes and its most significant nearby SNP
		top genes are partitioned into two categories, candidate and non-candidate 
	"""
	def cmpGeneSNPDistance(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, output_fname_prefix, \
						min_score=3.5, no_of_top_lines=None):
		
		sys.stderr.write("Comparing %s %s %s ...\n"%(call_method_id, phenotype_method_id, analysis_method_id))
		from GeneListRankTest import GeneListRankTest 
		candidate_gene_list = GeneListRankTest.dealWithCandidateGeneList(list_type_id)
		candidate_gene_set = set(candidate_gene_list)
		
		rm = Stock_250kDB.ResultsByGene.results_method
		rbg = Stock_250kDB.ResultsByGene.query.filter(rm.has(call_method_id=call_method_id)).\
				filter(rm.has(phenotype_method_id=phenotype_method_id)).filter(rm.has(analysis_method_id=analysis_method_id)).\
				filter_by(min_distance=40000).filter_by(get_closest=0).first()
		
		result_fname = getattr(rbg, 'filename', None)
		if result_fname is None or not os.path.isfile(result_fname):
			sys.stderr.write("%s doesn't exist.\n"%result_fname)
			return None
		import csv
		from pymodule import getColName2IndexFromHeader
		reader = csv.reader(open(result_fname), delimiter='\t')
		col_name2index = getColName2IndexFromHeader(reader.next())
		candidate_gene_snp_dist_ls = []
		non_candidate_gene_snp_dist_ls = []
		counter = 0
		no_of_candidate_genes = 0
		no_of_touches = 0
		for row in reader:
			gene_id = int(row[col_name2index['gene_id']])
			score = float(row[col_name2index['score']])
			disp_pos = int(row[col_name2index['disp_pos']])
			snps_id = int(row[col_name2index['snps_id']])
			snps_context = Stock_250kDB.SnpsContext.query.filter_by(snps_id=snps_id).filter_by(gene_id=gene_id).first()
			if snps_context.disp_pos_comment=='touch':	#touch, distance shall be zero
				disp_pos = 0
				no_of_touches += 1
			if no_of_top_lines is None and score<min_score:
				break
			elif counter>no_of_top_lines:
				break
			if gene_id in candidate_gene_set:
				candidate_gene_snp_dist_ls.append(abs(disp_pos))
				no_of_candidate_genes += 1
			else:
				non_candidate_gene_snp_dist_ls.append(abs(disp_pos))
			counter += 1
			if counter%5000==0:
				sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, no_of_touches))
		sys.stderr.write("%s\t%s\t%s\n"%('\x08'*40, counter, no_of_touches))
		
		import pylab
		pylab.clf()
		n1 = pylab.hist(candidate_gene_snp_dist_ls, min(max(10, no_of_candidate_genes/10),20), alpha=0.4, normed=1)
		n2 = pylab.hist(non_candidate_gene_snp_dist_ls, 20, alpha=0.4, normed=1, facecolor='r')
		pylab.title("Histogram of SNP-gene dist Top %s %s"%(no_of_top_lines, rbg.short_name))
		pylab.legend([n1[2][0], n2[2][0]], ['%s candidate genes'%no_of_candidate_genes, '%s non-candidate genes'%(counter-no_of_candidate_genes)])
		pylab.xlabel("gene-SNP distance")
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		#pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		sys.stderr.write("Done.\n")
		
		
	"""
call_method_id = 32
phenotype_method_id = 1
analysis_method_id = 7
list_type_id = 145
no_of_top_lines = 10000
output_fname_prefix = '/tmp/SNP_gene_dist_call_method_id_%s_phenotype_%s_analysis_%s_list_type_%s_Top%s'%\
			(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, no_of_top_lines)
cmpGeneSNPDistance(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, output_fname_prefix, no_of_top_lines=no_of_top_lines)
for phenotype_method_id in range(1,8):
	for analysis_method_id in [1,7]:
		output_fname_prefix = '/tmp/SNP_gene_dist_call_method_id_%s_phenotype_%s_analysis_%s_list_type_%s_Top%s'%\
			(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, no_of_top_lines)
		cmpGeneSNPDistance(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, output_fname_prefix, no_of_top_lines=no_of_top_lines)		
	"""
	
	"""
	2009-4-29
		check the number of SNPs per gene and see whether it's different between candidate and non-candidate 
			top genes are partitioned into two categories, candidate and non-candidate 
	"""
	def checkNoOfSNPsPerGene(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, snps_context_wrapper, output_fname_prefix, \
						min_score=3.5, no_of_top_lines=None):
		sys.stderr.write("Checking %s %s %s ...\n"%(call_method_id, phenotype_method_id, analysis_method_id))
		from GeneListRankTest import GeneListRankTest 
		candidate_gene_list = GeneListRankTest.dealWithCandidateGeneList(list_type_id)
		candidate_gene_set = set(candidate_gene_list)
		
		rm = Stock_250kDB.ResultsMethod.query.filter_by(call_method_id=call_method_id).\
			filter_by(phenotype_method_id=phenotype_method_id).filter_by(analysis_method_id=analysis_method_id).\
			filter_by(results_method_type_id=1).first()
		
		result_fname = getattr(rm, 'filename', None)
		if result_fname is None or not os.path.isfile(result_fname):
			sys.stderr.write("%s doesn't exist.\n"%result_fname)
			return None
		import csv
		from pymodule import PassingData
		genome_wide_result = GeneListRankTest.getResultMethodContent(rm, min_MAF=0.0)
		
		candidate_gene_snp_dist_ls = []
		non_candidate_gene_snp_dist_ls = []
		counter = 0
		no_of_candidate_genes = 0
		no_of_non_candidate_genes = 0
		if genome_wide_result is None:
			return None
		candidate_gene_id2no_of_snps = {}
		non_candidate_gene_id2no_of_snps = {}
		counter = 0
		for data_obj in genome_wide_result.data_obj_ls:
			score = data_obj.value
			snps_context_matrix = snps_context_wrapper.returnGeneLs(data_obj.chromosome, data_obj.position)
			if score<min_score:
				continue
			for snps_context in snps_context_matrix:
				snps_id, disp_pos, gene_id = snps_context
				if gene_id in candidate_gene_set:
					dc = candidate_gene_id2no_of_snps
					no_of_candidate_genes += 1
				else:
					dc = non_candidate_gene_id2no_of_snps
					no_of_non_candidate_genes += 1	
				if gene_id not in dc:
					dc[gene_id] = 0
				dc[gene_id] += 1
			counter += 1
			if counter%5000==0:
				sys.stderr.write("%s\t%s\t%s\t%s"%('\x08'*40, counter, no_of_candidate_genes, no_of_non_candidate_genes))
		sys.stderr.write("%s\t%s\t%s\t%s\n"%('\x08'*40, counter, no_of_candidate_genes, no_of_non_candidate_genes))
		import pylab
		pylab.clf()
		
		#recalculate two counters below because previous calculations are not unique gene
		no_of_candidate_genes = len(candidate_gene_id2no_of_snps)
		no_of_non_candidate_genes = len(non_candidate_gene_id2no_of_snps)
		n1 = pylab.hist(candidate_gene_id2no_of_snps.values(), min(max(10, no_of_candidate_genes/20),20), alpha=0.4, normed=1)
		n2 = pylab.hist(non_candidate_gene_id2no_of_snps.values(), 20, alpha=0.4, normed=1, facecolor='r')
		pylab.title("Histogram of #SNPs per gene (score>=%s) %s"%(min_score, rm.short_name))
		pylab.legend([n1[2][0], n2[2][0]], ['%s candidate genes'%no_of_candidate_genes, '%s non-candidate genes'%(no_of_non_candidate_genes)])
		pylab.xlabel("Number of SNPs per gene")
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		#pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		sys.stderr.write("Done.\n")
	
	
	"""
call_method_id = 32
phenotype_method_id = 1
analysis_method_id = 7
list_type_id = 145
min_score = 3
output_fname_prefix = '/tmp/no_of_SNPs_per_gene_call_method_id_%s_phenotype_%s_analysis_%s_list_type_%s_min_score_%s'%\
			(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, min_score)

from GeneListRankTest import GeneListRankTest, SnpsContextWrapper
min_distance = 40000
get_closest = 0
snps_context_picklef = '/Network/Data/250k/tmp-yh/snps_context/snps_context_g0_m40000'
snps_context_wrapper = GeneListRankTest.dealWithSnpsContextWrapper(snps_context_picklef, min_distance, get_closest)

checkNoOfSNPsPerGene(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, snps_context_wrapper, output_fname_prefix, min_score=min_score)
for phenotype_method_id in range(1,8):
	for analysis_method_id in [1,7]:
		if analysis_method_id==1:
			min_score = 5
		else:
			min_score = 3
		output_fname_prefix = '/tmp/no_of_SNPs_per_gene_call_method_id_%s_phenotype_%s_analysis_%s_list_type_%s_min_score_%s'%\
			(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, min_score)
		checkNoOfSNPsPerGene(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, snps_context_wrapper, output_fname_prefix, min_score=min_score)
	"""
	
	def simulatePvalue(curs, snps_table, output_fname):
		"""
		2008-01-30
			simulate uniform p-values for snps selected from db
		"""
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		curs.execute("select chromosome, position from %s order by chromosome, position"%snps_table)
		rows = curs.fetchall()
		import random
		for row in rows:
			chromosome, position = row
			writer.writerow([chromosome, position, random.expovariate(1)])
		del writer
	
	"""
	snps_table = 'snps_250k'
	output_fname = '/tmp/simulate.pvalue'
	simulatePvalue(curs, snps_table, output_fname)
	"""
	
	def minusLogPvalue(input_fname, output_fname):
		"""
		2008-02-14
			take log of pvalue in input_fname
		"""
		import csv
		reader = csv.reader(open(input_fname), delimiter='\t')
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		import math
		for row in reader:
			chromosome, position, pvalue = row
			pvalue = float(pvalue)
			if pvalue>0:
				pvalue = -math.log(float(pvalue))
				writer.writerow([chromosome, position, pvalue])
		del writer, reader
	
	def exponentialMinusLogPvalue(input_fname, output_fname):
		"""
		2008-02-14
			take log of pvalue in input_fname
		"""
		import csv
		reader = csv.reader(open(input_fname), delimiter='\t')
		reader.next()	#skip the header
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		import math
		for row in reader:
			chromosome, position, pvalue = row
			pvalue = float(pvalue)
			if pvalue>0:
				pvalue = math.exp(-float(pvalue))
				writer.writerow([chromosome, position, pvalue])
		del writer, reader
	
	@classmethod
	def estimateFDRofGWA(cls, results_id, output_fname_prefix, top_p_value_cutoff=1, lambda_gap=0.1):
		"""
		2008-11-04
			try estimating FDR on our genome association results
		"""
		from Stock_250kDB import Stock_250kDB, ResultsMethod
		rm = ResultsMethod.get(results_id)
		from GeneListRankTest import GeneListRankTest
		from pymodule import PassingData
		pd = PassingData(do_log10_transformation=False)
		gwr = GeneListRankTest.getResultMethodContent(rm, min_MAF=0, pdata=pd)
		pvalue_ls = [data_obj.value for data_obj in gwr.data_obj_ls]
		pvalue_ls.sort()
		from transfac.src.AnalyzeTRANSFACHits import AnalyzeTRANSFACHits
		AnalyzeTRANSFACHits_ins = AnalyzeTRANSFACHits()
		AnalyzeTRANSFACHits_ins.remove_top_p_values(pvalue_ls, top_p_value_cutoff)
		figure_fname = '%s_p_value_hist.png'%output_fname_prefix
		AnalyzeTRANSFACHits_ins.draw_pvalue_histogram(pvalue_ls, figure_fname)
		figure_fname = '%s_pi0Tolambda.png'%output_fname_prefix
		lambda_list, pi0_list = AnalyzeTRANSFACHits_ins.calculate_pi0_list(pvalue_ls, figure_fname, \
																		top_p_value_cutoff=top_p_value_cutoff, lambda_gap=lambda_gap)
		
		estimated_pi0 = AnalyzeTRANSFACHits_ins.estimate_pi0(lambda_list, pi0_list)
		
		AnalyzeTRANSFACHits_ins.cal_q_value_list(pvalue_ls, estimated_pi0, top_p_value_cutoff, output_fname_prefix)
		
	"""
	results_id = 2318
	output_fname_prefix = '/tmp/2318_FDR'
	GWA.estimateFDRofGWA(results_id, output_fname_prefix)
	"""
	
	def outputGWABetaPvalue(input_fname, output_fname, pos_index=1, need_beta=True, min_value_cutoff=None, do_log10_transformation=False):
		"""
		2008-01-04 upgrade it to output a specified beta or its pvalue
		2008-12-18 take a genome-wide-result file, replace score with beta1 and output into a new file
		"""
		from pymodule import getGenomeWideResultFromFile
		genome_wide_result = getGenomeWideResultFromFile(input_fname, min_value_cutoff, do_log10_transformation)
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = ['chromosome', 'position', 'score', 'MAF', 'MAC', 'genotype_var_perc', 'beta0']
		writer.writerow(header)
		for data_obj in genome_wide_result.data_obj_ls:
			value = None
			if data_obj.comment:
				beta_and_pvalue_ls = data_obj.comment.split(',')	#a list of betas are all separated by ','
				if len(beta_and_pvalue_ls)>pos_index:
					beta, pvalue = beta_and_pvalue_ls[pos_index].split(':')	#beta and its pvalue are separted by ':'
					if need_beta:
						value = abs(float(beta)*10)	#increase it 10 fold to match pvalue, also take absolute value
					else:
						value = abs(float(pvalue))
			if value is not None:
				row = [data_obj.chromosome, data_obj.position, value, data_obj.maf, data_obj.mac, data_obj.genotype_var_perc] + data_obj.extra_col_ls
				writer.writerow(row)
		del writer
	
	"""
	input_fname = '/Network/Data/250k/db/results/type_1/2898_results.tsv'	#LM_with_PC12 on LD
	output_fname = '/tmp/LM_with_PC12_on_LD_beta1.tsv'
	outputGWABetaPvalue(input_fname, output_fname)
	
	input_fname = '/Network/Data/250k/db/results/type_1/2736_results.tsv'	#Emma on LD
	output_fname = '/tmp/Emma_on_LD_beta1.tsv'
	outputGWABetaPvalue(input_fname, output_fname)
	
	#2008-01-04 check the gene-environ interaction pvalue
	input_fname = '/Network/Data/250k/tmp-yh/association_results/lm/call_method_17_y6_pheno_2_LD+V_1_LD.tsv'
	output_fname = '/Network/Data/250k/tmp-yh/association_results/lm/call_method_17_y6_pvalue_by_int_pheno_2_LD+V_1_LD.tsv'
	outputGWABetaPvalue(input_fname, output_fname, pos_index=3, need_beta=False)
	
	pheno_pair_ls = ['2_LD+V_1_LD', '1_LD_3_SD', '1_LD_4_SD+V', '2_LD+V_3_SD', '2_LD+V_4_SD+V', '4_SD+V_3_SD', '5_FT_10C_6_FT_16C', '5_FT_10C_7_FT_22C', '6_FT_16C_7_FT_22C', '39_JIC0W_42_JIC8W']
	common_input_fname = '/Network/Data/250k/tmp-yh/association_results/lm/call_method_17_y6'
	for pheno_pair in pheno_pair_ls:
		input_fname = '%s_pheno_%s.tsv'%(common_input_fname, pheno_pair)
		output_fname = '%s_pvalue_by_int_pheno_%s.tsv'%(common_input_fname, pheno_pair)
		outputGWABetaPvalue(input_fname, output_fname, pos_index=3, need_beta=False)
	"""
	
	@classmethod
	def subtractTwoGWAS(cls, input_fname1, input_fname2):
		"""
		2010-3-9
			return a new genome_wide_result = (-log(Pvalue) from input_fname1) - (-log(Pvalue) from input_fname2)
		"""
		from pymodule import getGenomeWideResultFromFile
		gwa1 = getGenomeWideResultFromFile(input_fname1, min_value_cutoff=None, do_log10_transformation=True,\
										construct_chr_pos2index=True)
		gwa2 = getGenomeWideResultFromFile(input_fname2, min_value_cutoff=None, do_log10_transformation=True,\
										construct_chr_pos2index=True)
		
		sys.stderr.write("Subtracting two GWASs ...")
		chr_pos_ls1 = gwa1.chr_pos2index.keys()
		chr_pos_ls2 = gwa2.chr_pos2index.keys()
		chr_pos_set = set(chr_pos_ls1)&set(chr_pos_ls2)	#intersection of all SNPs
		
		total_chr_pos_ls = list(chr_pos_set)
		total_chr_pos_ls.sort()
		no_of_total_snps = len(chr_pos_set)
		from pymodule import GenomeWideResult
		gwr = GenomeWideResult(name="%s-%s"%(gwa1.name, gwa2.name), construct_chr_pos2index=True, \
							construct_data_obj_id2index=True)
		for i in range(int(no_of_total_snps)):
			#for chr, pos in chr_pos_set:
			chr, pos = total_chr_pos_ls[i]
			data_obj1 = gwa1.get_data_obj_by_chr_pos(chr, pos)
			data_obj2 = gwa2.get_data_obj_by_chr_pos(chr, pos)
			data_obj1.value = data_obj1.value - data_obj2.value
			gwr.add_one_data_obj(data_obj1)
		sys.stderr.write("Done.\n")
		return gwr
	
	"""
	# 2010-3-9
	input_fname1 = '/Network/Data/250k/db/results/type_1/4023_results.tsv'	# LM_cofactor_G. orontii conidiophore_32
	input_fname2 = '/Network/Data/250k/db/results/type_1/4116_results.tsv'	# LM_G. orontii conidiophore_32
	gwr = GWA.subtractTwoGWAS(input_fname1, input_fname2)
	output_fname_prefix = '/tmp/LM_cofactor-LM_G. orontii conidiophore_32'
	GWA.drawGWANicer(db_250k, gwr, output_fname_prefix, min_value=None, ylim_type=2)
	
	input_fname1 = '/Network/Data/250k/db/results/type_1/4024_results.tsv'	# EMMA_cofactor_G. orontii conidiophore_32
	input_fname2 = '/Network/Data/250k/db/results/type_1/3860_results.tsv'	# EmmaTrans_Mildew sensitivity_32
	gwr = GWA.subtractTwoGWAS(input_fname1, input_fname2)
	output_fname_prefix = '/tmp/EMMA_cofactor-EMMA_G. orontii conidiophore_32'
	GWA.drawGWANicer(db_250k, gwr, output_fname_prefix, min_value=None, ylim_type=2)
	
	"""
	
	@classmethod
	def drawGWANicer(cls, db, genome_wide_result, output_fname_prefix, min_value=2.5, need_svg=False, ylim_type=1):
		"""
		2010-3-9
			if min_value is None, no filter.
			add argument ylim_type:
				1: ylim = ax.get_ylim(); ax.set_ylim([0, ylim[1]])
				2: ax.set_ylim([min_y, max_y])
		2008-1-11 draw nicer genome wide plots
		"""
		chr2xy_ls = {}
		for data_obj in genome_wide_result.data_obj_ls:
			if min_value and data_obj.value<min_value:	#2010-3-9
				continue
			chr = data_obj.chromosome
			if chr not in chr2xy_ls:
				chr2xy_ls[chr] = [[],[]]
			chr2xy_ls[chr][0].append(data_obj.position)
			chr2xy_ls[chr][1].append(data_obj.value)
		
		from variation.src.common import get_chr_id2size, get_chr_id2cumu_size
		chr_id_int2size = get_chr_id2size(db.metadata.bind)
		chr_id2cumu_size, chr_gap, chr_id_ls = get_chr_id2cumu_size(chr_id_int2size, chr_gap=0)
		
		import pylab
		pylab.clf()
		fig = pylab.figure(figsize=(10,2))
		#ax = pylab.axes()
		ax = fig.gca()
		import numpy
		chr_ls = chr2xy_ls.keys()
		chr_ls.sort()
		max_y = None
		min_y = None
		for chr in chr_ls:
			xy_ls = chr2xy_ls[chr]
			x_ls = numpy.array(xy_ls[0])
			x_ls += chr_id2cumu_size[chr]-chr_id_int2size[chr]
			if xy_ls:
				if max_y is None:
					max_y = max(xy_ls[1])
				else:
					max_y = max(max_y, max(xy_ls[1]))
				if min_y is None:
					min_y = min(xy_ls[1])
				else:
					min_y = min(min_y, min(xy_ls[1]))
				ax.plot(x_ls, xy_ls[1], '.', markeredgewidth=0, markersize=5, alpha=0.8)
		
		#separate each chromosome
		#for chr in chr_ls[:-1]:
		#	print chr
		#	ax.axvline(chr_id2cumu_size[chr], linestyle='--', color='k', linewidth=0.8)
		
		
		#draw the bonferroni line
		bonferroni_value = -math.log10(0.01/len(genome_wide_result.data_obj_ls))
		ax.axhline(bonferroni_value, linestyle='--', color='k', linewidth=0.8)
		
		#ax.set_ylabel("-log(P-value)")
		#ax.set_xlabel('Chromosomal Position')
		#ax.set_xlim([0, chr_id2cumu_size[chr_ls[-1]]])
		
		if ylim_type==1:
			ylim = ax.get_ylim()
			ax.set_ylim([0, ylim[1]])
		elif ylim_type==2:
			ax.set_ylim([min_y, max_y])
		
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		if need_svg:
			pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
	
	
	"""
		input_fname = '/Network/Data/250k/db/results/type_1/3018_results.tsv'	#KW on LD, call method 22
		from pymodule import getGenomeWideResultFromFile
		genome_wide_result = getGenomeWideResultFromFile(input_fname, min_value_cutoff=None, do_log10_transformation=True)
		output_fname_prefix = '/tmp/call_22_KW_on_LD'
		GWA.drawGWANicer(db_250k, genome_wide_result, output_fname_prefix)
		
		input_fname = '/Network/Data/250k/db/results/type_1/3025_results.tsv'	#Emma on LD, call method 22
		from pymodule import getGenomeWideResultFromFile
		genome_wide_result2 = getGenomeWideResultFromFile(input_fname, min_value_cutoff=None, do_log10_transformation=True)
		output_fname_prefix = '/tmp/call_22_Emma_on_LD'
		GWA.drawGWANicer(db_250k, genome_wide_result2, output_fname_prefix, min_value=1)
		
		input_fname = '/Network/Data/250k/db/results/type_1/3024_results.tsv'	#KW on FT_22C, call method 22
		from pymodule import getGenomeWideResultFromFile
		genome_wide_result3 = getGenomeWideResultFromFile(input_fname, min_value_cutoff=None, do_log10_transformation=True)
		output_fname_prefix = '/tmp/call_22_KW_on_7_FT_22C'
		GWA.drawGWANicer(db_250k, genome_wide_result3, output_fname_prefix)
		
		input_fname = '/Network/Data/250k/db/results/type_1/3031_results.tsv'	#Emma on FT_22C, call method 22
		from pymodule import getGenomeWideResultFromFile
		genome_wide_result4 = getGenomeWideResultFromFile(input_fname, min_value_cutoff=None, do_log10_transformation=True)
		output_fname_prefix = '/tmp/call_22_Emma_on_7_FT_22C'
		GWA.drawGWANicer(db_250k, genome_wide_result4, output_fname_prefix, min_value=1)
		
		input_fname_ls = ['/Network/Data/250k/db/results/type_1/729_results.tsv',\
						'/Network/Data/250k/db/results/type_1/3554_results.tsv',\
						'/Network/Data/250k/db/results/type_1/3881_results.tsv']
		output_fname_prefix_ls = [os.path.expanduser('~/doc/20090930EckerLabVisit/figures/call_7_KW_on_Germ10'),\
								os.path.expanduser('~/doc/20090930EckerLabVisit/figures/call_32_KW_on_Germ10'),\
								os.path.expanduser('~/doc/20090930EckerLabVisit/figures/call_43_KW_on_Germ10')]
		from pymodule import getGenomeWideResultFromFile
		for i in range(len(input_fname_ls)):
			input_fname = input_fname_ls[i]
			output_fname_prefix = output_fname_prefix_ls[i]
			genome_wide_result2 = getGenomeWideResultFromFile(input_fname, min_value_cutoff=None, do_log10_transformation=True)
			GWA.drawGWANicer(db_250k, genome_wide_result2, output_fname_prefix, min_value=1)
		
		# 2010-3-10
		input_fname1 = '/Network/Data/250k/db/results/type_1/4023_results.tsv'	# LM_cofactor_G. orontii conidiophore_32
		input_fname2 = '/Network/Data/250k/db/results/type_1/4116_results.tsv'	# LM_G. orontii conidiophore_32
		gwr = GWA.subtractTwoGWAS(input_fname1, input_fname2)
		output_fname_prefix = '/tmp/LM_cofactor-LM_G. orontii conidiophore_32'
		GWA.drawGWANicer(db_250k, gwr, output_fname_prefix, min_value=None, ylim_type=2)
		
		input_fname1 = '/Network/Data/250k/db/results/type_1/4024_results.tsv'	# EMMA_cofactor_G. orontii conidiophore_32
		input_fname2 = '/Network/Data/250k/db/results/type_1/3860_results.tsv'	# EmmaTrans_Mildew sensitivity_32
		gwr = GWA.subtractTwoGWAS(input_fname1, input_fname2)
		output_fname_prefix = '/tmp/EMMA_cofactor-EMMA_G. orontii conidiophore_32'
		GWA.drawGWANicer(db_250k, gwr, output_fname_prefix, min_value=None, ylim_type=2)
		
		#2010-10-23
		call_method_id=57
		query = Stock_250kDB.ResultsMethod.query.filter_by(call_method_id=call_method_id)
		output_dir = os.path.expanduser('~/doc/compbiophd/figures/deletionGWAS/')
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		from pymodule import getGenomeWideResultFromFile, PassingData
		pdata = PassingData(min_MAF=0.1)
		for row in query:
			genome_wide_result = getGenomeWideResultFromFile(row.filename, min_value_cutoff=None, \
															do_log10_transformation=True, pdata=pdata)
			output_fname_prefix = os.path.join(output_dir, 'deletionGWASPlot_c%s_p%s_a%s'%(row.call_method_id, \
																	row.phenotype_method_id, row.analysis_method_id))
			GWA.drawGWANicer(db_250k, genome_wide_result, output_fname_prefix, min_value=None, ylim_type=2)
		sys.exit(0)
	"""
	
	@classmethod
	def drawGWAPvalueHist(cls, db, genome_wide_result, output_fname_prefix, min_value=2.5, need_svg=False, ylim_type=1):
		"""
		2010-10-8
		"""
		pvalue_ls = []
		for data_obj in genome_wide_result.data_obj_ls:
			pvalue_ls.append(data_obj.value)
		
		import pylab
		pylab.clf()
		fig = pylab.figure(figsize=(10,8))
		#ax = pylab.axes()
		ax = fig.gca()
		
		ax.hist(pvalue_ls, 50, alpha=0.8)
		
		#separate each chromosome
		#for chr in chr_ls[:-1]:
		#	print chr
		#	ax.axvline(chr_id2cumu_size[chr], linestyle='--', color='k', linewidth=0.8)
		
		
		#ax.set_ylabel("-log(P-value)")
		#ax.set_xlabel('Chromosomal Position')
		#ax.set_xlim([0, chr_id2cumu_size[chr_ls[-1]]])
		
		
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		if need_svg:
			pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
	"""
		#2010-10-8
		rm = Stock_250kDB.ResultsMethod.query.filter_by(call_method_id=32).filter_by(phenotype_method_id=5).filter_by(analysis_method_id=1).first()
		from pymodule import getGenomeWideResultFromFile
		genome_wide_result = getGenomeWideResultFromFile(rm.filename, min_value_cutoff=None, do_log10_transformation=False)
		output_fname_prefix = '/tmp/pvalue_hist_call_32_KW_on_LD'
		GWA.drawGWAPvalueHist(db_250k, genome_wide_result, output_fname_prefix)
		sys.exit(0)
	"""
	def plot_maf_vs_pvalue(maf_vector, input_fname, do_log10_transformation=True):
		"""
		2008-06-27 read in pvalues from a file
		"""
		from GenomeBrowser import GenomeBrowser
		genome_wide_result = GenomeBrowser.getGenomeWideResultFromFile(input_fname, do_log10_transformation=do_log10_transformation)
		pvalue_ls = [genome_wide_result.data_obj_ls[i].value for i in range(len(genome_wide_result.data_obj_ls))]
		import pylab
		pylab.clf()
		pylab.plot(maf_vector, pvalue_ls, '.')
		pylab.show()
	
	"""
pvalue_fname = '/Network/Data/250k/db/results/type_1/394_results.tsv'
plot_maf_vs_pvalue(maf_vector, pvalue_fname)
for i in range(389, 758):
	pvalue_fname = '/Network/Data/250k/db/results/type_1/%s_results.tsv'%i
	plot_maf_vs_pvalue(maf_vector, pvalue_fname)
	"""
	
	
	@classmethod
	def plotCmpEnrichmentOfTwoAnalysisMethods(cls, db, output_fname_prefix):
		"""
		2009-4-17
			plot the enrichment ratios in each cell of the 2X2 partition between two analysis methods across all phenotypes
		"""
		import matplotlib
		import pylab, numpy
		pylab.clf()
		#Stock_250kDB.CmpEnrichmentOfTwoAnalysisMethods
		rows = db.metadata.bind.execute("select c.*, r.phenotype_method_id, r.call_method_id from %s c, %s r, %s p where c.results_id1=r.id and \
								p.id=r.phenotype_method_id and p.biology_category_id=1 order by r.phenotype_method_id, c.r1_min_score"%\
								(Stock_250kDB.CmpEnrichmentOfTwoAnalysisMethods.table.name, Stock_250kDB.ResultsMethod.table.name,\
								Stock_250kDB.PhenotypeMethod.table.name))
		axe_height = 0.2
		ax1 = pylab.axes([0.1, 0.1, 0.8, axe_height], frameon=False)
		axe_gap = 0.02
		ax2 = pylab.axes([0.1, 0.1+axe_height+axe_gap, 0.8, axe_height], frameon=False, sharey=ax1)
		ax3 = pylab.axes([0.1, 0.1+axe_height*2+axe_gap*2, 0.8, axe_height], frameon=False, sharey=ax1)
		ax4 = pylab.axes([0.1, 0.1+axe_height*3+axe_gap*3, 0.8, axe_height], frameon=False, sharey=ax1)
		non_significant_enrich_ratio_ls = []
		emma_enrich_ratio_ls = []
		kw_enrich_ratio_ls = []
		double_enrich_ratio_ls = []
		xlabel_ls = []
		for row in rows:
			if row.enrichment_ratio is None:
				enrichment_ratio = 0
			else:
				enrichment_ratio = row.enrichment_ratio
			if row.r1_max_score is None and row.r2_max_score is None:
				pm = Stock_250kDB.PhenotypeMethod.get(row.phenotype_method_id)
				xlabel_ls.append(pm.short_name)
				double_enrich_ratio_ls.append(enrichment_ratio)
			elif row.r1_max_score is None and row.r2_max_score is not None:
				emma_enrich_ratio_ls.append(enrichment_ratio)
			elif row.r1_max_score is not None and row.r2_max_score is None:
				kw_enrich_ratio_ls.append(enrichment_ratio)
			else:
				non_significant_enrich_ratio_ls.append(enrichment_ratio)
		N = len(xlabel_ls)
		ind = numpy.arange(N)    # the x locations for the groups
		width = 0.35       # the width of the bars: can also be len(x) sequenc
		p1 = ax1.bar(ind, double_enrich_ratio_ls, width)
		ax1.set_ylabel('double')
		p2 = ax2.bar(ind, emma_enrich_ratio_ls, width)
		ax2.set_ylabel('Emma')
		p3 = ax3.bar(ind, kw_enrich_ratio_ls, width)
		ax3.set_ylabel('KW')
		p4 = ax4.bar(ind, non_significant_enrich_ratio_ls, width)
		
		ax1.set_xticks(ind+width/2.,  xlabel_ls)
		ax1.set_ylim([0,5])
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)

	"""
output_fname_prefix = '/tmp/enrich_ratio'
plotCmpEnrichmentOfTwoAnalysisMethods(db_250k, output_fname_prefix)
	"""
	
	@classmethod
	def generateCommonSQLClause(cls, biology_category_id=1, threshold_type_id=12, type_id=56):
		"""
		2009-5-2
			used by plot3DCmpEnrichmentOfTwoAnalysisMethods() + cmpRankAtThresholdOfTwoAnalysisMethods()
		"""
		common_sql_sentence = "from %s c, %s r, %s p where c.results_id1=r.id and c.type_id=%s and \
							p.id=r.phenotype_method_id and p.biology_category_id=%s and c.threshold_type_id=%s order by \
							r.phenotype_method_id, c.results_id1, c.results_id2, c.r1_min_score"%\
							(Stock_250kDB.CmpEnrichmentOfTwoAnalysisMethods.table.name, Stock_250kDB.ResultsMethod.table.name,\
							Stock_250kDB.PhenotypeMethod.table.name, type_id, biology_category_id, threshold_type_id)
		return common_sql_sentence
	
	@classmethod
	def plot3DCmpEnrichmentOfTwoAnalysisMethods(cls, db, biology_category_id=1, threshold_type_id=12, type_id=56, \
											min_non_candidate_sample_size=20):
		"""
		2009-4-17
			3D-version plot of plotCmpEnrichmentOfTwoAnalysisMethods 
		"""
		import matplotlib
		import pylab, numpy
		pylab.clf()
		
		#Stock_250kDB.CmpEnrichmentOfTwoAnalysisMethods
		common_sql_sentence = cls.generateCommonSQLClause(biology_category_id, threshold_type_id, type_id)
		
		rows = db.metadata.bind.execute("select distinct r.phenotype_method_id %s"%common_sql_sentence)
		no_of_phenotypes = rows.rowcount*2
		
		rows = db.metadata.bind.execute("select c.*, r.phenotype_method_id, r.call_method_id %s"%common_sql_sentence)

		x, y = numpy.mgrid[0:2*no_of_phenotypes:1, 0:4:1]	#added a gap of 1 column between two phenotypes. one phenotype occupies two rows & two columns.
		
		#remove the gap in x & y
		needed_index_ls = []
		for i in range(0, no_of_phenotypes):
			needed_index_ls.append(2*i)
			#needed_index_ls.append(3*i+1)
			#y[3*i+1][1]=2
		x = x[needed_index_ls]
		y = y[needed_index_ls]
		enrichment_matrix = numpy.ones(x.shape, numpy.float)
		
		non_significant_enrich_ratio_ls = []
		emma_enrich_ratio_ls = []
		kw_enrich_ratio_ls = []
		double_enrich_ratio_ls = []
		xlabel_ls = []
		i = 0	#
		old_results_id_tuple = None
		for row in rows:
			results_id_tuple = (row.results_id1, row.results_id2)
			if old_results_id_tuple == None:
				old_results_id_tuple = results_id_tuple
				pm = Stock_250kDB.PhenotypeMethod.get(row.phenotype_method_id)
				xlabel_ls.append(pm.short_name)
			elif old_results_id_tuple!=results_id_tuple:
				i += 1
				old_results_id_tuple = results_id_tuple
				pm = Stock_250kDB.PhenotypeMethod.get(row.phenotype_method_id)
				xlabel_ls.append(pm.short_name)
			
			if row.enrichment_ratio is None or row.non_candidate_sample_size<=min_non_candidate_sample_size:
				enrichment_ratio = 1
			else:
				enrichment_ratio = row.enrichment_ratio
				
			if row.r1_max_score is None and row.r2_max_score is None:
				enrichment_matrix[i, 3] = enrichment_ratio
			elif row.r1_max_score is None and row.r2_max_score is not None:
				enrichment_matrix[i, 2] = enrichment_ratio
				#emma_enrich_ratio_ls.append(enrichment_ratio)
			elif row.r1_max_score is not None and row.r2_max_score is None:
				enrichment_matrix[i, 1] = enrichment_ratio
				#kw_enrich_ratio_ls.append(enrichment_ratio)
			else:
				enrichment_matrix[i, 0] = enrichment_ratio
				#non_significant_enrich_ratio_ls.append(enrichment_ratio)
		
		from enthought.mayavi import mlab
		mlab.clf()
		bar = mlab.barchart(x, y , enrichment_matrix, lateral_scale=0.8, opacity=1.)
		#mlab.ylabel("KW")
		#mlab.xlabel("Emma")
		#mlab.zlabel("Enrichment Ratio")
		from pymodule.DrawMatrix import get_font 
		font = get_font()
		
		for i in range(len(xlabel_ls)):
			label = xlabel_ls[i]
			char_width, char_height = font.getsize(label)	#W is the the biggest(widest)
			
			mlab.text(2*i, 0, label, z=0, width=char_width/1500.)	#min(0.0075*len(label), 0.04))
		
		s = numpy.ones(x.shape, numpy.int)
		surf = mlab.surf(x, y, s, opacity=0.6, extent=[-1, 2*no_of_phenotypes, -1, 4, 1.0, 1.0])
		
	
	"""
	biology_category_id = 1
	threshold_type_id = 1
	type_id = 57
	min_non_candidate_sample_size = 10
	output_fname_prefix = '/tmp/enrich_ratio_biology_category_id_%s_threshold_type_id_%s_type_id_%s'%(biology_category_id, threshold_type_id, type_id)
	#GWA.cmpRankAtThresholdOfTwoAnalysisMethods(db_250k, output_fname_prefix, biology_category_id, threshold_type_id, type_id)
	GWA.plot3DCmpEnrichmentOfTwoAnalysisMethods(db_250k, biology_category_id, threshold_type_id, type_id, min_non_candidate_sample_size)
	"""
	
	
	@classmethod
	def cmpRankAtThresholdOfTwoAnalysisMethods(cls, db, output_fname_prefix, biology_category_id=1, threshold_type_id=12, type_id=56):
		"""
		2009-5-1
		GWA based on two analysis methods are selected at different cutoffs to do enrichment.
			this function looks at the pvalue rank at that cutoff (=number of SNPs above that cutoff) in different phenotypes.
			to see how different two analysis methods differ in terms of that.
			
			an output table of all values allows further invesitgation of whether the number of SNPs above that cutoff
			has any correlation with enrichment ratio or other correlations.
		"""
		sys.stderr.write("Comparing ...")
		import matplotlib
		import pylab, numpy
		pylab.clf()
		
		#Stock_250kDB.CmpEnrichmentOfTwoAnalysisMethods
		common_sql_sentence = cls.generateCommonSQLClause(biology_category_id, threshold_type_id, type_id)
		
		rows = db.metadata.bind.execute("select distinct r.phenotype_method_id %s"%common_sql_sentence)
		no_of_phenotypes = rows.rowcount
		
		rows = db.metadata.bind.execute("select c.*, r.phenotype_method_id, r.call_method_id %s"%common_sql_sentence)
		
		x, y = numpy.mgrid[0:3*no_of_phenotypes:1, 0:2:1]	#added a gap of 1 column between two phenotypes. one phenotype occupies two rows & two columns.
		
		non_significant_enrich_ratio_ls = []
		emma_enrich_ratio_ls = []
		kw_enrich_ratio_ls = []
		double_enrich_ratio_ls = []
		
		non_significant_sample_size_ls = []
		emma_sample_size_ls = []
		kw_sample_size_ls = []
		double_sample_size_ls = []
		xlabel_ls = []
		i = 0	#
		old_results_id_tuple = None
		r1_threshold = None
		r2_threshold = None
		for row in rows:
			results_id_tuple = (row.results_id1, row.results_id2)
			if old_results_id_tuple == None:
				old_results_id_tuple = results_id_tuple
			elif old_results_id_tuple!=results_id_tuple:
				i += 1
				old_results_id_tuple = results_id_tuple
			if row.candidate_sample_size is None:
				candidate_sample_size = 0
			else:
				candidate_sample_size = row.candidate_sample_size
				
			if row.non_candidate_sample_size is None:
				non_candidate_sample_size = 0
			else:
				non_candidate_sample_size = row.non_candidate_sample_size
			if row.enrichment_ratio is None:
				enrichment_ratio = 1
			else:
				enrichment_ratio = row.enrichment_ratio
			
			#if enrichment_ratio==0:
			#	enrichment_ratio==0.01
			sample_size = candidate_sample_size+non_candidate_sample_size
			if row.r1_max_score is None and row.r2_max_score is None:
				pm = Stock_250kDB.PhenotypeMethod.get(row.phenotype_method_id)
				xlabel_ls.append(pm.short_name)
				double_sample_size_ls.append(sample_size)
				double_enrich_ratio_ls.append(enrichment_ratio)
			elif row.r1_max_score is None and row.r2_max_score is not None:
				emma_sample_size_ls.append(sample_size)
				emma_enrich_ratio_ls.append(enrichment_ratio)
			elif row.r1_max_score is not None and row.r2_max_score is None:
				kw_sample_size_ls.append(sample_size)
				kw_enrich_ratio_ls.append(enrichment_ratio)
			else:
				r1_threshold = row.r1_max_score
				r2_threshold = row.r2_max_score
				non_significant_sample_size_ls.append(sample_size)
				non_significant_enrich_ratio_ls.append(enrichment_ratio)
				#enrichment_matrix[2*i, 0] = enrichment_ratio
		import numpy
		double_sample_size_ls = numpy.array(double_sample_size_ls, numpy.float)
		emma_sample_size_ls = numpy.array(emma_sample_size_ls, numpy.float)
		kw_sample_size_ls = numpy.array(kw_sample_size_ls, numpy.float)
		real_emma_sample_size_ls = emma_sample_size_ls + double_sample_size_ls
		real_kw_sample_size_ls = kw_sample_size_ls + double_sample_size_ls
		real_emma_div_kw_sample_size_ls = real_emma_sample_size_ls/real_kw_sample_size_ls
		emma_div_kw_sample_size_ls = emma_sample_size_ls/kw_sample_size_ls
		emma_div_kw_enrichment_ratio_ls = numpy.array(emma_enrich_ratio_ls)/numpy.array(kw_enrich_ratio_ls)
		
		#output a table for xy plot investigation
		import csv
		writer = csv.writer(open('%s.tsv'%output_fname_prefix, 'w'), delimiter='\t')
		header = ['phenotype_name', 'phenotype_name', \
				'double_sample_size', 'emma_only_sample_size', 'kw_only_sample_size', \
				'emma_div_kw_sample_size', \
				'real_emma_sample_size', 'real_kw_sample_size', 'real_emma_div_kw_sample_size', \
				'double_enrich_ratio', 'emma_enrich_ratio', 'kw_enrich_ratio', 'emma_div_kw_enrichment_ratio']
		writer.writerow(header)
		for i in range(len(xlabel_ls)):
			row = [xlabel_ls[i], xlabel_ls[i], double_sample_size_ls[i], emma_sample_size_ls[i], kw_sample_size_ls[i], emma_div_kw_sample_size_ls[i],\
				real_emma_sample_size_ls[i], real_kw_sample_size_ls[i], real_emma_div_kw_sample_size_ls[i], 
				double_enrich_ratio_ls[i],\
				emma_enrich_ratio_ls[i], kw_enrich_ratio_ls[i], emma_div_kw_enrichment_ratio_ls[i]]
			writer.writerow(row)
		del writer
		import pylab
		pylab.clf()
		n1 = pylab.hist(kw_sample_size_ls, 10, alpha=0.4)
		n2 = pylab.hist(emma_sample_size_ls, 10, alpha=0.4, facecolor='r')
		pylab.legend([n1[2][0], n2[2][0]], ['#SNPs above %s in KW'%r2_threshold, '#SNPs above %s in Emma'%r1_threshold])
		pylab.xlabel("#SNPs above a threshold")
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		#pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		sys.stderr.write("Done.\n")
	
	"""
biology_category_id = 1
output_fname_prefix = '/tmp/NoOfSNPsAboveThreshold_Flowering'
GWA.cmpRankAtThresholdOfTwoAnalysisMethods(db_250k, output_fname_prefix, biology_category_id=biology_category_id)
	"""
	
	def drawRandomEffectAndResidual(association_output, chromosome=None, position=None):
		"""
		2009-3-24
			association_output is special output by Association.py with y-xb trailling in each row
			get (u+\epsilon) = y-xb, in order to check its distribution
		"""
		import csv, pylab
		reader = csv.reader(open(association_output), delimiter='\t')
		for row in reader:
			chr = int(row[0])
			pos = int(row[1])
			if chromosome is not None and position is not None and (chr!=chromosome or pos!=position):
				continue
			random_effect_residual_ls = row[8:]
			random_effect_residual_ls = map(float, random_effect_residual_ls)
			pylab.clf()
			pylab.title('SNP %s %s mean=%s'%(row[0], row[1], pylab.mean(random_effect_residual_ls)))
			pylab.hist(random_effect_residual_ls, 20)
			pylab.show()
			to_continue = raw_input("Conitnue? (Y/n)")
			if to_continue=='n' or to_continue=='N':
				break
			_chromosome = raw_input("SNP chromosome: (%s)"%(chromosome))
			if _chromosome:
				chromosome = int(_chromosome)
			_position = raw_input("SNP position: (%s)"%(position))
			if _position:
				position = int(_position)
		del reader

	"""
association_output = os.path.expanduser('~/panfs/250k/association_results/call_method_29_y3_pheno_41_JIC4W.tsv')
drawRandomEffectAndResidual(association_output)
	"""
	
	@classmethod
	def contrastPvalueFromTwoGWA(cls, input_fname1, input_fname2, output_fname_prefix, list_type_id=132):
		"""
		2009-7-8
			contrast pvalue from two methods. each dot is a SNP. x-axis is one method. y-axis is the other.
			output the figure to output_fname_prefix
			
			- for input_fname1,  do_log10_transformation=False
			- for input_fname2,  do_log10_transformation=False
		"""
		
		from pymodule import PassingData
		from pymodule.SNP import getGenomeWideResultFromFile
		import Stock_250kDB
		from PlotCmpTwoAnalysisMethods import PlotCmpTwoAnalysisMethods
		from GeneListRankTest import GeneListRankTest, SnpsContextWrapper
		
		param_data = PassingData(need_the_value=1, construct_chr_pos2index=True)
		gwar1 = getGenomeWideResultFromFile(input_fname1, do_log10_transformation=False, pdata=param_data)
		gwar2 = getGenomeWideResultFromFile(input_fname2, do_log10_transformation=False, pdata=param_data)
		
		genome_wide_result_ls = [gwar1, gwar2]
		
		#snps_context_wrapper is not to be used,
		snps_context_picklef = '/Network/Data/250k/tmp-yh/snps_context/snps_context_g0_m5000'
		snps_context_wrapper = PlotCmpTwoAnalysisMethods.dealWithSnpsContextWrapper(snps_context_picklef, min_distance=5000, get_closest=0)
		candidate_gene_set = PlotCmpTwoAnalysisMethods.dealWithCandidateGeneList(list_type_id, return_set=True)
		pvalue_matching_data = PlotCmpTwoAnalysisMethods.matchPvaluesFromTwoResults(genome_wide_result_ls, snps_context_wrapper, \
																				candidate_gene_set)
		import pylab
		pylab.clf()
		"""
		pylab.subplots_adjust(left=0.08, right=0.92,bottom = 0.05)
		#calculate the number of rows needed according to how many score_rank_data, always two-column
		no_of_rows = 1
		no_of_cols = 1
		rm1, rm2 = rm_ls
		ax_scatter_pvalue = pylab.subplot(no_of_rows, no_of_cols, 1, frameon=False)
		"""
		ax = pylab.axes([0.1, 0.1, 0.8,0.8], frameon=False)
		rm1 = Stock_250kDB.ResultsMethod.get(3847)
		rm2 = Stock_250kDB.ResultsMethod.get(3847)
		rm_ls = [rm1, rm2]
		PlotCmpTwoAnalysisMethods.plot_scatter(ax, pvalue_matching_data, rm_ls, data_type=2)
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		
	"""
from GeneListRankTest import GeneListRankTest, SnpsContextWrapper	#otherwise,  snps_context_wrapper = cPickle.load(picklef) won't work
input_fname1 = '/Network/Data/250k/db/results/type_1/3847_results.tsv' #RF_Emwa1_32
input_fname2 = '/Network/Data/250k/tmp-atarone/Adnane new analyses/One node 20K trees/RF_nT20K_10_Emwa1.imp'
output_fname_prefix = '/tmp/RF_10k_vs_20k'
GWA.contrastPvalueFromTwoGWA(input_fname1, input_fname2, output_fname_prefix, list_type_id=132)
	"""
	
	@classmethod
	def cofactorLM(cls, genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, cofactors=[],\
					cofactor_phenotype_id_ls=[], report=1, run_type=1, genotype_cofactor_interaction=False, debug=False):
		"""
		2010-1-29
			add argument genotype_cofactor_interaction to test the interation between SNP and the first cofactor
		2010-1-11
			add argument cofactor_phenotype_id_ls to have the functionality of treating some phenotypes as cofactors
		2010-3-16
			one phenotype at a time:
				1. create a new SNP matrix which includes accessions whose phenotypes (this phenotype + cofactor_phenotype_id_ls) are non-NA
				2. one SNP at a time
					1. add cofactor SNP matrix if cofactors are present
					2. add cofactor phenotype matrix if cofactor_phenotype_id_ls exist
					3. run association
			
			parameter run_type is passed to Association.linear_model()
				run_type 1: pure linear model by python
				run_type 2: EMMA
				run_type 3: pure linear model by R (Field test shows run_type 3 is same as 1.)
				run_type 4: EMMAX
			
			start_snp and end_snp are on the same chromosome
		"""
		sys.stderr.write("Running association (pure linear model or EMMA) with cofactor ... \n")
		from Association import Association
		
		start_chr, start_pos = start_snp.split('_')[:2]
		start_chr = int(start_chr)
		start_pos = int(start_pos)
		
		stop_chr, stop_pos = stop_snp.split('_')[:2]
		stop_chr = int(stop_chr)
		stop_pos = int(stop_pos)
		
		
		test_type = 3 #Emma
		all_phenotype_id_ls = list(set(phenotype_method_id_ls+cofactor_phenotype_id_ls))	# used to remove accessions with NA-phenotype
		initData = Association.readInData(phenotype_fname, genotype_fname, eigen_vector_fname=None, \
										phenotype_method_id_ls=all_phenotype_id_ls, test_type=test_type)
		
		from PlotGroupOfSNPs import PlotGroupOfSNPs
		#which_phenotype_index_ls = initData.which_phenotype_ls
		which_phenotype_index_ls = PlotGroupOfSNPs.findOutWhichPhenotypeColumn(initData.phenData, set(phenotype_method_id_ls))
		
		environment_matrix = None
		data_matrix = initData.snpData.data_matrix
		min_data_point = 3
		
		# 2010-1-11 create a cofactor_phenotype_index_ls
		cofactor_phenotype_index_ls = PlotGroupOfSNPs.findOutWhichPhenotypeColumn(initData.phenData, set(cofactor_phenotype_id_ls))
		sys.stderr.write("%s phenotype cofactors.\n"%(len(cofactor_phenotype_id_ls)))
		
		#create an index list of cofactor SNPs
		cofactors_indices = []
		cofactors_set = set(cofactors)
		for i in range(len(initData.snpData.col_id_ls)):
			col_id = initData.snpData.col_id_ls[i]
			if col_id in cofactors_set:
				cofactors_indices.append(i)
		sys.stderr.write("%s SNP cofactors.\n"%(len(cofactors_indices)))
		
		import numpy, rpy
		rpy.r.source(os.path.expanduser('~/script/variation/src/gwa/emma/R/emma.R'))
		for which_phenotype in which_phenotype_index_ls:
			phenotype_name = initData.phenData.col_id_ls[which_phenotype]
			phenotype_name = phenotype_name.replace('/', '_')	#'/' will be recognized as directory in output_fname
			output_fname='%s_pheno_%s.tsv'%(os.path.splitext(output_fname_prefix)[0], phenotype_name)	#make up a new name corresponding to this phenotype
			
			#create non NA phenotype
			phenotype_ls = initData.phenData.data_matrix[:, which_phenotype]
			non_phenotype_NA_row_index_ls = []
			non_NA_phenotype_ls = []
			for i in range(len(phenotype_ls)):
				# 2010-1-11 make sure no NA in the cofactor phenotype matrix
				this_row_has_NA_phenotype = False
				if cofactor_phenotype_index_ls:
					for phenotype_index in cofactor_phenotype_index_ls:
						if numpy.isnan(initData.phenData.data_matrix[i, phenotype_index]):
							this_row_has_NA_phenotype = True
							break
				
				if numpy.isnan(phenotype_ls[i]):
					this_row_has_NA_phenotype = True
				if not this_row_has_NA_phenotype:
					non_phenotype_NA_row_index_ls.append(i)
					non_NA_phenotype_ls.append(phenotype_ls[i])
			
			non_NA_phenotype_ar = numpy.array(non_NA_phenotype_ls)
			#new_data_matrix = data_matrix[non_phenotype_NA_row_index_ls,:]
			newSNPData = SNPData.keepRowsByRowIndex(initData.snpData, non_NA_phenotype_row_index_ls)
			new_data_matrix = newSNPData.data_matrix
			if run_type==2:
				kinship_matrix = newSNPData.get_kinship_matrix()
				eig_L = rpy.r.emma_eigen_L(None, kinship_matrix)	#to avoid repeating the computation of eig_L inside emma.REMLE
			elif run_type==4:
				kinship_matrix = newSNPData.get_kinship_matrix()
				preEMMAX_data = Association.preEMMAX(non_NA_phenotype_ls, kinship_matrix, debug=True)
				variance_matrix = preEMMAX_data.variance_matrix
				non_NA_phenotype_ar = preEMMAX_data.non_NA_phenotype_ar
				L_inverse = preEMMAX_data.L_inverse
				non_NA_phenotype_ls = non_NA_phenotype_ar
			else:
				kinship_matrix = None
				eig_L = None
			
			no_of_rows, no_of_cols = new_data_matrix.shape
			results = []
			counter = 0
			real_counter = 0
			
			#create the cofactor matrix based on new SNP data matrix
			if len(cofactor_phenotype_index_ls)>0:
				cofactor_phenotype_matrix = initData.phenData.data_matrix[non_phenotype_NA_row_index_ls,:][:,cofactor_phenotype_index_ls]
				# 2010-1-29 
				# initData.phenData.data_matrix[non_phenotype_NA_row_index_ls,cofactor_phenotype_index_ls] is 1D if cofactor_phenotype_index_ls is of length 1.
				# while [non_phenotype_NA_row_index_ls,:][:,cofactor_phenotype_index_ls] produces a 2D matrix.
			else:
				cofactor_phenotype_matrix = None
			
			#create the cofactor matrix based on new SNP data matrix
			if len(cofactors_indices)>0:
				cofactor_data_matrix = new_data_matrix[:, cofactors_indices]
			else:
				cofactor_data_matrix = None
			
			#do association
			for j in range(no_of_cols):
				col_id = initData.snpData.col_id_ls[j]
				if col_id not in cofactors_set:	#same SNP appearing twice would cause singular design matrix
					chr, pos = col_id.split('_')[:2]
					chr = int(chr)
					pos = int(pos)
					if chr>=start_chr and chr<=stop_chr and pos>=start_pos and pos<=stop_pos:
						genotype_matrix = new_data_matrix[:,j]
						if cofactor_data_matrix is not None:
							if len(genotype_matrix.shape)<2:
								genotype_matrix = genotype_matrix.reshape([len(genotype_matrix), 1])
							genotype_matrix = numpy.hstack((genotype_matrix, cofactor_data_matrix))
						if cofactor_phenotype_matrix is not None:
							if len(genotype_matrix.shape)<2:
								genotype_matrix = genotype_matrix.reshape([len(genotype_matrix), 1])
							genotype_matrix = numpy.hstack((genotype_matrix, cofactor_phenotype_matrix))
						if genotype_cofactor_interaction and len(genotype_matrix.shape)==2:	# 2010-1-29 genotype_matrix has to be 2D. 
							interaction_vector = genotype_matrix[:,0]*genotype_matrix[:,1]
							interaction_vector = interaction_vector.reshape([len(interaction_vector), 1])
							genotype_matrix = numpy.hstack((genotype_matrix, interaction_vector))
						if run_type==4:
							pdata = Association.linear_model(genotype_matrix, non_NA_phenotype_ls, min_data_point, snp_index=j, \
									kinship_matrix=None, eig_L=None, run_type=4, counting_and_NA_checking=False,\
									variance_matrix=variance_matrix, lower_triangular_cholesky_inverse=L_inverse)
						else:
							pdata = Association.linear_model(genotype_matrix, non_NA_phenotype_ls, min_data_point, snp_index=j, \
												kinship_matrix=kinship_matrix, eig_L=eig_L, run_type=run_type)
						if pdata is not None:
							results.append(pdata)
							real_counter += 1
				
				counter += 1
				if report and counter%2000==0:
					sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
			if report:
				sys.stderr.write("%s\t%s\t%s\n"%('\x08'*40, counter, real_counter))
			
			#output
			Association.output_lm_results(results, initData.snpData.col_id_ls, output_fname)
		
		sys.stderr.write("Done.\n")
	
	"""
	genotype_fname = '/Network/Data/250k/tmp-yh/250k_data/call_method_17_test.tsv'
	phenotype_fname = '/Network/Data/250k/tmp-yh/phenotype.tsv'
	phenotype_method_id_ls = [285]
	output_fname_prefix = '/tmp/cofactorLM.tsv'
	start_snp = '1_2005921'
	stop_snp = '1_20054408'
	cofactors = ['1_3832974']
	cofactor_phenotype_id_ls = [77]
	GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
		cofactors=cofactors, cofactor_phenotype_id_ls=cofactor_phenotype_id_ls)
	
	# FLC expression association with two FRI deletion alleles as cofactor
	genotype_fname = '/Network/Data/250k/db/dataset/call_method_32.tsv'
	phenotype_fname = '/Network/Data/250k/tmp-yh/phenotype.tsv'
	phenotype_fname = '/tmp/phenotype_43_FLC.tsv'
	phenotype_method_id_ls = [43]
	start_snp = '4_1'
	stop_snp = '4_1750000'
					
	run_type = 2	# EMMA
	cofactors = ['4_268809']
	output_fname_prefix = '/tmp/Runtype_%s_SNP_%s_%s_Cofactor_%s.tsv'%(run_type, start_snp, stop_snp, cofactors[0])
	GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
					cofactors=cofactors, run_type=run_type)
	
	cofactors = ['4_269962']
	output_fname_prefix = '/tmp/Runtype_%s_SNP_%s_%s_Cofactor_%s.tsv'%(run_type, start_snp, stop_snp, cofactors[0])
	GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
					cofactors=cofactors, run_type=run_type)
	
	cofactors = ['4_268809', '4_269962']
	output_fname_prefix = '/tmp/Runtype_%s_SNP_%s_%s_Cofactor_%s_%s.tsv'%(run_type, start_snp, stop_snp, cofactors[0], cofactors[1])
	GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
					cofactors=cofactors, run_type=run_type)
	
	
	run_type = 1
	cofactors = ['4_268809']
	output_fname_prefix = '/tmp/Runtype_%s_SNP_%s_%s_Cofactor_%s.tsv'%(run_type, start_snp, stop_snp, cofactors[0])
	GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
					cofactors=cofactors, run_type=run_type)
	
	cofactors = ['4_269962']
	output_fname_prefix = '/tmp/Runtype_%s_SNP_%s_%s_Cofactor_%s.tsv'%(run_type, start_snp, stop_snp, cofactors[0])
	GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
					cofactors=cofactors, run_type=run_type)
	
	cofactors = ['4_268809', '4_269962']
	output_fname_prefix = '/tmp/Runtype_%s_SNP_%s_%s_Cofactor_%s_%s.tsv'%(run_type, start_snp, stop_snp, cofactors[0], cofactors[1])
	GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
					cofactors=cofactors, run_type=run_type)
	
	run_type =3
	output_fname_prefix = '/tmp/Runtype_%s_SNP_%s_%s.tsv'%(run_type, start_snp, stop_snp)
	GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
					cofactors=[], run_type=run_type)
	
	genotype_fname = '/Network/Data/250k/tmp-yh/250k_data/call_method_17_test.tsv'
	phenotype_fname = '/Network/Data/250k/tmp-yh/phenotype.tsv'
	phenotype_method_id_ls = [285]
	output_fname_prefix = '/tmp/cofactorLM.tsv'
	start_snp = '1_2005921'
	stop_snp = '1_20054408'
	cofactors = ['1_3832974']
	
	# 2010-1-25 G. oronti with LES as cofactor
	genotype_fname = '/Network/Data/250k/db/dataset/call_method_32.tsv'
	phenotype_fname = '/Network/Data/250k/tmp-yh/phenotype.tsv'
	phenotype_method_id_ls = [285]
	start_snp = '1_1'
	stop_snp = '5_60000000'
	cofactor_phenotype_id_ls = [77]
	for run_type in [1,2]:
		output_fname_prefix = os.path.expanduser('~/Runtype_%s_Cofactor_Phenotype_%s.tsv'%(run_type, cofactor_phenotype_id_ls[0]))
		GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
			cofactors=[], cofactor_phenotype_id_ls=cofactor_phenotype_id_ls, run_type=run_type)
	for run_type in [1,2]:
		output_fname_prefix = os.path.expanduser('~/Runtype_%s_Cofactor_Phenotype_%s_with_interaction.tsv'%(run_type, cofactor_phenotype_id_ls[0]))
		GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
			cofactors=[], cofactor_phenotype_id_ls=cofactor_phenotype_id_ls, run_type=run_type, genotype_cofactor_interaction=True)
	
	# 2010-3-16 FT22 with FLC and/or FRI as cofactor
	genotype_fname = '/Network/Data/250k/db/dataset/call_method_32.tsv'
	phenotype_fname = '/Network/Data/250k/tmp-yh/phenotype/phenotype.tsv'
	phenotype_method_id_ls = [7]
	start_snp = '1_1'
	stop_snp = '5_60000000'
	for cofactor_phenotype_id in [43, 44]:
		cofactor_phenotype_id_ls = [cofactor_phenotype_id]
		for run_type in [4, 1]:
			output_fname_prefix = os.path.expanduser('~/tmp/Runtype_%s_Cofactor_Phenotype_%s.tsv'%(run_type, cofactor_phenotype_id_ls[0]))
			GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
				cofactors=[], cofactor_phenotype_id_ls=cofactor_phenotype_id_ls, run_type=run_type)
	
	cofactor_phenotype_id_ls = [43, 44]	# FT22 with FLC & FRI as cofactor
	for run_type in [4, 1]:
		output_fname_prefix = os.path.expanduser('~/tmp/Runtype_%s_Cofactor_Phenotype_%s_%s.tsv'%(run_type, cofactor_phenotype_id_ls[0], cofactor_phenotype_id_ls[1]))
		GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
			cofactors=[], cofactor_phenotype_id_ls=cofactor_phenotype_id_ls, run_type=run_type)
	
	
	# genome-wide association with two FRI deletion alleles as cofactor
	call_method_id = 32
	genotype_fname = '/Network/Data/250k/db/dataset/call_method_%s.tsv'%call_method_id
	phenotype_fname = '/Network/Data/250k/tmp-yh/phenotype/phenotype.tsv'
	phenotype_method_id_ls = [43, 44]	# 5, 6, 7, 43, 44
	start_snp = '1_1'
	stop_snp = '5_60000000'
					
	for run_type in [4, 1]:
		cofactors = ['4_268809']
		output_fname_prefix = os.path.expanduser('~/tmp/Call_%s_Runtype_%s_SNP_%s_%s_Cofactor_%s.tsv'%\
												(call_method_id, run_type, start_snp, stop_snp, cofactors[0]))
		GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
					cofactors=cofactors, run_type=run_type)
	for run_type in [4, 1]:
		cofactors = ['4_269962']
		output_fname_prefix = os.path.expanduser('~/tmp/Call_%s_Runtype_%s_SNP_%s_%s_Cofactor_%s.tsv'%\
												(call_method_id, run_type, start_snp, stop_snp, cofactors[0]))
		GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
						cofactors=cofactors, run_type=run_type)
	
	for run_type in [4, 1]:
		cofactors = ['4_268809', '4_269962']
		output_fname_prefix = os.path.expanduser('~/tmp/Call_%s_Runtype_%s_SNP_%s_%s_Cofactor_%s_%s.tsv'%\
												(call_method_id, run_type, start_snp, stop_snp, cofactors[0], cofactors[1]))
		GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
					cofactors=cofactors, run_type=run_type)
	
	# 2010-3-16 FT22 with FLC and/or FRI as cofactor
	genotype_fname = '/Network/Data/250k/db/dataset/call_method_32.tsv'
	phenotype_fname = '/Network/Data/250k/tmp-yh/phenotype/phenotype.tsv'
	phenotype_method_id_ls = [7]
	start_snp = '1_1'
	stop_snp = '5_60000000'
	for cofactor_phenotype_id in [43, 44]:
		cofactor_phenotype_id_ls = [cofactor_phenotype_id]
		for run_type in [4, 1]:
			output_fname_prefix = os.path.expanduser('~/tmp/Runtype_%s_Cofactor_Phenotype_%s.tsv'%(run_type, cofactor_phenotype_id_ls[0]))
			GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
				cofactors=[], cofactor_phenotype_id_ls=cofactor_phenotype_id_ls, run_type=run_type)
	
	cofactor_phenotype_id_ls = [43, 44]	# FT22 with FLC & FRI as cofactor
	for run_type in [4, 1]:
		output_fname_prefix = os.path.expanduser('~/tmp/Runtype_%s_Cofactor_Phenotype_%s_%s.tsv'%(run_type, cofactor_phenotype_id_ls[0], cofactor_phenotype_id_ls[1]))
		GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
			cofactors=[], cofactor_phenotype_id_ls=cofactor_phenotype_id_ls, run_type=run_type)
	
	"""
	@classmethod
	def investigateBjarniEMMAX(cls, genotype_fname, phenotype_fname, output_fname, kinship_fname, phenotypeColName='ft_spspring', \
					cofactors=[], snp_id_to_be_included_ls=[],\
					run_type=1, genotype_cofactor_interaction=False, run_genome_scan=True):
		"""
		2010-9-6
			try to see why Bjarni's EMMAX and mine are different.
			
			test with justin's phenotype data, kinship matrix, because kinship is already calculated.
			
			argument run_type:
				1: pure_linear_model
				2: emma
				3: pure_linear_model via R
				4: generalized least square with specified variance matrix
				
				5: bjarni's emmax
			
			1. read in the genotype file
			2. read in the phenotype file
			3. read in the kinship matrix.
			4. remove accessions in phenotype file which are not in kinship matirx file
			5. create a incidence matrix mapping each accession in phenotype to the kinship matrix
			6. extract genotype data contained in genotype file which are phenotyped
		"""
		from pymodule.SNP import SNPData
		from Association import Association
		import numpy, math
		from pymodule.utils import addExtraToFilenamePrefix
		
		if os.path.isfile(output_fname):
			sys.stderr.write("Error: File %s already exits.\n"%output_fname)
			return None
		
		JBData = SNPData(input_fname=phenotype_fname, turn_into_array=1, ignore_2nd_column=1, \
							data_starting_col=2, turn_into_integer=False)
		
		if kinship_fname:
			kinship_output_fname = kinship_fname
		else:
			kinship_output_fname = os.path.splitext(genotype_fname)[0]+'_kinship.tsv'
		
		
		Z = None	# to check whether Z would be computed or not
		
		kinshipData, JBData, Z = JBDataGWA.getExpandedKinship(kinship_output_fname, genotype_fname, JBData)
		phenotype_ls = JBDataGWA.getPhenotypeLsOutOfJBData(JBData, logPhenotype=False, phenotypeColName=phenotypeColName)
		
		#2010-3-31 create a snpData with only selected columns
		selectSNPData = SNPData(input_fname=genotype_fname, turn_into_array=1, ignore_2nd_column=1, \
						col_id_key_set=set(snp_id_to_be_included_ls))
		selectSNPData, allele_index2allele_ls = selectSNPData.convert2Binary(row_id_as_major_allele="6909")	# Col-0 as the 
		# expand the selectSNPData to include replicates
		if Z is None:
			JBData = JBData.removeRowsNotInTargetSNPData(selectSNPData)
			Z = JBDataGWA.createIndividualToLineIncidenceMatrix(JBData.row_id_ls, selectSNPData.row_id_ls)
		snp_matrix_with_replicates = numpy.dot(Z, selectSNPData.data_matrix)
		selectSNPDataWithReplicates = SNPData(row_id_ls=JBData.row_id_ls, col_id_ls=selectSNPData.col_id_ls, \
											data_matrix=snp_matrix_with_replicates)
		L_inverse = None
		non_NA_phenotype_ar = numpy.array(phenotype_ls)
		if run_type==5:
			import linear_models
			"""
			# 2010-9-7 check lls/dlls vs logdelta
			lmm = linear_models.LinearMixedModel(phenotype_ls)
			lmm.add_random_effect(kinship_matrix)
			K = lmm.random_effects[1][1]
			eig_L = lmm._get_eigen_L_(K)
			res = lmm.get_expedited_REMLE(eig_L=eig_L) #Get the variance estimates..
			log_deltas = res['log_deltas']
			lls = res['lls']
			dlls = res['dlls']
			import pylab
			pylab.clf()
			pylab.plot(log_deltas, lls)
			pylab.show()
			import pdb
			pdb.set_trace()
			
			pylab.clf()
			pylab.plot(log_deltas, dlls)
			pylab.show()
			"""
			res = linear_models.emmax(selectSNPDataWithReplicates.data_matrix.T, phenotype_ls, kinshipData.data_matrix, cofactors=None, \
							with_interactions=False,int_af_threshold=15)
			import csv
			writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
			writer.writerow(['snp_id', 'pvalue', 'beta'])
			for i , col_id in enumerate(selectSNPDataWithReplicates.col_id_ls):
				data_row = (col_id, res["ps"][i], res["betas"][i])
				writer.writerow(data_row)
			sys.stderr.write("Exit after run type 5.\n")
			return
		elif run_type==4:
			non_NA_genotype_ls = []
			preEMMAX_data = Association.preEMMAX(phenotype_ls, kinshipData.data_matrix, non_NA_genotype_ls=None, debug=True)
			variance_matrix = preEMMAX_data.variance_matrix
			non_NA_phenotype_ar = preEMMAX_data.non_NA_phenotype_ar
			L_inverse = preEMMAX_data.L_inverse
			
			#2010-9-7 plot the likelihood  and d(likelihood) curve
			one_emma_rs = preEMMAX_data.one_emma_rs
			logdelta = one_emma_rs.logdelta
			LL = one_emma_rs.LL
			dLL = one_emma_rs.dLL
			import pylab
			pylab.clf()
			pylab.plot(logdelta, LL)
			pylab.show()
			import pdb
			pdb.set_trace()
			pylab.clf()
			pylab.plot(logdelta, dLL)
			pylab.show()
			
		#else:
		#	sys.stderr.write("Exit.\n")
		#	return
		
		phenotype_variance = numpy.var(phenotype_ls)	# 2010-4-18 to calculate the variance explained without cholesky transformation
		#non_NA_phenotype_ar = numpy.dot(L_inverse, non_NA_phenotype_ar)	# numpy.dot and numpy.inner has subtle difference.
		
		extraVariateNameLs = ['S_square', 'var_perc','var_perc_real', 'vg', 've', 'delta', 'heritability',]	#2010-4-18 explicitly set the extra columns to be outputted
		writer = JBDataGWA.writeHeaderJBData(output_fname, base_formula=[], \
							GXE_environment_variate_name_ls=[],\
							special_interaction_snp_id_ls=[],\
							run_genome_scan=True, extraVariateNameLs=extraVariateNameLs)
		counter = 0
		real_counter = 0
		if run_genome_scan==1 or run_genome_scan==True:
			snpData = selectSNPDataWithReplicates
			Z = JBDataGWA.createIndividualToLineIncidenceMatrix(JBData.row_id_ls, snpData.row_id_ls)
			results = []
			no_of_cols = snpData.data_matrix.shape[1]
			snp_id_ls = snpData.col_id_ls
			for snp_id in snp_id_ls:
				j = snpData.col_id2col_index[snp_id]
				genotype_ls = snpData.data_matrix[:,j]
				genotype_matrix = genotype_ls.reshape([len(genotype_ls), 1])
				genotype_matrix = numpy.dot(Z, genotype_matrix)
				snp_id = snpData.col_id_ls[j]
				
				pdata = Association.linear_model(genotype_matrix, non_NA_phenotype_ar, min_data_point=3, \
												snp_index=snp_id, \
									kinship_matrix=kinshipData.data_matrix, eig_L=None, run_type=run_type, \
									counting_and_NA_checking=True,\
									variance_matrix=None, lower_triangular_cholesky_inverse=L_inverse)
				
				# counting_and_NA_checking=True to get MAF and MAC. doesn't matter if it's False.
				if pdata is not None:
					results.append(pdata)
					Association.output_multi_variate_lm_results([pdata], writer, run_genome_scan=run_genome_scan,\
															extraVariateNameLs=extraVariateNameLs)
					real_counter += 1
				counter += 1
				if counter%2000==0:
					sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
			sys.stderr.write("%s\t%s\t%s\n"%('\x08'*40, counter, real_counter))
			return
	"""
		#2010-9-6
		phenotype_genotype_fname = os.path.expanduser('~/Downloads/BLUP_381_spSpring.csv')
		
		genotype_fname_to_generate_kinship = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering20100820/call_method_49_core482_with_FRI_del_chr_order_one_time_impute_yu_format.tsv')
		kinship_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering20100820/call_method_49_core482_kinship.tsv')	
		interaction_snp_id_in_base_formula_ls = []
		snp_id_to_be_included_ls=['1_24345319', '1_3978063', '2_8516520', '3_9340928', \
								'4_1356197',  '4_158958', '4_268809', '4_269962', '4_387727', \
								'5_18620282', '5_25376551', '5_3188328']
		
		
		run_type = 5
		output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering20100820/investigateBjarniEMMAX_DTF16replicates_%s.tsv'%\
			(run_type))
		GWA.investigateBjarniEMMAX(genotype_fname_to_generate_kinship, phenotype_genotype_fname, \
								output_fname, kinship_fname,\
								snp_id_to_be_included_ls = snp_id_to_be_included_ls,\
								run_type=run_type)
		sys.exit(0)
		
		#2010-9-6
		phenotype_genotype_fname = os.path.expanduser('~/Downloads/BLUP_381_spSpring.csv')
		
		genotype_fname_to_generate_kinship = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering20100820/call_method_49_core482_with_FRI_del_chr_order_one_time_impute_yu_format.tsv')
		kinship_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering20100820/call_method_49_core482_kinship.tsv')	
		interaction_snp_id_in_base_formula_ls = []
		snp_id_to_be_included_ls=['1_24345319', '1_3978063', '2_8516520', '3_9340928', \
								'4_1356197',  '4_158958', '4_268809', '4_269962', '4_387727', \
								'5_18620282', '5_25376551', '5_3188328']
		
		
		run_type = 4
		output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering20100820/investigateBjarniEMMAX_DTF16replicates_%s.tsv'%\
			(run_type))
		GWA.investigateBjarniEMMAX(genotype_fname_to_generate_kinship, phenotype_genotype_fname, \
								output_fname, kinship_fname,\
								snp_id_to_be_included_ls = snp_id_to_be_included_ls,\
								run_type=run_type)
		sys.exit(0)
	"""
	@classmethod
	def linkEcotypeID2TargetEcotypeID(cls, db, JBd4DTF_fname, output_fname):
		"""
		2010-3-31
		
		"""
		sys.stderr.write("Link the ecotype id in d4DTF (JBLab flowering data) to target ecotypeid ... \n")
		from pymodule.SNP import SNPData
		JBData = SNPData(input_fname=JBd4DTF_fname, turn_into_array=1, ignore_2nd_column=1, \
						    data_starting_col=2, turn_into_integer=False)
		
		from common import get_ecotypeid2tg_ecotypeid
		ecotypeid2tg_ecotypeid = get_ecotypeid2tg_ecotypeid(db.metadata.bind)
		
		target_ecotype_id_ls = []
		no_of_total = len(JBData.row_id_ls)
		no_of_missed = 0
		for row_id in JBData.row_id_ls:
			ecotype_id = int(row_id)
			tg_ecotype_id =  ecotypeid2tg_ecotypeid.get(ecotype_id)
			if tg_ecotype_id:
				target_ecotype_id_ls.append(tg_ecotype_id)
			else:
				no_of_missed += 1
				target_ecotype_id_ls.append('NA')
		sys.stderr.write("%s out of %s rows miss target ecotype ids.\n"%(no_of_missed, no_of_total))
		
		JBData.row_id_ls = target_ecotype_id_ls
		JBData.strain_acc_list = target_ecotype_id_ls
		JBData.tofile(output_fname)
	
	"""
	# 2010-3-31 original JBLab's data compiled by Yan but with only 469 ecotypes
	JBd4DTF_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF.tsv')
	output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/d4DTF_tg_ecotypeid.tsv')
	GWA.linkEcotypeID2TargetEcotypeID(db_250k, JBd4DTF_fname, output_fname)
	
	# 2010-3-31 JBLab's phenotype data compiled by Yan with all ecotypes
	JBd4DTF_fname = os.path.expanduser('/tmp/DaysToFlower16replicates.csv')	# from google docs
	output_fname = os.path.expanduser('~/script/variation/data/JBLabSeasonFlowering/DaysToFlower16replicates_tg_ecotypeid.tsv')
	GWA.linkEcotypeID2TargetEcotypeID(db_250k, JBd4DTF_fname, output_fname)
	
	"""
	
	@classmethod
	def phenotypePCA(cls, db, phenotype_fname, output_fname_prefix):
		"""
		2009-9-1
			Run PCA on the phenotype matrix either ecotype-wise or phenotype-wise.
			
			output the data in a matrix fashion that the web MotionChartAppMCPanel app would recognize 
		"""
		from pymodule import read_data, SNPData
		import numpy, csv
		header_phen, strain_acc_list_phen, category_list_phen, data_matrix_phen = read_data(phenotype_fname, turn_into_integer=2, matrix_data_type=float)
		data_matrix_phen = numpy.array(data_matrix_phen)
		phenData = SNPData(header=header_phen, strain_acc_list=strain_acc_list_phen, category_list=category_list_phen, \
						data_matrix=data_matrix_phen)
		
		
		min_data_point = 4
		sys.stderr.write("Removing ecotypes with too few phenotype points (<%s) ..."%min_data_point)
		import pca_module
		rows_to_be_kept = []
		row_labels = []
		for i in range(len(phenData.row_id_ls)):
			row_id = phenData.row_id_ls[i]
			
			no_of_valid_points = sum(numpy.isfinite(phenData.data_matrix[i,:]))
			if no_of_valid_points>=min_data_point:
				rows_to_be_kept.append(row_id)
				row_labels.append('%s %s'%(row_id[1], row_id[0]))
		filteredPhenData = SNPData.keepRowsByRowID(phenData, rows_to_be_kept)
			#order in rows_to_be_kept is kept because it's in the same order as phenData.row_id_ls
		phenData = filteredPhenData
		sys.stderr.write("Done.\n")
		
		phenData_trans = SNPData(row_id_ls=phenData.col_id_ls, col_id_ls=phenData.row_id_ls, \
								data_matrix=numpy.transpose(phenData.data_matrix))
		
		sys.stderr.write("Carrying out phenotype-wise PCA ...")
		# phenotype-wise PCA
		import Stock_250kDB 
		phenotypePCA_fname = '%s_phenotype.tsv'%output_fname_prefix
		phenotypePCA_writer = csv.writer(open(phenotypePCA_fname, 'w'), delimiter='\t')
		
		import pca_module
		from pymodule.PCA import PCA
		#T, P, explained_var = pca_module.PCA_svd(phenData_trans.data_matrix, standardize=True)
		T, P, explained_var = PCA.eig(phenData_trans.data_matrix, normalize=False)	#normalize=True causes missing value in the covariance matrix
		# get the category information for each phenotype
		header = ['phenotype_label', 'category|string', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6']
		phenotypePCA_writer.writerow(header)
		for i in range(len(phenData_trans.row_id_ls)):
			row_id = phenData_trans.row_id_ls[i]
			row_id_tuple = row_id.split('_')
			phenotype_method_id=int(row_id_tuple[0])
			pm = Stock_250kDB.PhenotypeMethod.get(phenotype_method_id)
			if pm.biology_category:
				category = pm.biology_category.short_name
			else:
				category = 'other'
			data_row = [row_id, category] + list(T[i,0:6])
			phenotypePCA_writer.writerow(data_row)
		del phenotypePCA_writer
		sys.stderr.write("Done.\n")
		
		sys.stderr.write("Carrying out ecotype-wise PCA ... \n")
		# ecotype-wise PCA
		
		# run PCA
		#T, P, explained_var = pca_module.PCA_svd(phenData.data_matrix, standardize=True)	#SVD doesn't converge
		T, P, explained_var = PCA.eig(phenData.data_matrix, normalize=False)	#normalize=True gives rise to missing value in the covariance matrix
		from common import getEcotypeInfo
		ecotype_info = getEcotypeInfo(db)
		
		# output
		ecotypePCA_fname = '%s_ecotype.tsv'%output_fname_prefix
		ecotypePCA_writer = csv.writer(open(ecotypePCA_fname, 'w'), delimiter='\t')
		header = ['ecotype_label', 'region|string', 'country|string', 'lat', 'lon', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6']
		ecotypePCA_writer.writerow(header)
		for i in range(len(row_labels)):
			row_label = row_labels[i]
			ecotype_id = int(rows_to_be_kept[i][0])
			ecotype_obj = ecotype_info.ecotype_id2ecotype_obj.get(ecotype_id)
			if ecotype_obj is not None:
				region = ecotype_obj.region
				country = ecotype_obj.country
				lat = ecotype_obj.latitude
				lon = ecotype_obj.longitude
			else:
				region, country, lat, long = None, None, None, None
			data_row = [row_label, region, country, lat, lon] + list(T[i,0:6])
			ecotypePCA_writer.writerow(data_row)
		del ecotypePCA_writer
		sys.stderr.write("Done.\n")
	
	"""
	phenotype_fname = os.path.expanduser('~/mnt/panfs/250k/phenotype20090902.tsv')
	output_fname_prefix = '/tmp/phenotypePCA'
	GWA.phenotypePCA(db_250k, phenotype_fname, output_fname_prefix)
	
	phenotype_fname = os.path.expanduser('~/mnt/panfs/250k/phenotype20090902_cypress_col_1_9.tsv')	#~/mnt/panfs/250k/phenotype20090902_cypress.tsv')
	output_fname_prefix = '/tmp/phenotypePCA_cypress_col_1_9'
	GWA.phenotypePCA(db_250k, phenotype_fname, output_fname_prefix)
	
	phenotype_fname = os.path.expanduser('~/mnt/panfs/250k/phenotype20090902.tsv')
	output_fname_prefix = '/tmp/phenotypePCA'
	GWA.phenotypePCA(db_250k, phenotype_fname, output_fname_prefix)
	"""
	
	@classmethod
	def plotEnrichmentRatioVsCutoff(cls, db_250k, output_fname_prefix, phenotype_method_id_ls=[9,10,11,12,13], \
								analysis_method_id_ls=[1], list_type_id=149, type_id=2, stat_on_x_axis_type=1,\
								call_method_id=32, min_distance=20000):
		"""
		2009-10-22
			Fetch enrichment data from candidate_gene_top_snp_test_rm, plot the enrichment ratio against the cutoff or no of top SNPs.
				for one analysis method, and one/multiple phenotypes
			
			stat_on_x_axis_type=1: min_score
			stat_on_x_axis_type=2: no_of_top_snps
			
			
			The alternative way to get data matrix thru raw sql:
			
			select r.phenotype_method_id, p.short_name, r.analysis_method_id, newt.* from
				(select c.results_id, c.list_type_id, c.pvalue, c.candidate_sample_size, c.candidate_gw_size, 
				c.non_candidate_sample_size, c.non_candidate_gw_size, c.no_of_top_snps, c.max_score, c.min_score 
				from candidate_gene_top_snp_test_rm c where c.type_id=2 and c.results_id in (select id from 
				results_method where call_method_id =32 and phenotype_method_id >=9 and phenotype_method_id <=13 
				and analysis_method_id in (1,6)) and c.list_type_id=149 and c.min_distance=20000) as newt , 
				results_method r, phenotype_method p where r.id=newt.results_id and r.phenotype_method_id=p.id 
				order by phenotype_method_id, analysis_method_id, no_of_top_snps;
		"""
		sys.stderr.write("Plotting enrichment ratio vs cutoff ...")
		from Stock_250kDB import CandidateGeneTopSNPTestRM, ResultsMethod
		import math
		rows = ResultsMethod.query.filter(ResultsMethod.analysis_method_id.in_(analysis_method_id_ls)).\
			filter_by(call_method_id=call_method_id).\
			filter(ResultsMethod.phenotype_method_id.in_(phenotype_method_id_ls))
		results_id_ls = []
		for row in rows:
			results_id_ls.append(row.id)
		
		rows = CandidateGeneTopSNPTestRM.query.filter_by(type_id=type_id).filter_by(list_type_id=list_type_id).\
			filter_by(min_distance=min_distance).\
			filter(CandidateGeneTopSNPTestRM.results_id.in_(results_id_ls)).order_by(CandidateGeneTopSNPTestRM.results_id).\
			order_by(CandidateGeneTopSNPTestRM.no_of_top_snps)
		phenotype_id2name = {}
		phenotype_id2xy_ls = {}
		for row in rows:
			phenotype = row.result.phenotype_method
			phenotype_id = phenotype.id
			if phenotype_id not in phenotype_id2name:
				phenotype_id2name[phenotype_id] = phenotype.short_name
				phenotype_id2xy_ls[phenotype_id] = [[], []]
			candidate_ratio = row.candidate_sample_size/float(row.candidate_gw_size)
			non_candidate_ratio = row.non_candidate_sample_size/float(row.non_candidate_gw_size)
			if non_candidate_ratio >0:
				enrichment_ratio = candidate_ratio/non_candidate_ratio
				if stat_on_x_axis_type==1:
					stat_on_x_axis = row.min_score
					analysis_method = row.result.analysis_method
					if analysis_method.smaller_score_more_significant==1:
						stat_on_x_axis = math.pow(10, -stat_on_x_axis)	# the score is -log(pvalue), need to convert it back
				elif stat_on_x_axis_type==2:
					stat_on_x_axis = row.no_of_top_snps
				phenotype_id2xy_ls[phenotype_id][0].append(stat_on_x_axis)
				phenotype_id2xy_ls[phenotype_id][1].append(enrichment_ratio)
		
		import pylab
		pylab.clf()
		phenotype_id_ls = phenotype_id2xy_ls.keys()
		phenotype_id_ls.sort()
		legend_ls = []
		for phenotype_id in phenotype_id_ls:
			x_ls, y_ls = phenotype_id2xy_ls[phenotype_id]
			phenotype_name = phenotype_id2name[phenotype_id]
			legend_ls.append(phenotype_name)
			pylab.semilogx(x_ls, y_ls)
			
			#pylab.loglog(x_ls, y_ls, basex=10)
		pylab.ylabel("Enrichment Ratio")
		pylab.legend(legend_ls)
		if output_fname_prefix:
			pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
			pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		sys.stderr.write("Done.\n")
	
	"""
	output_fname_prefix = '/tmp/diseaseResistance_Wilcoxon'
	GWA.plotEnrichmentRatioVsCutoff(db_250k, output_fname_prefix, phenotype_method_id_ls=[9,10,11,12,13], \
								analysis_method_id_ls=[1], stat_on_x_axis_type=1)
	
	output_fname_prefix = '/tmp/diseaseResistance_Wilcoxon_no_of_top_SNPs_on_X'
	GWA.plotEnrichmentRatioVsCutoff(db_250k, output_fname_prefix, phenotype_method_id_ls=[9,10,11,12,13], \
								analysis_method_id_ls=[1], stat_on_x_axis_type=2)
	
	output_fname_prefix = '/tmp/diseaseResistance_RF'
	GWA.plotEnrichmentRatioVsCutoff(db_250k, output_fname_prefix, phenotype_method_id_ls=[9,10,11,12,13], \
								analysis_method_id_ls=[6], stat_on_x_axis_type=2)
	
	"""
	
	@classmethod
	def predictByRpart(cls, data_matrix, phenotype_ls, col_id_ls=None, output_fname=None):
		"""
		2009-11-18
			call rpart (=randomForest) to predict phenotype_ls based on data_matrix
		"""
				
		import rpy
		from rpy import r
		r.library('rpart')
		rpart_cp = 0.01
		rpy.set_default_mode(rpy.NO_CONVERSION)
			
		pre_data_frame_dict = {}
		for i in range(len(col_id_ls)):
			col_id = col_id_ls[i]
			pre_data_frame_dict[col_id] = rpy.r.as_factor(data_matrix[:,i])	# all SNPs are treated as factor
		
		pre_data_frame_dict["phenotype"] = phenotype_ls
		data_frame = r.as_data_frame(pre_data_frame_dict)
		fit = r.rpart(r("phenotype~."), data=data_frame, method="anova", control=r.rpart_control(cp=rpart_cp))
			#,\
			#	parms=r.list(prior=prior_prob, loss=r.matrix(loss_matrix) ) )
		r.postscript('%s_rsq.ps'%output_fname)
		a=r.rsq_rpart(fit)
		r.dev_off()
		r.postscript('%s_tree.ps'%output_fname)
		a=r.plot(fit)
		r.text(fit, use_n=rpy.r.TRUE)
		r.dev_off()
		r.postscript('%s_cp.ps'%output_fname)
		a=r.plotcp(fit)
		r.dev_off()
	
	@classmethod
	def predictBySVM(cls, data_matrix, phenotype_ls, col_id_ls=None, output_fname=None):
		"""
		2009-11-18
			use SVM (support vector machine) to predict phenotype_ls based on data_matrix
			package python-libsvm is based on http://www.csie.ntu.edu.tw/~cjlin/libsvm/
		"""
		from svm import svm_problem, svm_parameter, svm_model, cross_validation, LINEAR, POLY, RBF
		import numpy
		problem = svm_problem(phenotype_ls, data_matrix)
		size = len(phenotype_ls)
		
		total_correct = 0.
		kernels2kname = {LINEAR:'linear', POLY:'polynomial', RBF:'rbf'}
		kernels2kname = {POLY:'polynomial'}	# 2009- 11-18 only polynomial
		import sys
		mute_device = open('/dev/null', 'w')
		
		param = svm_parameter(C = 10)	#,nr_weight = 2,weight_label = [1,0],weight = [10,1])
		C_ls = [100, 10, 1, 0.1, 0.01]	# tried 1000, 10, 1, 0.1 ,0.01
		#C_ls = [100]
		#gamma_ls = [100, 10, 1, 0.1, 0.0001, 0]
		gamma_ls = [0]	# test shows that gamma affects rbf quite a bit
		nr_fold = 10	# 10 fold cross validation
		for k, kname in kernels2kname.iteritems():
			for C in C_ls:
				for gamma in gamma_ls:
					param.gamma = gamma
					param.C = C
					param.kernel_type = k;
					sys.stderr = mute_device
					sys.stdout = mute_device
					target = cross_validation(problem, param, nr_fold)
					#model = svm_model(problem,param)
					sys.stderr = sys.__stderr__
					sys.stdout = sys.__stdout__
					errors = 0.
					for i in range(size):
							#prediction = model.predict(data_matrix[i])
							prediction = target[i]
							#probability = model.predict_probability
							errors += (prediction-phenotype_ls[i])*(prediction-phenotype_ls[i])
							"""
							if (phenotype_ls[i] != prediction):
									errors = errors + 1
							"""
					errors = errors/size
					print "##########################################"
					print " kernel %s, C %s, gamma %s: mse = %s" % (kname, C, gamma, errors)
					print "##########################################"
	
	@classmethod
	def predictPhenotypeBasedOnGenotype(cls, genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, cofactors=[],\
									report=1, run_type=1):
		"""
		2009-11-18
			predict phenotype based on genotype matrix using machine-learning approaches (rpart, svm)
		"""
		sys.stderr.write("Predicting phenotype based on genotype ... \n")
		from Association import Association
		eigen_vector_fname = ''
		
		test_type = 3 #Emma
		initData = Association.readInData(phenotype_fname, genotype_fname, eigen_vector_fname, phenotype_method_id_ls, test_type=test_type)
		
		which_phenotype_index_ls = initData.which_phenotype_ls
		which_phenotype_index_ls.sort()
		environment_matrix = None
		data_matrix = initData.snpData.data_matrix
		min_data_point = 3
		
		#create an index list of cofactor SNPs
		cofactors_indices = []
		cofactors_set = set(cofactors)
		for i in range(len(initData.snpData.col_id_ls)):
			col_id = initData.snpData.col_id_ls[i]
			if col_id in cofactors_set:
				cofactors_indices.append(i)
		import numpy
		for which_phenotype in which_phenotype_index_ls:
			phenotype_name = initData.phenData.col_id_ls[which_phenotype]
			phenotype_name = phenotype_name.replace('/', '_')	#'/' will be recognized as directory in output_fname
			print phenotype_name
			output_fname='%s_pheno_%s'%(os.path.splitext(output_fname_prefix)[0], phenotype_name)	#make up a new name corresponding to this phenotype
			
			#create non NA phenotype
			phenotype_ls = initData.phenData.data_matrix[:, which_phenotype]
			non_phenotype_NA_row_index_ls = []
			non_NA_phenotype_ls = []
			for i in range(len(phenotype_ls)):
				if not numpy.isnan(phenotype_ls[i]):
					non_phenotype_NA_row_index_ls.append(i)
					non_NA_phenotype_ls.append(phenotype_ls[i])
			
			non_NA_phenotype_ar = numpy.array(non_NA_phenotype_ls)
			new_data_matrix = data_matrix[non_phenotype_NA_row_index_ls,:]
			if run_type==1:
				#new_data_matrix = new_data_matrix.tolist()
				cls.predictByRpart(new_data_matrix, non_NA_phenotype_ls, col_id_ls=initData.snpData.col_id_ls, output_fname=output_fname)			
			elif run_type==2:
				cls.predictBySVM(new_data_matrix, non_NA_phenotype_ls, col_id_ls=initData.snpData.col_id_ls, output_fname=output_fname)
			else:
				sys.stderr.write("Run type %s not supported.\n"%run_type)
		sys.stderr.write("Done.\n")
	"""
	genotype_fname = "/Network/Data/250k/tmp-yh/250k_data/call_method_17_test.tsv"
	genotype_fname = "/Network/Data/250k/db/dataset/call_method_32.tsv"
	phenotype_fname = '/Network/Data/250k/tmp-yh/phenotype20091113.tsv'
	phenotype_method_id_ls = [1,2]
	phenotype_method_id_ls = range(1,8)
	output_fname_prefix = '/tmp/predictPhenotype'
	run_type =2 
	GWA.predictPhenotypeBasedOnGenotype(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, run_type=run_type)
	"""
	
	
class FileFormatExchange(object):
	
	@classmethod
	def combineMultiFastaIntoOne(cls, input_dir, numberChrOrder, nonNumberChrOrder, output_fname):
		"""
		2011-1-31
			recover a function that has been lost.
			
			Given a directory of chromosome sequences downloaded from NCBI, find out which chromosome is in which file,
				order them in a pre-specified way (number + X Y chl MT) that is same as IGV, and output in a 
				multi-record fasta file.
				
		"""
		
	
	def convertChiamoOutput(chiamo_infname, chiamo_outfname, ref_250k_infname, output_fname, posterior_min=0.95):
		"""
		2008-03-18
			check chiamo output
		"""
		from variation.src.common import nt2number
		import csv
		import numpy
		
		reader = csv.reader(open(ref_250k_infname), delimiter='\t')
		SNP_acc_ls = reader.next()[2:]
		del reader
		
		SNP_id2allele_ls = {}
		SNP_id2info = {}
		for SNP_acc in SNP_acc_ls:
			chr, pos, alleleA, alleleB = SNP_acc.split('_')
			SNP_id = '%s_%s'%(chr, pos)
			SNP_id2allele_ls[SNP_id] = [alleleA, alleleB]
			SNP_id2info[SNP_id] = [SNP_acc, len(SNP_id2allele_ls)-1]	#[name, index]
		
		reader = csv.reader(open(chiamo_infname), delimiter='\t')
		strain_id_dup_ls = reader.next()[5:]
		strain_id_ls = []
		for i in range(0, len(strain_id_dup_ls), 2):
			strain_id_ls.append(strain_id_dup_ls[i][:-2])
		del reader
		
		sys.stderr.write("Reading chiamo output ...\n")
		reader = csv.reader(open(chiamo_outfname), delimiter=' ')
		snp_acc_ls = []
		data_matrix = numpy.zeros([len(strain_id_ls), len(SNP_id2allele_ls)], numpy.integer)
		counter = 0
		for row in reader:
			SNP_id = row[0]
			if SNP_id not in SNP_id2allele_ls:
				continue
			allele_ls = SNP_id2allele_ls[SNP_id]
			snp_acc, col_index = SNP_id2info[SNP_id]
			snp_acc_ls.append(snp_acc)
			for i in range(5, len(row)-1, 3):
				posterior_ls = [float(row[i]), float(row[i+1]), float(row[i+2])]	#posterior for 3 classes
				max_posterior_index = numpy.argmax(posterior_ls)
				if posterior_ls[max_posterior_index]>=posterior_min:
					if max_posterior_index==2:	#heterozygous
						call = '%s%s'%(allele_ls[0], allele_ls[1])
					else:
						call = allele_ls[max_posterior_index]
					row_index = (i-5)/3
					try:
						data_matrix[row_index][col_index] = nt2number[call]
					except:
						print SNP_id, snp_acc, col_index, allele_ls, row_index, strain_id_ls[row_index], call
						return
			counter += 1
			if counter%2000==0:
				sys.stderr.write("%s%s"%('\x08'*40, counter))
					
		del reader
		sys.stderr.write("Done.\n")
		
		from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header = ['ecotypeid', 'ecotypeid'] + snp_acc_ls
		FilterStrainSNPMatrix_instance.write_data_matrix(data_matrix, output_fname, header, strain_id_ls, [1]*len(strain_id_ls))
		sys.stderr.write("Finished.\n")
	
	"""
	chiamo_infname = os.path.expanduser('~/script/affy/250k_test/yanli8-29-07_chiamo_14SNPs.in')
	chiamo_outfname = os.path.expanduser('~/script/affy/250k_test/yanli8-29-07_chiamo.out_0_mcmc')
	ref_250k_infname = os.path.expanduser('~/script/variation/genotyping/250ksnp/data/data_250k.tsv')
	output_fname = os.path.expanduser('~/script/affy/250k_test/yanli8-29-07_chiamo_out.tsv')
	
	convertChiamoOutput(chiamo_infname, chiamo_outfname, ref_250k_infname, output_fname, posterior_min=0.95)
	"""
	
	

	@classmethod
	def removeRowsBasedOnSNPAllele(input_fname, output_fname, SNP_label, allele='-'):
		"""
		2008-08-04 investigate whether throwing off some rows help to increase significance
		"""
		from pymodule import read_data, write_data_matrix, nt2number
		header, strain_acc_list, category_list, data_matrix = read_data(input_fname, turn_into_integer=1)
		snp_acc_ls = header[2:]
		snp_index = -1
		for i in range(len(snp_acc_ls)):
			if snp_acc_ls[i]==SNP_label:
				snp_index = i
		
		allele = nt2number[allele]
		rows_to_be_tossed_out = set()
		for i in range(len(data_matrix)):
			if data_matrix[i][snp_index]==allele:
				rows_to_be_tossed_out.add(i)
		
		write_data_matrix(data_matrix, output_fname, header, strain_acc_list, category_list, rows_to_be_tossed_out=rows_to_be_tossed_out)
	
	"""
input_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix.tsv')
output_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix_2.tsv')
SNP_label = '4_268809_0'
FileFormatExchange.removeRowsBasedOnSNPAllele(input_fname, output_fname, SNP_label, allele='-')
	"""
	
	@classmethod
	def removeSNPsWithMoreThan2Alleles(cls, input_fname, output_fname):
		"""
		2008-08-05
			NPUTE can't work with SNPs with >2 alleles
		"""
		from pymodule import SNPData
		snpData = SNPData(input_fname=input_fname, turn_into_integer=1, turn_into_array=1)
		newSNPData = snpData.removeSNPsWithMoreThan2Alleles(snpData)
		newSNPData.tofile(output_fname)
	
	"""
input_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix.tsv')
output_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix_only_2_alleles.tsv')
FileFormatExchange.removeSNPsWithMoreThan2Alleles(input_fname, output_fname)
	"""
	
	
	@classmethod
	def turnNPUTEOutputIntoYuFormat(cls, input_fname, output_fname):
		"""
		2008-08-05
			NPUTE output format is SNPXStrain by and large.
				1st and 2nd column are same as input's 1st row. 1st row is input's 1st column. 2nd row is input's 2nd column.
			
		"""
		from pymodule import SNPData
		snpData = SNPData(input_fname=input_fname, turn_into_integer=1, turn_into_array=1, double_header=1, ignore_2nd_column=1)
		snpData.col_id_ls = snpData.header[0][2:]	#take the first header
		
		from pymodule.SNP import transposeSNPData
		newSNPData = transposeSNPData(snpData)
		newSNPData.strain_acc_list = newSNPData.row_id_ls
		newSNPData.category_list = snpData.header[1][2:]
		newSNPData.tofile(output_fname)
	"""
input_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix.NPUTE.tsv')
output_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix.NPUTE.yh.tsv')
FileFormatExchange.turnNPUTEOutputIntoYuFormat(input_fname, output_fname)

input_fname = os.path.expanduser('~/script/variation/data/ACD6/snp_matrix.84_1530.impute.tsv')
output_fname = os.path.expanduser('~/script/variation/data/ACD6/snp_matrix.84_1530.impute.Yu.format.tsv')
FileFormatExchange.turnNPUTEOutputIntoYuFormat(input_fname, output_fname)
	"""
	

	@classmethod
	def outputSNPmatrixGivenRegion(cls, snpData, output_fname, chr, start_pos, stop_pos):
		"""
		2008-12-15 output a chromosome region of the SNP matrix
		"""
		sys.stderr.write("Outputting a selected region of SNP matrix ...")
		from pymodule import read_data, write_data_matrix, nt2number
		#header, strain_acc_list, category_list, data_matrix = read_data(input_fname, turn_into_integer=1)
		snp_acc_ls = snpData.col_id_ls
		cols_to_be_tossed_out = set()
		for i in range(len(snp_acc_ls)):
			chr_pos = snp_acc_ls[i]
			chr_pos = chr_pos.split('_')
			chr_pos = map(int, chr_pos)
			pos = chr_pos[1]
			if not (chr_pos[0]==chr and pos>=start_pos and pos<=stop_pos):
				cols_to_be_tossed_out.add(i)
		
		write_data_matrix(snpData.data_matrix, output_fname, snpData.header, snpData.strain_acc_list, snpData.category_list, cols_to_be_tossed_out=cols_to_be_tossed_out)
	
	"""
input_fname= '/Network/Data/250k/tmp-yh/call_method_17.tsv'
from pymodule import SNPData
snpData = SNPData(input_fname=input_fname, turn_into_integer=1, turn_into_array=1)
chr=4
start_pos=100000
stop_pos=700000
output_fname = '/Network/Data/250k/tmp-yh/250k_data/call_method_17_chr%s_%s_%s.tsv'%(chr, start_pos, stop_pos)
FileFormatExchange.outputSNPmatrixGivenRegion(snpData, output_fname, chr, start_pos, stop_pos)
	"""

	def reduce250kOnlyToCertainSNPs(cls, snp_id_ls):
		#2008-10-07 form a smaller 250k test dataset for PlotGroupOfSNPs.py, snp_id_ls is the top 200 snps from (KW,LD)
		
		from pymodule import SNPData, write_data_matrix
		input_fname = '/Network/Data/250k/tmp-yh/call_method_17.tsv'
		snpData = SNPData(input_fname=input_fname, turn_into_integer=1, turn_into_array=1)
		output_fname = '/Network/Data/250k/tmp-yh/call_method_17_test.tsv'
		from sets import Set
		good_snp_id_set = Set(snp_id_ls)
		col_index_to_be_tossed_set = Set()
		for snp_id, col_index in snpData.col_id2col_index.iteritems():
			if snp_id not in good_snp_id_set:
				col_index_to_be_tossed_set.add(col_index)
		write_data_matrix(snpData.data_matrix, output_fname, snpData.header, snpData.strain_acc_list, snpData.category_list, 
						rows_to_be_tossed_out=None, cols_to_be_tossed_out=col_index_to_be_tossed_set)
	"""
FileFormatExchange.reduce250kOnlyToCertainSNPs = classmethod(reduce250kOnlyToCertainSNPs)
	"""
	
	"""
	2007-03-28
	"""
	@classmethod
	def strip_2010_strain_info(cls, input_fname, output_fname):
		import csv
		reader = csv.reader(open(input_fname), delimiter='\t')
		#skip the 1st two sentences
		reader.next()
		reader.next()
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		for row in reader:
			for i in range(len(row)):
				row[i] = row[i].strip()
			writer.writerow(row)
		del reader, writer
	
	"""
FileFormatExchange.strip_2010_strain_info('./script/variation/data/2010/2010_strain_info.csv', './script/variation/data/2010/2010_strain_info_stripped.csv')
	"""
	
	"""
	2007-03-15
	"""
	def reformat_data_for_chris(curs, input_fname, output_fname, snp_locus_table):
		from common import number2nt
		
		sys.stderr.write("Getting snp_acc2pos ...")
		snp_acc2pos = {}
		curs.execute("select acc, chromosome, position from %s"%snp_locus_table)
		rows = curs.fetchall()
		for row in rows:
			snp_acc, chromosome, position = row
			snp_acc2pos[snp_acc] = [chromosome, position]
		sys.stderr.write("Done.\n")
		
		import csv
		reader = csv.reader(open(input_fname), delimiter='\t')
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = reader.next()
		new_header = [header[0]] + header[2:]
		writer.writerow(new_header)
		
		chr_ls = ['']
		pos_ls = ['']
		for snp_acc in new_header[1:]:
			chr_ls.append(snp_acc2pos[snp_acc][0])
			pos_ls.append(snp_acc2pos[snp_acc][1])
		writer.writerow(chr_ls)
		writer.writerow(pos_ls)
		
		for row in reader:
			new_row = [row[0]]
			for call in row[2:]:
				nt = number2nt[int(call)]
				if nt=='NA':	#chris wants 'N' instead of 'NA'
					nt = 'N'
				new_row.append(nt)
			writer.writerow(new_row)
		del writer
	
	"""
	hostname='zhoudb'
	dbname='graphdb'
	schema = 'dbsnp'
	conn, curs = db_connect(hostname, dbname, schema)
	reformat_data_for_chris(curs, '../data/justin_data_filtered.csv', '../data/justin_data_filtered_for_chris.csv', 'snp_locus')
	"""


	@classmethod
	def removeRowsBasedOnPhenotypeValue(input_fname, phenotype_fname, phenotype_method_id, output_fname, min_pheno_value=None, max_pheno_value=None):
		"""
		2009-4-13 investigate whether throwing off some rows help to increase significance.
			similar purpose to removeRowsBasedOnSNPAllele()
			throw away accessions whose phenotypes are not in the phenotype range specified
		"""
		from pymodule import read_data, write_data_matrix, SNPData
		header, strain_acc_list, category_list, data_matrix = read_data(input_fname, turn_into_integer=1)
		
		snpData = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list,\
						data_matrix=data_matrix)
		header_phen, strain_acc_list_phen, category_list_phen, data_matrix_phen = read_data(phenotype_fname, turn_into_integer=0)
		from variation.src.Association import Association
		data_matrix_phen = Association.get_phenotype_matrix_in_data_matrix_order(strain_acc_list, strain_acc_list_phen, data_matrix_phen)
		phenData = SNPData(header=header_phen, strain_acc_list=snpData.strain_acc_list, data_matrix=data_matrix_phen)
		
		from variation.src.PlotGroupOfSNPs import PlotGroupOfSNPs
		which_phenotype_ls = PlotGroupOfSNPs.findOutWhichPhenotypeColumn(phenData, set([phenotype_method_id]))
		which_phenotype = which_phenotype_ls[0]
		
		import numpy
		rows_to_be_tossed_out = set()
		for i in range(len(phenData.data_matrix)):
			phenotype_value = phenData.data_matrix[i][which_phenotype]
			if numpy.isnan(phenotype_value):
				rows_to_be_tossed_out.add(i)
				continue
			
			if min_pheno_value is not None and phenotype_value<min_pheno_value:
				rows_to_be_tossed_out.add(i)
				continue
			if max_pheno_value is not None and phenotype_value>max_pheno_value:
				rows_to_be_tossed_out.add(i)
				continue
		
		write_data_matrix(data_matrix, output_fname, header, strain_acc_list, category_list, rows_to_be_tossed_out=rows_to_be_tossed_out)
		
	"""
	input_fname='/Network/Data/250k/tmp-yh/250k_data/call_method_17_test.tsv'
	input_fname='/Network/Data/250k/tmp-yh/250k_data/call_method_29.tsv'
	phenotype_fname='/Network/Data/250k/tmp-yh/phenotype_no_transform.tsv'
	phenotype_method_id=44
	min_pheno_value = 0.45
	output_fname = '/tmp/call_method_29_pheno_id_%s_above_%s.tsv'%(phenotype_method_id, min_pheno_value)
	FileFormatExchange.removeRowsBasedOnPhenotypeValue(input_fname, phenotype_fname, phenotype_method_id, output_fname, min_pheno_value)
	"""
	
	@classmethod
	def reverseComplementFastaInput(cls, input_fname, output_fname):
		"""
		2009-3-25
			read in sequences from input_fname (in fasta format), reverse complement each sequence and output them 
		"""
		inf = open(input_fname)
		of = open(output_fname, 'w')
		from Bio import SeqIO
		for seq_record in SeqIO.parse(inf, "fasta") :
			seq_rc = seq_record.seq.reverse_complement()
			of.write('>%s\n'%seq_record.id)
			of.write('%s\n'%seq_rc.tostring())
	
	
	"""
	input_fname = os.path.expanduser('~/script/variation/doc/FLC/CarolineDeanFLC.seq')
	output_fname = os.path.expanduser('~/script/variation/doc/FLC/CarolineDeanFLC_rc.seq')
	AnalyzeSNPData.reverseComplementFastaInput(input_fname, output_fname)
	"""
	
	
	"""
	2009-3-20
		output fasta format in DrawSNPMatrix.py input format.
		also partition the input into blocks each with number of columns <= maxNumberColsPerBlock.
		too many columns render some parts of the image generated by DrawSNPMatrix.py unviewable.
	"""
	@classmethod
	def fasta2DrawSNPMatrixFormat(cls ,input_fname, output_fname, maxNumberColsPerBlock=2000):
		import csv, math
		blockIndex2writer = {}
		inf = open(input_fname)
		from Bio import SeqIO
		isHeaderOutputted = False
		header_row = []
		for seq_record in SeqIO.parse(inf, "fasta") :
			one_row = []
			no_of_blocks = int(math.ceil(len(seq_record.seq)/float(maxNumberColsPerBlock)))
			for i in range(len(seq_record.seq)):
				nt = seq_record.seq[i]
				if not isHeaderOutputted:
					header_row.append(i+1)
				one_row.append(nt)
			
			#output header
			if not isHeaderOutputted:
				for i in range(no_of_blocks):
					start_index = i*maxNumberColsPerBlock
					stop_index = (i+1)*maxNumberColsPerBlock
					row_for_this_block = header_row[start_index:stop_index]
					if i not in blockIndex2writer:
						output_fname_ls = os.path.splitext(output_fname)
						block_output_fname = '%s_%s_%s%s'%(output_fname_ls[0], start_index, stop_index, output_fname_ls[1])
						writer = csv.writer(open(block_output_fname, 'w'), delimiter='\t')
						blockIndex2writer[i] = writer
					blockIndex2writer[i].writerow(['name',1]+row_for_this_block)
				isHeaderOutputted = True
			
			#output data
			for i in range(no_of_blocks):
				start_index = i*maxNumberColsPerBlock
				stop_index = (i+1)*maxNumberColsPerBlock
				row_for_this_block = one_row[start_index:stop_index]
				blockIndex2writer[i].writerow([seq_record.id, 1] + row_for_this_block)
		
		del blockIndex2writer, inf
	
	"""
	input_fname = os.path.expanduser('~/script/variation/doc/FLC/CarolineDeanFLCAlign.fasta')
	output_fname = os.path.expanduser('~/script/variation/doc/FLC/CarolineDeanFLCAlign.tsv')
	AnalyzeSNPData.fasta2DrawSNPMatrixFormat(input_fname, output_fname)
	"""
	
	@classmethod
	def splitMultiRecordFastaFileIntoSingleRecordFastaFiles(cls, input_fname, output_dir):
		"""
		2010-12-3
			One fasta file contains >1 sequences. This function writes each sequence to an individual file in output_dir.
			The sequence title would be the filename in output_dir.
			
			SyMAP requires this kind of splitting for it to work on fast files.
		"""
		import os,sys
		sys.stderr.write("Outputting each sequence in %s into single file ...\n"%(input_fname))
		inf = open(input_fname, 'r')
		outf = None
		counter = 0
		real_counter = 0
		for line in inf:
			counter += 1
			if line[0]=='>':
				real_counter += 1
				title = line[1:].strip()
				if outf is not None:
					outf.close()
				output_fname = os.path.join(output_dir, '%s.seq'%title)
				outf = open(output_fname, 'w')
			if outf is not None:
				outf.write(line)
			if counter%10000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
		sys.stderr.write("Done.\n")
	
	"""
		#2010-12-3
		input_fname = os.path.expanduser("~/script/variation/data/CNV/TAIR9/tair9.fas")
		output_dir = os.path.expanduser("~/script/variation/bin/symap_3_3/data/pseudo/thaliana/sequence/pseudo/")
		FileFormatExchange.splitMultiRecordFastaFileIntoSingleRecordFastaFiles(input_fname, output_dir)
		sys.exit(0)
		
		input_fname = os.path.expanduser("~/script/variation/data/lyrata/Araly1_assembly_scaffolds.fasta")
		output_dir = os.path.expanduser("~/script/variation/bin/symap_3_3/data/pseudo/lyrata/sequence/pseudo/")
		FileFormatExchange.splitMultiRecordFastaFileIntoSingleRecordFastaFiles(input_fname, output_dir)
		sys.exit(0)
		
		#2010-12-3
		input_fname = os.path.expanduser("~/script/variation/data/lyrata/Araly1_assembly_scaffolds.fasta")
		output_dir = os.path.expanduser("~/script/variation/bin/symap_3_3/data/pseudo/lyrata/sequence/pseudo/")
		FileFormatExchange.splitMultiRecordFastaFileIntoSingleRecordFastaFiles(input_fname, output_dir)
		sys.exit(0)
		
		#2011-1-31
		input_fname = os.path.expanduser("/Network/Data/NCBI/hs_genome.fasta")
		output_dir = os.path.expanduser("~/script/variation/bin/symap_3_3/data/pseudo/hg19/sequence/pseudo/")
		FileFormatExchange.splitMultiRecordFastaFileIntoSingleRecordFastaFiles(input_fname, output_dir)
		
		input_fname = os.path.expanduser("/Network/Data/NCBI/mm_genome.fasta")
		output_dir = os.path.expanduser("~/script/variation/bin/symap_3_3/data/pseudo/macaque/sequence/pseudo/")
		FileFormatExchange.splitMultiRecordFastaFileIntoSingleRecordFastaFiles(input_fname, output_dir)
		
		sys.exit(0)
		
	"""
		
class Data250k(object):
	def __init__(self):
		"""
		2009-11-14
			keep the old name alive
		"""
		choose192OutOf250k = self.retainSNPDataFromCertainEcotypes
	
	"""
	2008-08-16
		check if the array_ids and ecotype_ids in call files generated by bjarni match the ones in db
	"""
	@classmethod
	def checkBjarniFile(cls, input_fname, curs):
		import csv
		reader = csv.reader(open(input_fname))
		array_id_ls = reader.next()[2:]
		array_id_ls = map(int, array_id_ls)
		print "%s arrays"%(len(array_id_ls))
		ecotype_id_ls = reader.next()[2:]
		ecotype_id_ls = map(int, ecotype_id_ls)
		print "%s ecotypes"%(len(ecotype_id_ls))
		for i in range(len(array_id_ls)):
			array_id = array_id_ls[i]
			ecotype_id = ecotype_id_ls[i]
			curs.execute("select tg_ecotypeid from stock.ecotypeid2tg_ecotypeid where ecotypeid=%s"%ecotype_id)
			rows = curs.fetchall()
			if rows:
				bjarni_ecotype_id = rows[0][0]
			else:
				sys.stderr.write( "ecotype_id %s has no tg_ecotypeid in ecotypeid2tg_ecotypeid.\n"%(ecotype_id))
				bjarni_ecotype_id = ecotype_id
			curs.execute("select maternal_ecotype_id from stock_250k.array_info where id=%s"%array_id)
			rows = curs.fetchall()
			maternal_ecotype_id = rows[0][0]
			if bjarni_ecotype_id!=maternal_ecotype_id:
				print i, array_id, bjarni_ecotype_id, maternal_ecotype_id
		del reader
	
	"""
	input_fname = '/Network/Data/250k/dataFreeze_080608/250K_f8_080608.csv'
	checkBjarniFile(input_fname, curs)
	"""
	
	@classmethod
	def recoverArrayID2ndCol(cls, input_fname, data_with_array_id_fname, output_fname):
		"""
		2008-01-04 recover array id (2nd column) of a StrainXSNP data file based on its old version
			the 1st column (strain id) of input_fname has to be unique,
			otherwise a random one among duplicates would get a possibly wrong array id assigned.
		"""
		from pymodule import read_data, SNPData
		header, strain_acc_list, category_list, data_matrix = read_data(input_fname)
		snpData = SNPData(header=header, strain_acc_list=strain_acc_list, data_matrix=data_matrix)
		#
		header2, strain_acc_list2, category_list2, data_matrix2 = read_data(data_with_array_id_fname)
		#put the array_id from category_list2 into corresponding position in category_list
		for i in range(len(strain_acc_list2)):
			ecotype_id = strain_acc_list2[i]
			array_id = category_list2[i]
			row_index = snpData.row_id2row_index.get(ecotype_id)	#ecotype_id might not be in input_fname
			if row_index is not None:
				category_list[row_index] = array_id
		
		newSNPData = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list,\
							data_matrix=data_matrix)
		newSNPData.tofile(output_fname)
	
	"""
	input_fname = os.path.expanduser('~/panfs/NPUTE_data/input/250k_l3_y.85_uniq_ecotype_20080919_3_FRI_del_no_array_id.tsv')
	data_with_array_id_fname = os.path.expanduser('~/panfs/NPUTE_data/input/250k_l3_y.85_uniq_ecotype_20080919.tsv')
	output_fname = os.path.expanduser('~/panfs/NPUTE_data/input/250k_l3_y.85_uniq_ecotype_20080919_3_FRI_del.tsv')
	Data250k.recoverArrayID2ndCol(input_fname, data_with_array_id_fname, output_fname)
	"""
	
	@classmethod
	def chooseFirst96(cls, input_fname, db, output_fname):
		"""
		2009-2-17
			pick 1st 96 2010 accessions out of 250k data
		"""
		from pymodule import SNPData
		snp_data = SNPData(input_fname=input_fname, turn_into_array=1, ignore_2nd_column=1)
		rows = db.metadata.bind.execute("select * from at.accession2tg_ecotypeid where accession_id<=96")
		rows_to_be_preserved = set()
		for row in rows:
			row_id = '%s'%row.ecotype_id
			if row_id in snp_data.row_id2row_index:
				rows_to_be_preserved.add(snp_data.row_id2row_index[row_id])
		rows_to_be_tossed_out = set(range(len(snp_data.row_id_ls))) - rows_to_be_preserved
		snp_data.tofile(output_fname, rows_to_be_tossed_out=rows_to_be_tossed_out)
	
	"""
	input_fname = os.path.expanduser('~/panfs/250k/call_method_29.tsv')
	output_fname = os.path.expanduser('~/panfs/250k/call_method_29_only96.tsv')
	Data250k.chooseFirst96(input_fname, db_250k, output_fname)
	"""
	@classmethod
	def filterGWAToRetainSegregatingSNPs(cls, db_250k, call_method_id, phenotype_method_id, analysis_method_id, snpData, \
										 ecotype_id1, ecotype_id2, output_fname):
		"""
		2009-3-12
			read in a GWA result, and its affiliated genotype matrix
			and retain SNPs in GWA result that are segregating between ecotypes
		"""
		gwar = db_250k.getGWA(call_method_id, phenotype_method_id, analysis_method_id)
		row_index1 = snpData.row_id2row_index[str(ecotype_id1)]
		row_index2 = snpData.row_id2row_index[str(ecotype_id2)]
		
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		writer.writerow(['chromosome', 'position', 'pvalue'])
		for data_obj in gwar.data_obj_ls:
			snp_id = '%s_%s'%(data_obj.chromosome, data_obj.position)
			col_index = snpData.col_id2col_index[snp_id]
			allele1 = snpData.data_matrix[row_index1][col_index]
			allele2 = snpData.data_matrix[row_index2][col_index]
			if allele1!=allele2:
				writer.writerow([data_obj.chromosome, data_obj.position, data_obj.value])
		del writer

	"""
	call_method_id = 29
	phenotype_method_id=1
	analysis_method_id=7
	ecotype_id1=6917	'Fab-2'
	ecotype_id2=6904	'Br-0'
	output_fname = '/tmp/gwa_call_%s_phenotype_%s_analysis_%s_segregating_in_%s_%s.tsv'%(call_method_id, phenotype_method_id, analysis_method_id, ecotype_id1, ecotype_id2)
	snpData = db_250k.getSNPMatrix(call_method_id)
	Data250k.filterGWAToRetainSegregatingSNPs(db_250k, call_method_id, phenotype_method_id, analysis_method_id, snpData, ecotype_id1, ecotype_id2, output_fname)
	"""
	
	@classmethod
	def getEcotypeSetFromFile(cls, ecotype_id_fname):
		"""
		arguments:
			ecotype_id_fname: coma-delimited file which contains ecotype-ids to be kept. Only first column is read-in.
		
		2010-2-22
			split out of retainSNPDataFromCertainEcotypes()	
			
		"""
		import csv
		sys.stderr.write("Getting ecotype IDs from %s ..."%ecotype_id_fname)
		reader = csv.reader(open(ecotype_id_fname), delimiter=',')
		reader.next()	#toss the header
		ecotype_set = set()
		for row in reader:
			ecotype_set.add(int(row[0]))
		del reader
		sys.stderr.write("%s ecotypes.\n"%len(ecotype_set))
		return ecotype_set
	
	@classmethod
	def getEcotypeSetFromDBAccessionSet(cls, db_250k=None, accessionSetID=None):
		"""
		2010-2-22
			
		"""
		sys.stderr.write("Getting ecotype IDs given accession set %s from db ..."%accessionSetID)
		import Stock_250kDB
		accession_set = Stock_250kDB.AccessionSet.get(accessionSetID)
		ecotype_set = set()
		for ecotype in accession_set.ecotype_ls:
			ecotype_set.add(ecotype.ecotype_id)
		sys.stderr.write("%s ecotypes.\n"%len(ecotype_set))
		return ecotype_set
		 
	
	@classmethod
	def retainSNPDataFromCertainEcotypes(cls, input_fname, ecotype_id_fname=None, output_fname=None, removeMonomorphicCols=False, \
										accessionSetID=None, db_250k=None):
		"""
		Function:
			1. retain 250k data (input_fname) whose ecotype IDs are contained in ecotype_id_fname (one-line-header)
			2. output the filtered data into the new output_fname.
		arguments:
			input_fname: original 250k data in Strain X SNP format, whose first column is regarded as ecotype-id.
			ecotype_id_fname: coma-delimited file which contains ecotype-ids to be kept. Only first column is read-in.
			output_fname: where the output will be, Strain X SNP format
			removeMonomorphicCols: if True, will remove monomorphic SNPs after removing un-wanted ecotypes
		2010-2-22
			add argument accessionSetID, if available, get the ecotype_set from db
			add argument db_250k. must be valid if accessionSetID is given.
		2010-2-21
			add argument removeMonomorphicCols
			remove un-wanted rows, rather than ignore them during output
		2009-11-14
			function renamed from choose192OutOf250k()
		2009-3-27
		"""
		from pymodule import SNPData
		if accessionSetID:
			ecotype_set = cls.getEcotypeSetFromDBAccessionSet(db_250k, accessionSetID)
		elif ecotype_id_fname:
			ecotype_set = cls.getEcotypeSetFromFile(ecotype_id_fname)
		else:
			sys.stderr.write("Neither ecotype_id_fname nor accessionSetID is given. aborted.\n")
			return
		
		snp_data = SNPData(input_fname=input_fname, turn_into_array=1)	#2nd column array id is also needed
		
		rows_to_be_tossed_out = set()	# 2010-2-21 not used anymore.
		row_ids_to_be_kept = set()	# 2010-2-21
		for row_id in snp_data.row_id_ls:
			ecotype_id = int(row_id[0])	#1st column is ecotype_id, 2nd is array id
			if ecotype_id not in ecotype_set:
				rows_to_be_tossed_out.add(snp_data.row_id2row_index[row_id])
			else:
				row_ids_to_be_kept.add(row_id)
		
		snpData1 = SNPData.keepRowsByRowID(snp_data, row_ids_to_be_kept)
		if removeMonomorphicCols:
			snpData1 = SNPData.removeMonomorphicCols(snpData1)
		snpData1.tofile(output_fname)
	
	"""
	input_fname = '/Network/Data/250k/db/dataset/call_method_29.tsv'
	fname_with_192_ids = '/tmp/192_accessions_031009.csv' # email from bjarni. coma delimited, 1st column is ecotype id. actually 199 accessions.
	output_fname = '/tmp/call_method_29_with_192_accession.tsv'
	Data250k.retainSNPDataFromCertainEcotypes(input_fname, fname_with_192_ids, output_fname)
	
	# 2010-2-15 Yan Li's stuff
	input_fname = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49.tsv')
	ecotype_id_fname = os.path.expanduser('~/mnt/banyan/tmp/core482_tg_ecotypeid.csv')
	output_fname = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482.tsv')
	Data250k.retainSNPDataFromCertainEcotypes(input_fname, ecotype_id_fname, output_fname, removeMonomorphicCols=True)
	
	# 2010-2-22 for suzi's swedish versus non-swe GWAS
	input_fname = os.path.expanduser('/Network/Data/250k/db/dataset/call_method_43.tsv')
	output_fname = os.path.expanduser('/Network/Data/250k/db/dataset/call_method_43_Suzi295SWE.tsv')
	Data250k.retainSNPDataFromCertainEcotypes(input_fname, output_fname=output_fname, removeMonomorphicCols=True, \
										accessionSetID=4, db_250k=db_250k)
	
	input_fname = os.path.expanduser('/Network/Data/250k/db/dataset/call_method_43.tsv')
	output_fname = os.path.expanduser('/Network/Data/250k/db/dataset/call_method_43_Suzi295NoneSWE.tsv')
	Data250k.retainSNPDataFromCertainEcotypes(input_fname, output_fname=output_fname, removeMonomorphicCols=True, \
										accessionSetID=5, db_250k=db_250k)							
	"""
	
	@classmethod
	def filterCallMethod1ToHaveSameArraysAsCallMethod2(cls, call_method_id1, call_method_id2, output_fname=None,\
											call_method1_fname=None):
		"""
		2011-4-25
			add argument call_method1_fname. if given, it replaces call_method1.filename
		2009-5-19
			 generate a new snp dataset based on call_method_id1 but with only the arrays that are in call_method_id2
		"""
		import Stock_250kDB
		call_method2 = Stock_250kDB.CallMethod.get(call_method_id2)
		from pymodule import SNPData
		if call_method1_fname:
			input_fname = call_method1_fname
		else:
			call_method1 = Stock_250kDB.CallMethod.get(call_method_id1)
			input_fname=call_method1.filename
		snpData1 = SNPData(input_fname=input_fname, turn_into_array=1)
		
		sys.stderr.write("Selecting overlapping arrays ...")
		array_id_in_call_method2_set = set()
		for call_info in call_method2.call_info_ls:
			array_id_in_call_method2_set.add(call_info.array_id)
		
		#print array_id_in_call_method2_set
		
		row_id_wanted_ls = []
		for row_id in snpData1.row_id_ls:
			array_id = int(row_id[1])
			if array_id in array_id_in_call_method2_set:
				row_id_wanted_ls.append(row_id)
		#print row_id_wanted_ls
		sys.stderr.write("Done.\n")
		
		newSNPData = SNPData.keepRowsByRowID(snpData1, row_id_wanted_ls)
		newSNPData.tofile(output_fname)
		del newSNPData
		del snpData1
	
	"""
		#2009
		call_method_id1 = 34
		call_method_id2 = 33
		output_fname = '/tmp/call_method_%s_arrays_from_%s.tsv'%(call_method_id1, call_method_id2)
		Data250k.filterCallMethod1ToHaveSameArraysAsCallMethod2(call_method_id1, call_method_id2, output_fname)
		sys.exit(0)
		
		
		#2011-4-25
		cnv_method_id = 20
		output_fname = os.path.expanduser('~/script/variation/data/CNV/NonOverlapCNVAsSNP_cnvMethod%s_snpID.tsv'%cnv_method_id)
		#CNV.outputNonOverlappingCNVAsSNP(db_250k, output_fname, cnv_method_id=cnv_method_id, cnv_type_id=1)
		
		call_method_id1 = None
		call_method_id2 = 32
		call_method1_fname = output_fname
		output_fname = os.path.expanduser('~/script/variation/data/CNV/cnvMethod%s_arrays_only_in_call%s.tsv'%\
				(cnv_method_id, call_method_id2))
		Data250k.filterCallMethod1ToHaveSameArraysAsCallMethod2(call_method_id1, call_method_id2, output_fname,
				call_method1_fname=call_method1_fname)
		sys.exit(0)
		
	"""
	
	@classmethod
	def addUnimputedArrayCallsToImputedCalls(cls, db_250k, imputed_call_method_id=54, unimputed_array_ids='1516-1602',\
											output_fname=None, minCallProbability=0.85):
		"""
		2010-7-2
			for Rhonda Meyer's stuff
		"""
		sys.stderr.write("Adding unimputed array calls to imputed ones ....")
		import Stock_250kDB, csv
		from pymodule import figureOutDelimiter, getListOutOfStr, SNPData
		from pymodule.SNP import nt2number
		imputed_call_method = Stock_250kDB.CallMethod.get(imputed_call_method_id)
		reader = csv.reader(open(imputed_call_method.filename), \
						delimiter=figureOutDelimiter(imputed_call_method.filename))
		
		sys.stderr.write("Constructing chr_pos2col_index ...")
		header = reader.next()
		snp_chr_pos_ls = header[2:]
		no_of_snps = len(snp_chr_pos_ls)
		chr_pos2col_index = {}
		for i in range(no_of_snps):
			chr_pos = snp_chr_pos_ls[i]
			chr_pos2col_index[chr_pos] = i
		sys.stderr.write("Done.\n")
		
		
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		writer.writerow(header)
		
		sys.stderr.write("Outputting unimputed calls  ....\n")
		unimputed_array_ids = getListOutOfStr(unimputed_array_ids, data_type=int)
		for i in xrange(len(unimputed_array_ids)):
			array_id = unimputed_array_ids[i]
			call_info = Stock_250kDB.CallInfo.query.filter_by(method_id=3).filter_by(array_id=array_id).first()
			sys.stderr.write("%s%s\tarray %s"%('\x08'*80, i, call_info.array_id))
			tmp_reader = csv.reader(open(call_info.filename), delimiter='\t')
			tmp_reader.next()
			snp_call_ls = [0]*no_of_snps
			for row in tmp_reader:
				snp_id = row[0].split('_')	# 1_657_C_T
				snp_id = '_'.join(snp_id[:2])	# remove _C_T
				snp_col_index = chr_pos2col_index.get(snp_id)
				if snp_col_index is not None:
					oligo_probability = float(row[2])
					if oligo_probability>=minCallProbability:
						call_integer = nt2number[row[1]]
					else:
						call_integer = 0
					snp_call_ls[snp_col_index] = call_integer
			data_row = ['', call_info.array_id,] + snp_call_ls
			writer.writerow(data_row)
		sys.stderr.write("Done.\n")
		
		sys.stderr.write("Outputting imputed calls of method %s ....\n" %(imputed_call_method_id))
		counter = 0
		for row in reader:
			sys.stderr.write("%s%s"%('\x08'*80, counter))
			writer.writerow(row)
			counter += 1
		del writer
		sys.stderr.write("Done.\n")
	
	"""
		#2010-7-2
		imputed_call_method_id=54
		output_fname = os.path.expanduser('~/script/variation/data/call_%s_with_RhondaUnimputedCalls.tsv'%(imputed_call_method_id))
		Data250k.addUnimputedArrayCallsToImputedCalls(db_250k, imputed_call_method_id=imputed_call_method_id, \
				unimputed_array_ids='1516-1602',\
				output_fname=output_fname, minCallProbability=0.85)
		sys.exit(0)
		
		#2010-7-2 add Rhonda Myers' unimputed arrays into call 54
		imputed_call_method_id=54
		output_fname = os.path.expanduser('~/script/variation/data/call_%s_with_RhondaUnimputedCalls.tsv'%(imputed_call_method_id))
		Data250k.addUnimputedArrayCallsToImputedCalls(db_250k, imputed_call_method_id=imputed_call_method_id, \
				unimputed_array_ids='1516-1602',\
				output_fname=output_fname, minCallProbability=0.85)
		sys.exit(0)
		
		#2010-10-10 add Martin Koornneef's un-imputed arrays into call 54
		imputed_call_method_id=54
		output_fname = os.path.expanduser('~/script/variation/data/call_%s_with_KoornneefUnimputedCalls.tsv'%(imputed_call_method_id))
		Data250k.addUnimputedArrayCallsToImputedCalls(db_250k, imputed_call_method_id=imputed_call_method_id, \
				unimputed_array_ids='1603-1619',\
				output_fname=output_fname, minCallProbability=0.85)
		sys.exit(0)
		
	"""
class BooleanInteraction(object):
	
	"""
	2008-08-05
		test boolean relationship between SNPs
	"""
	
	def returnTop2Allele(snp_allele2count):
		"""
		2008-08-06 remove redundant argument snp_allele_ls
		2008-08-05 in the descending order of count for each allele, assign index ascending
		"""
		snp_allele_count_ls = []
		snp_allele_ls = snp_allele2count.keys()
		for snp_allele in snp_allele_ls:
			snp_allele_count_ls.append(snp_allele2count[snp_allele])
		import numpy
		argsort_ls = numpy.argsort(snp_allele_count_ls)
		new_snp_allele2index = {}
		for i in [-1, -2]:
			snp_index = argsort_ls[i]	#-1 is index for biggest, -2 is next biggest
			new_snp_allele2index[snp_allele_ls[snp_index]] = -i-1
		return new_snp_allele2index
	
	@classmethod
	def booleanMergeSNPs(cls, input_fname, output_fname, SNP_label1, SNP_label2, operator_type=1):	#1 is and, 2 is or
		"""
		2008-08-05
			alleles not in the top 2 are taken as NA. major allele is coded as 0. minor allele is coded as 1.
		"""
		from pymodule import read_data, write_data_matrix, nt2number
		header, strain_acc_list, category_list, data_matrix = read_data(input_fname, turn_into_integer=1)
		
		snp_acc_ls = header[2:]
		snp_index1 = -1
		snp_index2 = -1
		for i in range(len(snp_acc_ls)):
			if snp_acc_ls[i]==SNP_label1:
				snp_index1 = i
			if snp_acc_ls[i]==SNP_label2:
				snp_index2 = i
		
		snp_allele2count1 = {}
		snp_allele2count2 = {}
		no_of_rows = len(data_matrix)
		
		for i in range(no_of_rows):
			snp1_allele = data_matrix[i][snp_index1]
			snp2_allele = data_matrix[i][snp_index2]
			if snp1_allele!=0:
				if snp1_allele not in snp_allele2count1:
					snp_allele2count1[snp1_allele] = 0
				snp_allele2count1[snp1_allele] += 1
			if snp2_allele!=0:
				if snp2_allele not in snp_allele2count2:
					snp_allele2count2[snp2_allele] = 0
				snp_allele2count2[snp2_allele] += 1
		print snp_allele2count1
		print snp_allele2count2
		snp_allele2index1 = returnTop2Allele(snp_allele2count1)
		
		snp_allele2index2 = returnTop2Allele(snp_allele2count2)
		print snp_allele2index1
		print snp_allele2index2
		
		no_of_cols = 1
		new_data_matrix = data_matrix	#replace the 1st SNP's data with the new boolean result
		for i in range(no_of_rows):
			snp1_allele = data_matrix[i][snp_index1]
			snp2_allele = data_matrix[i][snp_index2]
			if snp1_allele in snp_allele2index1:
				snp_code1 = snp_allele2index1[snp1_allele]
			else:
				snp_code1 = 0
			
			if snp2_allele in snp_allele2index2:
				snp_code2 = snp_allele2index2[snp2_allele]
			else:
				snp_code2 = 0
			if operator_type==1:
				if snp1_allele in snp_allele2index1 and snp2_allele in snp_allele2index2:
					new_data_matrix[i][snp_index1] = (snp_code1 and snp_code2) + 1
			elif operator_type==2:
				if snp1_allele in snp_allele2index1 or snp2_allele in snp_allele2index2:
					new_data_matrix[i][snp_index1] = (snp_code1 or snp_code2) + 1
		write_data_matrix(new_data_matrix, output_fname, header, strain_acc_list, category_list)
	
	"""
	input_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix.tsv')
	output_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix_1_or_2.tsv')
	SNP_label1 = '4_268809_0'
	SNP_label2 = '4_269962_8'
	booleanMergeSNPs(input_fname, output_fname, SNP_label1, SNP_label2, operator_type=2)
	input_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix_1_or_2.tsv')
	output_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix_1_or_2_or_3.tsv')
	SNP_label3 = '4_270712_0'
	booleanMergeSNPs(input_fname, output_fname, SNP_label1, SNP_label3, operator_type=2)
	"""
	
	
	@classmethod
	def handleMostSignificantOperatorPerSNPPair(cls, input_fname, afterMatchFunction, param_obj):
		"""
		2009-2-19
			a general function which calls afterMatchFunction to do a bit processing for the most significant operator of each SNP pair.
				(results from the same SNP pair are always cluttered together.)
		"""
		import csv
		reader = csv.reader(open(input_fname), delimiter='\t')
		reader.next()
		prev_snp2_id = None
		prev_snp1_id = None
		min_pvalue = None
		row_with_min_pvalue = None
		i = 0
		#input_fname_basename = os.path.splitext(os.path.basename(input_fname))[0]	#function calls this general function should supply this
		#param_obj.input_fname_basename = input_fname_basename
		counter = 0
		real_counter = 0
		for row in reader:
			snp1_id, gene1_id, snp2_id, gene2_id, bool_type, pvalue, count1, count2 = row[:8]
			counter += 1
			pvalue = float(pvalue)
			if prev_snp2_id is None:	#first snp2
				prev_snp1_id = snp1_id
				prev_snp2_id = snp2_id
				min_pvalue = pvalue
				row_with_min_pvalue = row
			elif snp2_id == prev_snp2_id and snp1_id==prev_snp1_id:	#same snp2
				i += 1
				if pvalue<min_pvalue:
					min_pvalue=pvalue
					row_with_min_pvalue = row
			else:	#new pairs
				i = 0
				prev_snp1_id = row_with_min_pvalue[0]
				afterMatchFunction(row_with_min_pvalue, param_obj)
				real_counter += 1
				prev_snp1_id = snp1_id
				prev_snp2_id = snp2_id
				min_pvalue = pvalue
				row_with_min_pvalue = row
				if real_counter%50000==0:
					sys.stderr.write("%s%s\t\t%s"%('\x08'*40, real_counter, counter))
	
	"""
	2009-1-25
		input_fname is output of MpiIntraGeneSNPPairAsso.py or MpiInterGeneSNPPairAsso.py or MpiInterSNPPairAsso.py
		
		separate one snp-pairwise association results into single-SNP GWA files according to the snp1_id.
	"""
	@classmethod
	def outputRow(cls, row, param_obj):
		snp_id2writer = getattr(param_obj, "snp_id2writer")
		output_dir = getattr(param_obj, "output_dir")
		input_fname_basename = getattr(param_obj, "input_fname_basename")
		header = getattr(param_obj, "header")
		
		snp1_id, gene1_id, snp2_id, gene2_id, bool_type, pvalue, count1, count2 = row[:8]
		if snp1_id not in snp_id2writer:
			output_fname = os.path.join(output_dir, '%s_vs_%s.tsv'%(input_fname_basename, snp1_id))
			writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
			writer.writerow(header)
			snp_id2writer[snp1_id] = writer
		writer = snp_id2writer[snp1_id]
		chr, pos = snp2_id.split('_')
		count1 = int(count1)
		count2 = int(count2)
		mac = min(count1, count2)
		maf = float(mac)/(count1+count2)
		row = [chr, pos, pvalue, maf, mac, bool_type]
		writer.writerow(row)
	
	@classmethod
	def outputBooleanPairIntoGWAFormat(cls, input_fname, output_dir, pos_index=1, need_beta=True, min_value_cutoff=None, do_log10_transformation=False):
		"""
		2009-1-25
			call outputRow()
		"""
		import csv
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		snp_id2writer = {}	#each SNP1_ID would have a separate file for results from all SNPs that tested with that SNP
		
		input_fname_basename = os.path.splitext(os.path.basename(input_fname))[0]
		from pymodule import PassingData
		param_obj = PassingData(output_dir=output_dir, input_fname_basename=input_fname_basename, snp_id2writer=snp_id2writer,\
							header=header)
		cls.handleMostSignificantOperatorPerSNPPair(input_fname, cls.outputRow, param_obj)
		del snp_id2writer
	
	
	"""
	input_fname = os.path.expanduser('~/panfs/250k/boolean_gene_vs_snp_DOG1/SNPpair_41_JIC4W.tsv')
	output_dir = os.path.expanduser('~/250k/boolean_gene_vs_snp_DOG1_in_gwr/')
	outputBooleanPairIntoGWAFormat(input_fname, output_dir)
	"""
	
	@classmethod
	def getPvaluePerOperator(cls, input_fname, no_of_lines_to_skip=0):
		"""
		2009-2-19
			take boolean KW output as input_fname
			return pvalue list for each operator, no_of_lines_to_skip is the average number of lines to skip after each line to control intake.
		"""
		import csv, sys, traceback, random
		operator2pvalue_ls = {}
		reader = csv.reader(open(input_fname), delimiter='\t')
		reader.next()
		prev_snp2_id = None
		prev_snp1_id = None
		min_pvalue = None
		row_with_min_pvalue = None
		
		random_no_pool = []
		if no_of_lines_to_skip>0:
			pop_pool = range(0,no_of_lines_to_skip*2)
			for i in range(100):
				random_no_pool.append(random.sample(pop_pool, 1)[0])
		
		counter = 0
		real_counter = 0
		for row in reader:
			snp1_id, gene1_id, snp2_id, gene2_id, bool_type, pvalue, count1, count2 = row[:8]
			bool_type = int(bool_type)
			pvalue = float(pvalue)
			counter += 1
			real_counter += 1
			if bool_type not in operator2pvalue_ls:
				operator2pvalue_ls[bool_type] = []
			operator2pvalue_ls[bool_type].append(pvalue)
			if real_counter%50000==0:
				sys.stderr.write("%s%s\t\t%s"%('\x08'*40, real_counter, counter))
			if no_of_lines_to_skip>0:
				random_no = random_no_pool[counter%len(random_no_pool)]
				for i in range(random_no):
					try:
						reader.next()
						counter += 1
					except:
						traceback.print_exc()
						sys.stderr.write('%s.\n'%repr(sys.exc_info()))
						break
		del reader
		return operator2pvalue_ls
	
	@classmethod
	def stuffPvalueIntoDict(cls, row, param_obj):
		"""
		2009-2-19
			called by getPvaluePerOperatorOnlyTopInSNPPair()
		"""
		operator2pvalue_ls = getattr(param_obj, "operator2pvalue_ls")
		snp1_id, gene1_id, snp2_id, gene2_id, bool_type, pvalue, count1, count2 = row[:8]
		bool_type = int(bool_type)
		if bool_type not in operator2pvalue_ls:
			operator2pvalue_ls[bool_type] = []
		pvalue = float(pvalue)
		operator2pvalue_ls[bool_type].append(pvalue)
	
	@classmethod
	def getPvaluePerOperatorOnlyTopInSNPPair(cls, input_fname, no_of_lines_to_skip=0):
		"""
		2009-2-19
			take boolean KW output as input_fname
			return pvalue list for each operator, no_of_lines_to_skip is the average number of lines to skip after each line to control intake.
			
			call stuffPvalueIntoDict()
		"""
		import csv, sys, traceback, random
		operator2pvalue_ls = {}
		input_fname_basename = os.path.splitext(os.path.basename(input_fname))[0]
		from pymodule import PassingData
		param_obj = PassingData(operator2pvalue_ls=operator2pvalue_ls, input_fname_basename=input_fname_basename)
		cls.handleMostSignificantOperatorPerSNPPair(input_fname, cls.stuffPvalueIntoDict, param_obj)
		return operator2pvalue_ls
	
	"""
	input_fname = os.path.expanduser('~/panfs/250k/IntraGeneSNPPair_BooleanKW_m10k_call29/SNPpair_7_FT_22C.tsv')
	output_fname_prefix = os.path.expanduser('%s_top_bool_type_pvalue_hist'%os.path.splitext(input_fname)[0])
	operator_top2pvalue_ls = BooleanInteraction.getPvaluePerOperatorOnlyTopInSNPPair(input_fname, no_of_lines_to_skip=0)
	BooleanInteraction.drawPvalueHist(operator_top2pvalue_ls, output_fname_prefix, minus_log_pvalue=True)
	"""
	
	@classmethod
	def drawPvalueHist(cls, operator2pvalue_ls, output_fname_prefix, minus_log_pvalue=False):
		"""
		2009-2-19
			draw histogram distribution for each each operator's pvalue list
		"""
		import pylab
		pylab.clf()
		from variation.src.PlotGenePairAssoResult import PlotGenePairAssoResult
		#bool_type2marker_and_name
		import random
		hist_handler_ls = []
		legend_ls = []
		bool_type_ls = operator2pvalue_ls.keys()
		bool_type_ls.sort()
		color_ls = ['r', 'g', 'c', 'b', 'y', 'k']
		log_func = lambda x: -math.log(x)
		for i in range(len(bool_type_ls)):
			pylab.clf()
			bool_type = bool_type_ls[i]
			pvalue_ls = operator2pvalue_ls[bool_type]
			
			pvalue_ls = random.sample(pvalue_ls, 10000)
			
			if minus_log_pvalue:
				xlabel = '-logPvalue'
				pvalue_ls = map(log_func, pvalue_ls)
			else:
				xlabel = 'pvalue'
			n1 = pylab.hist(pvalue_ls, 50, alpha=0.3, normed=1, facecolor=color_ls[i])
			hist_handler_ls.append(n1[2][0])
			marker, bool_type_name = PlotGenePairAssoResult.bool_type2marker_and_name[bool_type][:2]
			legend_ls.append(bool_type_name)
			#pylab.legend(hist_handler_ls, legend_ls)
			pylab.ylabel('frequency')
			
			pylab.title(bool_type_name)
			pylab.savefig('%s_%s_%s.png'%(output_fname_prefix, bool_type, bool_type_name), dpi=300)
	
	"""
	input_fname = os.path.expanduser('~/panfs/250k/IntraGeneSNPPair_BooleanKW_m10k_call29/SNPpair_7_FT_22C.tsv')
	output_fname_prefix = os.path.expanduser('%s_bool_type_pvalue_hist'%os.path.splitext(input_fname)[0])
	operator2pvalue_ls = BooleanInteraction.getPvaluePerOperator(input_fname, no_of_lines_to_skip=0)
	BooleanInteraction.drawPvalueHist(operator2pvalue_ls, output_fname_prefix, minus_log_pvalue=True)
	"""
	
	@classmethod
	def countBoolTypeTopForEachSNPPair(cls, row, param_obj):
		"""
		2009-2-19
		"""
		import math, sys, traceback
		operator_int_pvalue2count = getattr(param_obj, "operator_int_pvalue2count")
		int_pvalue2count = getattr(param_obj, "int_pvalue2count")
		snp1_id, gene1_id, snp2_id, gene2_id, bool_type, pvalue, count1, count2 = row[:8]
		bool_type = int(bool_type)
		try:
			minus_log_pvalue = -math.log(float(pvalue))
		except:
			traceback.print_exc()
			sys.stderr.write('%s. minus_log_pvalue=50.\n'%repr(sys.exc_info()))
			minus_log_pvalue = 50
		int_pvalue = int(minus_log_pvalue)
		if int_pvalue not in int_pvalue2count:
			int_pvalue2count[int_pvalue] = 0
		int_pvalue2count[int_pvalue] += 1
		
		operator_int_pvalue = (bool_type, int_pvalue)
		if operator_int_pvalue not in operator_int_pvalue2count:
			operator_int_pvalue2count[operator_int_pvalue] = 0
		operator_int_pvalue2count[operator_int_pvalue] += 1
	
	@classmethod
	def plot_operator_int_pvalue2count(cls, operator_int_pvalue2count, int_pvalue2count, output_fname_prefix, \
									bool_type_ls=[1,2,4,6,7], int_pvalue_step=2):
		"""
		2009-2-19
		"""
		#from variation.src.PlotGenePairAssoResult import PlotGenePairAssoResult
		#PlotGenePairAssoResult.bool_type2marker_and_name[bool_type]
		int_pvalue_ls = int_pvalue2count.keys()
		int_pvalue_ls.sort()
		
		#group int_pvalue into groups according to int_pvalue_step
		row_id_ls = []
		for i in range(len(int_pvalue_ls)):
			int_pvalue = int_pvalue_ls[i]
			row_id = int_pvalue/int_pvalue_step
			if len(row_id_ls)==0:
				row_id_ls.append(row_id)
			else:
				if row_id!=row_id_ls[-1]:	#only add when it's new
					row_id_ls.append(row_id)
		
		col_id_ls = bool_type_ls
		import numpy
		#col_id_ls = numpy.array(col_id_ls)
		no_of_rows = len(row_id_ls)
		no_of_cols = len(col_id_ls)
		data_matrix = numpy.zeros([no_of_rows, no_of_cols], numpy.int)
		from pymodule import SNPData
		countData = SNPData(row_id_ls=row_id_ls, col_id_ls=col_id_ls, data_matrix=data_matrix)
		
		for int_pvalue in int_pvalue_ls:
			row_id = int_pvalue/int_pvalue_step
			row_index = countData.row_id2row_index[row_id]
			for bool_type in col_id_ls:
				operator_int_pvalue = (bool_type, int_pvalue)
				if operator_int_pvalue in operator_int_pvalue2count:
					count = operator_int_pvalue2count[operator_int_pvalue]
				else:
					count = 0
				col_index = countData.col_id2col_index[bool_type]
				countData.data_matrix[row_index][col_index] += count
		
		perc_matrix = numpy.zeros([no_of_rows, no_of_cols], numpy.float)
		percData = SNPData(row_id_ls=row_id_ls, col_id_ls=col_id_ls, data_matrix=perc_matrix)
		row_sum_array = numpy.sum(countData.data_matrix, 1, numpy.float)
		for j in range(no_of_cols):
			percData.data_matrix[:,j] = numpy.sum(countData.data_matrix[:,range(j+1)], 1)/row_sum_array	#cumulative percentage
		
		import pylab
		pylab.clf()
		pylab.grid(True, alpha=0.3)
		from variation.src.PlotGenePairAssoResult import PlotGenePairAssoResult
		color_ls = ['r', 'g', 'c', 'b', 'y', 'k']
		plot_ls = []
		legend_ls = []
		row_id_ls = numpy.array(row_id_ls)
		for j in range(no_of_cols):
			bool_type = percData.col_id_ls[j]
			marker, bool_type_name = PlotGenePairAssoResult.bool_type2marker_and_name[bool_type][:2]
			x_ls = percData.data_matrix[:,j]
			y_ls = row_id_ls*int_pvalue_step
			a = pylab.plot(x_ls, y_ls, '.-', color=color_ls[j])
			plot_ls.append(a)
			legend_ls.append(bool_type_name)
		
		pylab.legend(plot_ls, legend_ls)
		pylab.ylabel('-logPvalue')
		pylab.xlabel('percentage')
		pylab.title('Percentage of each bool type at certain log pvalue range')
		pylab.savefig('%s_perc_bool_by_log_pvalue.png'%(output_fname_prefix), dpi=300)
			
	
	"""
operator_int_pvalue2count = {}
int_pvalue2count = {}
from pymodule import PassingData
param_obj = PassingData(operator_int_pvalue2count=operator_int_pvalue2count, int_pvalue2count=int_pvalue2count)
BooleanInteraction.handleMostSignificantOperatorPerSNPPair(input_fname, BooleanInteraction.countBoolTypeTopForEachSNPPair, param_obj)
BooleanInteraction.plot_operator_int_pvalue2count(operator_int_pvalue2count, int_pvalue2count, output_fname_prefix)
	"""


class FRIDeletion(object):
	"""
		in processing FRI deletion data from Shindo2005. check plone doc, /research/variation/log-2008-07.
	2008-08-05
	"""
	@classmethod
	def outputShindo2005(cls, input_fname, output_fname, which_type_of_id_to_output=1):
		"""
		2008-08-06
			in output. if which_type_of_id_to_output==1, output in ecotypeid; else output in accession_id
		"""
		from variation.src.AtDB import AtDB, Sequence, Alignment
		db = AtDB(hostname='localhost')
		from pymodule import read_data, write_data_matrix
		header, strain_acc_list, category_list, data_matrix = read_data(input_fname, turn_into_integer=0)
		
		no_of_rows = len(strain_acc_list)
		import numpy
		new_data_matrix = numpy.ones([no_of_rows, 2], numpy.int8)
		new_strain_acc_list = []
		new_category_list = []
		for i in range(len(strain_acc_list)):
			strain_acc = strain_acc_list[i].upper()
			#2008-08-06 for magnus's data.csv	no swedish letters.
			rows = db.metadata.bind.execute("select a2e.ecotype_id, a2e.nativename, a2e.accession_id from accession2tg_ecotypeid a2e, magnus_192_vs_accession m2a where upper(m2a.linename)='%s' and m2a.accession_id=a2e.accession_id"%\
									(strain_acc))
			
			#2008-08-06 for stuff extracted from Supplemental Figure 4A, Figure 4B. this query doesn't help to recognize those swedish letters
			#rows = db.metadata.bind.execute("select a2e.ecotype_id, a2e.nativename, a2e.accession_id from accession2tg_ecotypeid a2e where upper(a2e.accession_name)='%s'"%\
			#						(strain_acc))
			try:
				row = rows.fetchone()
				if which_type_of_id_to_output==1:
					new_strain_acc_list.append(row.ecotype_id)
				else:
					new_strain_acc_list.append(row.accession_id)
				new_category_list.append(row.nativename)
				deletion_code = int(category_list[i])
				deletion_code_index = deletion_code-2
				if deletion_code_index>=0:
					new_data_matrix[i][deletion_code_index] = 2	#allele 2
			except:
				print i, strain_acc
				new_strain_acc_list.append(strain_acc)
				new_category_list.append(strain_acc)
			
		new_header = header[:2] + ['4_268809_0', '4_269962_8']
		write_data_matrix(new_data_matrix, output_fname, new_header, new_strain_acc_list, new_category_list)
	
	"""
	input_fname = os.path.expanduser('~/script/variation/doc/FRI/Shindo2005_data.csv')
	output_fname = os.path.expanduser('~/script/variation/doc/FRI/Shindo2005_data_SNP.tsv')
	outputShindo2005(input_fname, output_fname)
	"""


class Data149_Haplotype(object):
	"""
	2008-09-10
		Alex sent me some files of 149SNP haplotype consensus seq and country info to get pairwise picture. convert it to my format.
	"""
	
	def get_haplotype_group_name2country_ls(file_country):
		import csv, os, sys
		sys.stderr.write("Reading in country info of haplotypes ...")
		reader = csv.reader(open(file_country), delimiter='\t')
		reader.next()
		haplotype_group_name2country_ls = {}
		for row in reader:
			g_name = row[0][:-1]
			if g_name in haplotype_group_name2country_ls:
				sys.stderr.write("Error: %s already existed in haplotype_group_name2country_ls.\n"%g_name)
			haplotype_group_name2country_ls[g_name] = row[1:]
		del reader
		sys.stderr.write("Done.\n")
		return haplotype_group_name2country_ls
		
	def processAlexHaplotypeFile(file_country, file2, output_fname):
		haplotype_group_name2country_ls = get_haplotype_group_name2country_ls(file_country)
		import csv, os, sys
		sys.stderr.write("Reading in consensus of haplotypes ...")
		from pymodule import nt2number
		reader = csv.reader(open(file2), delimiter='\t')
		reader.next()
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		i = 0
		data_row = []
		for row in reader:
			g_name = row[0]
			if g_name in haplotype_group_name2country_ls:	#2008-09-12 has to be in the country file
				country_str = '_'.join(haplotype_group_name2country_ls[g_name])
				label = '%s_%s_%s'%(g_name, country_str, row[-1])
				data_row.append(label)
				data_row.append(label)
				for j in range(1, len(row)-2):
					for k in range(len(row[j])):
						data_row.append(nt2number[row[j][k]])
				if i==0:	#need to write a header
					writer.writerow(['','']+range(1, len(data_row)-1))
				writer.writerow(data_row)
				data_row = []
				i += 1
		sys.stderr.write("Done.\n")
		del reader, writer
	
	
	def getMatrixOutofCrossMatchResult(curs, cross_match_outfile, file_country, output_fname, max_mismatch_rate=1):
		haplotype_group_name2country_ls = get_haplotype_group_name2country_ls(file_country)
		import csv, os, sys
		longitude_g_name_ls = []
		for g_name, country_ls in haplotype_group_name2country_ls.iteritems():
			if len(country_ls)==1:
				curs.execute("select latitude, longitude, abbr from stock.country where abbr='%s'"%(country_ls[0]))
				rows = curs.fetchall()
				longitude = rows[0][1]
				longitude_g_name_ls.append((longitude, country_ls[0], g_name))
			else:
				longitude_g_name_ls.append((-360, '', g_name))
		
		haplotype_group_name2index = {}
		longitude_g_name_ls.sort()
		prev_country = None
		no_of_countries = 0
		label_ls = []
		for i in range(len(longitude_g_name_ls)):
			country = longitude_g_name_ls[i][1]
			if prev_country==None:
				prev_country = country
			elif country!=prev_country:	#insert a separator
				no_of_countries += 1
				g_name = '-%s'%no_of_countries
				haplotype_group_name2index[g_name] = len(haplotype_group_name2index)
				label_ls.append('')
				prev_country = country
			g_name = longitude_g_name_ls[i][-1]
			haplotype_group_name2index[g_name] = len(haplotype_group_name2index)
			label_ls.append(g_name)	#change it to a fuller one later
		
		sys.stderr.write("Getting data matrix ... \n")
		import numpy
		data_matrix = numpy.zeros([len(haplotype_group_name2index), len(haplotype_group_name2index)], numpy.float)
		data_matrix[:,:] = -1	#mark everything as NA
		reader = csv.reader(open(cross_match_outfile), delimiter='\t')
		#figure out which variable is in which column
		header = reader.next()
		col_name2index = {}
		for i in range(len(header)):
			column_name = header[i]
			col_name2index[column_name] = i
		
		for row in reader:
			source_id = row[col_name2index['source_id']]
			source_g_name = source_id.split('_')[0]
			target_id = row[col_name2index['target_id']]
			target_g_name = target_id.split('_')[0]
			mismatch_rate = float(row[col_name2index['mismatch_rate']])
			no_of_mismatches = int(row[col_name2index['no_of_mismatches']])
			no_of_non_NA_pairs = int(row[col_name2index['no_of_non_NA_pairs']])
			row_index = haplotype_group_name2index[source_g_name]
			col_index = haplotype_group_name2index[target_g_name]
			label_ls[row_index] = source_id
			if no_of_non_NA_pairs>=20 and mismatch_rate<=max_mismatch_rate:
				data_matrix[row_index][col_index] = data_matrix[col_index][row_index] = mismatch_rate
		
		for g_name, boundary_index in haplotype_group_name2index.iteritems():
			if g_name[0]=='-':
				data_matrix[boundary_index,:] = -3
				data_matrix[:,boundary_index] = -3
		sys.stderr.write("Done.\n")
		from pymodule import write_data_matrix
		write_data_matrix(data_matrix, output_fname, ['','']+label_ls, label_ls, label_ls)
			
	"""
file_country = os.path.expanduser('~/script/variation/data/149SNP/haps.countries.tsv')
file2 = os.path.expanduser('~/script/variation/data/149SNP/consensus_seqs.HG(0.005.S).tsv')
output_fname = os.path.expanduser('~/script/variation/data/149SNP/haps_in_149.tsv')
processAlexHaplotypeFile(file_country, file2, output_fname)

file_country = os.path.expanduser('~/panfs/149CrossMatch/haps.countries.tsv')
cross_match_outfile = os.path.expanduser('~/panfs/149CrossMatch/alex_hap_cross_match.tsv')
output_fname = os.path.expanduser('~/panfs/149CrossMatch/alex_hap_cross_match_matrix.tsv')
max_mismatch_rate=1
getMatrixOutofCrossMatchResult(curs, cross_match_outfile, file_country, output_fname, max_mismatch_rate=1)

#2008-09-12
file_country = os.path.expanduser('~/script/variation/data/149SNP/haps.countries.x11.tsv')
output_fname = os.path.expanduser('~/script/variation/data/149SNP/haps_in_149_no_4_bad_plates.tsv')
processAlexHaplotypeFile(file_country, file2, output_fname)

file_country = os.path.expanduser('~/banyan_home/script/variation/data/149SNP/haps.countries.x11.tsv')
cross_match_outfile = os.path.expanduser('~/panfs/149CrossMatch/alex_hap_no_4_bad_plates_cross_match.tsv')
output_fname = os.path.expanduser('~/panfs/149CrossMatch/alex_hap_no_4_bad_plates_cross_match_matrix_a0.3.tsv')
getMatrixOutofCrossMatchResult(curs, cross_match_outfile, file_country, output_fname, max_mismatch_rate=0.3)
	"""
	
	"""
	2008-09-11
		found out sequenom group 149SNPs into 4 blocks (38, 37, 37, 37) and then genotype each block on their plate separately
		the SNPs are not in chromosome,position order. It's in id of table snps order.
		now partition the data file into 4 blocks accordingly.
	"""
	def partition149SNPDataInto4Blocks(input_fname, db_149, output_fname_prefix):
		from pymodule import read_data, SNPData, write_data_matrix
		header, strain_acc_list, category_list, data_matrix = read_data(input_fname)
		snpData = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list, data_matrix=data_matrix, turn_into_array=1)
		
		block1_snp_id_ls = []
		block2_snp_id_ls = []
		block3_snp_id_ls = []
		block4_snp_id_ls = []
		sys.stderr.write("Grouping SNPs ...")
		rows = db_149.metadata.bind.execute("select * from snps order by id")
		for row in rows:
			snp_id = '%s_%s'%(row.chromosome, row.position)
			if row.id<=38:
				block1_snp_id_ls.append(snp_id)
			elif row.id>38 and row.id<=75:
				block2_snp_id_ls.append(snp_id)
			elif row.id>75 and row.id<=112:
				block3_snp_id_ls.append(snp_id)
			else:
				block4_snp_id_ls.append(snp_id)
		sys.stderr.write("Done.\n")
		sys.stderr.write("Splitting matrix ... \n")
		no_of_rows = snpData.data_matrix.shape[0]
		import numpy
		block_snp_id_ls_ls = [block1_snp_id_ls, block2_snp_id_ls, block3_snp_id_ls, block4_snp_id_ls]
		for i in range(len(block_snp_id_ls_ls)):
			output_fname = '%s_%s.tsv'%(output_fname_prefix, i)
			no_of_snps_in_this_block = len(block_snp_id_ls_ls[i])
			d_matrix = numpy.zeros([no_of_rows, no_of_snps_in_this_block], numpy.int8)
			for j in range(no_of_snps_in_this_block):
				snp_id = block_snp_id_ls_ls[i][j]
				col_index = snpData.col_id2col_index[snp_id]
				d_matrix[:,j] = snpData.data_matrix[:,col_index]
			header = ['', ''] + block_snp_id_ls_ls[i]
			write_data_matrix(d_matrix, output_fname, header, strain_acc_list, category_list)
		sys.stderr.write("Done.\n")
	"""
input_fname = os.path.expanduser('~/panfs/149CrossMatch/stock_149SNP_y0000110101.tsv')
output_fname_prefix = os.path.expanduser('~/panfs/149CrossMatch/149SNPSequenomBlock')
partition149SNPDataInto4Blocks(input_fname, db_149, output_fname_prefix)
	"""
	
	def getDistanceMatrixOutofCrossMatchFile(input_fname, max_no_of_strains, strainid2index, strainid_ls, min_no_of_non_NA_pairs=10):
		import os, sys,csv
		sys.stderr.write("Getting distance matrix from %s ... "%input_fname)
		reader = csv.reader(open(input_fname), delimiter='\t')
		header = reader.next()
		col_name2index = {}
		for i in range(len(header)):
			column_name = header[i]
			col_name2index[column_name] = i
		import numpy
		data_matrix = numpy.zeros([max_no_of_strains, max_no_of_strains], numpy.float)
		data_matrix[:] = -1	#set default to NA
		i = 0
		for row in reader:
			strainid = int(row[col_name2index['strainid']])
			target_id = int(row[col_name2index['target_id']])
			mismatch_rate = float(row[col_name2index['mismatch_rate']])
			no_of_mismatches = int(row[col_name2index['no_of_mismatches']])
			no_of_non_NA_pairs = int(row[col_name2index['no_of_non_NA_pairs']])
			if no_of_non_NA_pairs>=min_no_of_non_NA_pairs:
				if strainid not in strainid2index:
					strainid2index[strainid] = len(strainid2index)
					strainid_ls.append(strainid)
				if target_id not in strainid2index:
					strainid2index[target_id] = len(strainid2index)
					strainid_ls.append(target_id)
				row_index = strainid2index[strainid]
				col_index = strainid2index[target_id]
				if row_index<max_no_of_strains and col_index<max_no_of_strains:
					data_matrix[row_index][col_index] = mismatch_rate
					data_matrix[col_index][row_index] = mismatch_rate
				else:
					sys.stderr.write("no of strains in this data exceeds max_no_of_strains %s.\n"%(max_no_of_strains))
		del reader
		sys.stderr.write("%s strains. Done.\n"%(len(strainid_ls)))
		return data_matrix
	
	def findStrainShowMsmatchRateVariationIn4Blocks(block_fname_ls, output_fname, max_no_of_strains=7000, min_no_of_non_NA_pairs=10, min_var=0.03):
		import os, sys,csv
		from pymodule import PassingData
		no_of_blocks = len(block_fname_ls)
		strainid2index = {}
		strainid_ls = []
		data_matrix_ls = []
		for i in range(no_of_blocks):
			data_matrix = getDistanceMatrixOutofCrossMatchFile(block_fname_ls[i], max_no_of_strains, strainid2index, strainid_ls, min_no_of_non_NA_pairs)
			data_matrix_ls.append(data_matrix)
		
		sys.stderr.write("Looking for pairs having mismatch_rate var >= %s ..."%(min_var))
		import rpy
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		no_of_strains = len(strainid_ls)
		for i in range(no_of_strains):
			for j in range(i+1, no_of_strains):
				mismatch_ls = []
				mismatch_non_NA_ls = []
				for data_matrix in data_matrix_ls:
					if data_matrix[i][j]!=-1:
						mismatch_non_NA_ls.append(data_matrix[i][j])
					mismatch_ls.append(data_matrix[i][j])
				mismatch_var = rpy.r.var(mismatch_non_NA_ls)
				if mismatch_var>=min_var:
					writer.writerow([strainid_ls[i], strainid_ls[j], mismatch_var]+mismatch_ls)
		del writer
	
	"""
block_fname_ls = []
for i in range(4):
	block_fname_ls.append(os.path.expanduser('~/panfs/149CrossMatch/149SNPSequenomBlock_%s_cross_match.tsv'%(i)))

output_fname = os.path.expanduser('~/panfs/149CrossMatch/149SNPSequenomBlock_high_var.tsv')
findStrainShowMsmatchRateVariationIn4Blocks(block_fname_ls, output_fname, max_no_of_strains=7000, min_no_of_non_NA_pairs=10, min_var=0.03)
	"""
	
class CmpDifferentData(object):
	"""
	2009-2-17 like QC
	"""
	
	"""
	2008-10-11
		window average of SNP mismatch rates.
		input_fname is output of TwoSNPData.output_col_id2NA_mismatch_rate_InGWRFormat() (invoked in Qc.py)
	"""
	@classmethod
	def average_SNP_mismatch_rate_within_window(cls, input_fname, output_fname, window_size=100000):
		"""
		"""
		from pymodule import figureOutDelimiter
		import csv
		delimiter = figureOutDelimiter(input_fname)
		reader = csv.reader(open(input_fname), delimiter=delimiter)
		
		chr_pos2mismatch_rate = {}
		for row in reader:
			chr, pos, mismatch_rate = row
			chr = int(chr)
			pos = int(pos)
			mismatch_rate = float(mismatch_rate)
			chr_pos2mismatch_rate[(chr, pos)] = mismatch_rate
		
		prev_chr_pos = None
		start_pos = None
		no_of_mismatches = 0.
		no_of_snps = 0
		chr_pos_ls = chr_pos2mismatch_rate.keys()
		chr_pos_ls.sort()
		chr_start_stop2mismatch_rate_ls = {}
		for chr_pos in chr_pos_ls:
			chr, pos = chr_pos
			if prev_chr_pos==None:
				prev_chr_pos = chr_pos
			if start_pos == None:
				start_pos = pos
			
			if chr==prev_chr_pos[0]:
				if pos-start_pos<=window_size:
					no_of_mismatches += chr_pos2mismatch_rate[chr_pos]
					no_of_snps += 1
				else:	#out of window
					chr_start_stop2mismatch_rate_ls[(prev_chr_pos[0], start_pos, prev_chr_pos[1])] = [no_of_mismatches/no_of_snps, no_of_snps]
					
					#reset
					start_pos = pos
					no_of_mismatches = chr_pos2mismatch_rate[chr_pos]
					no_of_snps = 1
			elif chr!=prev_chr_pos[0]:
				chr_start_stop2mismatch_rate_ls[(prev_chr_pos[0], start_pos, prev_chr_pos[1])] = [no_of_mismatches/no_of_snps, no_of_snps]
				
				
				#reset
				start_pos = pos
				no_of_mismatches = chr_pos2mismatch_rate[chr_pos]
				no_of_snps = 1
			
			prev_chr_pos=chr_pos
		
		writer = csv.writer(open(output_fname, 'w'), delimiter=delimiter)
		chr_start_stop_ls = chr_start_stop2mismatch_rate_ls.keys()
		chr_start_stop_ls.sort()
		for chr_start_stop in chr_start_stop_ls:
			mismatch_rate_ls = chr_start_stop2mismatch_rate_ls[chr_start_stop]
			chr, start_pos, stop_pos = chr_start_stop
			writer.writerow([chr, start_pos, mismatch_rate_ls[0], stop_pos, mismatch_rate_ls[1]])
		del writer
	
	"""
	input_fname = '/tmp/chromosmal_SNP_mismatch_139'
	output_fname = '%s_window_100k'%(input_fname)
	average_SNP_mismatch_rate_within_window(input_fname, output_fname, window_size=100000)
	"""


class DB250k(object):
	@classmethod
	def dumpFTGene2File(cls, db, output_fname, list_type_id=28):
		"""
		2008-11-25 dump all flowering time genes from table ft_gene into file. (for MpiInterGeneSNPPairAsso.py)
		"""
		#rows = db.metadata.bind.execute("select distinct f.gene_id, g.gene_symbol from genome.gene g, ft_gene f where f.gene_id=g.gene_id")
		rows = db.metadata.bind.execute("select distinct l.gene_id, g.gene_symbol from genome.gene g, candidate_gene_list l where l.gene_id=g.gene_id and l.list_type_id=%s order by gene_symbol"%(list_type_id))
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		for row in rows:
			writer.writerow([row.gene_id, row.gene_symbol])
		del writer
	
	"""
output_fname = '/tmp/ft_gene.tsv'
dumpFTGene2File(db_250k, output_fname)
	"""
	
	@classmethod
	def outputSNPCandidateGeneAssociation(cls, snps_context_picklef, list_type_id, output_fname):
		"""
		2008-11-13
		output SNP-gene association (from a file containing the pickled snps_context_wrapper) into a file
		for magnus
		"""
		import cPickle, csv
		from GeneListRankTest import GeneListRankTest
		candidate_gene_set = GeneListRankTest.dealWithCandidateGeneList(list_type_id, return_set=True)
		from Stock_250kDB import SnpsContext
		
		picklef = open(snps_context_picklef)
		snps_context_wrapper = cPickle.load(picklef)
		del picklef
		
		chr_pos_ls = snps_context_wrapper.chrpos2snps_id.keys()
		chr_pos_ls.sort()
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = ['gene_id', 'snps_id', 'chromosome', 'position', 'disp_pos', 'left_or_right', 'disp_pos_comment']
		writer.writerow(header)
		for chr, pos in chr_pos_ls:
			snps_context_matrix = snps_context_wrapper.returnGeneLs(chr, pos)
			assign_snp_candidate_gene = 0
			assign_snp_non_candidate_gene = 0
			for snps_context in snps_context_matrix:
				snps_id, disp_pos, gene_id = snps_context
				if gene_id in candidate_gene_set:
					snps_context = SnpsContext.query.filter_by(snps_id=snps_id).filter_by(gene_id=gene_id).first()
					left_or_right = getattr(snps_context, 'left_or_right', '')
					disp_pos_comment = getattr(snps_context, 'disp_pos_comment', '')
					writer.writerow([gene_id, snps_id, chr, pos, disp_pos, left_or_right, disp_pos_comment])
		del writer
	
	"""
snps_context_picklef = './mnt2/panfs/250k/snps_context_g0_m500'
list_type_id = 28
output_fname = '/tmp/snps_context_g0_m500_list28.tsv'
outputSNPCandidateGeneAssociation(snps_context_picklef, list_type_id, output_fname)
	"""
	
	
	"""
	2008-04-11
		this is a one-time 250k pipeline fix to avoid memory-hefty intensity re-output.
		
		form 2 lists of intensity matrix filenames based on the new and old array_info_table
		output them as two-column (old_fname, new_fname) into a output_fname.
		shell/file_batch_move.py reads the output_fname and handle name changing.
	
	2008-04-11
		turns out to be useless because the header in each intensity matrix file has array_id embedded.
		have to output each array into intensity matrix.
	"""
	@classmethod
	def output_intensity_fname(cls, curs, new_array_info_table, old_array_info_table, output_fname):
		"""
		"""
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		curs.execute("select n.id, o.id from %s n, %s o where n.original_filename=o.filename"%(new_array_info_table, old_array_info_table))
		rows = curs.fetchall()
		for row in rows:
			new_array_id, old_array_id = row
			new_fname = '%s_array_intensity.tsv'%new_array_id
			old_fname = '%s_array_intensity.tsv'%old_array_id
			writer.writerow([old_fname, new_fname])
		del writer
	
	"""
	new_array_info_table='stock_250k.array_info'
	old_array_info_table='stock_250k.array_info_2008_04_11'
	output_fname = '/tmp/intensity_fname.rename.tsv'
	output_intensity_fname(curs, new_array_info_table, old_array_info_table, output_fname)
	"""
	
	
	"""
	2008-05-31
		output results from stock_250k.results, totally db based, not reading results.filename
	"""
	@classmethod
	def outputResults(cls, db, results_method_id, output_fname):
		import csv
		import sqlalchemy as sql
		
		conn = db.connection	#establish the connection before referring db.tables (in case it hasn't been setup)
		sys.stderr.write("Getting marker_id2position ... ")
		marker_id2position = {}
		snps_table = db.tables['snps'].alias()
		results = conn.execute(sql.select([snps_table.c.id, snps_table.c.chromosome, snps_table.c.position, snps_table.c.end_position]))
		for row in results:
			marker_id2position[row.id] = (row.chromosome, row.position, row.end_position)
		del results
		sys.stderr.write("Done.\n")
		
		sys.stderr.write("Outputting results_method_id=%s ... "%results_method_id)
		results_table = db.tables['results'].alias()
		results = conn.execute(sql.select([results_table.c.snps_id, results_table.c.score], results_table.c.results_method_id==results_method_id,\
										  order_by=[results_table.c.snps_id]))
		
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		writer.writerow(['chromosome', 'position', 'score'])
		for row in results:
			chromosome, position, end_position = marker_id2position[row.snps_id]
			writer.writerow([chromosome, position, row.score])
		del writer
		sys.stderr.write("Done.\n")
	
	"""
	from variation.src.db import Stock_250kDatabase
	db = Stock_250kDatabase(username='nordborglab',
					   password='papaya', hostname='papaya.usc.edu', database='stock_250k')
	conn = db.connection	#establish the connection before referring db.tables (it needs to be setup)
	outputResults(db, 5, '/tmp/5_results.tsv')
	outputResults(db, 6, '/tmp/6_results.tsv')
	"""
	
	"""
	2008-07-21
	"""
	@classmethod
	def drawScoreHistogram(cls, curs, results_method_id, list_type_id=1, do_log10_transformation=True, min_or_max_func='min'):
		from Stock_250kDB import ResultsMethod, GeneList, ResultsByGene, CandidateGeneRankSumTestResult
		rm = ResultsMethod.get(results_method_id)
		db_rows = GeneList.query.filter_by(list_type_id=list_type_id)
		from sets import Set
		gene_set = Set()
		for gl in db_rows:
			gene_set.add(gl.gene_id)
		
		import math
		score_ls1 = []
		score_ls2 = []
		#db_rows = ResultsByGene.query.filter_by(results_method_id=results_method_id)
		curs.execute("select r.gene_id, %s(r.score) as score from results_by_gene r  where r.results_method_id=%s group by r.gene_id"%(min_or_max_func, results_method_id))
		db_rows = curs.fetchall()
		for row in db_rows:
			gene_id, score = row
			if do_log10_transformation:
				score = -math.log(score)
			if gene_id in gene_set:
				score_ls1.append(score)
			else:
				score_ls2.append(score)
		import pylab
		pylab.clf()
		pylab.title('results_method_id=%s, (%s on %s) by list type: %s.'%(rm.id, rm.analysis_method.short_name, rm.phenotype_method.short_name, gl.list_type.short_name))
		n1 = pylab.hist(score_ls1, 100, alpha=0.4, normed=1)
		n2 = pylab.hist(score_ls2, 100, alpha=0.4, normed=1, facecolor='r')
		pylab.legend([n1[2][0], n2[2][0]], ['candidate gene list', 'non-candidate gene list'])
		pylab.show()
		return score_ls1, score_ls2
	
	"""
	score_ls1, score_ls2 = drawScoreHistogram(curs, 23, 1)
	"""
	
	@classmethod
	def plotArrayMismatchVsMedianIntensity(cls, db, array_id_ls, take_log=False):
		"""
		2009-3-11		
			plot array mismatch rate vs median intensity of all probes
			
			array_id_ls is comma or dash-separated list of array ids
		"""
		import math
		if array_id_ls:
			from pymodule import getListOutOfStr
			array_id_ls = getListOutOfStr(array_id_ls, data_type=str)
		else:
			array_id_ls = range(500)
			array_id_ls = map(str, array_id_ls)
		rows = db.metadata.bind.execute("select mismatch_rate, median_intensity from stock_250k.view_qc where array_id in (%s)"%(','.join(array_id_ls)))
		y_ls = []
		x_ls = []
		for row in rows:
			if take_log:
				median_intensity = math.log10(row.median_intensity)
			else:
				median_intensity = row.median_intensity
			x_ls.append(median_intensity)
			y_ls.append(row.mismatch_rate)
		
		import pylab
		pylab.plot(x_ls, y_ls, '.', alpha=0.5)
		pylab.title('mismatch rate vs median_intensity of arrays')
		pylab.xlabel('median_intensity')
		pylab.ylabel('mismatch rate')
		pylab.show()
	
	"""
plotArrayMismatchVsMedianIntensity(db_250k, '616-1045')
	"""
	
	@classmethod
	def saveArrayQCIntoTableFile(cls, db, array_id_ls, output_fname, take_log=False):
		"""
		2009-3-13		
			output array QC information in a matrix format which pymodule/DataMatrixGuiXYProbe.py can recognize and do clickable scatterplots 
		"""
		import math
		if array_id_ls:
			from pymodule import getListOutOfStr
			array_id_ls = getListOutOfStr(array_id_ls, data_type=str)
		else:
			array_id_ls = range(500)
			array_id_ls = map(str, array_id_ls)
		rows = db.metadata.bind.execute("select * from stock_250k.view_qc where array_id in (%s)"%(','.join(array_id_ls)))
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = ['nativename', 'original_array_name', 'ecotype_id', 'array_id', 'mismatch_rate', 'median_intensity', \
					'QC_NA_rate', 'no_of_mismatches', 'no_of_non_NA_pairs', 'qc_method_id', 'call_info_id', 'array_original_filename']
		writer.writerow(header)
		for row in rows:
			if take_log:
				median_intensity = math.log10(row.median_intensity)
			else:
				median_intensity = row.median_intensity
			original_array_name = os.path.basename(row.array_original_filename)
			output_row = []
			for col_name in header:
				if col_name=='original_array_name':
					output_row.append(original_array_name)
				elif col_name=='median_intensity':
					output_row.append(median_intensity)
				else:
					output_row.append(getattr(row, col_name))
			writer.writerow(output_row)
		del writer
		
	"""
output_fname = '/tmp/CHLA_2009_01_QC.tsv'
saveArrayQCIntoTableFile(db_250k, '616-1045', output_fname)
	"""

	"""
	2009-4-13 function to check all 250k calls to see which is bad , which is good.
		for bad ones, further check cross_match results to see if mis-labelling happens
		
		check view_qc
		check qc_cross_match_table to see if any cross-labeling
			no_of_non_NA_pairs>=40
			mismatch_rate<2%
		
		quality meaning:
			0: bad
			1: good
			2: mis-labelled
			
	"""
	@classmethod
	def reportBadArrays(cls, db, call_method_id, qc_method_id, output_fname, max_mismatch_rate=0.1, min_no_of_non_NA_pairs=40, \
					max_mislabel_mismatch_rate=0.02,\
					view_qc_table='view_qc', qc_cross_match_table='qc_cross_match'):
		from variation.src.common import get_ecotypeid2nativename
		ecotypeid2nativename = get_ecotypeid2nativename(db.metadata.bind, ecotype_table='stock.ecotype')
		sys.stderr.write("Reporting 250k arrays ... \n")
		rows = db.metadata.bind.execute("select * from %s where qc_method_id=%s and call_method_id=%s order by nativename"%\
									(view_qc_table, qc_method_id, call_method_id))
		
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header_row = ['ecotype_id', 'nativename', 'array_original_filename', 'array_created', 'quality', 'true ecotype_id', 'true nativename']
		writer.writerow(header_row)
		results = []
		counter = 0
		no_of_bad_ones = 0
		no_of_mis_labelled_ones = 0
		
		good_ecotype_id_set = set()
		bad_ecotype_id_set = set()
		true_ecotype_id_set = set()
		
		for row in rows:
			quality = 0
			counter += 1
			true_ecotype_id_ls = []
			true_nativename_ls = []
			if row.no_of_non_NA_pairs>=min_no_of_non_NA_pairs and row.mismatch_rate<=max_mismatch_rate:
				quality = 1
			else:
				no_of_bad_ones += 1
				cross_match_rows = db.metadata.bind.execute("select * from %s where call_info_id=%s and no_of_non_NA_pairs>=%s order by mismatch_rate"%\
										(qc_cross_match_table, row.call_info_id, min_no_of_non_NA_pairs))
				for cross_match_row in cross_match_rows:
					if cross_match_row.mismatch_rate<=max_mislabel_mismatch_rate:
						if cross_match_row.vs_ecotype_id in ecotypeid2nativename:	#make sure ecotypeid is still alive
							true_ecotype_id_ls.append(cross_match_row.vs_ecotype_id)
							true_nativename_ls.append(ecotypeid2nativename[cross_match_row.vs_ecotype_id])
				if len(true_ecotype_id_ls)>0:
					quality = 2
					no_of_mis_labelled_ones += 1
			if quality==1:
				good_ecotype_id_set.add(row.ecotype_id)
			elif quality==0 or quality==2:
				bad_ecotype_id_set.add(row.ecotype_id)
			
			if quality==2:
				true_ecotype_id_set.update(set(true_ecotype_id_ls))
			
			true_ecotype_id_ls = map(str, true_ecotype_id_ls)	#in order to use ','.join()
			output_row = [row.ecotype_id, row.nativename, row.array_original_filename, row.array_created, quality, \
						','.join(true_ecotype_id_ls), ','.join(true_nativename_ls)]
			writer.writerow(output_row)
			if counter%100==0:
				sys.stderr.write("%s\t%s\t%s\t%s"%('\x08'*80, no_of_mis_labelled_ones, no_of_bad_ones, counter))
		del writer
		
		rescued_bad_ecotype_id_set = set()
		for ecotype_id in bad_ecotype_id_set:
			if ecotype_id in good_ecotype_id_set or ecotype_id in true_ecotype_id_set:
				rescued_bad_ecotype_id_set.add(ecotype_id)
		
		# assign quality 2 to entries in true_ecotype_id_set if the ecotype_id is not in good_ecotype_id_set
		true_ecotype_id_quality_ls = []
		for ecotype_id in true_ecotype_id_set:
			if ecotype_id not in good_ecotype_id_set:
				true_ecotype_id_quality_ls.append([ecotype_id, 2])
		# remove the rescued bad ones from bad_ecotype_id_set
		bad_ecotype_id_set.difference_update(rescued_bad_ecotype_id_set)
		
		output_fname = '%s_simple%s'%(os.path.splitext(output_fname)[0], os.path.splitext(output_fname)[1])
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header_row = ['ecotype_id', 'nativename', 'quality']
		writer.writerow(header_row)
		for ecotype_id in good_ecotype_id_set:
			writer.writerow([ecotype_id, ecotypeid2nativename[ecotype_id], 1])
		for ecotype_id, quality in true_ecotype_id_quality_ls:
			writer.writerow([ecotype_id, ecotypeid2nativename[ecotype_id], quality])
		
		for ecotype_id in bad_ecotype_id_set:
			writer.writerow([ecotype_id, ecotypeid2nativename[ecotype_id], 0])
		del writer
		sys.stderr.write("%s mis-labelled out of %s bad ones, out of %s in total. Done.\n"%(no_of_mis_labelled_ones, no_of_bad_ones, counter))
		
	"""
call_method_id = 3
qc_method_id = 9
output_fname = '/tmp/250k_good_bad_cross_label_arrays.tsv'
reportBadArrays(db_250k, call_method_id, qc_method_id, output_fname)
	"""
	
	@classmethod
	def outputProbeIntensityOfOneSNPAcrossArrays(cls, db_250k, intensity_matrix_fname, snpData, snp_id, output_fname):
		"""
		2009-5-19
			intensity_matrix_fname is output of affy/CelQuantileNorm/gtype_cel_to_pq (plone doc: log/microarray)
				format is SNP X array. 1st row is name for the array and allele, with array_id embedded, like 67_raw_data_A.
					one array occupies two columns, allele A and allele B.
					1st column is snp_id (chr_pos like  1_657).
			
		"""
		sys.stderr.write("Outputting probe intensity of one SNP ... ")
		
		array_id2ecotype_id_name_ls = {}
		rows = db_250k.metadata.bind.execute("select * from view_array")
		for row in rows:
			array_id2ecotype_id_name_ls[row.array_id] = [row.maternal_ecotype_id, row.maternal_nativename]
		
		import csv
		reader = csv.reader(open(intensity_matrix_fname), delimiter='\t')
		header = reader.next()
		array_id_allele_label_ls = header[1:]
		nativename_ls = []
		ecotype_id_ls = []
		array_id_ls = []
		
		for i in range(0, len(array_id_allele_label_ls), 2):	# every other column
			label = array_id_allele_label_ls[i]
			array_id = int(label.split('_')[0])
			array_id_ls.append(array_id)
			ecotype_id_nativename_ls = array_id2ecotype_id_name_ls.get(array_id)
			if ecotype_id_nativename_ls is not None:
				ecotype_id, nativename = ecotype_id_nativename_ls
			else:
				ecotype_id = 0
				nativename = 'None'
			ecotype_id_ls.append(ecotype_id)
			nativename_ls.append(nativename)
		
		writer = csv.writer(open(output_fname,'w'), delimiter='\t')
		writer.writerow(['nativename', 'ecotype_id', 'array_id', 'alleleA-intensity', 'alleleB-intensity', 'call'])
		for row in reader:
			if row[0]==snp_id:
				intensity_ls = row[1:]
				for i in range(0, len(intensity_ls), 2):
					alleleA_intensity = intensity_ls[i]
					alleleB_intensity = intensity_ls[i+1]
					array_index = i/2
					ecotype_id = ecotype_id_ls[array_index]
					snp_row_index = snpData.row_id2row_index['%s'%ecotype_id]
					snp_col_index = snpData.col_id2col_index[snp_id]
					snp_allele = snpData.data_matrix[snp_row_index][snp_col_index]
					output_row = [nativename_ls[array_index], ecotype_id_ls[array_index], array_id_ls[array_index], alleleA_intensity, \
								alleleB_intensity, snp_allele]
					writer.writerow(output_row)
				break
		del writer, reader
		sys.stderr.write("Done.\n")
		
	"""
intensity_matrix_fname = os.path.expanduser('~/script/affy/250k_test/call_method_32_arrays.QN.tsv')
snp_id = '2_13265102'
output_fname = '/tmp/%s_intensity_matrix.tsv'%snp_id

from pymodule import SNPData
snp_fname = '/Network/Data/250k/db/dataset/call_method_32.tsv' 
snpData = SNPData(input_fname=snp_fname, turn_into_array=1, ignore_2nd_column=1)
DB250k.outputProbeIntensityOfOneSNPAcrossArrays(db_250k, intensity_matrix_fname, snpData, snp_id, output_fname)
	"""
	
	@classmethod
	def cleanUpTablePhenotype(cls, db, make_replicate_continuous=False):
		"""
		2009-8-25
			add an option (make_replicate_continuous) to make the replicate variable continuous, like 1,2,5,7 => 1,2,3,4.
			The unique key in phenotype table must be removed before applying this option. 
		2009-8-12
			increase replicate from one of the identical entries based on (method_id, ecotype_id, replicate)
		"""
		sys.stderr.write("cleaning up table phenotype ...\n")
		db.session.begin()
		
		sys.stderr.write("Getting data from db ...")
		import Stock_250kDB
		
		method_id_ecotype_id2db_id_ls = {}
		method_id_ecotype_id_replicate2count = {}
		method_id_ecotype_id2max_replicate = {}	#2009-8-25 added to find out (method_id, ecotype_id)s that have the discontinuous replicate numbers
		
		rows = Stock_250kDB.Phenotype.query.all()
		for row in rows:
			method_id_ecotype_id = (row.method_id, row.ecotype_id)
			method_id_ecotype_id_replicate = (row.method_id, row.ecotype_id, row.replicate)
			if method_id_ecotype_id not in method_id_ecotype_id2db_id_ls:
				method_id_ecotype_id2db_id_ls[method_id_ecotype_id] = []
			method_id_ecotype_id2db_id_ls[method_id_ecotype_id].append(row.id)
			
			if method_id_ecotype_id_replicate not in method_id_ecotype_id_replicate2count:
				method_id_ecotype_id_replicate2count[method_id_ecotype_id_replicate] = 0
			method_id_ecotype_id_replicate2count[method_id_ecotype_id_replicate] += 1
			
			#2009-8-25
			if method_id_ecotype_id not in method_id_ecotype_id2max_replicate:
				method_id_ecotype_id2max_replicate[method_id_ecotype_id] = 0
			if row.replicate>method_id_ecotype_id2max_replicate[method_id_ecotype_id]:
				method_id_ecotype_id2max_replicate[method_id_ecotype_id] = row.replicate
		sys.stderr.write("Done.\n")
		
		method_id_ecotype_id_need_fix_set = set()
		for method_id_ecotype_id_replicate, count in method_id_ecotype_id_replicate2count.iteritems():
			if count>1:
				method_id_ecotype_id = (method_id_ecotype_id_replicate[0], method_id_ecotype_id_replicate[1])
				method_id_ecotype_id_need_fix_set.add(method_id_ecotype_id)
		
		#2009-8-25
		if make_replicate_continuous:
			for method_id_ecotype_id, max_replicate in method_id_ecotype_id2max_replicate.iteritems():
				if max_replicate>len(method_id_ecotype_id2db_id_ls[method_id_ecotype_id]):
					method_id_ecotype_id_need_fix_set.add(method_id_ecotype_id)
		
		sys.stderr.write("%s (method_id, ecotype_id) pairs to be fixed.\n"%len(method_id_ecotype_id_need_fix_set))
		
		#re-assign the replicate value in the order of db ids
		for method_id_ecotype_id in method_id_ecotype_id_need_fix_set:
			db_id_ls = method_id_ecotype_id2db_id_ls.get(method_id_ecotype_id)
			if db_id_ls is None:
				continue
			
			db_id_ls.sort()
			for i in range(len(db_id_ls)):
				db_id = db_id_ls[i]
				db_entry = Stock_250kDB.Phenotype.get(db_id)
				db_entry.replicate = i +1
				db.session.add(db_entry)
				db.session.flush()
		db.session.commit()
		sys.stderr.write("Done.\n")
	
	"""	
DB250k.cleanUpTablePhenotype(db_250k)

DB250k.cleanUpTablePhenotype(db_250k, make_replicate_continuous=True)
	"""
	
	@classmethod
	def updatePhenotypeAvgEntry(cls, db, db_entry, individual_value_ls):
		"""
		2009-8-13
			
		"""
		import numpy
		db_entry.value = numpy.average(individual_value_ls)
		db_entry.sample_size = len(individual_value_ls)
		if db_entry.sample_size>1:
			db_entry.stddev = numpy.std(individual_value_ls) 
		db.session.add(db_entry)
		
	
	@classmethod
	def cleanUpTablePhenotypeAvg(cls, db):
		"""
		2009-8-12
			enforce the unique constraint, (method_id, ecotype_id) 
			
		"""
		sys.stderr.write("Cleaning up table phenotype_avg ... \n")
		db.session.begin()
		
		sys.stderr.write("Getting data from db ...")
		import Stock_250kDB
		
		method_id_ecotype_id2db_id_ls = {}
		
		rows = Stock_250kDB.PhenotypeAvg.query.all()
		for row in rows:
			method_id_ecotype_id = (row.method_id, row.ecotype_id)
			if method_id_ecotype_id not in method_id_ecotype_id2db_id_ls:
				method_id_ecotype_id2db_id_ls[method_id_ecotype_id] = []
			method_id_ecotype_id2db_id_ls[method_id_ecotype_id].append(row.id)
		
		sys.stderr.write("Done.\n")
		import numpy
		no_of_replicate_entries = 0
		no_of_replicate_entries_but_no_individual_data = 0
		no_of_replicate_entries_with_individual_data = 0
		for method_id_ecotype_id, db_id_ls in method_id_ecotype_id2db_id_ls.iteritems():
			if len(db_id_ls)>1:
				no_of_replicate_entries  += 1
				method_id, ecotype_id = method_id_ecotype_id
				rows = Stock_250kDB.Phenotype.query.filter_by(method_id=method_id).filter_by(ecotype_id=ecotype_id)
				db_id_ls.sort()
				if rows.count()>0:
					no_of_replicate_entries_with_individual_data += 1
					for i in range(len(db_id_ls)-1):	#delete all but the last one
						db_id = db_id_ls[i] 
						db_entry = Stock_250kDB.PhenotypeAvg.get(db_id)
						
						db.session.delete(db_entry)
					#update the last PhenotypeAvg entry
					db_entry = Stock_250kDB.PhenotypeAvg.get(db_id_ls[-1])
					individual_value_ls = [row.value for row in rows]
					cls.updatePhenotypeAvgEntry(db, db_entry, individual_value_ls)
				else:
					no_of_replicate_entries_but_no_individual_data += 1
					individual_value_ls = []
					for i in range(len(db_id_ls)-1):	#insert into Stock_250kDB.Phenotype and delete all but the last one from Stock_250kDB.PhenotypeAvg
						db_id = db_id_ls[i] 
						db_entry = Stock_250kDB.PhenotypeAvg.get(db_id)
						individual_value_ls.append(db_entry.value)
						phenotype_entry = Stock_250kDB.Phenotype(method_id=db_entry.method_id, ecotype_id=db_entry.ecotype_id, \
																value=db_entry.value, replicate=i+1)
						db.session.add(phenotype_entry)
						db.session.delete(db_entry)
					#insert the last phenotype_avg into Stock_250kDB.Phenotype and update it in  phenotype_avg
					db_entry = Stock_250kDB.PhenotypeAvg.get(db_id_ls[-1])
					individual_value_ls.append(db_entry.value)
					phenotype_entry = Stock_250kDB.Phenotype(method_id=db_entry.method_id, ecotype_id=db_entry.ecotype_id, \
																value=db_entry.value, replicate=len(db_id_ls))
					db.session.add(phenotype_entry)
					cls.updatePhenotypeAvgEntry(db, db_entry, individual_value_ls)
					#insert data from db_id_ls into Stock_250kDB.Phenotype
		db.session.commit()
		sys.stderr.write("%s total replicate entries. %s has individual data. %s has no individual data. Done.\n"%\
						(no_of_replicate_entries, no_of_replicate_entries_with_individual_data, no_of_replicate_entries_but_no_individual_data))
		
	"""
DB250k.cleanUpTablePhenotypeAvg(db)
	"""
	
	@classmethod
	def updatePhenotypeAvgBasedOnPhenotype(cls, db, phenotype_condition=None):
		"""
		2009-8-25
			After values in table phenotype are updated, corresponding ones in table phenotype_avg shall be updated too.
		"""
		sys.stderr.write("Updating phenotype_avg entries based on the ones in table phenotype ...\n")
		db.session.begin()
		
		method_id_ecotype_id2individual_value_ls = {}
		rows = db.metadata.bind.execute("select method_id, ecotype_id, value from phenotype %s"%phenotype_condition)
		replicate_count = 0
		for row in rows:
			method_id_ecotype_id = (row.method_id, row.ecotype_id)
			if method_id_ecotype_id not in method_id_ecotype_id2individual_value_ls:
				method_id_ecotype_id2individual_value_ls[method_id_ecotype_id] = []
			method_id_ecotype_id2individual_value_ls[method_id_ecotype_id].append(row.value)
			replicate_count += 1
		
		sys.stderr.write("need to update %s phenotype_avg entries.\n"%(len(method_id_ecotype_id2individual_value_ls)))
		
		avg_count = 0
		new_avg_count = 0
		multi_avg_count = 0
		import Stock_250kDB
		for method_id_ecotype_id, individual_value_ls in method_id_ecotype_id2individual_value_ls.iteritems():
			
			method_id, ecotype_id = method_id_ecotype_id
			#rows = db.metadata.bind.execute("select ecotype_id, avg(value) as avg_value, stddev(value) as stddev, count(value) as sample_size, \
			#	method_id from phenotype where \
			#	method_id=%s and ecotype_id=%s group by method_id, ecotype_id"%(method_id, ecotype_id))
			rows = Stock_250kDB.PhenotypeAvg.query.filter_by(method_id=method_id).filter_by(ecotype_id=ecotype_id)
			if rows.count()==1:
				phenotype_avg_entry = rows.one()
			elif rows.count()==0:
				phenotype_avg_entry = Stock_250kDB.PhenotypeAvg(method_id=method_id, ecotype_id=ecotype_id)
				new_avg_count += 1
			else:
				multi_avg_count += 1
				sys.stderr.write("Error: method_id=%s, ecotype_id=%s has multiple phenotype_avg entries.\n"%(method_id, ecotype_id))
				phenotype_avg_entry = rows.first()
				sys.exit(3)
			DB250k.updatePhenotypeAvgEntry(db, phenotype_avg_entry, individual_value_ls)
			#phenotype_avg_entry.value = row.avg_value
			#phenotype_avg_entry.stddev = row.stddev
			#phenotype_avg_entry.sample_size = row.sample_size
			#db.session.add(phenotype_avg_entry)
			avg_count += 1
		db.session.commit()
		sys.stderr.write("%s(%s new, %s multi) average values from %s replicates were updated. Done.\n"%\
						(avg_count, new_avg_count, multi_avg_count, replicate_count))
	
	"""
DB250k.updatePhenotypeAvgBasedOnPhenotype(db_250k, phenotype_condition='where date_created>"2009-08-23" order by method_id, ecotype_id')

#update all phenotype_avg entries that have inidividual values in table phenotype 
DB250k.updatePhenotypeAvgBasedOnPhenotype(db_250k);
	"""	
	
	@classmethod
	def convertOldFormatCallFileIntoNewFormat(cls, db_250k, method_id=3, priorTAIRVersion=False):
		"""
		2011-2-24
			add argument priorTAIRVersion, if true, it means using Snps.tair8_chromosome,Snps.tair8_position
				rather than Snps.chromosome,Snps.position.
		2011-1-24
			use Stock_250kDB.Snps.id to replace chr_pos... in the call file.
			i.e.
				old format: 1_3102_A_G      G       0.985079453549666
				new format: 2       G       0.985079453549666
		"""
		chr_pos2db_id = db_250k.getSNPChrPos2ID(priorTAIRVersion=priorTAIRVersion)
		import Stock_250kDB, os, sys, csv
		sys.stderr.write("Converting old format call files from method %s into new format ... \n"%(method_id))
		query = Stock_250kDB.CallInfo.query.filter_by(method_id=method_id)
		counter = 0
		no_of_files = query.count()
		for row in query:
			sys.stderr.write("%d/%d:\t%s "%(counter+1, no_of_files, row.filename))
			reader = csv.reader(open(row.filename), delimiter='\t')
			header_row = reader.next()
			first_data_row = reader.next()
			SNP_id_ls = first_data_row[0].split('_')
			if len(SNP_id_ls)>1:	#
				convertData = True
			else:
				convertData = False
			del reader
			
			if convertData:
				no_of_lines = 0
				no_of_lines_after_conversion = 0
				sys.stderr.write("...")
				# backup the old file first
				oldFormatFname = '%s.old.tsv'%(os.path.splitext(row.filename)[0])
				os.rename(row.filename, oldFormatFname)
				
				# writing data into the old filename in new format
				reader = csv.reader(open(oldFormatFname), delimiter='\t')
				header_row = reader.next()
				writer = csv.writer(open(row.filename, 'w'), delimiter='\t')
				writer.writerow(header_row)
				for row in reader:
					no_of_lines += 1
					chr_pos = row[0].split('_')[:2]
					chr_pos = tuple(map(int, chr_pos))
					db_id = chr_pos2db_id.get(chr_pos)
					if db_id is not None:
						no_of_lines_after_conversion += 1
						new_row = [db_id] + row[1:]
						writer.writerow(new_row)
				del reader, writer
				sys.stderr.write("%s lines different after conversion.\n"%(no_of_lines_after_conversion-no_of_lines))
			else:
				sys.stderr.write("no conversion. (maybe converted already).\n")
			counter += 1
	
	"""
		# 2011-1-24
		# 
		DB250k.convertOldFormatCallFileIntoNewFormat(db_250k, method_id=3)
		sys.exit(0)
		
		# 2011-2-15
		# 	Watch: use papaya's TAIR8 db because the call info files contains calls in TAIR8 coordinates.
		for call_method_id in xrange(20,73):
			DB250k.convertOldFormatCallFileIntoNewFormat(db_250k, method_id=call_method_id)
		sys.exit(0)
		
	"""
	
	@classmethod
	def convertSNPDatasetLocusIDIntoDBID(cls, db, input_fname, output_fname):
		"""
		2011-2-24
			input_fname is Strain X SNP format.
			This function replaces SNP id in chr_pos format with db id. Snps.id.
			
			
		"""
		from pymodule import read_data
		chr_pos_set = set()
		chr_pos2snp_id = db.getSNPChrPos2ID(keyType=1)
		#chr_pos2snp_id = db.chr_pos2snp_id
		for chr_pos, db_id in chr_pos2snp_id.iteritems():
			chr, pos = chr_pos[:2]
			chr_pos = '%s_%s'%(chr, pos)
			chr_pos_set.add(chr_pos)
		
		header, strain_acc_list, category_list, data_matrix = read_data(input_fname, col_id_key_set=chr_pos_set)
		
		new_header = header[:2]
		for chr_pos in header[2:]:
			chr_pos = map(int, chr_pos.split('_'))
			chr_pos = tuple(chr_pos)
			db_id = chr_pos2snp_id[chr_pos]
			new_header.append(db_id)
		header = new_header
		from pymodule.SNP import write_data_matrix
		write_data_matrix(data_matrix, output_fname, header, strain_acc_list, category_list, rows_to_be_tossed_out=None, \
			cols_to_be_tossed_out=None, nt_alphabet=0, transform_to_numpy=0,\
			discard_all_NA_rows=0, strain_acc2other_info=None, delimiter='\t', )
	"""
		# 2011-2-24
		input_fname = "/Network/Data/250k/db/reference_dataset/2010_149_384_20091005.csv"
		output_fname = "/Network/Data/250k/db/reference_dataset/2010_149_384_20091005_db_id.tsv"
		DB250k.convertSNPDatasetLocusIDIntoDBID(db, input_fname, output_fname)
		sys.exit(2)
		
		# 2011-2-24 be careful which db (old TAIR8 or TAIR9) to use. this case is TAIR8 (papaya).
		input_fname = "/Network/Data/250k/db/reference_dataset/2010_149_384_20091005.csv"
		output_fname = "/Network/Data/250k/db/reference_dataset/2010_149_384_20091005_db_id.tsv"
		DB250k.convertSNPDatasetLocusIDIntoDBID(db_250k, input_fname, output_fname)
		sys.exit(2)
	"""
	
	@classmethod
	def convertOldFormatResultMethodFileIntoNewFormat(cls, db_250k, call_method_id=None, priorTAIRVersion=False):
		"""
		2011-5-4
			use Stock_250kDB.Snps.id to replace chr, pos ... in the files associated with table ResultsMethod.
			
			i.e.
				old format: 1	657	0.985079453549666	...
				new format: 1		0.985079453549666
			
			argument priorTAIRVersion, if true, it means using Snps.tair8_chromosome,Snps.tair8_position
				rather than Snps.chromosome,Snps.position.
		"""
		import Stock_250kDB, os, sys, csv
		from pymodule import SNP
		sys.stderr.write("Converting old format results_method files from method %s into new format ... \n"%(call_method_id))
		chr_pos2db_id = db_250k.getSNPChrPos2ID(priorTAIRVersion=priorTAIRVersion)
		
		#type 2 is segment-based recombination rate. only one result.
		# analysis method 13 is BooleanSNPPair. ignore
		TableClass = Stock_250kDB.ResultsMethod
		query = TableClass.query.filter(TableClass.results_method_type_id!=2).filter(TableClass.analysis_method_id!=13)
		if call_method_id:
			query = query.filter_by(call_method_id=call_method_id)
		counter = 0
		no_of_files = query.count()
		for row in query:
			sys.stderr.write("%d/%d:\t%s "%(counter+1, no_of_files, row.filename))
			reader = csv.reader(open(row.filename), delimiter='\t')
			first_data_row = reader.next()
			hasHeader = False
			if SNP.pa_has_characters.search(first_data_row[0]):	#first column of 1st row has character in it. this file contains header. 
				first_data_row = reader.next()
				hasHeader = True
			first_data_row = reader.next()
			if first_data_row[1] and first_data_row[1]!='0':	#2011-5-3 2nd column is something other than '0'. conversion is needed.
				convertData = True
			else:
				convertData = False
			del reader
			
			if convertData:
				no_of_lines = 0
				no_of_lines_after_conversion = 0
				sys.stderr.write("...")
				# backup the old file first
				oldFormatFname = '%s.old.tsv'%(os.path.splitext(row.filename)[0])
				os.rename(row.filename, oldFormatFname)
				
				# writing data into the old filename in new format
				reader = csv.reader(open(oldFormatFname), delimiter='\t')
				writer = csv.writer(open(row.filename, 'w'), delimiter='\t')
				if hasHeader:
					header_row = reader.next()
					writer.writerow(['snp_id', 'none']+ header_row[2:])
				for row in reader:
					no_of_lines += 1
					chr_pos = row[:2]
					chr_pos = (chr_pos[0], int(chr_pos[1]))	#chr is of type str. pos is of type int. in chr_pos2db_id
					db_id = chr_pos2db_id.get(chr_pos)
					if db_id is not None:
						no_of_lines_after_conversion += 1
						new_row = [db_id, ''] + row[2:]	#2nd column (pos) should be empty in new format.
						writer.writerow(new_row)
				del reader, writer
				sys.stderr.write("%s lines different after conversion.\n"%(no_of_lines_after_conversion-no_of_lines))
			else:
				sys.stderr.write("no conversion. (maybe converted already).\n")
			counter += 1
	
	"""
		# 2011-5-4
		# priorTAIRVersion=True if it's connected to banyan's TAIR9 db while the association results are based off TAIR8.
		DB250k.convertOldFormatResultMethodFileIntoNewFormat(db_250k, call_method_id=1, priorTAIRVersion=True)
		sys.exit(0)
		
		# 2011-2-15
		# 	Watch: use papaya's TAIR8 db because the call info files contains calls in TAIR8 coordinates.
		for call_method_id in xrange(20,73):
			DB250k.convertOldFormatCallFileIntoNewFormat(db_250k, call_method_id=call_method_id, priorTAIRVersion=True)
		sys.exit(0)
		
	"""
	
	@classmethod
	def convertRMIntoJson(cls, db_250k, call_method_id_ls_str="", analysis_method_id_ls_str="", \
						phenotype_method_id_ls_str="", no_of_top_snps = 10000, commit=False):
		"""
		2010-5-3
			
		"""
		if commit:
			pass
		else:
			# begin the session but without commit
			db_250k.session.begin()
		from pymodule import getListOutOfStr
		import Stock_250kDB
		from Stock_250kDB import ResultsMethodJson, ResultsMethod
		query = ResultsMethod.query
		call_method_id_ls = getListOutOfStr(call_method_id_ls_str, data_type=int)
		if call_method_id_ls:
			query = query.filter(ResultsMethod.call_method_id.in_(call_method_id_ls))
		analysis_method_id_ls = getListOutOfStr(analysis_method_id_ls_str, data_type=int)
		if analysis_method_id_ls:
			query = query.filter(ResultsMethod.analysis_method_id.in_(analysis_method_id_ls))
		phenotype_method_id_ls = getListOutOfStr(phenotype_method_id_ls_str, data_type=int)
		if phenotype_method_id_ls:
			query = query.filter(ResultsMethod.phenotype_method_id.in_(phenotype_method_id_ls))
		
		no_of_total = query.count()
		count = 0
		for rm in query:
			count += 1
			sys.stderr.write("\t %s/%s, %s-%s-%s "%(count, no_of_total, rm.call_method_id, \
												rm.analysis_method_id, rm.phenotype_method_id))
			if rm.analysis_method.min_maf is not None:
				min_MAF = rm.analysis_method.min_maf
			else:
				min_MAF = 0
			rm_json = ResultsMethodJson.query.filter_by(results_id=rm.id).\
				filter(ResultsMethodJson.min_MAF<=min_MAF+0.0001).\
				filter(ResultsMethodJson.min_MAF>=min_MAF-0.0001).\
				filter_by(no_of_top_snps=no_of_top_snps).first()
			if rm_json:
				sys.stderr.write("json for result %s (%s, %s, %s) already exists in db.\n"%\
								(rm.id, rm.call_method_id, rm.analysis_method_id, rm.phenotype_method_id))
				continue
			from common import getOneResultJsonData
			json_data = getOneResultJsonData(rm, min_MAF, no_of_top_snps)
			rm_json = ResultsMethodJson(min_MAF=min_MAF, no_of_top_snps=no_of_top_snps)
			rm_json.result = rm
			rm_json.json_data = json_data
			db_250k.session.add(rm_json)
			db_250k.session.flush()
			sys.stderr.write("\n")
		if commit:
			db_250k.session.flush()
			db_250k.session.commit()
		else:
			pass
			#db_250k.session.commit()
	
	"""
	DB250k.convertRMIntoJson(db_250k, call_method_id_ls_str='54', analysis_method_id_ls_str='1,4,7,47,49', \
		phenotype_method_id_ls_str='409-603', commit=True)
	
	DB250k.convertRMIntoJson(db_250k, call_method_id_ls_str='29,32,43-56', commit=True)
	
	"""
	
	@classmethod
	def updateProbeChrPosToTAIR9(cls, db_250k, input_fname, min_no_of_matches=25):
		"""
		2010-6-20
			input_fname is output of MpiBlast.py
				mpiexec ~/script/variation/Dazhe/src/MpiBlast.py -o ~/script/variation/data/CNV/stock_250k_probes_vs_TAIR9.tsv 
				-p ~/script/variation/data/CNV/stock_250k_probes.csv -d ~/script/variation/data/CNV/TAIR9/tair9.fas -m 24
			
			If one probe has two perfect-matches in TAIR9, its chromosome, position would be set to NULL.
		"""
		from pymodule import getColName2IndexFromHeader, PassingData, figureOutDelimiter
		sys.stderr.write("Reading from %s ... \n"%input_fname)
		counter = 0
		real_counter = 0
		import csv, re
		reader = csv.reader(open(input_fname), delimiter=figureOutDelimiter(input_fname))
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header)
		chr_pattern = re.compile(r'.*Chr(\d)')
		probe_id2new_chr_pos_ls = {}
		
		for row in reader:
			chr = int(row[col_name2index['Chromosome']])
			pos = int(row[col_name2index['Position']])
			probe_id = int(row[col_name2index['Probe_ID']])
			tair9_label = row[col_name2index['Alignment_title']]	#Alignment_title is like "gnl|BL_ORD_ID|0 Chr1"
			
			new_chr = int(chr_pattern.search(tair9_label).group(1))
			
			no_of_matches = int(row[col_name2index['Number_matches']])
			alignment_start = int(row[col_name2index['Alignment_start']])
			alignment_stop = int(row[col_name2index['Alignment_stop']])
			new_pos = (alignment_start+alignment_stop)/2
			query_start = int(row[col_name2index['query_start']])
			if no_of_matches>=min_no_of_matches:
				if probe_id not in probe_id2new_chr_pos_ls:
					probe_id2new_chr_pos_ls[probe_id] = []
				probe_id2new_chr_pos_ls[probe_id].append((new_chr, new_pos))
				real_counter += 1
			counter += 1
			if counter%5000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
		sys.stderr.write("%s%s\t%s\n"%('\x08'*80, counter, real_counter))
		sys.stderr.write("%s probes, with %s matches from %s lines.\n"%(len(probe_id2new_chr_pos_ls), real_counter, counter))
		del reader
		
		sys.stderr.write("Updating database ... \n")
		import Stock_250kDB
		i = 0
		block_size = 1000
		no_of_new_chr_pos_updates = 0
		TableClass = Stock_250kDB.Probes
		rows = TableClass.query.offset(i).limit(block_size)
		while rows.count()!=0:
			for row in rows:
				new_chr_pos_ls = probe_id2new_chr_pos_ls.get(row.id)
				if new_chr_pos_ls is None:
					new_chr_pos = (None, None)
				elif len(new_chr_pos_ls)==1:
					new_chr_pos = new_chr_pos_ls[0]
				else:	#>1 matches in TAIR9
					new_chr_pos = (None, None)
				
				new_chr = new_chr_pos[0]
				new_pos = new_chr_pos[1]
				if new_chr!=row.chromosome or new_pos!=row.position:
					row.chromosome = new_chr
					row.position = new_pos
					db_250k.session.add(row)
					db_250k.session.flush()
					no_of_new_chr_pos_updates += 1
				
				i += 1
			if i%5000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, i, no_of_new_chr_pos_updates))
			rows = TableClass.query.offset(i).limit(block_size)
		sys.stderr.write("%s%s\t%s\n"%('\x08'*80, i, no_of_new_chr_pos_updates))
		sys.stderr.write("Out of %s, %s entries updated to new chr & pos, the rest set to Null. Done.\n"%(i, no_of_new_chr_pos_updates))
	"""
	input_fname = os.path.expanduser('~/mnt/hpc-cmb/script/variation/data/CNV/stock_250k_probes_vs_TAIR9.tsv')
	DB250k.updateProbeChrPosToTAIR9(db_250k, input_fname, min_no_of_matches=25)
	
	#2010-6-26 the hpc-cmb result is incomplete due to early MPI exit.
	input_fname = os.path.expanduser('~/script/variation/data/CNV/stock_250k_probes_vs_TAIR9_20100625.tsv')
	DB250k.updateProbeChrPosToTAIR9(db_250k, input_fname, min_no_of_matches=25)
	sys.exit(0)
	
	# 2010-6-26 After function above has finished, run this sql manually in a mysql client. 
		update snps s, probes p set s.chromosome=p.chromosome, s.position=p.position where s.id=p.snps_id;
	"""
		
	@classmethod
	def updateSNPIncludeAfterQC(cls, db_250k, call_method_id=72):
		"""
		2010-10-19
			update snps.include_after_qc based on a final snp dataset (call_method_id).
			
			call info file from this call method has Snps.id as first column, rather than chr_pos (old format).
		"""
		from pymodule import getColName2IndexFromHeader, PassingData, figureOutDelimiter
		sys.stderr.write("updating snps.include_after_qc based on one random call info from call method %s ... \n"%call_method_id)
		import Stock_250kDB
		call_info = Stock_250kDB.CallInfo.query.filter_by(method_id=call_method_id).first()
		counter = 0
		real_counter = 0
		import csv
		reader = csv.reader(open(call_info.filename), delimiter=figureOutDelimiter(call_info.filename))
		header = reader.next()
		for row in reader:
			snp_id = int(row[0])
			snp = Stock_250kDB.Snps.get(snp_id)
			counter += 1
			if snp:
				snp.include_after_qc = 1
				db_250k.session.add(snp)
				real_counter += 1
			if counter%5000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
		db_250k.session.flush()
		sys.stderr.write("%s%s\t%s\n"%('\x08'*80, counter, real_counter))
		sys.stderr.write("%s SNPs out of %s in the file have include_after_qc set. Done.\n"%(real_counter , counter))
		db_250k.session.commit()
	
	"""
	#2010-10-19
	DB250k.updateSNPIncludeAfterQC(db_250k)
	sys.exit(0)
	
	"""
	
	
class CNV(object):
	class Lyrata(object):
		"""
		2010-8-2
		"""
		def __init__(self):
			pass
		
		@classmethod
		def getLyrataNormalANDDeletionFromQuanFileInTAIR9(cls, input_fname, output_fname=None):
			"""
			2010-8-2
			"""
			sys.stderr.write("Getting lyrata sequences or deletions in TAIR9 coordinates ... \n")
			
			inf = open(input_fname)
			chromosome = None
			position = 0
			normalSegmentStart = None
			previousLetter = None
			normalSegmentLs = []
			lineNo = 0
			for line in inf:
				lineNo += 1
				line = line.strip()
				if line[0]=='>':
					if normalSegmentStart is not None and position>normalSegmentStart:
						#there's a normal segment in the end of the previous chromosome
						normal_segment = [chromosome, normalSegmentStart, position]
						normalSegmentLs.append(normal_segment)
					chromosome = line[4]
					position = 0
					normalSegmentStart = None
					previousLetter = None
				else:
					no_of_letters = len(line)
					for i in xrange(no_of_letters):
						position += 1
						letter = line[i]
						if previousLetter is None:
							previousLetter = letter
						elif previousLetter!='-' and letter=='-':	#from normal to deletion
							if normalSegmentStart is None:	#beginning of the chromosome
								normalSegmentStart = 1
							normal_segment = [chromosome, normalSegmentStart, position-1]
							normalSegmentLs.append(normal_segment)
							normalSegmentStart = None
						elif previousLetter=='-' and letter!='-':	#from deletion to normal
							normalSegmentStart = position
						previousLetter = letter
				if lineNo%10000==0:
					sys.stderr.write('%s%s\t%s\t%s'%('\x08'*80, 'Chr%s'%chromosome, lineNo, len(normalSegmentLs)))
			sys.stderr.write('%s%s\t%s\t%s\n'%('\x08'*80, 'Chr%s'%chromosome, lineNo, len(normalSegmentLs)))
			if normalSegmentStart is not None and position>normalSegmentStart:
				#there's a normal segment in the end of the previous chromosome
				normal_segment = [chromosome, normalSegmentStart, position]
				normalSegmentLs.append(normal_segment)
			
			sys.stderr.write("Outputting result ...")
			import csv
			writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
			
			header_row = ['start_probe_id', 'start_chr_pos', 'stop_probe_id', 'stop_chr_pos', 'no_of_probes', \
						'length', 'copy_number', 'size_difference']
			writer.writerow(header_row)
			for segment in normalSegmentLs:
				chromosome, normalSegmentStart, stop = segment[:3]
				data_row = ['NA', '%s_%s'%(chromosome, normalSegmentStart), 'NA', '%s_%s'%(chromosome, stop),\
						'NA', stop-normalSegmentStart+1, 1, 0]
				writer.writerow(data_row)
			del writer
			sys.stderr.write("Done.\n")
		
		"""
		#2010-8-2
		input_fname = os.path.expanduser('~/script/variation/data/lyrata/at_ancestor_t9.fa')
		output_fname = os.path.expanduser('~/script/variation/data/lyrata/at_ancestor_t9.normaSegment.tsv')
		CNV.Lyrata.getLyrataNormalANDDeletionFromQuanFileInTAIR9(input_fname, output_fname=output_fname)
		sys.exit(0)
		
		#2010-8-2
		input_fname = os.path.expanduser('~/script/variation/data/lyrata/at_ancestor_t9.fa')
		output_fname = os.path.expanduser('~/script/variation/data/lyrata/at_ancestor_t9.normaSegment.tsv')
		CNV.Lyrata.getLyrataNormalANDDeletionFromQuanFileInTAIR9(input_fname, output_fname=output_fname)
		sys.exit(0)
		
		"""
	class HMMSegmentation(object):
		"""
		2010-6-3
			test whether HMM would be feasible.
		"""
		@classmethod
		def outputViterbiPath(cls, writer, array_id, ecotype_id, probe_id_ls, chr_pos_ls, path, base_start_index=0):
			"""
			2010-6-3
			"""
			path = path[0]
			no_of_probes = len(path)
			old_state = None
			state_start_index = None
			state_stop_index = None
			for i in range(no_of_probes):
				if state_start_index is None:
					state_start_index = base_start_index + i
				if old_state is not None and path[i]!=old_state:
					state_stop_index = base_start_index + i - 1
					start_chr_pos = chr_pos_ls[state_start_index]
					stop_chr_pos = chr_pos_ls[state_stop_index]	#prevous one
					segment_len = int(stop_chr_pos[1]) - int(start_chr_pos[1]) + 25
					start_probe_id = probe_id_ls[state_start_index]
					stop_probe_id = probe_id_ls[state_stop_index]
					no_of_probes = state_stop_index - state_start_index + 1
					new_row = [ecotype_id, array_id, '_'.join(start_chr_pos), '_'.join(stop_chr_pos), no_of_probes, old_state, start_probe_id, \
							stop_probe_id, segment_len]
					writer.writerow(new_row)
					
					# reset the state start index
					state_start_index = base_start_index + i	# no -1
				old_state = path[i]
			if state_stop_index is None or state_stop_index<state_start_index:	# last segment isn't outputted
				state_stop_index = base_start_index + no_of_probes - 1
				
				start_chr_pos = chr_pos_ls[state_start_index]
				stop_chr_pos = chr_pos_ls[state_stop_index]	#prevous one
				segment_len = int(stop_chr_pos[1]) - int(start_chr_pos[1]) + 25
				start_probe_id = probe_id_ls[state_start_index]
				stop_probe_id = probe_id_ls[state_stop_index]
				no_of_probes = state_stop_index - state_start_index + 1
				new_row = [ecotype_id, array_id, '_'.join(start_chr_pos), '_'.join(stop_chr_pos), no_of_probes, old_state, start_probe_id, \
						stop_probe_id, segment_len]
				writer.writerow(new_row)
		
		@classmethod
		def viterbiPath(cls, input_fname, output_fname):
			"""
			2010-6-3
			"""
			import ghmm
			# example code for a continuous HMM with gaussian emissions
			
			F = ghmm.Float()  # emission domain of this model
			
			A = [[0.995, 0.0025, 0.0025], [0.0025, 0.995, 0.0025], [0.0025, 0.0025, 0.995]]   # transition matrix
			B = [[-0.4, 0.0025], [0.0, 0.0002], [0.4,0.002]]   # parameters of emission distributions in pairs of (mu, sigma)
			
			# Parameters of emission distributions 
			# Interpretation of B matrix for the mixture case (Example with three states and three components each):
			#  B = [ 
			#      [ ["mu11","mu12","mu13"],["sig11","sig12","sig13"],["w11","w12","w13"]   ],
			#      [  ["mu21","mu22","mu23"],["sig21","sig22","sig23"],["w21","w22","w23"]  ],
			#      [  ["mu31","mu32","mu33"],["sig31","sig32","sig33"],["w31","w32","w33"]  ],
			#      ]
			
			B = [[ [-0.4,0.0],[0.0025, 0.0002], [0.8,0.2]],
			     [ [-0.4,0.0,0.3],[0.0025, 0.0002, 0.002], [0.05,0.9,0.05]],
			     [ [0.0,0.3,0.6],[0.0002,0.002,0.002], [0.1,0.5,0.4]]]
			pi = [0.2, 0.7, 0.1]   # initial probabilities per state
			model = ghmm.HMMFromMatrices(F, ghmm.GaussianMixtureDistribution(F), A, B, pi)
			# re-normalize model parameters
			model.normalize()
			print model
			
			from pymodule import figureOutDelimiter
			import csv
			
			from pymodule.CNV import ArrayXProbeFileWrapper, get_chr2start_stop_index
			inputWrapper = ArrayXProbeFileWrapper(input_fname)
			
			reader = inputWrapper.reader
			probe_id_ls = inputWrapper.probe_id_ls
			chr_pos_ls = inputWrapper.chr_pos_ls
			
			
			chr2start_stop_index = get_chr2start_stop_index(chr_pos_ls)
			
			output_dir = os.path.split(output_fname)[0]
			if not os.path.isdir(output_dir):
				os.makedirs(output_dir)
			
			writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
			from RunGADA import RunGADA
			RunGADA.output_header(writer)
			
			import numpy
			counter = 0
			for row in reader:
				array_id = int(row[0])
				ecotype_id = row[1]
				intensity_ls = row[2:]
				intensity_array = map(float, intensity_ls)
				for chr, start_stop_index in chr2start_stop_index.iteritems():
					start_index, stop_index = start_stop_index[:2]
					observations = ghmm.EmissionSequence(ghmm.Float(), intensity_array[start_index:stop_index+1])
					
					#model.baumWelch(observations)
					#print "model after training by data from array %s e-id %s, chr %s"%(array_id, ecotype_id, chr)
					#print model
					path = model.viterbi(observations)
					cls.outputViterbiPath(writer, array_id, ecotype_id, probe_id_ls, chr_pos_ls, path, start_index)
	
	"""
	input_fname = os.path.expanduser('~/script/variation/data/CNV/call_48_Col-Ler-WABQN_b200_j100_lts.tsv')
	from pymodule.utils import addExtraToFilenamePrefix, addExtraLsToFilenamePrefix
	output_fname = addExtraToFilenamePrefix(input_fname, '3StateHMM')
	CNV.HMMSegmentation.viterbiPath(input_fname, output_fname)
	sys.exit(0)
	"""
	
	@classmethod
	def drawCNVAmpHist(cls, snpData, output_dir, max_counter=None, start_pos=None, stop_pos=None, chromosome_chosen=None):
		"""
		2009-2-20
			add arguments start_pos, stop_pos
		2008-12-12 draw histogram of CNV amplitudes probe by probe. input file is the amplitude output of RunGADA.py.
		"""
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		import pylab
		counter = 0
		for col_id in snpData.col_id_ls:
			col_id_split = col_id.split('_')
			col_id_split = map(int, col_id_split)
			chromosome, position = col_id_split
			
			#filter according to chromosome
			if chromosome_chosen is not None and chromosome!=chromosome_chosen:
				continue
			#filter according to position
			if start_pos is not None and stop_pos is not None and (position<start_pos-12 or position>stop_pos+12):
				continue
			
			col_index = snpData.col_id2col_index[col_id]
			sys.stderr.write("%s\t%s"%('\x08'*20, counter))
			output_fname = os.path.join(output_dir, '%s_amp_hist.png'%col_id)
			amp_ls = snpData.data_matrix[:,col_index]
			pylab.clf()
			pylab.hist(amp_ls, 40, alpha=0.6)
			pylab.title(col_id)
			pylab.xlabel('amplitude')
			pylab.ylabel('frequency')
			pylab.xlim([-1,1])
			pylab.savefig(output_fname, dpi=200)
			counter += 1
			if max_counter and counter>max_counter:
				break
	"""
		input_fname = '/Network/Data/250k/tmp-yh/CNV/call_method_17_CNV_array_intensity_chr4_line_no_888148_1107622_norm_GADA_out_amp.tsv'
		from pymodule import SNPData
		snpData = SNPData(input_fname=input_fname, turn_into_integer=1, turn_into_array=1, ignore_2nd_column=1, matrix_data_type=float)
		output_dir = '/Network/Data/250k/tmp-yh/CNV/amp_hist/'
		CNV.drawCNVAmpHist(snpData, output_dir, max_counter=1000)
	"""
	
	@classmethod
	def outputCNVMatrixIntoMemmapFormat(cls, input_fname, output_fname):
		"""
		2009-9-27
			input format is output of CNVNormalize.py & DB_250k2Array.py
			output is binary, probe X sample. test if matlab sees it as a memmapfile.
		"""
		from variation.src.CNVNormalize import CNVNormalize
		data_matrix, probe_id_ls, chr_pos_ls, header = CNVNormalize.get_input(input_fname)
		outf = open(output_fname, 'wb')
		import numpy
		data_matrix = numpy.transpose(data_matrix)	#somehow, without this, matlab will read it as transposed to the original matrix 
		data_matrix.tofile(outf, format='%.8f')	# format seems to not matter.
		"""
		no_of_rows, no_of_cols = data_matrix.shape
		for i in range(no_of_rows):
			for j in range(no_of_cols):
				outf.write(struct.pack('d', data_matrix[i,j])
		"""
		del outf
	"""
	input_fname = os.path.expanduser('~/panfs/250k/CNV/call_method_17_CNV_array_inensity_norm_chr4_n101.tsv')
	output_fname = os.path.expanduser('~/panfs/250k/CNV/call_method_17_CNV_array_intensity_norm_chr4_n101.memmap')
	CNV.outputCNVMatrixIntoMemmapFormat(input_fname, output_fname)
	
	for i in range(2,6):
		input_fname = os.path.expanduser('~/panfs/250k/CNV/call_43_CNV_norm_intensity_chr%s.tsv'%i)
		output_fname = os.path.expanduser('~/panfs/250k/CNV/call_43_CNV_norm_intensity_chr%s.memmap'%i)
		CNV.outputCNVMatrixIntoMemmapFormat(input_fname, output_fname)
	"""
	
	@classmethod
	def outputSelectedArraysInGWAFormat(cls, input_fname, output_dir, array_id_ls):
		"""
		2009-10-11
			input_fname is DB_250k2Array.py's output.
			output the intensity from specified arrays, one array one file to inspect genome-wide pattern:
				chromosome	position	amplitude 
		"""
		sys.stderr.write("Outputting selected arrays in GWA format ...\n")
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		
		import csv
		from pymodule import figureOutDelimiter
		reader = csv.reader(open(input_fname), delimiter=figureOutDelimiter(input_fname))
		
		
		array_id_set = set(array_id_ls)
		header = reader.next()
		array_index2writer = {}
		for i in range(1, len(header)-2):
			array_id = int(header[i])
			if array_id in array_id_set:
				output_fname = os.path.join(output_dir, 'array_%s_CNV.tsv'%array_id)
				writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
				array_index2writer[i] = writer
		
		counter = 0
		for row in reader:
			chr, pos = row[-2:]
			for i in range(1, len(row)-2):
				if i in array_index2writer:
					new_row = [chr, pos, row[i]]
					array_index2writer[i].writerow(new_row)
			counter += 1
			if counter%10000==0:
				sys.stderr.write("%s%s"%('\x08'*40, counter))
		del reader, array_index2writer
		sys.stderr.write("Done.\n")
	
	
	"""
	Col_array_id_ls = [1, 2, 43, 139, 145]	# Col-0 arrays
	Ler_array_id_ls = [3, 4, 41, 150, 151,]
	array_id_ls = Col_array_id_ls + Ler_array_id_ls
	
	input_fname = os.path.expanduser('~/mnt2/panfs/250k/CNV/call_48_CNV_QNormalize.tsv')
	output_dir = os.path.expanduser('~/mnt/banyan/tmp/call_48_CNV_QNormalize')
	CNV.outputSelectedArraysInGWAFormat(input_fname, output_dir, array_id_ls)
	
	
	input_fname = os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity.tsv')
	output_dir = os.path.expanduser('~/mnt/banyan/tmp/call_48_CNV')
	CNV.outputSelectedArraysInGWAFormat(input_fname, output_dir, array_id_ls)
	
	input_fname = os.path.expanduser('~/mnt2/panfs/250k/CNV/call_43_CNV_intensity.tsv')
	output_dir = os.path.expanduser('~/mnt/banyan/tmp/call_43_CNV')
	CNV.outputSelectedArraysInGWAFormat(input_fname, output_dir, array_id_ls)
	
	
	for i in range(1,6):
		input_fname = os.path.expanduser('~/mnt2/panfs/250k/CNV/call_43_CNV_norm_intensity_chr%s.tsv'%i)
		output_dir = os.path.expanduser('~/mnt/banyan/tmp/call_43_CNV_norm_intensity_chr%s'%i)
		CNV.outputSelectedArraysInGWAFormat(input_fname, output_dir, array_id_ls)
	"""
	
	@classmethod
	def outputMedianIntensityOfSelectedArraysInGWAFormat(cls, input_fname, output_fname, array_id_ls):
		"""
		2009-10-27
			similar to outputSelectedArraysInGWAFormat() but take median of probes across all arrays
			input_fname is DB_250k2Array.py's output.
			output the intensity from specified arrays, one array one file to inspect genome-wide pattern:
				chromosome	position	amplitude 
		"""
		sys.stderr.write("Outputting selected arrays in GWA format ...\n")
		
		import csv, numpy
		from pymodule import figureOutDelimiter
		reader = csv.reader(open(input_fname), delimiter=figureOutDelimiter(input_fname))
		
		
		array_id_set = set(array_id_ls)
		header = reader.next()
		array_index2writer = {}
		for i in range(1, len(header)-2):
			array_id = int(header[i])
			if array_id in array_id_set:
				#output_fname = os.path.join(output_dir, 'array_%s_CNV.tsv'%array_id)
				#writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
				array_index2writer[i] = None
		
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		counter = 0
		for row in reader:
			chr, pos = row[-2:]
			intensity_ls = []
			for i in range(1, len(row)-2):
				if i in array_index2writer:
					intensity_ls.append(float(row[i]))
			median = numpy.median(intensity_ls)
			new_row = [chr, pos, median]
			writer.writerow(new_row)
			counter += 1
			if counter%10000==0:
				sys.stderr.write("%s%s"%('\x08'*40, counter))
		del reader, array_index2writer, writer
		sys.stderr.write("Done.\n")
		
	"""
	Col_array_id_ls = [1, 2, 43, 139, 145]	# Col-0 arrays
	Ler_array_id_ls = [3, 4, 41, 150, 151,]
	
	input_fname = os.path.expanduser('~/mnt2/panfs/250k/CNV/call_48_CNV_QNormalize.tsv')
	output_fname = os.path.expanduser('~/mnt/banyan/tmp/call_48_CNV_QNormalize_Col_median.tsv')
	CNV.outputMedianIntensityOfSelectedArraysInGWAFormat(input_fname, output_fname, Col_array_id_ls)
	output_fname = os.path.expanduser('~/mnt/banyan/tmp/call_48_CNV_QNormalize_Ler_median.tsv')
	CNV.outputMedianIntensityOfSelectedArraysInGWAFormat(input_fname, output_fname, Ler_array_id_ls)
	
	"""
	
	@classmethod
	def pickCNVsForValidation(cls, db_250k, cnv_method_id=20, cnv_type_id=1, no_of_loci=100):
		"""
		2011-4-27
			select deletions based on frequency, size, probability
			output:
				accession name, deletion location, size, sequence deleted, 2-3kb flanking sequence
				
				which base in the flanking sequence is a SNP (either 250k or perlegen)
				
				avoid region with excessive polymorphism
		"""
		sys.stderr.write("Picking %s type %s CNVs from method %s for validation ..."%(no_of_loci, cnv_type_id, cnv_method_id))
		import Stock_250kDB
		query = Stock_250kDB.CNVArrayCall.query.filter_by(cnv_type_id=cnv_type_id).filter_by(cnv_method_id=cnv_method_id)
		
		sys.stderr.write("Done.\n")
		
	
	@classmethod
	def putCNVIntoSnps(cls, db_250k, cnv_method_id=None):
		"""
		2011-4-22
			put CNVs from particular cnv_method_id into table Snps
				so that results_method could just point to call_method, rather than cnv_method
		"""
		import Stock_250kDB
		sys.stderr.write("Putting cnvs from method %s into db ...\n"%(cnv_method_id))
		query = Stock_250kDB.CNV.query.filter_by(cnv_method_id=cnv_method_id)
		counter = 0
		real_counter = 0
		for row in query:
			counter += 1
			db_entry = db_250k.getSNP(chromosome=row.chromosome, start=row.start, stop=row.stop)
			if db_entry.id is None:
				real_counter += 1
			if counter%5000==0:
				sys.stderr.write("%s\t%s\t%s"%('\x08'*80, counter, real_counter))
		db_250k.session.flush()
		db_250k.session.commit()
		sys.stderr.write("%s\t%s\t%s.\n"%('\x08'*80, counter, real_counter))
	
	"""
		#2011-4-22
		cnv_method_id=20
		CNV.putCNVIntoSnps(db_250k, cnv_method_id=cnv_method_id)
		sys.exit(0)
		
		
	"""
	
	@classmethod
	def outputNonOverlappingCNVAsSNP(cls, db_250k, output_fname, cnv_method_id=None, cnv_type_id=None):
		"""
		2011-4-24
			replace CNV.id with Snps.id (same chr,start,stop).
			you can run CNV.putCNVIntoSnps() to make sure every CNV from that method is in db. 
		2011-2-17
			use CNV.id as column ID
		2010-8-7
			in matrix, normal is represented as 0. deletion or other chosen cnv type -> 1.
		"""
		sys.stderr.write("Outputting non-overlapping CNVs as SNPs for method %s, type %s ...\n"%(cnv_method_id, \
																			cnv_type_id))
		col_id_ls = []
		col_id2index = {}
		cnv_id2index = {}
		import Stock_250kDB
		query = Stock_250kDB.CNV.query.filter_by(cnv_method_id=cnv_method_id).filter_by(cnv_type_id=cnv_type_id).\
			order_by(Stock_250kDB.CNV.chromosome).order_by(Stock_250kDB.CNV.start).order_by(Stock_250kDB.CNV.stop)
		
		for row in query:
			#col_id = '%s_%s_%s'%(row.chromosome, row.start, row.stop)
			db_entry = db_250k.getSNP(chromosome=row.chromosome, start=row.start, stop=row.stop)
			if db_entry.id is None:
				db_entry.session.flush()
			col_id = db_entry.id		#2011-2-17 use CNV.id as column ID
			col_id_ls.append(col_id)
			col_id2index[col_id] = len(col_id2index)
			cnv_id2index[row.id] = col_id2index[col_id]
		sys.stderr.write("%s unique CNVs.\n"%(len(col_id2index)))
		
		row_id_ls = []
		row_id2index = {}
		query = Stock_250kDB.CNVArrayCall.query.filter_by(cnv_method_id=cnv_method_id)
		
		i = 0
		block_size = 10000
		real_counter = 0
		
		rows = query.offset(i).limit(block_size)
		session = db_250k.session
		data_matrix = []
		while rows.count()!=0:
			for row in rows:
				row_id = (row.array.maternal_ecotype_id, row.array_id)
				if row_id not in row_id2index:
					row_id2index[row_id] = len(row_id2index)
					row_id_ls.append(row_id)
					data_row = [0]*len(col_id_ls)
					data_matrix.append(data_row)
				row_index = row_id2index[row_id]
				data_row = data_matrix[row_index]
				col_index = cnv_id2index[row.cnv_id]
				data_row[col_index] = 1
				i += 1
			if i%5000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, i, real_counter))
			rows = query.offset(i).limit(block_size)
		sys.stderr.write("%s\t%s\t%s Done.\n"%('\x08'*80, i, real_counter))
		
		from pymodule.SNP import write_data_matrix
		header = ['ecotype_id', 'array_id'] + col_id_ls
		strain_acc_list = [row[0] for row in row_id_ls]
		category_list = [row[1] for row in row_id_ls]
		write_data_matrix(data_matrix, output_fname, header, strain_acc_list, category_list, rows_to_be_tossed_out=None, \
			cols_to_be_tossed_out=None, nt_alphabet=0, transform_to_numpy=0,\
			discard_all_NA_rows=0, strain_acc2other_info=None, delimiter='\t', )
	"""
		#2010-8-7
		CNV.outputNonOverlappingCNVAsSNP(db_250k, output_fname, cnv_method_id=cnv_method_id, cnv_type_id=1)
		sys.exit(0)
		
		#2010-8-7
		cnv_method_id = 20
		output_fname = os.path.expanduser('~/script/variation/data/CNV/NonOverlapCNVAsSNP_cnvMethod%s.tsv'%cnv_method_id)
		CNV.outputNonOverlappingCNVAsSNP(db_250k, output_fname, cnv_method_id=cnv_method_id, cnv_type_id=1)
		sys.exit(0)
		
	"""
		
	@classmethod
	def outputCNVFunctionalCompositionData(cls, db_250k, cnv_method_id=20, cnv_type_id=1):
		"""
		2010-10-24
		"""
		import Stock_250kDB
		from variation.src.DrawSNPRegion import DrawSNPRegion
		def dealWithGeneAnnotation():
			gene_annotation_picklef = '/Network/Data/250k/tmp-yh/at_gene_model_pickelf'
			DrawSNPRegion_ins = DrawSNPRegion(db_user=db_250k.username, db_passwd=db_250k.password, hostname=db_250k.hostname, \
										database=db_250k.database,\
										input_fname='/tmp/dumb', output_dir='/tmp', debug=0)
			gene_annotation = DrawSNPRegion_ins.dealWithGeneAnnotation(gene_annotation_picklef, cls_with_db_args=DrawSNPRegion_ins)
			return gene_annotation
		gene_annotation = dealWithGeneAnnotation()
		
		gene_type_id2gene_set = {}
		gene_id2type_id = {}
		for gene_id, model in gene_annotation.gene_id2model.iteritems():
			type_id = model.type_id
			if type_id==9 or type_id==11:	#turn TE gene into TE type
				type_id=11
			elif type_id==1:
				type_id=type_id
			else:
				type_id = 12
			if type_id not in gene_type_id2gene_set:
				gene_type_id2gene_set[type_id] = set()
			gene_type_id2gene_set[type_id].add(gene_id)
			gene_id2type_id[gene_id] = type_id
		
		sys.stderr.write("Outputting the number of genes in each type ... \n")
		for gene_type_id, gene_set in gene_type_id2gene_set.iteritems():
			print gene_type_id, len(gene_set)
		
		sys.stderr.write("Getting cnv.id from CNV ...")
		Table = Stock_250kDB.CNV
		query = Table.query.filter_by(cnv_method_id=cnv_method_id).filter_by(cnv_type_id=cnv_type_id)
		query = query.filter(Table.frequency<0.1)
		total_cnv_id_set = set()
		for row in query:
			total_cnv_id_set.add(row.id)
		sys.stderr.write("%s total cnvs under this condition Done.\n"%(len(total_cnv_id_set)))
		
		sys.stderr.write("Getting stuff from CNVContext ...")
		Table = Stock_250kDB.CNVContext
		query = Table.query.filter(Table.cnv.has(cnv_method_id=cnv_method_id)).filter(Table.cnv.has(cnv_type_id=cnv_type_id)).\
			filter(Table.cnv.has(Stock_250kDB.CNV.frequency<0.1))
		query = query.filter(Table.overlap_length>50)
		cnv_in_gene_context_id_set = set()
		gene_type_id2cnv_id_set = {}
		gene_type_id2gene_within_cnv_set = {}
		for row in query:
			cnv_in_gene_context_id_set.add(row.cnv_id)
			type_id = gene_id2type_id.get(row.gene_id)
			if type_id not in gene_type_id2cnv_id_set:
				gene_type_id2gene_within_cnv_set[type_id] = set()
				gene_type_id2cnv_id_set[type_id] = set()
			gene_type_id2cnv_id_set[type_id].add(row.cnv_id)
			gene_type_id2gene_within_cnv_set[type_id].add(row.gene_id)
		sys.stderr.write("Done.\n")
		
		print "Number of distinct cnv ids in any gene:", len(cnv_in_gene_context_id_set)
		print "Number of distinct cnv ids in each type of genes."
		for gene_type_id, cnv_id_set in gene_type_id2cnv_id_set.iteritems():
			print gene_type_id, len(cnv_id_set)
		
		print "Number of genes of each type in all these cnvs"
		for gene_type_id, gene_within_cnv_set in gene_type_id2gene_within_cnv_set.iteritems():
			print gene_type_id, len(gene_within_cnv_set)
		
	"""
		#2010-10-24
		CNV.outputCNVFunctionalCompositionData(db_250k, cnv_method_id=20, cnv_type_id=1)
		sys.exit(0)
	"""
	
	@classmethod
	def compareCNVSegmentsAgainstQCHandler(cls, input_fname_ls, ecotype_id2cnv_qc_call_data=None, afterMatchFunction=None, \
										param_obj=None, deletion_cutoff=None, max_boundary_diff=10000, \
										max_diff_perc=0.10, min_no_of_probes=5, count_embedded_segment_as_match=False, \
										min_reciprocal_overlap=0.6, functor_after_span=None, judgeTrueDeletionFunctor=None,\
										report=True):
		"""
		2010-6-17
			rename argument function_handler to afterMatchFunction
			add judgeTrueDeletionFunctor, default is isDeletionTrueBasedOnOtherDeletionData
		2010-1-26
			value of the ecotype_id2cnv_qc_call_data dictionary is a RBDict (RBTree dictionary) structure.
		2009-12-8
			add argument min_reciprocal_overlap
		2009-11-4
			a general handler to compare CNV segments from input_fname_ls with cnv_qc_call_data.
			Upon a match between a CNV segment from input_fname_ls and cnv_qc_call_data, function_handler would be called with param_obj as argument.
			
			If deletion_cutoff is None, all segments who have matches in ecotype_id2cnv_qc_call_data would be considered in function_handler().
			If deletion_cutoff is not None (some float), only those segments whose amplitude is below this value would be considered in function_handler().
		"""
		import fileinput
		from pymodule import getColName2IndexFromHeader, PassingData
		from pymodule.CNV import is_reciprocal_overlap, CNVSegmentBinarySearchTreeKey
		if ecotype_id2cnv_qc_call_data is None:
			no_of_ecotypes = 0
		else:
			no_of_ecotypes = len(ecotype_id2cnv_qc_call_data)
		sys.stderr.write("Running function on RunGADA.py output: %s for %s ecotypes ... \n"%\
						(repr(input_fname_ls), no_of_ecotypes))
		amp_ls = []
		array_id2array = {}
		counter = 0
		real_counter = 0
		no_of_deletions = 0
		no_of_valid_deletions = 0
		input_handler = fileinput.input(input_fname_ls)
		header = input_handler.readline().strip().split('\t')
		col_name2index = getColName2IndexFromHeader(header)
		ecotype_id2span_data = getattr(param_obj, 'ecotype_id2span_data', None)
		median_col_index = col_name2index.get('median')
		for line in input_handler:
			if line.find("array_id")!=-1:
				continue
			line = line.strip()
			row = line.split('\t')
			ecotype_id_idx = col_name2index.get('ecotype_id', col_name2index.get('array_id'))
			cnv_ecotype_id = int(row[ecotype_id_idx])
			array_id = int(row[col_name2index.get('array_id')])
			#row[ecotype_id_idx] = cnv_ecotype_id
			counter += 1
			if ecotype_id2cnv_qc_call_data is None or \
					(ecotype_id2cnv_qc_call_data is not None and cnv_ecotype_id in ecotype_id2cnv_qc_call_data):	# array is in CNVQCDat
				start_probe = row[col_name2index['start_probe']].split('_')	# split chr_pos
				start_probe = map(int, start_probe)
				start_probe_id = row[col_name2index['start_probe_id']]
				stop_probe = row[col_name2index['end_probe']].split('_')
				stop_probe = map(int, stop_probe)
				stop_probe_id = row[col_name2index['end_probe_id']]
				no_of_probes = int(row[col_name2index['length']])
				if no_of_probes<min_no_of_probes:
					continue
				amplitude = float(row[col_name2index['amplitude']])
				segment_chromosome = start_probe[0]
				if start_probe[0]!=stop_probe[0]:	#spurious. on different chromosomes.
					continue
				segment_start_pos = start_probe[1]-12
				segment_stop_pos = stop_probe[1]+12
				segment_length = abs(segment_stop_pos-segment_start_pos+1)
				if deletion_cutoff is not None and amplitude>deletion_cutoff:
					continue
				if median_col_index is not None:
					median_intensity = float(row[median_col_index])
				else:
					median_intensity = None
				
				cnv_segment_obj = PassingData(ecotype_id=cnv_ecotype_id, start_probe=start_probe, stop_probe=stop_probe,\
											no_of_probes=no_of_probes, amplitude=amplitude, segment_length=segment_length,\
											segment_chromosome=segment_chromosome, array_id=array_id,\
											start_probe_id=start_probe_id, stop_probe_id=stop_probe_id,\
											segment_start_pos=segment_start_pos, segment_stop_pos=segment_stop_pos,\
											median_intensity=median_intensity)
				cnvSegmentKey = CNVSegmentBinarySearchTreeKey(chromosome=segment_chromosome, span_ls=[segment_start_pos, segment_stop_pos],\
															min_reciprocal_overlap=min_reciprocal_overlap)
				
				no_of_deletions+=1
				if ecotype_id2span_data:
					span_data = ecotype_id2span_data.get(cnv_ecotype_id)
				else:
					span_data = None
				if span_data and cnvSegmentKey not in span_data:	# 2010-1-28 skip this CNV if span_data is available and it doesn't include this CNV. (to get good FPR)
					continue
				
				if span_data and functor_after_span:	# 2010-2-10
					functor_after_span(cnv_segment_obj, param_obj)
				
				#2010-6-17
				if judgeTrueDeletionFunctor is None:
					judgeTrueDeletionFunctor = cls.isDeletionTrueBasedOnOtherDeletionData
				
				param_obj.cnvSegmentKey = cnvSegmentKey
				judgeData = judgeTrueDeletionFunctor(cnv_segment_obj, param_obj)
				param_obj.isDeletionTrue = judgeData.isDeletionTrue
				if judgeData and judgeData.isDeletionTrue:
					no_of_valid_deletions += 1
					afterMatchFunction(param_obj, cnv_segment_obj, cnv_qc_call=getattr(judgeData, 'cnv_qc_call', None))
				else:
					# this stage, call the handler if it wants to record the number of deletions.
					afterMatchFunction(param_obj, cnv_segment_obj, cnv_qc_call=None)
				"""
				for cnv_qc_call in cnv_qc_call_data:
					qc_chromosome, qc_start, qc_stop = cnv_qc_call[:3]
					cnv_qc_call_id = cnv_qc_call[-1]
					valid_match = False
					if qc_chromosome==segment_chromosome:
						boundary_diff1 = abs(segment_start_pos-qc_start)
						boundary_diff2 = abs(segment_stop_pos-qc_stop)
						diff1_perc = boundary_diff1/float(segment_length)
						diff2_perc = boundary_diff2/float(segment_length)
						
						is_overlap = is_reciprocal_overlap([segment_start_pos, segment_stop_pos], [qc_start, qc_stop], \
													min_reciprocal_overlap=min_reciprocal_overlap)
						
						if is_overlap:
							no_of_valid_deletions += 1
							valid_match = True
						
						#if boundary_diff1<=max_boundary_diff and boundary_diff2<=max_boundary_diff and diff1_perc<=max_diff_perc and \
						#diff2_perc<=max_diff_perc:
						#	no_of_valid_deletions += 1
						#	valid_match = True
						#elif count_embedded_segment_as_match and segment_start_pos>=qc_start and segment_stop_pos<=qc_stop:	#the segment doesn't match the criteria but very small and within
						#	no_of_valid_deletions += 1
						#	valid_match = True
						
						if valid_match:
							afterMatchFunction(param_obj, cnv_segment_obj, cnv_qc_call, )
					elif qc_chromosome>segment_chromosome:
						break
				"""
			if report and counter%10000==0:
				sys.stderr.write('%s%s\t%s\t%s'%('\x08'*80, counter, no_of_deletions, no_of_valid_deletions))
		setattr(param_obj, "no_of_deletions", no_of_deletions)
		setattr(param_obj, "no_of_valid_deletions", no_of_valid_deletions)
		sys.stderr.write("\n")
	
	@classmethod
	def comparePeaksFromTwoAssociationResults(cls, db_250k, result1_id=None, result1_peak_type_id=None,\
											result2_id=None, result2_peak_type_id=None, result2_peak_ext_dist=0):
		"""
		2011-4-22
			peak data is from table ResultPeak.
			The function tells how many peaks from result1 are not found, and how many are found
				in result2.
		"""
		sys.stderr.write("Comparing peaks from result %s against result %s ...\n"%(result1_id, result2_id))
		import Stock_250kDB
		from pymodule.CNV import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
		from pymodule.RBTree import RBDict
		
		sys.stderr.write("Constructing RBDict for peaks from result %s ..."%(result2_id))
		result2_peakRBDict = RBDict()
		query = Stock_250kDB.ResultPeak.query.filter_by(result_id=result2_id).filter_by(result_peak_type_id=result2_peak_type_id)
		counter = 0
		real_counter = 0
		for row in query:
			counter += 1
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=row.chromosome, \
							span_ls=[max(1, row.start - result2_peak_ext_dist), row.stop + result2_peak_ext_dist], \
							min_reciprocal_overlap=1,)
							#2010-8-17 overlapping keys are regarded as separate instances as long as they are not identical.
			if segmentKey not in result2_peakRBDict:
				result2_peakRBDict[segmentKey] = []
			result2_peakRBDict[segmentKey].append(row)
		sys.stderr.write("%s peaks. Done.\n"%counter)
		
		compareIns = CNVCompare(min_reciprocal_overlap=0.0000001)	#any overlap is an overlap
		query = Stock_250kDB.ResultPeak.query.filter_by(result_id=result1_id).filter_by(result_peak_type_id=result1_peak_type_id)
		no_of_peaks_not_in_result2 = 0
		overlap_ls = []
		score_ls = []
		counter = 0
		for row in query:
			counter += 1
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=row.chromosome, \
							span_ls=[row.start, row.stop], \
							min_reciprocal_overlap=0.0000001, )	#min_reciprocal_overlap doesn't matter here.
				# it's decided by compareIns.
			node_ls = []
			result2_peakRBDict.findNodes(segmentKey, node_ls=node_ls, compareIns=compareIns)
			total_perc_overlapped_by_result2 = 0.
			for node in node_ls:
				overlap1, overlap2, overlap_length, overlap_start_pos, overlap_stop_pos = get_overlap_ratio(segmentKey.span_ls, \
										[node.key.start, node.key.stop])[:5]
				total_perc_overlapped_by_result2 += overlap1
			score_ls.append(row.peak_score)
			if total_perc_overlapped_by_result2==0:
				no_of_peaks_not_in_result2 += 1
				overlap_ls.append(-0.5)
			else:
				overlap_ls.append(total_perc_overlapped_by_result2)
		
		sys.stderr.write("%s/%s peaks in result %s not found in result %s. Done.\n"%(no_of_peaks_not_in_result2,\
												counter, result1_id, result2_id))
		from pymodule import yh_matplotlib
		yh_matplotlib.drawHist(overlap_ls, title='Histogram of overlap rate of %s peaks'%(len(overlap_ls)), \
							xlabel_1D='overlap rate',\
							outputFname='/tmp/hist_of_result%s_overlap_rate_in_result%s.png'%(result1_id, result2_id), \
							min_no_of_data_points=10, needLog=False)
		outputFname = '/tmp/hist_of_result%s_overlap_rate_in_result%s_2D.png'%(result1_id, result2_id)
		C_ls = [1]*len(overlap_ls)
		colorBarLabel='log(count)'
		reduce_C_function = yh_matplotlib.logSum
		yh_matplotlib.drawHexbin(overlap_ls, score_ls, C_ls, fig_fname=outputFname, gridsize=10, \
								title='%s peaks in result %s vs %s'%(counter, result1_id, result2_id), \
								xlabel = 'overlap rate', \
								ylabel = 'association-score',\
								colorBarLabel=colorBarLabel, reduce_C_function= reduce_C_function)
	"""
		# 2011-4-22
		result1_id = 4634	#cnv_20_LD_KW
		result1_peak_type_id = 1
		result2_id = 3395	#call_32_LD_KW
		result2_peak_type_id = 2
		CNV.comparePeaksFromTwoAssociationResults(db_250k, result1_id=result1_id, result1_peak_type_id=result1_peak_type_id,\
										result2_id=result2_id, result2_peak_type_id=result2_peak_type_id, result2_peak_ext_dist=0)
		sys.exit(0)
		
		# 2011-4-22
		result1_id = 3395	#call_32_LD_KW
		result1_peak_type_id = 2 #min_score=5
		result2_id = 3396	#call_32_LDV_KW
		result2_peak_type_id = 2	#min_score=5
		CNV.comparePeaksFromTwoAssociationResults(db_250k, result1_id=result1_id, result1_peak_type_id=result1_peak_type_id,\
										result2_id=result2_id, result2_peak_type_id=result2_peak_type_id, result2_peak_ext_dist=0)
		sys.exit(0)
		
		result1_id_ls = [4634, 4635, 4636, 4637]	#cnv_20_LD_KW, LDV, SD, SDV
		result1_peak_type_id = 1 #min_score=4
		result2_id_ls = [3395, 3396, 3397, 3398]	#call_32_LD_KW, LDV, SD, SDV
		result2_id = 3396	#call_32_LDV_KW
		result2_peak_type_id = 2	#min_score=5
		for i in xrange(len(result1_id_ls)):
			result1_id = result1_id_ls[i]
			result2_id = result2_id_ls[i]
		CNV.comparePeaksFromTwoAssociationResults(db_250k, result1_id=result1_id, result1_peak_type_id=result1_peak_type_id,\						result2_id=result2_id, result2_peak_type_id=result2_peak_type_id, result2_peak_ext_dist=0)
		sys.exit(0)
	"""
	
	@classmethod
	def putOneCNVSegmentIntoDBFunctor(cls, param_obj, cnv_segment_obj=None, cnv_qc_call=None):
		"""
		2010-7-1
			moved to CNVPredictDeletionBySVM.saveSegmentObj() of CNVPredictionBySVM.py
		2010-3-16
		"""
		from CNVPredictDeletionBySVM import CNVPredictDeletionBySVM
		CNVPredictDeletionBySVM.saveSegmentObj(param_obj, cnv_segment_obj)
	
	@classmethod
	def putCNVSegmentIntoDB(cls, db_250k, input_fname_ls, ecotype_id_set=None, cnv_type_id=1, cnv_method_id=6,\
								deletion_cutoff=-0.33, \
								min_no_of_probes=5, report=True):
		"""
		2010-3-16
			put any GADA segment whose amplitude is below the deletion_cutoff into db
		"""
		from pymodule import PassingData
		if ecotype_id_set is not None:
			ecotype_id2cnv_qc_call_data = {}
			for ecotype_id in ecotype_id_set:
				ecotype_id2cnv_qc_call_data[ecotype_id] = None
		else:
			ecotype_id2cnv_qc_call_data = None	# no filtering. all ecotypes into db.
		param_obj = PassingData(no_of_valid_deletions=0, session=db_250k.session, cnv_type_id=cnv_type_id, \
							cnv_method_id=cnv_method_id, no_of_total=0, no_of_into_db=0, report=report)
		cls.compareCNVSegmentsAgainstQCHandler(input_fname_ls, ecotype_id2cnv_qc_call_data=ecotype_id2cnv_qc_call_data, 
											afterMatchFunction=cls.putOneCNVSegmentIntoDBFunctor, param_obj=param_obj, \
											deletion_cutoff=deletion_cutoff, min_no_of_probes=min_no_of_probes, \
											min_reciprocal_overlap=None, report=report, \
											functor_after_span=None)
		ratio = param_obj.no_of_into_db/float(param_obj.no_of_total)
		sys.stderr.write("%s out of %s (%s) into db.\n"%(param_obj.no_of_into_db, param_obj.no_of_total, ratio))
		db_250k.session.flush()
		db_250k.session.commit()
	
	"""
		2010-3-16
		input_fname_ls = []
		for i in range(1,6):
			input_fname_ls.append(os.path.expanduser('~/mnt/panfs/250k/CNV/call_48_CNV_QNormalize_chr%s.GADA_A0.5T4M5.tsv'%i))
		
		deletion_cutoff = 2.3
		min_no_of_probes = 10
		CNV.putCNVSegmentIntoDB(db_250k, input_fname_ls, ecotype_id_set=set([6932, 6909, 8215, 6911, 6977, 6962]), cnv_type_id=1, cnv_method_id=6,\
									deletion_cutoff=deletion_cutoff, \
									min_no_of_probes=min_no_of_probes)
		
		# 2010-6-27
		input_fname_ls = [os.path.expanduser('~/script/variation/data/CNV/call_53_WABQN_b100_j50_lts_98_arrays_GADA_A0.5T4.0M5.tsv')]
		deletion_cutoff = -0.1
		min_no_of_probes = 5
		CNV.putCNVSegmentIntoDB(db_250k, input_fname_ls, ecotype_id_set=None, cnv_type_id=1, cnv_method_id=6,\
									deletion_cutoff=deletion_cutoff, \
									min_no_of_probes=min_no_of_probes)
		sys.exit(0)
		
		# 2010-6-28
		input_fname_ls = [os.path.expanduser('~/script/variation/data/CNV/call_53_WABQN_b100_j50_lts_98_arrays_2Flanking_GADA_A0.5T8.0M5.tsv')]
		deletion_cutoff = -0.1
		min_no_of_probes = 5
		cnv_method_id = 8
		CNV.putCNVSegmentIntoDB(db_250k, input_fname_ls, ecotype_id_set=set([6932]), cnv_type_id=1, cnv_method_id=cnv_method_id,\
									deletion_cutoff=deletion_cutoff, \
									min_no_of_probes=min_no_of_probes)
		sys.exit(0)
		
		# 2010-7-16 for the TAIR9 "-z banyan.usc.edu"
		input_fname_ls = [os.path.expanduser('~/mnt/hpc-cmb/script/variation/data/CNV/call_53_WABQN_b100_j50_lts_TAIR9_4Flanking_GADA_A0.5T12.0M5.tsv')]
		deletion_cutoff = -0.1
		min_no_of_probes = 5
		cnv_method_id = 7
		CNV.putCNVSegmentIntoDB(db_250k, input_fname_ls, ecotype_id_set=set([6932]), cnv_type_id=1, cnv_method_id=cnv_method_id,\
									deletion_cutoff=deletion_cutoff, \
									min_no_of_probes=min_no_of_probes)
		sys.exit(0)
		
		# 2010-7-19 "-z papaya.usc.edu"
		input_fname_ls = [os.path.expanduser('~/mnt/hpc-cmb/script/variation/data/CNV/call_53_WABQN_b100_j50_lts_4Flanking_GADA_A0.5T12.0M5.tsv')]
		deletion_cutoff = -0.1
		min_no_of_probes = 5
		cnv_method_id = 11
		CNV.putCNVSegmentIntoDB(db_250k, input_fname_ls, ecotype_id_set=set([6932]), cnv_type_id=1, cnv_method_id=cnv_method_id,\
									deletion_cutoff=deletion_cutoff, \
									min_no_of_probes=min_no_of_probes)
		sys.exit(0)
		
		# 2010-7-16 for the TAIR9 "-z banyan.usc.edu"
		ecotype_id_ls = CNV.getEcotypeIDListInCNVQCAccessionGivenDataSourceID(data_source_id=13)
		input_fname_ls = [os.path.expanduser('~/mnt/hpc-cmb/script/variation/data/CNV/call_53_WABQN_b100_j50_lts_TAIR9_4Flanking_GADA_A0.5T12.0M5.tsv')]
		deletion_cutoff = -0.09
		min_no_of_probes = 5
		cnv_method_id = 7
		CNV.putCNVSegmentIntoDB(db_250k, input_fname_ls, ecotype_id_set=set(ecotype_id_ls), \
							cnv_type_id=1, cnv_method_id=cnv_method_id,\
							deletion_cutoff=deletion_cutoff, \
							min_no_of_probes=min_no_of_probes)
		sys.exit(0)
	"""
	
	@classmethod
	def getTAIR8ToTAIR9_in_chr_pos(cls, db_250k, probeType=2, debug=False):
		"""
		2010-8-3
			return a map between tair8 and tair9
		"""
		sys.stderr.write("Getting TAIR8 to TAIR9 translation for probeType %s ...\n"%probeType)
		tair8_chr_pos2tair9_chr_pos = {}
		import Stock_250kDB
		probes_table = Stock_250kDB.Probes.table.name
		
		sql_query = 'select id, snps_id, chromosome, position, tair8_chromosome, tair8_position from %s '%(probes_table)
		where_sql_ls = []
		order_by_sql = ""
		if probeType==2 or probeType==4:
			where_sql_ls.append("snps_id is null and chromosome is not null and direction is not null and Tair9Copy=1")
			order_by_sql = "order by chromosome, position"
			# 2010-5-5
		elif probeType==1:	# SNP Probes
			where_sql_ls.append("snps_id is not null")
		elif probeType == 5:	# QC probes
			where_sql_ls.append("snps_id is null and chromosome is null and position is null")
		else:	#2010-5-5 get all probes
			order_by_sql = "order by chromosome, position"
		
		if where_sql_ls:
			where_sql = "where %s"%(' and '.join(where_sql_ls))
		else:
			where_sql = ""
		
		rows = db_250k.metadata.bind.execute("%s %s %s"%(sql_query, where_sql, order_by_sql))
		
		counter = 0
		for row in rows:
			probe_id = row.id
			snps_id = row.snps_id
			chromosome = row.chromosome
			position = row.position
			tair8_chr_pos = (row.tair8_chromosome, row.tair8_position)
			tair9_chr_pos = (row.chromosome, row.position)
			
			tair8_chr_pos2tair9_chr_pos[tair8_chr_pos] = tair9_chr_pos
			
			counter += 1
			if counter%5000==0:
				sys.stderr.write('%s%s'%('\x08'*80, counter,))
				if debug:
					break
		del rows
		sys.stderr.write("\t %s probes. Done.\n"%(counter))
		return tair8_chr_pos2tair9_chr_pos
		
		
	@classmethod
	def putSebastianCallIntoDB(cls, db_250k, input_dir, minProbabilityToCallDeletion=0.8,\
							cnv_method_id=22, cnv_type_id=1, debug=False):
		"""
		2010-8-3
			File is space-separated matrix.
			
			first line is header for array-id.
			2nd line and onwards look like:
				analysis_renorm_r30_m10_P0.001_c1.0.txt:CNV 1 type 2 start 1034 position 76629 length 17 end position 78133 frequency 0.184 \
				Carrier status: 0.00 0.00 0.00 0.00 0.00 0.14 0.00 0.00 0.00 0.00 0.00 0
		"""
		sys.stderr.write("Putting Sebastian calls into db ...\n")
		
		tair8_chr_pos2tair9_chr_pos = cls.getTAIR8ToTAIR9_in_chr_pos(db_250k, probeType=2, debug=debug)
		
		import csv, Stock_250kDB
		input_fname_ls = os.listdir(input_dir)
		import re
		cnvPosPattern = re.compile(r'position (\d+) length (\d+) end position (\d+)')
		session = db_250k.session
		counter = 0
		real_counter = 0
		no_of_cnvs_skipped = 0
		for input_fname in input_fname_ls:
			if input_fname[-3:]=='txt':
				chromosome = int(input_fname[8])	#'summary_4.txt'
				sys.stderr.write('Chr %s ... \n'%chromosome)
				input_fname = os.path.join(input_dir, input_fname)
				inf = open(input_fname)
				array_id_ls = inf.readline().strip().split()
				array_id_ls = map(int, array_id_ls)
				for line in inf:
					counter  += 1
					line = line.strip()
					part_ls = line.split(':')	#split into 3 parts
					cnvPosPattern_searchResult = cnvPosPattern.search(part_ls[1])
					start = int(cnvPosPattern_searchResult.group(1))
					#probe position is central
					tair9_chr_start = tair8_chr_pos2tair9_chr_pos.get((chromosome, start))
					if tair9_chr_start is None or tair9_chr_start[1] is None:	#skip if it's not found
						no_of_cnvs_skipped += 1
						continue
					start = tair9_chr_start[1]-12
					
					no_of_probes = int(cnvPosPattern_searchResult.group(2))
					stop = int(cnvPosPattern_searchResult.group(3))
					tair9_chr_stop = tair8_chr_pos2tair9_chr_pos.get((chromosome,stop))
					if tair9_chr_stop is None or tair9_chr_stop[1] is None:
						no_of_cnvs_skipped += 1
						continue
					
					stop = tair9_chr_stop[1] + 12
					
					probability_ls = part_ls[2].strip().split()
					probability_ls = map(float, probability_ls)
					
					cnv = Stock_250kDB.CNV(chromosome=chromosome, start=start, \
										stop=stop, no_of_probes_covered=no_of_probes, \
										cnv_method_id=cnv_method_id, cnv_type_id=cnv_type_id)
					session.add(cnv)
					for i in xrange(len(array_id_ls)):
						array_id = array_id_ls[i]
						probability = probability_ls[i]
						if probability>=minProbabilityToCallDeletion:
							cnv_array_call = Stock_250kDB.CNVArrayCall(array_id=array_id, cnv_method_id=cnv_method_id,\
														score=probability)
							cnv_array_call.cnv = cnv
							session.add(cnv_array_call)
							real_counter += 1
							if real_counter%10000==0:
								session.flush()
								session.expunge_all()
					if counter%5000==0:
						sys.stderr.write('%s%s\t%s\t%s'%('\x08'*80, counter, real_counter, no_of_cnvs_skipped))
		sys.stderr.write('%s%s\t%s\t%s\n'%('\x08'*80, counter, real_counter, no_of_cnvs_skipped))
		session.flush()
		session.expunge_all()
	
	"""
		#2010-8-3
		input_dir = os.path.expanduser()
		CNV.putSebastianCallIntoDB(db_250k, input_dir, minProbabilityToCallDeletion=0.8,\
							cnv_method_id=22, cnv_type_id=1)
		sys.exit(0)
		
		#2010-8-3
		input_dir = os.path.expanduser('~/script/variation/data/CNV/SebastianHMMSummaryCalls/')
		CNV.putSebastianCallIntoDB(db_250k, input_dir, minProbabilityToCallDeletion=0.8,\
							cnv_method_id=21, cnv_type_id=1, debug=debug)
		sys.exit(0)
		
	"""
	
	@classmethod
	def putCNVSNPLDIntoDB(cls, db_250k, input_fname, method_id=1):
		"""
		2010-10-10
		"""
		sys.stderr.write("Putting CNV-SNP LD into db ...\n")
		from pymodule.utils import getColName2IndexFromHeader
		import csv
		reader = csv.reader(open(input_fname), delimiter='\t')
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header)
		cnv_id2max_LD = {}
		counter = 0
		real_counter = 0
		for line in reader:
			cnv_id = int(line[col_name2index['snp2']])
			r2_value = float(line[col_name2index['r2']])
			if cnv_id not in cnv_id2max_LD:
				cnv_id2max_LD[cnv_id] = r2_value
				real_counter += 1
			if r2_value>cnv_id2max_LD[cnv_id]:
				cnv_id2max_LD[cnv_id] = r2_value
			counter += 1
			if counter%10000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
		del reader
		sys.stderr.write("%s CNVs with LD.\n"%(len(cnv_id2max_LD)))
		
		
		import Stock_250kDB
		method = Stock_250kDB.GenomeWideResultMethod.get(method_id)
		if method is None:
			method = Stock_250kDB.GenomeWideResultMethod(id=method_id)
			db_250k.session.add(method)
		
		counter = 0
		real_counter = 0
		for cnv_id, max_LD in cnv_id2max_LD.iteritems():
			cnv = Stock_250kDB.CNV.get(cnv_id)
			marker = Stock_250kDB.GenomeMarker.query.filter_by(chromosome=cnv.chromosome).filter_by(start=cnv.start).\
				filter_by(stop=cnv.stop).first()
			if marker is None:
				marker = Stock_250kDB.GenomeMarker(chromosome=cnv.chromosome, start=cnv.start, stop=cnv.stop)
				db_250k.session.add(marker)
			result = Stock_250kDB.GenomeWideResult(value = max_LD)
			result.marker = marker
			result.method = method
			db_250k.session.add(result)
			counter += 1
			if counter%5000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
				db_250k.session.flush()
		sys.stderr.write("Done.\n")
		db_250k.session.flush()
	
	"""
		# 2010-10-10
		input_fname = os.path.expanduser("~/cnvMethod20_vs_callMethod32_LD.tsv")
		CNV.putCNVSNPLDIntoDB(db_250k, input_fname, method_id=1)
		sys.exit(0)
		
		# 2010-10-10
		input_fname = os.path.expanduser("~/cnvMethod22_vs_callMethod32_LD.tsv")
		CNV.putCNVSNPLDIntoDB(db_250k, input_fname, method_id=2)
		sys.exit(0)
	"""
	
	@classmethod
	def mergeTwoCNVQCSegments(cls, segment, next_segment):
		"""
		2010-7-19
			called by getCNVQCDataFromDB()
			
			segment = (row.chromosome, row.start, row.stop, \
					row.size_affected, row.no_of_probes_covered, row.copy_number, row.id, row.ecotype_id)
		"""
		new_segment_stop = max(segment[2], next_segment[2])
		
		# to get the no_of_probes_covered for the merged segment
		if segment[4] is not None and next_segment[4] is not None:
			overlapping_len = float(segment[2]-next_segment[1])	# it could be negative, which would end up increasing the no_of_probes_covered
			import numpy
			# estimate the number of probes in the overlapping region for each segment, take the average
			no_of_probes_in_the_overlapping_of_s1 = int(round(segment[4]*overlapping_len/segment[3]))
			no_of_probes_in_the_overlapping_of_s2 = int(round(next_segment[4]*overlapping_len/next_segment[3]))
			no_of_probes_in_the_overlapping = numpy.mean([no_of_probes_in_the_overlapping_of_s1, no_of_probes_in_the_overlapping_of_s2])
			no_of_probes_covered = segment[4] + next_segment[4] - no_of_probes_in_the_overlapping
		else:
			no_of_probes_covered = None
		
		new_segment = [segment[0], segment[1], new_segment_stop, new_segment_stop-segment[1]+1, no_of_probes_covered] + \
					segment[5:]
		return new_segment
	
	@classmethod
	def getCNVQCDataFromDB(cls, db_250k=None, data_source_id=1, ecotype_id=None, cnv_type_id=None, \
						min_QC_segment_size=None, min_no_of_probes=None, min_reciprocal_overlap=0.6, \
						cmpfn=cmp, mergeOverlappingSegments=True, parseNumOfReadsFromComment=False, \
						cnv_method_id=None, ecotype_id_ls=None):
		"""
		2010-8-4
			fetch the accession-id & ecotype-id info in advance to avoid table-joining involving potentially
				large-data operation (busting the disk space left by creating a huge temporary file in /tmp)
		2010-7-22
			add argument cnv_method_id, ecotype_id_ls to restrict QC
				ecotype_id_ls will be used only if ecotype_id is None.
		2010-7-21
			add argument parseNumOfReadsFromComment.
				If True, number of reads will be parsed out of CNVQCCall.comment and put into segment as score.
		2010-7-19
			add argument mergeOverlappingSegments.
				If True, segments that are overlapping or next to each other will be merged before being put into RBDict.
				
		2010-6-27
			add an overlooked argument db_250k, which was a global variable before.
		2010-2-9
			add argument cmpfn
		2010-1-26
			replace the list structure of cnv_qc_call_data in ecotype_id2cnv_qc_call_data with binary_tree structure
		2009-12-9
			add no_of_probes_covered into returning data
			add cnv_type_id
		2009-11-4
			get CNV QC data from database
		"""
		sys.stderr.write("Getting CNV QC data ... \n")
		import Stock_250kDB
		
		if ecotype_id is not None and not ecotype_id_ls:
			ecotype_id_ls = [ecotype_id]
		
		accession_id_ls = cls.getAccessionIDLsFromCNVQCAccessionGivenEcotypeIDLs(db_250k, ecotype_id_ls, \
									data_source_id=data_source_id)
		accession_id2ecotype_id = cls.getAccessionID2EcotypeIDFromCNVQCAccessionGivenDataSourceID(db_250k, \
									data_source_id=data_source_id)
		
		"""
		sql_string = "select a.ecotype_id, c.chromosome, c.start, c.stop, c.size_affected, c.no_of_probes_covered, \
				c.copy_number, c.id, c.score, c.comment from %s c,\
				%s a where c.accession_id=a.id and a.data_source_id=%s \
				and (c.chromosome=c.stop_chromosome or c.stop_chromosome is null)"%\
				(Stock_250kDB.CNVQCCall.table.name, Stock_250kDB.CNVQCAccession.table.name, data_source_id)
		"""
		str_accession_id_ls = map(str, accession_id_ls)
		sql_string = "select c.accession_id, c.chromosome, c.start, c.stop, c.size_affected, c.no_of_probes_covered, \
				c.copy_number, c.id, c.score, c.comment from %s c where c.accession_id in (%s) \
				and (c.chromosome=c.stop_chromosome or c.stop_chromosome is null)"%\
				(Stock_250kDB.CNVQCCall.table.name, ','.join(str_accession_id_ls))
		if cnv_type_id is not None:
			sql_string += " and c.cnv_type_id=%s"%cnv_type_id
		"""
		if ecotype_id is not None:
			sql_string += " and a.ecotype_id=%s"%ecotype_id
		elif ecotype_id_ls:
			str_ecotype_id_ls = map(str, ecotype_id_ls)
			sql_string += " and a.ecotype_id in (%s)"%','.join(str_ecotype_id_ls)
		"""
		
		if min_no_of_probes is not None:
			sql_string += " and c.no_of_probes_covered>=%s"%min_no_of_probes
		if min_QC_segment_size is not None:
			sql_string += " and (c.stop-c.start+1)>=%s"%min_QC_segment_size
		if cnv_method_id is not None:
			sql_string += " and c.cnv_method_id=%s"%cnv_method_id
		
		if mergeOverlappingSegments:	# 2010-7-19	#assume mapping from accession to ecotype_id is one-to-one mapping for one data source 
			sql_string += " order by accession_id, chromosome, start, stop"
		else:
			sql_string += " order by RAND()"	# 2010-1-26 random ordering to optimize the binary_tree, not needed for RBDict.
		
		rows = db_250k.metadata.bind.execute(sql_string)
		counter = 0
		ecotype_id2cnv_qc_call_data = {}
		#from pymodule.BinarySearchTree import binary_tree	#2010-1-26 replace the list structure of cnv_qc_call_data \
															# in ecotype_id2cnv_qc_call_data with binary_tree structure
		from pymodule.RBTree import RBDict	# 2010-1-26 RBDict is more efficiency than binary_tree.
		from pymodule.CNV import CNVSegmentBinarySearchTreeKey
		
		# 2010-7-19
		ecotype_id2segment_ls = {}
		no_of_original_segments = 0
		
		# 2010-7-20 temporary. extract num_Reads out  of comment
		import re
		num_Reads_re = re.compile(r'num_Reads (\d+),')
		
		for row in rows:	#segments are sorted in chromosomal order
			no_of_original_segments += 1
			ecotype_id = accession_id2ecotype_id.get(row.accession_id)
			if ecotype_id is None:	#2010-8-4 ignore acessions with no ecotype id
				continue
			if ecotype_id not in ecotype_id2segment_ls:
				ecotype_id2segment_ls[ecotype_id] = []
			
			if parseNumOfReadsFromComment:
				num_Reads_search_result = num_Reads_re.search(row.comment)
				if num_Reads_search_result:
					num_Reads = int(num_Reads_search_result.group(1))
				else:
					num_Reads = -1
				score = num_Reads
			else:
				score = row.score
			segment = [row.chromosome, row.start, row.stop, \
					row.size_affected, row.no_of_probes_covered, row.copy_number, row.id, score, ecotype_id]
			
			segment_ls = ecotype_id2segment_ls.get(ecotype_id)
			if len(segment_ls)>0:
				previous_segment = segment_ls[-1]
			else:
				previous_segment = None
			
			if previous_segment is not None and mergeOverlappingSegments and previous_segment[0]==segment[0] and previous_segment[2]>=segment[1]-1:
				previous_segment = cls.mergeTwoCNVQCSegments(previous_segment, segment)
				segment_ls[-1] = previous_segment
			else:
				segment_ls.append(segment)
			ecotype_id2segment_ls[ecotype_id] = segment_ls
			if counter%5000==0:
				sys.stderr.write('%s%s'%('\x08'*80, counter,))
		sys.stderr.write('%s%s\n'%('\x08'*80, counter,))
		
		ecotype_id2data = {}
		from pymodule import PassingData
		counter = 0
		for ecotype_id, segment_ls in ecotype_id2segment_ls.iteritems():
			for segment in segment_ls:
				ecotype_id = segment[-1]
				chromosome, start, stop = segment[:3]
				if ecotype_id not in ecotype_id2cnv_qc_call_data:
					#ecotype_id2cnv_qc_call_data[row.ecotype_id] = binary_tree()
					ecotype_id2cnv_qc_call_data[ecotype_id] = RBDict(cmpfn=cmpfn)
					ecotype_id2data[ecotype_id] = PassingData(length=0, length_affected=0, no_of_segments=0)
				ecotype_id2data[ecotype_id].length += abs(stop-start)+1
				ecotype_id2data[ecotype_id].length_affected += segment[3]
				ecotype_id2data[ecotype_id].no_of_segments += 1
				
				segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=chromosome, span_ls=[start, stop], \
														min_reciprocal_overlap=min_reciprocal_overlap)
				ecotype_id2cnv_qc_call_data[ecotype_id][segmentKey] = segment
				
				#cnv_qc_call_data = ecotype_id2cnv_qc_call_data[row.ecotype_id]
				#cnv_qc_call_data.append((row.chromosome, row.start, row.stop, row.size_affected, row.no_of_probes_covered, row.copy_number, row.id))
				counter += 1
				if counter%5000==0:
					sys.stderr.write('%s%s'%('\x08'*80, counter,))
		sys.stderr.write('%s%s\n'%('\x08'*80, counter,))
		
		import math
		for ecotype_id, tree in ecotype_id2cnv_qc_call_data.iteritems():
			print "\tDepth of Ecotype %s's tree: %d" % (ecotype_id, tree.depth())
			print "\tOptimum Depth: %f (%d) (%f%% depth efficiency)" % (tree.optimumdepth(), math.ceil(tree.optimumdepth()),
													math.ceil(tree.optimumdepth()) / tree.depth())
			#cnv_qc_call_data.sort()
			#ecotype_id2cnv_qc_call_data[ecotype_id] = cnv_qc_call_data
		
		sys.stderr.write("\t%s cnv qc calls out of %s segments for %s ecotypes. Done.\n"%\
						(counter, no_of_original_segments, len(ecotype_id2cnv_qc_call_data)))
		return PassingData(ecotype_id2cnv_qc_call_data=ecotype_id2cnv_qc_call_data, ecotype_id2data=ecotype_id2data)
	
	
	@classmethod
	def getLerContigSpanDataFromDB(cls, data_source_id=1, ecotype_id=None, \
								min_QC_segment_size=None, min_no_of_probes=None, min_reciprocal_overlap=0.6):
		"""
		2010-3-18
			moved the trunk to Stock_250kDB.getSequenceFragmentRefPosFromDBInRBDict() in Stock_250kDB.py
		2010-1-28
			This data provides information regarding which parts of Col genome are actually covered by the Ler contigs,
			which is key to get an accurate FPR for the tiling-array deletions because you'll know for sure this deletion
			is a false positive only when the deletion falls into the Ler contig coverage.
			
		"""
		sys.stderr.write("Getting Ler Contig Span data ... \n")
		return db_250k.getSequenceFragmentRefPosFromDBInRBDict(data_source_id=data_source_id, ecotype_id=ecotype_id, \
								min_QC_segment_size=min_QC_segment_size, min_no_of_probes=min_no_of_probes,\
								min_reciprocal_overlap=min_reciprocal_overlap)
	
	@classmethod
	def checkOverlapOfTwoCNVSegmentRBTree(cls, tree1, tree2):
		"""
		2010-1-28
			check how overlapping the two trees are
		"""
		sys.stderr.write("Checking how overlapping two trees are ...")
		no_of_items_in_tree1 = len(tree1)
		no_of_items_in_tree2 = len(tree2)
		no_of_overlapping = 0
		for nodeKey in tree1:
			if nodeKey in tree2:
				no_of_overlapping += 1
		overlapping_ratio1 = no_of_overlapping/float(no_of_items_in_tree1)
		overlapping_ratio2 = no_of_overlapping/float(no_of_items_in_tree2)
		sys.stderr.write("%s (%s) out of tree1 exist in tree2 (%s).\n"%(no_of_overlapping, overlapping_ratio1, overlapping_ratio2))
	
	@classmethod
	def countMatchedDeletionsFunctor(cls, param_obj, cnv_segment_obj=None, cnv_qc_call=None):
		"""
		used in countNoOfCNVDeletionsMatchQC() to be passed to compareCNVSegmentsAgainstQCHandler()
		
		2009-12-9
			store qc data in param_obj.array_id2qc_data
		2009-11-4
			a functor to be called in 
			
			cnv_segment_obj = PassingData(ecotype_id=cnv_ecotype_id, start_probe=start_probe, stop_probe=stop_probe,\
											no_of_probes=no_of_probes, amplitude=amplitude, segment_length=segment_length,\
											segment_chromosome=segment_chromosome, array_id=array_id,\
											start_probe_id=start_probe_id, stop_probe_id=stop_probe_id,\
											segment_start_pos=segment_start_pos, segment_stop_pos=segment_stop_pos)
			cnvSegmentKey = CNVSegmentBinarySearchTreeKey(chromosome=segment_chromosome, span_ls=[segment_start_pos, segment_stop_pos],\
															min_reciprocal_overlap=min_reciprocal_overlap)
		"""
		from pymodule import PassingData
		if not hasattr(param_obj, 'no_of_valid_deletions'):
			setattr(param_obj, 'no_of_valid_deletions', 0)
		if not hasattr(param_obj, "array_id2qc_data"):
			param_obj.array_id2qc_data = {}
		if not hasattr(param_obj, "array_id2no_of_probes2qc_data"):	# for FPR per no_of_probes
			param_obj.array_id2no_of_probes2qc_data = {}
		if not hasattr(param_obj, "array_id2qc_no_of_probes2qc_data"):	# for FNR per no_of_probes
			param_obj.array_id2qc_no_of_probes2qc_data = {}
		
		array_id = cnv_segment_obj.array_id
		no_of_probes = cnv_segment_obj.no_of_probes
		if array_id not in param_obj.array_id2qc_data:
			param_obj.array_id2qc_data[array_id] = PassingData(ecotype_id=cnv_segment_obj.ecotype_id, \
															no_of_valid_deletions=0,\
															no_of_deletions=0,\
															cnv_qc_call_id_set=set())
			param_obj.array_id2no_of_probes2qc_data[array_id] = {}
			param_obj.array_id2qc_no_of_probes2qc_data[array_id] = {}
		if no_of_probes not in param_obj.array_id2no_of_probes2qc_data[array_id]:
			param_obj.array_id2no_of_probes2qc_data[array_id][no_of_probes] = PassingData(ecotype_id=cnv_segment_obj.ecotype_id, \
															no_of_valid_deletions=0,\
															no_of_deletions=0,\
															cnv_qc_call_id_set=set())
		# 2010-1-26 increase the no_of_deletions counter no matter whether there's a corresponding cnv_qc_call or not
		param_obj.array_id2qc_data[array_id].no_of_deletions += 1
		param_obj.array_id2no_of_probes2qc_data[array_id][no_of_probes].no_of_deletions += 1
		if cnv_qc_call is not None:
			qc_chromosome, qc_start, qc_stop = cnv_qc_call[:3]
			cnv_qc_call_id = cnv_qc_call[-1]
			param_obj.array_id2qc_data[array_id].cnv_qc_call_id_set.add(cnv_qc_call_id)
			param_obj.array_id2qc_data[array_id].no_of_valid_deletions += 1
			
			qc_no_of_probes = cnv_qc_call[4]
			if qc_no_of_probes not in param_obj.array_id2qc_no_of_probes2qc_data[array_id]:
				param_obj.array_id2qc_no_of_probes2qc_data[array_id][qc_no_of_probes] = PassingData(ecotype_id=cnv_segment_obj.ecotype_id, \
															cnv_qc_call_id_set=set())
				param_obj.array_id2qc_no_of_probes2qc_data[array_id][qc_no_of_probes].cnv_qc_call_id_set.add(cnv_qc_call_id)
			
			param_obj.array_id2no_of_probes2qc_data[array_id][no_of_probes].cnv_qc_call_id_set.add(cnv_qc_call_id)
			param_obj.array_id2no_of_probes2qc_data[array_id][no_of_probes].no_of_valid_deletions += 1

	@classmethod
	def countDeletionsMatchedByNormalCoverage(cls, param_obj, cnv_segment_obj=None, **keywords):
		"""
		2010-6-17
			similar to countMatchedDeletionsFunctor, but 
			
			used in countNoOfCNVDeletionsMatchQC() to be passed to compareCNVSegmentsAgainstQCHandler()
			
			
			cnv_segment_obj = PassingData(ecotype_id=cnv_ecotype_id, start_probe=start_probe, stop_probe=stop_probe,\
											no_of_probes=no_of_probes, amplitude=amplitude, segment_length=segment_length,\
											segment_chromosome=segment_chromosome, array_id=array_id,\
											start_probe_id=start_probe_id, stop_probe_id=stop_probe_id,\
											segment_start_pos=segment_start_pos, segment_stop_pos=segment_stop_pos)
			cnvSegmentKey = CNVSegmentBinarySearchTreeKey(chromosome=segment_chromosome, span_ls=[segment_start_pos, segment_stop_pos],\
															min_reciprocal_overlap=min_reciprocal_overlap)
		"""
		isDeletionTrue = param_obj.isDeletionTrue
		
		from pymodule import PassingData
		if not hasattr(param_obj, 'no_of_valid_deletions'):
			setattr(param_obj, 'no_of_valid_deletions', 0)
		if not hasattr(param_obj, "array_id2qc_data"):
			param_obj.array_id2qc_data = {}
		if not hasattr(param_obj, "array_id2no_of_probes2qc_data"):	# for FPR per no_of_probes
			param_obj.array_id2no_of_probes2qc_data = {}
		if not hasattr(param_obj, "array_id2qc_no_of_probes2qc_data"):	# for FNR per no_of_probes
			param_obj.array_id2qc_no_of_probes2qc_data = {}
		
		array_id = cnv_segment_obj.array_id
		no_of_probes = int(cnv_segment_obj.no_of_probes/float(cnv_segment_obj.segment_length/1000.))	#2010-6-17 number of probes per kb
		if array_id not in param_obj.array_id2qc_data:
			param_obj.array_id2qc_data[array_id] = PassingData(ecotype_id=cnv_segment_obj.ecotype_id, \
															no_of_valid_deletions=0,\
															no_of_deletions=0,\
															cnv_qc_call_id_set=set())
			param_obj.array_id2no_of_probes2qc_data[array_id] = {}
			param_obj.array_id2qc_no_of_probes2qc_data[array_id] = {}
		if no_of_probes not in param_obj.array_id2no_of_probes2qc_data[array_id]:
			param_obj.array_id2no_of_probes2qc_data[array_id][no_of_probes] = PassingData(ecotype_id=cnv_segment_obj.ecotype_id, \
															no_of_valid_deletions=0,\
															no_of_deletions=0,\
															cnv_qc_call_id_set=set())
		# 2010-1-26 increase the no_of_deletions counter no matter whether there's a corresponding cnv_qc_call or not
		param_obj.array_id2qc_data[array_id].no_of_deletions += 1
		param_obj.array_id2no_of_probes2qc_data[array_id][no_of_probes].no_of_deletions += 1
		if isDeletionTrue:
			cnv_qc_call_id = len(param_obj.array_id2qc_data[array_id].cnv_qc_call_id_set)
			param_obj.array_id2qc_data[array_id].cnv_qc_call_id_set.add(cnv_qc_call_id)
			param_obj.array_id2qc_data[array_id].no_of_valid_deletions += 1
			
			param_obj.array_id2no_of_probes2qc_data[array_id][no_of_probes].cnv_qc_call_id_set.add(cnv_qc_call_id)
			param_obj.array_id2no_of_probes2qc_data[array_id][no_of_probes].no_of_valid_deletions += 1
			
	@classmethod
	def aggregateEnoughDataPoints(cls, no_of_data_points_ls=[], lists_to_be_summed=[], \
								lists_to_be_sampled=[], lists_to_be_concatenated=[],\
								min_no_of_data_points=20):
		"""
		2010-2-9
			Given a range of indices (aggregating_start_index:aggregating_stop_index) which amass enough data points,
			each list in following "lists" will undergo different operations:
				lists_to_be_summed: new-element = sum(ls[aggregating_start_index:aggregating_stop_index])
				lists_to_be_sampled: new-element = ls[aggregating_start_index]
				lists_to_be_concatenated: new-ls=[];for k in range(aggregating_start_index, aggregating_stop_index):;new-ls.extend(ls[k]);
											new-element = new-ls
		2010-2-9
			add argument lists_to_be_concatenated
		2010-1-28
			it's a function to ensure a statistic having enough data points to be meaningful.
			no_of_data_points_ls is used to check from which to which needs aggregation (summarization).
		"""
		n = len(no_of_data_points_ls)
		if n==0:	# 2010-2-10
			return None
		aggregating_index_ls = [0]	# If aggregating_index_ls = [0,2,4, ...], no_of_data_points_ls = [sum(no_of_data_points_ls[0:2]), sum(no_of_data_points_ls[2:4], ...)
		for i in range(1, n+1):
			last_index = aggregating_index_ls[-1]
			no_of_data_points = sum(no_of_data_points_ls[last_index:i])	# from 
			if no_of_data_points>=min_no_of_data_points:	# there's possibility that the last batch never makes to the min_no_of_data_points and gets tossed.
				aggregating_index_ls.append(i)
		
		if len(aggregating_index_ls)==1 and n>0:	#2010-2-9 it's either no_of_data_points_ls has only one item in it or sum(no_of_data_points_ls)<min_no_of_data_points
			aggregating_index_ls.append(n)
			
		new_no_of_data_points_ls = []
		
		no_of_indices = len(aggregating_index_ls)
		new_lists_to_be_summed = []
		new_lists_to_be_sampled = []
		new_lists_to_be_concatenated = []
		for i in range(no_of_indices-1):
			aggregating_start_index = aggregating_index_ls[i]
			aggregating_stop_index = aggregating_index_ls[i+1]
			new_no_of_data_points_ls.append(sum(no_of_data_points_ls[aggregating_start_index:aggregating_stop_index]))
			for j in range(len(lists_to_be_summed)):
				ls = lists_to_be_summed[j]
				if i==0:	# first time to add 
					new_lists_to_be_summed.append([])
				new_lists_to_be_summed[j].append(sum(ls[aggregating_start_index:aggregating_stop_index]))	# sum the range
			for j in range(len(lists_to_be_sampled)):
				ls = lists_to_be_sampled[j]
				if i==0:	# first time to add 
					new_lists_to_be_sampled.append([])
				new_lists_to_be_sampled[j].append(ls[aggregating_start_index])	# take the first one in the list
			for j in range(len(lists_to_be_concatenated)):
				ls = lists_to_be_concatenated[j]
				if i==0:	# first time to add 
					new_lists_to_be_concatenated.append([])
				new_ls = []
				for k in range(aggregating_start_index, aggregating_stop_index):
					new_ls.extend(ls[k])	# extend the list
				new_lists_to_be_concatenated[j].append(new_ls)
		from pymodule import PassingData
		return_data = PassingData(no_of_data_points_ls=new_no_of_data_points_ls, lists_to_be_summed=new_lists_to_be_summed, \
					lists_to_be_sampled=new_lists_to_be_sampled, lists_to_be_concatenated=new_lists_to_be_concatenated)
		return return_data
	
	@classmethod
	def outputFalseNegativeRate(cls, param_obj):
		"""
		2009-12-9
			calculate FNR for each class with same number of probes
		2009-11-4
		"""
		for array_id, qc_data in param_obj.array_id2qc_data.iteritems():
			no_of_QCCalls_matched = len(qc_data.cnv_qc_call_id_set)
			no_of_total_QCCalls = len(param_obj.ecotype_id2cnv_qc_call_data[qc_data.ecotype_id])
			false_negative_rate = (no_of_total_QCCalls-no_of_QCCalls_matched)/float(no_of_total_QCCalls)
			sys.stderr.write("Array %s false negative rate: %s/%s(%s).\n"%(array_id, \
																		no_of_total_QCCalls-no_of_QCCalls_matched,\
																		no_of_total_QCCalls, false_negative_rate))
			
			if getattr(param_obj, 'array_id2qc_no_of_probes2qc_data', None):
				qc_no_of_probes2qc_data = param_obj.array_id2qc_no_of_probes2qc_data[array_id]
				no_of_probes_ls = qc_no_of_probes2qc_data.keys()
				no_of_probes_ls.sort()
				no_of_total_QCCalls_ls = []
				no_of_FalseNegatives_ls = []
				for no_of_probes in no_of_probes_ls:
					qc_data = qc_no_of_probes2qc_data[no_of_probes]
					no_of_QCCalls_matched = len(qc_data.cnv_qc_call_id_set)
					no_of_total_QCCalls = len(param_obj.ecotype_id2qc_no_of_probes2cnv_qc_call_id_set[qc_data.ecotype_id][no_of_probes])
					false_negative_rate = (no_of_total_QCCalls-no_of_QCCalls_matched)/float(no_of_total_QCCalls)
					no_of_total_QCCalls_ls.append(no_of_total_QCCalls)
					no_of_FalseNegatives_ls.append(no_of_total_QCCalls-no_of_QCCalls_matched)
					sys.stderr.write("\t%s\t%s\t%s\t%s\n"%(no_of_probes, \
														no_of_total_QCCalls-no_of_QCCalls_matched,\
														no_of_total_QCCalls, false_negative_rate))
	
				if getattr(param_obj, 'output_fname_prefix', None):
					output_fname = '%s_array_%s_FNR.png'%(param_obj.output_fname_prefix, array_id)
					aggregateReturnData = cls.aggregateEnoughDataPoints(no_of_total_QCCalls_ls, \
																	[no_of_FalseNegatives_ls], 
																	[no_of_probes_ls])
					new_no_of_total_QCCalls_ls = aggregateReturnData.no_of_data_points_ls
					new_no_of_FalseNegatives_ls = aggregateReturnData.lists_to_be_summed[0]
					new_no_of_probes_ls = aggregateReturnData.lists_to_be_sampled[0]
					divide_functor = lambda x: x[1]/float(x[0])
					FNR_ls = map(divide_functor, zip(new_no_of_total_QCCalls_ls, new_no_of_FalseNegatives_ls))
					import pylab
					pylab.clf()
					pylab.plot(new_no_of_probes_ls, FNR_ls, '.')
					pylab.title('array %s'%array_id)
					pylab.xlabel('No of Probes')
					pylab.ylabel('FNR')
					pylab.savefig(output_fname, dpi=300)
	@classmethod
	def outputFalsePositiveRate(cls, param_obj):
		"""
		2009-12-9
			calculate FNR for each class with same number of probes
		2009-11-4
		"""
		for array_id, qc_data in param_obj.array_id2qc_data.iteritems():
			no_of_valid_deletions = qc_data.no_of_valid_deletions
			no_of_deletions = qc_data.no_of_deletions
			no_of_non_valid_deletions = no_of_deletions-no_of_valid_deletions
			FDR = no_of_non_valid_deletions/float(no_of_deletions)
			sys.stderr.write("Array %s false discovery rate: %s/%s(%s).\n"%(array_id, \
							no_of_non_valid_deletions, no_of_deletions, FDR))
			if getattr(param_obj, 'array_id2no_of_probes2qc_data', None):
				no_of_probes2qc_data = param_obj.array_id2no_of_probes2qc_data[array_id]
				output_fname_prefix = getattr(param_obj, 'output_fname_prefix', None)
				if output_fname_prefix:
					output_fname_prefix = '%s_array_%s'%(output_fname_prefix, array_id)
					
					array_id2label = getattr(param_obj, 'array_id2label', None)
					if array_id2label:
						title = array_id2label.get(array_id)
						title += ' FDR %.3f %s objs.'%(FDR, sum(no_of_deletions_ls))
					else:
						title = None
					if not title:
						title = 'array %s'%array_id
						
					cls.plotOneFDR(no_of_probes2qc_data, output_fname_prefix, min_no_of_data_points=getattr(param_obj, 'min_no_of_data_points', 30), \
								fig_title=title, \
								xlabel=getattr(param_obj, 'xlabel', 'No of Probes'), ylabel = getattr(param_obj, 'ylabel', 'FDR'))
	
	@classmethod
	def outputCNVSegmentObj(cls, cnv_segment_obj, param_obj):
		"""
		2010-2-10
			output cnv_segment_obj from CNV.compareCNVSegmentsAgainstQCHandler() in the RunGADA.py output format.
		"""
		cnv_segment_output_fname_prefix = getattr(param_obj, 'cnv_segment_output_fname_prefix', '')
		
		if not cnv_segment_output_fname_prefix:	# file name is not, not output 
			return
		
		if not hasattr(param_obj, 'cnv_segment_output_dict'):
			param_obj.cnv_segment_output_dict = {}
		if cnv_segment_obj.no_of_probes not in param_obj.cnv_segment_output_dict:
			from RunGADA import RunGADA
			import csv
			output_fname = '%s_%s_probes.tsv'%(cnv_segment_output_fname_prefix, cnv_segment_obj.no_of_probes)
			cnv_segment_output_writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
			RunGADA.output_header(cnv_segment_output_writer)
			param_obj.cnv_segment_output_dict[cnv_segment_obj.no_of_probes] = cnv_segment_output_writer
		cnv_segment_output_writer = param_obj.cnv_segment_output_dict.get(cnv_segment_obj.no_of_probes)
		start_probe = map(str, cnv_segment_obj.start_probe)
		start_probe = '_'.join(start_probe)
		stop_probe = map(str, cnv_segment_obj.stop_probe)
		stop_probe = '_'.join(stop_probe)
		new_row = [cnv_segment_obj.ecotype_id, cnv_segment_obj.array_id, start_probe, stop_probe, cnv_segment_obj.segment_length, \
				cnv_segment_obj.amplitude, start_probe, stop_probe]
		cnv_segment_output_writer.writerow(new_row)
	
	@classmethod
	def isDeletionTrueBasedOnOtherDeletionData(cls, cnv_segment_obj, param_obj):
		"""
		2010-6-17
			be used in countNoOfCNVDeletionsMatchQC() to pass it to compareCNVSegmentsAgainstQCHandler()
		
			cnv_segment_obj = PassingData(ecotype_id=cnv_ecotype_id, start_probe=start_probe, stop_probe=stop_probe,\
											no_of_probes=no_of_probes, amplitude=amplitude, segment_length=segment_length,\
											segment_chromosome=segment_chromosome, array_id=array_id,\
											start_probe_id=start_probe_id, stop_probe_id=stop_probe_id,\
											segment_start_pos=segment_start_pos, segment_stop_pos=segment_stop_pos)
			cnvSegmentKey = CNVSegmentBinarySearchTreeKey(chromosome=segment_chromosome, span_ls=[segment_start_pos, segment_stop_pos],\
															min_reciprocal_overlap=min_reciprocal_overlap)
		"""
		from pymodule import PassingData
		ecotype_id = cnv_segment_obj.ecotype_id
		cnvSegmentKey = param_obj.cnvSegmentKey
		ecotype_id2cnv_qc_call_data = getattr(param_obj, "ecotype_id2cnv_qc_call_data", None)
		
		if ecotype_id2cnv_qc_call_data is not None and ecotype_id in ecotype_id2cnv_qc_call_data:
			cnv_qc_call_data = ecotype_id2cnv_qc_call_data.get(ecotype_id)
			if cnv_qc_call_data:
				cnv_qc_call = cnv_qc_call_data.get(cnvSegmentKey)
			else:
				cnv_qc_call = None
			if cnv_qc_call:
				return PassingData(isDeletionTrue=True, cnv_qc_call=cnv_qc_call)
			else:
				return PassingData(isDeletionTrue=False, cnv_qc_call=None)
		elif ecotype_id2cnv_qc_call_data is None:
			return PassingData(isDeletionTrue=False, cnv_qc_call=None)
	
	@classmethod
	def isDeletionTrueBasedOnNormalCoverage(cls, cnv_segment_obj=None, param_obj=None, ecotype_id=None):
		"""
		2010-6-27
			add argument ecotype_id to be invoked in updateCNVCallPercUnCoveredByLerContig()
		2010-6-17
			be used in countNoOfCNVDeletionsMatchQC() to pass it to compareCNVSegmentsAgainstQCHandler()
		
			cnv_segment_obj = PassingData(ecotype_id=cnv_ecotype_id, start_probe=start_probe, stop_probe=stop_probe,\
											no_of_probes=no_of_probes, amplitude=amplitude, segment_length=segment_length,\
											segment_chromosome=segment_chromosome, array_id=array_id,\
											start_probe_id=start_probe_id, stop_probe_id=stop_probe_id,\
											segment_start_pos=segment_start_pos, segment_stop_pos=segment_stop_pos)
			cnvSegmentKey = CNVSegmentBinarySearchTreeKey(chromosome=segment_chromosome, span_ls=[segment_start_pos, segment_stop_pos],\
															min_reciprocal_overlap=min_reciprocal_overlap)
		"""
		#cnvSegmentKey = CNVSegmentBinarySearchTreeKey(chromosome=cnv_segment_obj.segment_chromosome, \
		#									span_ls=[cnv_segment_obj.segment_start_pos, cnv_segment_obj.segment_stop_pos],\
		#									min_reciprocal_overlap=param_obj.min_reciprocal_overlap)
		from pymodule import PassingData
		cnvSegmentKey = param_obj.cnvSegmentKey
		ecotype_id2normal_data = param_obj.ecotype_id2normal_data
		maxNormalPercInDeletion = param_obj.maxNormalPercInDeletion
		
		if ecotype_id is not None:
			ecotype_id = ecotype_id
		elif cnv_segment_obj is not None:
			ecotype_id = getattr(cnv_segment_obj, 'ecotype_id', None)
		
		if ecotype_id is None:
			return None
		normal_data_RBDict = ecotype_id2normal_data.get(ecotype_id)
		node_ls = []
		normal_data_RBDict.findNodes(cnvSegmentKey, node_ls)
		
		from pymodule.CNV import get_overlap_ratio
		normalPerc = 0.0
		overlap_ratio_ls = []
		for node in node_ls:
			key = node.key
			overlap1, overlap2 = get_overlap_ratio(cnvSegmentKey.span_ls, key.span_ls)[:2]
			normalPerc += overlap1
			overlap_ratio_ls.append((overlap1, overlap2))
		
		if normalPerc<=maxNormalPercInDeletion:
			return PassingData(isDeletionTrue=True, normalPerc=normalPerc, overlap_ratio_ls=overlap_ratio_ls)
		else:
			return PassingData(isDeletionTrue=False, normalPerc=normalPerc, overlap_ratio_ls=overlap_ratio_ls)
		
	
	@classmethod
	def isDeletionTrueBasedOnDeletedFraction(cls, cnv_segment_obj=None, param_obj=None, ecotype_id=None):
		"""
		2010-7-20
			opposite of isDeletionTrueBasedOnNormalCoverage()
		
			cnv_segment_obj = PassingData(ecotype_id=cnv_ecotype_id, start_probe=start_probe, stop_probe=stop_probe,\
											no_of_probes=no_of_probes, amplitude=amplitude, segment_length=segment_length,\
											segment_chromosome=segment_chromosome, array_id=array_id,\
											start_probe_id=start_probe_id, stop_probe_id=stop_probe_id,\
											segment_start_pos=segment_start_pos, segment_stop_pos=segment_stop_pos)
			cnvSegmentKey = CNVSegmentBinarySearchTreeKey(chromosome=segment_chromosome, span_ls=[segment_start_pos, segment_stop_pos],\
															min_reciprocal_overlap=min_reciprocal_overlap)
		"""
		#cnvSegmentKey = CNVSegmentBinarySearchTreeKey(chromosome=cnv_segment_obj.segment_chromosome, \
		#									span_ls=[cnv_segment_obj.segment_start_pos, cnv_segment_obj.segment_stop_pos],\
		#									min_reciprocal_overlap=param_obj.min_reciprocal_overlap)
		from pymodule import PassingData
		cnvSegmentKey = param_obj.cnvSegmentKey
		ecotype_id2qc_data = param_obj.ecotype_id2qc_data
		minDeletedFraction = param_obj.minDeletedFraction
		
		if ecotype_id is not None:
			ecotype_id = ecotype_id
		elif cnv_segment_obj is not None:
			ecotype_id = getattr(cnv_segment_obj, 'ecotype_id', None)
		
		if ecotype_id is None:
			return None
		qc_data_RBDict = ecotype_id2qc_data.get(ecotype_id)
		node_ls = []
		qc_data_RBDict.findNodes(cnvSegmentKey, node_ls)
		
		from pymodule.CNV import get_overlap_ratio
		deletedFraction = 0.0
		overlap_ratio_ls = []
		for node in node_ls:
			key = node.key
			overlap1, overlap2 = get_overlap_ratio(cnvSegmentKey.span_ls, key.span_ls)[:2]
			deletedFraction += overlap1
			overlap_ratio_ls.append((overlap1, overlap2))
		
		if deletedFraction>=minDeletedFraction:
			return PassingData(isDeletionTrue=True, deletedFraction=deletedFraction, overlap_ratio_ls=overlap_ratio_ls)
		else:
			return PassingData(isDeletionTrue=False, deletedFraction=deletedFraction, overlap_ratio_ls=overlap_ratio_ls)
		
	
	@classmethod
	def countNoOfCNVDeletionsMatchQC(cls, db_250k, input_fname_ls, output_fname_prefix=None, ecotype_id=6909, \
								data_source_id=3, cnv_type_id=1, \
								min_QC_segment_size=200, deletion_cutoff=-0.33, max_boundary_diff=10000, \
								max_diff_perc=0.10, min_no_of_probes=5,\
								count_embedded_segment_as_match=True, min_reciprocal_overlap=0.0000001,\
								maxNormalPercInDeletion=0.4,\
								outputDeletionWithinSpan=False, judgeTrueDeletionFunctor=None, afterMatchFunction=None):
		"""
		2010-6-17
			add argument maxNormalPercInDeletion
				judgeTrueDeletionFunctor
					options:
						isDeletionTrueBasedOnNormalCoverage()
						isDeletionTrueBasedOnOtherDeletionData()
					default is cls.isDeletionTrueBasedOnNormalCoverage
				afterMatchFunction:
					options:
						countDeletionsMatchedByNormalCoverage()
						countMatchedDeletionsFunctor()
			
			Argument min_reciprocal_overlap is no longer used to check whether one CNV deletion matches the QC deletion.
				Now it's used to check if any normal data from ecotype_id2normal_data has any overlapping with a CNV deletion.
				So it's just barely above zero.
		2010-2-10
			add argument outputDeletionWithinSpan, output_fname_prefix
		2010-1-26
			pass min_reciprocal_overlap to cls.getCNVQCDataFromDB()
			count_embedded_segment_as_match is meaningless from this on.
		2009-12-9
			calculate FNR for each class with same number of probes
		2009-10-29
			for all CNV deletions, check how many are in QC dataset.
		"""
		from pymodule.CNV import CNVSegmentBinarySearchTreeKey, leftWithinRightAlsoEqualCmp, rightWithinLeftAlsoEqualCmp
		
		#ecotype_id2cnv_qc_call_data = cls.getCNVQCDataFromDB(data_source_id, ecotype_id, cnv_type_id=cnv_type_id, \
		#							min_QC_segment_size=min_QC_segment_size, \
		#							min_no_of_probes=min_no_of_probes,\
		#							min_reciprocal_overlap=min_reciprocal_overlap)
		ecotype_id2cnv_qc_call_data = None
		
		"""
		2010-4
		#ecotype_id2span_data = cls.getLerContigSpanDataFromDB(db_250k, data_source_id=11, ecotype_id=ecotype_id, \
		#													min_QC_segment_size=min_QC_segment_size, min_no_of_probes=min_no_of_probes, \
		#													min_reciprocal_overlap=min_reciprocal_overlap)	# 2010-1-28 special data_source_id 11
		"""
		
		"""
		2010-6-8 new span data from CNV.LerContig.discoverLerDeletionDuplication()
		"""
		#ecotype_id2span_data = cls.getCNVQCDataFromDB(data_source_id=11, ecotype_id=ecotype_id, cnv_type_id=6, \
		#											min_QC_segment_size=None, min_no_of_probes=5,\
		#											min_reciprocal_overlap=0.8, cmpfn=leftWithinRightAlsoEqualCmp)
						# 2010-6-8 special data_source_id 11, set min_reciprocal_overlap to 0.8 ().
						# leftWithinRightAlsoEqualCmp is key to judge whether a segment is in a span or not
		ecotype_id2span_data=None
		
		"""2010-6-17
			This data_source_id 15 is nucmer, maxmatch coords results. Overlapping fragments have all been merged.
			min_reciprocal_overlap=0.0001 because if there's any overlapping between CNV segment and stored nodes,
				that node should be returned.
			ecotype_id2normal_data is used to disqaulify a CNV deletion if it has too many normal data within.
		"""
		cnvQCData = cls.getCNVQCDataFromDB(db_250k, data_source_id=15, ecotype_id=ecotype_id, cnv_type_id=None, \
													min_QC_segment_size=None, min_no_of_probes=None,\
													min_reciprocal_overlap=min_reciprocal_overlap, \
													cmpfn=rightWithinLeftAlsoEqualCmp)
		ecotype_id2normal_data = cnvQCData.ecotype_id2cnv_qc_call_data
		
		if ecotype_id2cnv_qc_call_data and ecotype_id2span_data:
			for ecotype_id, cnv_qc_call_data in ecotype_id2cnv_qc_call_data.iteritems():
				if ecotype_id in ecotype_id2span_data:
					span_data = ecotype_id2span_data.get(ecotype_id)
					sys.stderr.write("Ecotype %s\n"%ecotype_id)
					cls.checkOverlapOfTwoCNVSegmentRBTree(cnv_qc_call_data, span_data)
				
		from pymodule import PassingData
		output_fname_tag = 'ecotype_%s_data_source_%s_deletion_cutoff_%s_min_no_of_probes_%s_min_reciprocal_overlap_%s_maxNormalPercInDeletion%s'%\
			(ecotype_id, data_source_id, deletion_cutoff, min_no_of_probes, min_reciprocal_overlap, maxNormalPercInDeletion)
		param_obj = PassingData(no_of_valid_deletions=0, cnv_qc_call_id_set=set(), array_id2qc_data={}, \
							ecotype_id2span_data = ecotype_id2span_data, ecotype_id2normal_data=ecotype_id2normal_data, \
							cnv_segment_output_fname_prefix = output_fname_prefix+'deletions_within_span_%s'%output_fname_tag,\
							maxNormalPercInDeletion=maxNormalPercInDeletion, min_reciprocal_overlap=min_reciprocal_overlap,\
							minDeletedFraction=1-maxNormalPercInDeletion)
		if outputDeletionWithinSpan:
			functor_after_span = cls.outputCNVSegmentObj
		else:
			functor_after_span = None
		if judgeTrueDeletionFunctor is None:
			judgeTrueDeletionFunctor=cls.isDeletionTrueBasedOnNormalCoverage
		
		cls.compareCNVSegmentsAgainstQCHandler(input_fname_ls, ecotype_id2cnv_qc_call_data=ecotype_id2normal_data, \
						param_obj=param_obj, \
						deletion_cutoff=deletion_cutoff, max_boundary_diff=max_boundary_diff, \
						max_diff_perc=max_diff_perc, min_no_of_probes=min_no_of_probes, \
						afterMatchFunction=afterMatchFunction,\
						count_embedded_segment_as_match=count_embedded_segment_as_match, \
						min_reciprocal_overlap=min_reciprocal_overlap, report=False, \
						functor_after_span=functor_after_span,\
						judgeTrueDeletionFunctor=judgeTrueDeletionFunctor)
		sys.stderr.write("For ecotype_id %s, data_source_id %s, min_QC_segment_size %s, deletion_cutoff: %s, \n\
						min_no_of_probes: %s, min_reciprocal_overlap: %s, maxNormalPercInDeletion %s.\n"%\
						(ecotype_id, data_source_id, min_QC_segment_size, deletion_cutoff, \
						min_no_of_probes, min_reciprocal_overlap, maxNormalPercInDeletion))
		param_obj.ecotype_id2cnv_qc_call_data = ecotype_id2cnv_qc_call_data
		
		"""
		param_obj.ecotype_id2qc_no_of_probes2cnv_qc_call_id_set = {}
		for ecotype_id, cnv_qc_call_data in ecotype_id2cnv_qc_call_data.iteritems():
			qc_no_of_probes2cnv_qc_call_id_set = {}
			for cnv_qc_call in  cnv_qc_call_data:
				no_of_probes_covered = cnv_qc_call[4]
				cnv_qc_call_id = cnv_qc_call[-1]
				if no_of_probes_covered not in qc_no_of_probes2cnv_qc_call_id_set:
					qc_no_of_probes2cnv_qc_call_id_set[no_of_probes_covered] = set()
				qc_no_of_probes2cnv_qc_call_id_set[no_of_probes_covered].add(cnv_qc_call_id)
			param_obj.ecotype_id2qc_no_of_probes2cnv_qc_call_id_set[ecotype_id] = qc_no_of_probes2cnv_qc_call_id_set
		"""
		
		param_obj.output_fname_prefix = output_fname_prefix+output_fname_tag
		cls.outputFalsePositiveRate(param_obj)
		#cls.outputFalseNegativeRate(param_obj)
		
		
	"""
	# 2009-11-4
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s.GADA_A0.5T4M5.tsv'%i))
	CNV.countNoOfCNVDeletionsMatchQC(db_250k, input_fname_ls, ecotype_id=8215, data_source_id=3, cnv_type_id=1, \
								min_QC_segment_size=200, deletion_cutoff=-0.33, max_boundary_diff=10000, max_diff_perc=0.10)
	
	# 2009-11-5
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s.GADA_A0.5T4M5.tsv'%i))
	#input_fname_ls = [os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr1.GADA_A0.5T4M5.tsv')]
	ecotype_id_data_source_id_ls = [(6911, 5), (8215, 3), (9169, 3), (6962, 3)]	# Cvi from Bob, Fei-0 from Schneeberger
	min_QC_segment_size = 5
	min_no_of_probes = 5
	count_embedded_segment_as_match = True
	for ecotype_id, data_source_id in ecotype_id_data_source_id_ls:
		for deletion_cutoff in [-0.33, -0.5]:
			for min_reciprocal_overlap in [0.4, 0.6, 0.8]:
				CNV.countNoOfCNVDeletionsMatchQC(db_250k, input_fname_ls, ecotype_id=ecotype_id, data_source_id=data_source_id, \
										cnv_type_id=1,\
										min_QC_segment_size=min_QC_segment_size, deletion_cutoff=deletion_cutoff, \
										min_no_of_probes=min_no_of_probes,\
										count_embedded_segment_as_match=count_embedded_segment_as_match,\
										min_reciprocal_overlap=min_reciprocal_overlap)
	
	
	# 2009-12-09 use min_reciprocal_overlap
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s.GADA_A0.5T4M5.tsv'%i))
	#input_fname_ls = [os.path.expanduser('~/mnt/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr1.GADA_A0.5T4M5.tsv')]
	#input_fname_ls = []
	#for i in range(1,6):
	#	input_fname_ls.append(os.path.expanduser('~/mnt2/panfs//250k/CNV/GADA_output/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s_A2.0T8.0M5.tsv'%i))
	
	ecotype_id_data_source_id_ls = [(6932, 7)]	# two different types of Ler-1 deletions from Ler contigs
	min_QC_segment_size = 5
	min_no_of_probes = 5
	count_embedded_segment_as_match = True
	for ecotype_id, data_source_id in ecotype_id_data_source_id_ls:
		for deletion_cutoff in [-0.33, -0.5]:
			for min_reciprocal_overlap in [0.4, 0.8]:
				CNV.countNoOfCNVDeletionsMatchQC(db_250k, input_fname_ls, \
					output_fname_prefix=os.path.expanduser('~/tmp/LerContig/call_method_48_CNV_intensity_QNorm_sub_ref.GADA_A0.5T4M5'),\
					ecotype_id=ecotype_id, data_source_id=data_source_id, \
										cnv_type_id=1,\
										min_QC_segment_size=min_QC_segment_size, deletion_cutoff=deletion_cutoff, \
										min_no_of_probes=min_no_of_probes,\
										count_embedded_segment_as_match=count_embedded_segment_as_match,\
										min_reciprocal_overlap=min_reciprocal_overlap, outputDeletionWithinSpan=True)
	
	# 2010-3-4 check histogram separation after wave correction
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt/panfs/250k/CNV/ColLerCviVanFeiSha/call_48_CNV_QNormalize_chr%s_intensity_2.7_3.1_span_2000_probes_wave_corrected.GADA_A0.5T4M5.tsv'%i))
	ecotype_id_data_source_id_ls = [(6932, 7)]	# two different types of Ler-1 deletions from Ler contigs
	min_QC_segment_size = 5
	min_no_of_probes = 5
	count_embedded_segment_as_match = True
	for ecotype_id, data_source_id in ecotype_id_data_source_id_ls:
		for deletion_cutoff in [-0.33, -0.5]:
			for min_reciprocal_overlap in [0.4, 0.8]:
				CNV.countNoOfCNVDeletionsMatchQC(db_250k, input_fname_ls, \
					output_fname_prefix=os.path.expanduser('~/tmp/LerContig/call_48_CNV_QNorm_intensity_wave_corrected_2.7_3.1_span_2000_probes.GADA_A0.5T4M5'),\
					ecotype_id=ecotype_id, data_source_id=data_source_id, \
										cnv_type_id=1,\
										min_QC_segment_size=min_QC_segment_size, deletion_cutoff=deletion_cutoff, \
										min_no_of_probes=min_no_of_probes,\
										count_embedded_segment_as_match=count_embedded_segment_as_match,\
										min_reciprocal_overlap=min_reciprocal_overlap, outputDeletionWithinSpan=True)
	
		#input_fname_ls = []
		#for i in range(1,6):
		#	input_fname_ls.append(os.path.expanduser('~/mnt/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s.GADA_A0.5T4M5.tsv'%i))
		#input_fname_ls = [os.path.expanduser('~/mnt/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr1.GADA_A0.5T4M5.tsv')]
		#input_fname_ls = []
		#for i in range(1,6):
		#	input_fname_ls.append(os.path.expanduser('~/mnt2/panfs//250k/CNV/GADA_output/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s_A2.0T8.0M5.tsv'%i))
		# 2010-6-25 use CNV.isDeletionTrueBasedOnNormalCoverage() and CNV.countDeletionsMatchedByNormalCoverage()
		common_prefix = '~/script/variation/data/CNV/call_53_WABQN_b100_j50_lts_98_arrays_GADA2_A0.5T4.0M5'
		input_fname = os.path.expanduser('%s.tsv'%common_prefix)
		input_fname_ls = [input_fname]
		ecotype_id_data_source_id_ls = [(6932, 13)]	# use LerDeletionInOnly1Contig
		min_QC_segment_size = None
		min_no_of_probes = 5
		count_embedded_segment_as_match = True
		maxNormalPercInDeletion=0.4
		min_reciprocal_overlap = 0.000000001
		output_fname_prefix=os.path.expanduser('~/tmp/LerContig/%s_NoOfProbesPKB_'%(os.path.basename(common_prefix)))
		for ecotype_id, data_source_id in ecotype_id_data_source_id_ls:
			for deletion_cutoff in [-0.33, -0.5]:
				for maxNormalPercInDeletion in [0.2, 0.4, 0.6]:	#, 0.8]:
					CNV.countNoOfCNVDeletionsMatchQC(db_250k, input_fname_ls, \
							output_fname_prefix=output_fname_prefix,\
							ecotype_id=ecotype_id, data_source_id=data_source_id, \
							cnv_type_id=1,\
							min_QC_segment_size=min_QC_segment_size, deletion_cutoff=deletion_cutoff, \
							min_no_of_probes=min_no_of_probes,\
							count_embedded_segment_as_match=count_embedded_segment_as_match,\
							min_reciprocal_overlap=min_reciprocal_overlap, outputDeletionWithinSpan=True,\
							maxNormalPercInDeletion=maxNormalPercInDeletion,\
							judgeTrueDeletionFunctor=CNV.isDeletionTrueBasedOnNormalCoverage,\
							afterMatchFunction=CNV.countDeletionsMatchedByNormalCoverage)
		sys.exit(0)
	"""
	
	@classmethod
	def updateCNVCallPercUnCoveredByLerContig(cls, db_250k, ecotype_id=6932, \
								qc_data_source_id=15, qc_cnv_method_id=1, qc_cnv_type_id=None, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=None, run_type=1,\
								ecotype_id_ls=None):
		"""
		2010-7-31
			run_type:
				1. update CNVCall.percUnCoveredByLerContig
				2. update CNVCall.fractionDeletedInPECoverageData
				3. update CNV.fractionNotCoveredByLyrata
					ecotype_id has to be -1
				4. update CNVArrayCall.fractionDeletedInPECoverageData
		2010-7-28
			rename cnv_type_id to qc_cnv_type_id
			
		2010-7-23
			add argument qc_cnv_method_id (to restrict CNVQCCall)
			add run_type
		2010-7-13
			add argument cnv_method_id to filter input.
				update CNVCall.percUnCoveredByLerContig only for the specified cnv_method_id
		2010-6-27
			update CNVCall.percUnCoveredByLerContig in db
		"""
		from pymodule.CNV import CNVSegmentBinarySearchTreeKey, leftWithinRightAlsoEqualCmp, rightWithinLeftAlsoEqualCmp
		
		"""2010-6-17
			This data_source_id 15 is nucmer, maxmatch coords results. Overlapping fragments have all been merged.
			min_reciprocal_overlap=0.0001 because if there's any overlapping between CNV segment and stored nodes,
				that node should be returned.
			ecotype_id2normal_data is used to disqaulify a CNV deletion if it has too many normal data within.
		"""
		cnvQCData = cls.getCNVQCDataFromDB(db_250k, data_source_id=qc_data_source_id, ecotype_id=ecotype_id, \
												cnv_type_id=qc_cnv_type_id, cnv_method_id=qc_cnv_method_id,\
												min_QC_segment_size=None, min_no_of_probes=None,\
												min_reciprocal_overlap=min_reciprocal_overlap, \
												cmpfn=rightWithinLeftAlsoEqualCmp)
		ecotype_id2normal_data = cnvQCData.ecotype_id2cnv_qc_call_data
		
		from pymodule import PassingData
		param_obj = PassingData(no_of_valid_deletions=0, cnv_qc_call_id_set=set(), array_id2qc_data={}, \
							ecotype_id2normal_data=ecotype_id2normal_data, ecotype_id2qc_data=ecotype_id2normal_data, \
							maxNormalPercInDeletion=0.5, min_reciprocal_overlap=min_reciprocal_overlap,\
							minDeletedFraction=0.6)
		# maxNormalPercInDeletion=0.5 is a random number. CNV.isDeletionTrueBasedOnNormalCoverage() needs it.
		# minDeletedFraction is just a place holder because CNV.isDeletionTrueBasedOnDeletedFraction() needs it.
		sys.stderr.write("Updating CNVCall.percUnCoveredByLerContig (or sth similar ) for method %s ecotype %s ... \n"%\
						(cnv_method_id, ecotype_id))
		i = 0
		block_size = 10000
		real_counter = 0
		import Stock_250kDB
		if run_type==1 or run_type==2:
			TableClass = Stock_250kDB.CNVCall
		elif run_type==3:
			TableClass = Stock_250kDB.CNV
		elif run_type==4:
			TableClass = Stock_250kDB.CNVArrayCall
		else:
			sys.stderr.write("run_type %s not  supported yet.\n"%run_type)
			return
		query = TableClass.query
		if cnv_method_id is not None:	#2010-7-13
			query = query.filter_by(cnv_method_id=cnv_method_id)
		if ecotype_id is not None and run_type!=3:
			query = query.filter(TableClass.array.has(maternal_ecotype_id=ecotype_id))
		elif ecotype_id_ls and run_type!=3:
			query = query.filter(TableClass.array.has(Stock_250kDB.ArrayInfo.maternal_ecotype_id.in_(ecotype_id_ls)))
		
		rows = query.offset(i).limit(block_size)
		session = db_250k.session
		while rows.count()!=0:
			for row in rows:
				if run_type==3:
					ecotype_id = -1	#fake the ID to be same as lyrata's ecotype_id
				else:
					ecotype_id = row.array.maternal_ecotype_id
				if ecotype_id in ecotype_id2normal_data:
					if run_type==4:
						dbObjWithChr = row.cnv
					else:
						dbObjWithChr = row
					cnvSegmentKey = CNVSegmentBinarySearchTreeKey(chromosome=dbObjWithChr.chromosome, \
												span_ls=[dbObjWithChr.start, dbObjWithChr.stop],\
												min_reciprocal_overlap=min_reciprocal_overlap)
					param_obj.cnvSegmentKey = cnvSegmentKey
					if run_type==1:
						judgeData = cls.isDeletionTrueBasedOnNormalCoverage(param_obj=param_obj, ecotype_id=ecotype_id)
					elif run_type==2 or run_type==4:
						judgeData = cls.isDeletionTrueBasedOnDeletedFraction(param_obj=param_obj, ecotype_id=ecotype_id)
					elif run_type==3:
						judgeData = cls.isDeletionTrueBasedOnNormalCoverage(param_obj=param_obj, ecotype_id=ecotype_id)
					if judgeData:
						if run_type==1:
							percUnCoveredByLerContig = 1-judgeData.normalPerc
							if row.percUnCoveredByLerContig != percUnCoveredByLerContig:
								row.percUnCoveredByLerContig = percUnCoveredByLerContig
								session.add(row)
								real_counter += 1
						elif run_type==2 or run_type==4:
							if row.fractionDeletedInPECoverageData != judgeData.deletedFraction:
								row.fractionDeletedInPECoverageData = judgeData.deletedFraction
								session.add(row)
								real_counter += 1
						elif run_type==3:
							fractionNotCoveredByLyrata = 1 - judgeData.normalPerc
							if row.fractionNotCoveredByLyrata != fractionNotCoveredByLyrata:
								row.fractionNotCoveredByLyrata = fractionNotCoveredByLyrata
								session.add(row)
								real_counter += 1
						
				i += 1
			if real_counter%5000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, i, real_counter))
				session.flush()
				#session.expunge_all()
			rows = query.offset(i).limit(block_size)
		sys.stderr.write("%s\t%s\t%s Done.\n"%('\x08'*80, i, real_counter))
		session.flush()
		try:
			session.commit()
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
		
	"""
		# 2010-6-28 update everything in CNVCall
		CNV.updateCNVCallPercUnCoveredByLerContig(db_250k, ecotype_id=6932, \
								qc_data_source_id=15, \
								min_reciprocal_overlap=0.0000001)
		sys.exit(0)
		
		#2010-7-12 update only CNVCall whose cnv_method_id is 10
		CNV.updateCNVCallPercUnCoveredByLerContig(db_250k, ecotype_id=6932, \
								qc_data_source_id=15, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=10)
		sys.exit(0)
		
		#2010-7-19 update only CNVCall from specific cnv_method_id
		CNV.updateCNVCallPercUnCoveredByLerContig(db_250k, ecotype_id=6932, \
								qc_data_source_id=15, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=12)
		sys.exit(0)
		
		#2010-7-16 for the TAIR9 "-z banyan.usc.edu"
		CNV.updateCNVCallPercUnCoveredByLerContig(db_250k, ecotype_id=6932, \
								qc_data_source_id=14, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=8)
		sys.exit(0)
		
		#2010-7-23 for the TAIR9 "-z banyan.usc.edu". ecotype_id=None would exhaust hard disk storage in /tmp/ due to multiple joins.
		# Run CNV.updateCNVCallFractionDeletedInPECoverageData() instead.
		CNV.updateCNVCallPercUnCoveredByLerContig(db_250k, ecotype_id=None, \
								qc_data_source_id=13, qc_cnv_method_id=9, qc_cnv_type_id=1, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=7, run_type=2)
		sys.exit(0)
		
		#2010-7-28 update CNVCall.percUnCoveredByLerContig for the TAIR9 "-z banyan.usc.edu", cnv_method_id 11
		CNV.updateCNVCallPercUnCoveredByLerContig(db_250k, ecotype_id=6932, \
								qc_data_source_id=14, qc_cnv_type_id=None, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=11)
		sys.exit(0)
		
		#2010-7-31 update CNV.fractionNotCoveredByLyrata for the TAIR9 "-z banyan.usc.edu", cnv_method_id 18
		CNV.updateCNVCallPercUnCoveredByLerContig(db_250k, ecotype_id=-1, \
								qc_data_source_id=15, qc_cnv_method_id=1, qc_cnv_type_id=None, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=18, run_type=3)
		sys.exit(0)
		
		#2010-7-31 update CNV.fractionNotCoveredByLyrata for the TAIR9 "-z banyan.usc.edu", cnv_method_id 18
		CNV.updateCNVCallPercUnCoveredByLerContig(db_250k, ecotype_id=-1, \
								qc_data_source_id=15, qc_cnv_method_id=1, qc_cnv_type_id=None, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=18, run_type=3)
		sys.exit(0)
	"""
	
	@classmethod
	def updateCNVFractionNotCoveredByLyrata(cls, db_250k, ecotype_id=-1, \
								qc_data_source_id=15, qc_cnv_method_id=1, qc_cnv_type_id=None, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=None, \
								ecotype_id_ls=None):
		"""
		2010-8-2
			call CNV.updateCNVCallPercUnCoveredByLerContig() with run_type=3
		"""
		CNV.updateCNVCallPercUnCoveredByLerContig(db_250k, ecotype_id=-1, \
								qc_data_source_id=qc_data_source_id, qc_cnv_method_id=qc_cnv_method_id,\
								qc_cnv_type_id=None, \
								min_reciprocal_overlap=min_reciprocal_overlap, cnv_method_id=cnv_method_id,\
								run_type=3)
	"""
		#2010-8-3
		CNV.updateCNVFractionNotCoveredByLyrata(db_250k, ecotype_id=-1, \
								qc_data_source_id=16, qc_cnv_method_id=1, qc_cnv_type_id=None, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=18)
		sys.exit(0)
		
		#2010-8-3
		CNV.updateCNVFractionNotCoveredByLyrata(db_250k, ecotype_id=-1, \
								qc_data_source_id=16, qc_cnv_method_id=1, qc_cnv_type_id=None, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=20)
		sys.exit(0)
		
		#2010-8-3	qc_data_source_id 16 is too much for mysql to handle, exhaust disk usage in /tmp
		CNV.updateCNVFractionNotCoveredByLyrata(db_250k, ecotype_id=-1, \
								qc_data_source_id=16, qc_cnv_method_id=1, qc_cnv_type_id=None, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=22)
		sys.exit(0)
	"""
		
	@classmethod
	def updateCNVCallFractionDeletedInPECoverageData(cls, db_250k, ecotype_id=6932, \
								qc_data_source_id=15, qc_cnv_method_id=1, qc_cnv_type_id=1, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=None):
		"""
		2010-7-13
			add argument cnv_method_id to filter input.
				update CNVCall.percUnCoveredByLerContig only for the specified cnv_method_id
		2010-6-27
			update CNVCall.percUnCoveredByLerContig in db
		"""
		if ecotype_id is not None:
			ecotype_id_ls = [ecotype_id]
		else:
			ecotype_id_ls = cls.getEcotypeIDListInCNVQCAccessionGivenDataSourceID(qc_data_source_id)
		for ecotype_id in ecotype_id_ls:
			#2010-7-23 for the TAIR9 "-z banyan.usc.edu"
			CNV.updateCNVCallPercUnCoveredByLerContig(db_250k, ecotype_id=ecotype_id, \
								qc_data_source_id=qc_data_source_id, qc_cnv_method_id=qc_cnv_method_id,\
								qc_cnv_type_id=qc_cnv_type_id, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=cnv_method_id, run_type=2)
	"""
		#2010-7-23 for the TAIR9 "-z banyan.usc.edu"
		CNV.updateCNVCallFractionDeletedInPECoverageData(db_250k, ecotype_id=None, \
								qc_data_source_id=13, qc_cnv_method_id=9, qc_cnv_type_id=1, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=7)
		sys.exit(0)
		
		#2010-7-23 for the TAIR9 "-z banyan.usc.edu"
		CNV.updateCNVCallFractionDeletedInPECoverageData(db_250k, ecotype_id=None, \
								qc_data_source_id=13, qc_cnv_method_id=9, qc_cnv_type_id=1, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=7)
		sys.exit(0)
		
		#2010-8-6 for the TAIR9 "-z banyan.usc.edu"
		for cnv_method_id in [27]:
			CNV.updateCNVCallFractionDeletedInPECoverageData(db_250k, ecotype_id=None, \
								qc_data_source_id=13, qc_cnv_method_id=9, qc_cnv_type_id=1, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=cnv_method_id)
		sys.exit(0)
		
	"""
	
	@classmethod
	def updateCNVArrayCallFractionDeletedInPECoverageData(cls, db_250k, ecotype_id=None, \
								qc_data_source_id=13, qc_cnv_method_id=9, qc_cnv_type_id=1, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=None, cnv_method_id_ls=[]):
		"""
		# 2010-8-4
			add argument cnv_method_id_ls
		2010-08-1
			similar to CNV.updateCNVCallFractionDeletedInPECoverageData()
			but updates CNVArrayCall.fractionDeletedInPECoverageData
		"""
		
		
		if ecotype_id is not None:
			ecotype_id_ls = [ecotype_id]
		else:
			ecotype_id_ls = cls.getEcotypeIDListInCNVQCAccessionGivenDataSourceID(qc_data_source_id)
		# 2010-8-4
		if not cnv_method_id_ls:	#empty or None
			cnv_method_id_ls = []
		if cnv_method_id is not None:
			cnv_method_id_ls.append(cnv_method_id)
		
		for cnv_method_id in cnv_method_id_ls:
			for ecotype_id in ecotype_id_ls:
				#2010-7-23 for the TAIR9 "-z banyan.usc.edu"
				CNV.updateCNVCallPercUnCoveredByLerContig(db_250k, ecotype_id=ecotype_id, \
									qc_data_source_id=qc_data_source_id, qc_cnv_method_id=qc_cnv_method_id,\
									qc_cnv_type_id=qc_cnv_type_id, \
									min_reciprocal_overlap=0.0000001, cnv_method_id=cnv_method_id, run_type=4)
	"""
		#2010-8-1 for the TAIR9 "-z banyan.usc.edu"
		CNV.updateCNVArrayCallFractionDeletedInPECoverageData(db_250k, ecotype_id=None, \
								qc_data_source_id=13, qc_cnv_method_id=9, qc_cnv_type_id=1, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=18)
		sys.exit(0)
		
		#2010-8-1 for the TAIR9 "-z banyan.usc.edu"
		CNV.updateCNVArrayCallFractionDeletedInPECoverageData(db_250k, ecotype_id=None, \
								qc_data_source_id=13, qc_cnv_method_id=9, qc_cnv_type_id=1, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=19)
		sys.exit(0)
		
		#2010-8-1 for the TAIR9 "-z banyan.usc.edu"
		CNV.updateCNVArrayCallFractionDeletedInPECoverageData(db_250k, ecotype_id=None, \
								qc_data_source_id=13, qc_cnv_method_id=9, qc_cnv_type_id=1, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=22)
		sys.exit(0)
		
	"""
	
	class ProcessCNVForFrequencyUpdate(object):
		"""
		2010-8-2
			used by updateCNVFrequency to be plugged into traverseCNV()
		"""
		def __init__(self, **keywords):
			"""
			keywords shall include db_250k, no_of_total_arrays
			"""
			self.real_counter = 0
			
			# 2010-7-29
			for keyword, value in keywords.iteritems():
				setattr(self, keyword, value)
		
		def run(self, row, param_obj=None):
			"""
			"""
			array_id_set = set()
			for cnv_array_call in row.cnv_array_call_ls:
				array_id_set.add(cnv_array_call.array_id)
			frequency = len(array_id_set)/float(self.no_of_total_arrays)
			if row.frequency != frequency:
				row.frequency = frequency
				self.db_250k.session.add(row)
				self.real_counter += 1
				if self.real_counter%5000==0:
					session = self.db_250k.session
					session.flush()
					"""
					#session.expunge_all()	#2010-8-3 this shouldn't be executed. and session.expire_on_commit = False
					otherwise, error like this will pop out.
					
					sqlalchemy.orm.exc.DetachedInstanceError: Parent instance <CNV at 0x4c53510> is not bound to a Session;
					"""
	
	@classmethod
	def updateCNVFrequency(cls, db_250k, cnv_method_id=None, run_type=1, cnv_method_id_ls=[]):
		"""
		2010-8-1
			CNV.frequency was set wrongly because of a bug in CNVMergeAcrossArrays.py
		"""
		# 2010-8-4
		if not cnv_method_id_ls:	#empty or None
			cnv_method_id_ls = []
		if cnv_method_id is not None:
			cnv_method_id_ls.append(cnv_method_id)
		
		for cnv_method_id in cnv_method_id_ls:
			sys.stderr.write("Updating CNV.frequency for method %s... \n"%cnv_method_id)
			import Stock_250kDB
			
			rows = db_250k.metadata.bind.execute("select distinct array_id from cnv_array_call where cnv_method_id=%s"%\
											cnv_method_id)
			no_of_total_arrays = rows.rowcount	#not count(). different from direct Table query.
			
			session = db_250k.session
			session.expire_on_commit = False
			processClass = cls.ProcessCNVForFrequencyUpdate
			processor = processClass(no_of_total_arrays=no_of_total_arrays, db_250k=db_250k)
			
			cls.traverseCNV(db_250k, cnv_method_id=cnv_method_id, processClassIns=processor, \
						param_obj=None)
			
			session.flush()
			try:
				session.commit()
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
	
	"""
	#2010-8-1 "-z banyan"
	CNV.updateCNVFrequency(db_250k, cnv_method_id=18, run_type=1)
	CNV.updateCNVFrequency(db_250k, cnv_method_id=19, run_type=1)
	sys.exit(0)
	"""
	
	@classmethod
	def traverseCNV(cls, db_250k, cnv_method_id=None, processClassIns=None, param_obj=None):
		"""
		2010-8-2
		"""
		sys.stderr.write("Traversing CNV for method %s ... \n"%cnv_method_id)
		import Stock_250kDB
		
		i = 0
		block_size = 10000
		real_counter = 0
		import Stock_250kDB
		TableClass = Stock_250kDB.CNV
		query = TableClass.query
		if cnv_method_id is not None:	#2010-7-13
			query = query.filter_by(cnv_method_id=cnv_method_id)
		
		rows = query.offset(i).limit(block_size)
		session = db_250k.session
		while rows.count()!=0:
			for row in rows:
				processClassIns.run(row, param_obj)
				i += 1
			if i%5000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, i, processClassIns.real_counter))
			rows = query.offset(i).limit(block_size)
		sys.stderr.write("%s\t%s\t%s Done.\n"%('\x08'*80, i, processClassIns.real_counter))
	
	class ProcessCNVForProbeInfoUpdate(object):
		"""
		2010-8-2
			used by CNV.updateCNVProbeInfo() to be plugged into traverseCNV()
		"""
		def __init__(self, **keywords):
			"""
			keywords shall include db_250k, probeRBDict, min_reciprocal_overlap
			"""
			self.real_counter = 0
			
			# 2010-7-29
			for keyword, value in keywords.iteritems():
				setattr(self, keyword, value)
		
		def run(self, row, param_obj=None):
			"""
			row is Stock_250kDB.CNV
			"""
			from pymodule.CNV import CNVSegmentBinarySearchTreeKey
			from pymodule.RBTree import RBDict
			
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=row.chromosome, span_ls=[row.start, row.stop], \
									min_reciprocal_overlap=self.min_reciprocal_overlap, )
			probeNodes = []
			self.probeRBDict.findNodes(segmentKey, probeNodes)
			no_of_probes_covered = len(probeNodes)
			
			chr_pos_RBNode_ls = []
			for i in range(no_of_probes_covered):
				probeRBNode = probeNodes[i]
				chr_pos_RBNode_ls.append((probeRBNode.key.chromosome, probeRBNode.key.start, probeRBNode.key.stop, \
										probeRBNode.key.probe_id))
			chr_pos_RBNode_ls.sort()
			
			row.no_of_probes_covered = no_of_probes_covered
			if len(chr_pos_RBNode_ls)>0:
				row.start_probe_id = chr_pos_RBNode_ls[0][-1]
				row.stop_probe_id = chr_pos_RBNode_ls[-1][-1]
			self.db_250k.session.add(row)
			self.real_counter += 1
			if self.real_counter%5000==0:
				session = self.db_250k.session
				session.flush()	#this would cause session to become invalid (Instance error ...)
				#session.expunge_all()
	
	@classmethod
	def getProbeRBDict(cls, db_250k, probeType=2, min_reciprocal_overlap=0.000001):
		"""
		2010-8-2
			construct a RBDict with each node representing a probe. span_ls = [position, position].
			probeType (same as run_type of DB_250k2Array.get_probes())
				1: SNP probes
				2/4: CNV probes 	(no need to specify function argument snps. Tair9Copy=1)
				5: QC probes
				else: all probes (SNP + CNV + QC)
			
			example used in CNV.updateCNVProbeInfo()
		"""
		sys.stderr.write("Constructing probe RBDict for probeType %s ... \n"%probeType)
		
		import Stock_250kDB
		from DB_250k2Array import DB_250k2Array
		
		#DB_250k2Array.debug = True	#2010-8-2 temporary debug purpose
		
		probes, xy_ls, chr_pos_ls, probe_id_ls = DB_250k2Array.get_probes(db_250k.metadata.bind, \
										Stock_250kDB.Probes.table.name, snps=None, run_type=2, \
										constructChrPos2index=False, constructXY2index=False, need_xy_ls=False,\
										need_probes=False)[:4]
		from pymodule.CNV import CNVSegmentBinarySearchTreeKey
		from pymodule.RBTree import RBDict
		probeRBDict = RBDict()
		no_of_probes = len(probe_id_ls)
		no_of_probes_ejected = 0
		for i in xrange(no_of_probes):
			probe_id = probe_id_ls[i]
			chr, pos = chr_pos_ls[i]
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=chr, span_ls=[pos, pos], \
									min_reciprocal_overlap=min_reciprocal_overlap, probe_id=probe_id)
			if segmentKey not in probeRBDict:
				probeRBDict[segmentKey] = None
			else:
				no_of_probes_ejected += 1
			if i%5000==0:
				sys.stderr.write("%s%s\t%s ejected"%('\x08'*80, i, no_of_probes_ejected))
		del chr_pos_ls, probe_id_ls
		sys.stderr.write("%s%s\t%s ejected. %s items in probeRBDict\n"%('\x08'*80, i, no_of_probes_ejected, len(probeRBDict)))
		return probeRBDict
	
	@classmethod
	def updateCNVProbeInfo(cls, db_250k, cnv_method_id=None, cnv_method_id_ls=[]):
		"""
		2010-8-2
			update start_probe_id, stop_probe_id, no_of_probes_covered on entries in table CNV
		"""
		# 2010-8-4
		if not cnv_method_id_ls:	#empty or None
			cnv_method_id_ls = []
		if cnv_method_id is not None:
			cnv_method_id_ls.append(cnv_method_id)
		
		sys.stderr.write("Updating start_probe_id, stop_probe_id, no_of_probes_covered of table CNV method %s ...\n"%\
						repr(cnv_method_id_ls))
		import Stock_250kDB
		probeRBDict = cls.getProbeRBDict(db_250k, probeType=2, min_reciprocal_overlap=0.000001)
		
		
		session = db_250k.session
		for cnv_method_id in cnv_method_id_ls:
			session.expire_on_commit = False
			processClass = cls.ProcessCNVForProbeInfoUpdate
			processor = processClass(probeRBDict=probeRBDict, db_250k=db_250k, min_reciprocal_overlap=0.000001)
			
			cls.traverseCNV(db_250k, cnv_method_id=cnv_method_id, processClassIns=processor, \
						param_obj=None)
			
			session.flush()
			try:
				session.commit()
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
	
	"""
		#2010-8-2
		CNV.updateCNVProbeInfo(db_250k, cnv_method_id_ls=[18, 19])
		sys.exit(0)
		
		#2010-8-1 "-z banyan"
		#CNV.updateCNVFrequency(db_250k, cnv_method_id=22, run_type=1)
		CNV.updateCNVProbeInfo(db_250k, cnv_method_id_ls=[22])
		sys.exit(0)
		
		#2010-8-2
		CNV.updateCNVProbeInfo(db_250k, cnv_method_id_ls=[18, 19])
		sys.exit(0)
		
	"""
	
	@classmethod
	def traverseCNVQCCall(cls, db_250k, cnv_method_id=8, processClassIns=None, param_obj=None):
		"""
		2010-11-24
			a common function which traverses through entries in table CNVQCCall.
		"""
		i = 0
		block_size = 10000
		real_counter = 0
		import Stock_250kDB
		TableClass = Stock_250kDB.CNVQCCall
		query = TableClass.query
		if cnv_method_id is not None:	#2010-7-13
			query = query.filter_by(cnv_method_id=cnv_method_id)
		
		rows = query.offset(i).limit(block_size)
		session = db_250k.session
		while rows.count()!=0:
			for row in rows:
				processClassIns.run(row, param_obj)
				i += 1
			if i%5000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, i, processClassIns.real_counter))
			rows = query.offset(i).limit(block_size)
		sys.stderr.write("%s\t%s\t%s Done.\n"%('\x08'*80, i, processClassIns.real_counter))
	
	class ProcessCNVQCCallForProbeInfoUpdate(object):
		"""
		2010-8-2
			used by CNV.updateCNVProbeInfo() to be plugged into traverseCNV()
		"""
		def __init__(self, **keywords):
			"""
			keywords shall include db_250k, probeRBDict, min_reciprocal_overlap
			"""
			self.real_counter = 0
			
			# 2010-7-29
			for keyword, value in keywords.iteritems():
				setattr(self, keyword, value)
		
		def run(self, row, param_obj=None):
			"""
			row is Stock_250kDB.CNV
			"""
			from pymodule.CNV import CNVSegmentBinarySearchTreeKey
			from pymodule.RBTree import RBDict
			
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=row.chromosome, span_ls=[row.start, row.stop], \
									min_reciprocal_overlap=self.min_reciprocal_overlap, )
			probeNodes = []
			self.probeRBDict.findNodes(segmentKey, probeNodes)
			no_of_probes_covered = len(probeNodes)
			
			row.no_of_probes_covered = no_of_probes_covered
			"""
			chr_pos_RBNode_ls = []
			for i in range(no_of_probes_covered):
				probeRBNode = probeNodes[i]
				chr_pos_RBNode_ls.append((probeRBNode.key.chromosome, probeRBNode.key.start, probeRBNode.key.stop, \
										probeRBNode.key.probe_id))
			chr_pos_RBNode_ls.sort()
			
			if len(chr_pos_RBNode_ls)>0:
				row.start_probe_id = chr_pos_RBNode_ls[0][-1]
				row.stop_probe_id = chr_pos_RBNode_ls[-1][-1]
			"""
			
			self.db_250k.session.add(row)
			self.real_counter += 1
			if self.real_counter%5000==0:
				session = self.db_250k.session
				session.flush()	#this would cause session to become invalid (Instance error ...)
				#session.expunge_all()
	
	@classmethod
	def updateCNVQCCallProbeInfo(cls, db_250k, cnv_method_id=None, cnv_method_id_ls=[]):
		"""
		2010-11-24
			update no_of_probes_covered on entries in table CNVQCCall
		"""
		# 2010-8-4
		if not cnv_method_id_ls:	#empty or None
			cnv_method_id_ls = []
		if cnv_method_id is not None:
			cnv_method_id_ls.append(cnv_method_id)
		
		sys.stderr.write("Updating no_of_probes_covered of entries in CNVQCCall, cnv method %s ...\n"%\
						repr(cnv_method_id_ls))
		import Stock_250kDB
		probeRBDict = cls.getProbeRBDict(db_250k, probeType=2, min_reciprocal_overlap=0.000001)
		
		session = db_250k.session
		session.expire_on_commit = False
		processClass = cls.ProcessCNVQCCallForProbeInfoUpdate
		
		for cnv_method_id in cnv_method_id_ls:
			processor = processClass(probeRBDict=probeRBDict, db_250k=db_250k, min_reciprocal_overlap=0.000001)
			
			cls.traverseCNVQCCall(db_250k, cnv_method_id=cnv_method_id, processClassIns=processor, \
						param_obj=None)
			
			session.flush()
			try:
				session.commit()
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
	
	"""
		#2010-11-24
		CNV.updateCNVQCCallProbeInfo(db_250k, cnv_method_id_ls=[9])
		sys.exit(0)
		
		#2010-11-24
		CNV.updateCNVQCCallProbeInfo(db_250k, cnv_method_id_ls=[31])
		sys.exit(0)
		
	"""
	
	@classmethod
	def plotCNVCallFDRVsNoOfProbesBasedOnPercUnCoveredByLerContig(cls, db_250k, output_dir=None, \
											ecotype_id=6932, minPercUnCoveredByLerContig=0.6, \
											cnv_method_id=6, useProbeDensity=False, run_type=1, minScore=None):
		"""
		2010-8-6
			add argument run_type
				1. update CNVCall.percUnCoveredByLerContig
				2. update CNVCall.fractionDeletedInPECoverageData
			add argument minScore
		2010-7-28
			add argument output_dir
		2010-7-19
			add argument useProbeDensity.
				If true, the number-of-probes variable will be replaced by number-of-probes-per-kb.
		2010-7-18
			This function is written to inspect FDR-vs-no of probes curve for the predicted CNV calls. 
			
			FDR = False discovery rate = FP/(TP+FP)
			
			This function pulls all entries from CNVCall (no amplitude cutoff), then counts one as true deletion
				or false deletion based on minPercUnCoveredByLerContig.
			
			The argument cnv_method_id is intended to be the method encompassing the CNV calls that were predicted
				by CNVPredictDeletionBySVM.py although it could be any method.
			
			
		"""
		sys.stderr.write("Plotting FDR vs no of probes for CNVCall method-id=%s ... \n"%cnv_method_id)
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		from pymodule import PassingData
		param_obj = PassingData(no_of_valid_deletions=0, no_of_deletions=0, array_id2qc_data={}, \
				array_id2no_of_probes2qc_data = {}, \
				output_fname_prefix = os.path.join(output_dir, 'FDRVsNoOfProbes_cnvMethod%s_minNotCoveredFraction%s_useProbeDensity%s_minScore%s_runType%s'%\
						(cnv_method_id, minPercUnCoveredByLerContig, useProbeDensity, minScore, run_type)),\
				minPercUnCoveredByLerContig=minPercUnCoveredByLerContig,\
				array_id2label={})
		
		from common import get_ecotypeid2nativename
		ecotype_id2nativename = get_ecotypeid2nativename(db_250k.metadata.bind)
		
		i = 0
		block_size = 10000
		real_counter = 0
		import Stock_250kDB
		TableClass = Stock_250kDB.CNVCall
		query = TableClass.query
		if ecotype_id is not None:	#2010-8-6
			query = query.filter(TableClass.array.has(maternal_ecotype_id=ecotype_id))
		if cnv_method_id is not None:	#2010-7-13
			query = query.filter_by(cnv_method_id=cnv_method_id)
		if minScore is not None:
			query = query.filter(TableClass.probability>=minScore)
		if run_type==1:
			query = query.filter(TableClass.percUnCoveredByLerContig!=None)
		elif run_type==2:
			query = query.filter(TableClass.fractionDeletedInPECoverageData!=None)
		
		rows = query.offset(i).limit(block_size)
		session = db_250k.session
		while rows.count()!=0:
			for row in rows:
				array_id = row.array_id
				ecotype_id = row.array.maternal_ecotype_id
				nativename = ecotype_id2nativename.get(ecotype_id)
				param_obj.no_of_deletions += 1
				if array_id not in param_obj.array_id2qc_data:
					param_obj.array_id2qc_data[array_id] = PassingData(no_of_deletions=0, no_of_valid_deletions=0,\
																	segmentSize_ls=[])
					param_obj.array_id2no_of_probes2qc_data[array_id] = {}
					param_obj.array_id2label[array_id] = 'aid %s eid %s, %s method %s delCut %s minSco %s'%(array_id, \
											ecotype_id, nativename, cnv_method_id, minPercUnCoveredByLerContig, minScore)
				
				size_affected = row.stop - row.start + 1
				param_obj.array_id2qc_data[array_id].no_of_deletions += 1
				param_obj.array_id2qc_data[array_id].segmentSize_ls.append(size_affected)
				
				if useProbeDensity:
					no_of_probes_covered = int(row.no_of_probes_covered*1000/float(row.size_affected))
				else:
					no_of_probes_covered = row.no_of_probes_covered
				
				if no_of_probes_covered not in param_obj.array_id2no_of_probes2qc_data[array_id]:
					param_obj.array_id2no_of_probes2qc_data[array_id][no_of_probes_covered] = \
						PassingData(no_of_deletions=0, no_of_valid_deletions=0)
				param_obj.array_id2no_of_probes2qc_data[array_id][no_of_probes_covered].no_of_deletions += 1
				
				if run_type==1:
					percUnCoveredByLerContig = row.percUnCoveredByLerContig
				elif run_type==2:
					percUnCoveredByLerContig = row.fractionDeletedInPECoverageData
				else:
					percUnCoveredByLerContig = row.percUnCoveredByLerContig
				if percUnCoveredByLerContig>=minPercUnCoveredByLerContig:
					param_obj.no_of_valid_deletions += 1
					param_obj.array_id2qc_data[array_id].no_of_valid_deletions += 1
					param_obj.array_id2no_of_probes2qc_data[array_id][no_of_probes_covered].no_of_valid_deletions += 1
					real_counter += 1
				i += 1
			if i%5000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, i, real_counter))
			rows = query.offset(i).limit(block_size)
		sys.stderr.write("%s\t%s\t%s Done.\n"%('\x08'*80, i, real_counter))
		
		for array_id, qc_data in param_obj.array_id2qc_data.iteritems():
			totalSizeInKB = int(sum(param_obj.array_id2qc_data[array_id].segmentSize_ls)/1000.0)
			param_obj.array_id2label[array_id] += ' %sKb'%totalSizeInKB
		
		param_obj.ylabel = 'FDR'
		if useProbeDensity:
			param_obj.xlabel = 'No of probes per kb'
		else:
			param_obj.xlabel = 'No of probes'
		
		cls.outputFalsePositiveRate(param_obj)
	"""
		# 2010-7-18
		output_dir = os.path.expanduser('~/script/variation/data/CNV/FDRVsNoOfProbes/')
		CNV.plotCNVCallFDRVsNoOfProbesBasedOnPercUnCoveredByLerContig(db_250k, output_dir=output_dir, \
											ecotype_id=6932, minPercUnCoveredByLerContig=0.6, \
											cnv_method_id=8, useProbeDensity=True)
		sys.exit(0)
		
		# 2010-7-18 -z banyan.usc.edu
		output_dir = os.path.expanduser('~/script/variation/data/CNV/FDRVsNoOfProbes/')
		CNV.plotCNVCallFDRVsNoOfProbesBasedOnPercUnCoveredByLerContig(db_250k, output_dir=output_dir, \
											ecotype_id=6932, minPercUnCoveredByLerContig=0.6, \
											cnv_method_id=11, useProbeDensity=False)
		sys.exit(0)
		
		# 2010-8-6 -z banyan.usc.edu
		output_dir = os.path.expanduser('~/script/variation/data/CNV/FDRVsNoOfProbes/')
		CNV.plotCNVCallFDRVsNoOfProbesBasedOnPercUnCoveredByLerContig(db_250k, output_dir=output_dir, \
											ecotype_id=None, minPercUnCoveredByLerContig=0.4, \
											cnv_method_id=16, useProbeDensity=False, run_type=2, minScore=None)
		sys.exit(0)
		
		# 2010-8-6 -z banyan.usc.edu
		output_dir = os.path.expanduser('~/script/variation/data/CNV/FDRVsNoOfProbes/')
		for cnv_method_id in [27]:
			CNV.plotCNVCallFDRVsNoOfProbesBasedOnPercUnCoveredByLerContig(db_250k, output_dir=output_dir, \
											ecotype_id=None, minPercUnCoveredByLerContig=0.4, \
											cnv_method_id=cnv_method_id, useProbeDensity=False, run_type=2, minScore=None)
		sys.exit(0)
	"""
	
	@classmethod
	def plotOneFDR(cls, x_data2qc_data, output_fname_prefix, min_no_of_data_points=30, fig_title='', xlabel='No of Probes',\
				ylabel = 'FDR'):
		"""
		2010-10-24
		"""
		no_of_probes_ls = x_data2qc_data.keys()
		no_of_probes_ls.sort()
		no_of_non_valid_deletions_ls = []
		no_of_deletions_ls = []
		for no_of_probes in no_of_probes_ls:
			qc_data = x_data2qc_data[no_of_probes]
			no_of_valid_deletions = qc_data.no_of_valid_deletions
			no_of_deletions = qc_data.no_of_deletions
			no_of_non_valid_deletions = no_of_deletions-no_of_valid_deletions
			no_of_non_valid_deletions_ls.append(no_of_non_valid_deletions)
			no_of_deletions_ls.append(no_of_deletions)
			#false_positive_rate = no_of_non_valid_deletions/float(no_of_deletions)
			
			#sys.stderr.write("\t%s\t%s\t%s\t%s\n"%(no_of_probes, \
			#				no_of_non_valid_deletions, no_of_deletions, false_positive_rate))
		
		if output_fname_prefix:
			sys.stderr.write("Drawing the FDR plot ...")
			output_fname = '%s_FDR.png'%(output_fname_prefix)
			aggregateReturnData = cls.aggregateEnoughDataPoints(no_of_deletions_ls, \
							[no_of_non_valid_deletions_ls], [no_of_probes_ls],\
							min_no_of_data_points=min_no_of_data_points)
			
			new_no_of_deletions_ls = aggregateReturnData.no_of_data_points_ls
			new_no_of_non_valid_deletions_ls = aggregateReturnData.lists_to_be_summed[0]
			new_no_of_probes_ls = aggregateReturnData.lists_to_be_sampled[0]
			divide_functor = lambda x: x[0]/float(x[1])
			FPR_ls = map(divide_functor, zip(new_no_of_non_valid_deletions_ls, new_no_of_deletions_ls))
			import pylab
			pylab.clf()
			pylab.plot(new_no_of_probes_ls, FPR_ls, '.')
			# make a title for the figure
			pylab.title(fig_title, size='small')
			pylab.xlabel(xlabel)
			pylab.ylabel(ylabel)
			pylab.savefig(output_fname, dpi=300)
		
	@classmethod
	def plotCNVArrayCallFDRVsFrequencyNoOfProbes(cls, db_250k, output_dir=None, \
											ecotype_id=None, ecotype_id_ls=None, qc_data_source_id=13, \
											minNotCoveredFraction=0.6, \
											cnv_method_id=18, xDataType=1, draw2D=1, minScore=None, ):
		"""
		2010-8-3
			add argument draw2D
				1: FDR vs frequency by x-data
				2: histogram (count) vs frequency by x-data
				
				x-data depends on xDataType.
					1 or 2: log10(no_of_probes_covered)
					3: probeDensity (#probes per kb)
			
			add argument xDataType for 1D FDR figure
				1: FDR vs frequency
				2: FDR vs no_of_probes_covered
				3: FDR vs probeDensity
		2010-8-1
			argument qc_data_source_id is used to figure out which arrays to look at.
			
			similar to plotCNVCallFDRVsNoOfProbesBasedOnPercUnCoveredByLerContig()
			but for table CNVArrayCall and x-axis is frequency of CNV.
		"""
		sys.stderr.write("Plotting FDR vs frequency for CNVArrayCall method-id=%s ... \n"%cnv_method_id)
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		from pymodule import PassingData
		output_fname_prefix = os.path.join(output_dir, 'FDRVsXDataType%s_cnvMethod%s_minNotCoveredFraction%s_minScore%s'%\
						(xDataType, cnv_method_id,minNotCoveredFraction, minScore))
		param_obj = PassingData(no_of_valid_deletions=0, no_of_deletions=0, array_id2qc_data={}, \
				array_id2no_of_probes2qc_data = {}, \
				output_fname_prefix = output_fname_prefix,\
				minPercUnCoveredByLerContig=minNotCoveredFraction,\
				min_no_of_data_points=60,\
				array_id2label={})
		
		if ecotype_id is None and not ecotype_id_ls:
			ecotype_id_ls = cls.getEcotypeIDListInCNVQCAccessionGivenDataSourceID(qc_data_source_id)
		
		from common import get_ecotypeid2nativename
		ecotype_id2nativename = get_ecotypeid2nativename(db_250k.metadata.bind)
		
		i = 0
		block_size = 10000
		real_counter = 0
		import Stock_250kDB
		TableClass = Stock_250kDB.CNVArrayCall
		query = TableClass.query.filter(TableClass.fractionDeletedInPECoverageData!=None)
		if cnv_method_id is not None:	#2010-7-13
			query = query.filter_by(cnv_method_id=cnv_method_id)
		if ecotype_id is not None:
			query = query.filter(TableClass.array.has(maternal_ecotype_id=ecotype_id))
		elif ecotype_id_ls:
			query = query.filter(TableClass.array.has(Stock_250kDB.ArrayInfo.maternal_ecotype_id.in_(ecotype_id_ls)))
		
		if minScore is not None:
			query = query.filter(TableClass.score>=minScore)
		
		rows = query.offset(i).limit(block_size)
		session = db_250k.session
		x_data2qc_data = {}	#2010-10-24 to get FDR for all data combined
		
		while rows.count()!=0:
			for row in rows:
				array_id = row.array_id
				param_obj.no_of_deletions += 1
				if array_id not in param_obj.array_id2qc_data:
					ecotype_id = row.array.maternal_ecotype_id
					param_obj.array_id2qc_data[array_id] = PassingData(no_of_deletions=0, no_of_valid_deletions=0,\
										no_of_probes_ls=[], frequency_ls =[], falseDiscoveryOrNot_ls =[],\
										segmentSize_ls=[])
					param_obj.array_id2no_of_probes2qc_data[array_id] = {}
					
					nativename = ecotype_id2nativename.get(ecotype_id)
					param_obj.array_id2label[array_id] = 'aid %s eid %s, %s method %s delCut %s minSco %s'%(array_id, \
											ecotype_id, nativename, cnv_method_id, minNotCoveredFraction, minScore)
				
				param_obj.array_id2qc_data[array_id].no_of_deletions += 1
				
				no_of_probes_covered = row.cnv.no_of_probes_covered
				size_affected = row.cnv.stop - row.cnv.start + 1
				param_obj.array_id2qc_data[array_id].frequency_ls.append(row.cnv.frequency)
				param_obj.array_id2qc_data[array_id].segmentSize_ls.append(size_affected)
				if xDataType==1:
					x_data = row.cnv.frequency
					param_obj.array_id2qc_data[array_id].no_of_probes_ls.append(math.log10(no_of_probes_covered))
				elif xDataType==2:
					x_data = no_of_probes_covered
					param_obj.array_id2qc_data[array_id].no_of_probes_ls.append(math.log10(no_of_probes_covered))
				elif xDataType==3:
					x_data = int(no_of_probes_covered*1000/float(size_affected))
					param_obj.array_id2qc_data[array_id].no_of_probes_ls.append(x_data)
				
				if x_data not in param_obj.array_id2no_of_probes2qc_data[array_id]:
					param_obj.array_id2no_of_probes2qc_data[array_id][x_data] = \
						PassingData(no_of_deletions=0, no_of_valid_deletions=0)
				param_obj.array_id2no_of_probes2qc_data[array_id][x_data].no_of_deletions += 1
				
				if x_data not in x_data2qc_data:
					x_data2qc_data[x_data] = PassingData(no_of_deletions=0, no_of_valid_deletions=0)
				x_data2qc_data[x_data].no_of_deletions += 1
				
				if row.fractionDeletedInPECoverageData>=minNotCoveredFraction:
					param_obj.array_id2qc_data[array_id].falseDiscoveryOrNot_ls.append(0)
					param_obj.no_of_valid_deletions += 1
					param_obj.array_id2qc_data[array_id].no_of_valid_deletions += 1
					param_obj.array_id2no_of_probes2qc_data[array_id][x_data].no_of_valid_deletions += 1
					real_counter += 1
					
					x_data2qc_data[x_data].no_of_valid_deletions += 1
				else:
					param_obj.array_id2qc_data[array_id].falseDiscoveryOrNot_ls.append(1)
				i += 1
			if i%5000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, i, real_counter))
			rows = query.offset(i).limit(block_size)
		sys.stderr.write("%s\t%s\t%s Done.\n"%('\x08'*80, i, real_counter))
		
		param_obj.ylabel = 'FDR'
		
		if xDataType==1:
			param_obj.xlabel = 'CNV frequency'
			xlabel_in_2D = 'log10(no of probes)'
		elif xDataType==2:
			param_obj.xlabel = 'No of probes'
			xlabel_in_2D = 'log10(no of probes)'
		elif xDataType==3:
			param_obj.xlabel = 'No of probes per kb'
			xlabel_in_2D = param_obj.xlabel
		
		for array_id, qc_data in param_obj.array_id2qc_data.iteritems():
			totalSizeInKB = int(sum(param_obj.array_id2qc_data[array_id].segmentSize_ls)/1000.0)
			param_obj.array_id2label[array_id] += ' %sKb'%totalSizeInKB
		
		#2010-10-24 to get FDR for all data combined
		cls.plotOneFDR(x_data2qc_data, output_fname_prefix, min_no_of_data_points=3000, fig_title='', xlabel=param_obj.xlabel,\
				ylabel = param_obj.ylabel)
		return	#2010-10-24 temporary. ignore the rest
		
		cls.outputFalsePositiveRate(param_obj)
		
		import numpy
		if draw2D:
			for array_id, qc_data in param_obj.array_id2qc_data.iteritems():
				if len(qc_data.no_of_probes_ls)>0:
					x_ls = qc_data.no_of_probes_ls
					y_ls = qc_data.frequency_ls
					if draw2D==2:
						C_ls = [1]*len(y_ls)
						reduce_C_function = CNV.logSum
						colorBarLabel='log10(count)'
					else:
						C_ls = qc_data.falseDiscoveryOrNot_ls
						reduce_C_function = numpy.mean
						colorBarLabel='FDR'
					#
					fig_fname = os.path.join(output_dir, \
						'FDRVsFrequencyByNoOfProbes_cnvMethod%s_minNotCoveredFraction%s_xDataType%s_minScore%s_2D%s_array%s.png'%\
						(cnv_method_id, minNotCoveredFraction, xDataType, minScore, draw2D, array_id))
					sys.stderr.write("Plotting FDR vs frequency X no-of-probes for array %s..."%array_id)
					cls.drawHexbin(x_ls, y_ls, C_ls, fig_fname=fig_fname, gridsize=20, \
							title='%s'%param_obj.array_id2label.get(array_id), \
							xlabel = xlabel_in_2D, \
							ylabel = 'frequency',\
							colorBarLabel=colorBarLabel, reduce_C_function= reduce_C_function)
					sys.stderr.write("Done.\n")
	"""
		# 2010-8-1 -z banyan.usc.edu
		output_dir = os.path.expanduser('~/script/variation/data/CNV/FDRVsFrequency/')
		CNV.plotCNVArrayCallFDRVsFrequencyNoOfProbes(db_250k, output_dir=output_dir, \
											ecotype_id=None, qc_data_source_id=13, minNotCoveredFraction=0.6, \
											cnv_method_id=18, xDataType=1)
		sys.exit(0)
		
		# 2010-8-1 -z banyan.usc.edu (1D FDR vs frequency, 2D FDR vs frequency by no-of-probes)
		output_dir = os.path.expanduser('~/script/variation/data/CNV/FDRVsFrequency/')
		CNV.plotCNVArrayCallFDRVsFrequencyNoOfProbes(db_250k, output_dir=output_dir, \
											ecotype_id=None, qc_data_source_id=13, minNotCoveredFraction=0.4, \
											cnv_method_id=21, xDataType=1, draw2D=True)
		sys.exit(0)
		
		# 2010-8-1 -z banyan.usc.edu
		output_dir = os.path.expanduser('~/script/variation/data/CNV/FDRVsFrequency/')
		minScore = 0.95
		cnv_method_id = 21
		CNV.plotCNVArrayCallFDRVsFrequencyNoOfProbes(db_250k, output_dir=output_dir, \
								ecotype_id=None, qc_data_source_id=13, minNotCoveredFraction=0.4, \
								cnv_method_id=cnv_method_id, xDataType=1, draw2D=True, minScore=minScore)
		# 2010-8-1 -z banyan.usc.edu (just FDR vs no-of-probes)
		CNV.plotCNVArrayCallFDRVsFrequencyNoOfProbes(db_250k, output_dir=output_dir, \
								ecotype_id=None, qc_data_source_id=13, minNotCoveredFraction=0.4, \
								cnv_method_id=cnv_method_id, xDataType=2, draw2D=False, minScore=minScore)
		# 2010-8-1 -z banyan.usc.edu
		CNV.plotCNVArrayCallFDRVsFrequencyNoOfProbes(db_250k, output_dir=output_dir, \
								ecotype_id=None, qc_data_source_id=13, minNotCoveredFraction=0.4, \
								cnv_method_id=cnv_method_id, xDataType=3, draw2D=True, minScore=minScore)
		sys.exit(0)
		
		# 2010-8-1 -z banyan.usc.edu
		output_dir = os.path.expanduser('~/script/variation/data/CNV/FDRVsFrequency/')
		minScore = None
		cnv_method_id = 21
		draw2D = 1
		minNotCoveredFraction = 0.8
		CNV.plotCNVArrayCallFDRVsFrequencyNoOfProbes(db_250k, output_dir=output_dir, \
								ecotype_id=None, qc_data_source_id=13, minNotCoveredFraction=minNotCoveredFraction, \
								cnv_method_id=cnv_method_id, xDataType=1, draw2D=draw2D, minScore=minScore)
		CNV.plotCNVArrayCallFDRVsFrequencyNoOfProbes(db_250k, output_dir=output_dir, \
								ecotype_id=None, qc_data_source_id=13, minNotCoveredFraction=minNotCoveredFraction, \
								cnv_method_id=cnv_method_id, xDataType=2, draw2D=False, minScore=minScore)
		sys.exit(0)
		
		# 2010-8-1 -z banyan.usc.edu
		minScore = None
		cnv_method_id = 18
		output_dir = os.path.expanduser('~/script/variation/data/CNV/FDRVsFrequency/')
		CNV.plotCNVArrayCallFDRVsFrequencyNoOfProbes(db_250k, output_dir=output_dir, \
								ecotype_id=None, qc_data_source_id=13, minNotCoveredFraction=0.4, \
								cnv_method_id=cnv_method_id, xDataType=1, draw2D=True, minScore=minScore)
		CNV.plotCNVArrayCallFDRVsFrequencyNoOfProbes(db_250k, output_dir=output_dir, \
								ecotype_id=None, qc_data_source_id=13, minNotCoveredFraction=0.4, \
								cnv_method_id=cnv_method_id, xDataType=2, draw2D=False, minScore=minScore)
		sys.exit(0)
		
		# 2010-8-1 -z banyan.usc.edu
		minScore = None
		cnv_method_id = 19
		output_dir = os.path.expanduser('~/script/variation/data/CNV/FDRVsFrequency/')
		CNV.plotCNVArrayCallFDRVsFrequencyNoOfProbes(db_250k, output_dir=output_dir, \
								ecotype_id=None, qc_data_source_id=13, minNotCoveredFraction=0.4, \
								cnv_method_id=cnv_method_id, xDataType=1, draw2D=True, minScore=minScore)
		CNV.plotCNVArrayCallFDRVsFrequencyNoOfProbes(db_250k, output_dir=output_dir, \
								ecotype_id=None, qc_data_source_id=13, minNotCoveredFraction=0.4, \
								cnv_method_id=cnv_method_id, xDataType=2, draw2D=False, minScore=minScore)
		sys.exit(0)
		
		# 2010-8-1 -z banyan.usc.edu
		minScore = None
		cnv_method_id = 22
		minNotCoveredFraction = 0.4
		output_dir = os.path.expanduser('~/script/variation/data/CNV/FDRVsFrequency/')
		CNV.plotCNVArrayCallFDRVsFrequencyNoOfProbes(db_250k, output_dir=output_dir, \
								ecotype_id=None, qc_data_source_id=13, minNotCoveredFraction=minNotCoveredFraction, \
								cnv_method_id=cnv_method_id, xDataType=1, draw2D=True, minScore=minScore)
		CNV.plotCNVArrayCallFDRVsFrequencyNoOfProbes(db_250k, output_dir=output_dir, \
								ecotype_id=None, qc_data_source_id=13, minNotCoveredFraction=minNotCoveredFraction, \
								cnv_method_id=cnv_method_id, xDataType=2, draw2D=False, minScore=minScore)
		sys.exit(0)
		
		
		
		
		#2010-8-1 "-z banyan"
		CNV.updateCNVProbeInfo(db_250k, cnv_method_id_ls=[23, 24])
		#2010-8-1 for the TAIR9 "-z banyan.usc.edu"
		for cnv_method_id in [23, 24]:
			CNV.updateCNVArrayCallFractionDeletedInPECoverageData(db_250k, ecotype_id=None, \
								qc_data_source_id=13, qc_cnv_method_id=9, qc_cnv_type_id=1, \
								min_reciprocal_overlap=0.0000001, cnv_method_id=cnv_method_id)
			minScore = None
			minNotCoveredFraction = 0.4
			output_dir = os.path.expanduser('~/script/variation/data/CNV/FDRVsFrequency/')
			CNV.plotCNVArrayCallFDRVsFrequencyNoOfProbes(db_250k, output_dir=output_dir, \
									ecotype_id=None, qc_data_source_id=13, minNotCoveredFraction=minNotCoveredFraction, \
									cnv_method_id=cnv_method_id, xDataType=1, draw2D=True, minScore=minScore)
			CNV.plotCNVArrayCallFDRVsFrequencyNoOfProbes(db_250k, output_dir=output_dir, \
									ecotype_id=None, qc_data_source_id=13, minNotCoveredFraction=minNotCoveredFraction, \
									cnv_method_id=cnv_method_id, xDataType=2, draw2D=False, minScore=minScore)
		sys.exit(0)
		
		# 2010-8-1 -z banyan.usc.edu
		output_dir = os.path.expanduser('~/script/variation/data/CNV/FDRVsFrequency/')
		CNV.plotCNVArrayCallFDRVsFrequencyNoOfProbes(db_250k, output_dir=output_dir, \
											ecotype_id=None, qc_data_source_id=13, minNotCoveredFraction=0.4, \
											cnv_method_id=19, xDataType=2)
		sys.exit(0)
	"""
	
	@classmethod
	def getEcotypeIDListInCNVQCAccessionGivenDataSourceID(cls, data_source_id=13):
		"""
		2010-7-20
			called by plotCNVCallVsQC_OverlappingRatio2DHistogram()
			
		"""
		import sys
		sys.stderr.write("Getting list of ecotype IDs for data source %s from CNVQCAccession ..."%data_source_id)
		# get a list of ecotype ids which are present in CNVQCAccession
		import Stock_250kDB
		query = Stock_250kDB.CNVQCAccession.query.filter_by(data_source_id=data_source_id)
		ecotype_id_ls = []
		for row in query:
			ecotype_id_ls.append(row.ecotype_id)
		sys.stderr.write("%s ecotypes.\n"%(len(ecotype_id_ls)))
		return ecotype_id_ls
	
	@classmethod
	def getAccessionID2EcotypeIDFromCNVQCAccessionGivenDataSourceID(cls, db_250k, data_source_id=13):
		"""
		2010-8-4
			get a map from accession_id to ecotype_id from CNVQCAccession.
			This function is called by CNV.getCNVQCDataFromDB().
		"""
		import sys
		sys.stderr.write("Getting accession_id2ecotype_id for data source %s from CNVQCAccession ..."%data_source_id)
		import Stock_250kDB
		query = Stock_250kDB.CNVQCAccession.query.filter_by(data_source_id=data_source_id)
		accession_id2ecotype_id = {}
		for row in query:
			accession_id2ecotype_id[row.id] = row.ecotype_id
		sys.stderr.write("%s accessions.\n"%(len(accession_id2ecotype_id)))
		return accession_id2ecotype_id
	
	@classmethod
	def getAccessionIDLsFromCNVQCAccessionGivenEcotypeIDLs(cls, db_250k, ecotype_id_ls=[], data_source_id=None):
		"""
		2010-8-4
			get a list of accession IDs given ecotype id list from CNVQCAccession.
			This function is called by CNV.getCNVQCDataFromDB().
		"""
		import sys
		sys.stderr.write("Getting the list of accession ids for %s ecotypes, data source %s from CNVQCAccession ..."%\
						(len(ecotype_id_ls), data_source_id))
		import Stock_250kDB
		query = Stock_250kDB.CNVQCAccession.query.filter(Stock_250kDB.CNVQCAccession.ecotype_id.in_(ecotype_id_ls))
		if data_source_id is not None:	# 2010-8-4
			query = query.filter_by(data_source_id=data_source_id)
		accession_id_ls = []
		for row in query:
			accession_id_ls.append(row.id)
		sys.stderr.write("%s accessions.\n"%(len(accession_id_ls)))
		return accession_id_ls
	
	@classmethod
	def plotCNVCallVsQC_OverlappingRatio2DHistogram(cls, db_250k, cnv_type_id=None, cnv_method_id=None, output_dir=None,\
									qc_data_source_id=13, qc_cnv_method_id=3, ecotype_id=None, ecotype_id_ls = None, \
									min_reciprocal_overlap=0.000000001, minDeletedFraction=0.6, \
									useProbeDensity=False):
		"""
		2010-8-6
			This program involves computing fractionDeletedInPECoverageData on the fly.
			CNV.plotCNVCallFDRVsNoOfProbesBasedOnPercUnCoveredByLerContig() fetches the pre-computed
				fractionDeletedInPECoverageData from db and plot FDR vs no-of-probes.
		2010-7-22
			add argument ecotype_id_ls to restrict the ecotypes to look at.
				ecotype_id_ls is only used when ecotype_id is None.
		2010-7-19
			draw a FDR vs no of probes (or probe density) plot for each array
			draw a 2D histogram of overlap ratio (ratio in CNVQCAccession vs ratio in CNVCall)
			
		"""
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		from pymodule.CNV import CNVSegmentBinarySearchTreeKey, leftWithinRightAlsoEqualCmp, rightWithinLeftAlsoEqualCmp
		
		cnvQCData = cls.getCNVQCDataFromDB(db_250k, data_source_id=qc_data_source_id, \
										ecotype_id=ecotype_id, ecotype_id_ls=ecotype_id_ls, cnv_type_id=cnv_type_id, \
										min_QC_segment_size=None, min_no_of_probes=None,\
										min_reciprocal_overlap=min_reciprocal_overlap, \
										cmpfn=rightWithinLeftAlsoEqualCmp,\
										cnv_method_id=qc_cnv_method_id)
		ecotype_id2qc_data = cnvQCData.ecotype_id2cnv_qc_call_data
		
		if ecotype_id is None and not ecotype_id_ls:
			ecotype_id_ls = cls.getEcotypeIDListInCNVQCAccessionGivenDataSourceID(qc_data_source_id)
		
		from pymodule import PassingData
		import Stock_250kDB
		i = 0
		block_size = 10000
		real_counter = 0
		TableClass = Stock_250kDB.CNVCall
		query = TableClass.query
		if cnv_method_id is not None:	#2010-7-13
			query = query.filter_by(cnv_method_id=cnv_method_id)
		if ecotype_id is not None:
			query = query.filter(TableClass.array.has(maternal_ecotype_id=ecotype_id))
		elif ecotype_id_ls:
			query = query.filter(TableClass.array.has(Stock_250kDB.ArrayInfo.maternal_ecotype_id.in_(ecotype_id_ls)))
		
		rows = query.offset(i).limit(block_size)
		session = db_250k.session
		
		minPercUnCoveredByLerContig = minDeletedFraction
		param_obj = PassingData(no_of_valid_deletions=0, no_of_deletions=0, array_id2qc_data={}, \
				array_id2no_of_probes2qc_data = {}, \
				output_fname_prefix = output_dir+'FDRVsNoOfProbes_QCDataSource%s_QCCNVMethod%s_CNVMethod%s_minDeletedFraction%s_useProbeDensity%s'%\
						(qc_data_source_id, qc_cnv_method_id, cnv_method_id, minDeletedFraction, useProbeDensity),\
				minPercUnCoveredByLerContig=minPercUnCoveredByLerContig,\
				ecotype_id2qc_data=ecotype_id2qc_data, minDeletedFraction=minDeletedFraction)
		
		array_id2overlap_ratio_ls = {}
		array_id2ecotype_id = {}
		while rows.count()!=0:
			for row in rows:
				if row.array.maternal_ecotype_id in ecotype_id2qc_data:
					qc_data = ecotype_id2qc_data.get(row.array.maternal_ecotype_id)
					cnvSegmentKey = CNVSegmentBinarySearchTreeKey(chromosome=row.chromosome, span_ls=[row.start, row.stop],\
												min_reciprocal_overlap=min_reciprocal_overlap)
					
					param_obj.cnvSegmentKey = cnvSegmentKey
					judgeData = cls.isDeletionTrueBasedOnDeletedFraction(param_obj=param_obj, \
														ecotype_id=row.array.maternal_ecotype_id)
					
					array_id = row.array_id
					param_obj.no_of_deletions += 1
					array_id2ecotype_id[array_id] = row.array.maternal_ecotype_id
					
					if array_id not in param_obj.array_id2qc_data:
						param_obj.array_id2qc_data[array_id] = PassingData(no_of_deletions=0, no_of_valid_deletions=0)
						param_obj.array_id2no_of_probes2qc_data[array_id] = {}
					
					param_obj.array_id2qc_data[array_id].no_of_deletions += 1
					
					if useProbeDensity:
						no_of_probes_covered = int(row.no_of_probes_covered*1000/float(row.size_affected))
					else:
						no_of_probes_covered = row.no_of_probes_covered
					
					if no_of_probes_covered not in param_obj.array_id2no_of_probes2qc_data[array_id]:
						param_obj.array_id2no_of_probes2qc_data[array_id][no_of_probes_covered] = \
							PassingData(no_of_deletions=0, no_of_valid_deletions=0)
					param_obj.array_id2no_of_probes2qc_data[array_id][no_of_probes_covered].no_of_deletions += 1
					
					if judgeData and judgeData.isDeletionTrue:
						if array_id not in array_id2overlap_ratio_ls:
							array_id2overlap_ratio_ls[array_id] = []
						
						array_id2overlap_ratio_ls[array_id] += judgeData.overlap_ratio_ls
						param_obj.no_of_valid_deletions += 1
						param_obj.array_id2qc_data[array_id].no_of_valid_deletions += 1
						param_obj.array_id2no_of_probes2qc_data[array_id][no_of_probes_covered].no_of_valid_deletions += 1
						real_counter += 1
				i += 1
			if real_counter%100==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, i, real_counter))
			rows = query.offset(i).limit(block_size)
		sys.stderr.write("%s\t%s\t%s Done.\n"%('\x08'*80, i, real_counter))
		
		# get label for each array
		from common import get_ecotypeid2nativename
		ecotype_id2nativename = get_ecotypeid2nativename(db_250k.metadata.bind)
		array_id2label = {}
		for array_id, ecotype_id in array_id2ecotype_id.iteritems():
			nativename = ecotype_id2nativename.get(ecotype_id)
			label = 'a-id %s, e-id %s, %s'%(array_id, ecotype_id, nativename)
			array_id2label[array_id] = label
		
		param_obj.ylabel = 'FDR'
		if useProbeDensity:
			param_obj.xlabel = 'No of probes per kb'
		else:
			param_obj.xlabel = 'No of probes'
		param_obj.array_id2label = array_id2label
		
		cls.outputFalsePositiveRate(param_obj)
		
		for array_id, overlap_ratio_ls in array_id2overlap_ratio_ls.iteritems():
			if len(overlap_ratio_ls)>0:
				x_ls = [row[0] for row in overlap_ratio_ls]
				y_ls = [row[1] for row in overlap_ratio_ls]
				C_ls = [1]*len(y_ls)
				fig_fname = output_dir+'OverlapRatio2DHistogram_QCDataSource%s_QCCNVMethod%s_CNVMethod%s_useProbeDensity%s_array%s.png'%\
								(qc_data_source_id, qc_cnv_method_id, cnv_method_id, useProbeDensity, array_id)
				cls.drawHexbin(x_ls, y_ls, C_ls, fig_fname=fig_fname, gridsize=20, \
						title='%s overlap-ratio 2D histogram'%array_id2label.get(array_id), \
						xlabel='overlapRatioInCNVCall', \
						ylabel='overlapRatioInCNVQCCall',\
						colorBarLabel='log10(count)', reduce_C_function=CNV.logSum)
			
	
	"""
		# 2010-7-20	QC against BreakDancer-derived deletions
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVCallVsQC_OverlapRatio2DHistogram/')
		# only for banyan.usc.edu's stock_250k db, which is in TAIR9.
		CNV.plotCNVCallVsQC_OverlappingRatio2DHistogram(db_250k, cnv_type_id=1, cnv_method_id=8, \
									output_dir=output_dir, \
									qc_data_source_id=13, qc_cnv_method_id=3, \
									min_reciprocal_overlap=0.000000001, minDeletedFraction=0.4, useProbeDensity=False)
		sys.exit(0)
		
		# 2010-7-20	QC against BreakDancer-derived deletions
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVCallVsQC_OverlapRatio2DHistogram/')
		# only for banyan.usc.edu's stock_250k db, which is in TAIR9.
		CNV.plotCNVCallVsQC_OverlappingRatio2DHistogram(db_250k, cnv_type_id=1, cnv_method_id=8, \
									output_dir=output_dir, \
									qc_data_source_id=13, qc_cnv_method_id=3, ecotype_id=None,\
									min_reciprocal_overlap=0.000000001, minDeletedFraction=0.6, useProbeDensity=True)
		sys.exit(0)
		
		# 2010-7-22 QC against PE-data coverage-derived deletions
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVCallVsQC_OverlapRatio2DHistogram/')
		# only for banyan.usc.edu's stock_250k db, which is in TAIR9.
		ecotype_id_ls = [8256, 6900, 6901, 5830, 5831]
		CNV.plotCNVCallVsQC_OverlappingRatio2DHistogram(db_250k, cnv_type_id=1, cnv_method_id=8, \
									output_dir=output_dir, \
									qc_data_source_id=13, qc_cnv_method_id=9, ecotype_id_ls=ecotype_id_ls,\
									min_reciprocal_overlap=0.000000001, minDeletedFraction=0.6, \
									useProbeDensity=False)
		sys.exit(0)
		
		# 2010-7-22 QC against PE-data coverage-derived deletions
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVCallVsQC_OverlapRatio2DHistogram/')
		# only for banyan.usc.edu's stock_250k db, which is in TAIR9.
		ecotype_id_ls = [8256, 6900, 6901, 5830, 5831]
		CNV.plotCNVCallVsQC_OverlappingRatio2DHistogram(db_250k, cnv_type_id=1, cnv_method_id=8, \
									output_dir=output_dir, \
									qc_data_source_id=13, qc_cnv_method_id=3, ecotype_id_ls=ecotype_id_ls,\
									min_reciprocal_overlap=0.000000001, minDeletedFraction=0.95, \
									useProbeDensity=False)
		sys.exit(0)
		
		# 2010-7-28 QC against PE-data coverage-derived deletions
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVCallVsQC_OverlapRatio2DHistogram/')
		# only for banyan.usc.edu's stock_250k db, which is in TAIR9.
		qc_data_source_id = 13
		ecotype_id_ls = CNV.getEcotypeIDListInCNVQCAccessionGivenDataSourceID(qc_data_source_id)
		for ecotype_id in ecotype_id_ls:
			CNV.plotCNVCallVsQC_OverlappingRatio2DHistogram(db_250k, cnv_type_id=1, cnv_method_id=11, \
									output_dir=output_dir, \
									qc_data_source_id=qc_data_source_id, qc_cnv_method_id=9, ecotype_id=ecotype_id, \
									min_reciprocal_overlap=0.000000001, minDeletedFraction=0.6, \
									useProbeDensity=False)
		sys.exit(0)
		
	"""
	
	@classmethod
	def outputCNVQCVsCNVCallStats(cls, db_250k, qc_data_source_id=None, cnv_type_id=1, cnv_method_id=8, \
						output_fname=None, \
						ecotype_id=None, min_no_of_probes=None, min_QC_segment_size=None):
		"""
		2010-7-20
			output tsv data for pymodule.DataMatrixGuiXYProbe.py to view.
		"""
		from pymodule.CNV import CNVSegmentBinarySearchTreeKey, leftWithinRightAlsoEqualCmp, rightWithinLeftAlsoEqualCmp
		cnvQCData = cls.getCNVQCDataFromDB(db_250k, data_source_id=qc_data_source_id, \
												ecotype_id=ecotype_id, cnv_type_id=cnv_type_id, \
												min_QC_segment_size=None, min_no_of_probes=None,\
												min_reciprocal_overlap=0.0000001, \
												cmpfn=rightWithinLeftAlsoEqualCmp)
		ecotype_id2qc_data = cnvQCData.ecotype_id2cnv_qc_call_data
		ecotype_id2data = cnvQCData.ecotype_id2data
		
		ecotype_id_ls = cls.getEcotypeIDListInCNVQCAccessionGivenDataSourceID(qc_data_source_id)
		
		from pymodule import PassingData
		import Stock_250kDB
		i = 0
		block_size = 10000
		real_counter = 0
		TableClass = Stock_250kDB.CNVCall
		query = TableClass.query
		if cnv_method_id is not None:	#2010-7-13
			query = query.filter_by(cnv_method_id=cnv_method_id)
		if ecotype_id is not None:
			query = query.filter(TableClass.array.has(maternal_ecotype_id=ecotype_id))
		else:
			query = query.filter(TableClass.array.has(Stock_250kDB.ArrayInfo.maternal_ecotype_id.in_(ecotype_id_ls)))
		
		rows = query.offset(i).limit(block_size)
		
		from common import get_ecotypeid2nativename
		ecotype_id2nativename = get_ecotypeid2nativename(db_250k.metadata.bind)
		
		array_id2data = {}
		while rows.count()!=0:
			for row in rows:
				if row.array.maternal_ecotype_id in ecotype_id2qc_data:
					array_id = row.array_id
					ecotype_id = row.array.maternal_ecotype_id
					if array_id not in array_id2data:
						nativename = ecotype_id2nativename.get(ecotype_id)
						array_id2data[array_id] = PassingData(no_of_deletions=0, deletion_size=0, ecotype_id=ecotype_id,\
															nativename=nativename)
					array_id2data[array_id].no_of_deletions += 1
					array_id2data[array_id].deletion_size += (row.stop-row.start+1.0)
					
					
					real_counter += 1
				i += 1
				if i%5000==0:
					sys.stderr.write("%s%s\t%s"%('\x08'*80, i, real_counter))
			rows = query.offset(i).limit(block_size)
		sys.stderr.write("%s\t%s\t%s Done.\n"%('\x08'*80, i, real_counter))
		
		sys.stderr.write("Outputting deletion stats ...")
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = ['array_id', 'nativename', 'ecotype_id', 'tiling_no_of_deletions', 'tiling_deletion_size', 'qc_no_of_deletions',\
				'qc_deletion_size', 'qc_length_total']
		writer.writerow(header)
		for array_id, stat_data in array_id2data.iteritems():
			qc_stat_data = ecotype_id2data.get(stat_data.ecotype_id)
			row = [array_id, stat_data.nativename, stat_data.ecotype_id, stat_data.no_of_deletions, stat_data.deletion_size,\
				qc_stat_data.no_of_segments, qc_stat_data.length_affected, qc_stat_data.length]
			writer.writerow(row)
		del writer
		sys.stderr.write("Done.\n")
	
	"""
		
		# 2010-7-20
		qc_data_source_id=13
		cnv_type_id=1;
		cnv_method_id=8
		output_fname = os.path.expanduser('~/script/variation/data/CNV/CNVQC_source%s_vs_CNVCall_Method%s_stats.tsv'%\
										(qc_data_source_id, cnv_method_id))
		CNV.outputCNVQCVsCNVCallStats(db_250k, qc_data_source_id=qc_data_source_id, cnv_type_id=cnv_type_id, \
									cnv_method_id=cnv_method_id, output_fname=output_fname,)
		sys.exit(0)
	"""
	
	@classmethod
	def inspectCNVCallProbeDensityAndAmplitudeAndNumberOfProbes(cls, db_250k, array_id=None, output_dir=None, \
											minPercUnCoveredByLerContig=0.6, gridsize=10, cnv_method_id=6, \
											replaceAmpWithMedianIntensity=False, deletedFractionType=1,\
											colorBarLabel='percentage uncovered by Contig Sequence', \
											replaceNoOfProbesWithProbeDensity=False):
		"""
		2010-7-24
			add argument deletedFractionType.
				1: percUnCoveredByLerContig
				2: fractionDeletedInPECoverageData
			add argument colorBarLabel and replaceNoOfProbesWithProbeDensity
		2010-6-29
			add argument replaceAmpWithMedianIntensity
		2010-6-27
			
		"""
		sys.stderr.write("Inspecting relationship between probe density, amplitude, number of probes, for each cnv call... \n")
		i = 0
		block_size = 5000
		real_counter = 0
		import Stock_250kDB
		TableClass = Stock_250kDB.CNVCall
		query = TableClass.query.filter_by(array_id=array_id).filter_by(cnv_method_id=cnv_method_id)
		rows = query.offset(i).limit(block_size)
		session = db_250k.session
		
		no_of_probes_ls = []
		probeDensityLs = []
		amplitude_ls = []
		percUnCoveredByLerContig_ls = []
		ecotype_id = None
		import math
		from pymodule.DrawMatrix import Value2Color
		import ImageColor
		c_ls = []
		while rows.count()!=0:
			for row in rows:
				ecotype_id = row.array.maternal_ecotype_id
				if deletedFractionType==1:
					deletedFraction = row.percUnCoveredByLerContig
				else:
					deletedFraction = row.fractionDeletedInPECoverageData
					
				if deletedFraction is not None:
					#x_ls.append(row.amplitude)
					no_of_probes_ls.append(math.log10(row.no_of_probes_covered))
					probeDensity = row.no_of_probes_covered*1000.0/(row.stop-row.start+1.0)
					probeDensityLs.append(probeDensity)
					if deletedFraction>=minPercUnCoveredByLerContig:
						isDeletionTrue = 1
						real_counter += 1
					else:
						isDeletionTrue = 0
					if replaceAmpWithMedianIntensity:
						amp = row.median_intensity
					else:
						amp = row.amplitude
					amplitude_ls.append(amp)
					hslColor = Value2Color.value2HSLcolor(deletedFraction, min_value=0., max_value=1.,\
											treat_above_max_as_NA=False)
					fc = ImageColor.getrgb(hslColor)
					fc = [color_value/255. for color_value in fc]	#matplotlib accepts rgb in [0-1] range
					c_ls.append(fc)
					percUnCoveredByLerContig_ls.append(deletedFraction)
				
				i += 1
			if i%5000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, i, real_counter))
			rows = query.offset(i).limit(block_size)
		sys.stderr.write("%s%s\t%s\n"%('\x08'*80, i, real_counter))
		
		import pylab
		import matplotlib
		from matplotlib import cm
		import numpy
		if replaceAmpWithMedianIntensity:
			xlabel = 'median intensity'
			xName = 'MedianIntensity'
		else:
			xlabel = 'mean intensity'
			xName = 'MeanIntensity'
		if replaceNoOfProbesWithProbeDensity:
			ylabel = 'no of probes per kb'
			yName = 'ProbeDensity'
		else:
			ylabel = 'no of probes'
			yName = 'NoOfProbes'
		
		if len(no_of_probes_ls)>100:
			c_ls = numpy.array(c_ls)
			sys.stderr.write(" gridSize: %s "%gridsize)
			if output_dir and not os.path.isdir(output_dir):	#2010-5-5 test if output_dir is something
				os.makedirs(output_dir)
			common_fig_fname = 'array_%s_ecotype_%s_cnvMethod%s.png'%(array_id, ecotype_id, cnv_method_id)
			from pymodule.utils import addExtraLsToFilenamePrefix, addExtraToFilenamePrefix
			fig_fname = addExtraToFilenamePrefix(common_fig_fname, '%sAsX_%sAsY_minNotCoveredFraction%s_gridsize_%s'%\
												(xName, yName, minPercUnCoveredByLerContig, gridsize))
			"""
			fig_fname = addExtraToFilenamePrefix(common_fig_fname, 'AmplitudeAsX_NoOfProbesAsY_minUncoverPerc%s_gridsize_%s'%\
												(minPercUnCoveredByLerContig, gridsize))
			"""
			fig_fname = os.path.join(output_dir, fig_fname)
			title = "array %s e-id %s minFrac %s gridsize %s"%(array_id, ecotype_id, minPercUnCoveredByLerContig, gridsize)
			
			"""
			ylabel = 'No of probes per kb'
			xlabel = 'No of probes'
			max_phenotype = numpy.nanmax(percUnCoveredByLerContig_ls)
			min_phenotype = numpy.nanmin(percUnCoveredByLerContig_ls)
			phenotype_gap = max_phenotype - min_phenotype
			phenotype_jitter = phenotype_gap/10.
			norm = matplotlib.colors.Normalize(vmin=min_phenotype-phenotype_jitter, vmax=max_phenotype+phenotype_jitter)
			
			import pylab
			from mpl_toolkits.mplot3d import Axes3D
			fig = pylab.figure()
			ax = Axes3D(fig)
			ax.scatter(no_of_probes_ls, probeDensityLs, amplitude_ls,c=c_ls)
			pylab.title(title)
			pylab.xlabel(xlabel)
			pylab.ylabel(ylabel)
			pylab.show()
			"""
			if replaceNoOfProbesWithProbeDensity:
				y_ls = probeDensityLs
			else:
				y_ls = no_of_probes_ls
			cls.drawHexbin(amplitude_ls, y_ls, percUnCoveredByLerContig_ls, fig_fname=fig_fname, \
				gridsize=gridsize, title=title, \
						xlabel=xlabel, \
						ylabel=ylabel, colorBarLabel=colorBarLabel,\
						reduce_C_function=numpy.mean)
			sys.stderr.write("\n")
	
	"""
		#2010-6-28
		output_dir = os.path.expanduser('~/script/variation/data/CNV/figures/callCutoff/')
		for array_id in [3,4,41,150,151]:
			for gridsize in [20, 40, 60, 100]:
				CNV.inspectCNVCallProbeDensityAndAmplitudeAndNumberOfProbes(db_250k, array_id=array_id, output_dir=output_dir, \
											minPercUnCoveredByLerContig=0.6, gridsize=gridsize)
		sys.exit(0)
		
		# replace mean intensity of probes in one segment (amplitude) with median intensity instead 
		output_dir = os.path.expanduser('~/script/variation/data/CNV/figures/callCutoff/')
		cnv_method_id = 9
		for array_id in [3,4,41,150,151]:
			for gridsize in [20, 40, 60, 100]:
				CNV.inspectCNVCallProbeDensityAndAmplitudeAndNumberOfProbes(db_250k, array_id=array_id, output_dir=output_dir, \
											minPercUnCoveredByLerContig=0.6, gridsize=gridsize, cnv_method_id=cnv_method_id,\
											replaceAmpWithMedianIntensity=True)
		sys.exit(0)
		
		#2010-7-24 banyan.usc.edu
		output_dir = os.path.expanduser('~/script/variation/data/CNV/figures/CNVCallLabelVsIntensityByNoOfProbes/')
		cnv_method_id = 7
		rows = db_250k.metadata.bind.execute("select distinct array_id from cnv_call where cnv_method_id=7 and percUnCoveredByLerContig is null order by array_id")
		array_id_ls = []
		for row in rows:
			array_id_ls.append(row.array_id)
		for array_id in array_id_ls:
			for gridsize in [20, 40, 60, 100]:
				CNV.inspectCNVCallProbeDensityAndAmplitudeAndNumberOfProbes(db_250k, array_id=array_id, output_dir=output_dir, \
											minPercUnCoveredByLerContig=0.6, gridsize=gridsize, cnv_method_id=cnv_method_id,\
											replaceAmpWithMedianIntensity=False, colorBarLabel='FractionNotCoveredByPEData', \
											deletedFraction=2)
		sys.exit(0)
		
		#2010-7-24 banyan.usc.edu
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVCallLabelVsIntensityByNoOfProbes/')
		cnv_method_id = 10
		rows = db_250k.metadata.bind.execute("select distinct array_id from cnv_call where cnv_method_id=%s and \
			percUnCoveredByLerContig is null order by array_id"%cnv_method_id)
		array_id_ls = []
		for row in rows:
			array_id_ls.append(row.array_id)
		print "%s arrays"%len(array_id_ls)
		for array_id in array_id_ls:
			for gridsize in [40, 60,]:
				CNV.inspectCNVCallProbeDensityAndAmplitudeAndNumberOfProbes(db_250k, array_id=array_id, output_dir=output_dir, \
									minPercUnCoveredByLerContig=0.6, gridsize=gridsize, cnv_method_id=cnv_method_id,\
									deletedFractionType=2, colorBarLabel='FractionNotCoveredByPEData',\
									replaceNoOfProbesWithProbeDensity=True)
		sys.exit(0)
		
		
	"""
	@classmethod
	def addAmplitudeFunctor(cls, param_obj, cnv_segment_obj, cnv_qc_call=None):
		"""
		2009-12-9
			adjust argument order and process only if cnv_qc_call is not None
		2009-11-4
		"""
		if not hasattr(param_obj, 'amp_ls'):
			setattr(param_obj, 'amp_ls', [])
		if cnv_qc_call is not None:
			qc_chromosome, qc_start, qc_stop = cnv_qc_call[:3]
			cnv_qc_call_id = cnv_qc_call[-1]
			param_obj.cnv_qc_call_id_set.add(cnv_qc_call_id)
			param_obj.amp_ls.append(cnv_segment_obj.amplitude)
	
	@classmethod
	def drawHistOfAmpOfValidatedDeletions(cls, db_250k, input_fname_ls, output_fname_prefix, data_source_id=1, cnv_type_id=1, \
								min_QC_segment_size=200, min_no_of_probes=5, max_boundary_diff=10000, \
								max_diff_perc=0.10, count_embedded_segment_as_match=True, min_reciprocal_overlap=0.6):
		"""
		2009-11-4
			draw histogram of amplitude of segments who are validated according to certain QC data
		"""
		cnvQCData = cls.getCNVQCDataFromDB(data_source_id, cnv_type_id=cnv_type_id, \
															min_QC_segment_size=min_QC_segment_size,\
															min_no_of_probes=min_no_of_probes, \
															min_reciprocal_overlap=min_reciprocal_overlap)
		ecotype_id2cnv_qc_call_data = cnvQCData.ecotype_id2cnv_qc_call_data
		from pymodule import PassingData
		param_obj = PassingData(amp_ls=[], cnv_qc_call_id_set=set())
		cls.compareCNVSegmentsAgainstQCHandler(input_fname_ls, ecotype_id2cnv_qc_call_data, cls.addAmplitudeFunctor, param_obj, \
											deletion_cutoff=None, max_boundary_diff=max_boundary_diff, max_diff_perc=max_diff_perc,\
											min_no_of_probes=min_no_of_probes, \
											count_embedded_segment_as_match=count_embedded_segment_as_match)

		sys.stderr.write("data_source_id %s, min_QC_segment_size %s, min_no_of_probes: %s, max_boundary_diff: %s, max_diff_perc: %s, count_embedded_segment_as_match: %s.\n"%\
						(data_source_id, min_QC_segment_size, min_no_of_probes, max_boundary_diff, max_diff_perc,\
						count_embedded_segment_as_match))
		param_obj.ecotype_id2cnv_qc_call_data = ecotype_id2cnv_qc_call_data
		cls.outputFalsePositiveRate(param_obj)
		cls.outputFalseNegativeRate(param_obj)
				
		sys.stderr.write("Drawing ...")
		import pylab
		pylab.clf()
		pylab.title("Histogram of amplitude of known deletions from source %s"%data_source_id)
		pylab.hist(param_obj.amp_ls, 20)
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		sys.stderr.write("Done.\n")
	
	"""
	# 2009-11-4
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s.GADA_A0.5T4M5.tsv'%i))
	#input_fname_ls = [os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr1.GADA_A0.5T4M5.tsv')]
	min_QC_segment_size = 5
	min_no_of_probes = 5
	count_embedded_segment_as_match = True
	data_source_id = 1
	for max_boundary_diff in [10000]:
		for max_diff_perc in [0.20, 0.3]:
			output_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/figures/Clark2007aDeletionsAmplitudeHistCall48_mxdiff_%s_mxperc_%s'%(max_boundary_diff, max_diff_perc))	
			CNV.drawHistOfAmpOfValidatedDeletions(db_250k, input_fname_ls, output_fname_prefix, data_source_id=data_source_id, \
									cnv_type_id=1, \
								min_QC_segment_size=min_QC_segment_size, \
								max_boundary_diff=max_boundary_diff, max_diff_perc=max_diff_perc, \
								min_no_of_probes=min_no_of_probes,\
								count_embedded_segment_as_match=count_embedded_segment_as_match)
	
	
	"""
	
	@classmethod
	def drawDeletedCNVProbeHistogram(cls, db_250k, input_fname_ls, CNV_qc_fname, output_fname_prefix, min_no_of_bp_deleted=20):
		"""
		2009-10-12
			input_fname_ls: a list of GADA output files, containing segmentation results and corresponding amplitudes
			for the probes that are deleted according to CNV_qc_fname, look at the histogram of their amplitudes outputted from GADA
		"""
		sys.stderr.write("Drawing deleted CNV probe intensity histogram ... \n")
		import fileinput
		
		from pymodule import SNPData
		CNVQCData = SNPData(input_fname=CNV_qc_fname, turn_into_array=1, ignore_2nd_column=1, data_starting_col=2, \
						turn_into_integer=False)
		new_col_id2col_index = {}
		new_col_id_ls = []
		for col_id in CNVQCData.col_id_ls:
			new_col_id = int(col_id.split('_')[0])	# it's ecotypeid, as in "6899_Bay-0".
			new_col_id_ls.append(new_col_id)
			new_col_id2col_index[new_col_id] = CNVQCData.col_id2col_index[col_id]
		CNVQCData.col_id_ls = new_col_id_ls
		CNVQCData.col_id2col_index = new_col_id2col_index
		
		from Stock_250kDB import Probes, ArrayInfo
		new_row_id2row_index = {}
		new_row_id_ls = []
		for row_id in CNVQCData.row_id_ls:
			probe_id = int(row_id)
			probe = Probes.get(probe_id)
			new_row_id = (probe.chromosome, probe.position)
			new_row_id_ls.append(new_row_id)
			new_row_id2row_index[new_row_id] = CNVQCData.row_id2row_index[row_id]
		CNVQCData.row_id_ls = new_row_id_ls
		CNVQCData.row_id2row_index = new_row_id2row_index
		
		sys.stderr.write("Getting probe amplitude from %s ... \n"%repr(input_fname_ls))
		amp_ls = []
		array_id2array = {}
		counter = 0
		real_counter = 0
		for line in fileinput.input(input_fname_ls):
			if line.find("array_id")==0:
				continue
			line = line.strip()
			row = line.split('\t')
			array_id = int(row[0])
			if array_id not in array_id2array:
				array = ArrayInfo.get(array_id)
				array_id2array[array_id] = array
			array = array_id2array[array_id]
			counter += 1
			if array.maternal_ecotype_id in CNVQCData.col_id2col_index:	# array is in CNVQCData
				start_probe = row[1].split('_')	# split chr_pos
				start_probe = map(int, start_probe)
				stop_probe = row[2].split('_')
				stop_probe = map(int, stop_probe)
				amplitude = float(row[4])
				col_index = CNVQCData.col_id2col_index[array.maternal_ecotype_id]
				for i in range(len(CNVQCData.row_id_ls)):
					row_id = CNVQCData.row_id_ls[i]
					chr, pos = row_id
					no_of_bp_deleted = CNVQCData.data_matrix[i][col_index]
					if no_of_bp_deleted!='?':	# not NA
						no_of_bp_deleted = int(no_of_bp_deleted)
						if no_of_bp_deleted>=min_no_of_bp_deleted and chr==start_probe[0] and pos>=start_probe[1] and pos<=stop_probe[1]:
							amp_ls.append(amplitude)
							real_counter += 1
			if counter%10000==0:
				sys.stderr.write('%s%s\t%s'%('\x08'*80, counter, real_counter))
		sys.stderr.write("Done.\n")
		
		import pylab
		pylab.clf()
		pylab.title("Histogram of intensity of probes with >=20 bases deleted")
		pylab.hist(amp_ls, 20)
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
	
	""""
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_43_CNV_norm_intensity_chr%s.GADAJRN_A0.5T4M5.tsv'%i))
	
	CNV_qc_fname = '/Network/Data/250K/tmp-yh/CNV/dazhe_CNV_probes_w_del.csv'
	output_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/figures/ProbeFromPerlegenCNVQC_AmplitudeHistCall43')
	CNV.drawDeletedCNVProbeHistogram(db_250k, input_fname_ls, CNV_qc_fname, output_fname_prefix, min_no_of_bp_deleted=20)
	
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_43_CNV_norm_intensity_chr%s.GADAJRN_A0.5T4M5.tsv'%i))
		CNV_qc_fname = '/Network/Data/250K/tmp-yh/CNV/dazhe_CNV_probes_w_del.csv'
		output_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/figures/ProbeFromPerlegenCNVQC_AmplitudeCall43Chr%sHist'%i)
		CNV.drawDeletedCNVProbeHistogram(db_250k, input_fname_ls, CNV_qc_fname, output_fname_prefix, min_no_of_bp_deleted=20)
	
	input_fname_ls = []
	for i in [2,4]:
		input_fname_ls = [os.path.expanduser('~/mnt2/panfs/250k/CNV/call_48_CNV_QNormalize_chr%s.GADAJRN_A0.5T4M10.tsv'%i)]
		CNV_qc_fname = '/Network/Data/250K/tmp-yh/CNV/dazhe_CNV_probes_w_del.csv'
		output_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/figures/ProbeFromPerlegenCNVQC_AmplitudeCall48Chr%sHist'%i)
		CNV.drawDeletedCNVProbeHistogram(db_250k, input_fname_ls, CNV_qc_fname, output_fname_prefix, min_no_of_bp_deleted=20)

	
	"""
	
	
	@classmethod
	def drawDeletedCNVDeletionSizeHist(cls, db_250k, input_fname_ls, output_fname_prefix, max_amplitude=-0.33):
		"""
		2009-10-12
			
		"""
		sys.stderr.write("Drawing deleted CNV probe intensity histogram ... \n")
		import fileinput
		
		from Stock_250kDB import Probes, ArrayInfo
		from pymodule import getColName2IndexFromHeader
		sys.stderr.write("Getting probe amplitude from %s ... \n"%repr(input_fname_ls))
		array_id2array = {}
		length_ls = []
		
		counter = 0
		real_counter = 0
		array_id2length_ls = {}
		input_handler = fileinput.input(input_fname_ls)
		header = input_handler.readline().strip().split('\t')
		col_name2index = getColName2IndexFromHeader(header)
		
		for line in input_handler:
			if line.find("array_id")!=-1:
				continue
			line = line.strip()
			row = line.split('\t')
			
			ecotype_id_idx = col_name2index.get('ecotype_id', col_name2index.get('array_id'))
			#cnv_ecotype_id = int(row[ecotype_id_idx])
			array_id = int(row[col_name2index.get('array_id')])
			if array_id not in array_id2array:
				array = ArrayInfo.get(array_id)
				array_id2array[array_id] = array
				array_id2length_ls[array_id] = []
			array = array_id2array[array_id]
			counter += 1			
			amplitude = float(row[col_name2index['amplitude']])
			if amplitude<=max_amplitude:
				start_probe = row[col_name2index['start_probe']].split('_')	# split chr_pos
				start_probe = map(int, start_probe)
				stop_probe = row[col_name2index['end_probe']].split('_')
				stop_probe = map(int, stop_probe)
				segment_chromosome = start_probe[0]
				segment_start_pos = start_probe[1]-12
				segment_stop_pos = stop_probe[1]+12
				segment_length = abs(segment_stop_pos-segment_start_pos+1)
				
				#length = stop_probe[1]-start_probe[1]+1
				
				length_ls.append(segment_length)
				array_id2length_ls[array_id].append(segment_length)
				real_counter += 1
			if counter%10000==0:
				sys.stderr.write('%s%s\t%s'%('\x08'*80, counter, real_counter))
		sys.stderr.write("Done.\n")
		
		import pylab
		for array_id, length_ls in array_id2length_ls.iteritems():
			if len(length_ls)>20:
				array = array_id2array[array_id]
				output_fname_prefix_array = '%s_array_%s_ecotype_%s'%(output_fname_prefix, array_id, array.maternal_ecotype_id)
				pylab.clf()
				pylab.title("Histogram of length of deletions with max_amplitude=%s"%max_amplitude)
				pylab.hist(length_ls, 40)
				pylab.savefig('%s.png'%output_fname_prefix_array, dpi=300)
				pylab.savefig('%s.svg'%output_fname_prefix_array, dpi=300)
	
	
	"""
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_43_CNV_norm_intensity_chr%s.GADAJRN_A0.5T4M5.tsv'%i))
	max_amplitude = -0.33
	output_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/figures/call_43_deletion_length_histogram_max_amp_%s'%max_amplitude)
	CNV.drawDeletedCNVDeletionSizeHist(db_250k, input_fname_ls, output_fname_prefix, max_amplitude=max_amplitude)
	
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_43_CNV_norm_intensity_chr%s.GADAJRN_A0.5T4M5.tsv'%i))
	
	max_amplitude = -0.33
	output_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/figures/call_43_deletion_length_histogram_max_amp_%s'%max_amplitude)
	CNV.drawDeletedCNVDeletionSizeHist(db_250k, input_fname_ls, output_fname_prefix, max_amplitude=max_amplitude)
	
	"""
	
	
	@classmethod
	def drawGADAOutputAmpHist(cls, db_250k,input_fname, output_fname_prefix, max_amplitude=None,\
											drawKernelDensity=False, gridsize=40):
		"""
		#2010-7-1
			if argument drawKernelDensity == 2:
				draw a 2D histogram, density by color
			add argument gridsize 40 for drawKernelDensity == 2
		2010-6-7
			It draws a histogram of "copy-number" (amplitude) of segments outputted by discoverLerDeletionDuplication (using GADA).
			
			Purpose of this function is to see whether there's a natural cutoff lying in a bi-modal distribution.
			
			argument max_amplitude is used to filter out segments whose amplitude is above this.
				purpose is to avoid saturation of histogram by normal fragments.
		"""
		import csv, math
		from pymodule import figureOutDelimiter, getColName2IndexFromHeader
		reader = csv.reader(open(input_fname), delimiter=figureOutDelimiter(input_fname))
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header)
		copy_number_col_index = col_name2index.get('amplitude')
		ecotype_id_col_index = col_name2index.get('ecotype_id')
		array_id_col_index = col_name2index.get('array_id')
		no_of_probes_col_index = col_name2index.get('length')
		
		array_id2ecotype_id_amplitude_ls = {}
		for row in reader:
			ecotype_id = int(row[ecotype_id_col_index])
			array_id = row[array_id_col_index]
			if array_id not in array_id2ecotype_id_amplitude_ls:
				array_id2ecotype_id_amplitude_ls[array_id] = [ecotype_id, [], []]
			amplitude = float(row[copy_number_col_index])
			if max_amplitude is not None and amplitude>max_amplitude:
				continue
			no_of_probes = int(row[no_of_probes_col_index])
			array_id2ecotype_id_amplitude_ls[array_id][1].append(amplitude)
			array_id2ecotype_id_amplitude_ls[array_id][2].append(math.log10(no_of_probes))
		
		from common import get_ecotypeid2nativename
		ecotype_id2nativename = get_ecotypeid2nativename(db_250k.metadata.bind)
		import Stock_250kDB
		
		import pylab
		for array_id, ecotype_id_amplitude_ls in array_id2ecotype_id_amplitude_ls.iteritems():
			pylab.clf()
			ecotype_id, amplitude_ls, no_of_probes_ls =  ecotype_id_amplitude_ls
			nativename = ecotype_id2nativename.get(ecotype_id)
			array = Stock_250kDB.ArrayInfo.get(array_id)
			
			title = 'Array %s, e %s, %s median %s'%(array_id, ecotype_id, nativename, array.median_intensity)
			pylab.title(title)
			xlabel = 'segment mean intensity'
			ylabel = 'density'
			if drawKernelDensity==True or drawKernelDensity==1:
				import statistics	# 2010-5-30 package from Michiel De Hoon
				y, x = statistics.pdf(amplitude_ls)
				pylab.plot(x, y, alpha=0.7)
			elif drawKernelDensity==2:	#2010-7-1	draw a 2D histogram, density by color
				import matplotlib.cm as cm
				pylab.hexbin(amplitude_ls, no_of_probes_ls, C=[1]*len(amplitude_ls), \
							gridsize=gridsize, reduce_C_function=CNV.logSum, cmap=cm.jet)
				ylabel = 'log10(no of probes)'
				cb = pylab.colorbar()
				cb.set_label("log10(count)")
			else:
				pylab.hist(amplitude_ls, 35)
			pylab.xlabel(xlabel)
			pylab.ylabel(ylabel)
			pylab.xlim([-1.5, 1])
			pylab.savefig('%s_e%s_a%s.png'%(output_fname_prefix, ecotype_id, array_id), dpi=300)
			sys.stderr.write("Done.\n")
	
	"""

		max_amplitude = 3	# 0.6
		drawKernelDensity = True
		common_prefix = '~/script/variation/data/CNV/call_53_WABQN_b100_j50_lts_98_arrays_GADA2_A0.5T4.0M5'
		input_fname = os.path.expanduser('%s.tsv'%common_prefix)
		if drawKernelDensity is False:
			common_prefix += '_hist'
		if max_amplitude is not None:
			common_prefix += '_maxAmp%.1f'%max_amplitude
		output_fname_prefix = os.path.expanduser('%s'%(common_prefix))
		CNV.drawGADAOutputAmpHist(db_250k, input_fname, output_fname_prefix, \
											max_amplitude=max_amplitude, drawKernelDensity=drawKernelDensity)
		sys.exit(0)


		#2010-7-1
		max_amplitude = 3	# 0.6
		drawKernelDensity = True
		common_prefix = 'call_53_WABQN_b100_j50_lts_2Flanking_GADA_A0.5T12.0M5'
		input_dir = os.path.expanduser('~/mnt/hpc-cmb/script/variation/data/CNV/')
		output_dir = os.path.expanduser('~/script/variation/data/CNV/GADAOutputAmpHist/')
		input_fname = os.path.join(input_dir, '%s.tsv'%common_prefix)
		if drawKernelDensity is False:
			common_prefix += '_hist'
		if max_amplitude is not None:
			common_prefix += '_maxAmp%.1f'%max_amplitude
		output_fname_prefix = os.path.join(output_dir, common_prefix)
		CNV.drawGADAOutputAmpHist(db_250k, input_fname, output_fname_prefix, \
											max_amplitude=max_amplitude, drawKernelDensity=drawKernelDensity)
		sys.exit(0)
		
		#2010-7-1
		max_amplitude = 3	# 0.6
		drawKernelDensity = 2
		common_prefix = 'call_53_WABQN_b100_j50_lts_2Flanking_GADA_A0.5T12.0M5'
		input_dir = os.path.expanduser('~/mnt/hpc-cmb/script/variation/data/CNV/')
		output_dir = os.path.expanduser('~/script/variation/data/CNV/GADAOutputAmpHist/')
		input_fname = os.path.join(input_dir, '%s.tsv'%common_prefix)
		if drawKernelDensity is False:
			common_prefix += '_hist'
		if max_amplitude is not None:
			common_prefix += '_maxAmp%.1f'%max_amplitude
		output_fname_prefix = os.path.join(output_dir, common_prefix)
		CNV.drawGADAOutputAmpHist(db_250k, input_fname, output_fname_prefix, \
								max_amplitude=max_amplitude, drawKernelDensity=drawKernelDensity, gridsize=40)
		sys.exit(0)
	"""
	@classmethod
	def traverseCNVCall(cls, db_250k, cnv_method_id=8, cnv_type_id=1, ecotype_id=None, ecotype_id_ls=None,\
					min_no_of_probes=None, min_segment_size=None, processClassIns=None, \
					order_by_string='order by array_id, chromosome, start, stop', param_obj=None):
		"""
		2010-7-29
			bugfix.
				start to use cnv_method_id to restrict query.
		2010-7-21
			a common function which traverses through all CNVCalls.
		"""
		from pymodule import PassingData
		import Stock_250kDB
		sql_string = "select a.maternal_ecotype_id as ecotype_id, c.array_id, c.chromosome, c.start, c.stop, c.size_affected, \
					c.no_of_probes_covered, c.amplitude, c.probability, c.percUnCoveredByLerContig, c.cnv_type_id, \
					c.cnv_method_id, c.comment from %s c, %s a \
					where c.array_id=a.id "%\
					(Stock_250kDB.CNVCall.table.name, Stock_250kDB.ArrayInfo.table.name)
					
		if cnv_type_id is not None:
			sql_string += " and c.cnv_type_id=%s"%cnv_type_id
		if cnv_method_id is not None:	#2010-7-29 bugfix
			sql_string += " and c.cnv_method_id=%s"%cnv_method_id
		
		if ecotype_id is not None:
			sql_string += " and a.maternal_ecotype_id=%s"%ecotype_id
		elif ecotype_id_ls:
			str_ecotype_id_ls = map(str, ecotype_id_ls)
			sql_string += ' and a.maternal_ecotype_id in (%s)'%(','.join(str_ecotype_id_ls))
		if min_no_of_probes is not None:
			sql_string += " and c.no_of_probes_covered>=%s"%min_no_of_probes
		if min_segment_size is not None:
			sql_string += " and (c.stop-c.start+1)>=%s"%min_segment_size
		
		#2010-7-30 offer custom order_by string
		sql_string += ' %s'%order_by_string
		
		#sql_string += " order by RAND()"	# 2010-1-26 random ordering to optimize the binary_tree, not needed for RBDict.
		
		rows = db_250k.metadata.bind.execute(sql_string)
		count = 0
		for row in rows:	#segments are sorted in chromosomal order
			count += 1
			processClassIns.run(row, param_obj)
			if count%5000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, count, processClassIns.real_counter))
		sys.stderr.write("%s%s\t%s Done.\n"%('\x08'*80, count, processClassIns.real_counter))
	
	class ProcessAdjacentCNVCallGapHistogram(object):
		"""
		2010-7-22
			used by drawCNVCallSizeHistogram which plugs it into traverseCNVCall()
			This class gets size of gaps between adjacent deletions on each chromosome for each array. 
		"""
		def __init__(self, ecotype_id2nativename=None, **keywords):
			"""
			2010-7-29
				add argument **keywords
			"""
			self.array_id2data = {}
			self.ecotype_id2nativename = ecotype_id2nativename
			self.real_counter = 0
			
			self.previous_row = None
			self.previous_array_id = None
			
			self.dataType = 1 # 2010-7-29
			for keyword, value in keywords.iteritems():
				setattr(self, keyword, value)
			
		def run(self, row, param_obj=None):
			"""
			This function assumes row comes in 
			"""
			from pymodule import PassingData
			import math
			array_id = row.array_id
			ecotype_id = row.ecotype_id
			if array_id not in self.array_id2data:
				nativename = self.ecotype_id2nativename.get(ecotype_id)
				self.array_id2data[array_id] = PassingData(x_ls=[], y_ls=[], ecotype_id=ecotype_id,\
													nativename=nativename)
			if self.previous_row is not None:
				if self.previous_array_id == row.array_id and self.previous_row.chromosome==row.chromosome:
					try:
						# abs() is used here because there is one case two deletions are overlapping.
						gap_length = abs(row.start - self.previous_row.stop + 1)
						gap_ratio1 = math.log10(float(gap_length)/row.size_affected)
						gap_ratio2 = math.log10(float(gap_length)/self.previous_row.size_affected)
						if self.dataType==1:
							data_x = gap_ratio1
							data_y = gap_ratio2
						elif self.dataType ==2:
							data_x = min(gap_ratio1, gap_ratio2)
							data_y = math.log10(gap_length)
					except:
						"""
						# two deletions below overlap with each other because the central positions of two probes are only 3-base apart. 
						#array_id: 65
						previous_row: (6918L, 65L, 2L, 5939417L, 5942609L, 3193L, 18L, -0.38029299999999999, 0.75223399999999996, None, 1L, 8L, None)
						row: (6918L, 65L, 2L, 5942588L, 5943887L, 1300L, 9L, -0.563971, 0.75422900000000004, None, 1L, 8L, None)
						"""
						print "array_id:", array_id
						print "previous_row:", self.previous_row
						print "row:", row 
						sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
						import traceback
						traceback.print_exc()
					self.array_id2data[array_id].x_ls.append(data_x)
					self.array_id2data[array_id].y_ls.append(data_y)
					#self.array_id2data[array_id].y_ls.append(abs(row.amplitude-self.previous_row.amplitude))
					self.real_counter += 1
			
			self.previous_row = row
			self.previous_array_id = row.array_id
	
	@classmethod
	def drawCNVCallGapHistogram(cls, db_250k, cnv_method_id=8, cnv_type_id=1, ecotype_id=None, qc_data_source_id=13,\
								output_dir=None, xlabel='log10(gap/segment2)', \
								ylabel_in_2D='log10(gap/segment1)', drawType=1, fileNamePrefix='CNVCallGapHist', dataType=1):
		"""
		2010-7-29
			add argument dataType
				1: gap_ratio1 in x_ls, gap_ratio2 in y_ls
				2: min(gap-ratio1, gap-ratio2) in x_ls, gap-length in y_ls
		2010-7-20
			drawType
				1: 1D histogram of gap ratio 1
				2: 2D histogram of both gap ratio 1 and gap ratio 2
		"""
		# for the gap ratio
		cls.drawCNVCallSizeHistogram(db_250k, cnv_method_id=cnv_method_id, cnv_type_id=cnv_type_id, \
								ecotype_id=ecotype_id, qc_data_source_id=qc_data_source_id,\
								output_dir=output_dir, drawType=drawType, xlabel=xlabel, \
								ylabel_in_2D=ylabel_in_2D, processClass=cls.ProcessAdjacentCNVCallGapHistogram,\
								xlim_in_1D=None,fileNamePrefix=fileNamePrefix, dataType=dataType)
		"""
		# for the gap size. need to change cls.ProcessAdjacentCNVCallGapHistogram.run()
		cls.drawCNVCallSizeHistogram(db_250k, cnv_method_id=cnv_method_id, cnv_type_id=cnv_type_id, \
								ecotype_id=ecotype_id, qc_data_source_id=qc_data_source_id,\
								output_dir=output_dir, drawType=drawType, xlabel='log10(gap size)', \
								ylabel_in_2D='amplitude delta', processClass=cls.ProcessAdjacentCNVCallGapHistogram,\
								xlim_in_1D=None)
		"""
	"""
		# 2010-7-20
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVCallGapSizeHist/')
		CNV.drawCNVCallGapHistogram(db_250k, cnv_method_id=8, cnv_type_id=1, ecotype_id=None, qc_data_source_id=13, \
								output_dir=output_dir, drawType=1, fileNamePrefix='CNVCallGapRatioHist', dataType=1)
		sys.exit(0)
		
		# 2010-7-29
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVCallGapRatioHist/')
		CNV.drawCNVCallGapHistogram(db_250k, cnv_method_id=11, cnv_type_id=1, ecotype_id=None, qc_data_source_id=13, \
								output_dir=output_dir, xlabel='min(gap_ratio1, gap_ratio2)', \
								ylabel_in_2D='gap_length', drawType=1, fileNamePrefix='CNVCallGapRatioVsGapLenHist', dataType=2)
		sys.exit(0)
		
		# 2010-7-20
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVCallGapRatioHist/')
		CNV.drawCNVCallGapHistogram(db_250k, cnv_method_id=11, cnv_type_id=1, ecotype_id=None, qc_data_source_id=13, \
								output_dir=output_dir, drawType=2, fileNamePrefix='CNVCallGapRatioHist', dataType=1)
		sys.exit(0)
		
		# 2010-7-20
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVCall_GapRatioHist/')
		CNV.drawCNVCallGapHistogram(db_250k, cnv_method_id=8, cnv_type_id=1, ecotype_id=None, qc_data_source_id=13, \
								output_dir=output_dir, drawType=1)
		sys.exit(0)
		
		
		# 2010-7-29
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVCallGapRatioHist/')
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		CNV.drawCNVCallGapHistogram(db_250k, cnv_method_id=17, cnv_type_id=1, ecotype_id=None, qc_data_source_id=13, \
								output_dir=output_dir, xlabel='min(gap_ratio1, gap_ratio2)', \
								ylabel_in_2D='gap_length', fileNamePrefix='CNVCallGapRatioVsGapLenHist', dataType=2)
		sys.exit(0)
		
	"""
	
	class ProcessAcrossArrayCNVCallOverlapHistogram(object):
		"""
		2010-7-22
			used by drawAcrossArrayCNVCallOverlapRatioHistogram() which plugs it into traverseCNVCall()
			This class gets size of gaps between adjacent deletions on each chromosome for each array. 
		"""
		def __init__(self, ecotype_id2nativename=None, **keywords):
			"""
			2010-7-29
				add argument **keywords
			"""
			self.array_id2data = {}
			self.ecotype_id2nativename = ecotype_id2nativename
			self.real_counter = 0
			
			self.previous_row = None
			self.previous_array_id = None
			
			self.dataType = 1 # 2010-7-29
			self.logRatio = False	#2010-7-30
			for keyword, value in keywords.iteritems():
				setattr(self, keyword, value)
			
			self.aheadSegmentLs = []
			self.x_ls = []
			self.y_ls = []
			
		def run(self, row, param_obj=None):
			"""
			This function assumes row comes in chromosomal order. all arrays mixed up.
			"""
			from pymodule import PassingData
			from pymodule.CNV import get_overlap_ratio
			import math
			array_id = row.array_id
			ecotype_id = row.ecotype_id
			if array_id not in self.array_id2data:
				nativename = self.ecotype_id2nativename.get(ecotype_id)
				self.array_id2data[array_id] = PassingData(x_ls=[], y_ls=[], ecotype_id=ecotype_id,\
													nativename=nativename)
			
			segment = [row.chromosome, row.start, row.stop, row.no_of_probes_covered, \
					row.size_affected,\
					row.amplitude, row.probability, row.array_id]
			self.aheadSegmentLs.append(segment)
			
			segmentOverlapIndex = 0	#this is the index in aheadSegmentLs for the new aheadSegmentLs to start
			no_of_segments = len(self.aheadSegmentLs)
			if no_of_segments>1:	#more than 1 segments
				#calcualte the overlap ratio between the current segments and all previous potentially overlapping segments
				for i in range(no_of_segments-1):
					aheadSegment = self.aheadSegmentLs[i]
					if aheadSegment[0]!=segment[0]:		#on different chromosomes, ignore
						segmentOverlapIndex = i+1
						continue
					overlap1, overlap2, overlap_length = get_overlap_ratio([aheadSegment[1], aheadSegment[2]], \
																[segment[1], segment[2]])[:3]
					if overlap1==0 or overlap2==0:
						segmentOverlapIndex = i+1
					elif overlap1<0 or overlap2<0:
						sys.stderr.write("aheadSegment %s.\n"%repr(aheadSegment))
						sys.stderr.write("segment %s.\n"%repr(segment))
						sys.stderr.write("Error: overlap1 %s or overlap2 %s is zero.\n"%(overlap1, overlap2)) 
					else:
						if self.logRatio:
							overlap1 = math.log10(overlap1)
							overlap2 = math.log10(overlap2)
						self.x_ls.append(overlap1)
						if self.dataType==1:
							self.y_ls.append(overlap2)
						elif self.dataType==2:
							self.y_ls.append(math.log10(overlap_length))
						self.real_counter += 1
				self.aheadSegmentLs = self.aheadSegmentLs[segmentOverlapIndex:]	#remove segments that are not overlapping with the current.
	
	@classmethod
	def drawAcrossArrayCNVCallOverlapRatioHistogram(cls, db_250k, cnv_method_id=8, cnv_type_id=1, ecotype_id=None, qc_data_source_id=13,\
								output_dir=None, drawType=1, xlabel='log10(overlap ratio1)', ylabel_in_2D='log10(overlap ratio2)',\
								processClass=None, xlim_in_1D=None, fileNamePrefix='AcrossArrayCNVCallOverlapHist', **keywords):
		"""
		2010-7-30
			dataType (hidden argument in **keywords)
				1: overlap ratio1 and overlap ratio2
				2: overlap ratio1 and overlap_length
			argument drawType is useless.
		"""
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		if qc_data_source_id is not None and ecotype_id is None:
			ecotype_id_ls = cls.getEcotypeIDListInCNVQCAccessionGivenDataSourceID(qc_data_source_id)
		else:
			ecotype_id_ls = None
		
		from pymodule import PassingData
		from common import get_ecotypeid2nativename
		ecotype_id2nativename = get_ecotypeid2nativename(db_250k.metadata.bind)
		import math
		if processClass is None:
			processClass = cls.ProcessAcrossArrayCNVCallOverlapHistogram
		processor = processClass(ecotype_id2nativename=ecotype_id2nativename, **keywords)
		
		order_by_string='order by chromosome, start, stop'
		param_obj = PassingData(order_by_string=order_by_string)
		
		cls.traverseCNVCall(db_250k, cnv_method_id=cnv_method_id, cnv_type_id=cnv_type_id, \
					ecotype_id=ecotype_id, ecotype_id_ls=ecotype_id_ls,\
					min_no_of_probes=None, min_segment_size=None, processClassIns=processor, \
					order_by_string=order_by_string, param_obj=param_obj)
		
		from pymodule.utils import addExtraLsToFilenamePrefix
		
		sys.stderr.write("drawing ...")
		x_ls = processor.x_ls
		y_ls = processor.y_ls
		
		import pylab
		pylab.clf()
		output_fname = os.path.join(output_dir, '%s_CNVMethod%s_CNVType%s_dataType%s_logRatio%s.png'%(fileNamePrefix, \
									cnv_method_id,	cnv_type_id, getattr(processor, 'dataType', 1), \
									getattr(processor, 'logRatio', False) ))
		title = '%s type %s objects.'%(len(x_ls), cnv_type_id)
		
		pylab.title(title)
		"""
		import statistics	# 2010-5-30 package from Michiel De Hoon
		y, x = statistics.pdf(x_ls)
		pylab.loglog(x, y, alpha=0.7)
		pylab.grid(True, alpha=0.6)
		"""
		pylab.hist(x_ls, 20, log=True)
		if xlim_in_1D:
			pylab.xlim(xlim_in_1D)
		
		pylab.xlabel(xlabel)
		pylab.ylabel('Count')
		#
		pylab.savefig(output_fname, dpi=300)
		
		pylab.clf()
		C_ls = [1]*len(y_ls)
		fig_fname = addExtraLsToFilenamePrefix(output_fname, ['2D'])
		cls.drawHexbin(x_ls, y_ls, C_ls, fig_fname=fig_fname, gridsize=20, \
				title=title, \
				xlabel=xlabel, \
				ylabel=ylabel_in_2D, \
				colorBarLabel='log10(count)', reduce_C_function=CNV.logSum)
		sys.stderr.write("Done.\n")
	
	"""
		# 2010-7-30
		output_dir = os.path.expanduser('~/script/variation/data/CNV/AcrossArrayCNVCallOverlapRatioHist/')
		CNV.drawAcrossArrayCNVCallOverlapRatioHistogram(db_250k, cnv_method_id=13, cnv_type_id=1, ecotype_id=None, \
				qc_data_source_id=13, \
				output_dir=output_dir, drawType=2, fileNamePrefix='AcrossArrayCNVCallOverlapRatioHist', dataType=1)
		sys.exit(0)
		
		# 2010-7-30
		output_dir = os.path.expanduser('~/script/variation/data/CNV/AcrossArrayCNVCallOverlapRatioHist/')
		CNV.drawAcrossArrayCNVCallOverlapRatioHistogram(db_250k, cnv_method_id=13, cnv_type_id=1, ecotype_id=None, \
				qc_data_source_id=13, xlabel='overlap ratio1', ylabel_in_2D='log10(overlap_length)',\
				output_dir=output_dir, drawType=2, fileNamePrefix='AcrossArrayCNVCallOverlapRatioHist', \
				dataType=2, logRatio=True)
		sys.exit(0)
		# 2010-7-30
		output_dir = os.path.expanduser('~/script/variation/data/CNV/AcrossArrayCNVCallOverlapRatioHist/')
		CNV.drawAcrossArrayCNVCallOverlapRatioHistogram(db_250k, cnv_method_id=16, cnv_type_id=1, ecotype_id=None, \
				qc_data_source_id=13, \
				output_dir=output_dir, drawType=2, fileNamePrefix='AcrossArrayCNVCallOverlapRatioHist', \
				dataType=1, logRatio=False)
		sys.exit(0)
		
	"""
	
	@classmethod
	def logSum(cls, ls):
		import math
		return math.log10(sum(ls))
	
	class ProcessCNVCallSizeHistogram(object):
		"""
		2010-7-21
			used by drawCNVCallSizeHistogram to be plugged into traverseCNVCall()
		"""
		def __init__(self, ecotype_id2nativename=None, **keywords):
			"""
			2010-7-29
				add argument **keywords
			"""
			self.array_id2data = {}
			self.ecotype_id2nativename = ecotype_id2nativename
			self.real_counter = 0
			
			# 2010-7-29
			for keyword, value in keywords.iteritems():
				setattr(self, keyword, value)
		
		def run(self, row, param_obj=None):
			"""
			"""
			from pymodule import PassingData
			import math
			array_id = row.array_id
			ecotype_id = row.ecotype_id
			if array_id not in self.array_id2data:
				nativename = self.ecotype_id2nativename.get(ecotype_id)
				self.array_id2data[array_id] = PassingData(x_ls=[], y_ls=[], ecotype_id=ecotype_id,\
													nativename=nativename)
			self.array_id2data[array_id].x_ls.append(math.log10(row.stop - row.start + 1))
			self.array_id2data[array_id].y_ls.append(row.probability)
			
			self.real_counter += 1
	
	@classmethod
	def drawCNVCallSizeHistogram(cls, db_250k, cnv_method_id=8, cnv_type_id=1, ecotype_id=None, qc_data_source_id=13,\
								output_dir=None, drawType=1, xlabel='log10(segment size)', ylabel_in_2D='score',\
								processClass=None, xlim_in_1D=[1.5, 6], fileNamePrefix='CNVCallSizeHist', **keywords):
		"""
		2010-7-29
			add argument fileNamePrefix, **keywords
		2010-7-21
			add argument xlabel to allow customization by drawCNVCallGapHistogram()
		2010-7-20
		"""
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		if qc_data_source_id is not None and ecotype_id is None:
			ecotype_id_ls = cls.getEcotypeIDListInCNVQCAccessionGivenDataSourceID(qc_data_source_id)
		else:
			ecotype_id_ls = None
		
		from pymodule import PassingData
		from common import get_ecotypeid2nativename
		ecotype_id2nativename = get_ecotypeid2nativename(db_250k.metadata.bind)
		import math
		if processClass is None:
			processClass = cls.ProcessCNVCallSizeHistogram
		processor = processClass(ecotype_id2nativename=ecotype_id2nativename, **keywords)
		
		cls.traverseCNVCall(db_250k, cnv_method_id=cnv_method_id, cnv_type_id=cnv_type_id, \
					ecotype_id=ecotype_id, ecotype_id_ls=ecotype_id_ls,\
					min_no_of_probes=None, min_segment_size=None, processClassIns=processor, \
					param_obj=None)
		from pymodule.utils import addExtraLsToFilenamePrefix
		
		for array_id, data in processor.array_id2data.iteritems():
			sys.stderr.write("array %s ..."%array_id)
			x_ls = data.x_ls
			y_ls = data.y_ls
			
			import pylab
			pylab.clf()
			output_fname = os.path.join(output_dir, '%s_CNVMethod%s_CNVType%s_a_id%s_e_id%s_%s.png'%(fileNamePrefix, \
										cnv_method_id, \
										cnv_type_id, array_id, data.ecotype_id, data.nativename))
			title = 'a-id %s, e-id %s, %s, %s type %s objects.'%(array_id, data.ecotype_id, data.nativename, \
														len(x_ls), cnv_type_id)
			
			pylab.title(title)
			"""
			import statistics	# 2010-5-30 package from Michiel De Hoon
			y, x = statistics.pdf(x_ls)
			pylab.loglog(x, y, alpha=0.7)
			pylab.grid(True, alpha=0.6)
			"""
			pylab.hist(x_ls, 20)
			if xlim_in_1D:
				pylab.xlim(xlim_in_1D)
			
			pylab.xlabel(xlabel)
			pylab.ylabel('Count')
			#
			pylab.savefig(output_fname, dpi=300)
			
			C_ls = [1]*len(y_ls)
			fig_fname = addExtraLsToFilenamePrefix(output_fname, ['2D'])
			cls.drawHexbin(x_ls, y_ls, C_ls, fig_fname=fig_fname, gridsize=20, \
					title=title, \
					xlabel=xlabel, \
					ylabel=ylabel_in_2D, \
					colorBarLabel='log10(count)', reduce_C_function=CNV.logSum)
			sys.stderr.write("Done.\n")
	"""
		# 2010-7-20
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVCall_SegmentSizeHist/')
		CNV.drawCNVCallSizeHistogram(db_250k, cnv_method_id=8, cnv_type_id=1, ecotype_id=None, qc_data_source_id=13, \
								output_dir=output_dir, drawType=1)
		sys.exit(0)
		
		
		# 2010-7-20
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVCall_SegmentSizeHist/')
		CNV.drawCNVCallSizeHistogram(db_250k, cnv_method_id=8, cnv_type_id=1, ecotype_id=None, qc_data_source_id=13, \
								output_dir=output_dir, drawType=2)
		sys.exit(0)
		
		# 2010-7-20
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVCallSegmentSizeHist/')
		CNV.drawCNVCallSizeHistogram(db_250k, cnv_method_id=8, cnv_type_id=1, ecotype_id=None, qc_data_source_id=13, \
								output_dir=output_dir, drawType=1)
		sys.exit(0)
		
	"""
	
	@classmethod
	def drawCNVSizeHistogram(cls, db_250k, cnv_method_id=8, cnv_type_id=1, \
								output_dir=None, drawType=1, xlabel='log(deletion length)', ylabel_in_2D='score',\
								processClass=None, xlim_in_1D=[1.5, 6], fileNamePrefix='CNVSizeHist', **keywords):
		"""
		2010-10-23
			draw histogram of size of cnvs in table CNV
		"""
		sys.stderr.write("Drawing CNV size histogram ...")
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		
		from pymodule import PassingData
		import math, numpy
		
		import Stock_250kDB
		query = Stock_250kDB.CNV.query.filter_by(cnv_method_id=cnv_method_id).filter_by(cnv_type_id=cnv_type_id)
		cnv_size_ls = []
		score_ls = []
		for row in query:
			cnv_size_ls.append(math.log10(row.stop-row.start+1))
			score_ls.append(row.score)
		median_size=  numpy.median(cnv_size_ls)
		sys.stderr.write("%s CNVs with median size %s.\n"%(len(cnv_size_ls), median_size))
		
		from pymodule.utils import addExtraLsToFilenamePrefix
		
		import pylab
		pylab.clf()
		output_fname = os.path.join(output_dir, '%s_CNVMethod%s_CNVType%s.png'%(fileNamePrefix, \
									cnv_method_id, cnv_type_id, ))
		title = ' %s method %s type %s deletions. median %s'%(len(cnv_size_ls), cnv_method_id, cnv_type_id, median_size)
		
		pylab.title(title)
		"""
		import statistics	# 2010-5-30 package from Michiel De Hoon
		y, x = statistics.pdf(x_ls)
		pylab.loglog(x, y, alpha=0.7)
		pylab.grid(True, alpha=0.6)
		"""
		pylab.hist(cnv_size_ls, 20,)
		if xlim_in_1D:
			pylab.xlim(xlim_in_1D)
		
		pylab.xlabel(xlabel)
		pylab.ylabel('Count')
		#
		pylab.savefig(output_fname, dpi=300)
		
		"""
		C_ls = [1]*len(y_ls)
		fig_fname = addExtraLsToFilenamePrefix(output_fname, ['2D'])
		cls.drawHexbin(x_ls, y_ls, C_ls, fig_fname=fig_fname, gridsize=20, \
				title=title, \
				xlabel=xlabel, \
				ylabel=ylabel_in_2D, \
				colorBarLabel='log10(count)', reduce_C_function=CNV.logSum)
		"""
		sys.stderr.write("Done.\n")
	"""
		# 2010-7-20
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVSegmentSizeHist/')
		CNV.drawCNVSizeHistogram(db_250k, cnv_method_id=22, cnv_type_id=1, \
								output_dir=output_dir, drawType=1, xlim_in_1D=None)
		sys.exit(0)
		
		# 2010-7-20
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVSegmentSizeHist/')
		CNV.drawCNVSizeHistogram(db_250k, cnv_method_id=19, cnv_type_id=1, \
								output_dir=output_dir, drawType=1, xlim_in_1D=None)
		sys.exit(0)
	"""
	
	@classmethod
	def drawCNVQCCallSizeHistogram(cls, db_250k, qc_data_source_id=13, cnv_type_id=1, ecotype_id=None, \
								output_dir=None, drawType=1, parseNumOfReadsFromComment=False):
		"""
		2010-7-20
			draw CNVQCCall segment size histogram
			
			drawType
				1: 1D histogram of segment size
				2: 2D histogram of segment size X score
					if parseNumOfReadsFromComment is True, score is replaced with number of reads.
		"""
		from pymodule.CNV import CNVSegmentBinarySearchTreeKey, leftWithinRightAlsoEqualCmp, rightWithinLeftAlsoEqualCmp
		cnvQCData = cls.getCNVQCDataFromDB(db_250k, data_source_id=qc_data_source_id, \
										ecotype_id=ecotype_id, cnv_type_id=cnv_type_id, \
										min_QC_segment_size=None, min_no_of_probes=None,\
										min_reciprocal_overlap=0.0000001, \
										cmpfn=rightWithinLeftAlsoEqualCmp, parseNumOfReadsFromComment=parseNumOfReadsFromComment)
		ecotype_id2qc_data = cnvQCData.ecotype_id2cnv_qc_call_data
		ecotype_id2data = cnvQCData.ecotype_id2data
		
		from common import get_ecotypeid2nativename
		ecotype_id2nativename = get_ecotypeid2nativename(db_250k.metadata.bind)
		
		x_ls = []
		y_ls = []
		import math
		def getSegmentSize(rb_node):
			x_ls.append(math.log10(rb_node.key.span_ls[1]-rb_node.key.span_ls[0]+1))
			y_ls.append(rb_node.value[7])
		
		
		from pymodule.utils import addExtraLsToFilenamePrefix
		
		for ecotype_id, qc_data in ecotype_id2qc_data.iteritems():
			sys.stderr.write("ecotype %s ..."%ecotype_id)
			x_ls = []
			y_ls = []
			qc_data.traverseTree(getSegmentSize)
			
			import pylab
			pylab.clf()
			nativename = ecotype_id2nativename.get(ecotype_id)
			output_fname = os.path.join(output_dir, 'dataSource%s_CNVType%s_e_id%s_%s.png'%(qc_data_source_id, \
																cnv_type_id, ecotype_id, nativename))
			title = 'e-id %s, %s, %s type %s segments.'%(ecotype_id, nativename, len(x_ls), cnv_type_id)
			
			if drawType==1:
				pylab.title(title)
				"""
				import statistics	# 2010-5-30 package from Michiel De Hoon
				y, x = statistics.pdf(x_ls)
				pylab.loglog(x, y, alpha=0.7)
				pylab.grid(True, alpha=0.6)
				"""
				pylab.hist(x_ls, 20)
				pylab.xlim([1.5, 6])
				
				pylab.xlabel('log10(segment size)')
				pylab.ylabel('Count')
				#
				pylab.savefig(output_fname, dpi=300)
			elif drawType==2:
				C_ls = [1]*len(y_ls)
				if parseNumOfReadsFromComment:
					fig_fname = addExtraLsToFilenamePrefix(output_fname, ['2D_SizeVsNumReads'])
					ylabel='num_Reads'
				else:
					fig_fname = addExtraLsToFilenamePrefix(output_fname, ['2D'])
					ylabel = 'score'
				cls.drawHexbin(x_ls, y_ls, C_ls, fig_fname=fig_fname, gridsize=20, \
						title=title, \
						xlabel='log10(segment size)', \
						ylabel=ylabel, \
						colorBarLabel='log10(count)', reduce_C_function=CNV.logSum)
			sys.stderr.write("Done.\n")
			
			
	
	"""
	
		# 2010-7-20
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVQCCall_SegmentSizeHist/')
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		CNV.drawCNVQCCallSizeHistogram(db_250k, qc_data_source_id=13, cnv_type_id=1, ecotype_id=None, \
								output_dir=output_dir)
		sys.exit(0)
	"""
	@classmethod
	def fillDeletionInCNVQCProbeCall(cls, db_250k, data_source_id=1, commit=False):
		"""
		2009-10-29
			construct CNVQCProbeCall based on CNVQCCall, only deletions
		"""
		sys.stderr.write("Filling in CNVQCProbeCall based on CNVQCCall ...\n")
		session = db_250k.session
		#session.begin()
		
		import Stock_250kDB
		from sqlalchemy import desc, asc
		qc_query = Stock_250kDB.CNVQCCall.query.filter_by(cnv_type_id=1)	# deletion only
		qc_call_data = []
		for qc_call in qc_query:
			qc_call_data.append((qc_call.chromosome, qc_call.start, qc_call.stop, qc_call.id, qc_call.accession_id))
		qc_call_data.sort()
		sys.stderr.write("%s CNVQCCall entries.\n"%(len(qc_call_data)))
		
		
		query = db_250k.metadata.bind.execute("select id, chromosome, position, xpos, ypos, allele, strand from %s where direction is not null order by chromosome, position"%(Stock_250kDB.Probes.table.name))
		# 2009-10-29 sqlalchemy takes up too much memory
		# query = Stock_250kDB.Probes.query.filter(Stock_250kDB.Probes.direction!=None).order_by(asc(Stock_250kDB.Probes.chromosome)).\
		#	order_by(asc(Stock_250kDB.Probes.position))
		starting_index_for_next_probe = 0
		no_of_probes = 0
		no_of_encounters = 0
		for probe in query:
			probe_start_pos = probe.position-12
			probe_stop_pos = probe.position+12
			no_of_segments_containing_this_probe = 0
			for i in range(starting_index_for_next_probe, len(qc_call_data)):
				chromosome, start, stop, qc_call_id, accession_id= qc_call_data[i]
				if chromosome==probe.chromosome and start<=probe_stop_pos and stop>=probe_start_pos:	# probe within
					no_of_segments_containing_this_probe += 1
					if no_of_segments_containing_this_probe==1:	# first segment encountering this probe, next probe should at least start from here.
						starting_index_for_next_probe = i
					target_position = probe_start_pos-start
					if start>=probe_start_pos:
						size_affected = 25-(start-probe_start_pos)
					elif stop<=probe_stop_pos:
						size_affected = 25-(probe_stop_pos-stop)
					else:
						size_affected = 25
					db_entry = Stock_250kDB.CNVQCProbeCall(accession_id=accession_id, probe_id=probe.id, chromosome=probe.chromosome, \
														position=probe.position, \
														size_affected=size_affected, target_position=target_position, \
														cnv_qc_call_id = qc_call_id, cnv_type_id=1)
					session.add(db_entry)
					session.flush()
					no_of_encounters += 1
				elif chromosome<probe.chromosome:	# skip
					continue
				elif chromosome>probe.chromosome or (chromosome==probe.chromosome and start>probe_stop_pos):	# all following segments is beyond this probe. exit
					break
			if no_of_segments_containing_this_probe==0:
				starting_index_for_next_probe = i-1		# this probe is not in any segment. Next probe should start from the (i-1)th segment,
											# which is the last segment ahead of the current probe.
			no_of_probes += 1
			if no_of_probes%5000==0:
				sys.stderr.write("%sNo of probes: %s\tNo of encounters: %s"%('\x08'*80, no_of_probes, no_of_encounters))
		#if commit:
		#	session.commit()
		sys.stderr.write("No of probes: %s\tNo of encounters: %s. Done.\n"%(no_of_probes, no_of_encounters))
	
	"""
	CNV.fillDeletionInCNVQCProbeCall(db_250k, commit=True)
	"""
	
	@classmethod
	def processColName2IndexForArrayInfo(cls, col_name2index):
		"""
		2009-11-20
		# construct ecotype_id2array_id_ls mapper
		"""
		import sys, os
		sys.stderr.write("Constructing ecotype_id2array_id_ls ...")
		from pymodule import PassingData
		import Stock_250kDB, re
		all_number_pattern = re.compile(r'^\d+$')
		array_id2index = {}
		ecotype_id2array_id_ls = {}
		array_label_ls = []
		for col_name in col_name2index:
			if all_number_pattern.search(col_name):
				array_id = int(col_name)
				array = Stock_250kDB.ArrayInfo.get(array_id)
				array_id2index[array_id] = col_name2index[col_name]
				ecotype_id = array.maternal_ecotype_id
				array_label_ls.append('e-%s %s'%(ecotype_id, array_id))
				if ecotype_id not in ecotype_id2array_id_ls:
					ecotype_id2array_id_ls[ecotype_id] = []
				ecotype_id2array_id_ls[ecotype_id].append(array_id)
		array_index_ls = array_id2index.values()
		array_index_ls.sort()
		min_array_index = min(array_index_ls)
		max_array_index = max(array_index_ls)
		returnData = PassingData(array_id2index=array_id2index, ecotype_id2array_id_ls=ecotype_id2array_id_ls,\
								min_array_index=min_array_index, max_array_index=max_array_index, array_label_ls=array_label_ls)
		sys.stderr.write("Done.\n")
		return returnData
	
	@classmethod
	def calculateProbeSTDDEV(cls, db_250k, input_fname_ls, output_fname, min_no_of_replicates=5):
		"""
		2009-11-19
			calculate stddev of one probe in multiple replicate arrays
		"""
		import fileinput
		from pymodule import getColName2IndexFromHeader, PassingData
		array_id2array = {}
		counter = 0
		real_counter = 0
		
		input_handler = fileinput.input(input_fname_ls)
		header = input_handler.readline().strip().split('\t')
		col_name2index = getColName2IndexFromHeader(header)
		
		returnData = cls.processColName2IndexForArrayInfo(col_name2index)
		ecotype_id2array_id_ls = returnData.ecotype_id2array_id_ls
		min_array_index = returnData.min_array_index
		max_array_index = returnData.max_array_index
		array_id2index = returnData.array_id2index
		
		# only to calculate stddev of probes of ecotypes who have >=min_no_of_replicates arrays.
		ecotype_id_to_check = []
		for ecotype_id, array_id_ls in ecotype_id2array_id_ls.iteritems():
			if len(array_id_ls)>=min_no_of_replicates:
				ecotype_id_to_check.append(ecotype_id)
		
		sys.stderr.write("Calculating stddev of probe intensity from replicate arrays from %s ... \n"%repr(input_fname_ls))
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = ['probe_id', 'chr_pos', 'ecotype_id', 'sample_size', 'stddev']
		writer.writerow(header)
		import numpy
		from annot.bin.codense.common import dict_map
		for line in input_handler:
			if line.find("probes_id")!=-1:
				continue
			line = line.strip()
			row = line.split('\t')
			probe_id_idx = col_name2index.get('probes_id')
			probe_id = int(row[probe_id_idx])
			chr = row[col_name2index.get('chromosome')]
			pos = row[col_name2index.get('position')]
			chr_pos = '%s_%s'%(chr, pos)
			total_probe_intensity_ls = row[min_array_index:max_array_index+1]
			total_probe_intensity_ls = map(float, total_probe_intensity_ls)
			total_probe_intensity_ls = numpy.array(total_probe_intensity_ls)
			for ecotype_id in ecotype_id_to_check:
				array_id_ls = ecotype_id2array_id_ls[ecotype_id]
				array_index_ls = dict_map(array_id2index, array_id_ls)
				probe_intensity_ls = total_probe_intensity_ls[array_index_ls]
				stddev = numpy.std(probe_intensity_ls)
				output_row = [probe_id, chr_pos, ecotype_id, len(probe_intensity_ls), stddev]
				writer.writerow(output_row)
			
			stddev = numpy.std(total_probe_intensity_ls)
			output_row = [probe_id, chr_pos, 'all', len(total_probe_intensity_ls), stddev]
			writer.writerow(output_row)
			
			counter += 1			
			if counter%10000==0:
				sys.stderr.write('%s%s'%('\x08'*80, counter))
		del writer
		sys.stderr.write("Done.\n")
	"""
	input_fname_ls = [os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity.tsv')]
	output_fname = '/Network/Data/250k/tmp-yh/CNV/call_method_48_CNV_probe_stddev.tsv'
	min_no_of_replicates = 5
	CNV.calculateProbeSTDDEV(db_250k, input_fname_ls, output_fname, min_no_of_replicates=min_no_of_replicates)
	
	
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s.tsv'%i))
	output_fname = '/Network/Data/250k/tmp-yh/CNV/call_method_48_CNV_QNorm_sub_ref_probe_stddev.tsv'
	min_no_of_replicates = 5
	CNV.calculateProbeSTDDEV(db_250k, input_fname_ls, output_fname, min_no_of_replicates=min_no_of_replicates)
	"""
	
	@classmethod
	def calculateProbeQuartilePerChromosome(cls, db_250k, output_fname_prefix=None):
		"""
		2009-11-25
			a newer and expanded version of this function is moved to DB_250k2Array.py
		2009-11-20
			calculate CNV probe intensity quartile (1st, median, 3rd) for each chromosome and store them in database
		"""
		import sys, os, numpy, math
		from scipy import stats	# for scoreatpercentile/percentileatscore to get quartiles
		import Stock_250kDB
		import rpy
		rpy.r.library('affy')
		from DB_250k2Array import DB_250k2Array
		session = db_250k.session
		
		probes, xy_ls, chr_pos_ls, probe_id_ls = DB_250k2Array.get_probes(db_250k.metadata.bind, Stock_250kDB.Probes.table.name, \
																		snps=None, run_type=2)
		
		
		sys.stderr.write("Getting probes into each chromosome ...")
		chr2xy_ls = {}
		chr2probe_id_ls = {}
		for i in range(len(xy_ls)):
			chr,pos = chr_pos_ls[i]
			if chr not in chr2xy_ls:
				chr2xy_ls[chr] = []
				chr2probe_id_ls[chr] = []	#initialize with the start_probe_id
			chr2xy_ls[chr].append(xy_ls[i])
			chr2probe_id_ls[chr].append(probe_id_ls[i])
		sys.stderr.write("Done.\n")
		
		
		rows = Stock_250kDB.ArrayInfo.query.all()
		no_of_arrays = rows.count()
		array_width = None
		counter = 0
		for row in rows:
			array_id = row.id
			filename = row.filename
			ecotype_id = row.maternal_ecotype_id
						
			sys.stderr.write("\t%d/%d: Extracting intensity from %s ... "%(counter+1, no_of_arrays, array_id))
			
			#read array by calling R
			array = rpy.r.read_affybatch(filenames=filename)
			intensity_array = rpy.r.intensity(array)	#return a lengthX1 2-Dimensional array.
			if array_width == None:
				intensity_array_size = len(intensity_array)
				array_width = int(math.sqrt(intensity_array_size))	#assume it's square array
			
			chr2intensity_ls = {}
			for chr, chr_xy_ls in chr2xy_ls.iteritems():
				chr2intensity_ls[chr] = []
				for xpos, ypos in chr_xy_ls:
					intensity_array_index = array_width*(array_width - xpos - 1) + ypos
					intensity = math.log10(intensity_array[intensity_array_index][0])
					chr2intensity_ls[chr].append(intensity)
				array_quartile = Stock_250kDB.ArrayQuartile(array_id=array_id, start_probe_id=chr2probe_id_ls[chr][0],\
														stop_probe_id=chr2probe_id_ls[chr][-1], \
														no_of_probes=len(chr2intensity_ls[chr]))
				array_quartile.minimum = numpy.min(chr2intensity_ls[chr])
				array_quartile.first_decile = stats.scoreatpercentile(chr2intensity_ls[chr], 10)
				array_quartile.lower_quartile = stats.scoreatpercentile(chr2intensity_ls[chr], 25)
				array_quartile.median = stats.scoreatpercentile(chr2intensity_ls[chr], 50)
				array_quartile.upper_quartile = stats.scoreatpercentile(chr2intensity_ls[chr], 75)
				array_quartile.last_decile = stats.scoreatpercentile(chr2intensity_ls[chr], 90)
				array_quartile.maximum = numpy.max(chr2intensity_ls[chr])
				
				session.add(array_quartile)
				
				# find and store the outliers
				IQR = array_quartile.upper_quartile - array_quartile.lower_quartile
				lower_whisker = array_quartile.lower_quartile - 1.5*IQR
				upper_whisker = array_quartile.upper_quartile + 1.5*IQR
				for i in range(len(chr2intensity_ls[chr])):
					intensity = chr2intensity_ls[chr][i]
					if intensity<lower_whisker or intensity>upper_whisker:
						probe_id = chr2probe_id_ls[chr][i]
						array_quartile_outlier = Stock_250kDB.ArrayQuartileOutlier(probe_id=probe_id, intensity=intensity)
						array_quartile_outlier.array_quartile = array_quartile
						session.add(array_quartile_outlier)
				session.flush()
			counter += 1
			sys.stderr.write("Done.\n")
		
		"""
		import fileinput
		from pymodule import getColName2IndexFromHeader, PassingData
		array_id2array = {}
		counter = 0
		real_counter = 0
		
		input_handler = fileinput.input(input_fname_ls)
		header = input_handler.readline().strip().split('\t')
		col_name2index = getColName2IndexFromHeader(header)
		
		returnData = cls.processColName2IndexForArrayInfo(col_name2index)
		min_array_index = returnData.min_array_index
		max_array_index = returnData.max_array_index
		
		output_dir = os.path.split(output_fname_prefix)[0]
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		
		sys.stderr.write("Drawing box plots of CNV arrays for each chromosome from %s ... \n"%repr(input_fname_ls))
		probe_intensity_matrix = []
		for i in range(len(returnData.array_label_ls)):
			probe_intensity_matrix.append([])
		import numpy, pylab
		from annot.bin.codense.common import dict_map
		prev_chromosome = None
		
		for line in input_handler:
			if line.find("probes_id")!=-1:
				continue
			line = line.strip()
			row = line.split('\t')
			probe_id_idx = col_name2index.get('probes_id')
			probe_id = int(row[probe_id_idx])
			chr = row[col_name2index.get('chromosome')]
			if prev_chromosome==None:
				prev_chromosome = chr
			#pos = row[col_name2index.get('position')]
			#chr_pos = '%s_%s'%(chr, pos)
			if (chr!=prev_chromosome and len(probe_intensity_matrix[0])>20 and counter>0) or (len(probe_intensity_matrix[0])%90000==0 and counter>0):
				sys.stderr.write("boxplot chr %s at %s ..."%(prev_chromosome, counter))
				for i in range(len(returnData.array_label_ls)):
					pylab.clf()
					pylab.boxplot(probe_intensity_matrix[i])
					pylab.title(returnData.array_label_ls[i])
					pylab.ylim([1,5])
					pylab.savefig('%s_%s_chr%s_%s.png'%(output_fname_prefix, returnData.array_label_ls[i], prev_chromosome, counter), dpi=300)
				sys.stderr.write("Done.\n")
				
				#clean up the matrix
				#del probe_intensity_matrix
				#probe_intensity_matrix = []
				for i in range(len(returnData.array_label_ls)):
					#probe_intensity_matrix.append([])
					probe_intensity_matrix[i] = []
			
			if chr!=prev_chromosome:
				prev_chromosome = chr
				for i in range(len(returnData.array_label_ls)):
					probe_intensity_matrix[i] = []
				
			total_probe_intensity_ls = row[min_array_index:max_array_index+1]
			total_probe_intensity_ls = map(float, total_probe_intensity_ls)
			for i in range(len(total_probe_intensity_ls)):
				probe_intensity_matrix[i].append(total_probe_intensity_ls[i])
			
			counter += 1			
			if counter%10000==0:
				sys.stderr.write('%s%s'%('\x08'*80, counter))
		sys.stderr.write("boxplot chr %s at %s ..."%(chr, counter))
		for i in range(len(returnData.array_label_ls)):
			pylab.clf()
			pylab.boxplot(probe_intensity_matrix[i])
			pylab.title(returnData.array_label_ls[i])
			pylab.ylim([1,5])
			pylab.savefig('%s_%s_chr%s_%s.png'%(output_fname_prefix, returnData.array_label_ls[i], prev_chromosome, counter), dpi=300)
		sys.stderr.write("Done.\n")
		"""
	
	"""
	output_fname_prefix = '/Network/Data/250k/tmp-yh/CNV/CNV_intensity_boxplot'
	CNV.calculateProbeQuartilePerChromosome(db_250k, output_fname_prefix)
	"""
	
	@classmethod
	def addDeletionIntoRBDict(cls, param_obj, cnv_segment_obj, cnv_qc_call=None):
		"""
		2010-6-1
		
		cnv_segment_obj = PassingData(ecotype_id=cnv_ecotype_id, start_probe=start_probe, stop_probe=stop_probe,\
											no_of_probes=no_of_probes, amplitude=amplitude, segment_length=segment_length,\
											segment_chromosome=segment_chromosome, array_id=array_id,\
											start_probe_id=start_probe_id, stop_probe_id=stop_probe_id,\
											segment_start_pos=segment_start_pos, segment_stop_pos=segment_stop_pos)
		"""
		from pymodule.CNV import CNVSegmentBinarySearchTreeKey, leftWithinRightAlsoEqualCmp
		deletionDict = param_obj.deletionDict
		segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=cnv_segment_obj.segment_chromosome, \
							span_ls=[cnv_segment_obj.segment_start_pos, cnv_segment_obj.segment_stop_pos], \
									min_reciprocal_overlap=param_obj.min_reciprocal_overlap)
		if segmentKey not in deletionDict:
			deletionDict[segmentKey] = []
		ls = deletionDict.get(segmentKey)
		try:
			ls.append(cnv_segment_obj)
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
	
	@classmethod
	def calDeletionFrequencySpectrum(cls, db_250k, input_fname_ls, output_fname_prefix, data_source_id=1, cnv_type_id=1, \
								min_QC_segment_size=200, min_no_of_probes=5, max_boundary_diff=10000, \
								max_diff_perc=0.10, count_embedded_segment_as_match=True, min_reciprocal_overlap=0.6,\
								deletion_cutoff=-0.4):
		"""
		2010-6-1
		"""
		
		from pymodule import PassingData
		from pymodule.RBTree import RBDict	# 2010-1-26 RBDict is more efficient than binary_tree.
		from pymodule.CNV import CNVSegmentBinarySearchTreeKey, leftWithinRightAlsoEqualCmp
		deletionDict = RBDict(cmpfn=leftWithinRightAlsoEqualCmp)
		
		param_obj = PassingData(deletionDict=deletionDict, min_reciprocal_overlap=min_reciprocal_overlap)
		ecotype_id2cnv_qc_call_data = None
		cls.compareCNVSegmentsAgainstQCHandler(input_fname_ls, ecotype_id2cnv_qc_call_data, cls.addDeletionIntoRBDict, param_obj, \
											deletion_cutoff=deletion_cutoff, max_boundary_diff=max_boundary_diff, \
											max_diff_perc=max_diff_perc,\
											min_no_of_probes=min_no_of_probes, \
											count_embedded_segment_as_match=count_embedded_segment_as_match)
		
		sys.stderr.write("data_source_id %s, min_QC_segment_size %s, min_no_of_probes: %s, max_boundary_diff: %s, max_diff_perc: %s, count_embedded_segment_as_match: %s.\n"%\
						(data_source_id, min_QC_segment_size, min_no_of_probes, max_boundary_diff, max_diff_perc,\
						count_embedded_segment_as_match))
		count_ls = [len(deletion_ls) for deletion_ls in deletionDict]
		sys.stderr.write("Drawing ...")
		import pylab
		pylab.clf()
		import statistics	# 2010-5-30 package from Michiel De Hoon
		y, x = statistics.pdf(count_ls)
		pylab.grid(True, alpha=0.3)
		pylab.title("density deletion count")
		pylab.plot(x, y, alpha=0.7)
		#pylab.title("Histogram of count of deletions")
		#pylab.hist(count_ls, 20)
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		sys.stderr.write("Done.\n")
	
	"""
	
	input_fname = os.path.expanduser('~/script/variation/data/CNV/call_53_WABQN_b100_j50_lts_98_arrays_GADA_output_A0.5T4.0M5.tsv')
	input_fname_ls = [input_fname]
	from pymodule.utils import addExtraToFilenamePrefix, addExtraLsToFilenamePrefix
	output_fname_prefix = os.path.splitext(addExtraToFilenamePrefix(input_fname, 'count_hist'))[0]
	CNV.calDeletionFrequencySpectrum(db_250k, input_fname_ls, output_fname_prefix, data_source_id=1, cnv_type_id=1, \
								min_QC_segment_size=200, min_no_of_probes=5, max_boundary_diff=10000, \
								max_diff_perc=0.10, count_embedded_segment_as_match=True, min_reciprocal_overlap=0.6,\
								deletion_cutoff=-0.4)
	sys.exit(0)
	
	"""
	
	class QuanLongPairedEndSolexaData(object):
		"""
		2010-7-22
			the class has functions dealing with Quan's PE data
		"""
		
		def __init__(self):
			pass
		
		@classmethod
		def discoverDeletionOrNonDeletionFromCoverageData(cls, input_dir=None, output_fname=None, max_coverage_deemed_missing=1.0,\
											copy_number=0, windowSize=100):
			"""
			2010-7-22
				any 100bp window whose coverage is below max_coverage_deemed_missing is regarded as deletion/missing.
				
				argument copy_number controls whether to discover deletions (=0) or non-deletions (>0).
			"""
			import os, sys, csv
			sys.stderr.write("Discovering deletions from Quan's PE coverage data ...\n")
			writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
			header = ['original_id','ecotype_id', 'chromosome', 'start', 'stop', 'size_affected', 'copy_number', 'mean_coverage', ]
			writer.writerow(header)
			for fname in os.listdir(input_dir):
				sys.stderr.write('%s ...'%fname)
				missing_start_pos = None
				last_missing_pos = None
				coverage_sum = 0.0
				input_fname = os.path.join(input_dir, fname)
				reader = csv.reader(open(input_fname), delimiter='\t')
				
				original_id, ecotype_id, chr = fname.split('.')[:3]
				ecotype_id = int(ecotype_id)
				count = 0
				for row in reader:
					position, coverage = row[:2]
					position = int(position)*windowSize	#coverage is caculated in 100-base window.
					coverage = float(coverage)
					
					if copy_number==0:
						condition = coverage<=max_coverage_deemed_missing
					elif copy_number>0:
						condition = coverage>max_coverage_deemed_missing
					
					if condition:
						if missing_start_pos is None:
							missing_start_pos = position - (windowSize-1)
						last_missing_pos = position
						coverage_sum += coverage*windowSize
					else:
						if last_missing_pos is not None and missing_start_pos is not None:
							segment_length = last_missing_pos-missing_start_pos+1
							mean_coverage = coverage_sum/float(segment_length)
							data_row = [original_id, ecotype_id, chr, missing_start_pos, last_missing_pos, segment_length, \
									copy_number, mean_coverage]
							writer.writerow(data_row)
							count += 1
						missing_start_pos = None
						last_missing_pos = None
						coverage_sum = 0.0
				# catch the deletion at the end of the chromosome
				if last_missing_pos is not None and missing_start_pos is not None:
					segment_length = last_missing_pos-missing_start_pos+1
					mean_coverage = coverage_sum/float(segment_length)
					data_row = [original_id, ecotype_id, chr, missing_start_pos, last_missing_pos,  segment_length, \
							copy_number, mean_coverage]
					writer.writerow(data_row)
					count += 1
				del reader
				sys.stderr.write(' %s deletions. Done.\n'%count)
			del writer
		"""
		# 2010-7-22 discover deletions out of coverage data
		input_dir = os.path.expanduser('~/script/variation/data/CNV/coverage4Yu/')
		output_fname = os.path.expanduser('~/script/variation/data/CNV/coverage4Yu_deletions.tsv')
		CNV.QuanLongPairedEndSolexaData.discoverDeletionOrNonDeletionFromCoverageData(input_dir=input_dir, \
										output_fname=output_fname, max_coverage_deemed_missing=1.0)
		sys.exit(0)
		
		
		# 2010-7-22 discover normal segments out of coverage data
		input_dir = os.path.expanduser('~/script/variation/data/CNV/coverage4Yu/')
		output_fname = os.path.expanduser('~/script/variation/data/CNV/coverage4Yu_NonDeletions.tsv')
		CNV.QuanLongPairedEndSolexaData.discoverDeletionOrNonDeletionFromCoverageData(input_dir=input_dir, \
										output_fname=output_fname, max_coverage_deemed_missing=1.0, copy_number=1)
		sys.exit(0)
		"""
		
		
		
		@classmethod
		def drawCoverageHistogram(cls, input_dir, output_dir, xlim_in_1D=None):
			"""
			2010-7-22
			"""
			import csv, os
			if not os.path.isdir(output_dir):
				os.makedirs(output_dir)
			
			for fname in os.listdir(input_dir):
				sys.stderr.write('%s ...'%fname)
				input_fname = os.path.join(input_dir, fname)
				reader=csv.reader(open(input_fname, 'r'), delimiter='\t')
				coverage_ls = []
				for row in reader:
					coverage = float(row[1])
					if coverage>0:
						coverage_ls.append(math.log10(coverage))
					else:
						coverage_ls.append(-3)
				
				import pylab
				pylab.clf()
				title = '%s'%os.path.splitext(fname)[0]
				xlabel = 'log10(coverage)'
				pylab.title(title)
				"""
				import statistics	# 2010-5-30 package from Michiel De Hoon
				y, x = statistics.pdf(x_ls)
				pylab.loglog(x, y, alpha=0.7)
				pylab.grid(True, alpha=0.6)
				"""
				pylab.hist(coverage_ls, 20, log=True)
				if xlim_in_1D:
					pylab.xlim(xlim_in_1D)
				
				pylab.xlabel(xlabel)
				pylab.ylabel('Count')
				output_fname = os.path.join(output_dir, '%s_coverage_hist.png'%os.path.splitext(fname)[0])
				pylab.savefig(output_fname, dpi=300)
				sys.stderr.write('Done.\n')
			
		"""
		input_dir = os.path.expanduser('~/script/variation/data/CNV/coverage4Yu/')
		output_dir = os.path.expanduser('~/script/variation/data/CNV/coverage4Yu_hist/')
		CNV.QuanLongPairedEndSolexaData.drawCoverageHistogram(input_dir, output_dir, xlim_in_1D=None)
		sys.exit(0)
		
		# 2010-7-22
		input_dir = os.path.expanduser('~/script/variation/data/CNV/coverage4Yu/')
		output_dir = os.path.expanduser('~/script/variation/data/CNV/coverage4Yu_hist/')
		CNV.QuanLongPairedEndSolexaData.drawCoverageHistogram(input_dir, output_dir, xlim_in_1D=None)
		sys.exit(0)
		"""
	
	class OneCNV(object):
		def __init__(self, max_boundary_diff=20000, max_diff_perc=0.2, min_reciprocal_overlap=0.6):
			self.chromosome = None
			self.start = None
			self.stop = None
			self.segment_length = None
			self.cnv_ls = []
			self.array_id_ls = []
			self.max_boundary_diff = max_boundary_diff
			self.max_diff_perc = max_diff_perc
			self.min_reciprocal_overlap = min_reciprocal_overlap
		
		def addFirstCNV(self, chromosome, start, stop, array_id=None):
			self.addOneCNV(chromosome, start, stop, array_id)
		
		def addOneCNV(self, chromosome, start, stop, array_id=None):
			"""
			"""
			if self.chromosome is not None and chromosome!=self.chromosome:
				return False
			
			if self.chromosome is None:
				self.chromosome = chromosome
			
			self.cnv_ls.append((start, stop))
			self.array_id_ls.append(array_id)
			self.adjustBoundary()
			return True
		
		def adjustBoundary(self):
			"""
			"""
			import numpy
			if len(self.cnv_ls)>0:
				self.start, self.stop = numpy.mean(self.cnv_ls, 0)
				self.segment_length = abs(self.stop-self.start)
		
		def __len__(self):
			return len(self.cnv_ls)
		
		def addNewCNV(self, chromosome, start, stop, array_id=None):
			"""
			"""
			if self.chromosome is None:
				self.addOneCNV(chromosome, start, stop, array_id)
			elif self.chromosome is not None and chromosome!=self.chromosome:
				return False
			else:
				"""
				boundary_diff1 = abs(start-self.start)
				boundary_diff2 = abs(stop-self.stop)
				diff1_perc = boundary_diff1/float(self.segment_length)
				diff2_perc = boundary_diff2/float(self.segment_length)
				if boundary_diff1<=self.max_boundary_diff and boundary_diff2<=self.max_boundary_diff and \
					diff1_perc<=self.max_diff_perc and diff2_perc<=self.max_diff_perc:
					self.addOneCNV(chromosome, start, stop, array_id)
				
				else:
					return False
				"""
				
				is_overlap = is_reciprocal_overlap([start, stop], [self.start, self.stop], \
													min_reciprocal_overlap=self.min_reciprocal_overlap)
				if is_overlap:
					self.addOneCNV(chromosome, start, stop, array_id)
				else:
					return False
			
	
	@classmethod
	def plotCNVOccurrenceInReplicatesHist(cls, db_250k, cnv_intensity_fname, GADA_output_fname_ls, output_fname_prefix, \
										min_no_of_replicates=5, min_no_of_probes=5,\
										deletion_cutoff=None, max_boundary_diff=10000, \
										max_diff_perc=0.10, min_reciprocal_overlap=0.6):
		"""
		2009-11-19
			calculate stddev of one probe in multiple replicate arrays
		"""
		import fileinput
		from pymodule import getColName2IndexFromHeader, PassingData
		array_id2array = {}
		counter = 0
		real_counter = 0
		
		input_handler = fileinput.input([cnv_intensity_fname])
		header = input_handler.readline().strip().split('\t')
		col_name2index = getColName2IndexFromHeader(header)
		
		returnData = cls.processColName2IndexForArrayInfo(col_name2index)
		ecotype_id2array_id_ls = returnData.ecotype_id2array_id_ls
		min_array_index = returnData.min_array_index
		max_array_index = returnData.max_array_index
		array_id2index = returnData.array_id2index
		
		# only to calculate stddev of probes of ecotypes who have >=min_no_of_replicates arrays.
		ecotype_id_to_check = []
		for ecotype_id, array_id_ls in ecotype_id2array_id_ls.iteritems():
			if len(array_id_ls)>=min_no_of_replicates:
				ecotype_id_to_check.append(ecotype_id)
		ecotype_id_to_check_set = set(ecotype_id_to_check)
		fileinput.close()	# must. otherwise next fileinput won't work.
		
		sys.stderr.write("Getting occurrence from replicate arrays from %s ... \n"%repr(GADA_output_fname_ls))
		counter = 0
		real_counter = 0
		no_of_deletions = 0
		input_handler = fileinput.input(GADA_output_fname_ls)
		header = input_handler.readline().strip().split('\t')
		col_name2index = getColName2IndexFromHeader(header)
		ecotype_id2cnv_ls = {}
		for line in input_handler:
			if line.find("array_id")!=-1:
				continue
			line = line.strip()
			row = line.split('\t')
			ecotype_id_idx = col_name2index.get('ecotype_id', col_name2index.get('array_id'))
			cnv_ecotype_id = int(row[ecotype_id_idx])
			array_id = int(row[col_name2index.get('array_id')])
			#row[ecotype_id_idx] = cnv_ecotype_id
			counter += 1
			if cnv_ecotype_id in ecotype_id_to_check_set:	# array is in CNVQCDat
				
				start_probe = row[col_name2index['start_probe']].split('_')	# split chr_pos
				start_probe = map(int, start_probe)
				stop_probe = row[col_name2index['end_probe']].split('_')
				stop_probe = map(int, stop_probe)
				no_of_probes = int(row[col_name2index['length']])
				if no_of_probes<min_no_of_probes:
					continue
				amplitude = float(row[col_name2index['amplitude']])
				segment_chromosome = start_probe[0]
				segment_start_pos = start_probe[1]-12
				segment_stop_pos = stop_probe[1]+12
				segment_length = abs(segment_stop_pos-segment_start_pos+1)
				if deletion_cutoff is not None and amplitude>deletion_cutoff:
					continue
				no_of_deletions+=1
				if cnv_ecotype_id not in ecotype_id2cnv_ls:
					ecotype_id2cnv_ls[cnv_ecotype_id] = []
				cnv_ls = ecotype_id2cnv_ls[cnv_ecotype_id]
				addToExistingCNV = False
				for one_cnv in cnv_ls:
					if one_cnv.addNewCNV(segment_chromosome, segment_start_pos, segment_stop_pos, array_id):
						addToExistingCNV = True
				if not addToExistingCNV:
					one_cnv = cls.OneCNV(max_boundary_diff =max_boundary_diff, max_diff_perc=max_diff_perc, \
										min_reciprocal_overlap=min_reciprocal_overlap)
					one_cnv.addFirstCNV(segment_chromosome, segment_start_pos, segment_stop_pos, array_id)
					cnv_ls.append(one_cnv)
				
			if counter%10000==0:
				sys.stderr.write('%s%s\t%s'%('\x08'*80, counter, no_of_deletions))
		
		for ecotype_id, cnv_ls in ecotype_id2cnv_ls.iteritems():
			occurrence_ls = map(len, cnv_ls)
			import pylab
			pylab.clf()
			no_of_replicates = len(ecotype_id2array_id_ls[ecotype_id])
			pylab.title('%s(%s reps), del-c %s, min #probes %s, max dist %s, perc %s, overlap %s'%(ecotype_id, no_of_replicates, deletion_cutoff, min_no_of_probes, \
																max_boundary_diff, max_diff_perc, min_reciprocal_overlap))
			pylab.hist(occurrence_ls, no_of_replicates)
			
			pylab.xlabel('occurrence of one deletion')
			pylab.ylabel('count')
			pylab.savefig('%s_e_%s.png'%(output_fname_prefix, ecotype_id))
		
		
		sys.stderr.write("Done.\n")
	
	"""
	
	cnv_intensity_fname = os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity.tsv')
	max_boundary_diff = 20000
	max_diff_perc = 0.2
	
	GADA_output_fname_ls = [os.path.expanduser('~/mnt2/panfs/250k/CNV/GADA_output/call_method_48_CNV_intensity_QNorm_sub_ref_chr1_A0.8T8.0M10.tsv')]
	deletion_cutoff = -0.33
	min_reciprocal_overlap = 0.6
	min_no_of_probes = 5
	output_fname_prefix = '/tmp/CNVOccurrence_call_48_chr1_A0.8T8.0M10_delCutoff%s_min_p%s_mdist%s_mperc%s_moverlap%s'%\
						(deletion_cutoff, min_no_of_probes, max_boundary_diff, max_diff_perc, min_reciprocal_overlap)
	CNV.plotCNVOccurrenceInReplicatesHist(db_250k, cnv_intensity_fname, GADA_output_fname_ls, output_fname_prefix, min_no_of_replicates=5, \
										min_no_of_probes=min_no_of_probes, \
										deletion_cutoff=deletion_cutoff, max_boundary_diff=max_boundary_diff, \
										max_diff_perc=max_diff_perc, min_reciprocal_overlap=min_reciprocal_overlap)
	
	
	cnv_intensity_fname = os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity.tsv')
	aAlpha_ls = [0.2, 0.5, 0.8, 1.0, 2.0]
	T_ls = [2.0,4.0,8.0,12.0,14.0]
	M_ls = [5,10]
	deletion_cutoff_ls = [-0.33, -0.4, -0.5, -0.6]
	
	max_boundary_diff = 20000
	max_diff_perc = 0.2
	
	for aAlpha in aAlpha_ls:
		for T in T_ls:
			for M in M_ls:
				for deletion_cutoff in deletion_cutoff_ls:
					for min_no_of_probes in [5,10,20,40]:
						for min_reciprocal_overlap in [0.4, 0.6, 0.8]:
							GADA_output_fname_ls = []
							for chr in range(1,6):
								fname = os.path.expanduser('~/mnt2/panfs/250k/CNV/GADA_output/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s_A%sT%sM%s.tsv'%(chr, aAlpha, T, M))
								if os.path.isfile(fname):
									GADA_output_fname_ls.append(fname)
							output_fname_prefix = '/Network/Data/250k/tmp-yh/CNV/CNVOccurrenceByOverlap/call_48_A%sT%sM%s_delCutoff%s_min_p%s_mdist%s_mperc%s_moverlap%s'%\
								(aAlpha, T, M, deletion_cutoff, min_no_of_probes, max_boundary_diff, max_diff_perc, min_reciprocal_overlap)
							CNV.plotCNVOccurrenceInReplicatesHist(db_250k, cnv_intensity_fname, GADA_output_fname_ls, \
																output_fname_prefix, \
																min_no_of_replicates=5, \
															min_no_of_probes=min_no_of_probes, deletion_cutoff=deletion_cutoff, \
															max_boundary_diff=max_boundary_diff, \
															max_diff_perc=max_diff_perc, min_reciprocal_overlap=min_reciprocal_overlap)
	
	"""
	
	@classmethod
	def drawHistOfDLRSpread(cls, db_250k, output_fname_prefix, call_method_id=None):
		"""
		2009-11-25
		"""
		import os, sys
		sys.stderr.write("Drawing histogram of DLRSpread for call %s ..."%call_method_id)
		import Stock_250kDB
		if call_method_id is not None:
			rows = Stock_250kDB.CallInfo.query.filter_by(method_id=call_method_id)
		else:
			rows = Stock_250kDB.ArrayQuartile.query.filter(Stock_250kDB.ArrayQuartile.dlrspread!=None)
		dlrspread_ls = []
		for row in rows:
			if call_method_id is not None:
				for array_quartile in row.array.array_quartile_ls: 
					if array_quartile.dlrspread is not None:
						dlrspread_ls.append(array_quartile.dlrspread)
			else:
				 dlrspread_ls.append(row.dlrspread)
		
		import pylab
		pylab.clf()
		if call_method_id is not None:
			title = 'arrays from call %s'%call_method_id
		else:
			title = 'all arrays'
		pylab.title(title)
		pylab.hist(dlrspread_ls, 30)
		pylab.xlabel('DLRSpread')
		pylab.ylabel('Count')
		pylab.xlim([0.1,0.5])
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		sys.stderr.write("Done.\n")
	
	"""
	output_fname_prefix = '/tmp/DLRSpreadHistAllArrays'
	CNV.drawHistOfDLRSpread(db_250k, output_fname_prefix)
	
	call_method_id = 48
	output_fname_prefix = '/tmp/DLRSpreadHistCall%sArrays'%call_method_id
	CNV.drawHistOfDLRSpread(db_250k, output_fname_prefix, call_method_id=call_method_id)
	
	output_fname_prefix = '/tmp/DLRSpreadHistAllArrays'
	CNV.drawHistOfDLRSpread(db_250k, output_fname_prefix)
	
	for call_method_id in (43, 32, 29, 48, 49):
		output_fname_prefix = '/tmp/DLRSpreadHistCall%sArrays'%call_method_id
		CNV.drawHistOfDLRSpread(db_250k, output_fname_prefix, call_method_id=call_method_id)
	
	"""
	
	@classmethod
	def drawOneArraySevenNumber(cls, title, output_fname_prefix, seven_number_2D_ls, ylim = None, no_of_outliers_ls=None,\
							median_bar_length = 0.8, quartile_box_width = 0.6, decile_bar_length = 0.4):
		"""
		2010-3-1
			called by drawArrayQuartile()
		"""
		import pylab
		pylab.clf()
		fig = pylab.figure()
		ax = fig.add_subplot(111)
		import matplotlib.patches as patches
		no_of_bars = len(seven_number_2D_ls)
		for i in range(no_of_bars):
			minimum, first_decile, lower_quartile, median, upper_quartile, last_decile, maximum = seven_number_2D_ls[i][:7]
			i += 1
			# quartile box
			box_lower_left_y = lower_quartile
			box_lower_left_x = i - quartile_box_width/2.0
			rect = patches.Rectangle((box_lower_left_x, box_lower_left_y), width=quartile_box_width, \
									height=abs(upper_quartile-lower_quartile),\
									transform=ax.transData, alpha=0.5)
			ax.add_patch(rect)
			#ax.axhspan(lower_quartile, upper_quartile, xmin=i-quartile_box_width/2.0, xmax=i+quartile_box_width/2.0, alpha=0.5,\
			#			transform=ax.transData)	# transform = ax.transData is a must. The default is ax.transAxes.
			# minimum to maximum verticle line
			ax.plot([i,i], [minimum, maximum], color='k')
			# median bar
			ax.plot([i-median_bar_length/2.0, i+median_bar_length/2.0], [median, median], color='r')
			# decile bar
			ax.plot([i-decile_bar_length/2.0, i+decile_bar_length/2.0], [first_decile, first_decile], color='b')
			ax.plot([i-decile_bar_length/2.0, i+decile_bar_length/2.0], [last_decile, last_decile], color='b')
			# if no_of_outliers_ls is not None and len(no_of_outliers_ls)>i: # 2010-3-2 deal with outlier counts later
		
		ax.set_xlim([0,no_of_bars+1])
		
		if ylim is not None and len(ylim)==2:
			ax.set_ylim(ylim)
		if title:
			ax.set_title(title)
		if output_fname_prefix:
			pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
	
	@classmethod
	def drawCNVFrequencyHist(cls, db_250k, output_dir=None, cnv_method_id=None, fileNamePrefix='CNVFrequency',\
							xlabel='frequency', ylabel_in_2D='log10(size)', xlim_in_1D=None, logHist=False,
							minFractionNotCoveredByLyrata=None, maxFractionNotCoveredByLyrata=None):
		"""
		2010-8-4
			add argument
				minFractionNotCoveredByLyrata: to get SFS/sequence frequency spectrum of derived insertions,
					usually >0.5
				maxFractionNotCoveredByLyrata: get SFS of derived deletions, usually <0.5
		2010-8-1
		
		"""
		sys.stderr.write("Getting data related to CNV frequency ...\n")
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		import Stock_250kDB, math
		query = Stock_250kDB.CNV.query.filter_by(cnv_method_id=cnv_method_id)
		extraFileSuffixLs = []
		if minFractionNotCoveredByLyrata is not None:
			query = query.filter(Stock_250kDB.CNV.fractionNotCoveredByLyrata>=minFractionNotCoveredByLyrata)
			extraFileSuffixLs.append('newInsertions_minFractionNotCoveredByLyrata%s'%minFractionNotCoveredByLyrata)
		if maxFractionNotCoveredByLyrata is not None:
			query = query.filter(Stock_250kDB.CNV.fractionNotCoveredByLyrata<=maxFractionNotCoveredByLyrata)
			extraFileSuffixLs.append('newDeletions_maxFractionNotCoveredByLyrata%s'%maxFractionNotCoveredByLyrata)
		x_ls = []
		y_ls = []
		for row in query:
			if minFractionNotCoveredByLyrata is not None:
				x_ls.append(1-row.frequency)
			elif maxFractionNotCoveredByLyrata is not None:
				x_ls.append(row.frequency)
			y_ls.append(math.log10(row.size_affected))
		sys.stderr.write("Done.\n")
		
		sys.stderr.write("Drawing ...")
		import pylab
		from pymodule.utils import addExtraLsToFilenamePrefix
		pylab.clf()
		output_fname = os.path.join(output_dir, '%s_CNVMethod%s.png'%(fileNamePrefix, \
									cnv_method_id,))
		output_fname = addExtraLsToFilenamePrefix(output_fname, extraFileSuffixLs)
		title = '%s method %s objects.'%(len(x_ls), cnv_method_id)
		
		pylab.title(title)
		"""
		import statistics	# 2010-5-30 package from Michiel De Hoon
		y, x = statistics.pdf(x_ls)
		pylab.loglog(x, y, alpha=0.7)
		pylab.grid(True, alpha=0.6)
		"""
		pylab.hist(x_ls, 30, log=logHist)
		if xlim_in_1D:
			pylab.xlim(xlim_in_1D)
		
		pylab.xlabel(xlabel)
		pylab.ylabel('Count')
		#
		if logHist:
			fig_fname = addExtraLsToFilenamePrefix(output_fname, ['logHist'])
		else:
			fig_fname = output_fname
		pylab.savefig(fig_fname, dpi=300)
		
		C_ls = [1]*len(y_ls)
		fig_fname = addExtraLsToFilenamePrefix(output_fname, ['2D'])
		cls.drawHexbin(x_ls, y_ls, C_ls, fig_fname=fig_fname, gridsize=20, \
				title=title, \
				xlabel=xlabel, \
				ylabel=ylabel_in_2D, \
				colorBarLabel='log10(count)', reduce_C_function=CNV.logSum)
		sys.stderr.write("Done.\n")
	
	"""
	
		#2010-8-1 -z banyan.usc.edu
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVFrequencyHist/')
		cnv_method_id=22
		CNV.drawCNVFrequencyHist(db_250k, output_dir, cnv_method_id=cnv_method_id, fileNamePrefix='CNVFrequency',\
							xlabel='frequency', ylabel_in_2D='log10(size)', xlim_in_1D=[0,1], logHist=False,\
							minFractionNotCoveredByLyrata=0.8, maxFractionNotCoveredByLyrata=None)
		CNV.drawCNVFrequencyHist(db_250k, output_dir, cnv_method_id=cnv_method_id, fileNamePrefix='CNVFrequency',\
							xlabel='frequency', ylabel_in_2D='log10(size)', xlim_in_1D=[0,1], logHist=False,\
							minFractionNotCoveredByLyrata=None, maxFractionNotCoveredByLyrata=0.2)
		sys.exit(0)
		
		
		#2010-8-1 -z banyan.usc.edu
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVFrequencyHist/')
		cnv_method_id=20
		CNV.drawCNVFrequencyHist(db_250k, output_dir, cnv_method_id=cnv_method_id, fileNamePrefix='CNVFrequency',\
							xlabel='frequency', ylabel_in_2D='log10(size)', xlim_in_1D=[0,1], logHist=True,\
							minFractionNotCoveredByLyrata=0.8, maxFractionNotCoveredByLyrata=None)
		CNV.drawCNVFrequencyHist(db_250k, output_dir, cnv_method_id=cnv_method_id, fileNamePrefix='CNVFrequency',\
							xlabel='frequency', ylabel_in_2D='log10(size)', xlim_in_1D=[0,1], logHist=True,\
							minFractionNotCoveredByLyrata=None, maxFractionNotCoveredByLyrata=0.2)
		sys.exit(0)
		
		#2010-8-1 -z banyan.usc.edu
		output_dir = os.path.expanduser('~/script/variation/data/CNV/CNVFrequencyHist/')
		cnv_method_id=22
		CNV.drawCNVFrequencyHist(db_250k, output_dir, cnv_method_id=cnv_method_id, fileNamePrefix='CNVFrequency',\
							xlabel='frequency', ylabel_in_2D='log10(size)', xlim_in_1D=[0,1], logHist=True,\
							minFractionNotCoveredByLyrata=0.8, maxFractionNotCoveredByLyrata=None)
		CNV.drawCNVFrequencyHist(db_250k, output_dir, cnv_method_id=cnv_method_id, fileNamePrefix='CNVFrequency',\
							xlabel='frequency', ylabel_in_2D='log10(size)', xlim_in_1D=[0,1], logHist=True,\
							minFractionNotCoveredByLyrata=None, maxFractionNotCoveredByLyrata=0.2)
		sys.exit(0)
	"""
	
	@classmethod
	def drawCNVLDHist(cls, db_250k, input_fname, output_dir=None, cnv_method_id=20, fileNamePrefix='CNVLD',\
							xlabel='LD', ylabel_in_2D='#SNPs per kb', xlim_in_1D=None, logHist=False,
							minFractionNotCoveredByLyrata=None, maxFractionNotCoveredByLyrata=None):
		"""
		2010-10-12
		
		"""
		sys.stderr.write("Getting data related to CNV LD ...\n")
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		import Stock_250kDB, math
		
		from pymodule.utils import getColName2IndexFromHeader
		from pymodule import PassingData
		import csv
		reader = csv.reader(open(input_fname), delimiter='\t')
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header)
		cnv_id2data = {}
		counter = 0
		real_counter = 0
		for line in reader:
			cnv_id = int(line[col_name2index['snp2']])
			r2_value = float(line[col_name2index['r2']])
			if cnv_id not in cnv_id2data:
				cnv_id2data[cnv_id] = PassingData(max_LD=r2_value, no_of_snps=0)
				real_counter += 1
			if r2_value>cnv_id2data[cnv_id].max_LD:
				cnv_id2data[cnv_id].max_LD = r2_value
			cnv_id2data[cnv_id].no_of_snps += 1
			counter += 1
			if counter%10000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
		del reader
		sys.stderr.write("\n%s CNVs with LD.\n"%(len(cnv_id2data)))
		
		
		import Stock_250kDB
		
		counter = 0
		real_counter = 0
		x_ls = []
		y_ls = []
		LD_ls = []
		deletionSizeLs = []
		for cnv_id, data in cnv_id2data.iteritems():
			cnv = Stock_250kDB.CNV.get(cnv_id)
			x_ls.append(cnv.frequency)
			deletionSize =  cnv.stop-cnv.start+1
			if deletionSize>0:
				deletionSizeLs.append(math.log10(deletionSize))
			else:
				deletionSizeLs.append(0)
			regionSize = 40000+deletionSize	#add deletionSize
			y_ls.append(data.no_of_snps/(regionSize/1000.))
			LD_ls.append(data.max_LD)
			counter += 1
			if counter%5000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
		sys.stderr.write("Done.\n")
		
		sys.stderr.write("Drawing ...")
		import pylab
		from pymodule.utils import addExtraLsToFilenamePrefix
		pylab.clf()
		output_fname = os.path.join(output_dir, '%s_CNVMethod%s.png'%(fileNamePrefix, \
									cnv_method_id,))
		extraFileSuffixLs = []
		output_fname = addExtraLsToFilenamePrefix(output_fname, extraFileSuffixLs)
		title = '%s method %s objects.'%(len(x_ls), cnv_method_id)
		
		pylab.title(title)
		"""
		import statistics	# 2010-5-30 package from Michiel De Hoon
		y, x = statistics.pdf(x_ls)
		pylab.loglog(x, y, alpha=0.7)
		pylab.grid(True, alpha=0.6)
		"""
		pylab.hist(x_ls, 30, log=logHist)
		if xlim_in_1D:
			pylab.xlim(xlim_in_1D)
		
		pylab.xlabel(xlabel)
		pylab.ylabel('Count')
		#
		if logHist:
			fig_fname = addExtraLsToFilenamePrefix(output_fname, ['logHist'])
		else:
			fig_fname = output_fname
		pylab.savefig(fig_fname, dpi=300)
		
		C_ls = LD_ls
		fig_fname = addExtraLsToFilenamePrefix(output_fname, ['deletionFrequency_vs_snpDensity_2D'])
		cls.drawHexbin(x_ls, y_ls, C_ls, fig_fname=fig_fname, gridsize=20, \
				title=title, \
				xlabel="deletion frequency", \
				ylabel=ylabel_in_2D, \
				colorBarLabel='median(LD)')
		sys.stderr.write("Done.\n")
		
		C_ls = LD_ls
		fig_fname = addExtraLsToFilenamePrefix(output_fname, ['deletionFrequency_vs_deletionSize_2D'])
		cls.drawHexbin(x_ls, deletionSizeLs, C_ls, fig_fname=fig_fname, gridsize=20, \
				title=title, \
				xlabel="deletion frequency", \
				ylabel="log10(deletion size)", \
				colorBarLabel='median(LD)')
		sys.stderr.write("Done.\n")
	
	"""
		#2010-8-1 -z banyan.usc.edu
		input_fname = os.path.expanduser("~/cnvMethod20_vs_callMethod32_LD.tsv")
		output_dir = os.path.expanduser('~/doc/compbiophd/figures/')
		cnv_method_id=20
		CNV.drawCNVLDHist(db_250k, input_fname, output_dir, cnv_method_id=cnv_method_id, fileNamePrefix='CNVLD',\
							xlabel='LD', ylabel_in_2D='#SNPs per kb', xlim_in_1D=[0,1], logHist=False,)
		sys.exit(0)
	"""
	
	@classmethod
	def drawArrayQuartile(cls, db_250k, output_dir, call_method_id=None):
		"""
		2010-3-1
			draw the array quartile data from table ArrayQuartile
				1. 7-number summary of intensity data
				2. 7-number summary of dlrspread data
		"""
		import os, sys
		from common import get_ecotypeid2nativename
		ecotypeid2nativename = get_ecotypeid2nativename(db_250k.metadata.bind)
		
		sys.stderr.write("Visualizing data from ArrayQuartile call method %s ...\n"%call_method_id)
		import Stock_250kDB
		
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		
		# finding the minimum and maximum to be set as ylim for array quartile & dlr
		rows = db_250k.metadata.bind.execute("select min(minimum) as q_min, max(maximum) as q_max, min(dlr_minimum) as dlr_min, \
											max(dlr_maximum) as dlr_max from %s"%\
											Stock_250kDB.ArrayQuartile.table.name)
		row = rows.fetchone()
		quartile_ylim = [row.q_min, row.q_max]
		dlr_ylim = [row.dlr_min, row.dlr_max]
		
		call_info_ls = Stock_250kDB.CallInfo.query
		if call_method_id is not None:
			call_info_ls = call_info_ls.filter_by(method_id=call_method_id)
		
		field_name_ls = ['minimum', 'first_decile', 'lower_quartile', 'median', 'upper_quartile', 'last_decile', 'maximum']
		no_of_total_arrays = call_info_ls.count()
		count = 0
		for call_info in call_info_ls:
			count += 1
			sys.stderr.write("Array %s, Ecotype %s, %s/%s ..."%(call_info.array.id, call_info.array.maternal_ecotype_id, count, no_of_total_arrays))
			rows = Stock_250kDB.ArrayQuartile.query.filter_by(array_id=call_info.array.id).order_by(Stock_250kDB.ArrayQuartile.start_probe_id)
			seven_number_2D_ls = []
			no_of_outliers_ls = []
			dlr_2D_ls = []
			for row in rows:
				seven_number_2D_ls.append([])
				dlr_2D_ls.append([])
				#no_of_outliers_ls.append(len(row.array_quartile_outlier_ls))
				for field_name in field_name_ls:
					seven_number_2D_ls[-1].append(getattr(row, field_name))
					dlr_2D_ls[-1].append(getattr(row, 'dlr_%s'%field_name))
			ecotype_id = call_info.array.maternal_ecotype_id
			nativename = ecotypeid2nativename.get(ecotype_id)
			title = 'ecotype_%s_%s_array_%s'%(nativename, ecotype_id, call_info.array.id)
			output_fname_prefix = os.path.join(output_dir, '%s_array_quartile'%title)
			cls.drawOneArraySevenNumber(title, output_fname_prefix, seven_number_2D_ls, ylim=quartile_ylim)
			output_fname_prefix = os.path.join(output_dir, '%s_dlr'%title)
			cls.drawOneArraySevenNumber(title, output_fname_prefix, dlr_2D_ls, ylim=dlr_ylim)
			sys.stderr.write("\n")
		sys.stderr.write("Done.\n")
	
	"""
	output_dir = '/Network/Data/250k/tmp-yh/CNV/ArrayIntensitySummary'
	CNV.drawArrayQuartile(db_250k, output_dir, call_method_id=48)
	
	"""
	class LerContig(object):
		@classmethod
		def drawLerContigSizeHist(cls, fasta_input_fname, output_fname_prefix):
			"""
			2009-12-7
				This function draws a histogram of size of Ler contigs.			
					Ler contigs were downloaded from http://www.arabidopsis.org/browse/Cereon/index.jsp, in fasta format.
					a png file with prefix as output_fname_prefix will be generated.
			"""
			import os, sys
			sys.stderr.write("Drawing histogram of Ler contig size ... ")
			inf = open(fasta_input_fname)
			from Bio import SeqIO
			contig_size_ls = []
			for seq_record in SeqIO.parse(inf, "fasta") :
				contig_size_ls.append(len(seq_record.seq))
				#of.write('>%s\n'%seq_record.id)
			
			import pylab
			pylab.clf()
			title = 'Histogram of Ler Contig Size'
			pylab.title(title)
			pylab.hist(contig_size_ls, 30)
			pylab.xlabel('contig size')
			pylab.ylabel('Count')
			pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
			sys.stderr.write("Done.\n")
			
		"""
		fasta_input_fname = os.path.expanduser('~/script/variation/data/CNV/Cereon_Ath_Ler.fasta')
		output_fname_prefix = '/tmp/Ler-contig-size-hist'
		CNV.drawLerContigSizeHist(fasta_input_fname, output_fname_prefix)
		"""
		
		@classmethod
		def putLerContigsIntoDB(cls, db, data_source, accession_name, fasta_input_fname, addSeqIntoDB=True):
			"""
			2010-11-12
				add argument addSeqIntoDB:
					indicating whether raw sequences shall be added into db
			2010-1-27
				put ler contigs sequences and ids into table Stock_250kDB.SequenceFragment
			"""
			import os, sys
			sys.stderr.write("Putting Ler contigs into db ... ")
			session = db.session
			import Stock_250kDB
			from CNVQCConstruct import CNVQCConstruct
			data_source_obj = CNVQCConstruct.getDBObj(session, Stock_250kDB.DataSource, data_source)
			acc_obj = CNVQCConstruct.getCNVQCAccessionObjFromDB(session, accession_name, data_source_obj)
			inf = open(fasta_input_fname)
			from Bio import SeqIO
			for seq_record in SeqIO.parse(inf, "fasta"):
				sequence_fragment = Stock_250kDB.SequenceFragment(short_name=seq_record.id, size=len(seq_record.seq),\
																description=seq_record.description)
				if addSeqIntoDB:
					sequence_fragment.sequence = seq_record.seq.tostring()
				sequence_fragment.accession = acc_obj
				session.add(sequence_fragment)
			session.flush()
			sys.stderr.write("Done.\n")
		
		"""
		# 2010-1-27
		fasta_input_fname = os.path.expanduser('~/script/variation/data/CNV/Cereon_Ath_Ler.fasta')
		data_source = 'LerContig'
		accession_name = 'Ler-1'
		CNV.LerContig.putLerContigsIntoDB(db_250k, data_source, accession_name, fasta_input_fname)
		sys.exit(0)
		
		# 2010-11-12
		fasta_input_fname = os.path.expanduser('~/script/variation/data/lyrata/Araly1_assembly_scaffolds.fasta')
		data_source = 'TAIR9LyrataNormalMaxMatch'
		accession_name = 'Lyrata'
		CNV.LerContig.putLerContigsIntoDB(db_250k, data_source, accession_name, fasta_input_fname, addSeqIntoDB=False)
		sys.exit(0)
		
		# 2010-11-12
		fasta_input_fname = os.path.expanduser('~/script/variation/data/lyrata/Araly1_assembly_scaffolds.fasta')
		data_source = 'TAIR9LyrataNormalMaxMatch'
		accession_name = 'Lyrata'
		CNV.LerContig.putLerContigsIntoDB(db_250k, data_source, accession_name, fasta_input_fname, addSeqIntoDB=False)
		sys.exit(0)
		
		"""
		
		@classmethod
		def outputLerContigSpanDataFromDBInRunGADAOutputFormat(cls, db, output_fname):
			"""
			2010-2-10
			"""
			from RunGADA import RunGADA
			import csv
			writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
			RunGADA.output_header(writer)
			rows = db.metadata.bind.execute("select a.ecotype_id, sp.chromosome, sp.start, sp.stop, sp.stop-sp.start as length, sp.start_probe_id, sp.stop_probe_id \
					from sequence_fragment_ref_pos sp, sequence_fragment s, cnv_qc_accession a where sp.sequence_fragment_id=s.id and s.accession_id=a.id")
			for row in rows:
				start_probe = '%s_%s'%(row.chromosome, row.start+12)
				stop_probe = '%s_%s'%(row.chromosome, row.stop-12)
				new_row = [row.ecotype_id, 3, start_probe, stop_probe, row.stop-row.start, \
					-2, row.start_probe_id, row.stop_probe_id]
				writer.writerow(new_row)
			del writer
		
		"""
		#2010-2-10
		output_fname = '/tmp/LerContigSpan.tsv'
		CNV.outputLerContigSpanDataFromDBInRunGADAOutputFormat(db_250k, output_fname)
		"""
		
		@classmethod
		def discoverLerDeletionWithinOneContig(cls, start_index, stop_index, probe_pos2fragment_start, chr_pos_ls,
											total_probe_id_ls, probes, contig_id=None):
			"""
			2010-4-29
				split out of discoverLerDeletionDuplication(). now called by it.
				call deletion within one contig.
					order the matching probes within one contig in chromosomal order.
					any segment that has no matching probes but surrounded by two matching probes
			"""
			deletion_ls = []
			index_of_prev_probe_within_a_contig = None
			prev_copy_number = 0
			for i in range(start_index, stop_index+1):
				probe_pos = chr_pos_ls[i]
				copy_number = 0
				if probe_pos in probe_pos2fragment_start:
					copy_number = 1
				if i>start_index:	#from the 2nd probe
					previous_probe_pos = chr_pos_ls[i-1]
					
					prev_copy_number = 0
					if previous_probe_pos in probe_pos2fragment_start:
						prev_copy_number = 1
					
					if prev_copy_number>0 and copy_number==0:
						index_of_prev_probe_within_a_contig = i-1
					elif prev_copy_number==0 and copy_number>0:	# from deletion to non-deletion
						if index_of_prev_probe_within_a_contig is not None:	# found a potential deletion
							prev_non_deleted_probe_pos = chr_pos_ls[index_of_prev_probe_within_a_contig]
							prev_fragment_start = probe_pos2fragment_start.get(prev_non_deleted_probe_pos)
							current_fragment_start = probe_pos2fragment_start.get(probe_pos)
							size_difference = abs(probe_pos[1] - prev_non_deleted_probe_pos[1]) - abs(current_fragment_start-prev_fragment_start)
							
							deletion_start_probe_id = total_probe_id_ls[index_of_prev_probe_within_a_contig+1]
							deletion_start_probe = probes.get_one_probe(deletion_start_probe_id)
							deletion_start_chr_pos = '%s_%s'%(deletion_start_probe.chr, deletion_start_probe.pos)
							deletion_stop_probe_id = total_probe_id_ls[i-1]
							deletion_stop_probe = probes.get_one_probe(deletion_stop_probe_id)
							deletion_stop_chr_pos = '%s_%s'%(deletion_stop_probe.chr, deletion_stop_probe.pos)
							row = [deletion_start_probe_id, deletion_start_chr_pos, deletion_stop_probe_id, \
									deletion_stop_chr_pos, \
									i-index_of_prev_probe_within_a_contig-1, \
									deletion_stop_probe.pos-deletion_start_probe.pos, \
									prev_copy_number, contig_id, size_difference]
							deletion_ls.append(row)
						index_of_prev_probe_within_a_contig = i
					elif prev_copy_number>0 and copy_number>0:	# from non-deletion to non-deletion
						index_of_prev_probe_within_a_contig = i
			return deletion_ls
		
		@classmethod
		def findNearbyProbeWithFragmentStart(cls, probe_index, probe_pos2fragment_start, chr_pos_ls, step=0, maxDepth=10):
			"""
			2010-6-9
				pre-stop the loop if step is beyond maxDepth
			2010-4-29
				find the nearest probe which has fragment start
			"""
			probe_with_fragment_start = None
			if step>maxDepth:	#2010-6-9 pre-stop it
				return None, None, None
			negative_step = - step
			if negative_step == step:
				step_ls = [step]
			else:
				step_ls = [negative_step, step]
			for _step in step_ls:
				probe_pos = chr_pos_ls[probe_index+_step]
				if probe_pos in probe_pos2fragment_start:
					fragment_start = probe_pos2fragment_start.get(probe_pos)[0]
					return probe_pos, fragment_start, _step
			# increase the step by 1 if it's not successful
			return cls.findNearbyProbeWithFragmentStart(probe_index, probe_pos2fragment_start, chr_pos_ls, step=step+1)
			
		@classmethod
		def discoverLerDeletionWithinOneContigByGADA(cls, start_index, stop_index, probe_pos2fragment_start, chr_pos_ls,\
											total_probe_id_ls, probes, contig_id=None, \
											tmp_fname_prefix='/tmp/GADA',\
											IO_thru_file=False,\
											GADA_path=os.path.expanduser('~/script/variation/bin/GADA/GADA'), \
											aAlpha=0.5, TBackElim=2, MinSegLen=5, debug=False,
											probe_pos2contig_id_ls=None):
			"""
			2010-4-29
				called by discoverLerDeletionDuplication()
			"""
			tmp_input_fname = '%s_%s_input'%(tmp_fname_prefix, contig_id.replace(' ', '_'))
			tmp_output_fname = '%s_%s_output'%(tmp_fname_prefix, contig_id.replace(' ', '_'))
			probe_copy_number_ls = []
			for i in range(start_index, stop_index+1):
				probe_pos = chr_pos_ls[i]
				copy_number = 0
				"""
				if probe_pos in probe_pos2fragment_start:
					copy_number = probe_pos2fragment_start.get(probe_pos)[1]
				"""
				if probe_pos in probe_pos2contig_id_ls:
					copy_number = len(probe_pos2contig_id_ls.get(probe_pos))
				probe_copy_number_ls.append(float(copy_number))
			probe_copy_number_set = set(probe_copy_number_ls)
			if len(probe_copy_number_set)<=1:
				# 2010-04-30 only one type of copy number in the data. GADA will fail to run if no variance.
				# just return empty.
				return []
			if len(probe_copy_number_ls)<MinSegLen:	#2010-6-7 shorter than even the MinSegLen. ignore
				return []
			
			if IO_thru_file:
				from RunGADA import RunGADA
				input_file = RunGADA.prepareGADAinput(probe_copy_number_ls, tmp_input_fname, IO_thru_file=IO_thru_file)
				GADA_output = RunGADA._GADA(GADA_path, input_file, tmp_output_fname, aAlpha, \
										TBackElim, MinSegLen, IO_thru_file=IO_thru_file)
			else:
				import GADA	#2010-4-29 moved here to avoid import errors if this module is not compiled well
				ins = GADA.GADA(debug)
				GADA_output = ins.run(map(float, probe_copy_number_ls), aAlpha, TBackElim, MinSegLen)
			
			if type(GADA_output)==str:
				import cStringIO, csv
				GADA_output = cStringIO.StringIO(GADA_output)
				GADA_output = csv.reader(GADA_output, delimiter='\t')
			counter = 0
			deletion_ls = []
			for row in GADA_output:
				if type(row[0])==str and (row[0].find('#')==0 or row[0].find('Start')==0):	# ignore the comments
					#skip the first line of real data. it's header "Start   Stop    Lenght  Ampl    State"
					continue
				counter += 1
				probe1_index, probe2_index, no_of_probes, amplitude = row[:4]
				no_of_probes = int(no_of_probes)
				amplitude = float(amplitude)
				probe1_index = int(probe1_index) + start_index - 1
				probe2_index = int(probe2_index) + start_index - 1
				probe1 = chr_pos_ls[probe1_index]
				probe1_chr_pos = '%s_%s'%(probe1[0], probe1[1])
				probe2 = chr_pos_ls[probe2_index]
				probe2_chr_pos = '%s_%s'%(probe2[0], probe2[1])
				
				start_probe_with_fragment_pos, fragment_start, start_step = cls.findNearbyProbeWithFragmentStart(probe1_index,\
																probe_pos2fragment_start, chr_pos_ls, step=0)
				
				stop_probe_with_fragment_pos, fragment_stop, stop_step = cls.findNearbyProbeWithFragmentStart(probe2_index,\
																probe_pos2fragment_start, chr_pos_ls, step=0)
				if start_step is not None and stop_step is not None:
					start_probe_with_fragment_index = probe1_index + start_step
					stop_probe_with_fragment_index = probe2_index + stop_step
					step_tuple = (start_step, stop_step)
					if start_probe_with_fragment_index<=stop_probe_with_fragment_index and \
						start_probe_with_fragment_index<=probe2_index and \
						stop_probe_with_fragment_index>=probe1_index:
						
						size_difference = abs(stop_probe_with_fragment_pos[1]-start_probe_with_fragment_pos[1]) - abs(fragment_stop-fragment_start)
					else:
						size_difference = None
				else:
					step_tuple = (None, None)
					size_difference = None
				probe1_id = total_probe_id_ls[probe1_index]
				probe2_id = total_probe_id_ls[probe2_index]
				new_row = [probe1_id, probe1_chr_pos, probe2_id, probe2_chr_pos, \
						no_of_probes, probe2[1]-probe1[1], \
						amplitude, contig_id, size_difference, step_tuple]
				deletion_ls.append(new_row)
			return deletion_ls
		
		@classmethod
		def addSegmentIntoRefCoverTupleLs(cls, segment, refCoverTuple_ls):
			"""
			2010-6-8
			"""
			probe1_chr_pos = segment[1]
			probe1_chr_pos = probe1_chr_pos.split('_')
			probe1_chr_pos = map(int, probe1_chr_pos)
			probe2_chr_pos = segment[3]
			probe2_chr_pos = probe2_chr_pos.split('_')
			probe2_chr_pos = map(int, probe2_chr_pos)
			
			refCoverTuple_ls.append([probe1_chr_pos[0], probe1_chr_pos[1], probe2_chr_pos[1], segment])
		
		@classmethod
		def mergeTwoSegments(cls, segment, next_segment):
			"""
			2010-6-8
				segment = ['start_probe_id', 'start_chr_pos', 'stop_probe_id', 'stop_chr_pos', 'no_of_probes', \
						'length', 'copy_number', 'contig', 'size_difference']
			"""
			new_contig_id = '%s,%s'%(segment[7], next_segment[7])
			if segment[8] is not None and next_segment[8] is not None:
				size_difference = segment[8]+next_segment[8]
			else:
				size_difference = None
			import numpy
			new_segment = [segment[0], segment[1], next_segment[2], next_segment[3], segment[4]+next_segment[4],\
						segment[5]+next_segment[5], numpy.mean([segment[6],next_segment[6]]), new_contig_id,\
						size_difference,]
			return new_segment
			
		@classmethod
		def discoverLerDeletionDuplication(cls, db_250k, ler_blast_result_fname, output_fname_prefix, deletion_only=True, \
										min_no_of_matches=25, max_deletion_copy_amp=0.1, min_coverage_copy_amp=0.5,\
										min_dup_copy_amp=1.9, maxDeletionLength=50000, IO_thru_file=False, \
										tmp_fname_prefix='/tmp/GADA'):
			"""
			2010-6-9
				ToDo: A portion of the final deletions  are overlapping or too close and need to be merged.
				But this function is abandoned in favor of discoverLerDeletionAndSpanFromMummerCoordsOutput().
					Because the probe-blast scheme can only focus on 25-exact-matches (computationally too intensive to
					get down the number of matches). 
			2010-6-7
				max_deletion_copy_amp: deleltion's GADA amp has to be smaller than this.
				min_coverage_copy_amp: any ref segment which has GADA amp above this is counted as covered.
				min_dup_copy_amp: ref segment's GADA amp has to be above this to be called duplication. (doesn't care 2 or more copies)
				maxDeletionLength: discovered-deletion has to be shorter than this.
				
				A bona-fide deletion will be called only when
					1. its amp is max_deletion_copy_amp.
					2. its length <= maxDeletionLength
					3. flanked by two segments whose amp are above min_coverage_copy_amp
				
				A bona-fide deletion is also counted as a covered segment apart from segments whose amp are above min_coverage_copy_amp.
			2010-4-28
				output contig ids in deletion_only=True mode
			2009-12-7
				ler_blast_result_fname is the output of blasting all CNV probes against Ler contigs
					http://www.arabidopsis.org/browse/Cereon/index.jsp.
				Two functions:
					1. deletion_only=True. make sure the deletion is covered by the sequencing.
						one naive criteria is if the boundary (two adjacent non-deleted probes) is within the same contig, then yes.
					2. deletion_only=False, detect copy number changes. If two adjacent probes have different number of contigs, 
						then it's a copy number change point.
			"""
			from pymodule import getColName2IndexFromHeader, PassingData, figureOutDelimiter
			sys.stderr.write("Reading from %s ... \n"%ler_blast_result_fname)
			counter = 0
			real_counter = 0
			import csv
			reader = csv.reader(open(ler_blast_result_fname), delimiter=figureOutDelimiter(ler_blast_result_fname))
			header = reader.next()
			col_name2index = getColName2IndexFromHeader(header)
			probe_id2contig_id_ls = {}
			contig_id2probe_pos_ls = {}
			for row in reader:
				chr = int(row[col_name2index['Chromosome']])
				pos = int(row[col_name2index['Position']])
				probe_pos = (chr, pos)
				contig_label = row[col_name2index['Alignment_title']]	#original looks like "gnl|BL_ORD_ID|58400 ATL8C37289 ATL7C108502_1"
				probe_id = int(row[col_name2index['Probe_ID']])
				no_of_matches = int(row[col_name2index['Number_matches']])
				alignment_start = int(row[col_name2index['Alignment_start']])
				query_start = int(row[col_name2index['query_start']])
				if no_of_matches>=min_no_of_matches:
					contig_id = ' '.join(contig_label.split()[1:])
					if contig_id not in contig_id2probe_pos_ls:
						contig_id2probe_pos_ls[contig_id] = []
					contig_id2probe_pos_ls[contig_id].append((probe_pos, alignment_start))
					probe_key = probe_pos	# 2010-6-7 previously it was probe_id
					if probe_key not in probe_id2contig_id_ls:
						probe_id2contig_id_ls[probe_key] = []
					probe_id2contig_id_ls[probe_key].append(contig_id)
			sys.stderr.write("%s probes in %s contigs Done.\n"%(len(probe_id2contig_id_ls), len(contig_id2probe_pos_ls)))
			del reader
			
			import Stock_250kDB
			from DB_250k2Array import DB_250k2Array
			probes, xy_ls, chr_pos_ls, total_probe_id_ls = DB_250k2Array.get_probes(db_250k.metadata.bind, \
										Stock_250kDB.Probes.table.name, snps=None, run_type=2)
			
			returnData = DB_250k2Array.organizeProbesIntoChromosome(xy_ls, chr_pos_ls, total_probe_id_ls)
			chr2xy_ls = returnData.chr2xy_ls
			chr2probe_id_ls = returnData.chr2probe_id_ls
			chr_pos2index = returnData.chr_pos2index
			
			delDupWriter = csv.writer(open('%s_delDup.tsv'%output_fname_prefix, 'w'), delimiter='\t')
			refCoveredWriter =  csv.writer(open('%s_ref_covered.tsv'%output_fname_prefix, 'w'), delimiter='\t')
			
			header_row = ['start_probe_id', 'start_chr_pos', 'stop_probe_id', 'stop_chr_pos', 'no_of_probes', \
						'length', 'copy_number', 'contig', 'size_difference']
			for writer in [delDupWriter, refCoveredWriter]:
				writer.writerow(header_row)
			
			sys.stderr.write("Discovering deletions ...\n")
			counter = 0
			real_counter = 0
			MinSegLen = 5
			refCoverTuple_ls = []
			if deletion_only==True:	#2010-4-29 new way of detecting deletion covered by Ler contigs
				for contig_id, probe_pos_ls in contig_id2probe_pos_ls.iteritems():
					probe_pos_ls.sort()
					probe_pos2fragment_start = {}
					chr2start_stop_probe_ls = {}
					for probe_pos, fragment_start in probe_pos_ls:
						chr, pos = probe_pos
						if chr not in chr2start_stop_probe_ls:
							chr2start_stop_probe_ls[chr] = []
						chr2start_stop_probe_ls[chr].append(probe_pos)
						probe_pos2fragment_start[probe_pos] = [fragment_start, len(probe_id2contig_id_ls[probe_pos])]
					for chr, start_stop_probe_ls in chr2start_stop_probe_ls.iteritems():
						if len(start_stop_probe_ls)<2:
							continue
						start_probe_pos = start_stop_probe_ls[0]
						stop_probe_pos = start_stop_probe_ls[-1]
						start_index = chr_pos2index.get(start_probe_pos)
						stop_index = chr_pos2index.get(stop_probe_pos)
						if start_index is not None and stop_index is not None:
							"""
							deletion_ls = cls.discoverLerDeletionWithinOneContig(start_index, stop_index, \
											probe_pos2fragment_start, chr_pos_ls,
											total_probe_id_ls, probes, contig_id=contig_id)
							"""
							
							segment_ls = cls.discoverLerDeletionWithinOneContigByGADA(start_index, stop_index, \
											probe_pos2fragment_start, chr_pos_ls,
											total_probe_id_ls, probes, contig_id=contig_id, \
											tmp_fname_prefix='%s_%s_%s'%(tmp_fname_prefix, contig_id.replace(' ', '_'), chr),\
											IO_thru_file=IO_thru_file, \
											aAlpha=0.5, TBackElim=2, MinSegLen=MinSegLen, \
											probe_pos2contig_id_ls=probe_id2contig_id_ls)
							
							no_of_segments = len(segment_ls)
							real_counter += no_of_segments
							ampIndex = 6
							for i in xrange(no_of_segments):
								segment = segment_ls[i]
								amp = segment[ampIndex]
								segmentLen = segment[5]
								if amp>=min_coverage_copy_amp:
									cls.addSegmentIntoRefCoverTupleLs(segment, refCoverTuple_ls)
									#refCoverTuple_ls.append(segment)
									if amp>=min_dup_copy_amp:
										segment[ampIndex] = 2
										delDupWriter.writerow(segment)
								elif amp<=max_deletion_copy_amp and segmentLen<=maxDeletionLength:
									if i>0 and i<no_of_segments-1:	# not in the beginning and not in the end
										prevSegmentAmp = segment_ls[i-1][ampIndex]
										nextSegmentAmp = segment_ls[i+1][ampIndex]
										if prevSegmentAmp>=min_coverage_copy_amp and nextSegmentAmp>=min_coverage_copy_amp:
											cls.addSegmentIntoRefCoverTupleLs(segment, refCoverTuple_ls)
											
											segment[ampIndex] = 0	#modify the copy number variable to 0
											delDupWriter.writerow(segment)
					counter += 1
					if counter%1000==0:
						sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
			del delDupWriter
			sys.stderr.write("Done.\n")
			
			sys.stderr.write("Outputting reference coverage ...")
			refCoverTuple_ls.sort()
			
			refCoverTuple = refCoverTuple_ls.pop(0)	# start from beginning
			while refCoverTuple_ls:
				refCoverTuple_chr = refCoverTuple[0]
				refCoverTuple_end = refCoverTuple[2]
				refCoverTuple_segment = refCoverTuple[3]
				
				next_refCoverTuple = refCoverTuple_ls.pop(0)	# start from beginning
				next_refCoverTuple_chr = next_refCoverTuple[0]
				next_refCoverTuple_start = next_refCoverTuple[1]
				next_refCoverTuple_end = next_refCoverTuple[2]
				next_refCoverTuple_segment = next_refCoverTuple[3]
				if next_refCoverTuple_chr==refCoverTuple_chr and next_refCoverTuple_start<=refCoverTuple_end:
					#merge the two
					refCoverTuple[2] = next_refCoverTuple_end
					refCoverTuple[3] = cls.mergeTwoSegments(refCoverTuple_segment, next_refCoverTuple_segment)
				else:
					
					# output it and replace the current refCoverTuple
					refCoveredWriter.writerow(refCoverTuple[3])
					
					refCoverTuple = next_refCoverTuple
			# don't forget to output the last refCoverTuple
			refCoveredWriter.writerow(refCoverTuple[3])
			del refCoveredWriter
			sys.stderr.write("Done.\n")
			return
			
			# 2010-4-29 old code
			for chr, probe_id_ls in chr2probe_id_ls.iteritems():
				no_of_probes = len(probe_id_ls)
				index_of_prev_probe_within_a_contig = None
				index_of_prev_probe_with_a_different_copy_number = None
				for i in range(no_of_probes):
					probe_id = probe_id_ls[i]
					contig_id_ls = probe_id2contig_id_ls.get(probe_id,[])
					copy_number = len(contig_id_ls)
					if i==0:
						index_of_prev_probe_with_a_different_copy_number = -1	# set before the first probe (index=0)
					else:
						prev_probe_contig_id_ls = probe_id2contig_id_ls.get(probe_id_ls[i-1],[])
						prev_copy_number = len(prev_probe_contig_id_ls)
						if not deletion_only:
							if copy_number != prev_copy_number:	# a change point of copy number
								if index_of_prev_probe_with_a_different_copy_number is not None:
									start_probe_id = probe_id_ls[index_of_prev_probe_with_a_different_copy_number+1]
									start_probe = probes.get_one_probe(start_probe_id)
									start_chr_pos = '%s_%s'%(start_probe.chr, start_probe.pos)
									stop_probe_id = probe_id_ls[i-1]
									stop_probe = probes.get_one_probe(stop_probe_id)
									stop_chr_pos = '%s_%s'%(stop_probe.chr, stop_probe.pos)
									row = [start_probe_id, start_chr_pos, stop_probe_id, stop_chr_pos, \
											i-index_of_prev_probe_with_a_different_copy_number-1, \
											stop_probe.pos-start_probe.pos, prev_copy_number]
									writer.writerow(row)
									real_counter += 1
								index_of_prev_probe_with_a_different_copy_number = i-1
						else:	# look for deleton only. The only difference from above is make sure the deletion is covered by the sequencing.
							#one naive criteria is if the boundary is within the same contig, then yes.
							if prev_copy_number>0 and copy_number==0:	# from non-deletion to deletion
								index_of_prev_probe_within_a_contig = i-1
							elif prev_copy_number==0 and copy_number>0:	# from deletion to non-deletion
								if index_of_prev_probe_within_a_contig is not None:	# found a potential deletion
									current_contig_id_set = set(contig_id_ls)
									prev_non_deleted_probe_id = probe_id_ls[index_of_prev_probe_within_a_contig]
									prev_non_deleted_probe_contig_id_ls = probe_id2contig_id_ls.get(prev_non_deleted_probe_id, [])
									prev_non_deleted_probe_contig_id_set = set(prev_non_deleted_probe_contig_id_ls)
									shared_contig_id_ls = list(prev_non_deleted_probe_contig_id_set&current_contig_id_set)
									if len(shared_contig_id_ls)>0:	#share at least one contig. deletion confirmed
										deletion_start_probe_id = probe_id_ls[index_of_prev_probe_within_a_contig+1]
										deletion_start_probe = probes.get_one_probe(deletion_start_probe_id)
										deletion_start_chr_pos = '%s_%s'%(deletion_start_probe.chr, deletion_start_probe.pos)
										deletion_stop_probe_id = probe_id_ls[i-1]
										deletion_stop_probe = probes.get_one_probe(deletion_stop_probe_id)
										deletion_stop_chr_pos = '%s_%s'%(deletion_stop_probe.chr, deletion_stop_probe.pos)
										row = [deletion_start_probe_id, deletion_start_chr_pos, deletion_stop_probe_id, \
											deletion_stop_chr_pos, \
											i-index_of_prev_probe_within_a_contig-1, deletion_stop_probe.pos-deletion_start_probe.pos, \
											prev_copy_number, ','.join(shared_contig_id_ls)]
										writer.writerow(row)
										real_counter += 1
								index_of_prev_probe_within_a_contig = i
							elif prev_copy_number>0 and copy_number>0:	# from non-deletion to non-deletion
								index_of_prev_probe_within_a_contig = i
					
					counter += 1
					if counter%10000==0:
						sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
				
				# don't forget the last segment if it's not in the deletion_only mode.
				if not deletion_only and index_of_prev_probe_with_a_different_copy_number is not None:
					start_probe_id = probe_id_ls[index_of_prev_probe_with_a_different_copy_number+1]
					start_probe = probes.get_one_probe(start_probe_id)
					start_chr_pos = '%s_%s'%(start_probe.chr, start_probe.pos)
					stop_probe_id = probe_id_ls[i]	# watch: not i-1.
					stop_probe = probes.get_one_probe(stop_probe_id)
					stop_chr_pos = '%s_%s'%(stop_probe.chr, stop_probe.pos)
					row = [start_probe_id, start_chr_pos, stop_probe_id, stop_chr_pos, \
							i-index_of_prev_probe_with_a_different_copy_number, stop_probe.pos-start_probe.pos, copy_number]	# watch no -1, and it's copy_number
					writer.writerow(row)
					real_counter += 1
			sys.stderr.write("Done.\n")
			del writer
		"""
	ler_blast_result_fname = '/Network/Data/250k/tmp-dazhe/ler_raw_CNV_QC.csv'
	ler_blast_result_fname = os.path.expanduser('~/mnt/hpc-cmb/script/variation/data/CNV/cnv_probe_blast_against_ler_contig.tsv')
	output_fname = '/tmp/Ler-deletions.tsv'
	CNV.LerContig.discoverLerDeletionDuplication(db_250k, ler_blast_result_fname, output_fname, deletion_only=True)
	
	ler_blast_result_fname = '/Network/Data/250k/tmp-dazhe/ler_raw_CNV_QC.csv'
	output_fname = '/tmp/Ler-copy-number.tsv'
	CNV.LerContig.discoverLerDeletionDuplication(db_250k, ler_blast_result_fname, output_fname, deletion_only=False)
	
	ler_blast_result_fname = '/Network/Data/250k/tmp-dazhe/tair9_raw.csv'
	output_fname = '/tmp/Col-copy-number.tsv'
	CNV.LerContig.discoverLerDeletionDuplication(db_250k, ler_blast_result_fname, output_fname, deletion_only=False)
	
	ler_blast_result_fname = os.path.expanduser('~/mnt/hpc-cmb/script/variation/data/CNV/cnv_probe_blast_against_ler_contig.tsv')
	output_fname = '/tmp/Ler-deletions.tsv'
	CNV.LerContig.discoverLerDeletionDuplication(db_250k, ler_blast_result_fname, output_fname, deletion_only=True)
	sys.exit(0)
	
	max_deletion_copy_amp = 0.1
	ler_blast_result_fname = os.path.expanduser('~/mnt/hpc-cmb/script/variation/data/CNV/cnv_probe_blast_against_ler_contig.tsv')
	output_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/CNVProbeBlastAgainstLerContig_maxDelAmp%s_bugfix'%(max_deletion_copy_amp))
	CNV.LerContig.discoverLerDeletionDuplication(db_250k, ler_blast_result_fname, output_fname_prefix, \
												deletion_only=True, max_deletion_copy_amp=max_deletion_copy_amp, \
												IO_thru_file=False)
	sys.exit(0)
	
	
		"""
		@classmethod
		def drawLerDelDupByGADAOutputCopyNumberHist(cls, discoverLerDelDupByGADA_output_fname, output_fname, max_amplitude=None,\
												drawKernelDensity=False):
			"""
			2010-6-7
				It draws a histogram of "copy-number" (amplitude) of segments outputted by discoverLerDeletionDuplication (using GADA).
				
				Purpose of this function is to see whether there's a natural cutoff lying in a bi-modal distribution.
				
				argument max_amplitude is used to filter out segments whose amplitude is above this.
					purpose is to avoid saturation of histogram by normal fragments.
			"""
			import csv
			input_fname = discoverLerDelDupByGADA_output_fname
			from pymodule import figureOutDelimiter, getColName2IndexFromHeader
			reader = csv.reader(open(input_fname), delimiter=figureOutDelimiter(input_fname))
			header = reader.next()
			col_name2index = getColName2IndexFromHeader(header)
			amplitude_ls = []
			for row in reader:
				copy_number_col_index = col_name2index.get('copy_number')
				amplitude = float(row[copy_number_col_index])
				if max_amplitude is not None and amplitude>max_amplitude:
					continue
				amplitude_ls.append(amplitude)
			
			import pylab
			pylab.clf()
			title = 'Histogram of copy number amplitude (GADA)'
			pylab.title(title)
			if drawKernelDensity:
				import statistics	# 2010-5-30 package from Michiel De Hoon
				y, x = statistics.pdf(amplitude_ls)
				pylab.plot(x, y, alpha=0.7)
			else:
				pylab.hist(amplitude_ls, 35)
			pylab.ylabel('density')
			pylab.xlabel('LerContig probe copy number amplitude')
			pylab.savefig(output_fname, dpi=300)
			sys.stderr.write("Done.\n")
		
		"""
	
	max_amplitude = 3	# 0.6
	drawKernelDensity = True
	common_prefix = '~/script/variation/data/CNV/Ler-deletions-within-one-contig-by-GADAA0.5T2M5'
	common_prefix = '~/script/variation/data/CNV/Ler-deletions-within-1-contig-dup-by-GADAA0.5T2M5'
	discoverLerDelDupByGADA_output_fname = os.path.expanduser('%s.tsv'%common_prefix)
	if drawKernelDensity is False:
		common_prefix += '_hist'
	if max_amplitude is not None:
		common_prefix += '_maxAmp%.1f'%max_amplitude
	output_fname = os.path.expanduser('%s.png'%(common_prefix))
	CNV.LerContig.drawLerDelDupByGADAOutputCopyNumberHist(discoverLerDelDupByGADA_output_fname, output_fname, \
										max_amplitude=max_amplitude, drawKernelDensity=drawKernelDensity)
	sys.exit(0)
	
		"""
		
		@classmethod
		def discoverLerContigSpanOverCol(cls, ler_blast_result_fname, output_fname, min_no_of_matches=25,
										max_delta_ratio=0.4, max_length_delta=50000):
			"""
			2010-1-27
				discover the span of Ler Contigs in terms of Col reference genome
				
				The algorithm:
					for each contig, get all the probes that perfectly match some part of it.
						sort the probes in chromosomal order
						set the start-probe to be the 1st probe of all probes matched to the contig
							for stop-probe in [last probe, ..., 2nd probe]:	# in the reverse order
								if [start-probe, stop-probe] and [alignment-start, alignment-stop] meet the condition:
									col-span of this Ler contig = the chromosomal positions of [start-probe, stop-probe]
				
				The condition is either the length delta is <= max_length_delta, or either of the two delta ratios
					(delta/length1, delta/length2) <= max_delta_ratio.
				
			"""
			from pymodule import getColName2IndexFromHeader, PassingData, figureOutDelimiter
			sys.stderr.write("Reading from %s ... \n"%ler_blast_result_fname)
			counter = 0
			real_counter = 0
			import csv
			reader = csv.reader(open(ler_blast_result_fname), delimiter=figureOutDelimiter(ler_blast_result_fname))
			header = reader.next()
			col_name2index = getColName2IndexFromHeader(header)
			contig_id2probe_chr_pos_ls = {}
			for row in reader:
				contig_label = row[col_name2index['Alignment_title']]
				probe_id = int(row[col_name2index['Probe_ID']])
				no_of_matches = int(row[col_name2index['Number_matches']])
				if no_of_matches>=min_no_of_matches:
					contig_id = contig_label.split()[1]	# 2010-1-28, take "ATL8C9990" in the case of "ATL8C9990 ATL7C121_1"
					if contig_id not in contig_id2probe_chr_pos_ls:
						contig_id2probe_chr_pos_ls[contig_id] = []
					chr = int(row[col_name2index['Chromosome']])
					pos = int(row[col_name2index['Position']])
					alignment_start = int(row[col_name2index['Alignment_start']])
					probe_chr_pos = (chr, pos)
					probe_id = int(probe_id)
					contig_id2probe_chr_pos_ls[contig_id].append((probe_chr_pos, probe_id, alignment_start))
			del reader
			sys.stderr.write("Done.\n")
			
			sys.stderr.write("Discovering Ler Contig Span over Col ...\n")		
			import csv, Stock_250kDB
			writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
			header_row = ['start_probe_id', 'start_chr_pos', 'stop_probe_id', 'stop_chr_pos', 'no_of_probes', 'length', 'copy_number', \
						'contig_id', 'contig_start', 'contig_stop', 'size_difference']
			writer.writerow(header_row)
			counter = 0
			for contig_id, probe_chr_pos_ls in contig_id2probe_chr_pos_ls.iteritems():
				sys.stderr.write('%sContig %s, %s'%('\x08'*80, counter, contig_id))
				probe_chr_pos_ls.sort()	# sorted according to probe position on Col genome
				n = len(probe_chr_pos_ls)
				span_start_probe_index = 0
				while span_start_probe_index<n:
					old_span_start_probe_index = span_start_probe_index
					stop_index_candidate_ls = range(span_start_probe_index+1, n)
					stop_index_candidate_ls.reverse()	# starting from the end to cover as much as it can
					for i in stop_index_candidate_ls:
						start_probe_chr_pos, start_probe_id, start_alignment_pos = probe_chr_pos_ls[span_start_probe_index]
						stop_probe_chr_pos, stop_probe_id, stop_alignment_pos = probe_chr_pos_ls[i]
						start_chr, start_pos = start_probe_chr_pos
						stop_chr, stop_pos = stop_probe_chr_pos
						if start_chr == stop_chr:
							col_length = abs(stop_pos-start_pos)
							ler_length = abs(stop_alignment_pos-start_alignment_pos)
							length_delta = col_length-ler_length
							if col_length>0:
								col_ratio = abs(length_delta)/float(col_length)
							else:
								continue	# 2010-1-28 ignore if col_length = 0
							if ler_length>0:
								ler_ratio = abs(length_delta)/float(ler_length)
							else:
								continue	# 2010-1-28 ignore if ler_length = 0
							if abs(length_delta)<=max_length_delta or col_ratio<=max_delta_ratio or ler_ratio<=max_delta_ratio:
								no_of_probes = i-span_start_probe_index+1
								row = [start_probe_id, '%s_%s'%(start_chr, start_pos-12), stop_probe_id, '%s_%s'%(stop_chr, stop_pos+12), no_of_probes, stop_pos-start_pos+25, \
									'', contig_id, start_alignment_pos, stop_alignment_pos, length_delta]
								writer.writerow(row)
								span_start_probe_index = i+1	# set the new starting index
								break
					if old_span_start_probe_index == span_start_probe_index:	# no change in the for loop above. or nothing found to meet the min_reciprocal_overlap
						span_start_probe_index += 1	# auto increment it by 1
				counter += 1
			del writer
			sys.stderr.write("Done.\n")
		
		"""
		ler_blast_result_fname = '/Network/Data/250k/tmp-dazhe/ler_raw_CNV_QC.csv'
		max_delta_ratio = 0.4
		for max_length_delta in range(1,7):		# range(1,11)+[12,15,20]:
			max_length_delta = max_length_delta*10000
			output_fname = '/tmp/Ler-span-over-Col-mdr%s-mld%s.tsv'%(max_delta_ratio, max_length_delta)
			CNV.discoverLerContigSpanOverCol(ler_blast_result_fname, output_fname, min_no_of_matches=25,
											max_delta_ratio=max_delta_ratio, max_length_delta=max_length_delta)
											
		"""
		
		@classmethod
		def mergeTwoMummerCoords(cls, ref_pos, next_ref_pos):
			"""
			2010-6-10
				ref_pos = [chr, ref_start, ref_stop, input_data_index, contig_id_set, query_start, query_stop, identity_perc]
				
				used by discoverLerDeletionAndSpanFromMummerCoordsOutput()
			"""
			len1 = abs(ref_pos[2] - ref_pos[1] + 1)
			len2 = abs(next_ref_pos[2] - next_ref_pos[1] + 1)
			
			#take the maximum of two stops
			ref_pos[2] = max(ref_pos[2], next_ref_pos[2])
			
			
			# combine two contig_id_set
			ref_pos[4] = ref_pos[4]|next_ref_pos[4]
			
			for i in [5,6]:	#reset the query_start, query_stop
				ref_pos[i] = "NA"
			
			# recalculate the identity_perc
			ref_pos[7] = (ref_pos[7]*len1 + next_ref_pos[7]*len2)/float(len1+len2)
			
			return ref_pos
		
		@classmethod
		def getGapBetweenTwoMummerCoords(cls, ref_pos, next_ref_pos):
			"""
			2010-6-10
				ref_pos = [chr, ref_start, ref_stop, input_data_index, contig_id_set, query_start, query_stop, identity_perc]
				
				use by discoverLerDeletionAndSpanFromMummerCoordsOutput()
				
			"""
			chr = ref_pos[0]
			deletion_start = ref_pos[2]+1
			deletion_stop = next_ref_pos[1]-1
			deletion_length = deletion_stop - deletion_start + 1
			#combine two contig_id_set
			contig_id_set = ref_pos[4]|next_ref_pos[4]
			return [chr, deletion_start, deletion_stop, 'NA', contig_id_set, 'NA', "NA", 0]
		
		@classmethod
		def calGapRatioOfTwoSegments(cls, segment, next_segment):
			"""
			2010-6-11
			"""
			current_chr = segment[0]
			current_start = segment[1]
			current_end = segment[2]
				
			next_chr = next_segment[0]
			next_start = next_segment[1]
			next_end = next_segment[2]
			
			gap_len = None
			mergedLength = None
			if next_chr==current_chr:
				gap_len = next_start - current_end - 1
				len1 = current_end - current_start + 1
				len2 = next_end - next_start + 1
				mergedLength = next_end - current_start + 1
				perc1 = gap_len/float(len1)
				perc2 = gap_len/float(len2)
				gap_ratio = min(perc1, perc2)
			else:
				gap_ratio = None
			return (gap_len, gap_ratio, mergedLength,)
		
		@classmethod
		def isTwoSegmentsMergable(cls, segment, next_segment, max_reciprocal_gap_ratio=0.1, max_gap_len=10000,
								maxDeletionLength=50000):
			"""
			2010-6-11
			"""
			current_chr = segment[0]
			next_chr = next_segment[0]
			
			mergable = False
			gap_ratio = None
			if next_chr==current_chr:
				gap_len, gap_ratio, mergedLength = cls.calGapRatioOfTwoSegments(segment, next_segment)[:3]
				if mergedLength<=maxDeletionLength:
					if gap_len<=0:
						#+1 because if they are next to each other, merge them as well
						mergable = True
						gap_ratio = 0
					elif gap_ratio<=max_reciprocal_gap_ratio and gap_len<=max_gap_len:
						mergable = True
			return (gap_ratio, mergable)
		
		@classmethod
		def findChromosomeNeighbor(cls, i, segment_ls, direction=1, max_gap_len=10000):
			"""
			2010-6-11
				finding neighbors of segment_ls[i] on the chromosome in certain direction within max_gap_len.
				max direction
					-1: go left
					+1: go right
			"""
			segment = segment_ls[i]
			current_chr = segment[0]
			current_start = segment[1]
			current_end = segment[2]
			
			count = 0
			next_segment = None
			no_of_segments = len(segment_ls)
			neighbor_index_ls = []
			while 1:
				count += 1
				new_index = i + count*direction
				if new_index<0 or new_index>=no_of_segments:
					next_segment = None
					break
				next_segment = segment_ls[new_index]
				if next_segment:
					next_chr = next_segment[0]
					next_start = next_segment[1]
					next_end = next_segment[2]
					if next_chr!=current_chr:	#beyond this chromosome, no neighbor. break
						next_segment = None
						break
					
					if direction>0:	#to the right 
						gap_len = (next_start - current_end) - 1
					else:	#to the left
						gap_len = current_start - next_end - 1
					if gap_len>max_gap_len:	# too far. stop here
						next_segment = None
						break
					neighbor_index_ls.append(new_index)
			return neighbor_index_ls
			
			
		@classmethod
		def mergeOverlappingORCloseSegmentsByGraph(cls, segment_ls, max_reciprocal_gap_ratio=0.1, max_gap_len=100, \
										mergeFunc=None, maxDeletionLength=50000):
			"""
			2010-6-11
				a graph-based counterpart to mergeOverlappingORCloseMummerCoords()
				
				each segment =
					[chr, ref_start, ref_stop, input_data_index, contig_id_set, query_start, query_stop, identity_perc]
			"""
			sys.stderr.write("Merging overlapping or close-by %s segments by graph ... \n"%(len(segment_ls)))
			if mergeFunc is None:
				mergeFunc=cls.mergeTwoMummerCoords
			segment_ls.sort()
			no_of_segments = len(segment_ls)
			import networkx as nx
			G=nx.Graph()
			
			from heapq import heappush, heappop
			edgeHeap = []
			sys.stderr.write("Constructing graph ... \n")
			counter  = 0
			real_counter = 0
			for i in xrange(no_of_segments):
				segment = segment_ls[i]
				current_chr = segment[0]
				current_start = segment[1]
				current_end = segment[2]
				G.add_node(i)
				for j in xrange(i+1, no_of_segments):
					next_segment = segment_ls[j]
					
					next_chr = next_segment[0]
					next_start = next_segment[1]
					next_end = next_segment[2]
					gap_ratio, mergable = cls.isTwoSegmentsMergable(segment, next_segment, \
										max_reciprocal_gap_ratio=max_reciprocal_gap_ratio, max_gap_len=max_gap_len,\
										maxDeletionLength=maxDeletionLength)[:2]
					if mergable:
						G.add_edge(i, j, gap_ratio)
						#heappush(edgeHeap, (gap_ratio, i,j))
						real_counter += 1
				counter += 1
				if counter%1000==0:
					sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
			sys.stderr.write(" . %s nodes, %s edges. Done.\n"%(G.number_of_nodes(), G.number_of_edges()))
			
			sys.stderr.write("Greedy graph algorithm to merge segments ...\n")
			counter = 0
			real_counter = 0
			import random
			while G.number_of_edges()>0:
				no_of_edges = G.number_of_edges()
				edge_index = random.randint(1, no_of_edges)
				#edge_index = 1
				i, j = G.edges()[edge_index-1]
				try:
					neighbor_set = set(G.neighbors(i))
					neighbor_set |= set(G.neighbors(j))
					
					neighbor_set |= set(cls.findChromosomeNeighbor(i, segment_ls, direction=-1, max_gap_len=max_gap_len))
					
					neighbor_set |= set(cls.findChromosomeNeighbor(i, segment_ls, direction=+1, max_gap_len=max_gap_len))
					
					neighbor_set |= set(cls.findChromosomeNeighbor(j, segment_ls, direction=-1, max_gap_len=max_gap_len))
					
					neighbor_set |= set(cls.findChromosomeNeighbor(j, segment_ls, direction=+1, max_gap_len=max_gap_len))
					
					segment = segment_ls[i]
					next_segment = segment_ls[j]
					
					segment = mergeFunc(segment, next_segment)
					segment_ls[i] = segment
					segment_ls[j] = None
					
					# remove the two old nodes and affiliated edges
					G.remove_node(i)
					G.remove_node(j)
					
					# add one node back
					G.add_node(i)
					
					#i or j could be added when they are next to each other on chromosome
					neighbor_set.remove(i)
					neighbor_set.remove(j)
					
					for node in neighbor_set:
						if i<node:
							first_segment_index = i
							snd_segment_index = node
						elif i>node:
							first_segment_index = node
							snd_segment_index = i
						else:
							continue
						segment = segment_ls[first_segment_index]
						next_segment = segment_ls[snd_segment_index]
						gap_ratio, mergable = cls.isTwoSegmentsMergable(segment, next_segment, \
										max_reciprocal_gap_ratio=max_reciprocal_gap_ratio, max_gap_len=max_gap_len,\
										maxDeletionLength=maxDeletionLength)[:2]
						if mergable:
							G.add_edge(first_segment_index, snd_segment_index, gap_ratio)
							#heappush(edgeHeap, (gap_ratio, first_segment_index, snd_segment_index))
				except:
					sys.stderr.write('Except type at counter=%s, real_counter=%s: %s\n'%(counter, real_counter, repr(sys.exc_info())))
					import traceback
					traceback.print_exc()
					for eH in edgeHeap:
						print eH
					print "%s nodes:"%(G.number_of_nodes())
					for node in G.nodes():
						print node
					print "%s edges:"%(G.number_of_edges())
					for edge in G.edges():
						print edge
					raise
				counter  += 1
				if counter%1000==0:
					sys.stderr.write("%s %s nodes %s edges"%('\x08'*80, G.number_of_nodes(), G.number_of_edges()))
			sys.stderr.write("%s %s nodes %s edges. Done.\n"%('\x08'*80, G.number_of_nodes(), G.number_of_edges()))
			
			return_ls = []
			for segment in segment_ls:
				if segment is not None:
					return_ls.append(segment)
			sys.stderr.write("shrinked to %s segments. Done.\n"%(len(return_ls)))
			return return_ls
		
		@classmethod
		def mergeOverlappingORCloseMummerCoords(cls, refCoverTuple_ls, max_reciprocal_gap_ratio=0.1, max_gap_len=100):
			"""
			2010-6-11
				need to repeatedly to do routine until no more segments could be merged.
			2010-6-10
				
			"""
			no_of_old_tuples = len(refCoverTuple_ls)
			no_of_new_tuples = 0
			round = 0
			while no_of_new_tuples!=no_of_old_tuples:
				no_of_old_tuples = len(refCoverTuple_ls)
				round += 1
				sys.stderr.write("Round %s Merging overlapping or close-by %s mummer coords ..."%(round, len(refCoverTuple_ls)))
				refCoverTuple_ls.sort()
				return_ls = []
				refCoverTuple = refCoverTuple_ls.pop(0)	# start from beginning
				while refCoverTuple_ls:
					refCoverTuple_chr = refCoverTuple[0]
					refCoverTuple_start = refCoverTuple[1]
					refCoverTuple_end = refCoverTuple[2]
					
					next_refCoverTuple = refCoverTuple_ls.pop(0)	# start from beginning
					next_refCoverTuple_chr = next_refCoverTuple[0]
					next_refCoverTuple_start = next_refCoverTuple[1]
					next_refCoverTuple_end = next_refCoverTuple[2]
					mergable = False
					if next_refCoverTuple_chr==refCoverTuple_chr:
						if next_refCoverTuple_start<=refCoverTuple_end+1:
							#+1 because if they are next to each other, merge them as well
							mergable = True
						else:
							try:
								gap_len = next_refCoverTuple_start - refCoverTuple_end - 1
								len1 = refCoverTuple_end - refCoverTuple_start + 1
								len2 = next_refCoverTuple_end - next_refCoverTuple_start + 1
								perc1 = gap_len/float(len1)
								perc2 = gap_len/float(len2)
								if gap_len<=max_gap_len and perc1<=max_reciprocal_gap_ratio and perc2<=max_reciprocal_gap_ratio:
									mergable = True
							except:
								import traceback
								sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
								traceback.print_exc()
								
					if mergable:
						refCoverTuple = cls.mergeTwoMummerCoords(refCoverTuple, next_refCoverTuple)
					else:
						return_ls.append(refCoverTuple)
						refCoverTuple = next_refCoverTuple
				# don't forget the last refCoverTuple
				return_ls.append(refCoverTuple)
				sys.stderr.write("shrinked to %s segments. Done.\n"%(len(return_ls)))
				no_of_new_tuples = len(return_ls)
				refCoverTuple_ls = return_ls
			return refCoverTuple_ls
		
		@classmethod
		def outputMummerCoords(cls, writer, ref_pos_ls, copy_number='NA'):
			"""
			2010-6-10
				In ref_pos_ls, each ref_pos
					= [chr, ref_start, ref_stop, input_data_index, contig_id_set, query_start, query_stop, identity_perc]
				output is
			
			header_row = ['start_probe_id', 'start_chr_pos', 'stop_probe_id', 'stop_chr_pos', 'no_of_probes', \
						'length', 'copy_number', 'contig', 'size_difference']	
			"""
			sys.stderr.write("Outputting %s coords ..."%len(ref_pos_ls))
			for ref_pos in ref_pos_ls:
				chr, ref_start, ref_stop, input_data_index, contig_id_set, query_start, query_stop, identity_perc = ref_pos[:8]
				start_chr_pos = '%s_%s'%(chr, ref_start)
				stop_chr_pos = '%s_%s'%(chr, ref_stop)
				contig_id_ls = list(contig_id_set)
				contig_id_ls.sort()
				contig_id = ','.join(contig_id_ls)
				length = ref_stop - ref_start + 1
				size_difference = -int(length*(1-identity_perc))
				data_row = ['NA', start_chr_pos, 'NA', stop_chr_pos, 'NA', length, copy_number, contig_id, size_difference]
				writer.writerow(data_row)
			sys.stderr.write("Done.\n")
		
		@classmethod
		def discoverLerDeletionAndSpanFromMummerCoordsOutput(cls, nucmer_coords_fname, output_fname_prefix, \
											min_identity_perc=0.9, maxDeletionLength=50000,\
											max_reciprocal_gap_ratio=0.1, max_gap_len=10000, debug=False):
			"""
			2010-6-7
				The idea is similar to discoverLerDeletionDuplication().
				Instead of take a probe as a unit, it takes every nucmer output segment as a unit.
				not implemented yet.
				
				nucmer_coords_fname is already sorted in reference's order 
			"""
			import csv
			from pymodule import getColName2IndexFromHeader, PassingData, figureOutDelimiter
			from pymodule.utils import getColName2IndexFromHeader
			# skip 4 lines.
			reader = csv.reader(open(nucmer_coords_fname), delimiter=figureOutDelimiter(nucmer_coords_fname))
			for i in range(3):
				reader.next()
			
			col_name2index = getColName2IndexFromHeader(reader.next())
			
			sys.stderr.write("Reading in data ...")
			ref_start_col_index = col_name2index.get('[S1]')
			ref_stop_col_index = col_name2index.get('[E1]')
			query_start_col_index = col_name2index.get('[S2]')
			query_stop_col_index = col_name2index.get('[E2]')
			
			identity_col_index = col_name2index.get('[% IDY]')
			ref_query_tag_col_index = col_name2index.get('[TAGS]')
			
			contig_id2chr2input_index_ls = {}
			count = 0
			input_data_ls = []
			for row in reader:
				ref_start = int(row[ref_start_col_index])
				ref_stop = int(row[ref_stop_col_index])
				query_start = int(row[query_start_col_index])
				query_stop = int(row[query_stop_col_index])
				
				identity_perc = float(row[identity_col_index])/100.
				if identity_perc<min_identity_perc:
					continue
				ref_query_tag_ls = row[ref_query_tag_col_index:]
				chr = ref_query_tag_ls[0]
				chr = int(chr[3])
				contig_id = ref_query_tag_ls[1]
				"""
				if contig_id not in contig_id2chr2input_index_ls:
					contig_id2chr2input_index_ls[contig_id] = {}
				if chr not in contig_id2chr2input_index_ls[contig_id]:
					contig_id2chr2input_index_ls[contig_id][chr] = []
				"""
				contig_id_set = set([contig_id])
				data = [chr, ref_start, ref_stop, count, contig_id_set, query_start, query_stop, identity_perc]
				#contig_id2chr2input_index_ls[contig_id][chr].append(count)
				input_data_ls.append(data)
				count += 1
			sys.stderr.write("%s segments. Done.\n"%(len(input_data_ls)))
			
			sys.stderr.write("Discovering deletions ...\n")
			counter = 0
			real_counter = 0
			MinSegLen = 5
			refCoverTuple_ls = []
			deletion_ls = []
			ref_pos = input_data_ls[0]
			no_of_segments = len(input_data_ls)
			for i in range(1, no_of_segments):
				current_ref_start = ref_pos[1]
				current_ref_stop = ref_pos[2]
				next_ref_pos = input_data_ls[i]
				next_ref_start = next_ref_pos[1]
				next_ref_stop = next_ref_pos[2]
				
				if ref_pos[0]==next_ref_pos[0] and next_ref_start<=current_ref_stop+1:
					#+1 because if they are next to each other, merge them as well
					#merge the two
					ref_pos = cls.mergeTwoMummerCoords(ref_pos, next_ref_pos)
				elif ref_pos[0]!=next_ref_pos[0]:	# different chromosome, skip to next_ref_pos
					refCoverTuple_ls.append(ref_pos)
					ref_pos = next_ref_pos
				else:	#there's gap between two
					refCoverTuple_ls.append(ref_pos)
					deletion  = cls.getGapBetweenTwoMummerCoords(ref_pos, next_ref_pos)
					deletion_length = deletion[2]-deletion[1] + 1
					if deletion_length<=maxDeletionLength and deletion_length>0:
						real_counter += 1
						deletion_ls.append(deletion)
					elif deletion_length<=0:
						sys.stderr.write("Error: Deletion length %s <=0.\n"%(deletion_length))
						print "ref_pos:", ref_pos
						print "next_ref_pos:", next_ref_pos
						sys.exit(1)
					ref_pos = next_ref_pos
				counter += 1
				if debug and counter==2000:	# in debug mode. break now
					break
				if counter%1000==0:
					sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
			#don't forget the last ref_pos
			refCoverTuple_ls.append(ref_pos)
			sys.stderr.write("%s ref coverage tuples, %s deletions. Done.\n"%(len(refCoverTuple_ls), len(deletion_ls)))
			
			"""
			sys.stderr.write("Discovering deletions ...\n")
			counter = 0
			real_counter = 0
			MinSegLen = 5
			refCoverTuple_ls = []
			deletion_ls = []
			for contig_id, chr2input_index_ls in contig_id2chr2input_index_ls.iteritems():
				for chr, input_index_ls in chr2input_index_ls.iteritems():
					no_of_ref_pos = len(input_index_ls)
					current_input_index = input_index_ls[0]
					ref_pos = input_data_ls[current_input_index]
					
					for i in range(1, no_of_ref_pos):
						current_input_index = ref_pos[3]	# ref_pos might have changed
						next_input_index = input_index_ls[i]
						for j in range(current_input_index+1, next_input_index+1):	# fill in between
							current_ref_start = ref_pos[1]
							current_ref_stop = ref_pos[2]
							next_ref_pos = input_data_ls[j]
							next_ref_start = next_ref_pos[1]
							next_ref_stop = next_ref_pos[2]
							
							if ref_pos[0]==next_ref_pos[0] and next_ref_start<=current_ref_stop+1:
								#+1 because if they are next to each other, merge them as well
								#merge the two
								ref_pos = cls.mergeTwoMummerCoords(ref_pos, next_ref_pos)
							elif ref_pos[0]!=next_ref_pos[0]:
								print "ref_pos:", ref_pos
								print "next_ref_pos:", next_ref_pos
								sys.stderr.write("Error different chr within same contig.\n")
								sys.exit(1)
							else:	#there's gap between two
								refCoverTuple_ls.append(ref_pos)
								deletion  = cls.getGapBetweenTwoMummerCoords(ref_pos, next_ref_pos)
								deletion_length = deletion[2]-deletion[1] + 1
								if deletion_length<=maxDeletionLength and deletion_length>0:
									deletion_ls.append(deletion)
								elif deletion_length<=0:
									sys.stderr.write("Error: Deletion length %s <=0.\n"%(deletion_length))
									print "ref_pos:", ref_pos
									print "next_ref_pos:", next_ref_pos
									sys.exit(1)
								ref_pos = next_ref_pos
					#don't forget the last ref_pos
					refCoverTuple_ls.append(ref_pos)
					real_counter += 1
				counter += 1
				if debug and counter==50:	# in debug mode. break now
					break
				if counter%1000==0:
					sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
			sys.stderr.write("%s ref coverage tuples, %s deletions. Done.\n"%(len(refCoverTuple_ls), len(deletion_ls)))
			"""
			
			sys.stderr.write("Outputting reference coverage ...")
			delDupWriter = csv.writer(open('%s_delDup.tsv'%output_fname_prefix, 'w'), delimiter='\t')
			refCoveredWriter =  csv.writer(open('%s_ref_covered.tsv'%output_fname_prefix, 'w'), delimiter='\t')
			
			header_row = ['start_probe_id', 'start_chr_pos', 'stop_probe_id', 'stop_chr_pos', 'no_of_probes', \
						'length', 'copy_number', 'contig', 'size_difference']
			for writer in [delDupWriter, refCoveredWriter]:
				writer.writerow(header_row)
			
			"""
			deletion_ls = cls.mergeOverlappingORCloseSegmentsByGraph(deletion_ls, \
							max_reciprocal_gap_ratio=max_reciprocal_gap_ratio, max_gap_len=max_gap_len,\
							maxDeletionLength=maxDeletionLength)
			"""
			cls.outputMummerCoords(delDupWriter, deletion_ls, copy_number=0)
			#refCoverTuple_ls += deletion_ls
			"""
			refCoverTuple_ls = cls.mergeOverlappingORCloseSegmentsByGraph(refCoverTuple_ls, \
							max_reciprocal_gap_ratio=max_reciprocal_gap_ratio, max_gap_len=max_gap_len,
							maxDeletionLength=maxDeletionLength)
			"""
			cls.outputMummerCoords(refCoveredWriter, refCoverTuple_ls, copy_number=1)
			del delDupWriter, refCoveredWriter
			sys.stderr.write("Done.\n")
		"""
		common_prefix = '~/script/variation/data/CNV/mummer/TAIR8_LerContig_maxmatch_maxgap40_mincluster60'
		#common_prefix = '~/script/variation/data/CNV/mummer/TAIR8_LerContig_mum_maxgap40k_mincluster60'
		common_prefix = '~/script/variation/data/CNV/mummer/TAIR8_LerContig_mum_maxgap40_mincluster60'
		nucmer_coords_fname = os.path.expanduser('%s.coords'%(common_prefix))
		max_reciprocal_gap_ratio = 0.1
		max_gap_len = 10000
		min_identity_perc = 0.9
		output_fname_prefix = os.path.expanduser('%s_min_identity_perc%s_max_reciprocal_gap_ratio%s_max_gap_len%s_randomEdge'%\
									(common_prefix, min_identity_perc, max_reciprocal_gap_ratio, max_gap_len))
		CNV.LerContig.discoverLerDeletionAndSpanFromMummerCoordsOutput(nucmer_coords_fname, output_fname_prefix, \
							min_identity_perc=min_identity_perc, maxDeletionLength=50000, \
							max_reciprocal_gap_ratio=max_reciprocal_gap_ratio, max_gap_len=max_gap_len, debug=debug)
		sys.exit(0)
		
		
		# 2010-7-15 TAIR9
		common_prefix = '~/script/variation/data/CNV/mummer/TAIR9_LerContig_maxmatch_maxgap40_mincluster60'
		common_prefix = '~/script/variation/data/CNV/mummer/TAIR9_Lyrata_maxmatch_maxgap40_mincluster60'
		nucmer_coords_fname = os.path.expanduser('%s.coords'%(common_prefix))
		max_reciprocal_gap_ratio = 0.1
		max_gap_len = 10000
		min_identity_perc = 0.9
		output_fname_prefix = os.path.expanduser('%s_min_identity_perc%s'%\
									(common_prefix, min_identity_perc))
		CNV.LerContig.discoverLerDeletionAndSpanFromMummerCoordsOutput(nucmer_coords_fname, output_fname_prefix, \
							min_identity_perc=min_identity_perc, maxDeletionLength=50000, \
							max_reciprocal_gap_ratio=max_reciprocal_gap_ratio, max_gap_len=max_gap_len, debug=self.debug)
		sys.exit(0)
		
		"""
		
		@classmethod
		def investigateNucmerCoordsOutput(cls, nucmer_coords_fname, output_fname):
			"""
			2010-6-9
				try to understand what parameter -maxgap (40 vs 40000) means.
				1. does it mean the ref and query match could have the length difference close to maxgap?
				
			"""
			import csv
			reader = csv.reader(open(nucmer_coords_fname), delimiter='\t')
			
			from pymodule.utils import getColName2IndexFromHeader
			# skip 4 lines
			for i in range(3):
				reader.next()
			
			col_name2index = getColName2IndexFromHeader(reader.next())
			len1_col_index = col_name2index.get('[LEN 1]')
			len2_col_index = col_name2index.get('[LEN 2]')
			identity_col_index = col_name2index.get('[% IDY]')
			ref_query_tag_col_index = col_name2index.get('[TAGS]')
			writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
			for row in reader:
				len1 = int(row[len1_col_index])
				len2 = int(row[len2_col_index])
				identity_perc = float(row[identity_col_index])
				ref_query_tag_ls = row[ref_query_tag_col_index:]
				len_delta = len1-len2
				new_row = row[:4] + [len_delta, len1, len2, identity_perc] + ref_query_tag_ls
				writer.writerow(new_row)
		"""
	common_prefix = '~/script/variation/data/CNV/mummer/TAIR8_LerContig_maxmatch_maxgap40_mincluster60'
	nucmer_coords_fname = os.path.expanduser('%s.coords'%(common_prefix))
	output_fname = os.path.expanduser('%s.coords_delta_len'%(common_prefix))
	CNV.LerContig.investigateNucmerCoordsOutput(nucmer_coords_fname, output_fname)
	sys.exit(0)
	
	common_prefix = '~/script/variation/data/CNV/mummer/TAIR8_LerContig_maxmatch_maxgap40k_mincluster60'
	nucmer_coords_fname = os.path.expanduser('%s.coords'%(common_prefix))
	output_fname = os.path.expanduser('%s.coords_delta_len'%(common_prefix))
	CNV.LerContig.investigateNucmerCoordsOutput(nucmer_coords_fname, output_fname)
	sys.exit(0)
		"""
		
	@classmethod
	def drawIntensityVsCopyNumber(cls, copy_number2intensity_ls, title='', output_fname_prefix=None, max_copy_number=16, \
										plot_type=1, min_no_of_data_points=10, max_no_of_data_points=100000):
		"""
		2010-2-9
			called by CNV.drawProbeIntensityVsTrueCopyNumber()
			argument plot_type
				1: box plot on a XY-plot in which X is copy number, Y is intensity
				2: intensity histograms stratified by copy_number
				3: kernel density function for intensity histograms. all in one plot.
			argument min_no_of_data_points
				aggregate copy number bins whose counts are below min_no_of_data_points
			argument max_no_of_data_points
				if more than this number, sample this number of points
			
		"""
		plot_type2name = {1: 'boxplot', 2: 'histogram'}
		sys.stderr.write("Drawing intensity vs probe copy number in %s for %s ...\n"%(plot_type2name.get(plot_type, ''), title))
		import pylab, random, numpy
		candidate_facecolor_ls = ['r','g', 'b', 'c', 'k']
		
		pylab.clf()
		
		if 0 in copy_number2intensity_ls and 1 in copy_number2intensity_ls:	# do a ks test between copy number 0 and 1
			copy_0_intensity_ls = copy_number2intensity_ls.get(0)
			copy_1_intensity_ls = copy_number2intensity_ls.get(1)
			import rpy
			ksTestResult = rpy.r.ks_test(copy_0_intensity_ls, copy_1_intensity_ls)
			ksPvalue= ksTestResult['p.value']
			ksStat = ksTestResult['statistic'].get('D')
			
			wilcoxTestResult = rpy.r.wilcox_test(copy_0_intensity_ls, copy_1_intensity_ls)
			wilcoxStat = wilcoxTestResult['statistic'].get('W')
		else:
			ksPvalue = None
			ksStat = None
			wilcoxStat = None
		
		if plot_type==1:
			ax1 = pylab.axes([0.05, 0.05, 0.9, 0.9], frameon=False)	# 2010-1-26 remove the frame
			ax1.grid(True, alpha=0.2)
			if title:
				ax1.title.set_text(title)
		elif plot_type==2:
			pylab.subplots_adjust(left=0.08, right=0.92, bottom = 0.08, top=0.92, wspace=0.2, hspace = 0.25)
		elif plot_type==3:
			pylab.grid(True, alpha=0.2)
		
		copy_number_ls = copy_number2intensity_ls.keys()
		copy_number_ls.sort()
			
		intensity_ls_ls = []
		no_of_probes_ls = []
		truncated_copy_number_ls = []
		for copy_number in copy_number_ls:
			if copy_number<=max_copy_number:
				truncated_copy_number_ls.append(copy_number)
				intensity_ls = copy_number2intensity_ls[copy_number]
				intensity_ls_ls.append(intensity_ls)
				no_of_probes_ls.append(len(intensity_ls))
		sys.stderr.write("\t %s different copy numbers shrunk to %s.\n"%(len(copy_number_ls), len(truncated_copy_number_ls)))
		
		aggregateReturnData = cls.aggregateEnoughDataPoints(no_of_data_points_ls = no_of_probes_ls, \
														lists_to_be_summed = [], \
														lists_to_be_sampled = [truncated_copy_number_ls], \
														lists_to_be_concatenated = [intensity_ls_ls], \
														min_no_of_data_points = min_no_of_data_points)
		
		new_no_of_probes_ls = aggregateReturnData.no_of_data_points_ls
		sys.stderr.write("\t further shrunk to %s.\n"%(len(new_no_of_probes_ls)))
		
		new_copy_number_ls = aggregateReturnData.lists_to_be_sampled[0]
		new_intensity_ls_ls = aggregateReturnData.lists_to_be_concatenated[0]
		
		hist_obj_ls = []
		legend_ls = []
		no_of_different_copy_numbers = len(new_copy_number_ls)
		no_of_vertical_axes = 2
		no_of_horizontal_axes = int(math.ceil(no_of_different_copy_numbers/float(no_of_vertical_axes)))
		first_axe = None
		for i in range(no_of_different_copy_numbers):
			copy_number = new_copy_number_ls[i]
			intensity_ls = new_intensity_ls_ls[i]
			intensity_mean = numpy.mean(intensity_ls)
			intensity_var = numpy.var(intensity_ls)
			no_of_probes = new_no_of_probes_ls[i]
			intensity_sigma = numpy.sqrt(intensity_var/(no_of_probes-1))
			if no_of_probes>max_no_of_data_points and plot_type!=3:	# plot_type 3 (density estimation) doesn't need sampling.
				intensity_ls = random.sample(intensity_ls, max_no_of_data_points)
			if plot_type==1:
				if len(intensity_ls)>=min_no_of_data_points:	# more than 10 values, draw a boxplot.
					ax1.boxplot(intensity_ls, positions=[copy_number])
				else:
					x_value_ls = [copy_number]*len(intensity_ls)
					ax1.plot(x_value_ls, intensity_ls, '.')
			elif plot_type==2:
				if len(intensity_ls)>=min_no_of_data_points:	# more than 10 values, draw a histogram.
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
					
					if first_axe is None:
						first_axe = pylab.subplot(no_of_horizontal_axes, no_of_vertical_axes, i+1, frameon=False)
					else:
						pylab.subplot(no_of_horizontal_axes, no_of_vertical_axes, i+1, frameon=False, sharex=first_axe)
					#facecolor_index = len(hist_obj_ls)%len(candidate_facecolor_ls)
					pylab.grid(True, alpha=0.3)
					if copy_number==0 and ksStat is not None:
						tmp_title = 'copy %s ks %.4f wlcx %.1f'%(copy_number, ksStat, wilcoxStat)
					else:
						tmp_title = '%s @ copy %s'%(title, copy_number)
					pylab.title(tmp_title)
					hist_obj = pylab.hist(intensity_ls, alpha=0.6)
					
					if i==no_of_different_copy_numbers-1:	#last axe
						pylab.ylabel('count')
						pylab.xlabel('intensity')
						#hist_obj_ls.append(hist_obj[2][0])
			elif plot_type==3:
				if len(intensity_ls)>=min_no_of_data_points:	# more than 10 values,
					import statistics	# 2010-5-30 package from Michiel De Hoon
					y, x = statistics.pdf(intensity_ls)
					pylab.plot(x, y, alpha=0.7)
					legend_ls.append('copy %s (%.3f, %.4f)'%(copy_number, intensity_mean, intensity_sigma))
					if copy_number==0:
						if ksStat is not None:
							tmp_title = '%s ks %.4f wlcx %.1f'%(title, ksStat, wilcoxStat)
							pylab.title(tmp_title)
						else:
							pylab.title(title)
					if i==no_of_different_copy_numbers-1:	#last axe
						pylab.ylabel('density')
						pylab.xlabel('intensity')
				
		if plot_type==1:
			ax1.set_xlim([-1, max_copy_number+1])
		elif plot_type==2:
			#pylab.legend(hist_obj_ls, legend_ls)
			pass
		elif plot_type==3:
			if legend_ls:
				pylab.legend(legend_ls)
			
		
		if output_fname_prefix:
			output_fname = '%s.png'%(output_fname_prefix)
			pylab.savefig(output_fname, dpi=300)
		sys.stderr.write("Done.\n")
	
	@classmethod
	def construct_array_id2copy_number2intensity_ls_fromProbeXArray(cls, input_fname_ls, ecotype_id2probe_id2copy_number, \
																probe_id_check_set, report=True):
		"""
		2010-5-27
			split out of drawProbeIntensityVsTrueCopyNumber()
		"""
		sys.stderr.write("Getting ProbeXArray intensity data from %s ...\n"%repr(input_fname_ls))
		import fileinput, csv
		from pymodule import getColName2IndexFromHeader, PassingData, figureOutDelimiter
		input_handler = fileinput.input(input_fname_ls)
		reader = csv.reader(input_handler, delimiter=figureOutDelimiter(input_fname_ls[0]))
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header)
		
		no_of_cols = len(header)-3
		
		import Stock_250kDB
		array_id_ls = header[1:-2]
		col_index_to_be_checked = set()
		array_id2ecotype_id = {}
		for array_id in array_id_ls:	# array_id is of type 'str'
			array_id = int(array_id)
			array = Stock_250kDB.ArrayInfo.get(array_id)
			if array.maternal_ecotype_id==array.paternal_ecotype_id and array.maternal_ecotype_id in ecotype_id2probe_id2copy_number:
				col_index_to_be_checked.add(col_name2index[str(array_id)])
				array_id2ecotype_id[array_id] = array.maternal_ecotype_id
		
		array_id2copy_number2intensity_ls = {}
		all_array_copy_number2intensity_ls = {}	# 2010-2-9 with all arrays together
		probe_id_index = col_name2index['probes_id']
		chr_index  = col_name2index['chromosome']
		pos_index =  col_name2index['position']
		counter = 0
		real_counter = 0
		for row in reader:
			if row[0].find("probes_id")!=-1:	# encounter header in one of the files of input_fname_ls
				continue
			probe_id = int(row[probe_id_index])
			if probe_id not in probe_id_check_set:	# save time
				continue
			#chr = int(row[chr_index])
			#pos = int(row[pos_index])
			for j in col_index_to_be_checked:
				array_id = int(header[j])
				ecotype_id = array_id2ecotype_id.get(array_id)
				if ecotype_id is None:
					continue
				probe_id2copy_number = ecotype_id2probe_id2copy_number.get(ecotype_id)
				if probe_id2copy_number is None:
					continue
				copy_number = probe_id2copy_number.get(probe_id)
				if copy_number is not None and row[j]!='nan':	# make sure intensity is not NA
					if array_id not in array_id2copy_number2intensity_ls:
						array_id2copy_number2intensity_ls[array_id] = {}
					if copy_number not in array_id2copy_number2intensity_ls[array_id]:
						array_id2copy_number2intensity_ls[array_id][copy_number] = []
					intensity = float(row[j])
					array_id2copy_number2intensity_ls[array_id][copy_number].append(intensity)
					
					# 2010-2-9
					if copy_number not in all_array_copy_number2intensity_ls:
						all_array_copy_number2intensity_ls[copy_number] = []
					all_array_copy_number2intensity_ls[copy_number].append(intensity)
					
					real_counter += 1
			counter += 1
			if counter%10000==0 and report:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
		
		del reader, input_handler
		sys.stderr.write("In total, %s different copy numbers encompass %s probes. Done.\n"%(len(all_array_copy_number2intensity_ls), \
																							real_counter))
		return all_array_copy_number2intensity_ls, array_id2copy_number2intensity_ls, array_id2ecotype_id
	
	@classmethod
	def constructProbeToCheckIndexLs(cls, probe_id_ls, probe_id_check_set):
		"""
		2010-5-28
			called by construct_array_id2copy_number2intensity_ls_fromArrayXProbe()
			construct this list to reduce the number of probes to check for each qualified array 
		"""
		probe_to_check_index_ls = []
		for i in xrange(len(probe_id_ls)):
			probe_id = probe_id_ls[i]
			if probe_id in probe_id_check_set:
				probe_to_check_index_ls.append(i)
		return probe_to_check_index_ls
		
	
	@classmethod
	def construct_array_id2copy_number2intensity_ls_fromArrayXProbe(cls, input_fname_ls, ecotype_id2probe_id2copy_number, \
																probe_id_check_set=None, report=True):
		"""
		2010-5-27
		"""
		sys.stderr.write("Getting ArrayXProbe intensity data from %s ...\n"%repr(input_fname_ls))
		import fileinput, csv
		from pymodule import getColName2IndexFromHeader, PassingData, figureOutDelimiter
		input_handler = fileinput.input(input_fname_ls)
		reader = csv.reader(input_handler, delimiter=figureOutDelimiter(input_fname_ls[0]))
		header = reader.next()
		probe_id_ls = header[2:]
		probe_id_ls = map(int, probe_id_ls)
		
		probe_to_check_index_ls = cls.constructProbeToCheckIndexLs(probe_id_ls, probe_id_check_set)
		
		chr_ls = reader.next()[2:]
		pos_ls = reader.next()[2:]
		
		array_id2copy_number2intensity_ls = {}
		all_array_copy_number2intensity_ls = {}	# 2010-2-9 with all arrays together
		array_id2ecotype_id = {}
		counter = 0
		real_counter = 0
		for row in reader:
			if row[0]=='':	#it's header in a new file. get a new list of probe ids
				probe_id_ls = row[2:]
				probe_id_ls = map(int, probe_id_ls)
				chr_ls = reader.next()[2:]
				pos_ls = reader.next()[2:]
				probe_to_check_index_ls = cls.constructProbeToCheckIndexLs(probe_id_ls, probe_id_check_set)
				continue
			array_id = int(row[0])
			ecotype_id = int(row[1])
			array_id2ecotype_id[array_id] = ecotype_id
			counter += 1
			
			if ecotype_id not in ecotype_id2probe_id2copy_number:
				continue
			
			probe_id2copy_number = ecotype_id2probe_id2copy_number.get(ecotype_id)
			intensity_ls = row[2:]
			for i in probe_to_check_index_ls:
				probe_id = probe_id_ls[i]
				if probe_id not in probe_id2copy_number:	# save time
					continue
				copy_number = probe_id2copy_number.get(probe_id)
				if intensity_ls[i]!='nan':	# make sure intensity is not NA
					if array_id not in array_id2copy_number2intensity_ls:
						array_id2copy_number2intensity_ls[array_id] = {}
					if copy_number not in array_id2copy_number2intensity_ls[array_id]:
						array_id2copy_number2intensity_ls[array_id][copy_number] = []
					intensity = float(intensity_ls[i])
					array_id2copy_number2intensity_ls[array_id][copy_number].append(intensity)
					
					# 2010-2-9
					if copy_number not in all_array_copy_number2intensity_ls:
						all_array_copy_number2intensity_ls[copy_number] = []
					all_array_copy_number2intensity_ls[copy_number].append(intensity)
					
					real_counter += 1
			if real_counter%5000==0 and report:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
		sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
		del reader, input_handler
		sys.stderr.write("In total, %s different copy numbers encompass %s probes. Done.\n"%(len(all_array_copy_number2intensity_ls), \
																							real_counter))
		return all_array_copy_number2intensity_ls, array_id2copy_number2intensity_ls, array_id2ecotype_id
	
	@classmethod
	def drawProbeIntensityVsTrueCopyNumber(cls, db_250k, input_fname_ls, output_fname_prefix, \
										ecotype_id=6909, data_source_id=3, cnv_type_id=None, \
										min_reciprocal_overlap=0.6, report=True, max_copy_number=16, \
										plot_type=1, min_no_of_data_points=10, max_no_of_data_points=10000,\
										min_no_of_probes=5, x_range=None, y_range=None):
		"""
		2010-2-9
			fix a bug: pass cmpfn=leftWithinRightAlsoEqualCmp to cls.getCNVQCDataFromDB() because here is not to check
				segment overlapping, rather it's to see which QC segment contains which probe.
			add argument plot_type
				1: box plot on a XY-plot in which X is copy number, Y is intensity
				2: intensity histograms stratified by copy_number
			add argument min_no_of_data_points
				aggregate copy number bins whose counts are below min_no_of_data_points
			add argument max_no_of_data_points
				if more than this number, sample this number of points
			split the drawing part into a standalone function, CNV.drawIntensityVsCopyNumber()
			
			if ecotype_id=None, the function will do for every ecotype available in ecotype_id2cnv_qc_call_data.
			
			add argument min_no_of_data_points to limit the QC data to be of certain length.
		2010-1-26
			value of the ecotype_id2cnv_qc_call_data dictionary is a RBDict (RBTree dictionary) structure.
			remove the frame from the final plot, add a grid instead.
			add argument max_copy_number to limit the plot only to probes whose copy numbers are <= this threshold.
		2009-12-9
			purpose is to check whether probe copy number based on Col/Ler blast results actually mean something.
			input_fname_ls is a list of filenames which contains CNV intensity, raw or normalized.
			
			argument min_reciprocal_overlap is not used.
			cnv_type_id = None means all cnv types.
		"""
		from pymodule.CNV import leftWithinRightAlsoEqualCmp
		cnvQCData = cls.getCNVQCDataFromDB(data_source_id, ecotype_id, cnv_type_id, \
															min_reciprocal_overlap=min_reciprocal_overlap,\
															cmpfn=leftWithinRightAlsoEqualCmp, min_no_of_probes=min_no_of_probes)
		ecotype_id2cnv_qc_call_data = cnvQCData.ecotype_id2cnv_qc_call_data
		import Stock_250kDB
		from DB_250k2Array import DB_250k2Array
		probes, xy_ls, chr_pos_ls, total_probe_id_ls = DB_250k2Array.get_probes(db_250k.metadata.bind, Stock_250kDB.Probes.table.name, \
																		snps=None, run_type=2, need_xy_ls=False, need_probes=False,\
																		x_range=x_range, y_range=y_range)
		
		sys.stderr.write("Establish the true copy number for each probe in each ecotype ...")	#2009 need to improve the running time by using CNVSegmentBinarySearchTreeKey from pymodule/CNV.py
		ecotype_id2probe_id2copy_number = {}
		probe_id_check_set = set()
		from pymodule.CNV import CNVSegmentBinarySearchTreeKey
		for i in range(len(total_probe_id_ls)):
			probe_id = total_probe_id_ls[i]
			probe_id_check_set.add(probe_id)
			chr, pos= chr_pos_ls[i]
			probeSegmentKey = CNVSegmentBinarySearchTreeKey(chromosome=chr, span_ls=[pos-12, pos+12], min_reciprocal_overlap=min_reciprocal_overlap)
			for ecotype_id, cnv_qc_call_data in ecotype_id2cnv_qc_call_data.iteritems():
				if ecotype_id not in ecotype_id2probe_id2copy_number:
					ecotype_id2probe_id2copy_number[ecotype_id] = {}
				
				cnv_qc_call = cnv_qc_call_data.get(probeSegmentKey)
				if cnv_qc_call:
					copy_number = cnv_qc_call[5]
					ecotype_id2probe_id2copy_number[ecotype_id][probe_id] = copy_number
		for ecotype_id, probe_id2copy_number in ecotype_id2probe_id2copy_number.iteritems():
			no_of_probes = len(probe_id2copy_number)
			sys.stderr.write("Ecotype %s has %s probes with copy_number assigned.\n"%(ecotype_id, no_of_probes))
		
		import fileinput, csv
		from pymodule import figureOutDelimiter
		input_handler = fileinput.input(input_fname_ls)
		reader = csv.reader(input_handler, delimiter=figureOutDelimiter(input_fname_ls[0]))
		header = reader.next()
		del reader
		input_handler.close()
		if header[0] == 'probes_id':
			all_array_copy_number2intensity_ls, array_id2copy_number2intensity_ls, array_id2ecotype_id = \
				cls.construct_array_id2copy_number2intensity_ls_fromProbeXArray(\
																input_fname_ls, ecotype_id2probe_id2copy_number, \
																probe_id_check_set, report=report)
		else:
			all_array_copy_number2intensity_ls, array_id2copy_number2intensity_ls, array_id2ecotype_id = \
				cls.construct_array_id2copy_number2intensity_ls_fromArrayXProbe(\
																input_fname_ls, ecotype_id2probe_id2copy_number, \
																probe_id_check_set, report=report)
		
		if len(all_array_copy_number2intensity_ls)>0:
			title = "Intensity Histogram For all arrays (source-id=%s)"%(data_source_id)
			output_fname_prefix_tmp = '%s_data_source_%s'%(output_fname_prefix, data_source_id)
			cls.drawIntensityVsCopyNumber(all_array_copy_number2intensity_ls, title=title, output_fname_prefix=output_fname_prefix_tmp, \
											max_copy_number=max_copy_number, plot_type=plot_type, \
											min_no_of_data_points=min_no_of_data_points, max_no_of_data_points=max_no_of_data_points)
		
		for array_id, copy_number2intensity_ls in array_id2copy_number2intensity_ls.iteritems():
			if len(copy_number2intensity_ls)>0:
				ecotype_id = array_id2ecotype_id.get(array_id)
				title = "Array %s E-id %s"%(array_id, ecotype_id)
				output_fname_prefix_tmp = '%s_array_%s_ecotype_%s_data_source_%s'%(output_fname_prefix, array_id, ecotype_id, data_source_id)
				cls.drawIntensityVsCopyNumber(copy_number2intensity_ls, title=title, output_fname_prefix=output_fname_prefix_tmp, \
											max_copy_number=max_copy_number, plot_type=plot_type, \
											min_no_of_data_points=min_no_of_data_points, max_no_of_data_points=max_no_of_data_points)
		return array_id2copy_number2intensity_ls
	"""
	# 2010-1-26 intensity vs copy number
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s.tsv'%i))
	
	input_fname_ls = [os.path.expanduser('~/panfs/250k/CNV/call_method_48_CNV_intensity.tsv')]
	output_fname_prefix = '/tmp/call_48_Col_intensity_vs_true_copy_number'
	CNV.drawProbeIntensityVsTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
										ecotype_id=6909, data_source_id=9, cnv_type_id=None)
	
	max_copy_number = 16
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s.tsv'%i))
	output_fname_prefix = os.path.expanduser('~/tmp/call_48_Col_intensity_QNorm_sub_ref_vs_true_copy_number_m%s'%max_copy_number)
	CNV.drawProbeIntensityVsTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
										ecotype_id=6909, data_source_id=9, cnv_type_id=None, max_copy_number=max_copy_number)
	
	input_fname_ls = [os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity.tsv')]
	output_fname_prefix = os.path.expanduser('~/tmp/call_48_Col_intensity_vs_true_copy_number_m%s'%max_copy_number)
	CNV.drawProbeIntensityVsTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
										ecotype_id=6909, data_source_id=9, cnv_type_id=None, max_copy_number=max_copy_number)
	
	input_fname_ls = [os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity.tsv')]
	output_fname_prefix = os.path.expanduser('~/tmp/call_48_Ler_intensity_vs_true_copy_number_m%s'%max_copy_number)
	CNV.drawProbeIntensityVsTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
										ecotype_id=6932, data_source_id=8, cnv_type_id=None, max_copy_number=max_copy_number)
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s.tsv'%i))
	output_fname_prefix = os.path.expanduser('~/tmp/call_48_Ler_intensity_QNorm_sub_ref_vs_true_copy_number_m%s'%max_copy_number)
	CNV.drawProbeIntensityVsTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
										ecotype_id=6932, data_source_id=8, cnv_type_id=None, max_copy_number=max_copy_number)
	
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_48_CNV_QNormalize_chr%s.tsv'%i))
	output_fname_prefix = os.path.expanduser('~/tmp/call_48_Ler_intensity_QNorm_vs_true_copy_number_m%s'%max_copy_number)
	CNV.drawProbeIntensityVsTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
										ecotype_id=6932, data_source_id=8, cnv_type_id=None, max_copy_number=max_copy_number)
	
	# 2010-2-9 draw intensity histograms conditioning on the copy number (plot_type=2)
	max_copy_number = 6
	input_fname_ls = [os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity.tsv')]
	output_fname_prefix = os.path.expanduser('~/tmp/call_48_Col_intensity_hist_by_copy_number_m%s'%max_copy_number)
	CNV.drawProbeIntensityVsTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
										ecotype_id=6909, data_source_id=9, cnv_type_id=None, max_copy_number=max_copy_number, \
										plot_type=2)
	
	input_fname_ls = [os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity.tsv')]
	output_fname_prefix = os.path.expanduser('~/tmp/call_48_Ler_intensity_hist_by_copy_number_m%s'%max_copy_number)
	CNV.drawProbeIntensityVsTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
										ecotype_id=6932, data_source_id=8, cnv_type_id=None, max_copy_number=max_copy_number, \
										plot_type=2)
	
	# 2010-2-10  intensity histograms conditioning on the copy number (plot_type=2) for Col
	for min_no_of_probes in [5, 10, 15]:
		input_fname_ls = []
		for i in range(1,6):
			input_fname_ls.append(os.path.expanduser('~/mnt/panfs/250k/CNV/call_48_CNV_QNormalize_chr%s.tsv'%i))
		output_fname_prefix = os.path.expanduser('~/tmp/call_48_Col_QNorm_intensity_hist_by_copy_number_max_copy%s_min_no_of_probes%s'%\
												(max_copy_number, min_no_of_probes))
		CNV.drawProbeIntensityVsTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
											ecotype_id=6909, data_source_id=9, cnv_type_id=None, max_copy_number=max_copy_number, \
											plot_type=2, min_no_of_probes=min_no_of_probes)
	# 2010-2-10 intensity histograms conditioning on the copy number (plot_type=2) for Ler
	for min_no_of_probes in [5, 10, 15]:
		input_fname_ls = []
		for i in range(1,6):
			input_fname_ls.append(os.path.expanduser('~/mnt/panfs/250k/CNV/call_48_CNV_QNormalize_chr%s.tsv'%i))
		output_fname_prefix = os.path.expanduser('~/tmp/call_48_Ler_QNorm_intensity_hist_by_copy_number_max_copy%s_min_no_of_probes%s'%\
												(max_copy_number, min_no_of_probes))
		CNV.drawProbeIntensityVsTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
											ecotype_id=6932, data_source_id=8, cnv_type_id=None, max_copy_number=max_copy_number, \
											plot_type=2, min_no_of_probes=min_no_of_probes)
	
	max_copy_number = 5
	# 2010-2-10 intensity histograms conditioning on the copy number (plot_type=2) for Clark2007a
	for min_no_of_probes in [5, 10, 15]:
		input_fname_ls = []
		for i in range(1,6):
			input_fname_ls.append(os.path.expanduser('~/mnt/panfs/250k/CNV/call_48_CNV_QNormalize_chr%s.tsv'%i))
		output_fname_prefix = os.path.expanduser('~/tmp/call_48_Clark2007a_Deletion_QNorm_intensity_hist_by_copy_number_max_copy%s_min_no_of_probes%s'%\
												(max_copy_number, min_no_of_probes))
		CNV.drawProbeIntensityVsTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
											ecotype_id=None, data_source_id=1, cnv_type_id=None, max_copy_number=max_copy_number, \
											plot_type=2, min_no_of_probes=min_no_of_probes)
	
	# 2010-2-10 probe intensity histograms for wave-corrected Col-0 intensities
	for min_no_of_probes in [5, 10, 15]:
		input_fname_ls = []
		for i in range(1,6):
			input_fname_ls.append(os.path.expanduser('~/mnt/panfs/250k/CNV/ColLerCviVanFeiSha/call_48_CNV_QNormalize_chr%s_wave_corrected.tsv'%i))
		output_fname_prefix = os.path.expanduser('~/tmp/call_48_CNV_QNorm_wave_corrected_hist_by_copy_number_max_copy%s_min_no_of_probes%s'%\
												(max_copy_number, min_no_of_probes))
		CNV.drawProbeIntensityVsTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
											ecotype_id=6909, data_source_id=9, cnv_type_id=None, max_copy_number=max_copy_number, \
											plot_type=2, min_no_of_probes=min_no_of_probes)
	
	# 2010-4-30 another wave correction for Col-0
	for min_no_of_probes  in [1,3,10,15]:
		input_fname_ls = []
		for i in range(1,6):
			input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/ColLerCviVanFeiSha/call_48_CNV_QNormalize_chr%s_intensity_1.5_4.0_wave_corrected.tsv'%i))
		output_fname_prefix = os.path.expanduser('~/tmp/call_48_CNV_QNorm_intensity_1.5_4.0_wave_corrected_hist_by_copy_number_max_copy%s_min_no_of_probes%s'%\
												(max_copy_number, min_no_of_probes))
		CNV.drawProbeIntensityVsTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
											ecotype_id=6909, data_source_id=9, cnv_type_id=None, max_copy_number=max_copy_number, \
											plot_type=2, min_no_of_probes=min_no_of_probes)
	
	# 2010-2-10 probe intensity histograms for wave-corrected Ler intensities
	for min_no_of_probes in [5, 10, 15]:
		for data_source_id in [7,8]:
			input_fname_ls = []
			for i in range(1,6):
				input_fname_ls.append(os.path.expanduser('~/mnt/panfs/250k/CNV/ColLerCviVanFeiSha/call_48_CNV_QNormalize_chr%s_wave_corrected.tsv'%i))
			output_fname_prefix = os.path.expanduser('~/tmp/call_48_CNV_QNorm_wave_corrected_hist_by_copy_number_max_copy%s_min_no_of_probes%s'%\
												(max_copy_number, min_no_of_probes))
			CNV.drawProbeIntensityVsTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
												ecotype_id=6932, data_source_id=data_source_id, cnv_type_id=None, \
												max_copy_number=max_copy_number, \
												plot_type=2, min_no_of_probes=min_no_of_probes)
	
	# 2010-2-10 probe intensity histograms for wave-corrected Ler intensities
	for min_no_of_probes in [5, 10, 15]:
		for data_source_id in [7,8]:
			input_fname_ls = []
			for i in range(1,6):
				input_fname_ls.append(os.path.expanduser('~/mnt/panfs/250k/CNV/ColLerCviVanFeiSha/call_48_CNV_QNormalize_chr%s_intensity_2.7_3.1_span_2000_probes_wave_corrected.tsv'%i))
			output_fname_prefix = os.path.expanduser('~/tmp/call_48_CNV_QNorm_wave_corrected_2.7_3.1_span_2000_probes_hist_by_copy_number_max_copy%s_min_no_of_probes%s'%\
												(max_copy_number, min_no_of_probes))
			CNV.drawProbeIntensityVsTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
												ecotype_id=6932, data_source_id=data_source_id, cnv_type_id=None, \
												max_copy_number=max_copy_number, \
												plot_type=2, min_no_of_probes=min_no_of_probes)
	
	# 2010-5-10 probe intensity histograms for subArrayQNorm or LeiLi's normalization
	max_copy_number = 5
	common_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/normalization/2010-05-07-quantile-100/medianOverAllRef/medianOverAllRef_probeType_2')
	common_fname_prefix = os.path.expanduser('~/mnt/panfs/250k/CNV/call_48_CNV_2_150_subArrayQNorm_b100_j25_m1666_obp0.85')
	for min_no_of_probes in [5, 10, 15]:
		for data_source_id in [7,8,12]:
			input_fname_ls = []
			input_fname_ls.append('%s.tsv'%common_fname_prefix)
			output_fname_prefix = '%s_hist_by_copy_number_max_copy%s_min_no_of_probes%s'%\
						(common_fname_prefix, max_copy_number, min_no_of_probes)
			CNV.drawProbeIntensityVsTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
												ecotype_id=6932, data_source_id=data_source_id, cnv_type_id=None, \
												max_copy_number=max_copy_number, \
												plot_type=2, min_no_of_probes=min_no_of_probes)
	sys.exit(0)
	
	# 2010-5-10 juxtapose probe intensity histograms from different data sources together (good copy 0 vs copy 1) 
	max_copy_number = 5
	common_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/normalization/2010-05-07-quantile-100/medianOverAllRef/medianOverAllRef_probeType_2')
	#common_fname_prefix = os.path.expanduser('~/mnt/panfs/250k/CNV/call_48_CNV_2_150_subArrayQNorm_b100_j25_m1666_obp0.85')
	#common_fname_prefix = os.path.expanduser('~/mnt/panfs/250k/CNV/call_48_CNV_2_150_subArrayQNorm_b100_j25_m1666_obp0.85')
	
	#common_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/normalization/medianOverRef-w50-lts/subArrayQuantileNorm_3_4_41_150_151_MedianAmong_Ref_1_2_43_139_145_probeType_2')
	common_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/normalization/medianOverRef-w75-qnorm/subArrayQuantileNorm_3_4_41_150_151_MedianAmong_Ref_1_2_139_145_probeType_2')
	common_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/call_48_Col-Ler-WABQN_b200_j100_lts')
	#common_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/call_48_WABQN_b200_j100_lts_2Flanking')
	y_range = None	#[1400, 1612]	#check small range to see if it works
	x_range = None	#[0, 200]
	
	min_no_of_probes = 5
	array_id2copy_number2intensity_ls_by_different_source = {}
	data_source_id_ls = [7,8]
	plot_type = 3	# kernel density estimation
	for data_source_id in data_source_id_ls:
		input_fname_ls = []
		input_fname_ls.append('%s.tsv'%common_fname_prefix)
		output_fname_prefix = '%s_hist_by_copy_number_max_copy%s_min_no_of_probes%s'%\
					(common_fname_prefix, max_copy_number, min_no_of_probes)
		
		#input_fname_ls = []
		#for i in range(1,6):
		#	input_fname_ls.append(os.path.expanduser('~/mnt/panfs/250k/CNV/call_48_CNV_QNormalize_chr%s.tsv'%i))
		#output_fname_prefix = os.path.expanduser('~/mnt/panfs/250k/CNV/call_48_CNV_QNormalize_hist_by_copy_number_max_copy%s_min_no_of_probes%s'%\
		#										(max_copy_number, min_no_of_probes))
		
		array_id2copy_number2intensity_ls = CNV.drawProbeIntensityVsTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
											ecotype_id=6932, data_source_id=data_source_id, cnv_type_id=None, \
											max_copy_number=max_copy_number, \
											plot_type=plot_type, min_no_of_probes=min_no_of_probes, x_range=x_range, y_range=y_range)
		if data_source_id==data_source_id_ls[1]:
			wanted_copy_number = 1
		elif data_source_id==data_source_id_ls[0]:
			wanted_copy_number = 0
		for array_id, copy_number2intensity_ls in array_id2copy_number2intensity_ls.iteritems():
			if array_id not in array_id2copy_number2intensity_ls_by_different_source:
				array_id2copy_number2intensity_ls_by_different_source[array_id] = {}
			array_id2copy_number2intensity_ls_by_different_source[array_id][wanted_copy_number] = copy_number2intensity_ls[wanted_copy_number]
	
	min_no_of_data_points=10
	max_no_of_data_points=10000
	for array_id, copy_number2intensity_ls in array_id2copy_number2intensity_ls_by_different_source.iteritems():
		if len(copy_number2intensity_ls)>0:
				title = "Array %s"%(array_id)
				output_fname_prefix_tmp = '%s_array_%s_source_%s_vs_%s'%(output_fname_prefix, array_id, data_source_id_ls[0], data_source_id_ls[1])
				CNV.drawIntensityVsCopyNumber(copy_number2intensity_ls, title=title, output_fname_prefix=output_fname_prefix_tmp, \
											max_copy_number=max_copy_number, plot_type=plot_type, \
											min_no_of_data_points=min_no_of_data_points, max_no_of_data_points=max_no_of_data_points)
	sys.exit(0)
	
	"""
	
	class exportCNVCallTraverseCNVCallProcessor(object):
		"""
		2010-9-10
			used by exportCNVCallData which plugs it into traverseCNVCall()
		"""
		def __init__(self, ecotype_id2nativename=None, **keywords):
			"""
			2010-7-29
				add argument **keywords
				which contains writer
			"""
			self.array_id2data = {}
			self.ecotype_id2nativename = ecotype_id2nativename
			self.real_counter = 0
			
			
			self.dataType = 1 # 2010-7-29
			for keyword, value in keywords.iteritems():
				setattr(self, keyword, value)
			
		def run(self, row, param_obj=None):
			"""
			This function assumes row comes in 
			"""
			array_id = row.array_id
			ecotype_id = row.ecotype_id
			nativename = self.ecotype_id2nativename.get(ecotype_id)
			data_row = [row.ecotype_id, nativename, array_id, row.chromosome, row.start, row.stop, row.stop-row.start+1,\
					row.no_of_probes_covered, row.amplitude, row.probability]
			self.writer.writerow(data_row)
			self.real_counter += 1
			if self.real_counter%5000==0:
				sys.stderr.write("%s%s"%('\x08'*80, self.real_counter))
			
	
	@classmethod
	def exportCNVCallData(cls, db_250k, cnv_method_id=16, output_fname=None, cnv_type_id=None, ecotype_id=None, \
						**keywords):
		"""
		2010-9-10
			export data from table CNVCall given call method.
		"""
		ecotype_id_ls = None
		
		from pymodule import PassingData
		from common import get_ecotypeid2nativename
		ecotype_id2nativename = get_ecotypeid2nativename(db_250k.metadata.bind)
		import math, csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = ['ecotype_id', 'nativename', 'array_id', 'chr', 'start', 'stop', 'size', 'no_of_probes', 'amplitude', 'probability']
		writer.writerow(header)
		
		processor = cls.exportCNVCallTraverseCNVCallProcessor(ecotype_id2nativename=ecotype_id2nativename, writer = writer, \
														**keywords)
		
		cls.traverseCNVCall(db_250k, cnv_method_id=cnv_method_id, cnv_type_id=cnv_type_id, \
					ecotype_id=ecotype_id, ecotype_id_ls=ecotype_id_ls,\
					min_no_of_probes=None, min_segment_size=None, processClassIns=processor, \
					param_obj=None)
		del writer
	
	"""
		#2010-9-10
		cnv_method_id = 16
		output_fname = '/tmp/cnv_method_%s.tsv'%cnv_method_id
		CNV.exportCNVCallData(db_250k, cnv_method_id=cnv_method_id, output_fname=output_fname)
		sys.exit(0)
	"""
	
	@classmethod
	def calculateArrayOutlierAlongChromosome(cls, input_fname,  output_fname_prefix, halfWindowSize=5, \
											drawHist=True):
		"""
		2010-5-9
			calculate the difference between intensity of one probe and the median of the surrounding 10 probes.
			do it in overlapping windows along the chromosome
				if one probe has multiple difference value, take the median of all for that probe
		"""
		sys.stderr.write("Calulating array outlier  for %s ...\n"%input_fname)
		from CNVNormalize import CNVNormalize
		import numpy
		from pymodule.utils import addExtraLsToFilenamePrefix, addExtraToFilenamePrefix
		data_matrix, probe_id_ls, chr_pos_ls, header = CNVNormalize.get_input(input_fname)
		array_id_ls = header[1:-2]
		no_of_arrays = len(array_id_ls)
		col_index_to_check_ls = range(no_of_arrays)
		
		no_of_rows, no_of_cols = data_matrix.shape
		
		no_of_arrays_to_fit = no_of_arrays
		
		outlier_data_matrix = numpy.zeros([no_of_rows, no_of_arrays_to_fit], numpy.float32)
		for k in range(no_of_arrays):
			j = col_index_to_check_ls[k]
			intensity_ls = data_matrix[:,j]
			array_id = array_id_ls[j]
			sys.stderr.write("array %s ..."%array_id)
			weight_ls = []
			for i in range(no_of_rows):
				left_start_index = max(0, i-halfWindowSize)
				left_stop_index = max(0, i-1)
				right_start_index = min(no_of_rows-1, i+1)
				right_stop_index = min(no_of_rows-1, i+halfWindowSize)
				neighbor_probe_index_ls = range(left_start_index, left_stop_index+1) + range(right_start_index, right_stop_index+1)
				neighbor_probe_intensity_ls = data_matrix[neighbor_probe_index_ls, j]
				neighbor_median = numpy.median(neighbor_probe_intensity_ls)
				delta = data_matrix[i,j] - neighbor_median
				outlier_data_matrix[i,k] = delta
			### plot a histogram here ...
			hist_fig_fname = addExtraToFilenamePrefix(output_fname_prefix, 'array_%s_w%s_hist'%(array_id, halfWindowSize))
			if drawHist:
				from PlotQCProbeIntensityHistogram import PlotQCProbeIntensityHistogram
				PlotQCProbeIntensityHistogram.plotHistogram(outlier_data_matrix[:,k], array_id, hist_fig_fname, \
							no_of_bins=50, \
							max_intensity=None, \
							xlim=None)
			sys.stderr.write("Done.\n")
		wave_output_fname = '%s_wave.tsv'%output_fname_prefix
		new_header = header
		CNVNormalize.output(outlier_data_matrix, probe_id_ls, chr_pos_ls, new_header, wave_output_fname)
		sys.stderr.write("Done.\n")
	"""
	# 2010-5-10
	for i in range(1,6):
		input_fname = os.path.expanduser('~/mnt/panfs/250k/CNV/call_48_CNV_intensity_chr%s_LerArrays.tsv'%i)
		output_fname_prefix = os.path.expanduser('~/mnt/panfs/250k/CNV/call_48_CNV_intensity_chr%s_outlier'%i)
		CNV.calculateArrayOutlierAlongChromosome(input_fname, output_fname_prefix, halfWindowSize=5)
	"""
	
	@classmethod
	def markCNVProbeTAIR9CopyNumber(cls, db_250k, tair9_blast_result_fname, min_no_of_matches=25):
		"""
		2010-6-23
			deal with the possibility that one probe has mutliple copies
		2010-5-23
			tair9_blast_result_fname is output of 
				~/scriptvariation/Dazhe/src/MpiBlast.py -o ~/script/variation/data/CNV/cnv_probe_blast_against_tair9.tsv -m 25
				-d ~/script/variation/data/CNV/TAIR9/tair9.fas
			
			This program marks Tair9Copy of probes in database to be 1 if the blast result shows it has a full copy in Tair9.
		"""
		from pymodule import getColName2IndexFromHeader, PassingData, figureOutDelimiter
		sys.stderr.write("Reading from %s ... \n"%tair9_blast_result_fname)
		counter = 0
		real_counter = 0
		session = db_250k.session
		import csv, Stock_250kDB
		reader = csv.reader(open(tair9_blast_result_fname), delimiter=figureOutDelimiter(tair9_blast_result_fname))
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header)
		probe_id2copy_number = {}
		for row in reader:
			counter += 1
			chr = int(row[col_name2index['Chromosome']])
			pos = int(row[col_name2index['Position']])
			probe_pos = (chr, pos)
			contig_label = row[col_name2index['Alignment_title']]
			probe_id = int(row[col_name2index['Probe_ID']])
			no_of_matches = int(row[col_name2index['Number_matches']])
			alignment_start = int(row[col_name2index['Alignment_start']])
			query_start = int(row[col_name2index['query_start']])
			if no_of_matches>=min_no_of_matches:
				if probe_id not in probe_id2copy_number:
					probe_id2copy_number[probe_id] = 0
				probe_id2copy_number[probe_id] += 1
				real_counter += 1
			if counter%5000==0:
				sys.stderr.write("%s\t%s\t%s"%('\x08'*80, counter, real_counter))
		del reader
		sys.stderr.write("%s\t%s\t%s perfect matches for %s distinct probes. Done.\n"%('\x08'*80, \
													counter, real_counter, len(probe_id2copy_number)))
		
		sys.stderr.write("Updating the database ...\n")
		i = 0
		block_size = 1000
		real_counter = 0
		TableClass = Stock_250kDB.Probes
		rows = TableClass.query.offset(i).limit(block_size)
		while rows.count()!=0:
			for row in rows:
				copy_number = probe_id2copy_number.get(row.id)
				if row.Tair9Copy!=copy_number:	# copy_number could be None, which sets Tair9Copy to NULL.
					row.Tair9Copy = copy_number
					session.add(row)
					real_counter += 1
				i += 1
			if i%5000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, i, real_counter))
				session.flush()
				session.expunge_all()
			rows = TableClass.query.offset(i).limit(block_size)
		session.flush()
		session.commit()
		sys.stderr.write("%s\t%s\t%s Done.\n"%('\x08'*80, i, real_counter))
	
	"""
	tair9_blast_result_fname = os.path.expanduser("~/script/variation/data/CNV/cnv_probe_blast_against_tair9.tsv")
	CNV.markCNVProbeTAIR9CopyNumber(db_250k, tair9_blast_result_fname)
	sys.exit(0)
	
	input_fname = os.path.expanduser('~/mnt/hpc-cmb/script/variation/data/CNV/stock_250k_probes_vs_TAIR9.tsv')
	CNV.markCNVProbeTAIR9CopyNumber(db_250k, input_fname)
	sys.exit(0)
	
	#2010-6-26 the hpc-cmb result is incomplete due to early MPI exit.
	input_fname = os.path.expanduser('~/script/variation/data/CNV/stock_250k_probes_vs_TAIR9_20100625.tsv')
	CNV.markCNVProbeTAIR9CopyNumber(db_250k, input_fname)
	sys.exit(0)
	"""
	
	@classmethod
	def plotXY(cls, x_ls, y_ls, title=None, output_fname=None, xlabel=None, ylabel=None, **keywords):
		"""
		2010-5-18
		"""
		import pylab
		pylab.clf()
		if xlabel:
			pylab.xlabel(xlabel)
		if ylabel:
			pylab.ylabel(ylabel)
		if title:
			pylab.title(title)
		pylab.plot(x_ls, y_ls, '.', **keywords)
		if output_fname:
			pylab.savefig(output_fname, dpi=300)
	
	class CNVPredictionBySVM(object):
		"""
			2010-6-30
		"""
		def __init__(self):
			"""
			2010-6-30
			"""
			pass
		
		@classmethod
		def outputDeletionStatForEachArrayData(cls,  db_250k, cnv_type_id=1, cnv_method_id=10, \
						replaceAmpWithMedianIntensity=False, output_fname=None, debug=False):
			"""
			2010-7-2
				output no_of_deletions, total_genome_size_deleted for each array, along with ecotype id and name
			"""
			i = 0
			block_size = 5000
			real_counter = 0
			import Stock_250kDB
			TableClass = Stock_250kDB.CNVCall
			query = TableClass.query.filter_by(cnv_method_id=cnv_method_id).filter_by(cnv_type_id=cnv_type_id)
			rows = query.offset(i).limit(block_size)
			session = db_250k.session
			array_id2data = {}
			
			from pymodule import PassingData
			from common import get_ecotypeid2nativename
			ecotype_id2nativename = get_ecotypeid2nativename(db_250k.metadata.bind)
			ecotype_id = None
			percUnCoveredByLerContig_ls = []
			feature_data = []
			class_label_ls = []
			c_ls = []
			while rows.count()!=0:
				for row in rows:
					ecotype_id = row.array.maternal_ecotype_id
					array_id = row.array_id
					if array_id not in array_id2data:
						nativename = ecotype_id2nativename.get(ecotype_id)
						array_id2data[array_id] = PassingData(no_of_deletions=0, deletion_size=0, ecotype_id=ecotype_id,\
															nativename=nativename)
					array_id2data[array_id].no_of_deletions += 1
					array_id2data[array_id].deletion_size += (row.stop-row.start+1.0)
					i += 1
				if i%5000==0:
					sys.stderr.write("%s%s\t%s"%('\x08'*80, i, real_counter))
					if debug:
						break
				rows = query.offset(i).limit(block_size)
			sys.stderr.write("%s%s\t%s\n"%('\x08'*80, i, real_counter))
			
			sys.stderr.write("Outputting deletion stats ...")
			import csv
			writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
			header = ['array_id', 'nativename', 'ecotype_id', 'no_of_deletions', 'deletion_size']
			writer.writerow(header)
			for array_id, stat_data in array_id2data.iteritems():
				row = [array_id, stat_data.nativename, stat_data.ecotype_id, stat_data.no_of_deletions, stat_data.deletion_size]
				writer.writerow(row)
			del writer
			sys.stderr.write("Done.\n")
		"""
		#2010-7-2
		output_fname = os.path.expanduser('~/script/variation/data/CNV/cnv_method_10_deletion_stat.tsv')
		CNV.CNVPredictionBySVM.outputDeletionStatForEachArrayData(db_250k, \
						cnv_type_id=1, cnv_method_id=10, \
						replaceAmpWithMedianIntensity=False, output_fname=output_fname)
			
		"""
		@classmethod
		def getCNVFeatureData(cls,  db_250k, array_id=None, \
						minPercUnCoveredByLerContig=0.6, cnv_method_id=6, \
						replaceAmpWithMedianIntensity=False):
			"""
			2010-7-1
				moved to CNVPredictDeletionBySVM.py
			2010-6-30
			"""
			from CNVPredictDeletionBySVM import CNVPredictDeletionBySVM
			"""
			2010-7-1 "from variation.src.CNVPredictDeletionBySVM import CNVPredictDeletionBySVM" doesn't work
				because  db_250k in this misc.py is directly imported from Stock_250kDB, rather than variation.src.Stock_250kDB.
			"""
			return CNVPredictDeletionBySVM.getCNVFeatureData(db_250k, array_id=array_id, \
						minPercUnCoveredByLerContig=minPercUnCoveredByLerContig, cnv_method_id=cnv_method_id, \
						replaceAmpWithMedianIntensity=replaceAmpWithMedianIntensity)
		
		@classmethod
		def getCNVPredictionErrorBySVMCrossValidation(cls, db_250k, output_dir=None, array_id=None, \
						minPercUnCoveredByLerContig=0.6, cnv_method_id=6, \
						replaceAmpWithMedianIntensity=False, gridsize=40):
			"""
			2010-6-30
			"""
			import sys
			cnvFeatureData = cls.getCNVFeatureData(db_250k, array_id=array_id, \
						minPercUnCoveredByLerContig=minPercUnCoveredByLerContig, cnv_method_id=cnv_method_id, \
						replaceAmpWithMedianIntensity=replaceAmpWithMedianIntensity)
			
			from svm import svm_parameter, LINEAR, POLY, RBF
			import numpy
			class_label_ls = cnvFeatureData.class_label_ls
			size = len(class_label_ls)
			if size<100:
				return
			sys.stderr.write("Starting the prediciton ...\n")
			total_correct = 0.
			kernels2kname = {LINEAR:'linear', POLY:'polynomial', RBF:'rbf'}
			kernels2kname = {RBF:'rbf'}	# 2009- 11-18 only polynomial
			mute_device = open('/dev/null', 'w')
			eps = 1e-2
			param = svm_parameter(C = 10, eps=eps)	#,nr_weight = 2,weight_label = [1,0],weight = [10,1])
			C_ls = [100, 10, 1, 0.1, 0.01]	# tried 1000, 10, 1, 0.1 ,0.01
			C_ls = [10]
			gamma_ls = [100, 10, 1, 0.1, 0.0001, 0]
			gamma_ls = [0.]	# test shows that gamma affects rbf quite a bit
			nr_fold = 10	# 10 fold cross validation
			for k, kname in kernels2kname.iteritems():
				for C in C_ls:
					for gamma in gamma_ls:
						param.gamma = gamma
						param.C = C
						param.kernel_type = k;
						param.eps = eps
						
						fitResult = cls.fitAndPredictBySVM(param, cnvFeatureData.feature_data, class_label_ls,\
											crossValidationNrFold=nr_fold)
						cls.drawSVMClassPrediction(output_dir, cnvFeatureData.feature_data, fitResult.prediction_ls, \
							gridsize=gridsize, \
							extraFileNameTag = 'ByCrossValidation', \
							extraTitle='ByCV, PPV %.3f TP %s TPR %.3f FPR %.3f'%(fitResult.precision, \
											fitResult.TP, \
											fitResult.sensitivity, fitResult.FPR),\
							class_label_name='segment label', \
							array_id=array_id, ecotype_id=cnvFeatureData.ecotype_id, cnv_method_id=cnv_method_id, \
							minPercUnCoveredByLerContig=minPercUnCoveredByLerContig)
		"""
		output_dir = os.path.expanduser('~/script/variation/data/CNV/figures/callPredictedBySVM/')
		CNV.CNVPredictionBySVM.getCNVPredictionErrorBySVMCrossValidation(db_250k, array_id=array_id, \
						minPercUnCoveredByLerContig=0.6, cnv_method_id=6, \
						replaceAmpWithMedianIntensity=False)
		
		cnv_method_id = 9
		for array_id in [3,4,41,150,151]:
			CNV.CNVPredictionBySVM.getCNVPredictionErrorBySVMCrossValidation(db_250k, array_id=array_id, \
					minPercUnCoveredByLerContig=0.6, cnv_method_id=cnv_method_id, \
					replaceAmpWithMedianIntensity=False)
		sys.exit(0)
		
		#2010-6-30
		output_dir = os.path.expanduser('~/script/variation/data/CNV/figures/callPredictedBySVM/')
		gridsize = 40
		cnv_method_id = 8
		for array_id in [3, 4, 41, 150, 151]:
			for minPercUnCoveredByLerContig in [0.4, 0.6, 0.8]:
				CNV.CNVPredictionBySVM.getCNVPredictionErrorBySVMCrossValidation(db_250k, output_dir=output_dir,\
						array_id=array_id, \
						minPercUnCoveredByLerContig=minPercUnCoveredByLerContig, cnv_method_id=cnv_method_id, \
						replaceAmpWithMedianIntensity=False, gridsize=gridsize)
		sys.exit(0)

		"""
		@classmethod
		def fitAndPredictBySVM(cls, param, training_feature_data, training_class_label_ls, test_feature_data=None, \
							test_class_label_ls=None, crossValidationNrFold=0):
			"""
			2010-6-30
				either do cross_validation using training_feature_data & training_class_label_ls
				
				or use training_* to train the model and then predict for the test_* data 
			"""
			import sys
			mute_device = open('/dev/null', 'w')
			sys.stderr = mute_device
			sys.stdout = mute_device
			from svm import svm_problem, svm_parameter, svm_model, cross_validation, LINEAR, POLY, RBF
			problem = svm_problem(training_class_label_ls, training_feature_data)
			isCrossValidationMode = False 
			if (crossValidationNrFold is None or crossValidationNrFold ==0) and (test_feature_data is not None \
				and test_class_label_ls is not None):
				model = svm_model(problem, param)
				class_label_ls = test_class_label_ls
				_prediction_ls = None
			elif crossValidationNrFold is not None and crossValidationNrFold >0:
				_prediction_ls = cross_validation(problem, param, crossValidationNrFold)
				class_label_ls = training_class_label_ls
				model = None
			else:
				return None
			sys.stderr = sys.__stderr__
			sys.stdout = sys.__stdout__
			
			errors = 0.
			no_of_cnvs = len(class_label_ls)
			prediction_ls = []
			probability_ls = []
			TP = 0	# True positive
			FP = 0
			FN = 0
			TN = 0
			for i in range(no_of_cnvs):
				probability = None
				if model is not None:
					if param.probability:
						prediction, probability = model.predict_probability(test_feature_data[i])
						probability_ls.append(probability[prediction])
					else:
						prediction = model.predict(test_feature_data[i])
				else:
					prediction = _prediction_ls[i]
				prediction_ls.append(prediction)
				"""
				# this is for continuous class label (mse estimate)
				errors += (prediction-class_label_ls[i])*(prediction-class_label_ls[i])
				"""
				if class_label_ls[i]==-1 and prediction==-1:
					TP += 1
				elif class_label_ls[i]==-1 and prediction==1:
					FN += 1
				elif class_label_ls[i]==1 and prediction==-1:
					FP += 1
				elif class_label_ls[i]==1 and prediction==1:
					TN += 1
				
				if (class_label_ls[i] != prediction):
					errors = errors + 1
			errors = errors/float(no_of_cnvs)
			precision = TP/float(FP+TP)
			sensitivity = TP/float(TP+FN)
			FPR = FP/float(FP+TN)
			print "##########################################"
			print " kernel %s, C %s, gamma %s eps %s: error rate = %s" % (param.kernel_type, param.C, param.gamma, param.eps, errors)
			print "##########################################"
			from pymodule import PassingData
			return PassingData(prediction_ls=prediction_ls, errors= errors, probability_ls=probability_ls,\
							TP=TP, FP=FP, FN=FN, TN=TN, precision=precision, sensitivity=sensitivity, FPR=FPR)
			
		@classmethod
		def applyOneArraySVMModelToAnotherArray(cls, db_250k, output_dir=None, array1_id=None, array2_id=None,\
						minPercUnCoveredByLerContig=0.6, cnv_method_id=6, \
						replaceAmpWithMedianIntensity=False, gridsize=40):
			"""
			2010-6-30
			"""
			array1FeatureData = cls.getCNVFeatureData(db_250k, array_id=array1_id, \
						minPercUnCoveredByLerContig=minPercUnCoveredByLerContig, cnv_method_id=cnv_method_id, \
						replaceAmpWithMedianIntensity=replaceAmpWithMedianIntensity)
			array2FeatureData = cls.getCNVFeatureData(db_250k, array_id=array2_id, \
						minPercUnCoveredByLerContig=minPercUnCoveredByLerContig, cnv_method_id=cnv_method_id, \
						replaceAmpWithMedianIntensity=replaceAmpWithMedianIntensity)
			if len(array1FeatureData.feature_data)<100 or len(array2FeatureData.feature_data)<100:
				return
			import sys
			from svm import svm_parameter, LINEAR, POLY, RBF
			
			total_correct = 0.
			kernels2kname = {LINEAR:'linear', POLY:'polynomial', RBF:'rbf'}
			kernels2kname = {POLY:'polynomial'}	# 2009- 11-18 only polynomial
			kernels2kname = {RBF:'rbf'}
			
			param = svm_parameter(C = 10, eps=1e-2, probability = 1)	#,nr_weight = 2,weight_label = [1,0],weight = [10,1])
			C_ls = [100, 10, 1, 0.1, 0.01]	# tried 1000, 10, 1, 0.1 ,0.01
			C_ls = [10]
			#gamma_ls = [100, 10, 1, 0.1, 0.0001, 0]
			gamma_ls = [0]	# test shows that gamma affects rbf quite a bit
			nr_fold = 10	# 10 fold cross validation
			for k, kname in kernels2kname.iteritems():
				for C in C_ls:
					for gamma in gamma_ls:
						#param.gamma = gamma
						param.C = C
						param.kernel_type = k;
						fitResult = cls.fitAndPredictBySVM(param, array1FeatureData.feature_data, \
											array1FeatureData.class_label_ls,\
											test_feature_data = array2FeatureData.feature_data, \
											test_class_label_ls = array2FeatureData.class_label_ls, \
											crossValidationNrFold=0)
						cls.drawSVMClassPrediction(output_dir, array2FeatureData.feature_data, fitResult.prediction_ls, \
								gridsize=gridsize, \
								extraFileNameTag = 'SVMModelFromArray%s'%array1_id, \
								extraTitle='M-a %s, PPV %.3f TP %s TPR %.3f FPR %.3f'%(array1_id, fitResult.precision, \
											fitResult.TP, \
											fitResult.sensitivity, fitResult.FPR),\
								class_label_name='segment label', \
								array_id=array2_id, ecotype_id=array2FeatureData.ecotype_id, cnv_method_id=cnv_method_id, \
								minPercUnCoveredByLerContig=minPercUnCoveredByLerContig)
						#fitResult.FP, fitResult.FN, fitResult.TN, 
		"""
		output_dir = os.path.expanduser('~/script/variation/data/CNV/figures/callPredictedBySVM/')
		cnv_method_id = 8
		gridsize = 40
		CNV.CNVPredictionBySVM.applyOneArraySVMModelToAnotherArray(db_250k, array1_id=3, array2_id=150,\
						minPercUnCoveredByLerContig=0.6, cnv_method_id=cnv_method_id, \
						replaceAmpWithMedianIntensity=False, output_dir=output_dir, gridsize=gridsize)
		sys.exit(0)
		
		# 2010-6-30
		output_dir = os.path.expanduser('~/script/variation/data/CNV/figures/callPredictedBySVM/')
		cnv_method_id = 8
		gridsize = 40
		for array1_id in [3,41,150,151]:
			for array2_id in [3,41,150,151]:
				for gridsize in [60, 100]:
					CNV.CNVPredictionBySVM.applyOneArraySVMModelToAnotherArray(db_250k, array1_id=array1_id, \
								array2_id=array2_id,\
								minPercUnCoveredByLerContig=0.6, cnv_method_id=cnv_method_id, \
								replaceAmpWithMedianIntensity=False, output_dir=output_dir, gridsize=gridsize)
		sys.exit(0)
		
		# 2010-6-30
		output_dir = os.path.expanduser('~/script/variation/data/CNV/figures/callPredictedBySVM/')
		cnv_method_id = 6
		gridsize = 40
		for array1_id in [3,41,150,151]:
			for array2_id in [3,41,150,151]:
				for gridsize in [40, 60, 100]:
					for minPercUnCoveredByLerContig in [0.4, 0.6, 0.8]:
						CNV.CNVPredictionBySVM.applyOneArraySVMModelToAnotherArray(db_250k, array1_id=array1_id, \
								array2_id=array2_id,\
								minPercUnCoveredByLerContig=minPercUnCoveredByLerContig, cnv_method_id=cnv_method_id, \
								replaceAmpWithMedianIntensity=False, output_dir=output_dir, gridsize=gridsize)
		sys.exit(0)
		
		
		"""
		@classmethod
		def drawSVMClassPrediction(cls, output_dir, feature_data, class_label_ls, \
								extraFileNameTag = '', extraTitle='',\
								feature_column_name_ls=['amplitude', 'no of probes', 'ProbeDensity'], \
								class_label_name='fraction uncovered by Contig Sequence', gridsize=40, \
								array_id=None, ecotype_id=None, cnv_method_id=None, minPercUnCoveredByLerContig=0.6):
			"""
			2010-7-1
				pass this to drawHexbin()
					reduce_C_function=numpy.mean
			2010-6-30
			"""
			import pylab
			import matplotlib
			from matplotlib import cm
			import numpy
			from pymodule.utils import addExtraLsToFilenamePrefix, addExtraToFilenamePrefix
			
			if len(class_label_ls)>100:
				feature_data = numpy.array(feature_data)
				sys.stderr.write(" gridSize: %s "%gridsize)
				if output_dir and not os.path.isdir(output_dir):	#2010-5-5 test if output_dir is something
					os.makedirs(output_dir)
				common_fig_fname = 'array_%s_ecotype_%s_cnvMethod%s.png'%(array_id, ecotype_id, cnv_method_id)
				if extraFileNameTag is not None:
					fig_fname = addExtraToFilenamePrefix(common_fig_fname, extraFileNameTag)
				"""
				fig_fname = addExtraToFilenamePrefix(common_fig_fname, 'AmplitudeAsX_ProbeDensityAsY_minUncoverPerc%s_gridsize_%s'%\
													(minPercUnCoveredByLerContig, gridsize))
				fig_fname = addExtraToFilenamePrefix(common_fig_fname, 'AmplitudeAsX_NoOfProbesAsY_minUncoverPerc%s_gridsize_%s'%\
													(minPercUnCoveredByLerContig, gridsize))
				"""
				fig_fname = addExtraToFilenamePrefix(fig_fname, 'MeanIntensityAsX_NoOfProbesAsY_minUncoverPerc%s_gridsize_%s'%\
													(minPercUnCoveredByLerContig, gridsize))
				fig_fname = os.path.join(output_dir, fig_fname)
				title = "array %s minPerc %s gsize %s %s"%(array_id, minPercUnCoveredByLerContig, gridsize, extraTitle)
				
				"""
				ylabel = 'No of probes per kb'
				xlabel = 'No of probes'
				max_phenotype = numpy.nanmax(percUnCoveredByLerContig_ls)
				min_phenotype = numpy.nanmin(percUnCoveredByLerContig_ls)
				phenotype_gap = max_phenotype - min_phenotype
				phenotype_jitter = phenotype_gap/10.
				norm = matplotlib.colors.Normalize(vmin=min_phenotype-phenotype_jitter, vmax=max_phenotype+phenotype_jitter)
				
				import pylab
				from mpl_toolkits.mplot3d import Axes3D
				fig = pylab.figure()
				ax = Axes3D(fig)
				ax.scatter(no_of_probes_ls, probeDensityLs, amplitude_ls,c=c_ls)
				pylab.title(title)
				pylab.xlabel(xlabel)
				pylab.ylabel(ylabel)
				pylab.show()
				"""
				"""
				cls.drawHexbin(amplitude_ls, probeDensityLs, percUnCoveredByLerContig_ls, fig_fname=fig_fname, \
					gridsize=gridsize, title=title, \
							xlabel='amplitude', \
							ylabel='no of probes per kb', colorBarLabel='percentage uncovered by Contig Sequence')
				"""
				CNV.drawHexbin(feature_data[:,0], feature_data[:, 1], class_label_ls, fig_fname=fig_fname, \
					gridsize=gridsize, title=title, \
					xlabel=feature_column_name_ls[0], \
					ylabel=feature_column_name_ls[1], colorBarLabel=class_label_name,\
					reduce_C_function=numpy.mean)
				sys.stderr.write("\n")
		
	class Normalize(object):
		@classmethod
		def outputKinshipVsIntensityCorrelation(cls, db_250k, kinship_output_fname=None, coefficient_fname=None, output_fname=None,\
											refEcotypeID='6909'):
			"""
			2010-5-25
				correlation is between one array and the median of all arrays
				
				investigate whether correlation is a good measure.
			"""
			from pymodule.SNP import SNPData
			from pymodule import figureOutDelimiter
			import numpy, os, sys, csv, Stock_250kDB
			kinshipData = SNPData(input_fname=kinship_output_fname, ignore_2nd_column=1, \
								matrix_data_type=float, turn_into_array=1)
			ecotype_id2kinship = {}
			row_index = kinshipData.row_id2row_index[refEcotypeID]
			no_of_cols = len(kinshipData.col_id_ls)
			for col_id, col_index in kinshipData.col_id2col_index.iteritems():
				if col_id!=refEcotypeID:
					ecotype_id2kinship[col_id] = kinshipData.data_matrix[row_index][col_index]
			
			kinship_ls = []
			coefficient_ls = []
			reader = csv.reader(open(coefficient_fname), delimiter=figureOutDelimiter(coefficient_fname))
			header = reader.next()
			for row in reader:
				coefficient = float(row[0])
				array_id = row[1]
				array = Stock_250kDB.ArrayInfo.get(array_id)
				ecotype_id = str(array.maternal_ecotype_id)
				if ecotype_id in ecotype_id2kinship:
					kinship_ls.append(ecotype_id2kinship.get(ecotype_id))
					coefficient_ls.append(coefficient)
			import pylab
			pylab.plot(kinship_ls, coefficient_ls , '.')
			pylab.xlabel('kinship distance to ecotype %s'%refEcotypeID)
			pylab.ylabel('correlation to median of all')
			pylab.savefig(output_fname, dpi=300)
		
		"""
	kinship_output_fname = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482_kinship.tsv')
	coefficient_fname = os.path.expanduser('~/Downloads/CorrelationBetweenOneArrayAndRefMedian.csv')
	refEcotypeID='6909'
	output_fname = os.path.expanduser('~/script/variation/data/CNV/figures/kinship_to_e%s_vs_correlation_coeff.png'%refEcotypeID)
	CNV.Normalize.outputKinshipVsIntensityCorrelation(db_250k, kinship_output_fname=kinship_output_fname, \
			coefficient_fname=coefficient_fname, output_fname=output_fname,\
											refEcotypeID=refEcotypeID)
	
	
	kinship_output_fname = os.path.expanduser('~/mnt/panfs/250k/dataset/call_method_49_core482_kinship.tsv')
	coefficient_fname = os.path.expanduser('~/Downloads/CorrelationBetweenOneArrayAndRefMedian.csv')
	refEcotypeID='6932'
	output_fname = os.path.expanduser('~/script/variation/data/CNV/figures/kinship_to_e%s_vs_correlation_coeff.png'%refEcotypeID)
	CNV.Normalize.outputKinshipVsIntensityCorrelation(db_250k, kinship_output_fname=kinship_output_fname, \
			coefficient_fname=coefficient_fname, output_fname=output_fname,\
											refEcotypeID=refEcotypeID)
	sys.exit(0)
		"""
		
		@classmethod
		def drawMedianIntensityVsIntensityCorrelation(cls, db_250k, coefficient_fname=None, output_fname=None,):
			"""
			2010-5-25
				correlation is between one array and the median of all arrays
				
				investigate whether correlation is a good measure.
			"""
			from pymodule.SNP import SNPData
			from pymodule import figureOutDelimiter
			import numpy, os, sys, csv, Stock_250kDB
			
			x_ls = []
			coefficient_ls = []
			reader = csv.reader(open(coefficient_fname), delimiter=figureOutDelimiter(coefficient_fname))
			header = reader.next()
			for row in reader:
				coefficient = float(row[0])
				array_id = row[1]
				array = Stock_250kDB.ArrayInfo.get(array_id)
				if array.median_intensity:
					x_ls.append(array.median_intensity)
					coefficient_ls.append(coefficient)
			import pylab
			pylab.plot(x_ls, coefficient_ls , '.')
			pylab.xlabel('median intensity')
			pylab.ylabel('correlation to median of all arrays')
			pylab.savefig(output_fname, dpi=300)
		"""
	coefficient_fname = os.path.expanduser('~/Downloads/CorrelationBetweenOneArrayAndRefMedian.csv')
	output_fname = os.path.expanduser('~/script/variation/data/CNV/figures/median_intensity_vs_correlation_coeff.png')
	CNV.Normalize.drawMedianIntensityVsIntensityCorrelation(db_250k, \
			coefficient_fname=coefficient_fname, output_fname=output_fname)
	sys.exit(0)
	
		"""
		@classmethod
		def extendLsByResampling(cls, data_ls, max_no_of_points):
			"""
			2010-5-19
				
			"""
			import random
			# now to make up data in the block if its number of data points is < max_no_of_points_in_one_block
			no_of_points = len(data_ls)
			if no_of_points<max_no_of_points:
				no_of_points_to_add = max_no_of_points - no_of_points
				no_of_copies = int(no_of_points_to_add/no_of_points)
				original_block_data = data_ls[:]	#[:] is to avoid reference. deep copy.
				data_ls.extend(original_block_data*no_of_copies)
				no_of_points_to_sample = no_of_points_to_add%no_of_points
				if no_of_points_to_sample>0:	# do a sampling
					data_ls.extend(random.sample(original_block_data, no_of_points_to_sample))
			return data_ls
		
		@classmethod
		def quantileNormalizeAllBlocks(cls, blockData_ls, blockDataCodedIndex_ls, probes, qnorm_data_matrix, array_index=None,
									returnBlockMatrix=False, target_probeXblock_matrix=None, common_fig_fname=None, \
									refArrayID=None, targetArrayID=None):
			"""
			2010-5-19
			"""
			import numpy
			from CNVNormalize import CNVNormalize
			data_matrix = numpy.array(blockData_ls, numpy.float32).transpose()	#tranpose() is a must.
			"""
			quantile_normalize() takes each row as a probe. each column is an array/block.
			"""
			CNVNormalize.quantile_normalize(data_matrix)
			
			# aggregate data points of one probe
			codedIndex2QNormDataLs = {}
			no_of_blocks = len(blockData_ls)
			for i in range(no_of_blocks):
				blockDataCodedIndex = blockDataCodedIndex_ls[i]
				no_of_probes_in_block = len(blockDataCodedIndex)
				if target_probeXblock_matrix is not None and common_fig_fname:
					# ref vs target
					fit_y_ls = cls.ltsFit(data_matrix[:no_of_probes_in_block,i], target_probeXblock_matrix[:no_of_probes_in_block,i])
					data_matrix[:no_of_probes_in_block,i] = target_probeXblock_matrix[:no_of_probes_in_block,i] - fit_y_ls
					"""
					from pymodule.utils import addExtraLsToFilenamePrefix, addExtraToFilenamePrefix
					ref_target_fig_fname = addExtraToFilenamePrefix(common_fig_fname, \
													'ref_target_block%s_after_qnorm_by_ref_%s'%(i, refArrayID))
					cls.plotXY(data_matrix[:no_of_probes_in_block,i], target_probeXblock_matrix[:no_of_probes_in_block,i], \
								title="block %s"%repr(blockDataCodedIndex[0]), output_fname=ref_target_fig_fname, \
								xlabel="array %s"%refArrayID, ylabel="array %s"%targetArrayID,)
					"""
				for j in range(no_of_probes_in_block):
					xy_pos = blockDataCodedIndex[j]
					if xy_pos not in codedIndex2QNormDataLs:
						codedIndex2QNormDataLs[xy_pos] = []
					intensity = data_matrix[j,i]
					codedIndex2QNormDataLs[xy_pos].append(intensity)
			# get the median of all qnormed intensity for one probe
			#codedIndex2Median = {}
			for xy_pos, data_ls in codedIndex2QNormDataLs.iteritems():
				row_index = probes.xy2index.get(xy_pos)
				#codedIndex2Median[xy_pos] = numpy.median(data_ls)
				qnorm_data_matrix[row_index,array_index] = numpy.median(data_ls)
			if returnBlockMatrix:
				return data_matrix
			else:
				del data_matrix
		
		
		@classmethod
		def generateBlocksWithinOneArray(cls,array_width, blockSize=50, jumpStep=10,
										x_range=None, y_range=None):
			"""
			2010-5-17
				split out of subArrayQuantileNormalize()
			"""
			sys.stderr.write("Generating block outlines...")
			blockBottomLeftXY_ls = []
			import math
			if x_range:
				start_x = x_range[0]
				array_x_max = min(array_width, x_range[1]+1)
				array_x_span = x_range[1]-x_range[0] + 1
			else:
				start_x = 0
				array_x_max = array_width
				array_x_span = array_width
			if y_range:
				_start_y = y_range[0]
				array_y_max = min(array_width, y_range[1]+1)
				array_y_span = y_range[1]-y_range[0] + 1
			else:
				_start_y = 0
				array_y_max = array_width
				array_y_span = array_width
			no_of_blocks_in_x = int(math.ceil((array_x_span-blockSize)/float(jumpStep)))+1
			no_of_blocks_in_y = int(math.ceil((array_y_span-blockSize)/float(jumpStep)))+1
			for i in range(no_of_blocks_in_x):
				if start_x+blockSize>array_x_max:	#reduce the x to maintain the same block size
					start_x = array_x_max-blockSize
				start_y = _start_y
				for j in range(no_of_blocks_in_y):
					if start_y+blockSize>array_y_max:	#reduce the y to maintain the same block size
						start_y = array_y_max-blockSize
					blockBottomLeftXY_ls.append((start_x,start_y))
					start_y += jumpStep
				start_x += jumpStep
			sys.stderr.write("%s blocks. Done.\n"%(len(blockBottomLeftXY_ls)))
			return blockBottomLeftXY_ls
		
		@classmethod
		def ltsFit(cls, x_ls, y_ls, fractionUsed=0.6):
			"""
			2010-5-21
				use ROOT to fit the y=a+bx
			"""
			from pymodule.algorithm import ltsFit
			return ltsFit(x_ls, y_ls, fractionUsed=fractionUsed)
		
		"""
	import numpy
	x_ls = numpy.array(range(100), numpy.float)
	y_ls = x_ls/2.
	for i in range(len(y_ls)):
		import random
		new_y = random.random()-0.5
		y_ls[i] += new_y
	
	# mess up some portion of y
	for i in range(5):
		import random
		new_y = random.random()
		new_y_index = random.sample(range(100),1)
		y_ls[new_y_index[0]] = new_y
	import numpy
	x_ls = numpy.array([ 2.64884758,  3.51235008,  2.83090925,  3.41229248,  3.01451969,\
    2.49899888,  3.69988108,  2.74896216,  3.05307841,  3.75705409,\
    3.08653784,  3.10703993,  3.61071348,  3.21285319,  2.91460752,\
    3.53737831,  3.06333303,  3.35391617,  3.43568516,  3.34429312,\
    3.31576061,  2.8007164 ,  2.73639655,  3.14690256,  3.10174704,\
    2.80888581,  2.72754121,  2.90064001,  3.19270658,  3.50596333,\
    2.61804676,  3.18127131,  3.27542663,  3.09586573], dtype=numpy.float)	# numpy.float32 is not supported by ROOT
	y_ls = numpy.array([ 2.52827311,  3.27265358,  2.36172366,  2.95760489,  2.50920248,\
    2.3443923 ,  3.23502254,  2.35410833,  2.50582743,  2.48501062,\
    2.82510138,  2.70799541,  2.43136382,  2.76342535,  2.45178652,\
    3.08224201,  2.26481771,  2.7387805 ,  3.23274207,  2.82769203,\
    2.25042009,  2.56702638,  2.4082365 ,  2.44793224,  2.65127802,\
    2.57460976,  2.43136382,  2.39005065,  2.70027065,  3.04452848,\
    2.28555727,  2.71933126,  2.6468935 ,  2.54157925], dtype=numpy.float)
    
	fitResult = CNV.Normalize.ltsFit(x_ls, y_ls)
	
	import pylab
	pylab.plot(x_ls, y_ls, '.')
	pylab.plot(x_ls, fitResult.fit_y_ls, '.')
	pylab.legend(['raw data','fitted'])
	pylab.show()
	sys.exit(0)
		
		"""
		
		@classmethod
		def subArrayQuantileNormalize(cls, db_250k, array_id_ls, output_fname_prefix, refArrayIDLs=[], output_dir=None, \
									blockSize=50, jumpStep=10,
				probeType=2, x_range=None, y_range=None, minNoOfProbesPerBlock=40, draw2D=True, gridsize=200, drawHist=True,
				outlierBlockPercentile=1.0):
			"""
			2010-5-18
				add argument refArrayIDLs to offer reference
			2010-5-9
				split the whole array into overlapping blocks, then do quantile normalization over all blocks
				
				# combine all blocks together into a 2D array, for probes with more than 1 value, take the median
				
				probeType=2 means only CNV probes
			"""
			
			from DB_250k2Array import DB_250k2Array, probes_class
			import Stock_250kDB, numpy, random, rpy
			from CNVNormalize import CNVNormalize
			rpy.r.library('affy')
			array_1 = Stock_250kDB.ArrayInfo.get(1)
			returnData = DB_250k2Array.getArrayWidth(array_1.filename)
			array_width = returnData.array_width
			
			blockBottomLeftXY_ls = cls.generateBlocksWithinOneArray(array_width, blockSize=blockSize, jumpStep=jumpStep,
										x_range=x_range, y_range=y_range)
			
			probes, xy_ls, chr_pos_ls, probe_id_ls = DB_250k2Array.get_probes(db_250k.metadata.bind, \
										Stock_250kDB.Probes.table.name, snps=None, run_type=probeType, \
										x_range=x_range,\
										y_range=y_range, constructChrPos2index=False, constructXY2index=True)
			
			sys.stderr.write("Generating xy indices for each block ...")
			no_of_blocks = len(blockBottomLeftXY_ls)
			xy_set = set(xy_ls)
			# construct a list which stores a list of (x,y) in each block
			new_xy_set = set()	# some probes will be excluded if they are situated in blocks with too few probes
			blockDataCodedIndex_ls = []
			for k in range(no_of_blocks):
				bottomLeftXY = blockBottomLeftXY_ls[k]
				start_x, start_y = bottomLeftXY
				blockDataCodedIndex = []
				for x in xrange(start_x, start_x+blockSize):
					for y in xrange(start_y, start_y+blockSize):
						xy_pos = (x, y)
						if xy_pos in xy_set:	#make sure it's only the probes that we want
							blockDataCodedIndex.append(xy_pos)
				if len(blockDataCodedIndex)>=minNoOfProbesPerBlock:	# skip blocks that are too small
					blockDataCodedIndex_ls.append(blockDataCodedIndex)
					new_xy_set |= set(blockDataCodedIndex)
			sys.stderr.write("%s probes in %s blocks. Done.\n"%(len(new_xy_set), len(blockDataCodedIndex_ls)))
			
			sys.stderr.write("Generating new probes ...")
			new_probe_ls = [(probes.getProbeGivenXY(xy[0], xy[1]))for xy in new_xy_set]
			del probes
			new_probe_ls.sort()		# in chr,pos order
			x_ls = []	# for hexbin drawing
			y_ls = []
			new_probes = probes_class(constructChrPos2index=False, constructXY2index=True)
			for probe in new_probe_ls:
				x_ls.append(probe.xpos)
				y_ls.append(probe.ypos)
				new_probes.addOneProbe(probe)
			sys.stderr.write("Done.\n")
			
			if output_dir and not os.path.isdir(output_dir):	#2010-5-5 test if output_dir is something
				os.makedirs(output_dir)
			
			from pymodule.utils import addExtraLsToFilenamePrefix, addExtraToFilenamePrefix
			sys.stderr.write("sub-array quantile-normalize ...\n")
			no_of_arrays = len(array_id_ls)
			max_no_of_points_in_one_block = max([len(row) for row in blockDataCodedIndex_ls])
			qnorm_data_matrix = numpy.zeros([len(new_probes), no_of_arrays], numpy.float32)
			qnorm_data_matrix[:,:] = numpy.nan
			ref_qnorm_data_matrix = numpy.zeros([len(new_probes), len(refArrayIDLs)], numpy.float32)
			ref_qnorm_data_matrix[:,:] = numpy.nan
			for k in range(no_of_arrays):
				array_id = array_id_ls[k]
				sys.stderr.write("%s/%s: array %s  "%(k+1, no_of_arrays, array_id))
				if len(refArrayIDLs)>k:
					refArrayID = refArrayIDLs[k]
					refArrayInfo = Stock_250kDB.ArrayInfo.get(int(refArrayID))
					refArray = rpy.r.read_affybatch(filenames=refArrayInfo.filename)
					ref_intensity_array = rpy.r.intensity(refArray)
				blockData_ls = []
				refBlockData_ls = []
				arrayInfo = Stock_250kDB.ArrayInfo.get(int(array_id))
				array = rpy.r.read_affybatch(filenames=arrayInfo.filename)
				intensity_array = rpy.r.intensity(array)
				blockDataSigma_ls = []
				no_of_blocks = len(blockDataCodedIndex_ls)	# no_of_blocks has changed
				
				common_fig_fname = 'array_%s_subArrayQNorm_b%s_j%s_m%s_obp%s.png'%\
					(array_id, blockSize, jumpStep, minNoOfProbesPerBlock, outlierBlockPercentile)
				if x_range:
					common_fig_fname = addExtraLsToFilenamePrefix(common_fig_fname, ['x_range']+x_range)
				if y_range:
					common_fig_fname = addExtraLsToFilenamePrefix(common_fig_fname, ['y_range']+y_range)
				
				for i in range(no_of_blocks):
					blockDataCodedIndex = blockDataCodedIndex_ls[i]
					blockData_ls.append([])
					refOneBlockIntensityLs = []
					for xy_pos in blockDataCodedIndex:
						xpos, ypos = xy_pos
						intensity_array_index = array_width*(array_width - xpos - 1) + ypos
						intensity = math.log10(intensity_array[intensity_array_index][0])
						refIntensity = math.log10(ref_intensity_array[intensity_array_index][0])
						blockData_ls[i].append(intensity)
						refOneBlockIntensityLs.append(refIntensity)
					# now to make up data in the block if its number of points is < max_no_of_points_in_one_block
					blockData_ls[i] = cls.extendLsByResampling(blockData_ls[i], max_no_of_points_in_one_block)
					no_of_probes_in_ref_block = len(refOneBlockIntensityLs)
					if no_of_probes_in_ref_block>0:
						refBlockData_ls.append(cls.extendLsByResampling(refOneBlockIntensityLs, max_no_of_points_in_one_block))
					
					sigmaEstimate = numpy.std(blockData_ls[i])/numpy.sqrt(max_no_of_points_in_one_block-1)
					meanEstimate = numpy.mean(blockData_ls[i])
					blockDataSigma_ls.append(sigmaEstimate/meanEstimate)
					if no_of_probes_in_ref_block>0 and i%3==0:
						ref_target_fig_fname = addExtraToFilenamePrefix(common_fig_fname, \
																	'ref_target_block%s_before_qnorm_by_ref_%s'%(i, refArrayID))
						ref_target_fig_fname = os.path.join(output_dir, ref_target_fig_fname)
						CNV.plotXY(refOneBlockIntensityLs[:no_of_probes_in_ref_block], blockData_ls[i][:no_of_probes_in_ref_block], \
								title="block %s"%repr(blockDataCodedIndex[0]), output_fname=ref_target_fig_fname, \
								xlabel="array %s"%refArrayID, ylabel="array %s"%array_id,)
				if outlierBlockPercentile<1.0 and outlierBlockPercentile>0.0:
					# only take blocks whose sigma is in bottom outlierBlockPercentile, i.e. 90%
					blockDataSigma_argsort_ls = numpy.argsort(blockDataSigma_ls)
					maxRank =  int(no_of_blocks*outlierBlockPercentile)
					new_blockData_ls = []
					new_blockDataCodedIndex_ls = []
					new_blockDataSigma_ls = []
					new_refBlockData_ls = []
					for i in xrange(maxRank):
						old_block_index = blockDataSigma_argsort_ls[i]
						new_blockDataSigma_ls.append(blockDataSigma_ls[old_block_index])
						new_blockData_ls.append(blockData_ls[old_block_index])
						if len(refBlockData_ls)>old_block_index:
							new_refBlockData_ls.append(refBlockData_ls[old_block_index])
						
						new_blockDataCodedIndex_ls.append(blockDataCodedIndex_ls[old_block_index])
				else:
					new_blockData_ls = blockData_ls
					new_blockDataCodedIndex_ls = blockDataCodedIndex_ls
					new_blockDataSigma_ls = blockDataSigma_ls
					new_refBlockData_ls = refBlockData_ls
				new_no_of_blocks = len(new_blockData_ls)
				sys.stderr.write(" %s blocks -> %s blocks "%(no_of_blocks, new_no_of_blocks))
				
				if drawHist:
					hist_fig_fname = addExtraToFilenamePrefix(common_fig_fname, 'blockDataSigmaHist')
					hist_fig_fname = os.path.join(output_dir, hist_fig_fname)
					from PlotQCProbeIntensityHistogram import PlotQCProbeIntensityHistogram
					PlotQCProbeIntensityHistogram.plotHistogram(new_blockDataSigma_ls, array_id, hist_fig_fname, \
								no_of_bins=min(new_no_of_blocks/4, 50), \
								max_intensity=None, \
								xlim=None)
				
				target_probeXblock_matrix = cls.quantileNormalizeAllBlocks(new_blockData_ls, new_blockDataCodedIndex_ls, new_probes, \
											qnorm_data_matrix, array_index=k, returnBlockMatrix=True)
				if new_refBlockData_ls:
					ref_target_common_fig_fname = os.path.join(output_dir, common_fig_fname)
					cls.quantileNormalizeAllBlocks(new_refBlockData_ls, new_blockDataCodedIndex_ls, new_probes, \
										ref_qnorm_data_matrix, array_index=k, target_probeXblock_matrix=target_probeXblock_matrix, \
										common_fig_fname=ref_target_common_fig_fname, \
										refArrayID=refArrayID, targetArrayID=array_id)
					qnorm_data_matrix = ref_qnorm_data_matrix
				if draw2D and output_dir:
					fig_fname = addExtraToFilenamePrefix(common_fig_fname, '2D_gridsize_%s'%(gridsize))
					fig_fname = os.path.join(output_dir, fig_fname)
					title = "intensity of array %s, blockSize %s, step %s. %s probes."%(array_id, blockSize, jumpStep, len(x_ls))
					CNV.drawHexbin(x_ls, y_ls, qnorm_data_matrix[:,k], fig_fname=fig_fname, gridsize=gridsize, title=title)
				
				sys.stderr.write(".\n")
			
			header =['probes_id'] + array_id_ls + ['chromosome', 'position']
			array_id_str_ls = map(str, array_id_ls)
			output_fname = '%s_%s_b%s_j%s_m%s_obp%s.tsv'%(output_fname_prefix, '_'.join(array_id_str_ls), blockSize, \
														jumpStep, minNoOfProbesPerBlock, outlierBlockPercentile)
			CNVNormalize.output(qnorm_data_matrix, new_probes.get_probe_id_ls(), new_probes.get_chr_pos_ls(), header, output_fname)
			
		
		"""
	# 2010-5-11
	output_fname_prefix = os.path.expanduser('~/mnt/panfs/250k/CNV/call_48_CNV_subArrayQNorm')
	output_dir = os.path.expanduser('~/script/variation/data/CNV/normalization/figures/')
	array_id_ls = [1, 2, 43, 139, 145, 3, 4, 41, 150, 151]	# 5 Col and 5 Ler arrays
	CNV.Normalize.subArrayQuantileNormalize(db_250k, array_id_ls, output_fname_prefix, output_dir=output_dir, blockSize=50, jumpStep=10,
			probeType=2, x_range=None, y_range=None, minNoOfProbesPerBlock=40, draw2D=True, gridsize=200)
	
	# 2010-5-11
	output_fname_prefix = os.path.expanduser('~/mnt/panfs/250k/CNV/call_48_CNV_subArrayQNorm')
	output_dir = os.path.expanduser('~/script/variation/data/CNV/normalization/figures/')
	array_id_ls = [1, 2, 43, 139, 145, 3, 4, 41, 150, 151]	# 5 Col and 5 Ler arrays
	array_id_ls= [41,4]
	y_range = None	#[1200, 1612]	#check small range to see if it works
	x_range = None	#[0, 100]
	
	for blockSize in [100,]:
		jumpStep = blockSize/4
		minNoOfProbesPerBlock=blockSize*blockSize/6
		CNV.Normalize.subArrayQuantileNormalize(db_250k, array_id_ls, output_fname_prefix, output_dir=output_dir, blockSize=blockSize,\
								jumpStep=jumpStep,\
			probeType=2, x_range=x_range, y_range=y_range, minNoOfProbesPerBlock=minNoOfProbesPerBlock, \
			draw2D=True, gridsize=200, outlierBlockPercentile=0.85)
	sys.exit(0)
	
	# 2010-5-11
	output_fname_prefix = os.path.expanduser('~/mnt/panfs/250k/CNV/call_48_CNV_subArrayQNorm')
	output_dir = os.path.expanduser('~/script/variation/data/CNV/normalization/figures/')
	array_id_ls = [1, 2, 43, 139, 145, 3, 4, 41, 150, 151]	# 5 Col and 5 Ler arrays
	array_id_ls= [4,41,]
	array_id_ls = [3, 4, 41, 150, 151]
	refArrayIDLs = [1, 2 ,43, 139, 145]
	
	array_id_ls = [4,41,150]
	refArrayIDLs = [1, 1, 1, ]
	y_range = [1200, 1612]	#check small range to see if it works
	x_range = [0, 200]
	
	for blockSize in [100,]:
		jumpStep = blockSize/4
		minNoOfProbesPerBlock=blockSize*blockSize/6
		CNV.Normalize.subArrayQuantileNormalize(db_250k, array_id_ls, output_fname_prefix, refArrayIDLs=refArrayIDLs,\
								output_dir=output_dir, blockSize=blockSize,\
								jumpStep=jumpStep,\
			probeType=2, x_range=x_range, y_range=y_range, minNoOfProbesPerBlock=minNoOfProbesPerBlock, \
			draw2D=True, gridsize=200, outlierBlockPercentile=1.0)
	sys.exit(0)
	
		"""
		
		@classmethod
		def drawIntensityPlotRefVsTarget(cls, db_250k, input_prefix=None, output_dir=None, \
									gridsize_ls=[200], probeType=2, x_range=None,\
									y_range=None, targetArrayIDLs = [3, 4, 41, 150, 151], \
									ref_id_ls = [1, 2 ,43, 139, 145]):
			"""
			2010-5-23
				defunct, replace by drawMultiArrayIntensity2DFromLeiLiNormalizationResult()
				dont' know why it is here.
			2010-5-7
				ref_id_ls is the corresponding reference array used in normalization for each in input_fname_ls
			"""
			
			from DB_250k2Array import DB_250k2Array
			import Stock_250kDB, numpy, rpy
			from CNVNormalize import CNVNormalize
			probes, xy_ls, chr_pos_ls, probe_id_ls = DB_250k2Array.get_probes(db_250k.metadata.bind, \
										Stock_250kDB.Probes.table.name, snps=None, run_type=probeType, \
										x_range=x_range,\
										y_range=y_range, constructChrPos2index=False, constructXY2index=True)
			
			if output_dir and not os.path.isdir(output_dir):	#2010-5-5 test if output_dir is something
				os.makedirs(output_dir)
			
			no_of_arrays = len(targetArrayIDLs)
			data_matrix = numpy.zeros([len(probes), no_of_arrays], numpy.float32)
			ref_id2intensity_ls = {}
			for i in range(no_of_arrays):
				targetArrayID = targetArrayIDLs[i]
				input_fname_ls = []
				ref_fname_ls = []
				for ref_id in ref_id_ls:
					input_fname = os.path.join(input_prefix, 't_%s_raw_data_By_Ref_r_%s_raw_data.cel'%(targetArrayID, ref_id))
					input_fname_ls.append(input_fname)
					
					ref_fname = os.path.join(input_prefix, 'r_%s_raw_data_By_Ref_r_%s_raw_data.cel'%(ref_id, ref_id))
					ref_fname_ls.append(ref_fname)
				intensityOfOneArray = CNV.drawArrayIntensity2DFromLeiLiNormalizationResult(input_fname_ls, output_dir, ref_fname_ls, \
															gridsize_ls=gridsize_ls, probeType=probeType, array_id=targetArrayID, \
															probes=probes,\
															x_range=x_range, y_range=y_range, ref_id2intensity_ls=ref_id2intensity_ls)
				data_matrix[:,i] = intensityOfOneArray
			header =['probes_id'] + targetArrayIDLs + ['chromosome', 'position']
			
			
			from pymodule.utils import addExtraLsToFilenamePrefix, addExtraToFilenamePrefix
			targetArrays = "_".join(map(str, targetArrayIDLs))
			refArrays = "_".join(map(str, ref_id_ls))
			output_fname = 'subArrayQuantileNorm_%s_MedianAmong_Ref_%s_probeType_%s.tsv'%(targetArrays, refArrays, probeType)
			if x_range:
				output_fname = addExtraLsToFilenamePrefix(output_fname, ['x_range']+x_range)
			if y_range:
				output_fname = addExtraLsToFilenamePrefix(output_fname, ['y_range']+y_range)
			output_fname = os.path.join(output_dir, output_fname)
			CNVNormalize.output(data_matrix, probe_id_ls, chr_pos_ls, header, output_fname)
		
		"""
		"""
		
		@classmethod
		def removeWaveEffect(cls, input_fname, output_fname_prefix, db_250k=None, no_of_probes_to_span=200, \
							min_intensity=2.5, max_intensity=3.3, ecotype_id_ls=[]):
			"""
			2010-2-9
				input_fname is CNVNormalize.py input/output format but only one chromosome.
				
				use loessFit from R package limma to fit the long-range wave effect along the chromosome
				
			"""
			sys.stderr.write("Removing wave effect for %s ...\n"%input_fname)
			from CNVNormalize import CNVNormalize
			data_matrix, probe_id_ls, chr_pos_ls, header = CNVNormalize.get_input(input_fname)
			array_id_ls = header[1:-2]
			array_id_ls = map(int, array_id_ls)
			no_of_arrays = len(array_id_ls)
			if db_250k and ecotype_id_ls:
				ecotype_id_set = set(ecotype_id_ls)
				col_index_to_check_ls = []
				for i in range(no_of_arrays):
					array_id = array_id_ls[i]
					array = Stock_250kDB.ArrayInfo.get(array_id)
					if array.maternal_ecotype_id==array.paternal_ecotype_id and array.maternal_ecotype_id in ecotype_id_set:
						col_index_to_check_ls.append(i)
			else:
				col_index_to_check_ls = range(no_of_arrays)
			
			no_of_arrays_to_fit = len(col_index_to_check_ls)
			sys.stderr.write("wave correction for %s arrays ... \n"%(no_of_arrays_to_fit))
			import rpy, numpy
			rpy.r.library('limma')
			
			position_ls = [int(row[1]) for row in chr_pos_ls]
			no_of_rows, no_of_cols = data_matrix.shape
			
			wave_data_matrix = numpy.zeros([no_of_rows, no_of_arrays_to_fit], numpy.float32)
			
			for k in range(no_of_arrays_to_fit):
				j = col_index_to_check_ls[k]
				intensity_ls = data_matrix[:,j]
				array_id = array_id_ls[j]
				sys.stderr.write("array %s ..."%array_id)
				weight_ls = []
				for i in range(no_of_rows):
					if data_matrix[i,j]<min_intensity or data_matrix[i,j]>max_intensity:
						weight = 0
					else:
						weight = 1
					weight_ls.append(weight)
				no_of_valid_points = sum(weight_ls)
				if no_of_valid_points<no_of_probes_to_span:
					sys.stderr.write("Only %s points for array %s within the range [%s, %s]. Ignored.\n"%(no_of_valid_points, \
																	array_id, min_intensity, max_intensity))
					continue
				loessFitSpan = no_of_probes_to_span/float(no_of_valid_points)
				
				lft = rpy.r.loessFit(intensity_ls, position_ls, weights=weight_ls, span=loessFitSpan)
				wave_data_matrix[:,k] = lft["fitted"]
				data_matrix[:,j] = data_matrix[:,j]-wave_data_matrix[:,k]
				sys.stderr.write("Done.\n")
			if no_of_arrays_to_fit==no_of_arrays:
				new_header = header
			else:
				new_header = header[:1]
				for j in col_index_to_check_ls:
					new_header.append(header[j+1])
				new_header.extend(header[-2:])
				data_matrix = data_matrix[:,col_index_to_check_ls]
			
			corrected_output_fname = '%s_wave_corrected.tsv'%output_fname_prefix
			CNVNormalize.output(data_matrix, probe_id_ls, chr_pos_ls, new_header, corrected_output_fname)
			wave_output_fname = '%s_wave.tsv'%output_fname_prefix
			CNVNormalize.output(wave_data_matrix, probe_id_ls, chr_pos_ls, new_header, wave_output_fname)
			sys.stderr.write("Done.\n")
		"""
	# 2010-2-10
	
	# test
	input_fname = '/Network/Data/250k/tmp-yh/CNV/call_method_17_CNV_array_intensity_norm_chr3_head_n10000.tsv'
	output_fname_prefix = '/tmp/call_method_17_CNV_array_intensity_norm_chr3_head_n10000'
	CNV.removeWaveEffect( input_fname, output_fname_prefix, no_of_probes_to_span=100, min_intensity=-2, max_intensity=+2)
	
	input_fname = '/Network/Data/250k/tmp-yh/CNV/call_method_17_CNV_array_intensity_norm_chr3_head_n10000.tsv'
	output_fname_prefix = '/tmp/call_method_17_CNV_array_intensity_norm_chr3_head_n10000'
	CNV.removeWaveEffect( input_fname, output_fname_prefix, db_250k, no_of_probes_to_span=100, min_intensity=-2, max_intensity=+2,\
						ecotype_id_ls=[6909, 6932, 6911, 6977, 8215, 6962])
	
	# real dataset
	for i in range(1,6):
		input_fname = os.path.expanduser('~/mnt/panfs/250k/CNV/call_48_CNV_QNormalize_chr%s.tsv'%i)
		output_fname_prefix = os.path.expanduser('~/mnt/panfs/250k/CNV/call_48_CNV_QNormalize_chr%s'%i)
		CNV.removeWaveEffect(input_fname, output_fname_prefix, db_250k, no_of_probes_to_span=200, min_intensity=2.5, max_intensity=3.3,
							ecotype_id_ls=[6909, 6932, 6911, 6977, 8215, 6962])
	
	min_intensity = 1.5
	max_intensity = 4.0
	for i in range(1,6):
		input_fname = os.path.expanduser('~/panfs/250k/CNV/call_48_CNV_QNormalize_chr%s.tsv'%i)
		output_fname_prefix = os.path.expanduser('~/panfs/250k/CNV/call_48_CNV_QNormalize_chr%s_intensity_%s_%s'%\
												(i, min_intensity, max_intensity))
		CNV.removeWaveEffect(input_fname, output_fname_prefix, db_250k, no_of_probes_to_span=200, min_intensity=min_intensity, \
							max_intensity=max_intensity,
							ecotype_id_ls=[6909, 6932, 6911, 6977, 8215, 6962])
	
	# 2010-3-2
	min_intensity = 2.7
	max_intensity = 3.1
	no_of_probes_to_span=2000
	for i in range(1,6):
		input_fname = os.path.expanduser('~/panfs/250k/CNV/call_48_CNV_QNormalize_chr%s.tsv'%i)
		output_fname_prefix = os.path.expanduser('~/panfs/250k/CNV/call_48_CNV_QNormalize_chr%s_intensity_%s_%s_span_%s_probes'%\
												(i, min_intensity, max_intensity, no_of_probes_to_span))
		CNV.removeWaveEffect(input_fname, output_fname_prefix, db_250k, no_of_probes_to_span=no_of_probes_to_span, min_intensity=min_intensity, \
							max_intensity=max_intensity)
	
		"""
		
		@classmethod
		def plotSpatialWaveEffect(cls, db_250k, input_fname_ls, output_dir):
			"""
			2010-3-4
				Take the waves (loess fitting results) and plot them for all arrays.
				
				input_fname_ls is the output from removeWaveEffect().
			"""
			import os,sys
			
			from common import get_ecotypeid2nativename
			ecotypeid2nativename = get_ecotypeid2nativename(db_250k.metadata.bind)
			
			#sys.stderr.write("Plotting spatial wave effect ...\n")
			if not os.path.isdir(output_dir):
				os.makedirs(output_dir)
			
			from CNVNormalize import CNVNormalize
			import Stock_250kDB
			import pylab
			
			for input_fname in input_fname_ls:
				data_matrix, probe_id_ls, chr_pos_ls, header = CNVNormalize.get_input(input_fname)
				array_id_ls = header[1:-2]
				array_id_ls = map(int, array_id_ls)
				no_of_arrays = len(array_id_ls)
				for i in range(no_of_arrays):
					chr = chr_pos_ls[0][0]
					array_id = array_id_ls[i]
					array = Stock_250kDB.ArrayInfo.get(array_id)
					fname = os.path.split(input_fname)[1]
					fname_prefix = os.path.splitext(fname)[0]
					ecotype_id = array.maternal_ecotype_id
					nativename = ecotypeid2nativename.get(ecotype_id)
					
					sys.stderr.write('Plotting wave effect for chr %s, array %s, ecotype %s %s ...'%\
									(chr, array_id, ecotype_id, nativename))
					
					output_fname = os.path.join(output_dir, \
											'%s_ecotype_%s_%s_array_%s.png'%(fname_prefix, ecotype_id, nativename, array_id))
					
					pos_ls = [row[1] for row in chr_pos_ls]
					pos_ls = map(int, pos_ls)
					
					pylab.clf()
					pylab.xlabel('chr %s position'%chr)
					pylab.ylabel('intensity')
					pylab.title("Spatial wave effect ecotype %s %s array %s"%(ecotype_id, nativename, array_id))
					pylab.plot(pos_ls, data_matrix[:, i])
					pylab.savefig(output_fname, dpi=300)
					pylab.clf()
					sys.stderr.write("\n")
			sys.stderr.write("Done.\n")
			
		"""
	# 2010-3-4
	
	min_intensity = 2.7
	max_intensity = 3.1
	no_of_probes_to_span=2000
	input_fname_ls = []
	for i in range(1,6):
		input_fname = os.path.expanduser('~/mnt/panfs/250k/CNV/ColLerCviVanFeiSha/call_48_CNV_QNormalize_chr%s_intensity_%s_%s_span_%s_probes_wave.tsv'%\
												(i, min_intensity, max_intensity, no_of_probes_to_span))
		input_fname_ls.append(input_fname)
	
	output_dir = '/Network/Data/250k/tmp-yh/CNV/WaveEffect'
	CNV.plotSpatialWaveEffect(db_250k, input_fname_ls, output_dir)
		"""
		
		@classmethod
		def medianSmooth(cls, input_fname, output_fname, noOfFlankingProbes=2):
			"""
			2010-5-28
				input_fname is Array X Probe format
			"""
			sys.stderr.write("Median-smoothing for %s ...\n"%input_fname)
			import csv
			from pymodule.CNV import ArrayXProbeFileWrapper
			inputWrapper = ArrayXProbeFileWrapper(input_fname)
			
			reader = inputWrapper.reader
			probe_id_ls = inputWrapper.probe_id_ls
			chr_pos_ls = inputWrapper.chr_pos_ls
			
			no_of_probes = len(probe_id_ls)
			chr2start_stop_index = {}
			old_chr = None
			for i in xrange(no_of_probes):
				chr = chr_pos_ls[i][0]
				if chr not in chr2start_stop_index:
					chr2start_stop_index[chr] = [i]
				if old_chr is not None and chr!=old_chr:
					chr2start_stop_index[old_chr].append(i-1)
				
				old_chr = chr
			chr2start_stop_index[old_chr].append(i)	#append the the stop index to the last chromosome
			
			output_dir = os.path.split(output_fname)[0]
			if not os.path.isdir(output_dir):
				os.makedirs(output_dir)
			writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
			from CNVNormalize import CNVNormalize
			CNVNormalize.writeArrayXProbeHeader(writer, probe_id_ls, chr_pos_ls)
			
			import numpy
			counter = 0
			for row in reader:
				array_id = int(row[0])
				ecotype_id = row[1]
				intensity_ls = row[2:]
				intensity_array = numpy.array(map(float, intensity_ls), dtype=numpy.float32)
				new_row = [array_id, ecotype_id]
				for i in range(no_of_probes):
					chr = chr_pos_ls[i][0]
					start_index, stop_index = chr2start_stop_index[chr]
					windowLeft = max(i-noOfFlankingProbes, start_index)
					windowRight = min(i+noOfFlankingProbes, stop_index)
					intensity = numpy.median(intensity_array[windowLeft:windowRight+1])
					new_row.append(intensity)
				writer.writerow(new_row)
				counter +=1
				sys.stderr.write("%s%s"%("\x08"*80, counter))
			del writer
			sys.stderr.write("Done.\n")
		
		"""
	#2010-5-30
	input_fname = os.path.expanduser('~/script/variation/data/CNV/call_48_WABQN_b200_j100_lts.tsv')
	noOfFlankingProbes = 2
	output_fname = os.path.expanduser('~/script/variation/data/CNV/call_48_WABQN_b200_j100_lts_%sFlanking.tsv'%noOfFlankingProbes)
	CNV.Normalize.medianSmooth(input_fname, output_fname, noOfFlankingProbes=noOfFlankingProbes)
	sys.exit(0)
	
	input_fname = os.path.expanduser('~/script/variation/data/CNV/call_48_WABQN_b200_j100_lts.tsv')
	noOfFlankingProbes = 4
	output_fname = os.path.expanduser('~/script/variation/data/CNV/call_48_WABQN_b200_j100_lts_%sFlanking.tsv'%noOfFlankingProbes)
	CNV.Normalize.medianSmooth(input_fname, output_fname, noOfFlankingProbes=noOfFlankingProbes)
	sys.exit(0)
	
	input_fname = os.path.expanduser('~/script/variation/data/CNV/call_53_WABQN_b100_j50_lts_98_arrays.tsv')
	noOfFlankingProbes = 2
	from pymodule.utils import addExtraToFilenamePrefix, addExtraLsToFilenamePrefix
	output_fname = addExtraToFilenamePrefix(input_fname, '%sFlanking'%noOfFlankingProbes)
	CNV.Normalize.medianSmooth(input_fname, output_fname, noOfFlankingProbes=noOfFlankingProbes)
	
	#input_fname = os.path.expanduser('~/script/variation/data/CNV/call_48_Col-Ler-WABQN_b200_j100_lts.tsv')
	noOfFlankingProbes = 4
	from pymodule.utils import addExtraToFilenamePrefix, addExtraLsToFilenamePrefix
	output_fname = addExtraToFilenamePrefix(input_fname, '%sFlanking'%noOfFlankingProbes)
	CNV.Normalize.medianSmooth(input_fname, output_fname, noOfFlankingProbes=noOfFlankingProbes)
	sys.exit(0)
	
	#2010-6-7
	input_fname = os.path.expanduser('~/script/variation/data/CNV/call_53_WABQN_b100_j50_lts.tsv')
	noOfFlankingProbes = 2
	from pymodule.utils import addExtraToFilenamePrefix, addExtraLsToFilenamePrefix
	output_fname = addExtraToFilenamePrefix(input_fname, '%sFlanking'%noOfFlankingProbes)
	CNV.Normalize.medianSmooth(input_fname, output_fname, noOfFlankingProbes=noOfFlankingProbes)
	sys.exit(0)
	
		
		#2010-6-30
		input_fname = os.path.expanduser('~/panfs/250k/CNV/call_53_WABQN_b100_j50_lts_TAIR9.tsv')
		input_fname = os.path.expanduser('~/script/variation/data/CNV/call_53_WABQN_b100_j50_lts.tsv')
		input_fname = os.path.expanduser('~/panfs/250k/CNV/call_53_WABQN_b100_j50_lts.tsv')
		noOfFlankingProbes = 4
		from pymodule.utils import addExtraToFilenamePrefix, addExtraLsToFilenamePrefix
		output_fname = addExtraToFilenamePrefix(input_fname, '%sFlanking'%noOfFlankingProbes)
		CNV.Normalize.medianSmooth(input_fname, output_fname, noOfFlankingProbes=noOfFlankingProbes)
		sys.exit(0)
		
		"""
	@classmethod
	def drawHexbin(cls, x_ls, y_ls, C_ls, fig_fname=None, gridsize=100, title=None, xlabel=None, ylabel=None,\
				colorBarLabel=None, reduce_C_function=numpy.median):
		"""
		2010-7-1
			add argument reduce_C_function()
		2010-6-28
			add argument xlabel & ylabel
		2010-5-11
		"""
		import pylab, numpy
		import matplotlib.cm as cm
		pylab.clf()
		pylab.hexbin(x_ls, y_ls, C=C_ls, gridsize=gridsize, reduce_C_function=reduce_C_function, cmap=cm.jet)
		pylab.axis([min(x_ls), max(x_ls), min(y_ls), max(y_ls)])
		if title is None:
			title = "gridsize %s, %s probes."%(gridsize, len(x_ls))
		pylab.title(title)
		if xlabel:
			pylab.xlabel(xlabel)
		if ylabel:
			pylab.ylabel(ylabel)
		cb = pylab.colorbar()
		if colorBarLabel:
			cb.set_label(colorBarLabel)
		if fig_fname:
			pylab.savefig(fig_fname, dpi=300)
	
	@classmethod
	def checkWhetherProbesAreRandomized(cls, db_250k, output_fname_prefix, \
							x_range=None, y_range=None, probeType=2):
		"""
		2010-5-16
			check whether probes are randomized on the 250k-tiling array
				draw the histogram of distance on the affy chip of two chromosomally-adjacent probes
			argument probeType, passed to DB_250k2Array.get_probes()'s run_type:
				1: SNP probes
				2/4: CNV probes
				5: QC probes
				else: all probes (SNP + CNV + QC)
		"""
		
		import Stock_250kDB
		from DB_250k2Array import DB_250k2Array
		probes, xy_ls, chr_pos_ls, probe_id_ls = DB_250k2Array.get_probes(db_250k.metadata.bind, \
									Stock_250kDB.Probes.table.name, snps=None, run_type=probeType, \
									x_range=x_range,\
									y_range=y_range)
		sys.stderr.write("Calculating adjacent-probe distance on the affy-chip ...")
		affy_chip_dist_ls = []
		no_of_probes = len(probes)
		import numpy
		for i in range(1, no_of_probes):
			prev_chr, prev_pos = chr_pos_ls[i-1][:2]
			chr, pos = chr_pos_ls[i][:2]
			if prev_chr==chr:
				prev_xy = numpy.array(xy_ls[i-1][:2])
				xy = numpy.array(xy_ls[i][:2])
				affy_chip_dist = numpy.linalg.norm(prev_xy-xy)
				affy_chip_dist_ls.append(affy_chip_dist)
		sys.stderr.write("%s distances. Done.\n"%(len(affy_chip_dist_ls)))
		
		from pymodule.utils import addExtraLsToFilenamePrefix, addExtraToFilenamePrefix
		common_fig_fname = '%s_probeType_%s.png'%(output_fname_prefix, probeType)
		if x_range:
			common_fig_fname = addExtraLsToFilenamePrefix(common_fig_fname, ['x_range']+x_range)
		if y_range:
			common_fig_fname = addExtraLsToFilenamePrefix(common_fig_fname, ['y_range']+y_range)
		
		sys.stderr.write("Drawing ...")
		import pylab
		pylab.clf()
		pylab.title("affy-chip distance histogram of adjacent probes. probeType %s"%probeType)
		pylab.hist(affy_chip_dist_ls, 50)
		pylab.savefig(common_fig_fname, dpi=300)
		sys.stderr.write("Done.\n")
	
	"""
	#2010-5-16
	y_range = None	#[1200, 1611]	#check small range to see if it works
	x_range = None	#[0, 400]
	output_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/figures/AffyChipDistHist')
	CNV.checkWhetherProbesAreRandomized(db_250k, output_fname_prefix, \
							x_range=x_range, y_range=y_range, probeType=2)
	sys.exit(0)
	"""
	
	@classmethod
	def drawArrayIntensity2DAndHist(cls, db_250k, curs=None, array_id_ls=None, output_dir=None, gridsize_ls=[50,100,200,400], \
							x_range=None, y_range=None, probeType=None, drawHist=True):
		"""
		2010-5-9
			argument probeType, passed to DB_250k2Array.get_probes()'s run_type:
				1: SNP probes
				2/4: CNV probes
				5: QC probes
				else: all probes (SNP + CNV + QC)
		2010-5-4
			array_id_ls is a list of ids in string
			2010-5-9 x_range and y_range are used to restrict data points on the 2D array
		"""
		
		import Stock_250kDB
		from DB_250k2Array import DB_250k2Array
		probes, xy_ls, chr_pos_ls, total_probe_id_ls = DB_250k2Array.get_probes(db_250k.metadata.bind, \
									Stock_250kDB.Probes.table.name, snps=None, run_type=probeType, \
									x_range=x_range,\
									y_range=y_range)
		# run_type=0 to get all probes
		
		if output_dir and not os.path.isdir(output_dir):	#2010-5-5 test if output_dir is something
			os.makedirs(output_dir)
		import numpy
		from pymodule.utils import addExtraLsToFilenamePrefix, addExtraToFilenamePrefix
		sys.stderr.write("Start to draw 2D array intensity image for probeType %s.\n"%(probeType))
		for array_id in array_id_ls:
			arrayIntensityData = DB_250k2Array.outputArray(db_250k.session, db_250k.metadata.bind, output_dir=None, \
							array_info_table=Stock_250kDB.ArrayInfo.table.name, snps=None, probes=probes,\
							array_id_ls=[array_id], xy_ls=xy_ls, \
							chr_pos_ls=chr_pos_ls, probes_id_ls=total_probe_id_ls,\
							call_method_id=None, run_type=2, \
							array_file_directory=None, outputCNVIntensity=False,\
							returnArrayIntensityData=True)
			no_of_rows, no_of_cols = arrayIntensityData.data_matrix.shape
			for j in xrange(no_of_cols):
				array_id = arrayIntensityData.col_id_ls[j]
				sys.stderr.write("\t %s/%s Array %s"%(j+1, no_of_cols, array_id))
				x_ls = []
				y_ls = []
				C_ls = []
				for i in xrange(no_of_rows):
					xy = arrayIntensityData.row_id_ls[i]
					if x_range and (xy[0]<x_range[0] or xy[0]>x_range[1]):
						continue
					if y_range and (xy[1]<y_range[0] or xy[1]>y_range[1]):
						continue
					x_ls.append(xy[0])
					y_ls.append(xy[1])
					C_ls.append(arrayIntensityData.data_matrix[i][j])
				
				common_fig_fname = 'array_%s_probeType_%s.png'%(array_id, probeType)
				if x_range:
					common_fig_fname = addExtraLsToFilenamePrefix(common_fig_fname, ['x_range']+x_range)
				if y_range:
					common_fig_fname = addExtraLsToFilenamePrefix(common_fig_fname, ['y_range']+y_range)
				hist_fig_fname = addExtraToFilenamePrefix(common_fig_fname, 'hist')
				hist_fig_fname = os.path.join(output_dir, hist_fig_fname)
				if drawHist:
					from PlotQCProbeIntensityHistogram import PlotQCProbeIntensityHistogram
					PlotQCProbeIntensityHistogram.plotHistogram(C_ls, array_id, hist_fig_fname, \
								no_of_bins=100, \
								max_intensity=None, \
								xlim=[1,5])
				for gridsize in gridsize_ls:
					sys.stderr.write(" gridSize: %s "%gridsize)
					fig_fname = addExtraToFilenamePrefix(common_fig_fname, '2D_gridsize_%s'%(gridsize))
					fig_fname = os.path.join(output_dir, fig_fname)
					title = "intensity of array %s, gridsize %s, %s probes."%(array_id, gridsize, len(x_ls))
					cls.drawHexbin(x_ls, y_ls, C_ls, fig_fname=fig_fname, gridsize=gridsize, title=title)
					sys.stderr.write("\n")
	
	"""
	array_id_ls = ["3", "4"]
	output_dir = '/tmp/'
	CNV.drawArrayIntensity2DAndHist(db_250k, curs, array_id_ls, output_dir)
	
	CNV.drawArrayIntensity2DAndHist(db_250k, curs, array_id_ls, output_dir, gridsize_ls=[100])
	
	# 2010-5-9 check where SNP probes are
	array_id_ls = ['1', '3']
	output_dir = os.path.expanduser('~/script/variation/data/CNV/normalization/2010-05-07-quantile-100/figures/')
	CNV.drawArrayIntensity2DAndHist(db_250k, curs, array_id_ls, output_dir, gridsize_ls=[200,400,1612], probeType=1)
	
	# 2010-5-9 check where CNV probes are
	array_id_ls = ['1', '3']
	output_dir = os.path.expanduser('~/script/variation/data/CNV/normalization/2010-05-07-quantile-100/figures/')
	CNV.drawArrayIntensity2DAndHist(db_250k, curs, array_id_ls, output_dir, gridsize_ls=[200,400,1612], probeType=2)
	
	# 2010-5-9 check how CNV probes are embeddd in the SNP probe panel
	array_id_ls = ['1', '3']
	output_dir = os.path.expanduser('~/script/variation/data/CNV/normalization/2010-05-07-quantile-100/figures/')
	block_size = 200
	no_of_blocks = int(math.ceil(1612/float(block_size)))
	for i in range(no_of_blocks):
		for j in range(no_of_blocks):
			x_range = [block_size*i, min(1612, block_size*(i+1))]
			y_range = [block_size*j, min(1612, block_size*(j+1))]
			CNV.drawArrayIntensity2DAndHist(db_250k, curs, array_id_ls, output_dir, x_range=x_range, y_range=y_range, \
				gridsize_ls=[block_size], probeType=2)
	sys.exit(0)
	
	# 2010-5-9 check where QC probes are
	array_id_ls = ['1', '2', '43', '139', '141', '142', '145', '1339']
	output_dir = os.path.expanduser('~/script/variation/data/CNV/normalization/2010-05-07-quantile-100/figures/')
	CNV.drawArrayIntensity2DAndHist(db_250k, curs, array_id_ls, output_dir, gridsize_ls=[1612], probeType=5)
	
	
	array_id_ls = ['1', '3', '4']
	output_dir = os.path.expanduser('~/script/variation/data/CNV/normalization/2010-05-07-quantile-100/figures/')
	gridsize_ls=[200,]
	# 2010-5-9 check where SNP probes are
	CNV.drawArrayIntensity2DAndHist(db_250k, curs, array_id_ls, output_dir, gridsize_ls=gridsize_ls, probeType=1)
	sys.exit(0)
	
	# 2010-5-9 check where QC probes are
	CNV.drawArrayIntensity2DAndHist(db_250k, curs, array_id_ls, output_dir, gridsize_ls=gridsize_ls, probeType=5)
	
	# 2010-5-9 check where CNV probes are
	CNV.drawArrayIntensity2DAndHist(db_250k, curs, array_id_ls, output_dir, gridsize_ls=gridsize_ls, probeType=2)
	
	# 2010-5-9 check how CNV probes are embeddd in the SNP-probe panel
	block_size = 200
	no_of_blocks = int(math.ceil(1612/float(block_size)))
	for i in range(no_of_blocks):
		for j in range(no_of_blocks):
			x_range = [block_size*i, min(1612, block_size*(i+1))]
			y_range = [block_size*j, min(1612, block_size*(j+1))]
			CNV.drawArrayIntensity2DAndHist(db_250k, curs, array_id_ls, output_dir, x_range=x_range, y_range=y_range, \
				gridsize_ls=[block_size], probeType=2)
	sys.exit(0)
		
		#2010-8-6 a few arrays whose number of deletions per chromosome varies greatly among chromosomes
		array_id_ls = ["182", "442", "559", "560", "562", "567", "692", "907", "988", "1019"]
		output_dir = os.path.expanduser('~/script/variation/data/CNV/ArrayIntesntiy2DPlot/')
		CNV.drawArrayIntensity2DAndHist(db_250k, array_id_ls=array_id_ls, output_dir=output_dir, gridsize_ls=[100])
		sys.exit(0)
		
	
	"""
	
	@classmethod
	def getIntensityDataOutOfAffyV3(cls, input_fname, probes=None, arrayWidth=1612):
		"""
		2010-5-8
			Lei Li's normalization program outputs data in AffyV3
		"""
		input_fname_basename = os.path.basename(input_fname)
		sys.stderr.write("Getting probe intensity from %s ..."%(input_fname_basename))
		inf = open(input_fname)
		import numpy
		if probes:
			x_ls = [0]*len(probes.probe_ls)
			y_ls = [0]*len(probes.probe_ls)
			C_ls = [numpy.nan]*len(probes.probe_ls)
		else:
			x_ls = []
			y_ls = []
			C_ls = []
		intensityRegion = False
		counter = 0 
		while 1:
			try:
				line = inf.next()
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
				break
			if line.find('[INTENSITY]')==0:
				intensityRegion = True
				# skip two lines
				inf.next()
				inf.next()
			elif line.find('[')==0 and intensityRegion is True:
				intensityRegion=False
				break
			
			if intensityRegion:
				if not line:
					continue
				row = line.split()
				if len(row)>=3:
					xy = row[:2]
					xy = map(int, xy)
					xy = (arrayWidth-1-xy[1], xy[0])	# 2010-5-11 xypos recorded in database has been turned 90 degree clockwise
					intensity = math.log10(float(row[2]))
					if probes:
						index = probes.xy2index.get(xy)
						if index is not None:
							x_ls[index] = xy[0]
							y_ls[index] = xy[1]
							C_ls[index] = intensity
							counter += 1
					else:
						x_ls.append(xy[0])
						y_ls.append(xy[1])
						C_ls.append(intensity)
						counter += 1
		sys.stderr.write(" %s probes.\n"%(counter))
		from pymodule import PassingData
		return PassingData(x_ls=x_ls, y_ls=y_ls, C_ls=C_ls)
	
	@classmethod
	def drawArrayIntensity2DFromLeiLiNormalizationResult(cls, input_fname_ls, output_dir, \
														ref_fname_ls=[], gridsize_ls=[50, 100, 400], \
														drawHist=False, probeType=2, array_id=None, probes=None,\
														x_range=None, y_range=None, ref_id2intensity_ls={}):
		"""
		2010-5-7
			ref_fname_ls is the corresponding reference array used in normalization for each in input_fname_ls
			
			draw 1. intensity histogram
				2. 2D array plot by hexbin
		"""
		
		import numpy
		import math
		x_ls = None
		y_ls = None
		C_array = []
		for i in range(len(input_fname_ls)):
			input_fname = input_fname_ls[i]
			input_fname_basename = os.path.basename(input_fname)
			intensityData = cls.getIntensityDataOutOfAffyV3(input_fname, probes=probes)
			x_ls = intensityData.x_ls
			y_ls = intensityData.y_ls
			C_ls = numpy.array(intensityData.C_ls)
			if len(ref_fname_ls)>i:
				ref_fname = ref_fname_ls[i]
				if ref_fname not in ref_id2intensity_ls:
					refIntensityData = cls.getIntensityDataOutOfAffyV3(ref_fname, probes=probes)
					ref_id2intensity_ls[ref_fname] = refIntensityData.C_ls
				refIntensityLs = ref_id2intensity_ls.get(ref_fname)
				C_ls = C_ls - numpy.array(refIntensityLs)	#subtract the reference
			C_array.append(C_ls)
		C_array = numpy.vstack(C_array)
		C_ls = numpy.median(C_array, axis=0)	#take the median of all
		if x_ls and y_ls:
			from pymodule.utils import addExtraLsToFilenamePrefix, addExtraToFilenamePrefix
			common_fig_fname = 'array_%s_probeType_%s.png'%(array_id, probeType)
			if x_range:
				common_fig_fname = addExtraLsToFilenamePrefix(common_fig_fname, ['x_range']+x_range)
			if y_range:
				common_fig_fname = addExtraLsToFilenamePrefix(common_fig_fname, ['y_range']+y_range)
			hist_fig_fname = addExtraToFilenamePrefix(common_fig_fname, 'hist')
			hist_fig_fname = os.path.join(output_dir, hist_fig_fname)
			if drawHist:
				from PlotQCProbeIntensityHistogram import PlotQCProbeIntensityHistogram
				PlotQCProbeIntensityHistogram.plotHistogram(C_ls, array_id, hist_fig_fname, \
							no_of_bins=50, \
							max_intensity=None, \
							xlim=None)
			for gridsize in gridsize_ls:
				sys.stderr.write(" gridSize: %s "%gridsize)
				fig_fname = addExtraToFilenamePrefix(common_fig_fname, '2D_gridsize_%s'%(gridsize))
				fig_fname = os.path.join(output_dir, fig_fname)
				title = "intensity of array %s, gridsize %s, %s probes."%(array_id, gridsize, len(x_ls))
				cls.drawHexbin(x_ls, y_ls, C_ls, fig_fname=fig_fname, gridsize=gridsize, title=title)
				sys.stderr.write("\n")
		return C_ls
	"""
	# 2010-5-7
	input_prefix = os.path.expanduser('~/script/variation/data/CNV/normalization/2010-05-07-quantile-100/')
	input_fname_ls = []
	for ref_id in [1, 2 ,43, 145, 139]:
		input_fname = os.path.join(input_prefix, 't_3_raw_data_By_Ref_r_%s_raw_data.cel'%(ref_id))
	output_dir = '/tmp/'
	CNV.drawArrayIntensity2DFromLeiLiNormalizationResult(input_fname_ls, output_dir)
	
	#2010-5-8 subtract the reference
	input_prefix = os.path.expanduser('~/script/variation/data/CNV/normalization/2010-05-07-quantile-100/')
	input_fname_ls = []
	ref_fname_ls = []
	output_dir = os.path.expanduser('~/script/variation/data/CNV/normalization/2010-05-07-quantile-100/')
	for targetArrayID in [4, 41, 150, 151]:
		for ref_id in [1, 2 ,43, 139, 145]:
			input_fname = os.path.join(input_prefix, 't_%s_raw_data_By_Ref_r_%s_raw_data.cel'%(targetArrayID, ref_id))
			input_fname_ls.append(input_fname)
			
			ref_fname = os.path.join(input_prefix, 'r_%s_raw_data_By_Ref_r_%s_raw_data.cel'%(ref_id, ref_id))
			ref_fname_ls.append(ref_fname)
		CNV.drawArrayIntensity2DFromLeiLiNormalizationResult(input_fname_ls, output_dir, ref_fname_ls)
	sys.exit(0)
	"""
	
	@classmethod
	def plotCNVCallGivenCNV(cls, CNV_ins, array_id2cnv_call_index_ls, extra_cnv_call_index_ls, output_dir=None, param_obj=None):
		"""
		2010-8-6
			called by CNV.drawMergedAcrossArraysCNVSpanningMultipleCNVCallsFromOneArray()
		"""
		sys.stderr.write("Plotting CNV %s (%s, %s-%s) ..."%(CNV_ins.id, CNV_ins.chromosome, CNV_ins.start, CNV_ins.stop))
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		if not hasattr(param_obj, 'baseYPos'):	#2010-8-10
			setattr(param_obj, 'baseYPos', 0)
		if not hasattr(param_obj, 'clearPlotEveryTime'):	#2010-8-10
			setattr(param_obj, 'clearPlotEveryTime', False)
		if not hasattr(param_obj, 'yticks'):	#2010-8-12
			setattr(param_obj, 'yticks', [])
		
		output_fname = os.path.join(output_dir, 'CNV_method%s_%s_%s_%s_%s.png'%(CNV_ins.cnv_method_id, CNV_ins.id, \
															CNV_ins.chromosome, CNV_ins.start, CNV_ins.stop))
		import pylab
		
		if param_obj.clearPlotEveryTime:	#2010-8-10
			pylab.clf()
		
		no_of_arrays = len(array_id2cnv_call_index_ls)
		title = 'CNV %s (%s, %s-%s), %s CNVCalls from %s arrays. method %s, type %s'%(CNV_ins.id, CNV_ins.chromosome, \
							CNV_ins.start, CNV_ins.stop, len(CNV_ins.cnv_call_ls), no_of_arrays,\
							CNV_ins.cnv_method_id, CNV_ins.cnv_type_id)
		pylab.title(title)
		ytick = 'CNV %s'%(CNV_ins.id)
		yticks = [ytick]
		param_obj.yticks.append(ytick)
		
		param_obj.baseYPos += 1
		pylab.plot([CNV_ins.start, CNV_ins.stop], [param_obj.baseYPos, param_obj.baseYPos], linewidth=5, alpha=0.3)	#extra thick
		
		#array_id_ls = array_id2cnv_call_index_ls.keys()
		#array_id_ls.sort()	#sorted in ecotype_id order
		extra_cnv_call_index_ls.sort()
		for i in xrange(len(extra_cnv_call_index_ls)):
			start, stop, ecotype_id, array_id, cnv_call_index = extra_cnv_call_index_ls[i][:5]
			#ecotype_id, array_id = array_id_ls[i][:2]
			
			#cnv_call_index_ls = array_id2cnv_call_index_ls.get(array_id_ls[i])
			#for cnv_call_index in cnv_call_index_ls:
			cnv_call = CNV_ins.cnv_call_ls[cnv_call_index]
			ytick = 'aid %s eid %s %s %s'%(array_id, ecotype_id, param_obj.ecotype_id2nativename.get(ecotype_id), cnv_call.id)
			param_obj.yticks.append(ytick)
			yticks.append(ytick)
			param_obj.baseYPos += 1
			pylab.plot([cnv_call.start, cnv_call.stop], [param_obj.baseYPos, param_obj.baseYPos])
		
		no_of_yticks = len(yticks)
		pylab.ylim([-10, param_obj.baseYPos+10])
		if param_obj.clearPlotEveryTime:	#2010-8-10 save it every time if this is true.
			yticks_index_ls = range(1, param_obj.baseYPos+1)
			pylab.yticks(yticks_index_ls, yticks, size=2)
			if hasattr(param_obj,'start') and hasattr(param_obj, 'stop'):
				pylab.axvspan(param_obj.start, param_obj.stop, facecolor='g', alpha=0.3)
			pylab.savefig(output_fname, dpi=400)
		sys.stderr.write("Done.\n")
	
	@classmethod
	def drawMergedAcrossArraysCNVSpanningMultipleCNVCallsFromOneArray(cls, db_250k, output_dir=None,\
											cnv_method_id=None, cnv_type_id=None):
		"""
		2010-8-6
			draw original CNVCall from which the merged CNV is derived from.
		"""
		from common import get_ecotypeid2nativename
		ecotype_id2nativename = get_ecotypeid2nativename(db_250k.metadata.bind)
		
		from pymodule import PassingData
		param_obj = PassingData(no_of_valid_deletions=0, no_of_deletions=0, \
				array_id2label={}, ecotype_id2nativename=ecotype_id2nativename,\
				clearPlotEveryTime=True, baseYPos=0)
		
		sys.stderr.write("Drawing Merged-across-array CNVs that span multiple CNV Calls from one array. method %s type %s ...\n"%\
						(cnv_method_id, cnv_type_id))
		
		
		import Stock_250kDB
		
		TableClass=Stock_250kDB.CNV
		query = TableClass.query.filter_by(cnv_method_id=cnv_method_id).filter_by(cnv_type_id=cnv_type_id)
		counter = 0 
		for row in query:
			counter += 1
			array_id2cnv_call_index_ls = {}	#the key is (ecotype_id, array_id) in order to sort in terms of ecotype_id
			extra_cnv_call_index_ls = []
			no_of_cnv_calls = 0
			for i in xrange(len(row.cnv_call_ls)):
				cnv_call = row.cnv_call_ls[i]
				extra_cnv_call_index_ls.append((cnv_call.start, cnv_call.stop, cnv_call.array.maternal_ecotype_id, cnv_call.array_id, i))
				key = (cnv_call.array.maternal_ecotype_id, cnv_call.array_id,)
				if key not in array_id2cnv_call_index_ls:
					array_id2cnv_call_index_ls[key] = []
				array_id2cnv_call_index_ls[key].append(i)
				no_of_cnv_calls += 1
			if len(array_id2cnv_call_index_ls)!=no_of_cnv_calls:
				param_obj.baseYPos = 0	# 2010-8-10 param_obj.clearPlotEveryTime=True. so reset the y-position to zero.
				cls.plotCNVCallGivenCNV(row, array_id2cnv_call_index_ls, extra_cnv_call_index_ls, output_dir=output_dir, param_obj=param_obj)
			if counter%5000==0:
				sys.stderr.write("%s%s"%('\x08'*80, counter))
	"""
		#2010-8-6
		output_dir = os.path.expanduser('~/script/variation/data/CNV/MergedAcrossArraysCNVSpanningMultipleCNVCallsFromOneArray/')
		CNV.drawMergedAcrossArraysCNVSpanningMultipleCNVCallsFromOneArray(db_250k, output_dir=output_dir,\
									cnv_method_id=22, cnv_type_id=1)
		sys.exit(0)
	"""
	
	@classmethod
	def drawMergedAcrossArraysCNVGivenPosition(cls, db_250k, output_dir=None,\
											cnv_method_id=None, cnv_type_id=None, chromosome=None, start=None, stop=None):
		"""
		2010-8-10
			draw original CNVCall from which the merged CNV is derived from.
		"""
		from common import get_ecotypeid2nativename
		ecotype_id2nativename = get_ecotypeid2nativename(db_250k.metadata.bind)
		
		from pymodule import PassingData
		param_obj = PassingData(no_of_valid_deletions=0, no_of_deletions=0, \
				array_id2label={}, ecotype_id2nativename=ecotype_id2nativename, \
				clearPlotEveryTime=False, baseYPos=0, start=start, stop=stop, chromosome=chromosome, yticks=[])
		
		sys.stderr.write("Drawing Merged-across-array CNVs that overlap with Chr %s (%s-%s). method %s type %s ...\n"%\
						(cnv_method_id, cnv_type_id, chromosome, start, stop))
		
		
		import Stock_250kDB
		
		TableClass=Stock_250kDB.CNV
		query = TableClass.query.filter_by(cnv_method_id=cnv_method_id).filter_by(cnv_type_id=cnv_type_id)
		
		query = db_250k.limitQueryByChrPosition(query, TableClass, chr=chromosome, start=start, stop=stop)
		
		counter = 0
		no_of_total_cnv_calls = 0
		array_id_set = set()
		for row in query:
			counter += 1
			array_id2cnv_call_index_ls = {}	#the key is (ecotype_id, array_id) in order to sort in terms of ecotype_id
			extra_cnv_call_index_ls = []
			no_of_cnv_calls = 0
			for i in xrange(len(row.cnv_call_ls)):
				cnv_call = row.cnv_call_ls[i]
				extra_cnv_call_index_ls.append((cnv_call.start, cnv_call.stop, cnv_call.array.maternal_ecotype_id, cnv_call.array_id, i))
				key = (cnv_call.array.maternal_ecotype_id, cnv_call.array_id,)
				if key not in array_id2cnv_call_index_ls:
					array_id2cnv_call_index_ls[key] = []
				array_id2cnv_call_index_ls[key].append(i)
				array_id_set.add(cnv_call.array_id)
				no_of_cnv_calls += 1
			no_of_total_cnv_calls += no_of_cnv_calls
			cls.plotCNVCallGivenCNV(row, array_id2cnv_call_index_ls, extra_cnv_call_index_ls, \
								output_dir=output_dir, param_obj=param_obj)
			if counter%5000==0:
				sys.stderr.write("%s%s"%('\x08'*80, counter))
		if not param_obj.clearPlotEveryTime:
			import pylab
			output_fname = os.path.join(output_dir, 'CNV_method%s_%s_%s_%s.png'%(cnv_method_id, \
															chromosome, start, stop))
			title = '%s CNVs covering %s CNVCalls in %s arrays. method %s type %s'%(counter, \
							no_of_total_cnv_calls, len(array_id_set), cnv_method_id, cnv_type_id)
			pylab.title(title)
			yticks_index_ls = range(1, param_obj.baseYPos+1)
			pylab.yticks(yticks_index_ls, param_obj.yticks, size=2)
			pylab.axvspan(start, stop, facecolor='g', alpha=0.3)
			pylab.savefig(output_fname, dpi=400)
	"""
		#2010-8-10
		output_dir = os.path.expanduser('~/script/variation/data/CNV/MergedAcrossArraysCNV/')
		CNV.drawMergedAcrossArraysCNVGivenPosition(db_250k, output_dir=output_dir,\
									cnv_method_id=23, cnv_type_id=1, chromosome=1, start=6159423, stop=6160879)
		sys.exit(0)
		
		#2010-8-10
		for cnv_method_id in [18, 19, 22]:
			output_dir = os.path.expanduser('~/script/variation/data/CNV/MergedAcrossArraysCNV/')
			CNV.drawMergedAcrossArraysCNVGivenPosition(db_250k, output_dir=output_dir,\
									cnv_method_id=cnv_method_id, cnv_type_id=1, chromosome=1, start=6156000, stop=6162000)
		sys.exit(0)
		
		#2010-8-10
		output_dir = os.path.expanduser('~/script/variation/data/CNV/MergedAcrossArraysCNV/')
		CNV.drawMergedAcrossArraysCNVGivenPosition(db_250k, output_dir=output_dir,\
									cnv_method_id=24, cnv_type_id=1, chromosome=2, start=7063000, stop=7065260)
		sys.exit(0)
	"""
	
	@classmethod
	def drawCNVCallWithinRegion(cls, db_250k, output_dir=None,\
							cnv_method_id=None, cnv_type_id=None, chromosome=None, start=None, stop=None):
		"""
		2010-8-10
			draw original CNVCall which falls into chromosome, (start-stop)
		"""
		from common import get_ecotypeid2nativename
		ecotype_id2nativename = get_ecotypeid2nativename(db_250k.metadata.bind)
		
		from pymodule import PassingData
		param_obj = PassingData(no_of_valid_deletions=0, no_of_deletions=0, \
				array_id2label={}, ecotype_id2nativename=ecotype_id2nativename, \
				clearPlotEveryTime=False, baseYPos=0, start=start, stop=stop, chromosome=chromosome)
		
		sys.stderr.write("Drawing CNVCalls that overlap with Chr %s (%s-%s). method %s type %s ...\n"%\
						(cnv_method_id, cnv_type_id, chromosome, start, stop))
		
		
		import Stock_250kDB
		
		TableClass=Stock_250kDB.CNVCall
		query = TableClass.query.filter_by(cnv_method_id=cnv_method_id).filter_by(cnv_type_id=cnv_type_id)
		
		query = db_250k.limitQueryByChrPosition(query, TableClass, chr=chromosome, start=start, stop=stop)
		
		query = query.order_by(TableClass.chromosome).order_by(TableClass.start).order_by(TableClass.stop)
		counter = 0
		no_of_total_cnv_calls = 0
		array_id_set = set()
		ecotype_id_set = set()
		yticks = []
		import pylab
		pylab.clf()
		for row in query:
			counter += 1
			no_of_total_cnv_calls += 1
			array_id_set.add(row.array_id)
			ecotype_id = row.array.maternal_ecotype_id
			ecotype_id_set.add(ecotype_id)
			ytick = 'aid %s eid %s %s'%(row.array_id, ecotype_id, param_obj.ecotype_id2nativename.get(ecotype_id))
			yticks.append('%s %s'%(ytick, row.id))
			pylab.plot([row.start, row.stop], [counter, counter])
			if counter%10==0:
				sys.stderr.write("%s%s"%('\x08'*80, counter))
		sys.stderr.write("%s%s\n"%('\x08'*80, counter))
		no_of_yticks = len(yticks)
		yticks_index_ls = range(1, counter+1)
		pylab.yticks(yticks_index_ls, yticks, size=2)
		pylab.ylim([-10, no_of_yticks+10])

		output_fname = os.path.join(output_dir, 'CNV_method%s_%s_%s_%s.png'%(cnv_method_id, \
														chromosome, start, stop))
		title = '%s CNVCalls in %s arrays (%s ecotypes). method %s type %s'%\
				(no_of_total_cnv_calls, len(array_id_set), len(ecotype_id_set), cnv_method_id, cnv_type_id)
		pylab.title(title)
		pylab.axvspan(start, stop, facecolor='g', alpha=0.3)
		pylab.savefig(output_fname, dpi=400)
	"""
		#2010-8-10
		output_dir = os.path.expanduser('~/script/variation/data/CNV/MergedAcrossArraysCNV/')
		CNV.drawCNVCallWithinRegion(db_250k, output_dir=output_dir,\
							cnv_method_id=23, cnv_type_id=1, chromosome=1, start=6159423, stop=6160879)
		sys.exit(0)
		
		#2010-8-12
		output_dir = os.path.expanduser('~/script/variation/data/CNV/MergedAcrossArraysCNV/')
		CNV.drawCNVCallWithinRegion(db_250k, output_dir=output_dir,\
									cnv_method_id=16, cnv_type_id=1, chromosome=1, start=540847, stop=562919)
		sys.exit(0)
	"""
	
	@classmethod
	def drawMultiArrayIntensity2DFromLeiLiNormalizationResult(cls, db_250k, input_prefix=None, output_dir=None, \
														gridsize_ls=[200], probeType=2, x_range=None,\
														y_range=None, targetArrayIDLs = [3, 4, 41, 150, 151], \
														ref_id_ls = [1, 2 ,43, 139, 145]):
		"""
		2010-5-7
			ref_fname_ls is the corresponding reference array used in normalization for each in input_fname_ls
		"""
		
		from DB_250k2Array import DB_250k2Array
		import Stock_250kDB, numpy, rpy
		from CNVNormalize import CNVNormalize
		probes, xy_ls, chr_pos_ls, probe_id_ls = DB_250k2Array.get_probes(db_250k.metadata.bind, \
									Stock_250kDB.Probes.table.name, snps=None, run_type=probeType, \
									x_range=x_range,\
									y_range=y_range, constructChrPos2index=False, constructXY2index=True,\
									need_xy_ls=False)
		
		if output_dir and not os.path.isdir(output_dir):	#2010-5-5 test if output_dir is something
			os.makedirs(output_dir)
		
		no_of_arrays = len(targetArrayIDLs)
		data_matrix = numpy.zeros([len(probe_id_ls), no_of_arrays], numpy.float32)
		ref_id2intensity_ls = {}
		for i in range(no_of_arrays):
			targetArrayID = targetArrayIDLs[i]
			input_fname_ls = []
			ref_fname_ls = []
			for ref_id in ref_id_ls:
				input_fname = os.path.join(input_prefix, 't_%s_raw_data_By_Ref_r_%s_raw_data.cel'%(targetArrayID, ref_id))
				input_fname_ls.append(input_fname)
				
				ref_fname = os.path.join(input_prefix, 'r_%s_raw_data_By_Ref_r_%s_raw_data.cel'%(ref_id, ref_id))
				ref_fname_ls.append(ref_fname)
			intensityOfOneArray = CNV.drawArrayIntensity2DFromLeiLiNormalizationResult(input_fname_ls, output_dir, ref_fname_ls, \
														gridsize_ls=gridsize_ls, probeType=probeType, array_id=targetArrayID, \
														probes=probes,\
														x_range=x_range, y_range=y_range, ref_id2intensity_ls=ref_id2intensity_ls)
			data_matrix[:,i] = intensityOfOneArray
		header =['probes_id'] + targetArrayIDLs + ['chromosome', 'position']
		
		
		from pymodule.utils import addExtraLsToFilenamePrefix, addExtraToFilenamePrefix
		targetArrays = "_".join(map(str, targetArrayIDLs))
		refArrays = "_".join(map(str, ref_id_ls))
		output_fname = 'subArrayQuantileNorm_%s_MedianAmong_Ref_%s_probeType_%s.tsv'%(targetArrays, refArrays, probeType)
		if x_range:
			output_fname = addExtraLsToFilenamePrefix(output_fname, ['x_range']+x_range)
		if y_range:
			output_fname = addExtraLsToFilenamePrefix(output_fname, ['y_range']+y_range)
		output_fname = os.path.join(output_dir, output_fname)
		CNVNormalize.output(data_matrix, probe_id_ls, chr_pos_ls, header, output_fname)
	
	"""
	input_prefix = os.path.expanduser('~/script/variation/data/CNV/normalization/2010-05-07-quantile-100/')
	output_dir = os.path.expanduser('~/script/variation/data/CNV/normalization/2010-05-07-quantile-100/medianOverRef1/')
	
	#2010-5-11 test
	y_range = [1200, 1611]	#check small range to see if it works
	x_range = [0, 100]
	CNV.drawMultiArrayIntensity2DFromLeiLiNormalizationResult(db_250k, input_prefix=input_prefix, output_dir=output_dir, \
									gridsize_ls=[200], probeType=2, x_range=x_range,\
									y_range=y_range)
	sys.exit(0)
	
	#2010-5-11 full array
	y_range = None	#[1200, 1611]	#check small range to see if it works
	x_range = None	#[0, 100]
	CNV.drawMultiArrayIntensity2DFromLeiLiNormalizationResult(db_250k, input_prefix=input_prefix, output_dir=output_dir, \
									gridsize_ls=[200], probeType=2, x_range=x_range,\
									y_range=y_range, ref_id_ls = [1, 2 , 139, 145])
	sys.exit(0)
	
	#2010-5-13 reference without array 43
	y_range = None	#[1200, 1611]	#check small range to see if it works
	x_range = None	#[0, 100]
	CNV.drawMultiArrayIntensity2DFromLeiLiNormalizationResult(db_250k, input_prefix=input_prefix, output_dir=output_dir, \
									gridsize_ls=[200], probeType=2, x_range=x_range,\
									y_range=y_range, ref_id_ls = [1, 2 , 139, 145])
	sys.exit(0)

	#2010-5-25
	input_prefix = '/media/FreeAgent Disk/huangyu/2010-05-13-quantile-70-40'
	output_dir = os.path.expanduser('~/script/variation/data/CNV/normalization/medianOverRef1/')
	y_range = None	#[1200, 1611]	#check small range to see if it works
	x_range = None	#[0, 100]
	CNV.drawMultiArrayIntensity2DFromLeiLiNormalizationResult(db_250k, input_prefix=input_prefix, output_dir=output_dir, \
									gridsize_ls=[], probeType=2, x_range=x_range,\
									y_range=y_range, ref_id_ls = [1, 2 , 139, 145])
	sys.exit(0)
	
	#2010-5-25
	input_prefix = '/media/FreeAgent Disk/huangyu/2010-05-10-quantile'
	output_dir = os.path.expanduser('~/script/variation/data/CNV/normalization/medianOverRef-w50/')
	y_range = None	#[1200, 1611]	#check small range to see if it works
	x_range = None	#[0, 100]
	CNV.drawMultiArrayIntensity2DFromLeiLiNormalizationResult(db_250k, input_prefix=input_prefix, output_dir=output_dir, \
									gridsize_ls=[], probeType=2, x_range=x_range,\
									y_range=y_range, ref_id_ls = [1, 2 , 43, 139, 145])
	sys.exit(0)
	
	
	
	"""
	
class AnalyzeSNPData(object):
	@classmethod
	def DrawStrainSNP_NA_PercHist(cls, data_matrix_fname, output_fname, need_savefig=0):
		"""
		2007-03-20
		"""
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(data_matrix_fname)
		import Numeric
		data_matrix = Numeric.array(data_matrix)
		strain_NA_perc_ls = []
		for i in range(data_matrix.shape[0]):
			strain_NA_perc_ls.append(sum(data_matrix[i,:]==0)/float(data_matrix.shape[1]))	#0 is NA
		SNP_NA_perc_ls = []
		for i in range(data_matrix.shape[1]):
			SNP_NA_perc_ls.append(sum(data_matrix[:,i]==0)/float(data_matrix.shape[0]))	#0 is NA
		import pylab,os
		base_fname = os.path.basename(data_matrix_fname)
		pylab.clf()
		pylab.hist(strain_NA_perc_ls, 20)
		pylab.title("%s Strain NA perc histogram"%base_fname)
		if need_savefig:
			pylab.savefig('%s_strain_NA_perc.eps'%output_fname, dpi=300)
			pylab.savefig('%s_strain_NA_perc.svg'%output_fname, dpi=300)
			pylab.savefig('%s_strain_NA_perc.png'%output_fname, dpi=300)
		pylab.show()
		
		pylab.clf()
		pylab.hist(SNP_NA_perc_ls, 20)
		pylab.title("%s SNP NA perc histogram"%base_fname)
		if need_savefig:
			pylab.savefig('%s_SNP_NA_perc.eps'%output_fname, dpi=300)
			pylab.savefig('%s_SNP_NA_perc.svg'%output_fname, dpi=300)
			pylab.savefig('%s_SNP_NA_perc.png'%output_fname, dpi=300)
		pylab.show()
	
	"""
	DrawStrainSNP_NA_PercHist('./script/variation/data/justin_data_y.csv', './script/variation/data/justin_data_y', 1)
	"""
	@classmethod
	def DrawStrain_Heterozygotes_PercHist(cls, data_matrix_fname, output_fname, need_savefig=0):
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(data_matrix_fname)
		import Numeric
		data_matrix = Numeric.array(data_matrix)
		strain_Hetero_perc_ls = []
		for i in range(data_matrix.shape[0]):
			strain_Hetero_perc_ls.append(sum(data_matrix[i,:]>4)/float(data_matrix.shape[1]))	#bigger than 4 is heterozygotes
		
		import pylab,os
		base_fname = os.path.basename(data_matrix_fname)
		pylab.clf()
		pylab.hist(strain_Hetero_perc_ls, 20)
		pylab.title("%s Strain Heterozygote perc histogram"%base_fname)
		if need_savefig:
			pylab.savefig('%s_strain_hz_perc.eps'%output_fname, dpi=300)
			pylab.savefig('%s_strain_hz_perc.svg'%output_fname, dpi=300)
			pylab.savefig('%s_strain_hz_perc.png'%output_fname, dpi=300)
		pylab.show()
		return strain_Hetero_perc_ls
	
	"""
	strain_Hetero_perc_ls = DrawStrain_Heterozygotes_PercHist('./script/variation/data/justin_data_y.csv', './script/variation/data/justin_data_y', 1)
	"""
	
	"""
	2007-09-24 increase the #bins of histogram to 40
	2007-03-21 draw histogram of pairwise accession genetic distance
	"""
	
	@classmethod
	def turn_row_id2pairwise_dist_into_ecotype_id_pair_dict(cls, row_id2pairwise_dist):
		"""
		2010-10-23
			called by outputPairwiseDistanceFromTwoDatasets()
		"""
		ecotype_id_pair2dist = {}
		for row_id, pairwise_dist in row_id2pairwise_dist.iteritems():
			row_id = int(row_id)
			for dist in pairwise_dist:
				mismatch_rate, row_id2, no_of_mismatches, no_of_non_NA_pairs = dist[:4]
				row_id2 = int(row_id2)
				key = (min(row_id, row_id2), max(row_id, row_id2))
				ecotype_id_pair2dist[key] = [mismatch_rate, no_of_mismatches, no_of_non_NA_pairs]
		return ecotype_id_pair2dist
	
	@classmethod
	def outputPairwiseDistanceFromTwoDatasets(cls, input_fname1, input_fname2, output_fname, ref_ecotype_id=6909):
		"""
		2010-10-23
			calculate the genetic distance from non-Col ecotypes to Col for each dataset and then output the two in juxtaposation 
		"""
		from pymodule import SNPData
		snpData1 = SNPData(input_fname=input_fname1, turn_into_array=1, ignore_2nd_column=1)
		snpData2 = SNPData(input_fname=input_fname2, turn_into_array=1, ignore_2nd_column=1)
		
		row_id2pairwise_dist1 = snpData1.calRowPairwiseDist(NA_set =set(['NA', 'N', -2, '|']), ref_row_id=str(ref_ecotype_id))
		ecotype_id_pair2dist1 = cls.turn_row_id2pairwise_dist_into_ecotype_id_pair_dict(row_id2pairwise_dist1)
		
		#row_id2pairwise_dist2 = snpData2.calRowPairwiseDist(NA_set =set(['NA', 'N', -2, '|']), ref_row_id=str(ref_ecotype_id))
		#ecotype_id_pair2dist2 = cls.turn_row_id2pairwise_dist_into_ecotype_id_pair_dict(row_id2pairwise_dist2)
		# CNV dataset doesn't have ref_ecotype_id 6909
		row_id2fractionData = snpData2.calFractionOfLociCarryingNonRefAllelePerRow(NA_set =set(['NA', 'N', -2, '|']), ref_allele=0)
		ecotype_id_pair2dist2 = {}
		for row_id, fractionData in row_id2fractionData.iteritems():
			row_id = int(row_id)
			key = (min(ref_ecotype_id, row_id), max(ref_ecotype_id, row_id))
			ecotype_id_pair2dist2[key] = fractionData
		
		key_overlap = set(ecotype_id_pair2dist1.keys()) & set(ecotype_id_pair2dist2.keys())
		
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter ='\t')
		header =  ['ecotype_id_pair', 'mismatch_rate_1', 'mismatch_rate_2', 'no_of_non_NA_pairs_1', 'no_of_non_NA_pairs_2']
		writer.writerow(header)
		for key in key_overlap:
			dist_1 = ecotype_id_pair2dist1.get(key)
			dist_2 = ecotype_id_pair2dist2.get(key)
			key = map(str, key)
			row = ['-'.join(key), dist_1[0], dist_2[0], dist_1[-1], dist_2[-1]]
			writer.writerow(row)
		del writer
	"""
		# 2010-10-23 test
		input_fname1 = '/Network/Data/250k/tmp-yh/250k_data/call_method_17_test.tsv'
		input_fname2 = '/Network/Data/250k/tmp-yh/250k_data/call_method_17_test.tsv'
		output_fname = '/tmp/pairwiseDistFromTwoDatasets.tsv'
		AnalyzeSNPData.outputPairwiseDistanceFromTwoDatasets(input_fname1, input_fname2, output_fname)
		sys.exit(0)
		
		
		# 2010-10-23
		input_fname1 = '/Network/Data/250k/db/dataset/call_method_54.tsv'
		input_fname2 = os.path.expanduser('~/mnt/panfs/250k/CNV/NonOverlapCNVAsSNP_cnvMethod20.tsv')
		output_fname = os.path.expanduser('~/pairwiseDist_call54_vs_cnv20.tsv')
		AnalyzeSNPData.outputPairwiseDistanceFromTwoDatasets(input_fname1, input_fname2, output_fname)
		sys.exit(0)
		
	"""
	
	
	
	@classmethod
	def DrawDistanceHistogram(cls, data_matrix_fname, output_fname, need_savefig=0):
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(data_matrix_fname)
		import Numeric
		distance_ls = []
		no_of_NA_pairs = 0
		no_of_strains = len(data_matrix)
		no_of_snps = len(data_matrix[0])
		distance_matrix = Numeric.zeros([no_of_strains, no_of_strains], Numeric.Float)
		for i in range(no_of_strains):
			for j in range(i+1, no_of_strains):
				no_of_valid_pairs = 0.0
				no_of_matching_pairs = 0.0
				for k in range(no_of_snps):
					if data_matrix[i][k]!=0 and data_matrix[j][k]!=0:
						no_of_valid_pairs += 1
						if data_matrix[i][k] == data_matrix[j][k]:
							no_of_matching_pairs += 1
				if no_of_valid_pairs!=0:
					distance = 1 - no_of_matching_pairs/no_of_valid_pairs
					distance_matrix[i,j] = distance_matrix[j,i] = distance
					distance_ls.append(distance)
				else:
					no_of_NA_pairs += 1
		print "out of %s pairs, %s are NA"%((no_of_strains*(no_of_strains-1))/2, no_of_NA_pairs)
		import pylab
		pylab.clf()
		pylab.hist(distance_ls, 40)
		pylab.title("Histogram of non-NA distances")
		if need_savefig:
			pylab.savefig('%s_distance_hist.eps'%output_fname, dpi=300)
			pylab.savefig('%s_distance_hist.svg'%output_fname, dpi=300)
			pylab.savefig('%s_distance_hist.png'%output_fname, dpi=300)
		pylab.show()
		return distance_matrix
	
	"""
	distance_matrix = DrawDistanceHistogram('./script/variation/data/justin_data_filtered.csv', './script/variation/data/justin_data_filtered' , 1)
	"""
	
	"""
	2007-09-17
		check allele frequency
	"""
	@classmethod
	def cal_maf_vector(cls, input_fname):
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
		import Numeric
		data_matrix = Numeric.array(data_matrix)
		from EstimateSelfingGeneration import EstimateSelfingGeneration
		EstimateSelfingGeneration_instance = EstimateSelfingGeneration()
		locus_allele_prob_vector = EstimateSelfingGeneration_instance.cal_locus_allele_prob_vector(data_matrix)
		maf_vector = Numeric.zeros(locus_allele_prob_vector.shape[0], Numeric.Float)	#the minor allele frequency vector
		for i in range(locus_allele_prob_vector.shape[0]):
			maf_vector[i] = min(locus_allele_prob_vector[i])
		import pylab
		pylab.hist(maf_vector, 10)
		return maf_vector
	"""
	input_fname = 'script/variation/stock20070829/data_row_na_col_na_bad_snps.tsv'
	input_fname = '/mnt/hpc-cmb/KW/input/250K_method_5_after_imputation_noRedundant_051908.tsv'
	maf_vector = cal_maf_vector(input_fname)
	import pylab
	pylab.clf()
	pylab.plot(range(len(maf_vector)), maf_vector)
	"""
	
	"""
	2007-09-23
		calculate the percentage of NAs in a data_matrix
	"""
	@classmethod
	def calculate_NA_perc(cls, input_fname):
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
		import Numeric
		data_matrix = Numeric.array(data_matrix)
		no_of_NAs = 0
		for i in range(data_matrix.shape[0]):
			for j in range(data_matrix.shape[1]):
				if data_matrix[i][j] == 0:
					no_of_NAs += 1
		no_of_totals = data_matrix.shape[0]*data_matrix.shape[1]
		print '%s/%s=%s'%(no_of_NAs, no_of_totals, float(no_of_NAs)/no_of_totals)
	
	"""
	calculate_NA_perc('.//script/variation/stock20070919/data.tsv')
	"""
	
	
	#2008-2-16
	#check the frequency of all allele-combinations between two SNPs
	@classmethod
	def get_allele_combo2freq(cls, snp_data, snp_id1, snp_id2):
		snp_index1 = snp_data.col_id2col_index[snp_id1]
		snp_index2 = snp_data.col_id2col_index[snp_id2]
		no_of_rows, no_of_cols = snp_data.data_matrix.shape
		allele_combo2freq = {}
		for i in range(no_of_rows):
			allele_combo = (snp_data.data_matrix[i][snp_index1], snp_data.data_matrix[i][snp_index2])
			if allele_combo not in allele_combo2freq:
				allele_combo2freq[allele_combo] = 0
			allele_combo2freq[allele_combo] += 1
		return allele_combo2freq
	"""
	from pymodule import SNPData
	snp_data = SNPData(input_fname='panfs/250k/call_method_29_binary.tsv', turn_into_array=1)
	get_allele_combo2freq(snp_data, '1_18234094', '5_18607474')
	"""
	
	@classmethod
	def getAlignmentMatrixFromFasta(cls, fasta_input_fname, output_fname, chromosome=1, start=1, pickPolymorphicColumns=True):
		"""
		2009-5-28
			add argument pickPolymorphicColumns
		2009-3-26
			generate alignment matrix out of a sequence file in fasta format
				(if ID includes an underscore, the part behind underscore is regarded as real accession id.)
			chromosome & start are the position of the 1st nucleotide in the alignment
		"""
		import os, sys, numpy
		sys.stderr.write("Getting alignment matrix from %s ..."%(fasta_input_fname))
		snp_pos_ls = []
		accession_id_ls = []
		name_ls = []
		data_matrix = []
		inf = open(fasta_input_fname)
		from Bio import SeqIO
		from DiscoverSNPFromAlignment import DiscoverSNPFromAlignment
		from pymodule import dict_map, nt2number, PassingData
		counter = 0
		for seq_record in SeqIO.parse(inf, "fasta"):
			if counter == 0:
				snp_pos_ls = DiscoverSNPFromAlignment.get_snp_pos_ls(seq_record.seq.tostring().upper(), chromosome, start)
			record_id_split = seq_record.id.split('_')
			if len(record_id_split)==2:
				record_id = record_id_split[1]
			else:
				record_id = record_id_split[0]
			accession_id_ls.append(record_id)
			name_ls.append(seq_record.id)
			data_row = dict_map(nt2number, seq_record.seq.tostring().upper())
			data_matrix.append(data_row)
			counter += 1
		data_matrix = numpy.array(data_matrix, numpy.int8)
		passingdata = PassingData(snp_pos_ls=snp_pos_ls, accession_id_ls=accession_id_ls, name_ls=name_ls, data_matrix=data_matrix)
		sys.stderr.write(' %s accessions, %s bases. Done.\n'%(len(accession_id_ls), len(snp_pos_ls)))
		
		if pickPolymorphicColumns:
			DiscoverSNPFromAlignment.pickPolymorphicColumns(passingdata)
		
		header = ['id', 'name']
		for snp_pos in passingdata.snp_pos_ls:
			header.append('%s_%s_%s'%snp_pos)
		passingdata.header = header
		return passingdata
	
	@classmethod
	def outputSNPsOutOfFastaAlignMatrixInYuFormat(cls, fasta_input_fname, output_fname, chromosome=1, start=1):
		"""
		2009-5-28
			call getAlignmentMatrixFromFasta(), then write the returned data into file
		"""
		passingdata = cls.getAlignmentMatrixFromFasta(fasta_input_fname, output_fname, chromosome, start)
		from pymodule import write_data_matrix
		write_data_matrix(passingdata.data_matrix, output_fname, passingdata.header, \
						passingdata.accession_id_ls, passingdata.name_ls)
		return passingdata
	
	
	"""
fasta_input_fname = os.path.expanduser('~/script/variation/doc/FLC/CarolineDeanFLC_rc_with_ref.align.fasta')
output_fname = os.path.expanduser('~/script/variation/doc/FLC/CarolineDeanFLC_rc_with_ref_snp_matrix.tsv')
snpData = AnalyzeSNPData.outputSNPsOutOfFastaAlignMatrixInYuFormat(fasta_input_fname, output_fname, chromosome=5, start=3170860)
	"""
	
	@classmethod
	def filterColumnsOfFasta(cls, fasta_input_fname, output_fname, chromosome=1, start=1, pickPolymorphicColumns=False):
		"""
		2009-5-28
			read the fasta_input_fname (in fasta format, if ID includes an underscore, the part behind underscore is regarded as real accession id.)
			purpose of this program is to remove columns full of 'N' or '-', which could cause trouble for alignment with other sequences 
			output in the same fasta format for further data fiddling
		"""
		passingdata = cls.getAlignmentMatrixFromFasta(fasta_input_fname, output_fname, chromosome, start, pickPolymorphicColumns=False)
		sys.stderr.write("Filtering all-NA or all-deletion columns ... ")
		from pymodule import nt2number, number2nt, dict_map, number2single_char_nt
		import numpy
		no_of_accessions = len(passingdata.accession_id_ls)
		all_del_col = [nt2number['-']]*no_of_accessions
		all_NA_col = [nt2number['N']]*no_of_accessions
		non_del_NA_col_indices = []
		all_NA_counter = 0
		all_del_counter = 0
		for j in range(passingdata.data_matrix.shape[1]):
			all_del_truth_vector = passingdata.data_matrix[:,j]==all_del_col
			all_NA_truth_vector = passingdata.data_matrix[:,j]==all_NA_col
			
			if all_del_truth_vector.all():
				all_del_counter += 1
			elif all_NA_truth_vector.all():
				all_NA_counter += 1
			else:
				non_del_NA_col_indices.append(j)
		
		
		outf = open(output_fname, 'w')
		for i in range(no_of_accessions):
			sequence_id = passingdata.name_ls[i]
			outf.write('>%s\n'%sequence_id)
			nt_ls = dict_map(number2single_char_nt, passingdata.data_matrix[i, non_del_NA_col_indices])
			outf.write('%s\n'%''.join(nt_ls))
		del outf
		sys.stderr.write("%s columns are all NA. %s columns are all del. Done.\n"%(all_NA_counter, all_del_counter))
	
	"""
fasta_input_fname = os.path.expanduser('~/script/variation/data/ACD6/ACD6.non_Kz10.sequenced.fasta')
output_fname = os.path.expanduser('~/script/variation/data/ACD6/ACD6.non_Kz10.sequenced.del_NA_trimed.fasta')
AnalyzeSNPData.filterColumnsOfFasta(fasta_input_fname, output_fname, chromosome=4, start=8293290)
	"""
	
	@classmethod
	def removeInsertionSNPsFromSNPMatrix(cls, input_fname, output_fname):
		"""
		2009-5-29
			The most comprehensive way to represent a SNP is 'chr_pos_offset'.
			However, some of my programs can only handle 'chr_pos'.
			This function removes those SNPs embedded in insertion (relative to reference genome).
			
		"""
		sys.stderr.write("Removing SNPs embedded in insertions ")
		from pymodule import SNPData
		snpData1 = SNPData(input_fname=input_fname, turn_into_array=1, ignore_2nd_column=1)
		col_id_to_be_kept_ls = []
		for col_id in snpData1.col_id_ls:
			col_id_tuple = col_id.split('_')
			offset = int(col_id_tuple[2])
			if offset==0:
				col_id_to_be_kept_ls.append(col_id)
		snpData2 = SNPData.keepColsByColID(snpData1, col_id_to_be_kept_ls)
		snpData2.tofile(output_fname)
	
	"""
input_fname = os.path.expanduser('~/script/variation/data/ACD6/ACD6.non_Kz10.sequenced.del_NA_trimed.with_ref.clustalx.snp_matrix.84_1530.tsv')
output_fname = os.path.expanduser('~/script/variation/data/ACD6/snp_matrix.84_1530.no_insert.tsv')
AnalyzeSNPData.removeInsertionSNPsFromSNPMatrix(input_fname, output_fname)
	"""
	
	@classmethod
	def filterSNPMatrixBeforeImputation(cls, input_fname, output_fname):
		"""
		2009-5-29
			1. convert het to NA
			2. remove SNPs with >25% missing calls
			3. remove SNPs with >2 alleles cuz NPUTE doesn't support >2 alleles.
			4. remove SNPs with MAF<4%
		"""
		sys.stderr.write("Removing SNPs in preparation for imputation ...")
		from pymodule import SNPData
		snpData1 = SNPData(input_fname=input_fname, turn_into_array=1, ignore_2nd_column=1)
		snpData2 = SNPData.convertHetero2NA(snpData1)
		snpData3 = SNPData.removeColsByNARate(snpData2, max_NA_rate=0.25)
		snpData4 = SNPData.removeSNPsWithMoreThan2Alleles(snpData3)
		snpData5 = SNPData.removeColsByMAF(snpData4, min_MAF=0.04)
		snpData5.tofile(output_fname)
		sys.stderr.write("Done.\n")
	
	"""
input_fname = os.path.expanduser('~/script/variation/data/ACD6/snp_matrix.84_1530.no_insert.tsv')
output_fname = os.path.expanduser('~/script/variation/data/ACD6/snp_matrix.84_1530.no_insert.no_het.maxNA0.25.noSNPsMoreThan2Alleles.minMAF0.04.tsv')
AnalyzeSNPData.filterSNPMatrixBeforeImputation(input_fname, output_fname)
	"""
	
	@classmethod
	def generateDatasetWithImputedCallsOnly(cls, unImputedFname, imputedFname, output_fname):
		"""
		2009-5-22
			generate a dataset which only contains imputed calls and masks everything else as NA (missing) 
				unImputedFname and imputedFname are of same shape, with same accessions & SNPs.
				The former contains the data before imputation. The latter contains the one after imputation. 
		"""
		sys.stderr.write("Calculating stats of imputed calls in the dataset ... \n")
		from pymodule import SNPData
		snpData1 = SNPData(input_fname=unImputedFname, turn_into_array=1)
		snpData2 = SNPData(input_fname=imputedFname, turn_into_array=1)
		row_id2stats = {}
		
		import numpy, copy
		newSnpData = SNPData(col_id_ls=copy.deepcopy(snpData1.col_id_ls), row_id_ls=snpData1.row_id_ls)
		newSnpData.data_matrix = numpy.zeros(snpData1.data_matrix.shape, numpy.int8)
		row_index = 0
		for row_id, row_index1 in snpData1.row_id2row_index.iteritems():
			row_index2 = snpData2.row_id2row_index[row_id]
			if row_id not in row_id2stats:
				row_id2stats[row_id] = 0
			for col_id, col_index1 in snpData1.col_id2col_index.iteritems():
				col_index2 = snpData2.col_id2col_index[col_id]
				allele_unimputed = snpData1.data_matrix[row_index1][col_index1]
				allele_imputed = snpData2.data_matrix[row_index2][col_index2]
				if (allele_unimputed==0 or allele_unimputed==-2) and allele_imputed!=0 and allele_imputed!=-2:
					#If before imputation it's NA (missing), keep the imputed call.
					row_id2stats[row_id] += 1
					newSnpData.data_matrix[row_index1][col_index1] = allele_imputed
		newSnpData.tofile(output_fname)
		sys.stderr.write("Done.\n")
		return row_id2stats
	
	"""
unImputedFname = '/Network/Data/250k/db/dataset/call_method_35.tsv'
imputedFname = '/Network/Data/250k/db/dataset/call_method_33.tsv'
output_fname = '/tmp/call_method_33_imputed_only.tsv'
row_id2no_of_imputed = AnalyzeSNPData.generateDatasetWithImputedCallsOnly(unImputedFname, imputedFname, output_fname)
	"""
	
	@classmethod
	def calculateNumberOfDifferingAllelesFromOneEcotypeInCallMethod(cls, db_250k, call_method_id, ecotype_id=6909):
		"""
		2010-9-29
			Paul Dervent's email regarding too many non-reference alleles disappearing
			 from version 3.02 (call 43) to 3.04 (call 54)
		"""
		sys.stderr.write("Calculating number alleles from non-%s ecotypes to %s ... \n"%(ecotype_id, ecotype_id))
		import Stock_250kDB
		cm = Stock_250kDB.CallMethod.get(call_method_id)
		from pymodule import SNPData, TwoSNPData, PassingData
		snpData1 = SNPData(input_fname=cm.filename, turn_into_array=1)
		no_of_differing_alleles = 0
		no_of_cols = len(snpData1.col_id2col_index)
		base_ecotype_row_index = None
		for row_id, row_index in snpData1.row_id2row_index.iteritems():
			row_ecotype_id = int(row_id[0])
			if row_ecotype_id==ecotype_id:
				base_ecotype_row_index = row_index
				break
		
		no_of_non_base_accessions = 0
		for row_id, row_index in snpData1.row_id2row_index.iteritems():
			row_ecotype_id = int(row_id[0])
			if row_ecotype_id!= ecotype_id:
				no_of_non_base_accessions += 1
				no_of_identities = sum(snpData1.data_matrix[row_index,:]==snpData1.data_matrix[base_ecotype_row_index,:])
				no_of_differing_alleles += (no_of_cols - no_of_identities)
		
		no_of_total_calls = len(snpData1.row_id_ls)*no_of_cols
		fraction_of_differing = no_of_differing_alleles/float(no_of_total_calls)
		print "%s (%.4f) alleles from %s rows different from ecotype %s."%(no_of_differing_alleles, fraction_of_differing, \
											no_of_non_base_accessions, ecotype_id)
	
	"""
		# 2010-9-29
		AnalyzeSNPData.calculateNumberOfDifferingAllelesFromOneEcotypeInCallMethod(db_250k, 54, ecotype_id=6909)
		sys.exit(0)
		
		# 2010-9-29
		AnalyzeSNPData.calculateNumberOfDifferingAllelesFromOneEcotypeInCallMethod(db_250k, 43, ecotype_id=6909)
		sys.exit(0)
		
	"""
	
	@classmethod
	def cmpTwoSNPDatasets(cls, inputFname1, inputFname2):
		"""
		2009-6-12
			compare two SNP datasets, report:
				#rows deleted/added
				#columns deleted/added
				#SNPs changed (from which to which)
		"""
		sys.stderr.write("Comparing two SNP datasets ... \n")
		from pymodule import SNPData, TwoSNPData, PassingData
		snpData1 = SNPData(input_fname=inputFname1, turn_into_array=1)
		snpData2 = SNPData(input_fname=inputFname2, turn_into_array=1)
		row_id_deleted = []
		row_id_added = []
		col_id_deleted = []
		col_id_added = []
		total_row_id_set = set(snpData1.row_id_ls) | set(snpData2.row_id_ls)
		total_col_id_set = set(snpData1.col_id_ls) | set(snpData2.col_id_ls)
		for row_id in total_row_id_set:
			if row_id not in snpData1.row_id2row_index:
				row_id_added.append(row_id)
			elif row_id not in snpData2.row_id2row_index:
				row_id_deleted.append(row_id)
		for col_id in total_col_id_set:
			if col_id not in snpData1.col_id2col_index:
				col_id_added.append(col_id)
			elif col_id not in snpData2.col_id2col_index:
				col_id_deleted.append(col_id)
		
		twoSNPData = TwoSNPData(SNPData1=snpData1, SNPData2=snpData2)
		diff_data = twoSNPData.get_diff_matrix()
		diff_matrix = diff_data[0]
		import numpy
		sth_index_ls = [1] + range(3, diff_matrix.shape[1])	#"deletion" + 10 calls
		print sth_index_ls
		print diff_matrix
		counter_NA_to_sth = numpy.sum(diff_matrix[[0,2],:][:, sth_index_ls])
		counter_sth_to_NA = numpy.sum(diff_matrix[sth_index_ls,:][:, [0,2]])
		sub_diff_matrix = diff_matrix[sth_index_ls,:][:, sth_index_ls]
		counter_sth_to_sth_diff = numpy.sum(sub_diff_matrix) - numpy.sum(numpy.diagonal(sub_diff_matrix))
		
		print "%s rows deleted."%len(row_id_deleted)
		print "%s rows added."%len(row_id_added)
		print "%s cols deleted."%len(col_id_deleted)
		print "%s cols added."%len(col_id_added)
		print "%s SNPs from NA to sth."%counter_NA_to_sth
		print "%s SNPs from sth to NA."%counter_sth_to_NA
		print "%s SNPs from sth to sth different."%counter_sth_to_sth_diff
		return_data = PassingData(row_id_deleted=row_id_deleted, row_id_added=row_id_added, col_id_deleted=col_id_deleted, col_id_added=col_id_added)
		return_data.counter_NA_to_sth = counter_NA_to_sth
		return_data.counter_sth_to_NA = counter_sth_to_NA
		return_data.counter_sth_to_sth_diff = counter_sth_to_sth_diff
		sys.stderr.write("Done.\n")
		return return_data
	
	"""
	inputFname1 = '/Network/Data/250k/db/dataset/call_method_35.tsv'
	inputFname2 = '/Network/Data/250k/db/dataset/call_method_33.tsv'
	return_data = AnalyzeSNPData.cmpTwoSNPDatasets(inputFname1, inputFname2)
	"""
	
	@classmethod
	def cmpOneRowToTheOther(cls, inputFname, row_id1, row_id2):
		"""
		2009-6-17
			compare SNP data of one accession to the other in the same dataset
		"""
		sys.stderr.write("Comparing one row to the other ... \n")
		from pymodule import SNPData, TwoSNPData, PassingData
		row_id_key_set = set([row_id1, row_id2])
		snpData = SNPData(input_fname=inputFname, turn_into_array=1, row_id_key_set=row_id_key_set)
		twoSNPData = TwoSNPData(SNPData1=snpData, SNPData2=snpData)
		print twoSNPData.cmpOneRow(row_id1, row_id2)
	
	"""
	inputFname = '/Network/Data/250k/db/dataset/call_method_29.tsv'
	row_id1 = ('6910', '62')
	row_id2 = ('8290', '181')
	AnalyzeSNPData.cmpOneRowToTheOther(inputFname, row_id1, row_id2)
		
	inputFname = os.path.expanduser('~/mnt2/panfs/NPUTE_data/input/250k_l3_y.85_20091208.tsv')
	row_id1 = ('7034', '1338')	# Blh-1 from Versailles plate
	#row_id2 = ('8265', '243')	# another Blh-1
	row_id2 = ('7035', '336')	# Blh-2
	AnalyzeSNPData.cmpOneRowToTheOther(inputFname, row_id1, row_id2)
	"""
	
	@classmethod
	def cmpAllDuplicatesOfOneEcotype(cls, inputFname, ecotype_id_ls):
		"""
		2009-12-11
			For each ecotype_id in ecotype_id_ls, compare mismatch-rates between duplicates
		"""
		sys.stderr.write("Comparing one row to the other ... \n")
		from pymodule import SNPData, TwoSNPData, PassingData
		ecotype_id_set = set(ecotype_id_ls)
		def row_id_hash_func(row_id):
			return int(row_id[0])
		snpData = SNPData(input_fname=inputFname, turn_into_array=1, row_id_key_set=ecotype_id_set, row_id_hash_func=row_id_hash_func)
		
		ecotype_id2row_id_to_check_ls = {}
		for row_id in snpData.row_id_ls:
			ecotype_id = int(row_id[0])
			if ecotype_id in ecotype_id_set:
				if ecotype_id not in ecotype_id2row_id_to_check_ls:
					ecotype_id2row_id_to_check_ls[ecotype_id] = []
				ecotype_id2row_id_to_check_ls[ecotype_id].append(row_id)
		twoSNPData = TwoSNPData(SNPData1=snpData, SNPData2=snpData)
		for ecotype_id, row_id_to_check_ls in ecotype_id2row_id_to_check_ls.iteritems():
			if len(row_id_to_check_ls)>1:
				print "ecotype_id: %s"%ecotype_id
				no_of_arrays = len(row_id_to_check_ls)
				for i in range(no_of_arrays):
					for j in range(i+1, no_of_arrays):
						row_id1 = row_id_to_check_ls[i]
						row_id2 = row_id_to_check_ls[j]
						print "row_id1 %s vs row_id2 %s"%(row_id1, row_id2)
						print twoSNPData.cmpOneRow(row_id1, row_id2)
	
	"""
	inputFname = os.path.expanduser('~/mnt2/panfs/NPUTE_data/input/250k_l3_y.85_20091208.tsv')
	ecotype_id_ls = [8297, 7317, 6910, 8274, 6911, 6905, 7034, 6909, 6962, 7373, 7270, 6983, 6899]
	AnalyzeSNPData.cmpAllDuplicatesOfOneEcotype(inputFname, ecotype_id_ls)
	
	"""
	
	@classmethod
	def linkEcotypeIDFromSuziPhenotype(cls, fname_with_ID, fname_with_phenotype, output_fname):
		"""
		2009-7-31
			she gave me two files
				one has phenotype data and accession names but with no ecotype ID
				2nd is a map from accession name to ecotype ID
		"""
		sys.stderr.write("Linking accession names to ecotype ID ... ")
		import csv
		inf_phenotype = csv.reader(open(fname_with_phenotype, 'r'), delimiter='\t')
		accession_name_ls = []
		#skip two lines
		inf_phenotype.next()
		inf_phenotype.next()
		for row in inf_phenotype:
			accession_name_ls.append(row[0])
		del inf_phenotype
		
		inf_with_ID = csv.reader(open(fname_with_ID), delimiter='\t')
		inf_with_ID.next()
		accession_name2ecotype_id = {}
		for row in inf_with_ID:
			ecotype_id = row[0]
			accession_name = row[5]
			accession_name2ecotype_id[accession_name] = ecotype_id
		del inf_with_ID
		print accession_name2ecotype_id
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		for accession_name in accession_name_ls:
			ecotype_id = accession_name2ecotype_id.get(accession_name)
			writer.writerow([accession_name, ecotype_id])
		del writer
		
		sys.stderr.write("Done.\n")
	
	"""
fname_with_ID = '/tmp/seed_batch3_out.txt'
fname_with_phenotype = '/tmp/ft_LN_other phenotypes_10_batch_3_complete_7_29th.txt'
output_fname = '/tmp/batch_3_name2id.tsv'
AnalyzeSNPData.linkEcotypeIDFromSuziPhenotype(fname_with_ID, fname_with_phenotype, output_fname)
	"""

class Fun(object):
	"""
	2010-5-20
		fun projects
	"""
	@classmethod
	def printFrameshiftTripletDict(cls, frameshift2triplet2count):
		"""
		"""
		for frameshift, triplet2count in frameshift2triplet2count.iteritems():
			print "frameshift:", frameshift
			count_triplet_ls = []
			for triplet, count in triplet2count.iteritems():
				count_triplet_ls.append((count, triplet))
			count_triplet_ls.sort()
			count_triplet_ls.reverse()
			for count, triplet in count_triplet_ls:
				triplet = triplet.replace(' ', '|')
				print '%s\t%s'%(count, triplet)
	
	@classmethod
	def findLetterFreqSpectrum(cls, name):
		"""
		"""
		letter2count = {}
		for i in range(len(name)):
			#letter = letter.upper()
			letter = name[i]
			letter = letter.upper()
			if letter not in letter2count:
				letter2count[letter] = 0
			letter2count[letter] += 1
		count_letter_ls = [(count, letter) for letter, count in letter2count.iteritems()]
		count_letter_ls.sort()
		count_letter_ls.reverse()
		return [row[0] for row in count_letter_ls]
		
	@classmethod
	def findDNATripletFreqSpectrum(cls, dna_seq):
		"""
		"""
		no_of_triplets = len(dna_seq)/3 + 2
		letter2count = {}
		for i in range(no_of_triplets):
			#letter = letter.upper()
			start_index = i*3
			if start_index<len(dna_seq):
				triplet = dna_seq[start_index:start_index+3]
				if len(triplet)!=3:
					continue
				if triplet not in letter2count:
					letter2count[triplet] = 0
				letter2count[triplet] += 1
		count_letter_ls = [(count, letter) for letter, count in letter2count.iteritems()]
		count_letter_ls.sort()
		count_letter_ls.reverse()
		return [row[0] for row in count_letter_ls]
	
	@classmethod
	def matchNameWithDNAByTriplet(cls, name, dna_seq):
		"""
		2010-5-21
		"""
		no_of_triplets = len(dna_seq)/3
		triplet2letter = {}
		translated_name = ''
		for i in range(no_of_triplets):
			start_index = i*3
			triplet = dna_seq[start_index:start_index+3]
			if triplet not in triplet2letter:
				triplet2letter[triplet] = name[i]
			translated_name += triplet2letter[triplet]
		if name==translated_name:
			return triplet2letter
		else:
			return None
		
	@classmethod
	def crackDGibson2010ScienceWaterMarkCode(cls):
		"""
		2010-5-20
			Dan's new science paper came out today. Chris forwarded me the email from laura.
			
			compare the frequency spectrum of triplets in DNA sequence and the frequency of letters in names.
			
			The stop codon corresponds to space (or should be).
			
			Match the DNA blocks separated by the stop-codon with names (first & last).
			
				Their length and frequency spectrum should match for each (DNA-block, name) pair.
				
			the beginning&trailing 6-frame stop-codon and restricting site have been removed from s1-s4.
			
			
		"""
		s1 = 'GTTCGAATATTTCTATAGCTGTACATATTGTAATGCTGATAACTAATACTGTGCGCTTGACTGTGATCCTGATAAATAACTTCTTCTGTAGGGTAGAGTTTTA\
TTTAAGGCTACTCACTGGTTGCAAACCAATGCCGTACATTACTAGCTTGATCCTTGGTCGGTCATTGGGGGATATCTCTTACTAATAGAGCGGCCTATCGCGTATTCTCGCCG\
GACCCCCCTCTCCCACACCAGCGGTGTAGCATCACCAAGAAAATGAGGGGAACGGATGAGGAACGAGTGGGGGCTCATTGCTGATCATAATGACTGTTTATATACTAATGC\
CGTCAACTGTTTGCTGTGATACTGTGCTTTCGAGGGCGGGAGATTCGTTTTTGACATACATAAATATCATGACAAAACAGCCGGTCATGACAAAACAGCCGGTCATAATAGAT\
TAGCCGGTGACTGTGAAACTAAAGCTACTAATGCCGTCAATAAATATGATAATAGCAACGGCACTGACTGTGAAACTAAAGCCGGCACTCATAATAGATTAGCCGGAGTCGT\
ATTCATAGCCGGTAGATATCACTATAAGGCCCAGGATCATGATGAACACAGCACCACGTCGTCGTCCGAGTTTTTTTGCTGCGACGTCTATACCACGGAAGCTGATCATAAAT\
AGTTTTTTTGCTGCGGCACTAGAGCCGGACAAGCACACTACGTTTGTAAATACATCGTTCCGAATTGTAAATAATTTAATTTCGTATTTAAATTATATGATCACTGGCTATAGTC\
TAGTGATAACTACAATAGCTAGCAATAAGTCATATATAACAATAGCTGAACCTGTGCTACATATCCGCTATACGGTAGATATCACTATAAGGCCCAGGACAATAGCTGAACTGA\
CGTCAGCAACTACGTTTAGCTTGACTGTGGTCGGTTTTTTTGCTGCGACGTCTATACGGAAGCTCATAACTATAAGAGCGGCACTAGAGCCGGCACACAAGCCGGCACAGT\
CGTATTCATAGCCGGCACTCATGACAAAACAGC'
		
		s2 = 'CAACTGGCAGCATAAAACATATAGAACTACCTGCTATAAGTGATACAACTGTTTTCATAGTAAAACATACAACGTTGCTGATAGTACTCCTAAGTGATAGCTT\
AGTGCGTTTAGCATATATTGTAGGCTTCATAATAAGTGATATTTTAGCTACGTAACTAAATAAACTAGCTATGACTGTACTCCTAAGTGATATTTTCATCCTTTGCAATACAATAA\
CTACTACATCAATAGTGCGTGATATGCCTGTGCTAGATATAGAACACATAACTACGTTTGCTGTTTTCAGTGATATGCTAGTTTCATCTATAGATATAGGCTGCTTAGATTCCCT\
ACTAGCTATTTCTGTAGGTGATATACGTCCATTGCATAAGTTAATGCATTTAACTAGCTGTGATACTATAGCATCCCCATTCCTAGTGCATATTTTCATCCTAGTGCTACGTGAT\
ATAATTGTACTAATGCCTGTAGATAATTTAATGCCTGGCTCGTTTGTAGGTGATAATTTAGTGCCTGTAAAACATATACCTGAGTGCTCGTTGCGTGATAGTTCGTTCATGCAT\
ATACAACTAGGCTGCTGTGATATGGTCACTGCCCTTACTGTGCTACATATTACTGCGAGGGGGATGACGTATAAACCTGTTGTAAGTGATATGACGTATATAACTACTAGTGA\
TATGACGTATAGGCTAGAACAACGTGATATGACGTATATGACTACTGTCCCAAACATCAGTGATATGACGTATACTATAATTTCTATAATAGTGATAAATAAACCTGGGCTAAA\
TACGTTCCTGAATACGTGGCATAAACCTGGGCTAACGAGGAATACCCATAGTTTAGCAATAAGCTATAGTTCGTCATTTTTAA'
		s3 = 'TTTAACCATATTTAAATATCATCCTGATTTTCACTGGCTCGTTGCGTGATATAGATTCTACTGTAGTGCTAGATAGTTCTGTACTAGGTGATACTATAGATTTC\
ATAGATAGCACTACTGGCTTCATGCTAGGCATCCCAATAGCTAGTGATAGTTTAGTGCATACAACGTCATGTGATACAACGTTGCTGGCTGTAGATACAACGTCGTATTCTGT\
AAGTGATACAATAGCTATTGCTGTGCATAGGCCTATAGTGGCTGTAACTAGTGATATCACGTAACAACCATATAAGTTAGATTTAATGCCCCTGACTGAACGCTCGTTGCGTG\
ATAGTTTAGGCTCGTTGCATACAACTGTGATTTTCATAAAACAACGTGATAATTTAGTGCTAGATAAGTTCCGCTTAGCAAGTGATAGTTTCCGCTTGACTGTGCATAGTTCGT\
TCATGCGCTCGTTGCGTGATAAACTAGGCAGCTTCACAACTGATAATTTAATTGCTGATATTGCTGGCTGTCTAGTGCTAGTGATCATAGTGCGTGATAGTTTAAGCTGCTCT\
GTTTTAGATATCACGTGCTTGATAATGAAACTAACTAGTGATACTACGTAGTTAACTATGAATAGGCCTACTGTAAATTCAATAGTGCGTGATATTGAACTAGATTCTGCAACTG\
CTAATATGCCGTGCTGCACGTTTGGTGATAGTTTAGCATGCTTCACTATAATAAATATGGTAGTTGTAACTACTGCGAATAGGGGGAGCTTAATAAATATGATCACTGTGCTAC\
GCTATATGCCGTTGAATATAGGCTATATGATCATAACATATATAGCTATAAGTGATAAGTTCCTGAATATAGGCTATATGATCATAACATATACAACTGTACTCATGAATAAGTT\
AACGAGGA'
		s4 = 'TTTCATTGCTGATCACTGTAGATATAGTGCATTCTATAAGTCGCTCCCACAGGCTAGTGCTGCGCACGTTTTTCAGTGATATTATCCTAGTGCTACATAACA\
TCATAGTGCGTGATAAACCTGATACAATAGGTGATATCATAGCAACTGAACTGACGTTGCATAGCTCAACTGTGATCAGTGATATAGATTCTGATACTATAGCAACGTTGCGT\
GATATTTTCACTACTGGCTTGACTGTAGTGCATATGATAGTACGTCTAACTAGCATAACTAGTGATAGTTATATTTCTATAGCTGTACATATTGTAATGCTGATAACTAGTGATA\
TAATCCAACTAGATAGTCCTGAACTGATCCCTATGCTAACTAGTGATAAACTAACTGATACATCGTTCCTGCTACGTGATAGCTTCACTGAGTTCCATACATCGTCGTGCTTAA\
ACATCAGTGATAACACTATAGAGTTCATAGATACTGCATTAACTAGTGATATGACTGCAAATAGCTTGACGTTTTGCAGTCTAAAACAACGTGATAATTCTGTAGTGCTAGATA\
CTATAGATTTCCTGCTAAGTGATAAGTCTACTGATTTACTAATGAATAGCTTGGTTTTGGCATACACTGTGCGCTGCACTGGTGATAGCTTTTCGTTGATGAATAATTTCCCTA\
GCACTGTGCGTGATATGCTAGATTCTGTAGATAGGCTAAATTCGTCTACGTTTGTAGGTGATAGTTTAGTTGCTGTAACTAATATTATCCCTGTGCCGTTGCTAAGCTGTGATA\
TCATAGTGCTGCTAGATATGATAAGCAAACTAATAGAGTCGAGGGGGAGTCTCATAGTGAATACTGATATTTTAGTGCTGCCGTTGAATAAGTTCCCTGAACATTGTGATACT\
GATATTTTAGTGCTGCCGTTGAATATCCTGCATTTAACTAGCTTGATAGTGCATTCGAGGAATACCCATACTACTGTTTTCATAGCTAATTATAGGCTAACATTGCCAATAGTGC'
		s1 = s1+s2+s3+s4
		print "length of watermark 1:", len(s1)
		authors = 'Daniel G. Gibson,1 John I. Glass,1 Carole Lartigue,1 Vladimir N. Noskov,1 Ray-Yuan Chuang,1 Mikkel A. \
Algire,1 Gwynedd A. Benders,2 Michael G. Montague,1 Li Ma,1 Monzia M. Moodie,1 Chuck Merryman,1Sanjay Vashee,1 \
Radha Krishnakumar,1 Nacyra Assad-Garcia,1 Cynthia Andrews-Pfannkoch,1 Evgeniya A. Denisova,1 Lei Young,1 Zhi-Qing Qi,1 \
Thomas H. Segall-Shapiro,1 Christopher H. Calvey,1 Prashanth P. Parmar,1 Clyde A. Hutchison III,2 Hamilton O. Smith,2 \
J. Craig Venter1,2*\
To live, to err, to fall, to triumph, to recreate life out of life. James Joyce See things not as they are, but as they might be. What I cannot built, I cannot understand.'
		authors = authors.upper() + "ABCDEFGHIJKLMNOPQRSTUVWXYZ abcdefghijklmnopqrstuvwxyz "
		
		authors_1 = authors.replace(",", " ")
		authors_1 = authors_1.replace("1", " ")
		authors_1 = authors_1.replace("2", " ")
		authors_1 = authors_1.replace("*", " ")
		author_list = authors_1.split()
		
		stop_codon_set = set(['TAG', 'TAA', 'TGA'])
		frameshift2stop_codon_counts = {}
		frameshift2triplet2count = {}
		frameshift2double_triplet2count = {}
		new_s1 = ""	#replace the stop codon with empty strings
		for frameshift in [0,]:	# 0,1,and 2 have been tested. 0 harbors the most stop codons.
			frameshift2stop_codon_counts[frameshift] = 0
			frameshift2triplet2count[frameshift] = {}
			no_of_triplets = len(s1)/3 + 2
			for i in range(no_of_triplets):
				start_index = i*3+frameshift
				if start_index<len(s1):
					triplet = s1[start_index:start_index+3]
					if len(triplet)!=3:
						continue
					if triplet in stop_codon_set:
						frameshift2stop_codon_counts[frameshift] += 1
						new_s1 += "   "
					else:
						new_s1 += triplet
					if triplet not in frameshift2triplet2count[frameshift]:
						frameshift2triplet2count[frameshift][triplet] = 0
					frameshift2triplet2count[frameshift][triplet] += 1
		print frameshift2stop_codon_counts
		cls.printFrameshiftTripletDict(frameshift2triplet2count)
		
		new_s1 = s1
		new_s1.replace('ATA', ' ')
		dna_block_ls = new_s1.split()
		author2dna_block_ls = {}
		author2triplet2letter_ls = {}
		for author in author_list:
			if len(author)<5:
				continue
			#for i in range(len(new_s1)):
			for dna_block in dna_block_ls:
				name_length = len(author)
				no_of_bases = name_length*3
				for i in range(len(dna_block)-no_of_bases+1):
					dna_seq = dna_block[i:i+no_of_bases]
					#name_freq_ls = cls.findLetterFreqSpectrum(author)
					#dna_triplet_freq_ls = cls.findDNATripletFreqSpectrum(dna_block)
					#if name_freq_ls==dna_triplet_freq_ls:
					triplet2letter = cls.matchNameWithDNAByTriplet(author, dna_seq)
					if triplet2letter is not None:
						if author not in author2dna_block_ls:
							author2dna_block_ls[author] = []
							author2triplet2letter_ls[author] = []
						author2dna_block_ls[author].append(dna_seq)
						author2triplet2letter_ls[author].append([triplet2letter, 0])
		for author, dna_block_ls in author2dna_block_ls.iteritems():
			print author, ":", dna_block_ls
		
		"""
		1. build a graph. each node is a triplet2letter. edge if two conflicting triplet2letter conflict with each other.
		2. eliminate all edges by a greedy algorithm. the node with most edges go first.
		"""
		
		"""
		#
		for frameshift in [0,]:	# 0,1,and 2 have been tested. 0 harbors the most stop codons.
			frameshift2double_triplet2count[frameshift] = {}
			no_of_triplets = len(new_s1)/3 + 2
			for i in range(no_of_triplets):
				start_index = i*3+frameshift
				if start_index<len(new_s1):
					triplet = new_s1[start_index:start_index+3]
					if len(triplet)!=3:
						continue
					next_triplet = new_s1[start_index+3:start_index+6]
					double_triplet = triplet + next_triplet
					if double_triplet not in frameshift2double_triplet2count[frameshift]:
						frameshift2double_triplet2count[frameshift][double_triplet] = 0
					frameshift2double_triplet2count[frameshift][double_triplet] += 1
		
		cls.printFrameshiftTripletDict(frameshift2double_triplet2count)
		"""
		
		
		letter2count = {}
		for i in range(len(authors)):
			#letter = letter.upper()
			letter = authors[i]
			letter = letter.upper()
			if letter not in letter2count:
				letter2count[letter] = 0
			letter2count[letter] += 1
		count_letter_ls = [(count, letter) for letter, count in letter2count.iteritems()]
		count_letter_ls.sort()
		count_letter_ls.reverse()
		for count, letter in count_letter_ls:
			letter = letter.replace(' ', '|')
			print '%s\t%s'%(count, letter)
		
		"""
		new_dict = {"ATA":"a", "GTT":"n"}
		no_of_triplets = len(new_s1)/3 + 2
		for i in range(no_of_triplets):
			start_index = i*3+frameshift
			if start_index<len(new_s1):
				triplet = new_s1[start_index:start_index+3]
				if len(triplet)!=3:
					continue
				alphabet = new_dict.get(triplet)
				if alphabet is not None:
					sys.stdout.write("-%s-"%alphabet)
				else:
					sys.stdout.write("%s"%triplet)
		print
		print
		
		
		new_dict = {"ATACAA":"an"}
		no_of_triplets = len(new_s1)/3 + 2
		for i in range(no_of_triplets):
			start_index = i*3+frameshift
			if start_index<len(new_s1):
				triplet = new_s1[start_index:start_index+3]
				if len(triplet)!=3:
					continue
				
				next_triplet = new_s1[start_index+3:start_index+6]
				double_triplet = triplet + next_triplet
				alphabet = new_dict.get(double_triplet)
				if alphabet is not None:
					sys.stdout.write("-%s-"%alphabet)
				else:
					sys.stdout.write("%s"%triplet)
					
		"""
	"""
	Fun.crackDGibson2010ScienceWaterMarkCode()
	sys.exit(0)
	"""
	
class DBGenome(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['genome', 'd', 1, 'database name', ],\
							('schema', 0, ): ['', 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2010-6-22
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		from pymodule import GenomeDB
		self.db = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		self.db.setup(create_tables=False)	#2010-6-22
		
	def filterValue(self, value, data_type=None, NA_str="-"):
		"""
		2010-6-22
			NCBI's gene_info.gz has '-' as NA.
		"""
		if value==NA_str:
			value = None
		if value is not None and data_type is not None:
			value = data_type(value)
		return value
	
	def fillTableGene(self, inputFname):
		"""
		2010-6-22
			inputFname could be gene_info.gz or not.
			downloaded from ftp://ftp.ncbi.nlm.nih.gov/gene/
			
			plone doc /research/annot/tf/genomedb/README has information.
		"""
		sys.stderr.write("Filling table gene with data from %s ...\n"%(inputFname))
		session = self.db.session
		import csv
		from pymodule import figureOutDelimiter, getColName2IndexFromHeader
		if inputFname[-2:]=='gz':	#turn on the gzip flag based on the filename extension
			import gzip
			inf = gzip.open(inputFname, 'r')
		else:
			inf = open(inputFname, 'r')
		reader = csv.reader(inf, delimiter='\t')
		header = reader.next()
		"""
		## 2010-6-22 header looks like this.
		#Format: tax_id GeneID Symbol LocusTag Synonyms dbXrefs chromosome map_location description
		type_of_gene Symbol_from_nomenclature_authority Full_name_from_nomenclature_authority
		Nomenclature_status Other_designations Modification_date (tab is used as a separator, pound sign - start of a comment)
		"""
		from pymodule import GenomeDB
		from datetime import datetime
		header = header[0].split(' ')[1:]
		col_name2index = getColName2IndexFromHeader(header)
		counter = 0
		real_counter = 0
		for row in reader:
			counter += 1
			if row[0][0]=='#':	#ignore comments
				continue
			gene_id = self.filterValue(row[col_name2index['GeneID']], int)
			gene = GenomeDB.Gene.get(gene_id)
			if gene:	# already in db
				continue
			gene = Gene()
			gene.tax_id = self.filterValue(row[col_name2index['tax_id']], int)
			gene.gene_id = gene_id
			gene.gene_symbol = self.filterValue(row[col_name2index['Symbol']])
			gene.locustag = self.filterValue(row[col_name2index['LocusTag']])
			gene.synonyms = self.filterValue(row[col_name2index['Synonyms']])
			gene.dbxrefs = self.filterValue(row[col_name2index['dbXrefs']])
			gene.chromosome = self.filterValue(row[col_name2index['chromosome']])
			gene.map_location = self.filterValue(row[col_name2index['map_location']])
			gene.description = self.filterValue(row[col_name2index['description']])
			gene.type_of_gene = self.filterValue(row[col_name2index['type_of_gene']])
			gene.symbol_from_nomenclature_authority = self.filterValue(row[col_name2index['Symbol_from_nomenclature_authority']])
			gene.full_name_from_nomenclature_authority = self.filterValue(row[col_name2index['Full_name_from_nomenclature_authority']])
			gene.nomenclature_status = self.filterValue(row[col_name2index['Nomenclature_status']])
			gene.other_designations = self.filterValue(row[col_name2index['Other_designations']])
			gene.modification_date = self.filterValue(row[col_name2index['Modification_date']])
			if gene.modification_date:
				gene.modification_date = datetime.strptime(gene.modification_date, '%Y%m%d')	#20080317
			real_counter += 1
			session.add(gene)
			if real_counter%5000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
				session.flush()
				session.expunge_all()	# 2010-6-22 sqlalchemy 0.6
		sys.stderr.write("%s%s\t%s\n"%('\x08'*80, counter, real_counter))
	
	"""
	hostname='localhost'
	dbname='genome'
	dbGenome = DBGenome(drivername=drivername, db_user=db_user,
					db_passwd=db_passwd, hostname=hostname, dbname=dbname, schema=schema)
	inputFname = os.path.expanduser('~/script/variation/data/genome/gene_info.gz')
	dbGenome.fillTableGene(inputFname)
	sys.exit(0)
	"""
	
	def fillTableGene2go(self, inputFname):
		"""
		2010-6-22
			inputFname is gene2go.gz, downloaded from ftp://ftp.ncbi.nlm.nih.gov/gene/
			
			plone doc /research/annot/tf/genomedb/README has information.
		"""
		session = self.db.session
		sys.stderr.write("Filling table gene2go with data from %s ...\n"%(inputFname))
		import csv
		from pymodule import figureOutDelimiter, getColName2IndexFromHeader
		if inputFname[-2:]=='gz':	#turn on the gzip flag based on the filename extension
			import gzip
			inf = gzip.open(inputFname, 'r')
		else:
			inf = open(inputFname, 'r')
		reader = csv.reader(inf, delimiter='\t')
		header = reader.next()
		"""
		## 2010-6-22 header looks like this.
#Format: tax_id GeneID GO_ID Evidence Qualifier GO_term PubMed Category (tab is used as a separator, pound sign - start of a comment)
		"""
		from pymodule import GenomeDB
		from datetime import datetime
		header = header[0].split(' ')[1:]
		col_name2index = getColName2IndexFromHeader(header)
		counter = 0
		real_counter = 0
		for row in reader:
			counter += 1
			if row[0][0]=='#':	#ignore comments
				continue
			db_entry = GenomeDB.Gene2go()
			db_entry.tax_id = self.filterValue(row[col_name2index['tax_id']], int)
			db_entry.gene_id = self.filterValue(row[col_name2index['GeneID']], int)
			db_entry.go_id = self.filterValue(row[col_name2index['GO_ID']])
			db_entry.evidence = self.filterValue(row[col_name2index['Evidence']])
			db_entry.go_qualifier = self.filterValue(row[col_name2index['Qualifier']])
			db_entry.go_description = self.filterValue(row[col_name2index['GO_term']])
			db_entry.pubmed_ids = self.filterValue(row[col_name2index['PubMed']])
			db_entry.category = self.filterValue(row[col_name2index['Category']])
			real_counter += 1
			session.add(db_entry)
			if real_counter%5000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
				session.flush()
				session.expunge_all()	# 2010-6-22 sqlalchemy 0.6
		sys.stderr.write("%s%s\t%s\n"%('\x08'*80, counter, real_counter))
	
	"""
		hostname='localhost'
		dbname='genome'
		dbGenome = DBGenome(drivername=drivername, db_user=db_user,
						db_passwd=db_passwd, hostname=hostname, dbname=dbname, schema=schema)
		inputFname = os.path.expanduser('~/script/variation/data/genome/gene2go.gz')
		dbGenome.fillTableGene2go(inputFname)
		sys.exit(0)
		
		hostname='localhost'
		dbname='genome'
		dbGenome = DBGenome(drivername=drivername, db_user=db_user,
						db_passwd=db_passwd, hostname=hostname, dbname=dbname, schema=schema)
		inputFname = os.path.expanduser('~/script/variation/data/genome/gene_info.gz')
		dbGenome.fillTableGene(inputFname)
		inputFname = os.path.expanduser('~/script/variation/data/genome/gene2go.gz')
		dbGenome.fillTableGene2go(inputFname)
		sys.exit(0)
	
	"""
	
	def getGeneFamily(self, session, family_name=None, family_type=None, super_family=None):
		"""
		2010-8-19
			used by addTEGenesAndFamilyInfo()
		"""
		from pymodule import GenomeDB
		family = GenomeDB.GeneFamily.query.filter_by(short_name=family_name).first()
		if not family:
			family = GenomeDB.GeneFamily(short_name=family_name, family_type=family_type)
			family.super_family = super_family
			session.add(family)
			session.flush()
		return family
	
	def addTEGenesAndFamilyInfo(self, TAIR_TE_URL='ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR9_genome_release/TAIR9_Transposable_Elements.txt', \
							tax_id=3702, type_of_gene ='TRANSPOSABLE_ELEMENT', family_type='TE'):
		"""
		2011-1-25
			change the inputFname to TAIR_TE_URL, which would be fetched on the fly.
		2010-8-19
			inputFname, TAIR9_Transposable_Elements.txt is downloaded from 
				ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR9_genome_release/TAIR9_Transposable_Elements.txt
				
			add genes into Gene, EntrezgeneMapping
			add family into GeneFamily, Gene2Family
			add new symbols into Gene_symbol2id
		"""
		sys.stderr.write("Adding TE genes/fragments and family info into db ...")
		from transfac.src.UpdateGenomeDB import UpdateGenomeDB
		inputFname, = UpdateGenomeDB.getInputFileList([TAIR_TE_URL], '/tmp/')
		
		import csv
		from pymodule import GenomeDB
		reader = csv.reader(open(inputFname), delimiter='\t')
		header = reader.next()
		no_of_total = 0
		no_of_into_db = 0
		no_of_genes_already_in_db = 0
		no_of_entrezgene_mappings_already_in_db = 0
		counter = 0
		session = self.db.session
		for row in reader:
			counter += 1
			transposon_name, orientation_is_5prime, Transposon_min_Start, Transposon_max_End,\
				Transposon_Family, Transposon_Super_Family = row[:6]
			chromosome = transposon_name[2]
			if orientation_is_5prime=='true':
				strand = '+1'
			else:
				strand = '-1'
			start = int(Transposon_min_Start)
			stop = int(Transposon_max_End)
			
			annot_assembly = GenomeDB.AnnotAssembly.query.filter_by(chromosome=chromosome).filter_by(tax_id=tax_id).first()
			if not annot_assembly:
				sys.stderr.write(" Warning: chromosome %s tax_id=%s not in db yet. all TEs from this chromosome are ignored.\n"%\
								(chromosome, tax_id))
				continue
			entrezgene_type = GenomeDB.EntrezgeneType.query.filter_by(type=type_of_gene).first()	#query the db to see if it exists or not
			if not entrezgene_type:
				entrezgene_type = GenomeDB.EntrezgeneType(type=type_of_gene)
				session.add(entrezgene_type)
				session.flush()
				no_of_into_db += 1
			
			gene = GenomeDB.Gene.query.filter_by(tax_id=tax_id).filter_by(chromosome=chromosome).\
				filter_by(strand=strand).filter_by(start=start).filter_by(stop=stop).first()
			if gene:
				if gene.gene_symbol!=transposon_name and gene.locustag!=transposon_name:
					gene_symbol2id = GenomeDB.Gene_symbol2id(tax_id=tax_id, gene_symbol=transposon_name, symbol_type='TE_name')
					gene_symbol2id.gene = gene
					session.add(gene_symbol2id)
					no_of_into_db += 1
				no_of_genes_already_in_db += 1
			else:
				gene = GenomeDB.Gene(gene_symbol=transposon_name, tax_id=tax_id, locustag=transposon_name, \
					chromosome=chromosome,\
					type_of_gene=type_of_gene, description='TE Family %s Super-Family %s'%(Transposon_Family, Transposon_Super_Family), \
					strand=strand, start=start, stop=stop)
				gene.genomic_annot_assembly = annot_assembly
				gene.entrezgene_type = entrezgene_type
				session.add(gene)
				no_of_into_db += 1
			
			"""
			# 2011-1-22	EntrezgeneMapping has been merged into Gene.
			if gene.gene_id:
				entrezgene_mapping = GenomeDB.EntrezgeneMapping.query.filter_by(gene_id=gene.gene_id).first()
			else:
				entrezgene_mapping = None
			if not entrezgene_mapping:
				entrezgene_mapping = GenomeDB.EntrezgeneMapping(tax_id=tax_id, chromosome=chromosome,\
												start=start, stop=stop, strand=strand,)
				entrezgene_mapping.genomic_annot_assembly = annot_assembly
				entrezgene_mapping.gene = gene
				entrezgene_mapping.entrezgene_type = entrezgene_type
				
				session.add(entrezgene_mapping)
				no_of_into_db += 1
			else:
				no_of_entrezgene_mappings_already_in_db += 1
			"""
			super_family = self.getGeneFamily(session, family_name=Transposon_Super_Family, family_type=family_type, \
									super_family=None)
			family = self.getGeneFamily(session, family_name=Transposon_Family, family_type=family_type, \
									super_family=super_family)
			gene2family = GenomeDB.Gene2Family(family=family)
			gene2family.gene = gene
			session.add(gene2family)
			no_of_into_db += 1
			
			if no_of_into_db>2000:
				session.flush()
				no_of_total += no_of_into_db
				no_of_into_db = 0
				sys.stderr.write('\t %s into db out of %s genes. %s gene(s) already_in_db. %s entrezgene_mapping(s) already_in_db.\n'%\
						(no_of_total, counter, no_of_genes_already_in_db,\
						no_of_entrezgene_mappings_already_in_db,))
		sys.stderr.write('\t %s into db out of %s genes. %s gene(s) already_in_db. %s entrezgene_mapping(s) already_in_db.\n'%\
						(no_of_total, counter, no_of_genes_already_in_db,\
						no_of_entrezgene_mappings_already_in_db,))
		session.flush()
		session.expunge_all()	# 2010-6-22 sqlalchemy 0.6
	"""
	
		# 2010-8-19
		hostname='banyan'
		dbname='genome_tair'
		dbGenome = DBGenome(drivername=self.drivername, db_user=self.db_user,
						db_passwd=self.db_passwd, hostname=hostname, dbname=dbname, schema=self.schema)
		#inputFname = os.path.expanduser('~/script/variation/data/TAIR9/TAIR9_genome_release/TAIR9_Transposable_Elements.txt')
		TAIR_TE_URL = 'ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR9_genome_release/TAIR9_Transposable_Elements.txt'
		TAIR_TE_URL = 'ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR9_Transposable_Elements.txt'
		dbGenome.addTEGenesAndFamilyInfo(TAIR_TE_URL, tax_id=3702, type_of_gene ='TRANSPOSABLE_ELEMENT')
		sys.exit(0)
		
		#2011-1-22
		hostname='banyan.usc.edu'
		dbname='genome'
		dbGenome = DBGenome(drivername=self.drivername, db_user=self.db_user,
						db_passwd=self.db_passwd, hostname=hostname, dbname=dbname, schema=self.schema)
		#inputFname = os.path.expanduser('~/script/variation/data/TAIR9/TAIR9_genome_release/TAIR9_Transposable_Elements.txt')
		TAIR_TE_URL = 'ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR9_genome_release/TAIR9_Transposable_Elements.txt'
		TAIR_TE_URL = 'ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR9_Transposable_Elements.txt'
		dbGenome.addTEGenesAndFamilyInfo(TAIR_TE_URL, tax_id=3702, type_of_gene ='TRANSPOSABLE_ELEMENT')
		sys.exit(0)
		
	"""
class Main(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('input_fname', 0, ): ['', 'i', 1, 'common input file.', ],\
							('output_fname', 0, ): ['', 'o', 1, 'common output file', ],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self, **keywords):
		"""
		2008-12-05
		This class is the entry to all others.
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def run(self,):
		if self.debug:	# 2010-4-18 enter debug mode "~/.../variation/misc.py -b"
			import pdb
			pdb.set_trace()
			debug = True
		else:
			debug =False
		
		
		import Stock_250kDB
		db_250k = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_250k.setup(create_tables=False)
		self.db_250k = db_250k
		
		#import MySQLdb
		#conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.db_user, passwd = self.db_passwd)
		#curs = conn.cursor()
		
		# 2011-5-4
		# priorTAIRVersion=True if it's connected to banyan's TAIR9 db because the association results are based off TAIR8.
		# call_method_id=None if you want the function to go through all possible files.
		DB250k.convertOldFormatResultMethodFileIntoNewFormat(db_250k, call_method_id=None, priorTAIRVersion=True)
		sys.exit(0)
		
		
		# 2011-4-22
		result1_id = 4634	#cnv_20_LD_KW
		result1_peak_type_id = 1 #min_score=4
		result2_id = 4635	#cnv_20_LDV_KW
		result2_peak_type_id = 1	#min_score=4
		CNV.comparePeaksFromTwoAssociationResults(db_250k, result1_id=result1_id, result1_peak_type_id=result1_peak_type_id,\
										result2_id=result2_id, result2_peak_type_id=result2_peak_type_id, result2_peak_ext_dist=0)
		sys.exit(0)


#2007-03-05 common codes to initiate database connection
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	#sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	#sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
	#sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/test/python')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/variation/src')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	#sys.path.insert(0, os.path.expanduser('~/lib/python'))
	#sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
	#sys.path.insert(0, os.path.join(os.path.expanduser('~/script/test/python')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/variation/src')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

#import matplotlib; matplotlib.use("Agg")	#to avoid popup and collapse in X11-disabled environment

"""
from codense.common import db_connect, form_schema_tables
hostname='dl324b-1'
dbname='yhdb'
schema = 'dbsnp'
hostname='zhoudb'
dbname='graphdb'
schema = 'dbsnp'
#conn, curs = db_connect(hostname, dbname, schema)

hostname='localhost'
dbname='stock20070829'
import MySQLdb
conn0 = MySQLdb.connect(db=dbname,host=hostname)
curs0 = conn0.cursor()

hostname='papaya.usc.edu'
dbname='stock_250k'
db_user='yh'
db_passwd = ''
import MySQLdb
conn = MySQLdb.connect(db=dbname,host=hostname, user=db_user, passwd=db_passwd)
curs = conn.cursor()

drivername='mysql'
schema = None
import Stock_250kDB
db_250k = Stock_250kDB.Stock_250kDB(drivername=drivername, username=db_user,
				password=db_passwd, hostname=hostname, database=dbname, schema=schema)
db_250k.setup(create_tables=False)

drivername='mysql'
hostname='papaya.usc.edu'
dbname='stock'
db_user='yh'
db_passwd = ''
schema = None
from StockDB import StockDB
db_149 = StockDB(drivername=drivername, username=db_user,
				password=db_passwd, hostname=hostname, database=dbname, schema=schema)
db_149.setup(create_tables=False)

"""

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = Main
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
