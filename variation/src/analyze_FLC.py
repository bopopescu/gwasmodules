"""
Analyze FLC locus 
"""

import sequences,dataParsers,snpsdata
import pdb

FLC_dict = {'202D': 6046, '202G': 6900, '309G': 8258, '210C': 7183, '206B': 6969, 
	'206C': 6826, '202H': 6901, 'Ull-2-5': 6974, '210H': 6945, '311G': 6064, 
	'311E': 8237, '303G': 7182, '208D': 6980, '311D': 8240, '309A': 8242, 
	'cut_3184000_31700': 6909,'cut_3184000_3170000': 6909, 'Lov1': 6043, '211H': 6929, '203F': None, 
	'203C': 6964, '308E': 8247, '308D': 8306, '308F': 8266, 'Var-2-6': 7517, 
	'211E': 6971, '310D': None, 'Col': 6909, '312E': 8284, '209C': 6897, 
	'306D': 7285, '204E': 6956}


def load_FLC_phenotypes():
	import csv
	import phenotypeData as pd
	"""
	Load the new 01/12/10 phenotype file
	"""
	filename = "/Users/bjarnivilhjalmsson/Projects/FLC_analysis/phenotype_data_011710.csv"
	print "Loading phenotype/accession file:",filename
	f = open(filename,"r")
	reader = csv.reader(f)
	phenotype_names = reader.next()[3:]
	#print phenotype_names
	accession_names = []
	accession_ID = []
	flc_id_to_ecotype_map = {}
	phenotypes=[[] for i in range(len(phenotype_names))]	#[phenotype_name][acc_id]
	for row in reader:
		accession_names.append(row[2].lower())
		accession_ID.append(row[0])
		for i, phen_val in enumerate(row[3:]):
			try:
				p_val = float(phen_val)
			except Exception:
				p_val = "NA"
			#print p_val
			phenotypes[i].append(p_val)
	f.close()
	#print accession_names
	acc_dict = pd._getAccessionToEcotypeIdDict_(accession_names)
	#acc_dict["cibc-5"] = 6908
	#acc_dict["pla-0"] = 8357
	#print acc_dict
	#accession_names.sort()
	new_phenotypes=[[] for i in range(len(phenotype_names))]
	ecotypes = []
	for acc in acc_dict:
		acc_i = accession_names.index(acc)
		ecotype = acc_dict[acc]
		flc_id_to_ecotype_map[accession_ID[acc_i]] = ecotype
		ecotypes.append(ecotype)
		for i in range(len(phenotype_names)):
			new_phenotypes[i].append(phenotypes[i][acc_i])
	#print new_phenotypes
	#print len(ecotypes)
	#return {"phenotypes":new_phenotypes,"phenotype_names":phenotype_names, "ecotypes":ecotypes}
	phenotypes = map(list,zip(*new_phenotypes))
	phend = pd.PhenotypeData(ecotypes, phenotype_names, phenotypes)
	phend.writeToFile("/tmp/FLC_phenotypes_011710.tsv", delimiter = "\t")
	flc_id_to_ecotype_map['Col']=6909
	flc_id_to_ecotype_map['Lov1']=6043
	flc_id_to_ecotype_map['Ull-2-5']=6974
	flc_id_to_ecotype_map['Var-2-6-Part']=7517
	print flc_id_to_ecotype_map
	return phend, flc_id_to_ecotype_map

def load_raw_flc_data():
	"""
	Load the new 01/12/10 sequences files
	"""
	import sequences
	phend, e_map = load_FLC_phenotypes()
	print e_map
	file_name_1 = "/Users/bjarnivilhjalmsson/Projects/FLC_analysis/FLC_new_seqs_Feb_2010.fas"
	sd = sequences.readFastaFile(file_name_1,e_map,split_name_by='_')
#	file_name_1 = "/Users/bjarnivilhjalmsson/Projects/FLC_analysis/FLC_first_batch.txt"
#	file_name_2 = "/Users/bjarnivilhjalmsson/Projects/FLC_analysis/FLC_second_batch.txt"
#	sd = sequences.readFastaFile(file_name_1,e_map)
#	sd.add(sequences.readFastaFile(file_name_2,e_map))
	sd.reverse_sequences()
	col_seq = sequences.readFastaSequence("/Users/bjarnivilhjalmsson/Projects/Data/Col_sequence_TAIR_071009/chr5_tair8.fas", 3170500, 3183000, "ref_col-0", 6909)
	sd.add_sequence(col_seq)
	return sd,phend

def _save_raw_flc_as_file(filename="/Users/bjarnivilhjalmsson/Projects/FLC_analysis/flc_seqs_040510.fasta"):
	sd,phend = load_raw_flc_data()
	sd.save_fasta_file(filename,data_id='raw')

def _process_raw_flc_(filename="/Users/bjarnivilhjalmsson/Projects/FLC_analysis/flc_seqs_aln_050310.fasta"):
	ref_seq_name = "raw_ref_col-0"
	ref_start = 3170501
	ref_chr = 5
	ad = sequences.readFastaAlignment(filename,
					ref_seq_name=ref_seq_name,ref_start=ref_start,ref_chr=ref_chr,
					alignment_type="muscle",ref_direction=1)
	ad.merge_ecotypes()
	#ad.save_fasta_file("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/flc_seqs_aln_merged_011810.fasta",data_id='')
	ad.save_fasta_file("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/flc_seqs_aln_merged_050410.fasta",data_id='')



def _impute_FLC_192_():
	phed =  pd.readPhenotypeFile("/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/FLC_phenotypes_011710.tsv") 

	d250k_file = env.home_dir+"Projects/Data/250k/250K_192_043009.csv"
	d250k_sd = dataParsers.parse_snp_data(d250k_file)
	d250k_sd.filter_accessions(phed.accessions)
	d250k_sd.filter_maf_snps(0.05)
	
	seq_snpsd = dataParsers.parseCSVData(data_dir+"/flc_seqs_aln_imputed_snps_012710.csv")
	seq_snpsd.onlyBinarySnps()
	
	d250k_sd.snpsDataList[4].compareWith(seq_snpsd)
	d250k_sd.snpsDataList[4].merge_data(seq_snpsd)
	

def _generate_250K_2010_FLC_data_(impute=True):
	"""
	Create a combined version of 
	250K, overlapping with the FLC phenotypes.
	Then merge with 2010 data (including indels).
	Then merge with FLC sequences.
	Impute missing SNPs.
	write to file.
	"""
	import phenotypeData as pd
	import env
	
	phed =  pd.readPhenotypeFile("/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/FLC_phenotypes_011710.tsv") 
	
	d2010_file = env.home_dir+"Projects/Data/2010/2010_imputed_012610.csv"
	d2010_sd = dataParsers.parse_snp_data(d2010_file,id="2010_data")
	d2010_sd.filter_accessions(phed.accessions)
	d2010_sd.filter_na_snps()
	d2010_sd.filter_maf_snps(0.05)

	#d250k_file = env.home_dir+"Projects/Data/250k/250K_t54.csv"
	d250k_file = env.home_dir+"Projects/Data/250k/250K_192_043009.csv"
	d250k_sd = dataParsers.parse_snp_data(d250k_file)
	d250k_sd.filter_accessions(phed.accessions)
	d250k_sd.filter_maf_snps(0.05)
	
	d250k_sd.merge_snps_data(d2010_sd)
	d250k_sd.filter_na_accessions()
	d250k_sd.filter_na_snps(0.7)	
	d250k_sd.filter_monomorphic_snps()
	

	ref_seq_name = "raw_ref_col-0"
	ref_start = 3170501
	ref_chr = 5
	seq_file = env.home_dir+"Projects/FLC_analysis/flc_seqs_aln_merged_050410.fasta"
	ad = sequences.readFastaAlignment(seq_file,ref_seq_name=ref_seq_name,ref_start=ref_start,
			ref_chr=ref_chr,alignment_type="muscle",ref_direction=1)
#	ref_start = 3170500
#	ad2 = sequences.readFastaAlignment(seq_file,ref_seq_name=ref_seq_name,ref_start=ref_start,
#			ref_chr=ref_chr,alignment_type="muscle",ref_direction=1)
#	ref_start = 3170502
#	ad3 = sequences.readFastaAlignment(seq_file,ref_seq_name=ref_seq_name,ref_start=ref_start,
#			ref_chr=ref_chr,alignment_type="muscle",ref_direction=1)
	pdb.set_trace()
	r = ad.get_snps(type=0)
	seq_snpsd1 = r['snpsd']
	seq_snpsd1.merge_data(r['indels'],error_threshold=0.0)

#	r2 = ad2.get_snps(type=0)
#	seq_snpsd2 = r2['snpsd']
#	seq_snpsd2.merge_data(r2['indels'],error_threshold=0.0)
#
#	r3 = ad3.get_snps(type=0)
#	seq_snpsd3 = r3['snpsd']
#	seq_snpsd3.merge_data(r3['indels'],error_threshold=0.0)
	
	
	print "Now merging data.."

	d250k_sd.snpsDataList[4].compareWith(seq_snpsd1)
#	d250k_sd.snpsDataList[4].compareWith(seq_snpsd2)
#	d250k_sd.snpsDataList[4].compareWith(seq_snpsd3)
	d250k_sd.snpsDataList[4].merge_data(seq_snpsd1,union_accessions=False)
	d250k_sd.filter_na_accessions()
	d250k_sd.filter_na_snps(0.7)	
	d250k_sd.filter_monomorphic_snps()	
	d250k_sd.snpsDataList[4].impute_data()
	d250k_sd.writeToFile("/tmp/test.csv")
	print "YEAH!"
	
	
	
	
def _insert_merged_data_in_db_():
	"""
	Ad hoc, to fix a problem..
	"""
	sd_t54 = dataParsers.parse_snp_data('/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t54.csv',filter=0.001)
	d = {}
	for eid,aid in zip(sd_t54.accessions,sd_t54.array_ids):
		d[eid]=aid
	sd_t57 = dataParsers.parse_snp_data('/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_merged_2010_250K_w_FLC_seq.csv')
	aids = [d[eid] for eid in sd_t57.accessions]
#	for sd in sd_t57.snpsDataList: 
#		sd.arrayIds=aids
	sd_t57.arrayIds=aids
	sd_t57.write_to_file_yu_format('/Users/bjarnivilhjalmsson/Projects/Data/250k/call_method_57.tsv')
	#sd_t57._generate_db_call_files_(call_method=57,array_ids=aids)
	

def _insert_markers_into_db_():
	import bisect,dbutils
	sd_192 = dataParsers.parse_snp_data('/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_192_043009.csv')
	sd_t57 = dataParsers.parse_snp_data('/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_merged_2010_250K_w_FLC_seq.csv')
	cpl_192 = sd_192.getChrPosList()
	cpsl_t57 = sd_t57.getChrPosSNPList()
	found_count = 0
        conn = dbutils.connect_to_papaya()
        cursor = conn.cursor()                            
	for c,p,snp in cpsl_t57:
		if c<5:
			continue
		i = bisect.bisect(cpl_192,(c,p))
		if not cpl_192[i-1]==(c,p):
  	                #Check if the SNP is in the DB.      
  	                alleles = list(set(snp))
	                sql_statement = "SELECT id FROM stock_250k.snps WHERE chromosome=%d AND position=%d AND \
	                		(allele1='%s' OR allele1='%s');"%(c,p,alleles[0],alleles[1])
	                print sql_statement
                	found = False
                        num_rows = int(cursor.execute(sql_statement))
                        row = cursor.fetchone()
                        if row:
                        	print row
                        	found=True
                        if not found:
                        	#Insert SNP in DB.
                        	snp_name = str(c)+"_"+str(p)+"_"+alleles[0]+"_"+alleles[1]
		                sql_statement = "INSERT INTO stock_250k.snps (name, chromosome, position, allele1, allele2)\
		                		 VALUES ('%s',%d,%d,'%s','%s');"%(snp_name,c,p,alleles[0],alleles[1])
		                print sql_statement
		                try:
		                	cursor.execute(sql_statement)
                                        print "Committing transaction (making changes permanent)."
                                        conn.commit()  			
                                except:
                                	print "insert failed... moving on"
		                	
        #Close connection
        cursor.close()
        conn.close()
 			
		


#def compareFLCand250kSNPs():
#	(FLC_positions,FLC_snps) = getSNPsFromSequences("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/FLC_muscle_072109.aln",ref_start_pos=3170000,ref_seq="cut_3184000_3170000",reversed=True)
#	import dataParsers
#	snpsd = dataParsers.parseCSVData("/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_f13_012609.csv")[4]
#	i = 0
#	snps_cut = []
#	positions_cut = []
#	while snpsd.positions[i]<3170000:
#		i += 1
#	while snpsd.positions[i]<=3184000:
#		snps_cut.append(snpsd.snps[i])
#		positions_cut.append(snpsd.positions[i])
#		i+=1
#	print positions_cut


def plot_flc_haplotype():
	import analyzeHaplotype as ah
	#res = readFastaFile("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/FLC_muscle_072109.aln")
	#seqs = res["sequences"]
	#seq_names = res["names"]
	(positions,aln_snps,seq_names) = getSNPsFromSequences("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/FLC_muscle_072109.aln",
							ref_start_pos=3170001,ref_seq="cut_3184000_3170000",reversed=True)
	#aln_snps = map(list,zip(*seqs))
	#seqs = reverse_sequences(seqs)
	i = seq_names.index("cut_3184000_3170000")
	seq_names[i] = "Col_TAIR8"	
	(flc_250k_snps,flc_snps,flc_250K_positions,accessions,flc_data_acc_map) = get_overlapping_snps_in_region()
	print flc_data_acc_map
	import phenotypeData as pd
	a_dict = pd._getEcotypeIdInfoDict_()
	new_accessions = []
	for acc in accessions:
		new_accessions.append(unicode(a_dict[int(acc)][0],'iso-8859-1'))
	accessions = new_accessions
	ah.plot_flc_haplotypes(aln_snps,positions=positions,accessions=seq_names,
			haplotypeFile="/Users/bjarnivilhjalmsson/tmp/aln_haplotype.pdf",
			treeFile="/Users/bjarnivilhjalmsson/tmp/aln_tree.pdf",acc_250k=flc_data_acc_map,
			flc_250K_positions=flc_250K_positions)
	
	#ah.plot_haplotypes(flc_250k_snps,positions=flc_250K_positions,accessions=accessions,haplotypeFile="/Users/bjarnivilhjalmsson/tmp/205k_FLC_haplotype.pdf",treeFile="/Users/bjarnivilhjalmsson/tmp/250k_FLC_tree.pdf")


#def load_sequences_in_dir(dir="/Users/bjarnivilhjalmsson/Projects/FLC_analysis/data_102509/Sequences"):
#	"""
#	Load fasta sequences in directory.
#	"""
#	import os
#	seq_files = os.listdir(dir)
#	seq_list = []
#	accessions = [] 
#	for s_file in seq_files:
#		#print s_file
#		acc = s_file[0:-4].split()[0]
#		accessions.append(acc.lower())
#		f = open(dir+"/"+s_file,"r")
#		lines = f.readlines()
##		i = 0
##		line = lines[i]
##		if line[0] == "^" :
##			i += 1
##			while line[0]!="^" and i < len(lines)-1:
##				print line[0:-1]
##				i += 1
##				line = lines[i]
#		seq = ""
#		for line in lines:
#			seq += line.strip()
#		
#		seq = seq.upper()
#		s = sequences.Sequence(seq,acc)
#		seq_list.append(s)
#		f.close()
#	#print sequences
#	#print accessions
#	print len(accessions), "sequences loaded from the following directory:",dir
#	return sequences.SequenceData(seq_list)



#def compare_col_seqs():
#	import pdb
#	region_start = 3170001
#	region_end = 3184001
#	#get ref. col seq.
#	ref_seq = cutFastaFile("/Users/bjarnivilhjalmsson/Projects/Data/Col_sequence_TAIR_071009/NCBI_chr5",region_start,region_end)
#	#print ref_seq
#	#get FLC col seq.
#	res = readFastaFile("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/FLC_muscle_072109.aln")
#	seqs = res["sequences"]
#	seqs = reverse_sequences(seqs)
#	seq_names = res["names"]
#	ref_index = seq_names.index("cut_3184000_3170000")
#	col_index = seq_names.index("Col")
#	aln_ref_seq = seqs[ref_index].replace("-","")
#	i = 0
#	while seqs[col_index][i] == "-":
#		i+= 1
#	aln_col_seq = "-"*i 
#	j = 1
#	while seqs[col_index][len(seqs[col_index])-j] == "-":
#		j+= 1
#	while  i < len(seqs[col_index])-j+1:
#		if seqs[col_index][i] != "-":
#			aln_col_seq += seqs[col_index][i]
#		i += 1
#	aln_col_seq += "-"*(j-1)
#	#get perlegen col seq.
#	snpsds = dataParsers.parseCSVData("/Users/bjarnivilhjalmsson/Projects/Data/perlegen/perlegen_011609.csv")
#	print len(snpsds)
#	snpsd = snpsds[4]
#	col_index = snpsd.accessions.index("6909")
#	print col_index
#	seq_perlegen = ""
#	seq_perlegen_positions = []
#	i = 0
#	pos = snpsd.positions[i]
#	while pos < region_start:
#		i += 1
#		pos = snpsd.positions[i]
#	print 'i:',i
#	print 'pos:',pos
#	while pos < region_end:
#		nt = snpsd.snps[i][col_index]
#		if nt =="NA":
#			nt = 'N'
#		seq_perlegen += nt
#		seq_perlegen_positions.append(pos)
#		i += 1
#		pos = snpsd.positions[i]
#	print 'i:',i
#	print 'pos:',pos
#	print "Found",len(seq_perlegen),"SNPs in perlegen data"
#	print seq_perlegen
#	print seq_perlegen_positions
#	col_perlegen_seq = ""
#	last_pos = region_start
#	for (nt,pos) in zip(seq_perlegen,seq_perlegen_positions):
#		print pos
#		col_perlegen_seq+="-"*(pos-last_pos-1)
#		col_perlegen_seq+=nt
#		last_pos = pos
#		
#	col_perlegen_seq+="-"*(region_end-last_pos)
#	pdb.set_trace()
#	
#	
#	
#	#get 250k col seq.
#	snpsds = dataParsers.parseCSVData("/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_f13_012609.csv")
#	print len(snpsds)
#	snpsd = snpsds[4]
#	col_index = snpsd.accessions.index("6909")
#	print col_index
#	seq_250k = ""
#	seq_250k_positions = []
#	i = 0
#	pos = snpsd.positions[i]
#	while pos < region_start:
#		i += 1
#		pos = snpsd.positions[i]
#	while pos < region_end:
#		seq_250k += snpsd.snps[i][col_index]
#		seq_250k_positions.append(pos)
#		i += 1
#		pos = snpsd.positions[i]
#	print "Found",len(seq_250k),"SNPs in 250K data"
#	col_250k_seq = ""
#	last_pos = region_start
#	for (nt,pos) in zip(seq_250k,seq_250k_positions):
#		col_250k_seq+="-"*(pos-last_pos-1)
#		col_250k_seq+=nt
#		last_pos = pos
#	col_250k_seq+="-"*(region_end-last_pos)
#	print "len(aln_col_seq)",len(aln_col_seq)
#	print "len(aln_ref_seq)",len(aln_ref_seq)
#	print "len(ref_seq)",len(ref_seq)
#	print "len(col_250k_seq)",len(col_250k_seq)
#	print "len(col_perlegen_seq)",len(col_250k_seq)
#	print seq_250k_positions
#	#print aln_col_seq[0:2000]
#	#print aln_ref_seq[0:2000]
#	#print ref_seq[0:2000]
#	print col_250k_seq[0:2000]
#	print col_perlegen_seq[0:2000]
#	pdb.set_trace()
#	
	

	
	#get 2010 col seq.
	#overlap_nts = get_overlaping_nts("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/FLC_muscle_072109.aln",6909,3170000,True)
	



def _run_():
	sd = load_sequences_in_dir()
	cons_sd = load_sequences_in_dir(dir="/Users/bjarnivilhjalmsson/Projects/FLC_analysis/data_102509/Consensus_Sequences")
	print sd
	sd.reverse_sequences()
	sd.add(cons_sd)
	
	#Retrieve 2010 sequences for the FLC region.
	chr = 5
	region_start = 3170001
	region_end = 3184001
	sd_list = sequences.get2010sequence_datas(chr,region_start,region_end)
	sd_list[0].reverse_sequences()
	a = sd_list[0]
	print len(a.sequences)
	a.strip()
	print len(a.sequences)

	#for i in range(1,len(sd_list)):
		#sd_list[i].reverse_sequences()
	#	a.extend_alignment_with(sd_list[i])
	sd_2010 = a.get_sequence_data()
	#sd.add(sd_2010)
	
	#Load COL reference sequence.
	col_sd = sequences.readFastaFile("/Users/bjarnivilhjalmsson/Projects/Data/Col_sequence_TAIR_071009/chr5_FLC.fasta")
	#col_sd.add(sd)
	col_sd.add(sd_2010)
	col_sd.generate_unique_names()
	col_sd.save_fasta_file("/tmp/FLC_with_2010_1.fasta")

	#Load 2010 data
#	snps_dataset_2010 = dataParsers.parse_snp_data("/Users/bjarnivilhjalmsson/Projects/Data/2010/2010_073009.csv")
#	snpsd = snps_dataset_2010.get_region_snpsd(chr,region_start,region_end)
#	print "2010 data has",len(snpsd.snps),"SNPs in the FLC region."
#	print snpsd.positions
#	
	
	#Compare with other data
	#	- Previous sequences
	#	- 2010
	#	- Perlegen
	#	- 250K
	#print res["sequences"]
	#save_fasta_file("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/data_102509/sequences_v2.fasta", sequences, accessions)


#def locating_2010_FLC_sequences():
#	col_seq = sequences.readFastaSequence("/Users/bjarnivilhjalmsson/Projects/Data/Col_sequence_TAIR_071009/chr5_tair8.fas",3170000,3184000,name="Col-0",ecotype=6909)
#	len(col_seq.seq)
#	col_sd = sequences.SequenceData([col_seq])
#	col_sd.add_data_id("ref_1")
#	col_seq = sequences.readFastaSequence("/Users/bjarnivilhjalmsson/Projects/Data/Col_sequence_TAIR_071009/chr5_tair8.fas",3170000,3184000,name="Col-0",ecotype=6909)
#	len(col_seq.seq)
#	col_2 = sequences.SequenceData([col_seq])
#	col_2.add_data_id("ref_2")
#	col_sd.add(col_2)
#	
#	#Loading 2010 data
#	a = sequences.get2010sequence_datas(5,3170000,3184000)
#	sd_2010 = a[0]
#	print sd_2010.ecotypes
#	sd_2010.filter_missing_sequences()
#	sd_2010.strip()
#	sd_2010.filter_accessions_calls()
#	for i in range(1,len(a)):
#		a[i].filter_missing_sequences()
#		a[i].strip()
#		a[i].filter_accessions_calls()
#		sd_2010.extend_alignment_with(a[i])
#	sd_2010 = sd_2010.get_sequence_data()
#	print len(sd_2010.sequences)
#	sd_2010.strip()
#	sd_2010.add_data_id("2010")
#	col_sd.add(sd_2010)
#	
#	#Old Fasta sequences.
#	sd_old = sequences.readFastaFile("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/FLC_sequences050209.fasta")
#	new_sequences = []
#	for s in sd_old.sequences:
#		s.ecotype = FLC_dict[s.seq_name] 
#		if s.ecotype:
#			new_sequences.append(s)
#	sd_old = sequences.SequenceData(new_sequences)
#	sd_old.reverse_sequences()
#	sd_old.add_data_id("old")
#	col_sd.add(sd_old)
#
##	#Short sequences
##	sd = load_sequences_in_dir()
##	sd._guess_ecotypes_()
##	sd.add_data_id("short")
##	col_sd.add(sd)
##	
#	cons_sd = load_sequences_in_dir(dir="/Users/bjarnivilhjalmsson/Projects/FLC_analysis/data_102509/Consensus_Sequences")
#	cons_sd.reverse_sequences()
#	cons_sd._guess_ecotypes_()
#	cons_sd.add_data_id("cons")
#	print cons_sd.ecotypes
#	col_sd.add(cons_sd)
#
#	col_sd.filter_missing_sequences()
#	col_sd.strip()
#	#col_sd.generate_unique_names()	
#	print col_sd.ecotypes
#	data_dir = "/Users/bjarnivilhjalmsson/Projects/FLC_analysis/data_102509/"
#	file_name = data_dir+"FLC_seqs2_120209.fasta"
#	col_sd.filter_missing_ecotypes()
#	print col_sd.ecotypes, len(col_sd.ecotypes),len(set(col_sd.ecotypes))
#	col_sd.save_fasta_file(file_name)
	

#def cut_and_merge_alignments():
#	data_dir = "/Users/bjarnivilhjalmsson/Projects/FLC_analysis/data_102509/"
#	ref_seq_name = "ref_1_Col-0"
#	ref_start = 3170001
#	ref_chr = 5
#	ad_2010 = sequences.readFastaAlignment(data_dir+"FLC_seqs2_120209.aln.fasta",ref_seq_name=ref_seq_name,ref_start=ref_start,
#				ref_chr=ref_chr,alignment_type="muscle",ref_direction=1)
#	
#	d250k_file = "/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t43_192.csv"
#	d250k_sd = dataParsers.parse_snp_data_region(d250k_file,5,3170001,3184000)
#	#ad_2010.compare_with_snps_data(d250k_sd)
#	
#	ad_short = sequences.readFastaAlignment(data_dir+"FLC_seqs1_120209.aln.fasta",ref_seq_name=ref_seq_name,ref_start=ref_start,
#				ref_chr=ref_chr,alignment_type="muscle",ref_direction=1)
#	#ad_short.compare_with_snps_data(d250k_sd)
#	
#	ad_2010.cut_alignment(3170000,3178000)
#	#ad_2010.save_fasta_file("/tmp/test.fasta")
#	ad_short.cut_alignment(3178001,3184000)
#	#ad_short.save_fasta_file("/tmp/test.fasta")
#	ad_2010.extend_alignment_with(ad_short,type=2)
#	ad_2010.save_fasta_file(data_dir+"FLC_full.aln.fasta")
#	import env
#	ad_2010.save_fasta_file(env.home_dir+"tmp/test.fasta")
#	



def get_alignment():
	data_dir = "/Users/bjarnivilhjalmsson/Projects/FLC_analysis/data_102509/"
	ref_seq_name = "ref_2_Col-0"
	ref_start = 3170001
	ref_chr = 5
	ad_2010 = sequences.readFastaAlignment(data_dir+"FLC_full.aln.fasta",ref_seq_name=ref_seq_name,ref_start=ref_start,
				ref_chr=ref_chr,alignment_type="muscle",ref_direction=1)
	d250k_file = "/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t43_192.csv"
	d250k_sd = dataParsers.parse_snp_data_region(d250k_file,5,3170001,3184000)
	ad_2010.compare_with_snps_data(d250k_sd)
	
	#m = ad_2010.get_similarity_matrix()
	#print m
	ad_2010.merge_ecotypes()
	ad_2010.save_fasta_file(data_dir+"FLC_full_merged.aln.fasta")
	ad_2010.compare_with_snps_data(d250k_sd)
	ad_2010 = sequences.readFastaAlignment(data_dir+"FLC_full_merged.aln.fasta",ref_seq_name=ref_seq_name,ref_start=ref_start,
				ref_chr=ref_chr,alignment_type="muscle",ref_direction=1)
	ad_2010.compare_with_snps_data(d250k_sd)



def analyzeSNPs():
	import KW,phenotype_parsers,phenotypeData
	import Emma
	result_id = "filtered_imputed"
	data_dir = "/Users/bjarnivilhjalmsson/Projects/FLC_analysis/"
	#ref_seq_name = "2010_Col-0"
	ref_seq_name = "raw_ref_col-0"
	ref_start = 3170501
	ref_chr = 5
	#ad_2010 = sequences.readFastaAlignment(data_dir+"FLC_full_edited_merged.aln.fasta",ref_seq_name=ref_seq_name,ref_start=ref_start,
	#			ref_chr=ref_chr,alignment_type="muscle",ref_direction=1)
	#ad_2010 = sequences.readFastaAlignment(data_dir+"FLC_full_merged.aln.fasta",ref_seq_name=ref_seq_name,ref_start=ref_start,
	#			ref_chr=ref_chr,alignment_type="muscle",ref_direction=1)
	#ad = sequences.readFastaAlignment(data_dir+"flc_seqs_aln_merged_011810.fasta",ref_seq_name=ref_seq_name,ref_start=ref_start,
	#		ref_chr=ref_chr,alignment_type="muscle",ref_direction=1)

	#r = ad.get_snps(type=1)
	#seq_snpsd = r['snpsd']
	#seq_snpsd = seq_snpsd.getSnpsData(missingVal='NA')
	#seq_snpsd.onlyBinarySnps()
	#i_snpsd = r['indels']
	#print indels
	#i_snpsd = i_snpsd.getSnpsData(missingVal='NA')
	#print zip(i_snpsd.positions, i_snpsd.snps)
	#print i_snpsd.accessionsl
	seq_snpsd = dataParsers.parseCSVData(data_dir+"/flc_seqs_aln_imputed_snps_012510.csv")[0]	
	seq_snpsd = seq_snpsd.getSnpsData(missingVal='NA')

#	d2010_file = "/Users/bjarnivilhjalmsson/Projects/Data/2010/2010_073009.csv"
	d2010_file = "/Users/bjarnivilhjalmsson/Projects/Data/2010/2010_imputed_012610.csv"
	d2010_sd = dataParsers.parse_snp_data(d2010_file,id="2010_data")
#	d2010_sd.filter_na_accessions()
	d2010_sd.filter_na_snps()
	d2010_sd.convert_2_binary()
	d2010_sd.filter_maf_snps(0.05)
	#kinship_2010 = Emma.calcKinship(d2010_sd.getSnps(0.05))
	d2010_sd = d2010_sd.get_region_snpsd(5,3140000,3220000)
	d2010_sd.remove_redundant_snps(w_missing=True)

	d250k_file = "/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_data_t43_081009.csv"
	snpsd = dataParsers.parse_snp_data(d250k_file)
	snpsd.filter_accessions(seq_snpsd.accessions)
	snpsd.convert_2_binary()
	snpsd.filter_maf_snps(0.05)
	#kinship_250k = Emma.calcKinship(snpsd.getSnps(0.02))
	
	snpsd = snpsd.get_region_snpsd(5,3140000,3220000)
	snpsd.remove_redundant_snps()
	
	seq_snpsd.remove_accessions(snpsd.accessions)
	seq_snpsd.snpsFilterRare(0.05)
	seq_snpsd.onlyBinarySnps()
	acc_map = []
	for i, acc in enumerate(seq_snpsd.accessions):
		acc_map.append((i,snpsd.accessions.index(acc)))
		
	seq_snpsd.orderAccessions(acc_map)
	seq_snpsd.remove_redundant_snps(w_missing=True)
	
	#snpsd.mergeDataUnion(d2010_sd,priority=2,unionType=3)
	#ad.compare_with_snps_data(snpsd) #Something missing here snpsd...?
	#i_snpsd = 
	#snpsd.mergeDataUnion(d250k_sd,unionType=3,verbose=True)


	#NOW PERFORM GWAS AND PLOT RESULT!!!!
	
	phend =  phenotypeData.readPhenotypeFile("/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/FLC_phenotypes_011710.tsv") 
	#phenotype_parsers.load_phentoype_file("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/data_102509/FLC_soil_data_102509.csv")
	results_colors = ['blue','green','red']
	#kinship_matrices = [kinship_250k,kinship_250k,kinship_2010]
	snpsds = [snpsd,seq_snpsd,d2010_sd]
	phenotypeIndices = phend.phenIds
	log_transforms = [1,2]
	import analyzePhenotype as ap
	import analyzeSNPResult as asr
	import copy
	
#	for i in phenotypeIndices:
#		#ap.drawHistogram(phend,i,pdfFile="/Users/bjarnivilhjalmsson/tmp/hist_"+str(phend.getPhenotypeName(i))+".pdf")
#		#if i in log_transforms:
#		phend.logTransform(i)
#		#print "log transforming"
#		results = []
#		filtered_sds=[]	
#		for sd,k in zip(snpsds,kinship_matrices):
#			new_sd = copy.deepcopy(sd)
#			res = Emma.run_emma_w_missing_data(new_sd,phend,i,5,k)
#			res.negLogTransform()
#			snps_indices_to_keep = res.filterMARF(minMaf=0.1)
#			print "Got",len(res.scores),len(res.positions),"p-values from Emma."
#			results.append(res)
#			#pvals = res.scores
#			#positions = res.positions
#			#pp = zip(pvals,positions)
#			#pp.sort()
#			#print pp 
#			#import plotResults as pr
#			#pr.plotResult(res,"/Users/bjarnivilhjalmsson/tmp/test.pdf")
#			new_sd.filter_snp_indices(snps_indices_to_keep)
#			filtered_sds.append(new_sd)
#		import regionPlotter as rp
#		reg_plotter = rp.RegionPlotter()
#		reg_plotter.plot_small_result(results,results_colors=results_colors,
#					pdf_file="/Users/bjarnivilhjalmsson/tmp/seqences_250k_"+result_id+"_emma_gwas_"+str(phend.getPhenotypeName(i))+".pdf")
#		for j,(r,sd) in enumerate(zip(results,filtered_sds)):
#			r_i = r.scores.index(max(r.scores))
#			phend.plot_marker_box_plot(i,sd,r_i,pdf_file="/Users/bjarnivilhjalmsson/tmp/box_plot_emma_"+str(phend.getPhenotypeName(i))+"_"+results_colors[j]+".pdf",marker_score=r.scores[r_i])
#			
	phend =  phenotypeData.readPhenotypeFile("/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/FLC_phenotypes_011710.tsv") #phenotype_parsers.load_phentoype_file("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/data_102509/FLC_soil_data_102509.csv")
	
	for i in phenotypeIndices:
		results = []	
		filtered_sds = []	
		for sd in snpsds:
			new_sd = copy.deepcopy(sd)
			res, f_sd= KW.run_kw(new_sd, phend, i,5)
			filtered_sds.append(f_sd)
			res.negLogTransform()
			print "Got",len(res.scores),len(res.positions),"p-values from KW."
			results.append(res)
			#pvals = res.scores
			#positions = res.positions
			#pp = zip(pvals,positions)
			#pp.sort()
			#print pp 
			#import plotResults as pr
			#pr.plotResult(res,"/Users/bjarnivilhjalmsson/tmp/test.pdf")
		import regionPlotter as rp
		reg_plotter = rp.RegionPlotter()
		reg_plotter.plot_small_result(results,results_colors=results_colors,
					pdf_file="/Users/bjarnivilhjalmsson/tmp/seqences_250k_"+result_id+"_gwas_"+str(phend.getPhenotypeName(i))+".pdf")
		for j,(r,sd) in enumerate(zip(results,filtered_sds)):
			if len(r.scores)!=len(sd.snps):
				print "Lengths not equal? %d, %d", (len(r.scores),len(sd.snps))
			r_i = r.scores.index(max(r.scores))
			phend.plot_marker_box_plot(i,sd,r_i,pdf_file="/Users/bjarnivilhjalmsson/tmp/box_plot_kw_"+str(phend.getPhenotypeName(i))+"_"+results_colors[j]+".pdf",marker_score=r.scores[r_i])
			
			



def compareSNPs():
	"""
	Compare SNPs with Perlegen and 250K data...
	"""
	data_dir = "/Users/bjarnivilhjalmsson/Projects/FLC_analysis/data_102509/"
	#ref_seq_name = "2010_Col-0"
	ref_seq_name = "ref_2_Col-0"
	ref_start = 3170000
	ref_chr = 5
	#ad = sequences.readFastaAlignment(data_dir+"FLC_full_edited_merged.aln.fasta",ref_seq_name=ref_seq_name,ref_start=ref_start-1,
	ad = sequences.readFastaAlignment(data_dir+"FLC_full_merged.aln.fasta",ref_seq_name=ref_seq_name,ref_start=ref_start-1,
				ref_chr=ref_chr,alignment_type="muscle",ref_direction=1)
	
	#Load perlegen data in region...
	chromosome = 5
	start_pos = 3170000
	end_pos = 3184000
	perlegen_file = "/Users/bjarnivilhjalmsson/Projects/Data/perlegen/perlegen_073009.csv"
	perlegen_sd = dataParsers.parse_snp_data_region(perlegen_file,chromosome,start_pos,end_pos)
	#m = perlegen_sd.get_mafs()
	for i in range(len(perlegen_sd.positions)-1):
		print perlegen_sd.positions[i],perlegen_sd.positions[i+1]-perlegen_sd.positions[i],perlegen_sd.snps[i].count('NA')
	print perlegen_sd.positions[-1] 
	
	perlg_col_i = perlegen_sd.accessions.index("6909")
	#seq_col_i = ad.seq_names.index("2010_Col-0")
	seq_col_i = ad.seq_names.index("ref_2_Col-0")
	for i in range(len(perlegen_sd.snps)):
		pos = perlegen_sd.positions[i]
		print "Position:",pos
		(alleles,a_i) = ad.get_allele_at_ref_pos(pos,window=[5,5])
		print "Index i:",a_i
		print "Perlg. Col allele:",perlegen_sd.snps[i][perlg_col_i]
		print "AD Col allele:",alleles[seq_col_i]

	#ad.reverse_sequences()
	r = ad.get_snps(type=1)
	snpsd = r['snpsd']  #Includes deletions as well...
	snpsd = snpsd.getSnpsData(missingVal='NA')
	snpsd.accessions = map(str,snpsd.accessions)
	#m = snpsd.get_mafs()
	#for i in range(len(snpsd.positions)-1):
	#	print snpsd.positions[i],snpsd.positions[i+1]-snpsd.positions[i], snpsd.snps[i].count('NA')
	#print snpsd.positions[-1] 

	perlegen_sd.compareWith(snpsd)
	
	

	#d250k_file = "/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t43_192.csv"
	#d250k_sd = dataParsers.parse_snp_data_region(d250k_file,chromosome,start_pos,end_pos)
	#d250k_sd.compareWith(snpsd)


def _read_tree_accession_file_():
	f = open("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/tree_accessions.csv")
	#import codecs
	#f = codecs.open('unicode.rst', encoding='utf-8')
	f.readline()
	lines = f.readlines()
	accessions = []
	for l in lines:
		line = l.split(',')
		accessions.append(line[2].lower())
	return accessions


def plot_local_tree():
	data_dir = "/Users/bjarnivilhjalmsson/Projects/FLC_analysis/"
	accs_to_keep=_read_tree_accession_file_()
	ref_seq_name = "raw_ref_col-0"
	ref_start = 3170501
	ref_end = 3183000
	ref_chr = 5
	intron_start = 3175600 
	intron_stop =  3179100
	#ad_2010 = sequences.readFastaAlignment(data_dir+"FLC_full_edited_merged.aln.fasta",ref_seq_name=ref_seq_name,ref_start=ref_start,
	#			ref_chr=ref_chr,alignment_type="muscle",ref_direction=1)
	#ad_2010 = sequences.readFastaAlignment(data_dir+"FLC_full_merged.aln.fasta",ref_seq_name=ref_seq_name,ref_start=ref_start,
	#			ref_chr=ref_chr,alignment_type="muscle",ref_direction=1)
	ad = sequences.readFastaAlignment(data_dir+"flc_seqs_aln_merged_011810.fasta",ref_seq_name=ref_seq_name,ref_start=ref_start,
			ref_chr=ref_chr,alignment_type="muscle",ref_direction=1)
	#ref_seq_name = "ref_2_Col-0"
	#ref_start = 3170001
	#ref_end = 3184000
	#ref_chr = 5
	#ad_2010 = sequences.readFastaAlignment(data_dir+"FLC_full_merged.aln.fasta",ref_seq_name=ref_seq_name,ref_start=ref_start,
	#			ref_chr=ref_chr,alignment_type="muscle",ref_direction=1)
	
	#r = ad_2010.get_snps(type=1,min_called_fraction=0.1)
	#seq_sd = r['snpsd'] #Raw SNPSs data
	#seq_sd.id = "Sequences"
	#seq_sd.remove_accessions(accs_to_keep,True)
	#seq_sd.filterMonoMorphicSnps()
	#print seq_sd.snps
	#snpsd = seq_sd.getSnpsData(missingVal='NA')
	r = ad.get_snps(type=1)
	seq_snpsd = r['snpsd']
	seq_snpsd.remove_accessions(accs_to_keep,True)
	seq_snpsd.filterMonoMorphicSnps()
	print len(seq_snpsd.snps)
	#i_snpsd = r['indels']
	
	#TREE and HAPLOTYPES
	import analyzeHaplotype as ah
	start_stop_list = [(3170500,3183000),(3172000,3181000),(3172000,3175000),(3175000,3178000),
			(3178000,3181000),(3176000,3181000),(intron_start,intron_stop)]
	for start,stop in start_stop_list:
		snpsd = seq_snpsd.get_region_snpsd(start,stop)
		tree_file="/Users/bjarnivilhjalmsson/tmp/aln_tree_"+str(start)+"_"+str(stop)+".pdf"
		ah.plot_tree(snpsd,tree_file,verbose=False)
	
	#250K
	d250k_file = "/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t43_192.csv"
	d250k = dataParsers.parse_snp_data(d250k_file)
	temp_d250k = snpsdata.RawSnpsData(snps=d250k.getSnps(0.05),accessions=d250k.accessions)
	tree_file="/Users/bjarnivilhjalmsson/tmp/250k_full_data_tree.pdf"
	ah.plot_tree(temp_d250k,tree_file,verbose=True)
	
	d250k_sd = d250k.get_region_snpsd(5,3140000,3220000)
	#d250k_sd = dataParsers.parse_snp_data_region(d250k_file,ref_chr,3140000,3220000,id="250K_data")
	start_stop_list = [(3140000,3220000),(3150000,3210000),(3170501,3183000),(3172000,3181000),
			(3172000,3175000),(3175000,3178000),(3178000,3181000),(3176000,3181000),(intron_start,intron_stop)]
	for start,stop in start_stop_list:
		snpsd = d250k_sd.get_region_snpsd(start,stop)
		tree_file="/Users/bjarnivilhjalmsson/tmp/250k_tree_"+str(start)+"_"+str(stop)+".pdf"
		ah.plot_tree(snpsd,tree_file,verbose=False)
	seq_snpsd.mergeDataUnion(d250k_sd,unionType=1,verbose=True)


	#2010 
	d2010_file = "/Users/bjarnivilhjalmsson/Projects/Data/2010/2010_073009.csv"
	d2010_sd = dataParsers.parse_snp_data_region(d2010_file,ref_chr,3140000,3220000,id="2010_data")
	d2010_sd.filterMissingSnps(50)
	d2010_sd._convert_to_tg_ecotypes_()
	d2010_sd.mergeDataUnion(d250k_sd,unionType=1,verbose=True)
	
	d250k_sd.remove_accessions(accs_to_keep,True)
	for start,stop in start_stop_list:
		snpsd = d250k_sd.get_region_snpsd(start,stop)
		tree_file="/Users/bjarnivilhjalmsson/tmp/250k_filtered_tree_"+str(start)+"_"+str(stop)+".pdf"
		ah.plot_tree(snpsd,tree_file,verbose=False)
	d250k = dataParsers.parse_snp_data(d250k_file)
	d250k.filter_accessions(accs_to_keep,True)
	d250k.filter_monomorphic_snps()
	snps=d250k.getSnps(0.05)
	temp_d250k = snpsdata.RawSnpsData(snps=snps,accessions=d250k.accessions)
	tree_file="/Users/bjarnivilhjalmsson/tmp/250k_full_data_filtered_tree.pdf"
	ah.plot_tree(temp_d250k,tree_file,verbose=False)
	
	
	#Perlegen
	perlegen_file = "/Users/bjarnivilhjalmsson/Projects/Data/perlegen/perlegen_073009.csv"
	perlegen_sd = dataParsers.parse_snp_data_region(perlegen_file,ref_chr,3140000,3220000,id="perlegen_data")
	perlegen_sd._convert_to_tg_ecotypes_()
	perlegen_sd.filterMissingSnps(10)
	d2010_sd.mergeDataUnion(perlegen_sd,priority=2,unionType=1,verbose=True)
	seq_snpsd.mergeDataUnion(perlegen_sd,priority=2,unionType=1,verbose=True)
	
	#250K, 2010, Perlegen TREE
	d2010_sd.filter_accessions_by_NAs(0.9)
	d2010_sd.filterMissingSnps(180)
	d2010_sd.filterMonoMorphicSnps()
	for start,stop in start_stop_list:
		snpsd = d250k_sd.get_region_snpsd(start,stop)
		tree_file="/Users/bjarnivilhjalmsson/tmp/250k_2010_perlegen_tree_"+str(start)+"_"+str(stop)+".pdf"
		ah.plot_tree(snpsd,tree_file,verbose=False)
	
	#250K, 2010, Sequences, Perlegen TREE
	seq_snpsd.filterMonoMorphicSnps()
	seq_snpsd.filter_accessions_by_NAs(0.9)
	seq_snpsd.filterMissingSnps(180)
	start_stop_list = [(3170500,3183000),(3172000,3181000),(3172000,3175000),(3175000,3178000),
			(3178000,3181000),(3176000,3181000),(intron_start,intron_stop)]
	for start,stop in start_stop_list:
		snpsd = seq_snpsd.get_region_snpsd(start,stop)
		tree_file="/Users/bjarnivilhjalmsson/tmp/Seq_250k_2010_perlegen_tree_"+str(start)+"_"+str(stop)+".pdf"
		ah.plot_tree(snpsd,tree_file,verbose=False)


def _plot_local_FLC_haplotype_():
	data_dir = "/Users/bjarnivilhjalmsson/Projects/FLC_analysis/"
	ref_seq_name = "raw_ref_col-0"
	ref_start = 3170501
	ref_chr = 5
	ad = sequences.readFastaAlignment(data_dir+"flc_seqs_aln_merged_011810.fasta",ref_seq_name=ref_seq_name,ref_start=ref_start,
			ref_chr=ref_chr,alignment_type="muscle",ref_direction=1)

	r = ad.get_snps(type=1)
	seq_snpsd = r['snpsd']
	seq_snpsd.filter_na_snps(max_na_rate=0.05)
#	seq_snpsd = dataParsers.parseCSVData(data_dir+"flc_seqs_aln_imputed_snps_012710.csv")[0]	
#	seq_snpsd = seq_snpsd.getSnpsData(missingVal='NA')
#	seq_snpsd.remove_accessions(map(str,ad.ecotypes))
	import phenotypeData as pd
	phend =  pd.readPhenotypeFile("/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/FLC_phenotypes_011710.tsv") #phenotype_parsers.load_phentoype_file("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/data_102509/FLC_soil_data_102509.csv")

	for i in phend.phenIds:
		phen_name = phend.getPhenotypeName(i)
		import analyzeHaplotype as ah
		ah.plot_haplotypes(seq_snpsd.snps,seq_snpsd.accessions,haplotypeFile="/Users/bjarnivilhjalmsson/tmp/flc_seq_haplotypes_old.pdf")
		for start,stop in [(3175500,3176000),(3176000,3176500),(3176500,3177000),(3177000,3177500),(3177500,3178000),(3178000,3178500),(3178500,3179000)]:
			ah.plot_local_haplotypes("/Users/bjarnivilhjalmsson/tmp/flc_seq_haplotypes_"+phen_name+"_"+str(start)+"_"+str(stop)+".pdf",
						seq_snpsd,start,stop,phenotypeData=phend,phen_id=i)






if __name__ == "__main__":
	#_run_()
	#locating_2010_FLC_sequences()
	#cut_and_merge_alignments()
	#get_alignment()
	#load_FLC_phenotypes()
	#plot_local_tree()
	#_save_raw_flc_as_file()
	#_impute_FLC_192_()
	#_plot_local_FLC_haplotype_()
	#analyzeSNPs()
	#_process_raw_flc_()
	#compareSNPs()
	#_generate_250K_2010_FLC_data_()
	_insert_merged_data_in_db_()
	#_insert_markers_into_db_()
	#_impute_FLC_192_()
	print "Done!"
	
	
	
	
