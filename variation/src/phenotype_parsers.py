"""
Contains helper functions for parsing various phenotype files.
"""
import csv
import phenotypeData as pd
from env import *

def load_phentoype_file(filename):
	"""
	Load a FLC type phenotype data file.
	"""
	print "Loading phenotype file:", filename
	f = open(filename, "r")
	reader = csv.reader(f)
	phenotype_names = reader.next()[3:]
	#print phenotype_names
	accession_names = []
	accession_ID = []
	phenotypes = [[] for i in range(len(phenotype_names))]	#[phenotype_name][acc_id]
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
	acc_dict["cibc-5"] = 6908
	acc_dict["pla-0"] = 8357
	#print acc_dict
	#accession_names.sort()
	new_phenotypes = [[] for i in range(len(phenotype_names))]
	ecotypes = []
	for acc in acc_dict:
		acc_i = accession_names.index(acc)
		ecotypes.append(acc_dict[acc])
		for i in range(len(phenotype_names)):
			new_phenotypes[i].append(phenotypes[i][acc_i])
	#print new_phenotypes
	#print len(ecotypes)
	#return {"phenotypes":new_phenotypes,"phenotype_names":phenotype_names, "ecotypes":ecotypes}
	phenotypes = map(list, zip(*new_phenotypes))
	phend = pd.PhenotypeData(ecotypes, phenotype_names, phenotypes)
	phend.writeToFile("/tmp/FLC_phenotypes_102809.tsv", delimiter="\t")
	return phend


def load_phentoype_file_Pecinka():
	accession_file = "/Users/bjarnivilhjalmsson/Projects/Ales_Pecinka/NatVar-AP-2010-Feb.csv"
	f = open(accession_file, "r")
	reader = csv.reader(f)
	accession_names = []
	accession_ID = []
	for row in reader:
		accession_names.append(row[1].split()[0].lower())
		accession_ID.append("CS" + row[0][1:])
	f.close()
	print accession_names
	acc_dict = pd._getAccessionToEcotypeIdDict_(accession_names + ["n13", "kno-10", "kno-10", "shahdara", "nd-1"])
	acc_dict["cibc-5"] = 6908
	acc_dict["wa-1"] = 6978
	acc_dict["gu-0"] = 6922
	acc_dict["cs22491"] = acc_dict["n13"]
	acc_dict["knox-10"] = acc_dict["kno-10"]
	acc_dict["knox-18"] = acc_dict["kno-10"]
	acc_dict["shakdara"] = acc_dict["shahdara"]
	acc_dict["wd-1"] = acc_dict["nd-1"]
	print acc_dict


	filename = "/Users/bjarnivilhjalmsson/Projects/Ales_Pecinka/NatVar-AP-2010-Feb_phen.csv"
	#phenotype_names = reader.next()[2:]
	phenotype_names = ["Absolute_root_growth", "Absolute_root_growth_sd", "Percentage_of_root_elongation", "Percentage_of_bent roots",
			   "Percentage_of_dead_plants", "Percentage_of_unaffected_plants", "Percentage_of_average_survival"]
	phenotype_indices = [1, 2, 5, 8, 11, 14, 17]
	phenotype_ecotypes = [0, 0, 4, 7, 10, 13, 16]
	print phenotype_names
	ecotype_ids = [[] for i in range(len(phenotype_names))]
	phenotypes = [[] for i in range(len(phenotype_names))]	#[phenotype_name][acc_id]
	f = open(filename, "r")
	reader = csv.reader(f)
	new_ecotype_ids = set()
	for row in reader:
		print row
		for i, (pi, ei) in enumerate(zip(phenotype_indices, phenotype_ecotypes)):
			if row[ei] != "":
				acc_name = (row[ei].split()[0]).lower()
				if acc_name in acc_dict:
					eid = acc_dict[(row[ei].split()[0]).lower()]
					new_ecotype_ids.add(eid)
					pv = float(row[pi])
					ecotype_ids[i].append(eid)
					phenotypes[i].append(pv)
				else:
					print "Wrong accession name?", acc_name


	new_phenotypes = []
	new_ecotype_ids = list(new_ecotype_ids)
	for i, phen_vals in enumerate(phenotypes):
		new_phen_vals = []
		for ei in new_ecotype_ids:
			if ei in ecotype_ids[i]:
				j = ecotype_ids[i].index(ei)
				new_phen_vals.append(phen_vals[j])
			else:
				new_phen_vals.append('NA')
		new_phenotypes.append(new_phen_vals)
	phenotypes = map(list, zip(*new_phenotypes))
	ecotypes = map(str, new_ecotype_ids)
	phed = pd.PhenotypeData(ecotypes, phenotype_names, phenotypes)
	phed.writeToFile("/Users/bjarnivilhjalmsson/Projects/Ales_Pecinka/phen_pecinka_170310.tsv", delimiter='\t')


def load_phentoype_file_wilczek():
	filename = "/Users/bjarnivilhjalmsson/Projects/Amity_Wilczek/PhenotypeDataWilczek.csv"
	f = open(filename, "r")
	reader = csv.reader(f)
	phenotype_names = reader.next()[2:]
	for i in range(len(phenotype_names)):
		phenotype_names[i] = phenotype_names[i].replace(" ", "_")
	print phenotype_names
	accession_names = []
	accession_ID = []
	for row in reader:
		accession_names.append(row[1].split()[0].lower())
		accession_ID.append(row[0])
	f.close()
	print accession_names
	acc_dict = pd._getAccessionToEcotypeIdDict_(accession_names)#+["n13","kno-10","kno-10","shahdara","nd-1"])
	acc_dict["cibc-5"] = 6908
	acc_dict["wa-1"] = 6978
	acc_dict["gu-0"] = 7149
	acc_dict['Rubezhnoe-1'] = 7323
	print len(acc_dict), acc_dict
	import env
	d250k_file = env.home_dir + "Projects/Data/250k/250K_t54.csv"
	import dataParsers
	d250k_sd = dataParsers.parse_snp_data(d250k_file)
	ecotypes = []
	key_file = "/Users/bjarnivilhjalmsson/Projects/Amity_Wilczek/unique_id_to_ecotype_id.csv"
	f = open(key_file, "w")
	f.write("unique_id, accession_name, ecotype_id, in_250k_data\n")
	for acc, acc_id in zip(accession_names, accession_ID):
		if not acc in acc_dict or acc_id == 'karl27' or acc_id == 'karl05':
			print "(%s, %s) is missing" % (acc, acc_id)
		else:
			ecotype = acc_dict[acc]
			ecotypes.append(ecotype)
			f.write("%s,%s,%s,%s\n" % (acc_id, acc, str(ecotype), str(str(ecotype) in d250k_sd.accessions)))
	f.close()

	#phenotype_names = reader.next()[2:]
	phenotype_indices = range(2, len(phenotype_names) + 2)
	phenotypes = []	#[acc_id][phenotype_name]
	f = open(filename, "r")
	reader = csv.reader(f)
	reader.next()

	for row in reader:
		#print row
		if row[1].split()[0].lower() in acc_dict:
			phen_vals = []
			for pv in row[2:]:
				if pv == "":
					pv = 'NA'
				else:
					pv = float(pv)
				phen_vals.append(pv)
			phenotypes.append(phen_vals)
		else:
			print "Missing:", row[1]

	phed = pd.PhenotypeData(ecotypes, phenotype_names, phenotypes)
	phed.writeToFile("/Users/bjarnivilhjalmsson/Projects/Amity_Wilczek/phen_wilzcek_050710.tsv", delimiter='\t')
	phed.writeToFile("/Users/bjarnivilhjalmsson/Projects/Amity_Wilczek/phen_wilzcek_050710.csv", delimiter=',')



def load_phentoype_file_dilkes():
	filename = env['phen_dir'] + 'dilkes_metabolites.csv'
	f = open(filename, "r")
	print f.next()
	header = f.next()
	phen_names = map(str.strip, header.split(',')[2:])
	print phen_names
	pids = range(len(phen_names))
	accessions = []
	phen_dict = {}
	for pid, name in zip(pids, phen_names):
		phen_dict[pid] = {'name':name, 'ecotypes':[], 'values':[]}
	for line in f:
		l = map(str.strip, line.split(','))
		accessions.append(l[1].lower())
		for i in pids:
			phen_dict[i]['values'].append(float(l[2 + i]))

	f.close()
	acc_dict = pd.get_250K_accession_to_ecotype_dict()
#	acc_dict['buckhorn'] = [0, 0, 0, 0, 7033]
	acc_dict['sakhdara'] = [0, 0, 0, 0, 6962]
	ecotypes = []
	for acc in accessions:
		if not acc in acc_dict:
			print "%s is missing in dictionary" % acc
		else:
			ecotype = acc_dict[acc][4]
			ecotypes.append(ecotype)

	for pid, name in zip(pids, phen_names):
		phen_dict[pid]['ecotypes'] = ecotypes

	phed = pd.phenotype_data(phen_dict)
	phed.write_to_file(env['phen_dir'] + 'b_dilkes_metabolites.csv')




def load_phentoype_file_nc_resistance():
	filename = "/Users/bjarnivilhjalmsson/Summary_results_330Arabidopsis_accessions.csv"
	f = open(filename, "r")
	line = map(str.strip, f.next().split(','))
	phenotype_names = line[-2:]
	print phenotype_names
	phenotypes = []
	accession_names = []
	full_accession_names = []
	for l in f:
		line = map(str.strip, l.split(','))
		accession_names.append(line[0].lower())
		full_accession_names.append(line[2].lower())
	f.close()
	print accession_names
	acc_dict = pd.get_accession_to_ecotype_id_dict(accession_names)#+["n13","kno-10","kno-10","shahdara","nd-1"])
#	acc_dict["cibc-5"] = 6908
#	acc_dict["wa-1"] = 6978
#	acc_dict["gu-0"] = 7149
#	acc_dict['Rubezhnoe-1'] = 7323
	print len(acc_dict), acc_dict
	import env
	d250k_file = env.home_dir + "Projects/Data/250k/250K_t54.csv.binary"
	import dataParsers
	d250k_sd = dataParsers.parse_binary_snp_data(d250k_file)
	ecotypes = []
	key_file = "/Users/bjarnivilhjalmsson/Projects/Amity_Wilczek/unique_id_to_ecotype_id.csv"
	f = open(key_file, "w")
	f.write("unique_id, accession_name, ecotype_id, in_250k_data\n")
	for acc, acc_id in zip(accession_names, full_accession_names):
		if not acc in acc_dict or acc_id == 'karl27' or acc_id == 'karl05':
			print "(%s, %s) is missing" % (acc, acc_id)
		else:
			ecotype = acc_dict[acc]
			ecotypes.append(ecotype)
			f.write("%s,%s,%s,%s\n" % (acc_id, acc, str(ecotype), str(str(ecotype) in d250k_sd.accessions)))

	#phenotype_names = reader.next()[2:]
	phenotypes = []	#[acc_id][phenotype_name]
	f = open(filename, "r")

	for l in f:
		line = map(str.strip, l.split(','))
		if line[0].lower() in acc_dict:
			phen_vals = []
			for pv in line[-2:]:
				if pv == "NA":
					pv = 'NA'
				else:
					pv = float(pv)
				phen_vals.append(pv)
			phenotypes.append(phen_vals)
		else:
			print "Missing:", line[0]

	phed = pd.PhenotypeData(ecotypes, phenotype_names, phenotypes)
	phed.insert_into_DB(growth_condition='Field', biology_category_id='2')
	phed.writeToFile("/Users/bjarnivilhjalmsson/Projects/Amity_Wilczek/nc14_resistance_091610.tsv", delimiter='\t')



def load_phentoype_file_nc_resistance_2():
	filename = "/Users/bjarnivilhjalmsson/Projects/Albugo_laibachii_nc14.csv"
	f = open(filename, "r")
	line = map(str.strip, f.next().split(','))
	phenotype_names = line[-1:]
	print phenotype_names
	phenotypes = []
	accession_names = []
	for l in f:
		line = map(str.strip, l.split(','))
		accession_names.append(line[0].lower())
	f.close()
	print accession_names
	acc_dict = pd.get_accession_to_ecotype_id_dict(accession_names)#+["n13","kno-10","kno-10","shahdara","nd-1"])
#	acc_dict["cibc-5"] = 6908
#	acc_dict["wa-1"] = 6978
#	acc_dict["gu-0"] = 7149
#	acc_dict['Rubezhnoe-1'] = 7323
	print len(acc_dict), acc_dict
	import env
	d250k_file = env.home_dir + "Projects/Data/250k/250K_t54.csv.binary"
	import dataParsers
	d250k_sd = dataParsers.parse_binary_snp_data(d250k_file)
	ecotypes = []
	for acc in accession_names:
		if not acc in acc_dict:
			print "%s is missing" % (acc)
		else:
			ecotype = acc_dict[acc]
			ecotypes.append(ecotype)

	phenotypes = []	#[acc_id][phenotype_name]
	f = open(filename, "r")

	for l in f:
		line = map(str.strip, l.split(','))
		if line[0].lower() in acc_dict:
			phen_vals = []
			for pv in line[-1:]:
				if pv == "NA":
					pv = 'NA'
				else:
					pv = float(pv)
				phen_vals.append(pv)
			phenotypes.append(phen_vals)
		else:
			print "Missing:", line[0]

	phed = pd.PhenotypeData(ecotypes, phenotype_names, phenotypes)
	phed.insert_into_DB(growth_condition='Field', biology_category_id='2')
	phed.writeToFile("/Users/bjarnivilhjalmsson/Projects/Amity_Wilczek/nc14_resistance_96accessions_092810.tsv", \
			delimiter='\t')





def load_phentoype_file_bergelsson():
	import env
	filename = "/Users/bjarnivilhjalmsson/Projects/Joy_Bergelsson/bergelsson_rosette_glucs.csv"
	f = open(filename, "r")
	reader = csv.reader(f)
	phenotype_names = reader.next()[2:]
	for i in range(len(phenotype_names)):
		phenotype_names[i] = phenotype_names[i].replace(" ", "_")
		phenotype_names[i] = 'jb_' + phenotype_names[i]
	print phenotype_names
	accession_names = []
	accession_ID = []
	for row in reader:
		accession_names.append(row[0].split()[0].lower())
		accession_ID.append(row[1])
	f.close()
	print accession_names
	#acc_dict = pd._getAccessionToEcotypeIdDict_(accession_names)#+["n13","kno-10","kno-10","shahdara","nd-1"])
	e_info_dict = pd._getEcotypeIdInfoDict_()
	ei_2_tgei = pd._getEcotype2TgEcotypeDict_()
	#print len(acc_dict),acc_dict
	ecotypes = []
        uncertain_list = []
	for acc, acc_id in zip(accession_names, accession_ID):
		#if not acc in acc_dict:
		if not int(acc_id) in ei_2_tgei:
			print "(%s, %s) is missing in dictionary" % (acc, acc_id)
			a_id = int(acc_id)
			if a_id in e_info_dict:
				e_info = e_info_dict[a_id]
				print "Guessing that it's:", e_info
			else:
				print "No good guess for it.  Look it up!!\n"
			#acc_dict[acc] = acc_id
			ecotypes.append(acc_id)
		else:
			#ecotype = acc_dict[acc]
			ecotype = ei_2_tgei[int(acc_id)]
			ecotypes.append(ecotype)
	phenotype_indices = range(2, len(phenotype_names) + 2)
	phenotypes = []	#[acc_id][phenotype_name]
	f = open(filename, "r")
	reader = csv.reader(f)
	reader.next()

	print len(set(accession_ID)), len(set(ecotypes))

	for row in reader:
		#print row
		#if row[0].split()[0].lower() in acc_dict:
			phen_vals = []
			for pv in row[2:]:
				if pv == "" or pv == 'NA':
					pv = 'NA'
				else:
					pv = float(pv)
				phen_vals.append(pv)
			if len(phen_vals) != len(phenotype_names):
				import pdb;
				pdb.set_trace()
			phenotypes.append(phen_vals)
		#else:
		#	print "Missing:",row[0]


	phed = pd.PhenotypeData(ecotypes, phenotype_names, phenotypes)
	phed.writeToFile("/Users/bjarnivilhjalmsson/Projects/Joy_Bergelsson/phen_bergelsson_051710.tsv", delimiter='\t')
	phed.writeToFile("/Users/bjarnivilhjalmsson/Projects/Joy_Bergelsson/phen_bergelsson_051710.csv", delimiter=',')






def load_phentoype_file_duszynska():
	fn1 = env['phen_dir'] + 'seed_size_2n.csv'
	fn2 = env['phen_dir'] + 'seed_size_3n_2x4.csv'
	fn3 = env['phen_dir'] + 'seed_size_3n_4x2.csv'
	fn4 = env['phen_dir'] + 'seed_size_spss.csv'
	fns = [fn1, fn2, fn3, fn4]
	phen_names = ['seed_size_2n', 'seed_size_3n_2x4', 'seed_size_3n_4x2', 'seed_size_spss']
	phen_dict = {}
	for i, pn in enumerate(phen_names):
		phen_dict[i + 1] = {'name':pn }
	accs_list = []
	for i, fn in enumerate(fns):
		f = open(fn, "r")
		acc_dict = pd.get_250K_accession_to_ecotype_dict()
		ecotypes = []
		values = []
		print f.next()
		for line in f:
			l = map(str.strip, line.split(','))
			acc = l[0].lower()
			if not acc in acc_dict:
				print "(%s) is missing in dictionary" % (acc)
			else:
				ecotypes.append(acc_dict[acc][4])
				values.append(float(l[1]))
		f.close()
		phen_dict[i + 1]['ecotypes'] = ecotypes
		phen_dict[i + 1]['values'] = values


	phed = pd.phenotype_data(phen_dict)
	phed.write_to_file(env['phen_dir'] + 'seed_size.csv')




def load_phentoype_file_riha():
	filename = env['phen_dir'] + 'telomere_lengths_192_raw.csv'
	f = open(filename, "r")
	phen_name = 'telomere_length'
	accession_names = []
	accession_ids = []
	parent_ids = []
	phen_vals = []
	print f.next()
	for line in f:
		l = map(str.strip, line.split(','))
		parent_ids.append(l[0])
		acc_l = l[1].split()
		acc_name = acc_l[0]
		if len(acc_l) > 1:
			acc_id = acc_l[1]
		else:
			acc_id = ''
		accession_names.append(acc_name.lower())
		accession_ids.append(acc_id)
		phen_vals.append(float(l[2]))

	f.close()
	print accession_names
	acc_dict = pd.get_250K_accession_to_ecotype_dict()
	acc_dict['buckhorn'] = [0, 0, 0, 0, 7033]
	acc_dict['shahdara'] = [0, 0, 0, 0, 6962]
	ecotypes = []
        uncertain_list = []
        new_phen_vals = []
	for acc, par_id, acc_id, phen_val in zip(accession_names, parent_ids, accession_ids, phen_vals):
		if not acc in acc_dict:
			print "(%s, %s, %s) is missing in dictionary" % (acc, par_id, acc_id)
		else:
			ecotype = acc_dict[acc][4]
			ecotypes.append(ecotype)
			new_phen_vals.append(phen_val)

	print len(set(accession_names)), len(set(ecotypes))
	phen_dict = {1:{'name':phen_name, 'ecotypes':ecotypes, 'values':new_phen_vals}}

	phed = pd.phenotype_data(phen_dict)
	phed.write_to_file(env['phen_dir'] + 'telomere_lengths_192.csv')





def add_phenotypes_to_db(phenotypes, phenotype_names, ecotypes, method_ids, method_descriptions=None, data_descriptions=None, host="papaya.usc.edu", user="bvilhjal", passwd="bjazz32"):
	"""
	Insert phenotypes into the DB
	
	Use the format from above...
	"""
	import MySQLdb
	print "Connecting to db, host=" + host
	if not user:
		import sys
		sys.stdout.write("Username: ")
		user = sys.stdin.readline().rstrip()
	if not passwd:
		import getpass
		passwd = getpass.getpass()
	try:
		conn = MySQLdb.connect (host=host, user=user, passwd=passwd, db="at")
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()
	#Retrieve the filenames
	print "Inserting data"

	for i in range(len(phenotype_names)):
		sql_statement = "INSERT INTO stock_250k.phenotype_method  (id, short_name, only_first_96, biology_category_id, method_description, data_description, data_type) VALUES (" + str(method_ids[i]) + ", '" + phenotype_names[i] + "', true, 1, '" + method_descriptions[i] + "', '" + data_descriptions[i] + "', 'quantitative')"
		print sql_statement
		numRows = int(cursor.execute(sql_statement))
		row = cursor.fetchone()
		if row:
			print row

		for j in range(0, len(ecotypes)):
			val = phenotypes[i][j]
			if val != "NA":
				sql_statement = "INSERT INTO stock_250k.phenotype_avg (ecotype_id, value, ready_for_publication, method_id, transformed_value) VALUES ( " + str(ecotypes[j]) + ", " + str(val) + ", 0, " + str(method_ids[i]) + ", " + str(val) + " )"
				print sql_statement
				numRows = int(cursor.execute(sql_statement))
				row = cursor.fetchone()
				if row:
					print row
	conn.commit()
	cursor.close ()
	conn.close ()

def _run_():
	pd = load_phentoype_file("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/data_102509/FLC_soil_data_102509.csv")
	phenotypes = pd["phenotypes"]
	phenotype_names = pd["phenotype_names"]
	new_phenotype_names = []
	for pn in phenotype_names:
		new_phenotype_names.append("FLC_v2_" + pn)
	phenotype_names = new_phenotype_names
	ecotypes = pd["ecotypes"]
	method_descriptions = ["without vernalization", "4 weeks after vernalization in cold", "10 days after 4 weeks of vernalization", "30 days after 4 weeks of vernalization", "Ratio between 4wT0 and NV", "Ratio between 4wT10 and NV", "Ratio between 4wT30 and 4wT0", "Ratio between 4wT30 and 4wT10"]
	data_descriptions = ["Stable soil growth condition and new Roch PCR machine." for i in range(len(phenotype_names))]
	add_phenotypes_to_db(phenotypes, phenotype_names, ecotypes, method_ids=range(352, 352 + len(phenotype_names)), method_descriptions=method_descriptions, data_descriptions=data_descriptions)

if __name__ == "__main__":
	#_run_()
	#load_phentoype_file("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/data_102509/FLC_soil_data_102509.csv")
	#load_phentoype_file_Pecinka()
	#load_phentoype_file_wilczek()
	load_phentoype_file_duszynska()
	print "Done!"


