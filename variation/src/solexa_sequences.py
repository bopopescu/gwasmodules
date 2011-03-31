"""
Solexa sequencing stuff.
"""
import phenotypeData as pd


def _read_seeds_files_(file_names = ["seeds_Nordborg_box1.csv","seeds_Nordborg_box2.csv","seeds_Nordborg_box3.csv","seeds_Nordborg_box4.csv"],
		file_dir = "/Users/bjarnivilhjalmsson/Projects/Solexa_sequencing/seeds/"):
	import csv
	ecotypes = []
	accessions = []
	for fn in file_names:
		f = file_dir+fn
		print "loading",f
		r = csv.reader(open(f,"r"), delimiter=',')
		for row in r:
			try:
				ecotypes.append(int(row[0]))
				accessions.append(row[1].strip())
			except Exception, err_str:
				print err_str
	l = zip(ecotypes,accessions)
	l.sort()
	print "Loaded", len(l), "accessions."
	return {"ecotypes":ecotypes,"accessions":accessions}


def _read_MPI_SALK_sequences_():
	file_dir = "/Users/bjarnivilhjalmsson/Projects/Solexa_sequencing/"
	file_names = ["MPI_1001_sequences.csv","salk_accessions_111209.csv"]
	import csv
	stock_parents = []
	accessions = []
	for fn in file_names:
		f = file_dir+fn
		print "loading",f
		r = csv.reader(open(f,"r"), delimiter=',')
		for row in r:
			if row[0].strip()!='Accession' :
				stock_parents.append(row[1].strip())
				accessions.append(row[0].strip())
	l = zip(accessions,stock_parents)
	s = set(l)
	l = list(s)
	l.sort()
	l = map(list,zip(*l))
	accessions = l[0]
	stock_parents = l[1]
	
	l = zip(accessions,stock_parents)
	for acc,sp in l:
		if accessions.count(acc)>1:
			print acc,":", sp
	print "Loaded", len(l), "accessions.... or",len(set(accessions))
	return {"stock_parents":stock_parents,"accessions":accessions}
	
	
#def _read_199_accessions_():
	
	
def _run_():
	l_192 = [6909,6977,100000,6906,8266,6897,6898,5837,6907,7438,6910,6913,6914,6918,6919,8214,
		6924,8424,6926,6928,6933,7520,7521,6936,7522,6937,6939,6900,6901,6908,6009,6915,6917,
		6920,6922,6923,6927,6929,6930,6940,6942,6943,7518,6946,8213,6951,6958,6959,7525,6961,
		6967,6973,6974,6976,7516,6979,6980,6982,6983,6985,6931,6043,6945,7519,7526,7523,6956,
		6960,7524,6963,6964,6965,6966,6969,6971,6975,7517,6978,6981,6984,6899,6903,6904,6905,
		6911,6916,8215,6921,6932,6046,6944,7515,7514,6962,6968,6972,6970,8329,7163,8258,8259,
		8290,7461,7323,8254,8270,8233,8285,6016,8423,8237,6040,6064,6957,8369,8247,8426,9058,
		8249,9057,6709,7000,7062,7460,7123,7147,7255,7275,8241,6988,8256,8264,8265,8231,8271,
		8274,8275,8420,8283,8284,6008,8422,8296,8297,8300,8235,8306,8310,8236,8311,8314,8239,
		8240,8323,8242,8325,8326,8222,8430,6042,8335,8343,6074,8351,8353,8354,7296,8365,8374,
		8376,8378,8412,8387,8389,6243,7306,7418,8312,8313,8334,8337,8357,8366,8411,8388,8395,
		7014,7081,8243,8245,7033,7064,7094,7424,7231,7282,7477,7346,8230]
	print len(l_192)
	
	acc_info = _read_seeds_files_()
	et_list = acc_info["ecotypes"]
	spi = pd._get_stock_parent_info_dict_()
	ei_dict = pd._getEcotypeIdInfoDict_()
	usc_ecotypes = _read_seeds_files_(file_names=["USC_seeds_info_from Joy.csv"],file_dir = "/Users/bjarnivilhjalmsson/Projects/Solexa_sequencing/")["ecotypes"]
	diff_set = set(et_list).difference(set(et_list).intersection(set(usc_ecotypes)))
	for e in diff_set:
		print ei_dict[e]
	ms_acc_info = _read_MPI_SALK_sequences_()
	ms_accessions = ms_acc_info["accessions"]
	ms_stock_parents = ms_acc_info["stock_parents"]	
	print len(spi)
	a2e = pd._getAccessionToEcotypeIdDict_(ms_accessions,)
	len(a2e)
	ecotypes = []
	for acc,sp in zip(ms_accessions,ms_stock_parents):
		if acc in a2e:
			ecotypes.append(int(a2e[acc]))
		elif sp in spi:
			ecotypes.append(int(spi[sp][0]))
		else:
			ecotypes.append(None)
			print acc,",", sp,": weren't found!!" 
	print len(ecotypes)-ecotypes.count(None),len(ecotypes)
	tg_e_dict = pd._getEcotype2TgEcotypeDict_()
	ms_e_set = set(ecotypes)
	print len(ms_e_set)
	elist = []
	for e in et_list:
		if e:
			elist.append(tg_e_dict[e])
		else:
			elist.append(None)
	
	e_set = set(elist)
	i_set = ms_e_set.intersection(e_set)
	print e_set
	print ms_e_set
	print i_set
	print ei_dict[list(i_set)[0]]
	
	ecotypes_192 = set(l_192)
	i_set1 = ecotypes_192.intersection(ms_e_set)
	i_set2 = ecotypes_192.intersection(e_set)
	s = ecotypes_192.intersection(ms_e_set.union(e_set))
	ds = ecotypes_192.difference(s)
	print len(i_set1), len(i_set2), len(s), len(ds)
	for e in ds:
		print ei_dict[e]
	
	
	f = open("/tmp/missing_gwas.csv","w")
	f.write("ecotype_id,accession_name,stock_parent,country_of_origin\n")
	for e in ds:
		f.write(str(e)+","+ei_dict[e][0]+","+ei_dict[e][1]+","+ei_dict[e][4]+"\n")
	f.close()
	
	
	
	

if __name__ == "__main__":
	_run_()
	print "Done!"
	
	