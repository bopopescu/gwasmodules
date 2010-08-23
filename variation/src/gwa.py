#!/usr/bin/env python2.5
"""
Usage: gwa.py [OPTIONS] [--id=id_text] SNPS_DATA_FILE PHENOTYPE_DATA_FILE [PHENOTYPE_INDEX]

Option:

	--id=...				ID string, used for output files.
	-d ..., --delim=...			default is ", "	  
	-m ..., --missingval=...		default is "NA"
	-h, --help				show this help
	--parallel=...				Run mapping methods on the cluster with standard parameters.  The argument is used for runid 
						as well as output files.  Note that if this option is used then no output file should be specified.
						If phenotype index is missing, then all phenotypes are used.
	--specific_methods=...			Apply specific methods, otherwise all available are applied:
						emma,emma_trans,kw,ft, etc.
	--use_existing_results			Use existing results when possible, to speed up.  (Otherwise existing files are overwritten)

	--region_plots=...			Include region plots for the top N (given num) peaks.
	--cand_genes_file=...			A file with a list of candidate genes.  (added to plots)
	
	--analysis_plots			Generate analysis plots, histograms, QQ-plots, test different transformations, etc.

	--addToDB				Adds the result(s) to the DB.
	--callMethodID=...			What data set is used. (Only applicable if result is added to DB.)
	--comment=...				Comment for DB. (Only applicable if result is added to DB.)
	
	--phenotype_w_db_id			Phenotypes have DB id as an prefix in their names.
	
	--memReq=...				Request memory (on cluster), otherwise use defaults
	--walltimeReq=...			Request time limit (on cluster), otherwise use defaults
	--proc_per_node=...			Request number of processors per node, default is 8. (Works only for EMMA.)
	--debug_filter=...			Filter SNPs randomly (speed-up for debugging)
								
						

Examples:
~/gwas_data$ python /home/GMI/bjarni.vilhjalmsson/Projects/py_src/gwa.py --id=../gwas_results/pi1 --specific_methods=kw --callMethodID=54 --phenotype_w_db_id --proc_per_node=8 250K_t54.csv phen_all_052010.tsv	
Description:
  Applies various GWA methods to to phenotypes.

  Methods include: Kruskal-Wallis, Fisher's Exact, Emma, etc.
"""
#from epydoc.docparser import str_prefixes
import sys, getopt, traceback
import os

import phenotypeData
import dataParsers
import snpsdata
import gwaResults
import util
import warnings
import gc
import multiprocessing as mp
from Queue import Empty
import time

from numpy import *
#import rpy2.robjects as robjects
#import rpy2.robjects.numpy2ri
#r_source = robjects.r('source')

#import rpy2.rpy_classic as rpy  
#rpy.set_default_mode(rpy.NO_CONVERSION)
#rpy.set_default_mode(rpy.CLASS_CONVERSION)
#Switching to rpy2

#HPC CLUSTER SETTING
#results_dir="/home/cmb-01/bvilhjal/results/"
#data_dir = '/home/cmb-01/bvilhjal/data/'
#scriptDir="/home/cmb-01/bvilhjal/Projects/Python-snps/"
#results_dir = '/home/cmbpanfs-01/bvilhjal/results/'
#data_dir = '/home/cmbpanfs-01/bvilhjal/data/'
#script_dir = '/home/cmbpanfs-01/bvilhjal/src/'

#LOCAL BETULA SETTING
#data_dir = '/tmp/'
#script_dir="/Users/bjarnivilhjalmsson/Projects/py_src/"

from env import *

analysis_methods_dict = {"kw":1,
			 "ft":2,
			 "emma":4,
			 "emma_trans":47,
			 "emma_trans_no":49}

def prepare_data(sd,phed,p_i,mapping_method):
	"""
	Coordinates phenotype and snps data for different mapping methods.
	"""
	num_outliers_removed = 0
	if "emma_trans"== mapping_method:
		phed.standard_transform(p_i)

	elif "emma_trans_no"== mapping_method:
		phed.standard_transform(p_i)
		num_outliers_removed = phed.naOutliers(p_i,iqrThreshold=3)
	sd.coordinate_w_phenotype_data(phed,p_i)
	return num_outliers_removed	
	
	

def get_perm_pvals(snps,phen_vals,mapping_method='kw',num_perm=100,snps_filter=0.05):
	import random
	if snps_filter<1.0:
		snps = random.sample(snps,int(snps_filter*len(snps)))
	pvals = []
	if mapping_method=='kw':
		for i in range(num_perm):
			random.shuffle(phen_vals)
			kw_res = util.kruskal_wallis(snps,phen_vals)
			pvals.extend(kw_res['ps'])	
	elif mapping_method=='ft':
		for i in range(num_perm):
			random.shuffle(phen_vals)
			pvals.extend(run_fet(snps, phen_vals))	
	return pvals
		
		

def _run_():
	if len(sys.argv)==1:
		print __doc__
		sys.exit(2)
	
	long_options_list=["id=", "delim=", "missingval=", "help", "parallel=", 
			   "addToDB", "callMethodID=", "comment=","memReq=",
			   "walltimeReq=","specific_methods=","analysis_plots",
			   "use_existing_results","region_plots=","cand_genes_file=",
			   "debug_filter=", "phenotype_w_db_id", "proc_per_node="]
	try:
		opts, args=getopt.getopt(sys.argv[1:], "d:m:h", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)	
	
	run_id=None
	delim=","
	missingVal="NA"
	help=0
	parallel=None
	addToDB=False
	callMethodID=None
	comment=""
	
	memReq = "11800mb"
	walltimeReq = "72:00:00"
	proc_per_node = 8
	
	specific_methods = ['kw','ft','emma',"emma_trans","emma_trans_no"]
	analysis_plots = False
	use_existing_results = False
	region_plots = 0
	cand_genes_file=None
	debug_filter=1
	phenotype_w_db_id = False

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help=1
			print __doc__
		elif opt in ("--id"):
			run_id=arg
		elif opt in ("--parallel"):
			parallel=arg
		elif opt in ("--addToDB"):
			addToDB=True
		elif opt in ("--callMethodID"):
			callMethodID=int(arg)
		elif opt in ("--comment"):
			comment=arg
		elif opt in ("-d", "--delim"):
			delim=arg
		elif opt in ("-m", "--missingval"):
			missingVal=arg
		elif opt in ("--memReq"):
			memReq=arg
		elif opt in ("--walltimeReq"):
			walltimeReq=arg
		elif opt in ("--proc_per_node"):
			proc_per_node=int(arg)
		elif opt in ("--specific_methods"):
			specific_methods = arg.split(",")
		elif opt in ("--analysis_plots"):
			analysis_plots = True
		elif opt in ("--use_existing_results"):
			use_existing_results = True
		elif opt in ("--region_plots"):
			region_plots = int(arg)
		elif opt in ("--cand_genes_file"):
			cand_genes_file = arg
		elif opt in ("--debug_filter"):
			debug_filter = float(arg)
		elif opt in ("--phenotype_w_db_id"):
			phenotype_w_db_id = True
		else:
			if help==0:
				print "Unkown option!!\n"
				print __doc__
			sys.exit(2)

	if len(args)<3 and not parallel:
		if help==0:
			print "Arguments are missing!!\n"
			print __doc__
		sys.exit(2)

	snpsDataFile=args[0]
	phenotypeDataFile=args[1]

	print "GWA runs are being set up with the following parameters:"
	print "run_id:",run_id
	print "phenotypeDataFile:",phenotypeDataFile
	print "snpsDataFile:",snpsDataFile
	print "parallel:",parallel
	print "walltimeReq:",walltimeReq
	print "proc_per_node:",proc_per_node
	print "memReq:",memReq	
	print "specific_methods:",specific_methods
	print "analysis_plots:",analysis_plots
	print "addToDB:",addToDB
	print "callMethodID:",callMethodID
	print "comment:",comment
	print "use_existing_results:",use_existing_results
	print "region_plots:",region_plots
	print "cand_genes_file:",cand_genes_file
	print "debug_filter:",debug_filter
	print 'phenotype_w_db_id:',phenotype_w_db_id

	def run_parallel(p_i,phed,mapping_method="analysis"):
		"""
		If no mapping_method, then analysis run is set up.
		"""

		if not p_i in phed.phenIds:
			print "Phenotype ID not found:%i"%p_i
			return 
		if phed.isBinary(p_i):
			if mapping_method in ["kw","emma_trans","emma_trans_no"]:
				return
			elif analysis_plots:
				specific_methods = ["ft","emma"]
		else:
			if mapping_method in ["ft"]:
				return
			elif analysis_plots:
				specific_methods = ["kw","emma", "emma_trans","emma_trans_no"]
			
		phenotype_name=phed.getPhenotypeName(p_i)
		#Cluster specific parameters
		print "Setting up a gwa run for phenotype:%s, pid=%d, using method:%s"%(phenotype_name,p_i,mapping_method)
		run_id=results_dir+parallel+"_pid"+str(p_i)+"_"+phenotype_name+"_"+mapping_method

		shstr = "#!/bin/csh\n"
		shstr += "#PBS -l walltime="+walltimeReq+"\n"
		if mapping_method in ['emma','emma_trans','emma_trans_no']:
			shstr += "#PBS -l nodes=1:ppn=%d \n"%proc_per_node
		shstr += "#PBS -l mem="+memReq+"\n"
		shstr +="#PBS -q cmb\n"
		
		shstr+="#PBS -N p"+str(p_i)+"_"+mapping_method+"_"+parallel+"_"+mapping_method+"\n"
		shstr+="(python "+script_dir+"gwa.py --id="+run_id+"  "
		shstr+="--region_plots="+str(region_plots)+"  "
		shstr+="--proc_per_node=%d  "%proc_per_node
		if use_existing_results:
			shstr+="--use_existing_results  "
		if analysis_plots:
			shstr+="--analysis_plots  --specific_methods="+",".join(specific_methods)+"  "	
		else:	
			shstr += " --specific_methods="+mapping_method+" "
		if cand_genes_file:
			shstr+="--cand_genes_file="+cand_genes_file+"  "
		if phenotype_w_db_id:
			shstr+="--phenotype_w_db_id  "	
		if debug_filter:
			shstr+="--debug_filter="+str(debug_filter)+"  "	
		if addToDB:
			shstr+="--addToDB "
		if callMethodID:
			shstr+="--callMethodID="+str(callMethodID)+"  "
		if comment:
			shstr+="--comment="+str(comment)+"  "
		shstr+=snpsDataFile+" "+phenotypeDataFile+" "+str(p_i)+" "
		shstr+="> "+run_id+"_job"+".out) >& "+run_id+"_job"+".err\n"

		f=open(parallel+".sh", 'w')
		f.write(shstr)
		f.close()

		#Execute qsub script
		os.system("qsub "+parallel+".sh ")

	#Load phenotype file
	phed=phenotypeData.readPhenotypeFile(phenotypeDataFile,with_db_ids=phenotype_w_db_id)  #Get Phenotype data 
	if parallel:
		if len(args)==3:
			t_pids = args[2].split(',')
			pids = []
			for s in t_pids:
				if '-' in s:
					pid_range = map(int,s.split('-'))
				        for pid in range(pid_range[0],pid_range[1]+1):
				        	pids.append(pid)
				else:
					pids.append(int(s))
				        
		if len(args)==2:  #phenotype index arguement is missing, hence all phenotypes are run/analyzed.
			pids = phed.phenIds				
		if analysis_plots:  #Running on the cluster..
			for p_i in pids:
				run_parallel(p_i,phed)
		else:
			for mapping_method in specific_methods:
				for p_i in pids:
					run_parallel(p_i,phed,mapping_method)
		return
	else:
		p_i=int(args[2])

	phen_is_binary = phed.isBinary(p_i)
	phenotype_name=phed.getPhenotypeName(p_i)
	print "Phenotype:%s, phenotype_id:%s"%(phenotype_name, p_i)

	if analysis_plots:
		print "\nAnalysing GWAs results jointly... QQ plots etc."
		
		#Genotype and phenotype data is only used for permutations.
		sd=dataParsers.parse_snp_data(snpsDataFile, format = 0, delimiter = delim, 
					      missingVal = missingVal,filter=debug_filter)
		prepare_data(sd,phed,p_i,'kw')
		snps = sd.getSnps()
		phen_vals = phed.getPhenVals(p_i)		

		print "Plotting accession phenotype map"
		sys.stdout.flush()
		accession_map_pdf_file =  run_id+"_acc_phen_map.pdf"
		phed.plot_accession_map(p_i,pdf_file=accession_map_pdf_file)
		
		#Load gwas results
		results = {}
		perm_pvals = None
		methods_found = []
		for mapping_method in specific_methods:
			try:
				l = run_id.split("_")[:-1]
				file_prefix = "_".join(l)+"_"+mapping_method
				file_name = file_prefix+".pvals"
				sys.stdout.write("Looking for file %s."%(file_name))
				if os.path.isfile(file_name):
					sys.stdout.write("..... found!\n")
					snps = sd.getSnps()
					res = gwaResults.Result(result_file=file_name,name=mapping_method+"_"+\
							        phenotype_name,snps=snps)
					pvals=True
					if mapping_method in ["emma","emma_trans","emma_trans_no"]:
						res.filter_attr("mafs",15)
					results[mapping_method]=res
						

				else:
					sys.stdout.write("..... not found.\n")
					sys.stdout.flush()
					file_name = file_prefix+".scores"
					sys.stdout.write("Looking for file %s."%(file_name))
					if os.path.isfile(file_name):
						sys.stdout.write("..... found!\n")
						snps = sd.getSnps()
						res = gwaResults.Result(result_file=file_name,name=mapping_method+\
								        "_"+phenotype_name, snps=snps)
						pvals=False
						results[mapping_method]=res
					else:
						sys.stdout.write("..... not found.\n")
				sys.stdout.flush()

			
				#Permutation tests for KW and FT..
				if mapping_method in ['kw','ft']:
					print "Retrieving permutations for %s"%(mapping_method)
					sys.stdout.flush()
					perm_pvals = get_perm_pvals(snps,phen_vals,mapping_method)
			except Exception, err_str:
				print err_str
				print "Failed loading file %s"%(file_name)
				sys.stdout.flush()
		if len(results.keys())>0:
			print "Drawing QQ plots for methods: %s"%(str(results.keys()))
			for k in results:
				print k, results[k].name
			sys.stdout.flush()
                        #Plotting the QQ plots.
			num_dots = 1000
			max_log_val = 5
			gwaResults.qq_plots(results, num_dots, max_log_val, file_prefix, method_types=results.keys(),
					    phen_name = phenotype_name, perm_pvalues=perm_pvals, is_binary = phen_is_binary)

	else:
	
		if len(specific_methods)!=1:
			print "More than one specific GWA methods cannot be applied unless running on the cluster."
			return
		mapping_method = specific_methods[0]
		
		#Load genotype file (in binary format)
		sd=dataParsers.parse_snp_data(snpsDataFile, format = 0, delimiter = delim, 
					      missingVal = missingVal,filter=debug_filter)
		num_outliers = prepare_data(sd,phed,p_i,mapping_method)
		if num_outliers==0 and mapping_method=="emma_trans_no":
			print "No outliers were removed, so it makes no sense to go on and perform GWA."
			return
		snps = sd.getSnps()
		phen_vals = phed.getPhenVals(p_i)		
	
		#Check whether result already exists.
		res = None
		if use_existing_results:
			print "\nRetrieving existing results now!"
			result_file = run_id+".pvals"
			if os.path.isfile(result_file):
				res = gwaResults.Result(result_file=result_file,name=mapping_method+"_"+phenotype_name,snps=snps)
				pvals=True
			else:
				result_file = run_id+".scores"
				if os.path.isfile(result_file):
					res = gwaResults.Result(result_file=result_file,name=mapping_method+"_"+phenotype_name, 
							        snps=snps)
					pvals=False
			if res:
				print "Found existing results.. (%s)"%(result_file)
			sys.stdout.flush()
		
		
		if not res: #If results weren't found in a file... then do GWA.
			gc.collect()
			print "\nStarting to apply method(s) now!"
			sys.stdout.flush()
			kwargs = {}
			additional_columns = []
			if "kw" == mapping_method:
				print "Applying KW to data."
				#Writing files
				#phed and phenotype
				
				if phen_is_binary:
					warnings.warn("Warning, applying KW to a binary phenotype")
		
				kw_res = util.kruskal_wallis(snps,phen_vals)
				pvals = kw_res['ps']
				kwargs['statistics'] = kw_res['ds']
				additional_columns.append('statistics')
	
		
			elif "ft" == mapping_method:
				pvals, or_est = run_fet(snps,phen_vals)
				kwargs['odds_ratio_est'] = or_est
				additional_columns.append('odds_ratio_est')
			
			else:  #Parametric tests below:		
				gc.collect()
				if mapping_method in ['emma','emma_trans','emma_trans_no']:
					print "Applying %s to data."%(mapping_method)
					sys.stdout.flush()
					sys.stdout.write("calculating the Kinship matrix K....")
					sys.stdout.flush()
					if callMethodID:
						k_file = data_dir+"kinship_matrix_cm"+str(callMethodID)+".pickled"
					else:
						k_file = None
					k = retrieve_kinship(accessions=sd.accessions,k_file=k_file,sd_file=snpsDataFile,
							 snps=snps,random_fraction=debug_filter)
					#k = calcKinship(snps)
					gc.collect()
					sys.stdout.write("Done!.\n")
					sys.stdout.write("Starting EMMA (fast)\n")
					sys.stdout.flush()
					res = run_emma_parallel(snps,phen_vals,k,num_proc=proc_per_node)
					#res = runEmma(snps,phen_vals,k)
					kwargs['genotype_var_perc'] = [val[0] for val in res['genotype_var_perc']]
					kwargs['beta0'] = [val[0] for val in res['beta0_est']]
					kwargs['beta1'] = [val[0] for val in res['beta1_est']]
					additional_columns.append('genotype_var_perc')
					additional_columns.append('beta0')
					additional_columns.append('beta1')
					pvals = [pval[0] for pval in res['ps']]
					sys.stdout.write("Done!\n")
					sys.stdout.flush()

			
			
			kwargs['correlations'] = calc_correlations(snps, phen_vals)
			additional_columns.append('correlations')
			 
			res = gwaResults.Result(scores=pvals, snps_data=sd, name=mapping_method+"_"+phenotype_name, **kwargs)
										
			if pvals:
			 	result_file=run_id+".pvals"
			else:
			 	result_file=run_id+".scores"	
			res.write_to_file(result_file,additional_columns)
		
		#add results to DB..
		if addToDB:
			if phenotype_w_db_id:
				db_pid = p_i
			else:
				db_pid = phed.get_db_pid(p_i)
			import results_2_db as rdb
			short_name=mapping_method+"_"+phenotype_name+"_call_method"+str(callMethodID)
			rdb.add_results_to_db(result_file,short_name,callMethodID,db_pid,analysis_methods_dict[mapping_method])
			

		#Load candidate genes from a file, if it is given
		cand_genes = None
		if cand_genes_file:
			cand_genes,tair_ids = gwaResults.load_cand_genes_file(cand_genes_file)
		else:
			cand_genes=None
			tair_ids=None
			
		
		print "Generating a GW plot."
		sys.stdout.flush()
		png_file = run_id+"_gwa_plot.png"
		png_file_max30 = run_id+"_gwa_plot_max30.png"
		if pvals:
			res.neg_log_trans()
			if mapping_method in ["kw","ft"]:
				res.plot_manhattan(png_file=png_file_max30,percentile=90,type="pvals",ylab="$-$log$_{10}(p)$", 
					       plot_bonferroni=True,cand_genes=cand_genes,max_score=30)
				res.plot_manhattan(png_file=png_file,percentile=90,type="pvals",ylab="$-$log$_{10}(p)$", 
					       plot_bonferroni=True,cand_genes=cand_genes)
			else:	
				res.filter_attr("mafs",15)
				res.plot_manhattan(png_file=png_file_max30,percentile=90,type="pvals",ylab="$-$log$_{10}(p)$", 
					       plot_bonferroni=True,cand_genes=cand_genes,max_score=30)
				res.plot_manhattan(png_file=png_file,percentile=90,type="pvals",ylab="$-$log$_{10}(p)$", 
					       plot_bonferroni=True,cand_genes=cand_genes)
		else:
			pass
		
		print "plotting histogram"
		if mapping_method in ["kw","ft"]:
			hist_png_file =  run_id+"_hist_raw.png"
			phed.plot_histogram(p_i,pngFile=hist_png_file)
		elif mapping_method=="emma_trans":
			hist_png_file = run_id+"_hist_log_trans.png"
			phed.plot_histogram(p_i,pngFile=hist_png_file)
		elif mapping_method=="emma_trans_no":
			hist_png_file = run_id+"_hist_log_trans_no.png"
			phed.plot_histogram(p_i,pngFile=hist_png_file)
			
			
	
		
		if region_plots:
			import regionPlotter as rp
			regions_results = res.get_top_region_results(region_plots)
			plotter = rp.RegionPlotter()
			print "Starting region plots..."
			for reg_res in regions_results:
				chromosome = reg_res.chromosomes[0]
				caption = phenotype_name+"_c"+str(chromosome)+"_"+mapping_method
				png_file = run_id+"_reg_plot_c"+str(chromosome)+"_s"+str(reg_res.positions[0])+"_e"+str(reg_res.positions[-1])+".png"
				tair_file = run_id+"_reg_plot_c"+str(chromosome)+"_s"+str(reg_res.positions[0])+"_e"+str(reg_res.positions[-1])+"_tair_info.txt"
				plotter.plot_small_result([reg_res],png_file=png_file,highlight_gene_ids=tair_ids,
							  caption=caption,tair_file=tair_file)
				
				#Plot Box-plot
				png_file = run_id+"_reg_plot_c"+str(chromosome)+"_s"+str(reg_res.positions[0])+"_e"+str(reg_res.positions[-1])+"_box_plot.png"
				(marker,score,chromosome,pos) = reg_res.get_max_snp()
				marker_accessions = sd.accessions
				phed.plot_marker_box_plot(p_i,marker=marker,marker_accessions=marker_accessions,png_file=png_file,
						     title="c"+str(chromosome)+"_p"+str(pos),marker_score=score,marker_missing_val=sd.missing_val)
				
				
			
def filter_k_for_accessions(k,indices_to_keep):
	new_k = zeros((len(indices_to_keep),len(indices_to_keep)))
	for i in range(len(indices_to_keep)):
		for j in range(len(indices_to_keep)):
			new_k[i,j]=k[indices_to_keep[i],indices_to_keep[j]]
	k = new_k
	return k

		
		
def retrieve_kinship(accessions=None,k_file=None,sd_file=None,snps=None,random_fraction=None):
	import pickle
	if os.path.isfile(k_file) and accessions:
		sys.stdout.write("Loading K.\n")
		sys.stdout.flush()
		f = open(k_file,'r')
		l = pickle.load(f)
		f.close()
		k = l[0]
		k_accessions = l[1]
		#Filter for used accessions
		indices_to_keep = []
		for i, acc in enumerate(k_accessions):
			if acc in accessions:
				indices_to_keep.append(i)		
		k = filter_k_for_accessions(k,indices_to_keep)
		
	elif sd_file and accessions: 
		#Else generate and save
		sys.stdout.write("Generating K.\n")
		sys.stdout.flush()
		sd = dataParsers.parse_snp_data(sd_file,format = 0,filter=random_fraction)
		snps = sd.getSnps()
		k = calcKinship(snps)
		if k_file:
			f = open(k_file,'w')
			pickle.dump([k,sd.accessions],f)
			f.close()
		indices_to_keep = []
		for i, acc in enumerate(sd.accessions):
			if acc in accessions:
				indices_to_keep.append(i)		
		k = filter_k_for_accessions(k,indices_to_keep)	
		
	elif snps:
		sys.stdout.write("Generating K.\n")
		sys.stdout.flush()
		if random_fraction:
			import random
			snps = random.sample(snps,int(random_fraction*len(snps)))			 
		k = calcKinship(snps)
	else:
		raise Exception
		
	return k
	

def calcKinship(snps):
	"""
	Requires EMMA to be installed.
	"""
	from rpy import r 
	a = array(snps)
	#r_source(script_dir+"emma_fast.R")
	#r_emma_kinship = robjects.r['emma.kinship']
	#return array(r_emma_kinship(a))
	r.source(script_dir+"emma_fast.R")
	return r.emma_kinship(a)

	
def _run_emma_mp_(in_que,out_que,confirm_que,num_splits=10):
	args = in_que.get(block=True)
	#r_source(script_dir+"emma_fast.R")
	#r_emma_REML_t = robjects.r['emma.REML.t']
	r.source(script_dir+"emma_fast.R")
	phenArray = array([args[2]])
	snps = args[1]
	print "Found data nr. %i, starting run.  pid=%d"%(args[0],os.getpid())
	sys.stdout.flush()
	confirm_que.put([args[0],"start confirmation!"],block=False)
	results = []
	for i in range(num_splits):
		snp_slice = snps[(len(snps)*i)/num_splits:(len(snps)*(i+1))/num_splits]
		snpsArray = array(snp_slice)
	        results.append(r.emma_REML_t(phenArray,snpsArray,args[3],[args[4]]))
		#results.append(r_emma_REML_t(phenArray,snpsArray,args[3],[args[4]]))
		if args[0]==0:
			print "About %d.0%% is done now."%(100*(i+1.0)/num_splits)
			sys.stdout.flush()
	new_res = {}
	for k in results[0]:
		new_res[k]=[]
	for nr in results:
		for k in nr:
			new_res[k].extend(list(nr[k]))	       
	out_que.put([args[0],new_res])

def run_emma_parallel(snps,phen_vals,k,num_proc=8):
	"""
	Run EMMA on multiple processors, using multiprocessing
	"""
	processes = []
	in_que = mp.Queue()
	out_que = mp.Queue()
	confirm_que = mp.Queue()
	import scipy 
	phen_var = scipy.var(phen_vals)
	
	#Populating the in_que
	for i in range(num_proc):
		snp_slice = snps[(len(snps)*i)/num_proc:(len(snps)*(i+1))/num_proc]
		print len(snp_slice)
		in_que.put([i,snp_slice,phen_vals,k,phen_var],block=False)
		p = mp.Process(target=_run_emma_mp_,args=(in_que, out_que, confirm_que))
		p.start()
		done=False
		fail_count=0
		while not done:
			try:
				confirm=confirm_que.get(block=True,timeout=90)
				done=True
			except:
				p.terminate()
				fail_count+= 1
				if fail_count>10:
					print "Giving up!!"
					for p in processes:
						p.terminate()
					raise Exception("Failed to start all processes.")
				else:
					print "Trying again with a new process!"
					in_que.put([i,snp_slice,phen_vals,k,phen_var],block=False)
					p = mp.Process(target=_run_emma_mp_,args=(in_que, out_que, confirm_que))
					p.start()
					
				
		print "ID=%i, Recieved %s"%(confirm[0],confirm[1])
		sys.stdout.flush()
		processes.append(p)
		
	for p in processes:
		print p, p.is_alive()
		
	results = []
	for i in range(num_proc):
		#import pdb;pdb.set_trace()
		if i>0:
			try:
				res = out_que.get(block=True,timeout=10000) #waits about 3 hours, if needed. 
				results.append(res)
			except Exception, err_str:
				print "The parent process didn't receive the results, after waiting over almost 3 hours."
			
		else:
			res = out_que.get()
			results.append(res)
	
	for p in processes:
		p.terminate()
		
	results.sort()
	new_res = {}
	for k in results[0][1]:
		new_res[k]=[]
	for nr in results:
		for k in nr[1]:
			new_res[k].extend(list(nr[1][k]))
		
	#import pdb;pdb.set_trace()
	return new_res
	#print "len(q):",len(q)
	

		
def runEmma(snps,phenValues,k):
	from rpy import r 
	#Assume that the accessions are ordered.
	#r_source("emma.R")
	r.source("emma.R")
	#r_emma_REML_t = robjects.r['emma.REML.t']

	phenArray = array([phenValues])
	snpsArray = array(snps)
	#res = r_emma_REML_t(phenArray,snpsArray,k)
	res = r.emma_REML_t(phenArray,snpsArray,k)
	return res




	


		
def run_fet(snps,phenotypeValues,verbose=False):
	from rpy import r 
	"""
	Fisher's exact test.
	"""
	print "Running Fisher's exact test on",len(snps),"snps, and",len(phenotypeValues),"phenotype values."
	pvals = []
	or_est = []  #Odds ratio estimates
	for snp in snps:
		res = r.fisher_test(phenotypeValues,snp)
		pvals.append(res["p.value"])
		or_est.append(res['estimate'])
	return pvals,or_est

def calc_correlations(snps,phen_vals):
	import scipy as sp
	corrs = []
	for snp in snps:
		corrs.append(sp.corrcoef(snp,phen_vals)[1,0])
	return corrs
	
	
	
if __name__=='__main__':
	_run_()
	print "Done!"
