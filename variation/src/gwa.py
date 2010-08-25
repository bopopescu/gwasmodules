#!/usr/bin/env python2.5
"""
Usage: gwa.py [OPTIONS] [--id=id_text] [PHENOTYPE_ID]

Option:

	--id=...				ID string, used for output files.
	-d ..., --delim=...			default is ", "	  
	-m ..., --missingval=...		default is "NA"
	-h, --help				show this help

	-t ..., --callMethodID=...		What data set is used.
	-r ..., --phen_file=...			Phenotype file, if left out then phenotypes are retireved from the DB (transformed values).
	-k ..., --kinship_file=...		Specify the path to a file containing the kinship matrix.  (otherwise default file is 
						used or it's generated.)

	-a ..., --specific_methods=...		Apply specific methods, otherwise all available are applied:
						lm,emma,emmax,kw,ft, etc.
	-b ..., --specific_transformations=... 	Apply a transformation to the data, default is none, other possibilities are 
						log_trans, box_cox_lambda (where lambda is a number)
	--remove_outliers=...			Should phenotype outliers be removed.  0 (no fence) is the default, else the outlier fence is 
						given in IQRs. (An int is required).
						
	-p ..., --parallel=...			Run mapping methods on the cluster with standard parameters.  The argument is used for runid 
						as well as output files.  Note that if this option is used then no output file should be specified.
						If phenotype index is missing, then all phenotypes are used.
						 
	--use_existing_results			Use existing results when possible, to speed up.  (Otherwise existing files are overwritten)

	--region_plots=...			Include region plots for the top N (given num) peaks.
	--cand_genes_file=...			A file with a list of candidate genes.  (added to plots)
	
	--analysis_plots			Generate analysis plots, histograms, QQ-plots, test different transformations, etc.

	--addToDB				Adds the result(s) to the DB.
	--only_add_2_db				Does not submit jobs, but only adds available result files to DB. (hack for usc hpc)
	
	--data_dir=...				Default is to look for the .gwa_config file in the home folder, and use the 
						specified directory there.
	
	--comment=...				Comment for DB. (Only applicable if result is added to DB.)
	
	--no_phenotype_ids			Phenotypes don't have DB id as an prefix in their names.  (The phenotype index should be used instead.)
	
	--memReq=...				Request memory (on cluster), otherwise use defaults 4GB for Kruskal-Wallis, 12GB for Emma.
	--walltimeReq=...			Request time limit (on cluster), otherwise use defaults
	--proc_per_node=...			Request number of processors per node, default is 8. (Works only for EMMA.)
	--debug_filter=...			Filter SNPs randomly (speed-up for debugging)
	
	--cofactor_chr_pos=...			A list of SNP (chromosome,positions) to be added as co-factors in the analysis.
	--cofactor_phen_id=...			A list of SNP positions to be added as co-factors in the analysis.
	--cofactor_file=...			A file specifying the cofactor.
	--cofactor_no_interact			Exclude interaction terms for cofactors.
						

Examples:
~/gwas_data$ python /home/GMI/bjarni.vilhjalmsson/Projects/py_src/gwa.py --id=../gwas_results/pi1 --specific_methods=kw --callMethodID=54 --phenotype_w_db_id --proc_per_node=8 250K_t54.csv phen_all_052010.tsv	
Description:
  Applies various GWA methods to to phenotypes.

  Methods include: Kruskal-Wallis, Fisher's Exact, EMMA, EMMAX, etc.
  
  If PHENOTYPE_DATA_FILE is left out, then papaya DB is searched and a PHENOTYPE_ID is required.
  
  
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
import time
#import pickle
import cPickle

import linear_models as lm

from numpy import *
#import rpy2.robjects as robjects
#import rpy2.robjects.numpy2ri
#r_source = robjects.r('source')

#import rpy2.rpy_classic as rpy  
#rpy.set_default_mode(rpy.NO_CONVERSION)
#rpy.set_default_mode(rpy.CLASS_CONVERSION)
#Switching to rpy2

from env import *
data_dir = env['data_dir']
script_dir = env['script_dir']
results_dir = env['results_dir']
script_dir = env['script_dir']


transformation_method_dict = {
			'none':1,
			'log_trans':2,
			'box_cox':3, 
			}


analysis_methods_dict = {"kw":1,
			 "ft":2,
			 "emma":4,
			 'lm':16,
			 "emmax":32,
			 }


def prepare_data(sd,phed,p_i,transformation_type,remove_outliers):
	"""
	Coordinates phenotype and snps data for different mapping methods.
	"""
	num_outliers_removed = 0
	if "log_trans"== transformation_type:
		print 'log transforming phenotypes..'
		phed.standard_transform(p_i)
	
	if remove_outliers:
		print 'removing outliers above IQR fence of',remove_outliers,'..'
		num_outliers_removed = phed.naOutliers(p_i,iqrThreshold=remove_outliers)
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
			   "walltimeReq=","specific_methods=",'specific_transformations=',
			   'remove_outliers=',"analysis_plots",'data_dir=',
			   "use_existing_results","region_plots=","cand_genes_file=",
			   "debug_filter=", "no_phenotype_ids", "proc_per_node=", 'phen_file=',
			   'kinship_file=', 'only_add_2_db', 'cofactor_chr_pos', 
			   'cofactor_phen_id', 'cofactor_file', 'cofactor_no_interact' ]
	try:
		opts, args=getopt.getopt(sys.argv[1:], "d:m:h:p:a:b:t:r:k:", long_options_list)

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
	
	memReq = "3800mb"
	walltimeReq = "12:00:00"
	proc_per_node = 8
	
	specific_methods = ['kw','ft','lm','emma','emmax']
	specific_transformations = ['none']
	remove_outliers = 0
	data_dir = env['data_dir']
	
	analysis_plots = False
	use_existing_results = False
	region_plots = 0
	cand_genes_file=None
	debug_filter=1
	no_phenotype_ids = False
	phen_file = None
	kinship_file=None
	only_add_2_db = False
	
	cofactor_chr_pos = None
	cofactor_phen_id = None
	cofactor_file = None
	cofactor_no_interact = False
	
	print ''

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help=1
			print __doc__
		elif opt =="--id":
			run_id=arg
		elif opt in ("-p", "--parallel"):
			parallel=arg
		elif opt =="--addToDB":
			addToDB=True
		elif opt in ('-t', "--callMethodID"):
			callMethodID=int(arg)
		elif opt in ('-r', "--phen_file"):
			phen_file=arg
		elif opt =="--comment":
			comment=arg
		elif opt in ("-d", "--delim"):
			delim=arg
		elif opt in ("-m", "--missingval"):
			missingVal=arg
		elif opt == "--memReq":
			memReq=arg
		elif opt == "--walltimeReq":
			walltimeReq=arg
		elif opt == "--proc_per_node":
			proc_per_node=int(arg)
		elif opt in ('-a', "--specific_methods"):
			specific_methods = arg.split(",")
		elif opt in ('-b', "--specific_transformations"):
			specific_transformations = arg.split(",")
		elif opt in ('-k', "--kinship_file"):
			kinship_file = arg
		elif opt == "--remove_outliers":
			remove_outliers = int(arg)
		elif opt == "--analysis_plots":
			analysis_plots = True
		elif opt == "--data_dir":
			data_dir = arg
		elif opt == "--use_existing_results":
			use_existing_results = True
		elif opt == "--region_plots":
			region_plots = int(arg)
		elif opt == "--cand_genes_file":
			cand_genes_file = arg
		elif opt == "--debug_filter":
			debug_filter = float(arg)
		elif opt == "--no_phenotype_ids":
			no_phenotype_ids = True
		elif opt == "--only_add_2_db":
			only_add_2_db = True
		else:
			if help==0:
				print "Unkown option!!\n"
				print __doc__
			sys.exit(2)

	if len(args)<1 and not parallel:
		if help==0:
			print "Arguments are missing!!\n"
			print __doc__
		sys.exit(2)

	print "GWA runs are being set up with the following parameters:"
	print "run_id:",run_id
	print "phen_file:",phen_file
	print "kinship_file:",kinship_file
	print "parallel:",parallel
	print "walltimeReq:",walltimeReq
	print "proc_per_node:",proc_per_node
	print "memReq:",memReq	
	print "specific_methods:",specific_methods
	print "specific_transformations:",specific_transformations
	print "remove_outliers:",remove_outliers
	print "analysis_plots:",analysis_plots
	print "addToDB:",addToDB
	print "callMethodID:",callMethodID
	print 'phen_file:',phen_file
	print "data_dir:",data_dir
	print "comment:",comment
	print "use_existing_results:",use_existing_results
	print "region_plots:",region_plots
	print "cand_genes_file:",cand_genes_file
	print "debug_filter:",debug_filter
	print 'no_phenotype_ids:',no_phenotype_ids
	print 'only_add_2_db:',only_add_2_db

	def run_parallel(p_i,phed,mapping_method="analysis",trans_method='none'):
		"""
		If no mapping_method, then analysis run is set up.
		"""

		if not p_i in phed.phenIds:
			print "Phenotype ID not found:%i"%p_i
			return 
		if phed.isBinary(p_i):
			if trans_method != 'none':
				return
			elif mapping_method in ["kw"]:
				return
			elif analysis_plots:
				specific_methods = ["ft",'lm',"emma",'emmax']
		else:
			if mapping_method in ["ft"]:
				return
			elif analysis_plots:
				specific_methods = ["kw",'lm',"emma",'emmax']
			
		phenotype_name=phed.getPhenotypeName(p_i)
		#Cluster specific parameters
		print "Setting up a gwa run for phenotype:%s, pid=%d, using method:%s, with transformation as:%s"\
			%(phenotype_name,p_i,mapping_method,trans_method)
		run_id=results_dir+parallel+"_pid"+str(p_i)+"_"+phenotype_name+"_"+mapping_method+"_"+trans_method
		
		#Add results to DB (a hack for papaya, to add results from main node to DB).
		if only_add_2_db:
			pval_file = run_id+".pvals"
			score_file = run_id+".scores"
			sys.stdout.write("Looking for files %s or %s."%(pval_file,score_file))
			result_file = None
			if os.path.isfile(pval_file):
				result_file = pval_file
			elif os.path.isfile(score_file):
				result_file = score_file
			if result_file:
				sys.stdout.write("..... found!\n")
				if no_phenotype_ids:
					db_pid = phed.get_db_pid(p_i)
				else:
					db_pid = p_i
	
				import results_2_db as rdb
				short_name="cm"+str(callMethodID)+"_pid"+str(db_pid)+"_"+phenotype_name+"_"+mapping_method+"_"+trans_method
				if remove_outliers:
					short_name+="_no"
				tm_id = transformation_method_dict[trans_method]
				rdb.add_results_to_db(result_file,short_name,callMethodID,db_pid,analysis_methods_dict[mapping_method], 
						tm_id, remove_outliers=remove_outliers)
				return
			else:
				sys.stdout.write("Result files not found!\n")
				sys.stdout.write("Setting up the run.\n")
				sys.stdout.flush()
		
		
		if remove_outliers:
			run_id+='_no'+str(remove_outliers)
		print "run_id:", run_id

		shstr = "#!/bin/csh\n"
		shstr+="#PBS -l walltime="+walltimeReq+"\n"
		if mapping_method in ['emma']:
			shstr+="#PBS -l nodes=1:ppn=%d \n"%proc_per_node
		shstr+="#PBS -l mem="+memReq+"\n"
		shstr+="#PBS -q cmb\n"
		
		shstr+="#PBS -N p"+str(p_i)+"_"+mapping_method+"_"+parallel+"_"+mapping_method+"_"+trans_method+"\n"
		shstr+="(python "+script_dir+"gwa.py --id="+run_id+"  "
		shstr+="--region_plots="+str(region_plots)+"  "
		shstr+="--proc_per_node=%d  "%proc_per_node
		shstr+="--data_dir=%s  "%data_dir
		if use_existing_results:
			shstr+="--use_existing_results  "
		if analysis_plots:
			shstr+="--analysis_plots  -a "+mapping_method+\
				"  -b "+trans_method+"  "
		else:	
			shstr += " -a "+mapping_method+" "
			shstr += " -b "+trans_method+" "

		if phen_file:
			shstr+="-r "+phen_file+"  "
		if kinship_file:
			shstr+="-k "+kinship_file+"  "			
		if cand_genes_file:
			shstr+="--cand_genes_file="+cand_genes_file+"  "
		if no_phenotype_ids:
			shstr+="--no_phenotype_ids  "
		if remove_outliers:
			shstr+='--remove_outliers='+str(remove_outliers)+'  '
		if debug_filter:
			shstr+="--debug_filter="+str(debug_filter)+"  "	
		if addToDB:
			shstr+="--addToDB "
		if callMethodID:
			shstr+="-t "+str(callMethodID)+"  "
		if comment:
			shstr+="--comment="+str(comment)+"  "
		
		shstr+=str(p_i)+" "  #phenotype ID
		shstr+="> "+run_id+"_job"+".out) >& "+run_id+"_job"+".err\n"
		#print '\n',shstr,'\n'
		f=open(parallel+".sh", 'w')
		f.write(shstr)
		f.close()

		#Execute qsub script
		os.system("qsub "+parallel+".sh ")
		
			
		
		
		

		
	#Load phenotype file
	if phen_file:
		print 'Loading phenotypes from file.'
		phed=phenotypeData.readPhenotypeFile(phen_file,with_db_ids=(not no_phenotype_ids))  #load phenotype file
	else:
		print 'Retrieving the phenotypes from the DB.'
		phed=phenotypeData.getPhenotypes()
	if parallel:
		if len(args)==0:  #phenotype index arguement is missing, hence all phenotypes are run/analyzed.
			if not phen_file:
				raise Exception('Phenotype file or phenotype ID is missing.')
			pids = phed.phenIds
		elif len(args)==1:
			t_pids = args[0].split(',')
			pids = []
			for s in t_pids:
				if '-' in s:
					pid_range = map(int,s.split('-'))
				        for pid in range(pid_range[0],pid_range[1]+1):
				        	pids.append(pid)
				else:
					pids.append(int(s))
		else:
			raise Exception('Too many arguments..')
		
		if analysis_plots:  #Running on the cluster..
			for p_i in pids:
				run_parallel(p_i,phed)
		else:
			for mapping_method in specific_methods:
				for trans_method in specific_transformations:
					for p_i in pids:
						run_parallel(p_i,phed,mapping_method,trans_method)
		return
	else:
		p_i=int(args[0])
	
	filter_accessions = phed.getNonNAEcotypes(p_i)
	phen_is_binary = phed.isBinary(p_i)
	phenotype_name=phed.getPhenotypeName(p_i)
	print "Phenotype:%s, phenotype_id:%s"%(phenotype_name, p_i)

	#SNPs data
	snps_data_file = data_dir+'250K_t'+str(callMethodID)+'.csv'

	if analysis_plots:
		print "\nAnalysing GWAs results jointly... QQ plots etc."
		
		#Genotype and phenotype data is only used for permutations.
		sd=dataParsers.parse_snp_data(snps_data_file , format = 0, delimiter = delim, 
					      missingVal = missingVal, filter = debug_filter,
					      filter_accessions=filter_accessions)
		prepare_data(sd,phed,p_i,trans_method,remove_outliers)
		snps = sd.getSnps()
		phen_vals = phed.getPhenVals(p_i)		

		try:
			print "Plotting accession phenotype map"
			sys.stdout.flush()
			accession_map_pdf_file =  run_id+"_acc_phen_map.pdf"
			phed.plot_accession_map(p_i,pdf_file=accession_map_pdf_file)
		except Exception, err_str:
			print 'Skipping accession-phenotype map... basemap package is probably not installed:',err_str
		
		#Load gwas results
		results = {}
		perm_pvals = None
		methods_found = []
		#Iterate over all possible combination of results
		for mapping_method in specific_methods:
			for trans_method in specific_transformations:
#				for remove_outliers in [True,False]:
				try:
					l = run_id.split("_")[:-1]
					#file_prefix = "_".join(l)+"_"+mapping_method
					file_prefix = "_".join(l)+"_"+mapping_method+'_'+trans_method
					if remove_outliers:
						file_prefix +='_no'
					file_name = file_prefix+".pvals"
					sys.stdout.write("Looking for file %s."%(file_name))
					if os.path.isfile(file_name):
						sys.stdout.write("..... found!\n")
						snps = sd.getSnps()
						res = gwaResults.Result(result_file=file_name,name=mapping_method+"_"+\
								        phenotype_name,snps=snps)
						pvals=True
						if mapping_method in ['lm','emma','emmax']:
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

	else: #If not analysis plots...
		assert len(specific_methods)==1, \
			"More than one specific GWA methods cannot be applied unless running on the cluster."
		mapping_method = specific_methods[0]
		trans_method = specific_transformations[0]
		
		result_name = phenotype_name+'_'+mapping_method+"_"+trans_method #(not really important)
		if remove_outliers:
			result_name += '_no'
		res = None


		
		#Check whether result already exists.
		if use_existing_results:
			if region_plots:
				sd=dataParsers.parse_snp_data(snps_data_file , format = 0, delimiter = delim, 
							      missingVal = missingVal, filter = debug_filter,
							      filter_accessions=filter_accessions)
				num_outliers = prepare_data(sd,phed,p_i,trans_method,remove_outliers)
				if remove_outliers:
					assert num_outliers!=0,"No outliers were removed, so it makes no sense to go on and perform GWA."
		
				snps = sd.getSnps()
			else:
				snps = None

			print "\nChecking for existing results."
			result_file = run_id+".pvals"
			if os.path.isfile(result_file):
				res = gwaResults.Result(result_file=result_file,name=result_name,snps=snps)
				pvals=True
			else:
				result_file = run_id+".scores"
				if os.path.isfile(result_file):
					res = gwaResults.Result(result_file=result_file,name=result_name,snps=snps)
					pvals=False
			if res:
				print "Found existing results.. (%s)"%(result_file)
			sys.stdout.flush()

	

		
		
		if not res: #If results weren't found in a file... then do GWA.
			gc.collect()
			#Do we need to calculate the K-matrix?
			if mapping_method in ['emma','emmax']:
				#Load genotype file (in binary format)
				sys.stdout.write("Retrieving the Kinship matrix K....")
				sys.stdout.flush()
				k_file = data_dir+"kinship_matrix_cm"+str(callMethodID)+".pickled"
				if not kinship_file and os.path.isfile(k_file): #Check if corresponding call_method_file is available
					kinship_file = k_file
				if kinship_file:   #Kinship file was supplied..
					sd=dataParsers.parse_snp_data(snps_data_file , format = 0, delimiter = delim, 
								      missingVal = missingVal, filter = debug_filter,
								      filter_accessions=filter_accessions)
					num_outliers = prepare_data(sd,phed,p_i,trans_method,remove_outliers)	
					k = lm.load_kinship_from_file(kinship_file,sd.accessions)
				else:
					sd=dataParsers.parse_snp_data(snps_data_file , format = 0, delimiter = delim, 
							      missingVal = missingVal, filter = debug_filter)
					print "No kinship file was found.  Generating kinship file:",k_file
					snps = sd.getSnps()		
					k_accessions = sd.accessions[:]
					if debug_filter:
						import random
						snps = random.sample(snps,int(debug_filter*len(snps)))			 
					k = lm.calc_kinship(snps)
					f = open(k_file,'w')
					cPickle.dump([k,sd.accessions],f)
					f.close()
					num_outliers = prepare_data(sd,phed,p_i,trans_method,remove_outliers)
					k = lm.filter_k_for_accessions(k, k_accessions, sd.accessions)
				sys.stdout.flush()
				gc.collect()
				sys.stdout.write("Done!.\n")
			else:
				sd=dataParsers.parse_snp_data(snps_data_file , format = 0, delimiter = delim, 
							      missingVal = missingVal, filter = debug_filter,
							      filter_accessions=filter_accessions)
				num_outliers = prepare_data(sd,phed,p_i,trans_method,remove_outliers)	
			
			snps = sd.getSnps()
			if remove_outliers:
				assert num_outliers!=0,"No outliers were removed, so it makes no sense to go on and perform GWA."
			phen_vals = phed.getPhenVals(p_i)		

			print "Applying %s to data."%(mapping_method)
			sys.stdout.flush()
			kwargs = {}
			additional_columns = []
			if "kw" == mapping_method:
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
				
				if mapping_method in ['emma']:
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
				elif mapping_method in ['emmax']:
					res = lm.emmax(snps,phen_vals,k)
					kwargs['genotype_var_perc'] = res['var_perc']
					kwargs['beta0'] = [val[0] for val in res['betas']]
					kwargs['beta1'] = [val[1] for val in res['betas']]
					additional_columns.append('genotype_var_perc')
					additional_columns.append('beta0')
					additional_columns.append('beta1')
					pvals = res['ps']
					sys.stdout.write("Done!\n")
					sys.stdout.flush()
					
				elif mapping_method in ['py_emma']:
					pass
					#CONNECT TO PYTHON EMMA
				elif mapping_method in ['lm']:
					res = lm.linear_model(snps,phen_vals)
					kwargs['genotype_var_perc'] = res['var_perc']
					kwargs['beta0'] = [val[0] for val in res['betas']]
					kwargs['beta1'] = [val[1] for val in res['betas']]
					additional_columns.append('genotype_var_perc')
					additional_columns.append('beta0')
					additional_columns.append('beta1')
					pvals = res['ps']
					sys.stdout.write("Done!\n")
					sys.stdout.flush()

			
			
			kwargs['correlations'] = calc_correlations(snps, phen_vals)
			additional_columns.append('correlations')
			 
			res = gwaResults.Result(scores=pvals, snps_data=sd, name=result_name, **kwargs)
										
			if pvals:
			 	result_file=run_id+".pvals"
			else:
			 	result_file=run_id+".scores"	
			res.write_to_file(result_file,additional_columns)
		
		#add results to DB..
		if addToDB:
			if no_phenotype_ids:
				db_pid = phed.get_db_pid(p_i)
			else:
				db_pid = p_i

			import results_2_db as rdb
			short_name="cm"+str(callMethodID)+"_pid"+str(db_pid)+"_"+phenotype_name+"_"+mapping_method+"_"+trans_method
			if remove_outliers:
				short_name+="_no"
			#rdb.add_results_to_db(result_file,short_name,callMethodID,db_pid,analysis_methods_dict[mapping_method])
			tm_id = transformation_method_dict[trans_method]
			rdb.add_results_to_db(result_file,short_name,callMethodID,db_pid,analysis_methods_dict[mapping_method], 
						tm_id, remove_outliers=remove_outliers)
			
			

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
		#png_file_max30 = run_id+"_gwa_plot_max30.png"
		if pvals:
			res.neg_log_trans()
			if mapping_method in ["kw","ft"]:
				#res.plot_manhattan(png_file=png_file_max30,percentile=90,type="pvals",ylab="$-$log$_{10}(p)$", 
				#	       plot_bonferroni=True,cand_genes=cand_genes,max_score=30)
				res.plot_manhattan(png_file=png_file,percentile=90,type="pvals",ylab="$-$log$_{10}(p)$", 
					       plot_bonferroni=True,cand_genes=cand_genes)
			else:	
				res.filter_attr("mafs",15)
				#res.plot_manhattan(png_file=png_file_max30,percentile=90,type="pvals",ylab="$-$log$_{10}(p)$", 
				#	       plot_bonferroni=True,cand_genes=cand_genes,max_score=30)
				res.plot_manhattan(png_file=png_file,percentile=90,type="pvals",ylab="$-$log$_{10}(p)$", 
					       plot_bonferroni=True,cand_genes=cand_genes)
		else:
			pass
		
		print "plotting histogram"
		hist_png_file =  run_id+"_hist_"+str(trans_method)+"_no"+str(remove_outliers)+".png"
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
	
	#run function ends here.


	
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
	r.source("emma.R")
	#r_emma_REML_t = robjects.r['emma.REML.t']

	phenArray = array([phenValues])
	snpsArray = array(snps)
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
