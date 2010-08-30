"""
Usage: ThreeLocusTest.py [OPTIONS] [-o OUT_FILE] SNPS_DATA_FILE [PHENOTYPE_INDEX]

Option:

	-o ..., --outputFile=...		Output files, one 'name'.rData file, and one 'name'.pvals file.
	-h, --help				show this help
	-a ..., --mapping_method=...		kw, emmax, lm
	--parallel=...				Run analysis on the cluster with standard parameters.  The arguement is used for runid 
						as well as output files.  Note that if this option is used then no output file should be specified.
	--parallelAll				Apply to all phenotypes.
	--numberPerRun=...			Number of SNPs (phenotypes) per node, default is 100
	--filter=...				Random fraction (given as parameters) of phenotypes will be used.
	--local=...				Do a local GWA with the given windowsize (bases).
	--pvalueThreshold=...			Discard all pvalues larger than the given one.  (Default is no threshold)
	
	--maf_filter=...			Generate the phenotypes on the fly, using SNPs with MAF greater than the given value.
	--latent_variable=...			Type of latent variable: swedish, north_swedish, northern.
	--latent_corr=...			Sample randomly from all real genomic markers which have correlation with 
						the latent allel, greater than the given one.
	--phenotype_model=... 			How to generate the phenotypes: additive, xor, or
	--score_file=...			Score file (e.g. selection scores), from which the top filter % SNPs will be used for the simulation. 			

	--noPvals				Don't write pvalue files.  (saves space and time) DEPRECATED!
	
	--summarizeRuns				Collect results and plot things (Does not generate pvalue files...) 
	--plot_all_models			Plot all phenotype models when summarizing results.
	
	--kinship_file=...			Necessary for running EMMAX.
	--phenotype_error=...			Percentage of variance, due to error.
	--kinship_error=...			Fraction of error, du to kinship term.
	

Examples:
	ThreeLocusTest.py -o outputFile  250K.csv phenotype_index 
	
Description:
	Applies Kruskal-Wallis to fake phenotype...
 
"""
#from epydoc.docparser import str_prefixes
import sys, getopt, traceback
import os, env
import phenotypeData
import dataParsers
import gc
import snpsdata
import plotResults 
import SecondStageAnalysis
import util
from rpy import r
import math
import time
import random
import pdb
import cPickle

#import AddResults2DB

resultDir="/home/cmbpanfs-01/bvilhjal/results/"
scriptDir="/home/cmbpanfs-01/bvilhjal/src/"
tempDir="/home/cmbpanfs-01/bvilhjal/tmp/"

def get_latent_snp(latent_variable,accessions,snpsdata=None):	
	f = open("/home/cmbpanfs-01/bvilhjal/data/eco_dict.pickle")
	ecotype_info_dict = cPickle.load(f)
	f.close()
	#ecotype_info_dict = phenotypeData._getEcotypeIdInfoDict_()
	latent_snp = []
	if latent_variable=="swedish":
		for acc in accessions:
			(native_name,stock_parent,latitude,longitude,country)=ecotype_info_dict[int(acc)]
			if country in ["SWE","FIN"]:
				latent_snp.append(1)
			else:
				latent_snp.append(0)
	elif latent_variable=="northern_swedish":
		for acc in accessions:
			(native_name,stock_parent,latitude,longitude,country)=ecotype_info_dict[int(acc)]
			if country in ["SWE","FIN"] and latitude>=60:
				latent_snp.append(1)
			else:
				latent_snp.append(0)
	elif latent_variable=="northern":
		for acc in accessions:
			(native_name,stock_parent,latitude,longitude,country)=ecotype_info_dict[int(acc)] 
			if latitude>50: #Where is the mean split?
				latent_snp.append(1)
			else:
				latent_snp.append(0)
	elif latent_variable=="random":
		for acc in accessions:
			if random.random()<0.5:
				latent_snp.append(1)
			else:
				latent_snp.append(0)
	elif latent_variable[0:2]=="pc": #Principle component
		pc_num = int(latent_variable[2:])
		pc_file = tempDir+"pc"+str(pc_num)+"_r0.1_pickle.dump"
		import os.path
		if os.path.isfile(pc_file):
			f = open(pc_file,'r')
			pc = cPickle.load(f)
			f.close()
		else:
			print "PC file wasn't found, calculating "+latent_variable+"."
			pc = snpsdata.get_pc()
			f = open(pc_file,'w')
			pc = cPickle.dump(pc,f)
			f.close()			
		
		mean_pc_val = sum(pc)/len(pc)
		for v in pc:
			if v>mean_pc_val:
				latent_snp.append(1)
			else:
				latent_snp.append(0)				
		
	else:
		latent_snp = [0]*len(accessions)
	print latent_snp
	print latent_snp.count(0) 
	return latent_snp


def _run_():
	if len(sys.argv)==1:
		print __doc__
		sys.exit(2)
	
	long_options_list=["outputFile=", "help", "parallel=", "numberPerRun=", "parallelAll", "filter=","local=",
			"pvalueThreshold=", "noPvals", "summarizeRuns", "maf_filter=", "latent_variable=",
			"phenotype_model=", "phenotypeFile=", "runId=", "score_file=","plot_all_models", "latent_corr=",
			'mapping_method=', 'kinship_file=', 'phenotype_error=', 'kinship_error=']
	try:
		opts, args=getopt.getopt(sys.argv[1:], "o:a:h", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
		
	delim=","
	missingVal="NA"
	withArrayIds=1
	outputFile=None
	help=0
	parallel=None
	parallelAll=False
	numberPerRun=100
	filter=None
	local=None
	pvalueThreshold = 1.01
	noPvals = False
	summarizeRuns = False
	maf_filter = 0 #No maf filter
	latent_variable = "northern"
	phenotype_model = None
	phenotypeFile = None
	runId = None
	score_file = None
	plot_all_models = False
	latent_corr=None
	mapping_method = 'kw'
	kinship_file = None
	phenotype_error = 0
	kinship_error = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help=1
			print __doc__
		elif opt in ("-o", "--outputFile"):
			outputFile=arg
		elif opt in ("--parallel"):
			parallel=arg
		elif opt in ("--parallelAll"):
			parallelAll=True
		elif opt in ("--numberPerRun"):
			numberPerRun=int(arg)
		elif opt in ("--pvalueThreshold"):
			pvalueThreshold=float(arg)
		elif opt in ("--filter"):
			filter=float(arg)
		elif opt in ("--local"):
			local=int(arg)
		elif opt in ("--noPvals"):
			noPvals=True
		elif opt in ("--runId"):
			runId=arg
		elif opt in ("--latent_corr"):
			latent_corr=float(arg)
		elif opt in ("--plot_all_models"):
			plot_all_models=True
		elif opt in ("--summarizeRuns"):
			summarizeRuns=True
		elif opt in ("--maf_filter"):
			maf_filter=float(arg)
		elif opt in ("--latent_variable"):
			latent_variable=arg
		elif opt in ("--score_file"):
			score_file=arg
		elif opt in ("--phenotypeFile"):
			phenotypeFile=arg
		elif opt in ("--phenotype_model"):
			phenotype_model=int(arg)
		elif opt in ('-a',"--mapping_method"):
			mapping_method=int(arg)
		elif opt in ("--kinship_file"):
			kinship_file=arg
		elif opt in ("--phenotype_error"):
			phenotype_error=arg
		elif opt in ("--kinship_error"):
			kinship_error=arg
		
		else:
			if help==0:
				print "Unkown option!!\n"
				print __doc__
			sys.exit(2)

	if len(args)<2 and not parallel:
		if help==0:
			print "Arguments are missing!!\n"
			print __doc__
		sys.exit(2)

	
	snpsDataFile=args[0]

	print "Three locus simulations are being set up with the following parameters:"
	print "snpsDataFile:",snpsDataFile
	print "parallel:",parallel
	print "parallelAll:",parallelAll
	print "numberPerRun:",numberPerRun
	print "filter:",filter
	print "local:",local
	print "noPvals:",noPvals
	print "plot_all_models:",noPvals
	print "summarizeRuns:",summarizeRuns
	print "mapping_method:",mapping_method
	print "maf_filter:",maf_filter
	print "latent_variable:",latent_variable 
	print "latent_corr:",latent_corr
	print "phenotype_model:",phenotype_model 
	print "score_file:",score_file
	print "kinship_file:",kinship_file
	print "phenotype_error:",phenotype_error
	print "kinship_error:",kinship_error
	print "runId:",runId

	def runParallel(phen_index,numberPerRun=0,runId=None,phenotypeFile=None,phenotype_model=None):
		#Cluster specific parameters
		if not runId:
			outputFile=resultDir+"TLS_"+parallel+"_"+str(phen_index)+"_"+str(numberPerRun)
		else: #Summarizing runs
			outputFile=resultDir+"TLS_"+parallel+"_"+str(phen_index)
		if phenotype_model and not summarizeRuns:
			outputFile = outputFile+"_pm"+str(phenotype_model)

		shstr="""#!/bin/csh
#PBS -q cmb
#PBS -l walltime=24:00:00
"""
		if runId:
			shstr += "#PBS -l mem=5000m \n"
		else:
			shstr += "#PBS -l mem=3000m \n"
		
		shstr+="#PBS -N TLS_"+str(phen_index)+"_"+parallel+"\n"
		shstr+="(python "+scriptDir+"ThreeLocusTest.py -o "+outputFile+" "
		shstr+=" --numberPerRun="+str(numberPerRun)+" "			
		if local:
			shstr+=" --local="+str(local)+" "			
		if noPvals:
			shstr+=" --noPvals "
		if latent_corr:
			shstr+=" --latent_corr="+str(latent_corr)+" "						
		if plot_all_models:
			shstr+=" --plot_all_models "			
		if maf_filter:
			shstr+=" --maf_filter="+str(maf_filter)+" "		
#		if parallel:
#			shstr+=" --parallel="+str(parallel)+" "		
		if phenotypeFile:
			shstr+=" --phenotypeFile="+str(phenotypeFile)+" "		
		if runId:
			shstr+=" --runId="+str(runId)+" "		
		if score_file:
			shstr+=" --score_file="+str(score_file)+" "		
		if maf_filter:
			shstr+=" --latent_variable="+str(latent_variable)+" "			
		if phenotype_model:
			shstr+=" --phenotype_model="+str(phenotype_model)+" "
		if kinship_file:
			shstr+=" --kinship_file="+str(kinship_file)+" "						
		if summarizeRuns:
			shstr+=" --summarizeRuns "			
		shstr+=" --pvalueThreshold=%f --mapping_method=%s  --phenotype_error=%f  --kinship_error=%f "\
			%(pvalueThreshold,mapping_method,phenotype_error,kinship_error)			
		shstr+=snpsDataFile+" "+str(phen_index)+" "
		shstr+="> "+outputFile+"_job"+".out) >& "+outputFile+"_job"+".err\n"
		#print shstr

		f=open(parallel+".sh", 'w')
		f.write(shstr)
		f.close()

		#Execute qsub script
		os.system("qsub "+parallel+".sh ")

	def summarizeAllRuns(runId, phenotype_model):
		"""
		Plot figures 
		"""
		print "\nWorking on phenotype model:",phenotype_model

		if local:
			sys.stdout.write("Summarizing local results \n")
			total_pos_list = []
			results = []
			r_list = []
			sys.stdout.write("Loading p-values\n")
			sys.stdout.flush()
			fail_count = 0
			phen_index = -numberPerRun
			while True:
				phen_index += numberPerRun
				sys.stdout.write(".")
				sys.stdout.flush()
				filename=resultDir+"TLS_"+runId+"_"+str(phen_index)+"_"+str(numberPerRun)+".pvals"
				try:
					f = open(filename,"r")
					l = cPickle.load(f)
					f.close()
					fail_count = 0
				except Exception, err_str:
					print "Had problems with file",filename
					print err_str
					fail_count += 1
					if fail_count >20:
						print "More than 20 consecutive file reads failed... Giving up!"
						break 
					else:
						continue
				total_pos_list += l[0]	
				results += l[1]	
				
			sys.stdout.write("\nP-values loaded\n")
			sys.stdout.write("Sorting into bins\n")
			sys.stdout.flush()
			max_dist = 100000
			n_bins = 100
			bins = [[] for i in range(0,n_bins)]
			other_bin = []
			for i in range(0,len(total_pos_list)):
				if i%1000==0:
					sys.stdout.write(".")
					sys.stdout.flush()
				chr_pos_pvals = results[i]
				for (chr,pos,pval) in chr_pos_pvals:
					
					d = abs(int(total_pos_list[i])-int(pos))
					if chr!=1 or d>=max_dist:
						other_bin.append(-math.log10(pval))
					else:
						bin_nr = (n_bins*d)/max_dist
						bins[bin_nr].append(-math.log10(pval))
			bin_pos = range(max_dist/(2*n_bins),max_dist+max_dist/(2*n_bins),max_dist/n_bins)			
			
			sys.stdout.write("\nRertieveing quartiles\n")
			sys.stdout.flush()
			min_pvals = []
			pvals_q1s = []
			pvals_q2s = []
			pvals_q3s = []
			pvals_q90s = []
			pvals_q99s = []
			
			for i in range(0,len(bins)):
				bin = bins[i]
				print i, len(bin)
				
				if len(bin)>0:
					bin.sort()
					min_pvals.append(bin[-1])
					pvals_q1s.append(bin[len(bin)/4])
					pvals_q2s.append(bin[len(bin)/2])
					pvals_q3s.append(bin[(3*len(bin))/4])
					pvals_q90s.append(bin[int(0.9*len(bin))])
					pvals_q99s.append(bin[int(0.99*len(bin))])
				else:
					min_pvals.append(0)
					pvals_q1s.append(0)
					pvals_q2s.append(0)
					pvals_q3s.append(0)
					pvals_q90s.append(0)
					pvals_q99s.append(0)
					

			sys.stdout.write("Plotting data\n")
			sys.stdout.flush()
			#Plot figure..
			plt.figure(figsize=(10,7))
			plt.axes([0.14,0.13,0.81,0.83])
			#plt.plot(bin_pos,min_pvals,label="Min. p-val")
			plt.plot(bin_pos,pvals_q1s,label="25% quantile")
			plt.plot(bin_pos,pvals_q2s,label="50% quantile")
			plt.plot(bin_pos,pvals_q3s,label="75% quantile")
			plt.plot(bin_pos,pvals_q90s,label="90% quantile")
			plt.plot(bin_pos,pvals_q99s,label="99% quantile")
			plt.ylabel("$-log_{10}(p)$")
			plt.xlabel("Distance from causative SNP (bases)")
			plt.legend()
			pdfFile = resultDir+"TLS_"+runId+".dist2.pdf"
			plt.savefig(pdfFile, format = "pdf")
			pngFile = resultDir+"TLS_"+runId+".dist2.png"
			plt.savefig(pngFile, format = "png")
			plt.clf()
			
						
		else: # then global, one model at a time.
			total_pos_list = []
			d_list = []
			r_list = []
			fail_count = 0
			phen_index = -numberPerRun
			loaded_phen_mafs = []
			max_pval_list = []
			obs_pval_list = []
			latent_pvals  = []
			latent_snps_ranks = []
			latent_distances_to_min_pval = []
			latent_loci_snp_chr_pos_mafs = []
			rank_statistics = []
			sign_fractions = []
			sign_statistics = []
			
				
			while True: #phen_index<len(phen_mafs)-numberPerRun:
				phen_index += numberPerRun
				sys.stdout.write(".")
				sys.stdout.flush()
				filename=resultDir+"TLS_"+runId+"_"+str(phen_index)+"_"+str(numberPerRun)+"_"+"pm"+str(phenotype_model)+".stats"
				print "Loading data from file:",filename
				try:
					f = open(filename,"r")
					p_d_r_stats = cPickle.load(f)
					f.close()
					fail_count = 0
				except Exception, err_str:
					#print "Had problems with file",filename
					#print err_str
					fail_count += 1
					if fail_count >20:
						print "More than 20 consecutive file reads failed... Giving up!"
						print "Had problems with file",filename
						print err_str
						break 
					else:
						continue

				try:
					total_pos_list += p_d_r_stats['phen_positions']
					d_list += p_d_r_stats['distances_to_min_pval']	
					r_list += p_d_r_stats['true_snps_ranks']				
					loaded_phen_mafs += p_d_r_stats['phen_mafs']
					obs_pval_list += p_d_r_stats['obs_pvals']
					max_pval_list += p_d_r_stats['min_pvals']
					sign_fractions += p_d_r_stats.get('sign_fractions')
					sign_statistics += p_d_r_stats.get('sign_statistics')
					
					latent_pvals += p_d_r_stats.get('latent_pvals',[])
					latent_snps_ranks += p_d_r_stats.get('latent_snps_ranks',[])
					latent_distances_to_min_pval += p_d_r_stats.get('latent_distances_to_min_pval',[])
					latent_loci_snp_chr_pos_mafs += p_d_r_stats.get('latent_loci_snp_chr_pos_mafs',[])
					rank_statistics += p_d_r_stats.get('rank_statistics',[])
					
						
					#print max_pval_list
					#print "Found",len(max_pval_list),"pvals."
				except Exception, err_str:
					print "Failed when loading dictionary..." 
					print err_str
					fail_count += 1
					if fail_count >20:
						print "More than 20 consecutive file reads failed... Giving up!"
						break 
			
			print "len(r_list), len(d_list), len(loaded_phen_mafs):",len(r_list), len(d_list), len(loaded_phen_mafs)
			print "len(sign_statistics):",len(sign_statistics)
			print "len(obs_pval_list), len(latent_pvals):", len(obs_pval_list), len(latent_pvals)
			print "Mean significant fraction:",sum(sign_fractions)/float(len(sign_fractions))


			
		
			#Setting figure legends, etc.
			show_y_label=False
			show_x_label=False
			show_legend=False
			if phenotype_model==1:
				fig_label='C'
				fig_label2='F'
				show_legend=True
			elif phenotype_model==2:
				fig_label='B'
				fig_label2='E'
				show_x_label=True
			elif phenotype_model==3:
				fig_label='A'
				fig_label2='D'
				show_y_label=True
			elif phenotype_model==4:
				fig_label=''
				fig_label2=''
				
				
			runId = runId+"_pm"+str(phenotype_model)
			filename = resultDir+"TLS_"+runId
			plot_max_dist_hist(sign_statistics,obs_pval_list,latent_pvals,filename,fig_label=fig_label,
					show_y_label=show_y_label,show_x_label=show_x_label,show_legend=show_legend)

			plot_rank_dist_hist(rank_statistics,obs_pval_list,latent_pvals,filename,fig_label=fig_label2,
					show_y_label=show_y_label,show_x_label=show_x_label,show_legend=show_legend)

			plot_dist_hist(d_list,loaded_phen_mafs,max_pval_list,filename,fig_label=fig_label2,
					show_y_label=show_y_label,show_x_label=show_x_label,show_legend=show_legend)
			plot_rank_hist(r_list,loaded_phen_mafs,max_pval_list,filename,fig_label=fig_label,
					show_y_label=show_y_label,show_x_label=show_x_label,show_legend=show_legend)
			if latent_pvals:
				filename = resultDir+"TLS_"+runId+"_latent"
				maf_list = map(list,zip(*latent_loci_snp_chr_pos_mafs))[-1]
				plot_dist_hist(latent_distances_to_min_pval,maf_list,max_pval_list,
						filename,fig_label=fig_label2,show_y_label=show_y_label,
						show_x_label=show_x_label,show_legend=show_legend)
				plot_rank_hist(latent_snps_ranks,maf_list,max_pval_list,
						filename,fig_label=fig_label,show_y_label=show_y_label,
						show_x_label=show_x_label,show_legend=show_legend)

				plot_rank_3D_hist(r_list,latent_snps_ranks,max_pval_list,filename)

				first_ranks = []
				first_mafs = []
				second_ranks = []
				second_mafs = []

				rank_origins = []
				dist_origins = []
				#Minimum distances, to either of the two causative loci.
				min_distances = []
				min_ranks = []
				min_rank_mafs = []
				min_dist_mafs = []
				for d1,d2,r1,r2,m1,m2 in zip(d_list,latent_distances_to_min_pval,r_list,latent_snps_ranks,loaded_phen_mafs,maf_list):
					ml = [m1,m2]
					if d2<0: 
						min_distances.append(d1)
						min_dist_mafs.append(m1)
						dist_origins.append(1)
					elif d1 <0:
						min_distances.append(d2)
						min_dist_mafs.append(m2)
						dist_origins.append(2)
					else:
						dl = [d1,d2]
						if d1<d2:
							dist_origins.append(1)						
						else:
							dist_origins.append(2)						
						min_d = min(dl)
						min_dist_mafs.append(ml[dl.index(min_d)])
						min_distances.append(min_d)
					rl = [r1,r2]
					if r1<r2:
						rank_origins.append(1)
						first_ranks.append(r1)
						first_mafs.append(m1)						
						second_ranks.append(r2)
						second_mafs.append(m2)
					else:
						rank_origins.append(2)						
						first_ranks.append(r2)
						first_mafs.append(m2)						
						second_ranks.append(r1)
						second_mafs.append(m1)
					min_r = min(rl)
					min_rank_mafs.append(ml[rl.index(min_r)])
					min_ranks.append(min_r)
				filename = resultDir+"TLS_"+runId+"_shortest"
				maf_list = map(list,zip(*latent_loci_snp_chr_pos_mafs))[-1]
				plot_dist_hist(min_distances,min_dist_mafs,max_pval_list,
						filename,fig_label=fig_label2,show_y_label=show_y_label,show_x_label=show_x_label,
						show_legend=show_legend,dist_origins=dist_origins)
				plot_rank_hist(min_ranks,min_rank_mafs,max_pval_list,
						filename,fig_label=fig_label,show_y_label=show_y_label,show_x_label=show_x_label,
						show_legend=show_legend,rank_origins=rank_origins)
					
				filename = resultDir+"TLS_"+runId+"_both_ranks"
				plot_both_ranks_hist(first_ranks,second_ranks,first_mafs,second_mafs,max_pval_list,
						filename,fig_label=fig_label,show_y_label=show_y_label,show_x_label=show_x_label,
						show_legend=show_legend,rank_origins=rank_origins)
	#End of summarizeAllRuns function.

	
	if phenotypeFile:
		phenotype_file = phenotypeFile
	else:
		print "Generating phenotype"
		phenotype_file = tempDir+parallel+".phen"

	if phenotype_model:  #If phenotype model is specified, then use that, otherwise run all..
		phenotype_models = [phenotype_model]
	else:
		phenotype_models = range(1,5)

		
	if parallelAll:
		#Generating the phenotype
#		if score_file:
#			snps_dataset = dataParsers.parse_snp_data(snpsDataFile,format=0)
#			#1. Load the score file...
#			(positions,phs_scores) = _read_phs_result_()
#			#2. Sort it by scores
#			score_pos_list = []
#			for i in range(0,5):
#				score_pos_list.extend(zip(phs_scores[i],positions[i],[i+1]*len(positions[i])))
#				score_pos_list.sort(reversed=True)
#			#3. choose the top %filter SNPs..
#			score_pos_list = score_pos_list[0:int(filter*len(score_pos_list))]
#			snps_to_keep_indices = []
#			for (score,pos,ch) in score_pos_list:
#				snps_to_keep_indices.append(snps_dataset)  #FINISH!!!


		snps_dataset = dataParsers.parse_binary_snp_data(snpsDataFile)
		if 0<maf_filter<=0.5:
			snps_dataset.filter_maf_snps(maf_filter)
		if latent_variable!='random_snp':
			latent_snp = get_latent_snp(latent_variable,snps_dataset.accessions,snps_dataset)
			if latent_corr:
				latent_snp_chr_pos_maf = snps_dataset.get_top_correlated_snp(latent_snp, latent_corr)
				print len(latent_snp_chr_pos_maf),"SNPs were found to have higher r^2 than",latent_corr,"with the latent loci."
		if latent_variable=='random_snp':
			latent_snp_chr_pos_maf = snps_dataset.get_all_snp_w_info()
			latent_corr = 1
		if filter:
			snps_dataset.sample_snps(filter)
		snps_list = snps_dataset.getSnps()
		anti_decoder = {1:0,0:1}
		chr_pos_list = snps_dataset.getChrPosList()
		mafs = snps_dataset.get_mafs()["marfs"]

		phen_dict = {}
		for phenotype_model in phenotype_models:
			phenotypes = []
			phen_positions = []
			phen_chr_pos = []
			phen_mafs = []
			latent_loci_snp_chr_pos_mafs = []
			for i, snp in enumerate(snps_list):
				(chr,pos) = chr_pos_list[i]
				maf = mafs[i]
				anti_snp = []
				for nt in snp:
					anti_snp.append(anti_decoder[nt]) 
				phenotype = []
				anti_phenotype = []
				if latent_corr:
					lsd = random.choice(latent_snp_chr_pos_maf)
					(latent_snp,latent_chr,latent_pos,latent_maf) = lsd
					while latent_snp==snp: #Make sure the two SNPs aren't identical.
						lsd = random.choice(latent_snp_chr_pos_maf)
						(latent_snp,latent_chr,latent_pos,latent_maf) = lsd
				if latent_variable=="random":
					latent_snp = []
					for acc in snps_dataset.accessions:
						if random.random()<0.5:
							latent_snp.append(1)
						else:
							latent_snp.append(0)								
				for j in range(len(snp)):
					if phenotype_model == 1:#xor
						#print "XOR"
						phen_val = int(snp[j]!=latent_snp[j])
						anti_phen_val = int(anti_snp[j]!=latent_snp[j])
					elif phenotype_model == 2:#or
						#print "OR/AND"
						phen_val = snp[j] or latent_snp[j] 
						anti_phen_val = anti_snp[j] or latent_snp[j] 
					elif phenotype_model == 3:#plus
						#print "+"
						phen_val = snp[j] + latent_snp[j] 
						anti_phen_val = anti_snp[j] + latent_snp[j] 
					elif phenotype_model == 4:#xor plus
						#print "XOR +0.5"
						phen_val = int(snp[j]!=latent_snp[j])+0.5*int(snp[j]==1 and latent_snp[j]==1)
						anti_phen_val = int(anti_snp[j]!=latent_snp[j])+0.5*int(anti_snp[j]==1 and latent_snp[j]==1)
					phenotype.append(phen_val)
					anti_phenotype.append(anti_phen_val)
				#Check whether phenotype is OK.
				#SHOULD THIS CHECK BE SKIPPED???
				if len(set(phenotype))>1:
					phenotypes.append(phenotype)
					#print phenotype
					phen_positions.append(pos)
					phen_chr_pos.append((chr,pos))
					phen_mafs.append(maf)
					if latent_corr:
						latent_loci_snp_chr_pos_mafs.append(lsd)
	
				else:
					print "Found problematic phenotype"
				if len(set(anti_phenotype))>1:
					phenotypes.append(anti_phenotype)		
					#print anti_phenotype
					phen_positions.append(pos)
					phen_chr_pos.append((chr,pos))		
					phen_mafs.append(maf)
					if latent_corr:
						latent_loci_snp_chr_pos_mafs.append(lsd)
				else:
					print "Found problematic anti-phenotype"
				if anti_phenotype == phenotype:
					print "Phenotype and anti-phenotype are the same!?!?!?"
				
			print "Phenotypes generated for phenotype model:",phenotype_model
			if latent_corr:
				d = {"phenotypes":phenotypes,"phen_positions":phen_positions,"phen_chr_pos":phen_chr_pos,
				"latent_snp":latent_snp, "phen_mafs":phen_mafs, 
				"latent_loci_snp_chr_pos_mafs":latent_loci_snp_chr_pos_mafs}
			else:
				d = {"phenotypes":phenotypes,"phen_positions":phen_positions,"phen_chr_pos":phen_chr_pos,
				"latent_snp":latent_snp, "phen_mafs":phen_mafs, "latent_loci_snp_chr_pos_mafs":[]}
			phen_dict[phenotype_model]=d
			#phenotype_models for loop ends.
		f = open(phenotype_file,"w")
		print "dumping phenotypes to file:",f
		cPickle.dump(phen_dict,f)
		f.close()
		
		
		
	elif not (summarizeRuns and parallel):	
		f = open(phenotype_file,"r")

		phen_dict = cPickle.load(f)
		print f
		print phen_dict.keys
		print phenotype_models[0]
		print phenotype_models
		d = phen_dict[phenotype_models[0]]
		phenotypes = d["phenotypes"]
		phen_positions = d["phen_positions"]
		phen_chr_pos = d["phen_chr_pos"]
		latent_snp = d["latent_snp"]
		phen_mafs = d["phen_mafs"]
		latent_loci_snp_chr_pos_mafs = d["latent_loci_snp_chr_pos_mafs"]
		f.close()
	else: #Submitting a summary run to the cluster!
		runParallel("summary",numberPerRun, runId=parallel, phenotypeFile=phenotype_file)
		return

			
	num_phenotypes = len(phenotypes) 
	print "num_phenotypes:",num_phenotypes
		
	
	if parallelAll and not summarizeRuns:  #Running on the cluster..
		phen_indices = range(0,num_phenotypes,numberPerRun)
		for phenotype_model in phenotype_models:
			print "Submitting jobs for phenotype model:",phenotype_model
			for phen_index in phen_indices:
				runParallel(phen_index,numberPerRun,phenotypeFile=phenotype_file,phenotype_model=phenotype_model)
		return

	elif summarizeRuns:		
		#summarizeAllRuns(phenotypes,phen_positions,phen_chr_pos,phen_mafs,latent_snp)
		
		for phenotype_model in phenotype_models:
			summarizeAllRuns(runId, phenotype_model)
		return
	else:
		phen_index=int(args[1])
	print "phen_index:",phen_index
	print "\nStarting simulation now!\n"
	print "output:",outputFile
	
	#Filter phenotypes for actual chosen SNPs:
	sys.stdout.write("Filtering phenotypes for this run.\n")
        sys.stdout.flush()
	new_phenotypes = []       
	lim = min(phen_index+numberPerRun,num_phenotypes)
	phen_positions = phen_positions[phen_index:lim]
	phen_chr_pos = phen_chr_pos[phen_index:lim]
	phen_mafs = phen_mafs[phen_index:lim]
	phenotypes = phenotypes[phen_index:lim]
	latent_loci_snp_chr_pos_mafs = latent_loci_snp_chr_pos_mafs[phen_index:lim]


	print "Loading SNPS dataset (again)" 
	snps_dataset = dataParsers.parse_snp_data(snpsDataFile,format=0)
	snps_positions = snps_dataset.getPositions()
	snps_list = snps_dataset.getSnps()
	chr_pos_list = snps_dataset.getChrPosList()
	#Run KW
	results = [] #[num_of_phen][num_of_snps]
	print "Running KW"
	if local:
		
		print "len(snps_list):",len(snps_list)
		def _get_local_snps_(phen_pos,local):
			snps = []
			curr_pos = snps_positions[0]
			i = 0
			while curr_pos < (phen_pos-local/2.0) and i < len(snps_positions)-1:
				i += 1
				curr_pos = snps_positions[i]
			prefix_num = i
			while curr_pos < (phen_pos+local/2.0) and i < len(snps_positions)-1:
				snps.append(snps_list[i])
				i += 1
				curr_pos = snps_positions[i]
			suffix_num = len(snps_positions)-i
			print len(snps),prefix_num,suffix_num
			return (snps,prefix_num,suffix_num)

		for i in range(0,len(phenotypes)):
			#Find local snps & fix results
			(snps,prefix_num,suffix_num) = _get_local_snps_(phen_positions[i],local)
			phenotype = phenotypes[i]
			print "KW is being applied to the "+str(i)+"'th phenotype (locally)."
			pvals = [2]*prefix_num
			pvals += util.kruskal_wallis(snps,phenotype)["ps"]
			pvals += [2]*suffix_num
			results.append(pvals)

	else: #Global run
		print "len(chr_pos_list):",len(chr_pos_list)
		print "len(snps_list):",len(snps_list)
		for i in range(0,len(phenotypes)):
			phenotype = phenotypes[i]
       			if mapping_method=='kw':
				sys.stdout.write("KW is being applied to the "+str(i)+"'th phenotype.\n")
	       			sys.stdout.flush()
			       	pvals = util.kruskal_wallis(snps_list,phenotype)["ps"]
			elif mapping_method=='emmax':
				pass
			elif mapping_method=='lm':
				pass
			results.append(pvals)
	print "len(results):",len(results)

	# ------------------------------------- SUMMARIZING RESULTS!!!!!! -------------------------------------
	#Calculating various statistics.
	#Distances to min p-val.
	distances_to_min_pval = [] 
	#Ranks of true SNPs
	true_snps_ranks = []
	min_pvals = []
	obs_pvals = []
	rank_statistics = []  #For the max ranks histogram!
	sign_statistics = []  #For the max significant histogram!
	sign_fractions = []
	if latent_corr:
		latent_distances_to_min_pval = [] 
		latent_snps_ranks = []
		latent_pvals = []
		second_dist_to_min_pval = []  
	print len(results),len(phen_mafs)
	for i in range(len(results)):
		sys.stdout.write("Summarizing results for the "+str(i)+"'th phenotype.\n")
		sys.stdout.flush()
		pvals = results[i]
		sign_count = 0
		tot_count = 0
		pval_chr_pos = []
		sign_pval_chr_pos = []	
		for j, pval in enumerate(pvals):
			tot_count += 1
			if math.log10(pval)<-6.6354837468149119:
				sign_count +=1
				sign_pval_chr_pos.append((pval,chr_pos_list[j][0],chr_pos_list[j][1]))
			pval_chr_pos.append((pval,chr_pos_list[j][0],chr_pos_list[j][1]))
		sign_fractions.append(sign_count/float(tot_count))
		pval_chr_pos.sort()
		
		min_pval = min(pvals)
		min_pvals.append(min_pval)
		j = pvals.index(min_pval)
		
		#chromosome and position of the causative SNP
		(phen_chr,phen_pos) = phen_chr_pos[i]
		
		#chromosome and position of the most significant SNP
		(chr,pos) = chr_pos_list[j]  
		
		if chr==phen_chr:
			distances_to_min_pval.append(abs(int(phen_pos)-int(pos)))
		else:
			distances_to_min_pval.append(-abs(int(phen_chr)-int(chr)))
		try:
			snp_i = chr_pos_list.index((phen_chr,phen_pos))
			obs_pval = pvals[snp_i]
			obs_pvals.append(obs_pval)
			r_list = util.getRanks(pvals)
			snp_rank = r_list[snp_i]
			true_snps_ranks.append(snp_rank)
		except Exception, err_str:
			print "Didn't find causative SNP:",phen_pos
			print err_str
			true_snps_ranks.append(-1)

		if latent_corr:
			(latent_snp,latent_chr,latent_pos,latent_maf) = latent_loci_snp_chr_pos_mafs[i]
			if chr==latent_chr:
				latent_distances_to_min_pval.append(abs(int(latent_pos)-int(pos)))
			else:
				latent_distances_to_min_pval.append(-abs(int(latent_chr)-int(chr)))
			try:
				latent_snp_i = chr_pos_list.index((latent_chr,latent_pos))
				obs_pval = pvals[latent_snp_i]
				latent_pvals.append(obs_pval)
				latent_rank = r_list[latent_snp_i]
				latent_snps_ranks.append(latent_rank)
			except Exception, err_str:
				print "Didn't find latent causative SNP:",phen_pos
				print err_str
				true_snps_ranks.append(-1)

			# Which is the closest one to the min. pval?
			if latent_distances_to_min_pval[-1]<0 and distances_to_min_pval[-1]<0:
				second_dist_to_min_pval.append(-3)  #Both causative loci on diff. chromosomes from the min pval.
			elif latent_distances_to_min_pval[-1]<0:
				second_chr=latent_chr
				second_pos=latent_pos				
			elif distances_to_min_pval[-1]<0:
				second_chr=phen_chr
				second_pos=phen_pos				
			else:
				second_dist_to_min_pval.append(-2)  #Both causative loci on same chromosome.
							
			# New min pval!!
#			elif chr != first_chr:
#				if chr != second_chr:
#					second_dist_to_min_pval.append(-1)
#				else:
#					second_dist_to_min_pval.append(abs(int(second_pos)-pos))
#			else:
#				
				


		#Calculating the dist to the farthest SNP with rank greater or equal to the second ranked causative SNP.
		if latent_corr:
			print latent_corr
			max_rank = int(math.ceil(max(latent_rank,snp_rank))+0.001)
		else:
			max_rank = int(math.ceil(snp_rank)+0.001)
			  
		im_pval_chr_pos = pval_chr_pos[:max_rank]
		max_dist = 0
		diff_chr = False
		for (im_pval,im_ch,im_pos) in im_pval_chr_pos:
			if latent_corr and im_ch == phen_chr and im_ch == latent_chr:
				max_dist = max(max_dist, min(abs(im_pos-phen_pos),abs(im_pos-latent_pos)))
			elif im_ch == phen_chr:
				max_dist = max(max_dist, abs(im_pos-phen_pos))
			elif latent_corr and im_ch == latent_chr:
				max_dist = max(max_dist, abs(im_pos-latent_pos))
			else:
				diff_chr = True
				break
		if diff_chr:
			rank_statistics.append(-1)
		else:
			rank_statistics.append(max_dist)
		
		#Calculating the max distance to a significant p-val.
		if sign_pval_chr_pos:
			max_dist = 0
			diff_chr = False
			for (s_pval,s_ch,s_pos) in sign_pval_chr_pos:
				if s_ch == phen_chr:
					if latent_corr and phen_chr == latent_chr:
						max_dist = max(max_dist, min(abs(s_pos-phen_pos),abs(phen_pos-latent_pos)))
					else:
						max_dist = max(max_dist, abs(s_pos-phen_pos))
				elif latent_corr and s_ch == latent_chr:
					max_dist = max(max_dist, abs(phen_pos-latent_pos))
				else:
					diff_chr = True
					break
			if diff_chr:
				sign_statistics.append(-1)
			else:
				sign_statistics.append(max_dist)
#					if max_dist == 0:
#						print "sign_pval_chr_pos:",sign_pval_chr_pos
#						print "pvals[snp_i]:",pvals[snp_i]
#						print "pvals[latent_snp_i]:",pvals[latent_snp_i]
#						print "phen_pos,phen_chr:",phen_pos,phen_chr
#						print "latent_pos,latent_chr:",latent_pos,latent_chr
		else:
			sign_statistics.append(-2)

		
				
			 		
			
			
			
		gc.collect()  #Calling garbage collector, in an attempt to clean up memory..
	print "Distances:", distances_to_min_pval
	print "Ranks:", true_snps_ranks
	print "Obs. pvals:", obs_pvals
	print "Min pvals:", min_pvals
	print "Significant fractions:",sign_fractions
	print "Significant pvals statistics:",sign_statistics
	if latent_corr:
		print "Distances:", latent_distances_to_min_pval
		print "Latent ranks:", latent_snps_ranks
		print "Latent obs. pvals:", latent_pvals
		print "Intermittent rank statistics:", rank_statistics

	
	#Discarding the missing p-values due to local association mapping. 
	#results = map(list,zip(*results))
	#print "len(results):",len(results), len(chr_pos)

	#Write results to file

	if not noPvals and local: 
		print "Writing p-values to file:"
		new_results = []
		for pvals in results: #For all phenotypes (results)
			chr_pos_pvals = []
			for i in range(0,len(chr_pos_list)): #For all SNPs
				pval = pvals[i]
				if pval!=2 and pval<=pvalueThreshold:
					(chr,pos)=chr_pos_list[i]
					chr_pos_pvals.append((chr,pos,pval))
			new_results.append(chr_pos_pvals)
		results = new_results	
		l = [phen_positions,results]
		pvalFile = outputFile+".pvals"
		print "Writing p-values to file:",pvalFile
		f = open(pvalFile,"w")
		cPickle.dump(l,f)	
		f.close()

	p_d_r_dict = {}
	if latent_corr:
		p_d_r_dict['latent_pvals'] = latent_pvals
		p_d_r_dict['latent_snps_ranks'] = latent_snps_ranks
		p_d_r_dict['latent_distances_to_min_pval'] = latent_distances_to_min_pval
		p_d_r_dict['latent_loci_snp_chr_pos_mafs'] = latent_loci_snp_chr_pos_mafs

	p_d_r_dict['phen_positions'] = phen_positions
	p_d_r_dict['distances_to_min_pval'] = distances_to_min_pval
	p_d_r_dict['true_snps_ranks'] = true_snps_ranks
	p_d_r_dict['phen_mafs'] = phen_mafs
	p_d_r_dict['obs_pvals'] = obs_pvals
	p_d_r_dict['min_pvals'] = min_pvals
	p_d_r_dict['sign_fractions'] = sign_fractions
	p_d_r_dict['rank_statistics'] = rank_statistics
	p_d_r_dict['sign_statistics'] = sign_statistics
	
	filename = outputFile+".stats"	
	print "Writing results to file:",filename
	f = open(filename,"w")
	cPickle.dump(p_d_r_dict,f)
	f.close



def plot_max_dist_hist(sign_statistics,snp_pvals,latent_pvals,filename,hist_color='#CCFF00',fig_label='',
		show_y_label=True,show_x_label=True,show_legend=True):
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt

#	plt.figure(figsize=(3.5,4.2))
#	plt.axes([0.17,0.28,0.79,0.69])

	group_names =["d$=$0","0$<$d$\leq$5", "5$<$d$\leq$50", "50$<$d", "Other chr."]
	if latent_pvals:
		sub_group_names = ["Neither significant","One significant","Both significant"]
	else:
		sub_group_names = ["Caus. not significant","Caus. significant"]
	group_counts = [[0]*len(sub_group_names) for i in range(len(group_names))]
	d_bin_thresholds = [0,1,5001,50001]
	valid_count = 0
	for i, smd in enumerate(sign_statistics):
		if smd!= -2:
			valid_count += 1
			if latent_pvals:
				sign_group = int(math.log(snp_pvals[i])<=-6.6354837468149119)+int(math.log(latent_pvals[i])<=-6.6354837468149119)
			else:
				sign_group = int(math.log(snp_pvals[i])<=-6.6354837468149119)
			if smd <0:
				group_counts[-1][sign_group]+=1
			elif smd>d_bin_thresholds[-1]:
				group_counts[-2][sign_group]+=1
			else:
				for i in range(len(d_bin_thresholds)-1):
					dthres1 = d_bin_thresholds[i]
					dthres2 = d_bin_thresholds[i+1]
					if dthres1-0.1<smd<dthres2:
						group_counts[i][sign_group]+=1
	
	group_counts = [[gc2/float(valid_count) for gc2 in gc1] for gc1 in group_counts]
	
	group_colors = ["#9999FF","#CCFF00","#33CC00"]#,"#3399FF","#FF9900"]
		
	
	
	
	plt.figure(figsize=(3.5,4.2))
	plt.axes([0.17,0.28,0.79,0.69])
	#Plotting the histogram!! (manually)
	for i, gc in enumerate(group_counts):
		for j in range(len(sub_group_names)):
			if i==len(group_counts)-1:			
				plt.bar(i,sum(gc[j:]),color=group_colors[j],label=sub_group_names[j])
			else:
				plt.bar(i,sum(gc[j:]),color=group_colors[j])

	minVal = 0
	maxVal = len(group_counts)-0.2
	x_range = maxVal - minVal
	y_max = max(group_counts)	 
	print group_counts,y_max
	plt.axis([minVal-0.035*x_range,maxVal+0.035*x_range,-0.035*1.18,1.035*1.18])
	
	x_ticks_pos = []
	for i in range(len(group_names)):
		x_ticks_pos.append(i+0.4)
	plt.xticks(x_ticks_pos,group_names,fontsize="medium",rotation=45,horizontalalignment='right')
	if show_y_label:
		plt.ylabel("Frequency",fontsize='large')
	plt.text(0.2,1.06,fig_label,fontsize=15)
	if show_x_label:
		plt.xlabel("Distance to nearest causative (kb)",fontsize='medium')
	if show_legend:
		fontProp = matplotlib.font_manager.FontProperties(size=11)
		if latent_pvals:
			plt.legend(loc = 9, numpoints = 1, handlelen = 0.04, markerscale = 1,prop=fontProp)
		else:
			plt.legend(loc = 1, numpoints = 1, handlelen = 0.04, markerscale = 1,prop=fontProp)
	pdfFile = filename+".max_sign_dist.pdf"
	pngFile = filename+".max_sign_dist.png"
	print "saving to figure:",pdfFile 
	plt.savefig(pdfFile, format = "pdf")
	print "saving to figure:",pngFile 
	plt.savefig(pngFile, format = "png")
	



def plot_rank_dist_hist(rank_statistics,snp_pvals,latent_pvals,filename,hist_color="#3399FF",fig_label='',
		show_y_label=True,show_x_label=True,show_legend=True):
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt

#	plt.figure(figsize=(3.5,4.2))
#	plt.axes([0.17,0.28,0.79,0.69])


	group_names =["d$=$0","0$<$d$\leq$5", "5$<$d$\leq$50", "50$<$d", "Other chr."]
	if latent_pvals:
		sub_group_names = ["Neither significant","One significant","Both significant"]
	else:
		sub_group_names = ["Caus. not significant","Caus. significant"]
	group_counts = [[0]*len(sub_group_names) for i in range(len(group_names))]
#	group_counts = [0,0,0,0,0]
	d_bin_thresholds = [0,1,5001,50001]
	valid_count = 0
	for i, smd in enumerate(rank_statistics):
		if smd!= -2:
			valid_count += 1
			if latent_pvals:
				sign_group = int(math.log(snp_pvals[i])<=-6.6354837468149119)+int(math.log(latent_pvals[i])<=-6.6354837468149119)
			else:
				sign_group = int(math.log(snp_pvals[i])<=-6.6354837468149119)
			if smd <0:
				group_counts[-1][sign_group]+=1
			elif smd>d_bin_thresholds[-1]:
				group_counts[-2][sign_group]+=1
			else:
				for i in range(len(d_bin_thresholds)-1):
					dthres1 = d_bin_thresholds[i]
					dthres2 = d_bin_thresholds[i+1]
					if dthres1-0.1<smd<dthres2:
						group_counts[i][sign_group]+=1
		
	group_counts = [[gc2/float(valid_count) for gc2 in gc1] for gc1 in group_counts]

	group_colors = ["#9999FF","#CCFF00","#33CC00"]#,"#3399FF","#FF9900"]
		
	
	print "Group counts:"
	for g_i, gn in enumerate(group_names):
		for sg_i, sgn in enumerate(sub_group_names):
			print gn+", "+sgn+": "+str(group_counts[g_i][sg_i])
			
	
	plt.figure(figsize=(3.5,4.2))
	plt.axes([0.17,0.28,0.79,0.69])
	#Plotting the histogram!! (manually)
	for i, gc in enumerate(group_counts):
		for j in range(len(sub_group_names)):
			if i==len(group_counts)-1:			
				plt.bar(i,sum(gc[j:]),color=group_colors[j],label=sub_group_names[j])
			else:
				plt.bar(i,sum(gc[j:]),color=group_colors[j])


	minVal = 0
	maxVal = len(group_counts)-0.2
	x_range = maxVal - minVal
	y_max = max(group_counts)	 
	print group_counts,y_max
	plt.axis([minVal-0.035*x_range,maxVal+0.035*x_range,-0.035*1.18,1.035*1.18])
	
	x_ticks_pos = []
	for i in range(len(group_names)):
		x_ticks_pos.append(i+0.4)
	plt.xticks(x_ticks_pos,group_names,fontsize="medium",rotation=45,horizontalalignment='right')
	if show_y_label:
		plt.ylabel("Frequency",fontsize='large')
	plt.text(0.2,1.06,fig_label,fontsize=15)
	if show_x_label:
		plt.xlabel("Distance to nearest causative (kb)",fontsize='medium')
	if show_legend:
		fontProp = matplotlib.font_manager.FontProperties(size=11)
		if latent_pvals:
			plt.legend(loc = 9, numpoints = 1, handlelen = 0.04, markerscale = 1,prop=fontProp)
		else:
			plt.legend(loc = 1, numpoints = 1, handlelen = 0.04, markerscale = 1,prop=fontProp)
	pdfFile = filename+".ranks_dist.pdf"
	pngFile = filename+".ranks_dist.png"
	print "saving to figure:",pdfFile 
	plt.savefig(pdfFile, format = "pdf")
	print "saving to figure:",pngFile 
	plt.savefig(pngFile, format = "png")
	



def plot_rank_hist(r_list,loaded_phen_mafs,max_pval_list,filename,fig_label='',
		show_y_label=True,show_x_label=True,show_legend=True,rank_origins=None):	
	"""
	Plot the histograms and other figures..
	"""
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	
	
	#Split data up by MAFs
	filtered_r_list = []
	filtered_origins = []
	maf_thresholds = [0.1,0.2,0.3,0.4,0.5]
	maf_r_list = [[] for i in range(len(maf_thresholds))]
	sign_count = 0
	for i, maf in enumerate(loaded_phen_mafs):
		maf_i = 0
		try:
			if math.log10(max_pval_list[i])<=-6.6354837468149119: #Filter out non sign. ones.
				maf_r_list[maf_i].append(r_list[i])			
				while maf_i <4 and maf > maf_thresholds[maf_i]:
					maf_i += 1
					maf_r_list[maf_i].append(r_list[i])			
				sign_count += 1
				filtered_r_list.append(r_list[i])
				if rank_origins:
					filtered_origins.append(rank_origins[i])
		except Exception, err_str:
			print "maf too big?:", maf, maf_i,", error str:", err_str
	print "Significant fraction was:",sign_count/float(len(loaded_phen_mafs))
	

	
	plt.figure(figsize=(3.5,4.2))
	plt.axes([0.17,0.28,0.79,0.69])
	#plt.axes([0.14,0.13,0.81,0.83])
	maf_thresholds = [0.1,0.2,0.3,0.4,0.5]

	#Binning (manually)
	r_bin_str =["r$=$1","1$<$r$\leq$10", "10$<$r$\leq$10$^{2}$", "10$^{2}<$r$\leq$10$^{3}$", "10$^{3}<$r"]
	r_bin_thresholds = [1,2,11,101,1001]
	#maf_colors = ["blue","green","red","cyan","purple"]
	maf_colors = ["#9999FF","#CCFF00","#33CC00","#3399FF","#FF9900"]
	origin_colors = ["#3399FF","#FF9900"]
	origin_names = ["First causative","Second causative"]
	max_bin_counts = []
	for i in range(len(r_bin_thresholds)-1):
		if rank_origins:
			c_count = 0
			l_count = 0
			for j, r in enumerate(filtered_r_list):
				
				if r_bin_thresholds[i]<=r<r_bin_thresholds[i+1]:
					c_count += 1
					if filtered_origins[j]==2:
						l_count += 1
			plt.bar(i,c_count/float(sign_count),color=origin_colors[0])
			plt.bar(i,l_count/float(sign_count),color=origin_colors[1])
			max_bin_counts.append(c_count/float(sign_count))
			max_bin_counts.append(l_count/float(sign_count))
		else:
			for maf_i in range(len(maf_thresholds)):
				count = 0
				for r in maf_r_list[maf_i]:
					if r_bin_thresholds[i]<=r<r_bin_thresholds[i+1]:
						count += 1
				plt.bar(i,count/float(sign_count),color=maf_colors[maf_i])
				if maf_i ==0:
					max_bin_counts.append(count/float(sign_count))
		
	if rank_origins:
		c_count = 0
		l_count = 0
		for j, r in enumerate(filtered_r_list):
			
			if r_bin_thresholds[-1]<r:
				c_count += 1
				if filtered_origins[j]==2:
					l_count += 1
		plt.bar(4,c_count/float(sign_count),color=origin_colors[0],label=origin_names[0])
		plt.bar(4,l_count/float(sign_count),color=origin_colors[1],label=origin_names[1])
		max_bin_counts.append(c_count/float(sign_count))
		max_bin_counts.append(l_count/float(sign_count))
	else:
		maf_thresholds.append(0)
		for maf_i in range(len(maf_thresholds)-1):
			count = 0
			for r in maf_r_list[maf_i]:
				if r_bin_thresholds[-1]<r: 
					count += 1
			plt.bar(len(r_bin_thresholds)-1,count/float(sign_count),color=maf_colors[maf_i],label=str(maf_thresholds[maf_i-1])+"$<$MAF$\leq$"+str(maf_thresholds[maf_i]))
			if maf_i ==0:
				max_bin_counts.append(count/float(sign_count))

	minVal = 0
	maxVal = len(r_bin_thresholds)-0.2
	x_range = maxVal - minVal
	y_max = max(max_bin_counts)	 
	print max_bin_counts,y_max
	#plt.axis([minVal-0.035*x_range,maxVal+0.035*x_range,-0.035*y_max,1.2*y_max])
	plt.axis([minVal-0.035*x_range,maxVal+0.035*x_range,-0.035*1.18,1.035*1.18])
	
	x_ticks_pos = []
	for i in range(len(r_bin_thresholds)):
		x_ticks_pos.append(i+0.4)
	plt.xticks(x_ticks_pos,r_bin_str,fontsize="medium",rotation=45,horizontalalignment='right')
	if show_y_label:
		plt.ylabel("Frequency",fontsize='large')
	plt.text(0.2,1.06,fig_label,fontsize=15)
	if show_x_label:
		if rank_origins:
			plt.xlabel("P-value rank",fontsize='medium')
		else:
			plt.xlabel("P-value rank of the causative SNP",fontsize='medium')
	if show_legend:
		fontProp = matplotlib.font_manager.FontProperties(size=11)
		if rank_origins:
			plt.legend(loc = 1, numpoints = 1, handlelen = 0.04, markerscale = 1,prop=fontProp)
		else:
			plt.legend(loc = 9, numpoints = 1, handlelen = 0.04, markerscale = 1,prop=fontProp)
	pdfFile = filename+".ranks.pdf"
	pngFile = filename+".ranks.png"
	print "saving to figure:",pdfFile 
	plt.savefig(pdfFile, format = "pdf")
	print "saving to figure:",pngFile 
	plt.savefig(pngFile, format = "png")
	plt.clf()



def plot_both_ranks_hist(first_ranks,second_ranks,first_mafs,second_mafs,max_pval_list,filename,fig_label='',
		show_y_label=True,show_x_label=True,show_legend=True,rank_origins=None):	
	"""
	Plot the histograms and other figures..
	"""
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	
	
	#Split data up by MAFs
	new_second_ranks = []
	new_second_mafs = []
	maf_thresholds = [0.1,0.2,0.3,0.4,0.5]
	first_maf_r_list = [[] for i in range(len(maf_thresholds))]
	sign_count = 0
	for i, maf in enumerate(first_mafs):
		maf_i = 0
		try:
			if math.log10(max_pval_list[i])<=-6.6354837468149119: #Filter out non sign. ones.
				first_maf_r_list[maf_i].append(first_ranks[i])			
				while maf_i <4 and maf > maf_thresholds[maf_i]:
					maf_i += 1
					first_maf_r_list[maf_i].append(first_ranks[i])			
				sign_count += 1
				new_second_ranks.append(second_ranks[i])
				new_second_mafs.append(second_mafs[i])
		except Exception, err_str:
			print "maf too big?:", maf, maf_i,", error str:", err_str
	print "Significant fraction was:",sign_count/float(len(first_ranks))
	second_mafs = new_second_mafs 
	second_ranks = new_second_ranks 
	
	second_maf_r_list = [[] for i in range(len(maf_thresholds))]
	for i, maf in enumerate(second_mafs):
		maf_i = 0
		try:
			second_maf_r_list[maf_i].append(second_ranks[i])			
			while maf_i <4 and maf > maf_thresholds[maf_i]:
				maf_i += 1
				second_maf_r_list[maf_i].append(second_ranks[i])			
		except Exception, err_str:
			print "maf too big?:", maf, maf_i,", error str:", err_str

	
	plt.figure(figsize=(3.5,4.2))
	plt.axes([0.17,0.28,0.79,0.69])
	#plt.axes([0.14,0.13,0.81,0.83])
	maf_thresholds = [0.1,0.2,0.3,0.4,0.5]

	#Binning (manually)
	r_bin_str =["r$=$1","1$<$r$\leq$10", "10$<$r$\leq$10$^{2}$", "10$^{2}<$r$\leq$10$^{3}$", "10$^{3}<$r"]
	r_bin_thresholds = [1,2,11,101,1001]
	maf_colors = ["#9999FF","#CCFF00","#33CC00","#3399FF","#FF9900"]
	max_bin_counts = []
	for i in range(len(r_bin_thresholds)-1):
		#for maf_i in range(len(maf_thresholds)):
		maf_i=0
		count = 0
		for r1 in first_maf_r_list[maf_i]:
			if r_bin_thresholds[i]<=r1<r_bin_thresholds[i+1]:
				count += 1
		plt.bar(i,count/float(sign_count),color=maf_colors[maf_i+3],width=0.4)
		max_bin_counts.append(count/float(sign_count))
		count = 0
		for r2 in second_maf_r_list[maf_i]:
			if r_bin_thresholds[i]<=r2<r_bin_thresholds[i+1]:
				count += 1
		plt.bar(i+0.4,count/float(sign_count),color=maf_colors[maf_i+4],width=0.4)
		max_bin_counts.append(count/float(sign_count))
		
	#maf_thresholds.append(0)
	#for maf_i in range(len(maf_thresholds)-1):
	maf_i=0
	count = 0
	for r in first_maf_r_list[maf_i]:
		if r_bin_thresholds[-1]<r: 
			count += 1
	plt.bar(len(r_bin_thresholds)-1,count/float(sign_count),color=maf_colors[maf_i+3],width=0.4,label="Top ranked causative")
	max_bin_counts.append(count/float(sign_count))
	count = 0
	for r in second_maf_r_list[maf_i]:
		if r_bin_thresholds[-1]<r: 
			count += 1
	#plt.bar(len(r_bin_thresholds)-1+0.4,count/float(sign_count),color=maf_colors[maf_i],width=0.4,label="$"+str(maf_thresholds[maf_i-1])+"< $MAF$\leq"+str(maf_thresholds[maf_i])+"$")
	plt.bar(len(r_bin_thresholds)-1+0.4,count/float(sign_count),color=maf_colors[maf_i+4],width=0.4,label="Other causative")
	max_bin_counts.append(count/float(sign_count))

	minVal = 0
	maxVal = len(r_bin_thresholds)-0.2
	x_range = maxVal - minVal
	y_max = max(max_bin_counts)	 
	print max_bin_counts,y_max
	#plt.axis([minVal-0.035*x_range,maxVal+0.035*x_range,-0.035*y_max,1.2*y_max])
	plt.axis([minVal-0.035*x_range,maxVal+0.035*x_range,-0.035*1.18,1.035*1.18])
	
	x_ticks_pos = []
	for i in range(len(r_bin_thresholds)):
		x_ticks_pos.append(i+0.4)
	plt.xticks(x_ticks_pos,r_bin_str,fontsize="medium",rotation=45,horizontalalignment='right')
	if show_y_label:
		plt.ylabel("Frequency",fontsize='large')
	plt.text(0.2,1.06,fig_label,fontsize=15)
	if show_x_label:
		if rank_origins:
			plt.xlabel("P-value rank",fontsize='medium')
		else:
			plt.xlabel("P-value rank of the causative SNP",fontsize='medium')
	if show_legend:
		fontProp = matplotlib.font_manager.FontProperties(size=11)
		plt.legend(loc = 1, numpoints = 1, handlelen = 0.04, markerscale = 1,prop=fontProp)
	pdfFile = filename+".ranks.pdf"
	pngFile = filename+".ranks.png"
	print "saving to figure:",pdfFile 
	plt.savefig(pdfFile, format = "pdf")
	print "saving to figure:",pngFile 
	plt.savefig(pngFile, format = "png")
	plt.clf()




def plot_dist_hist(d_list,loaded_phen_mafs,max_pval_list,filename,fig_label='',
		show_y_label=True,show_x_label=True,show_legend=True,dist_origins=None):	
	"""
	Plot the dist histogram.
	"""
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	
	filtered_d_list = []
	filtered_origins = []

	#Split data up by MAFs
	maf_thresholds = [0.1,0.2,0.3,0.4,0.5]
	maf_d_list = [[] for i in range(len(maf_thresholds))]
	sign_count = 0
	for i, maf in enumerate(loaded_phen_mafs):
		maf_i = 0
		try:
			if math.log10(max_pval_list[i])<=-6.6354837468149119:
				maf_d_list[maf_i].append(d_list[i])
				while maf_i <4 and maf > maf_thresholds[maf_i]:
					maf_i += 1
					maf_d_list[maf_i].append(d_list[i])
				sign_count += 1
				filtered_d_list.append(d_list[i])
				if dist_origins:
					filtered_origins.append(dist_origins[i])
		except Exception, err_str:
			print "maf too big?:", maf, maf_i,", error str:", err_str
	print "Significant fraction was:",sign_count/float(len(loaded_phen_mafs))

	
	plt.figure(figsize=(3.5,4.2))
	#plt.axes([0.17,0.26,0.81,0.73])
	plt.axes([0.17,0.28,0.79,0.69])

	#Binning (manually)
	#d_bin_str =["$ d=0 $","$ 0<d\le 2 \textrm{kb} $", "$ 2 \textrm{kb}<d \le 20\textrm{kb}$", "$ 20 \textrm{kb}<d \le 100\textrm{kb} $", "$ 100 \textrm{kb}<d $"]
	d_bin_str =["d$=$0","0$<$d$\leq$5", "5$<$d$\leq$50", "50$<$d", "Other chr."]
	d_bin_thresholds = [0,1,5001,50001]
	origin_colors = ["#3399FF","#FF9900"]
	origin_names = ["First causative","Second causative"]
	maf_colors = ["#9999FF","#CCFF00","#33CC00","#3399FF","#FF9900"]
	max_bin_counts = []
	for i in range(len(d_bin_thresholds)-1):
		if dist_origins:
			c_count = 0
			l_count = 0
			for j, d in enumerate(filtered_d_list):
				if d_bin_thresholds[i]<=d<d_bin_thresholds[i+1]:
					c_count += 1
					if filtered_origins[j]==2:
						l_count += 1
			plt.bar(i,c_count/float(len(d_list)),color=origin_colors[0])
			plt.bar(i,l_count/float(len(d_list)),color=origin_colors[1])
			max_bin_counts.append(c_count/float(len(d_list)))
			#max_bin_counts.append(l_count/float(len(d_list)))
		else:
			for maf_i in range(len(maf_thresholds)):
				count = 0
				for d in maf_d_list[maf_i]:
					if d_bin_thresholds[i]<=d<d_bin_thresholds[i+1]:
						count += 1
				plt.bar(i,count/float(len(d_list)),color=maf_colors[maf_i])
				if maf_i ==0:
					max_bin_counts.append(count/float(len(d_list)))
		
	#Drawing the bin for >50kb
	if dist_origins:
		c_count = 0
		l_count = 0
		for j, d in enumerate(filtered_d_list):
			if d_bin_thresholds[-1]<=d:
				c_count += 1
				if filtered_origins[j]==2:
					l_count += 1
		plt.bar(3,c_count/float(len(d_list)),color=origin_colors[0])
		plt.bar(3,l_count/float(len(d_list)),color=origin_colors[1])
		max_bin_counts.append(c_count/float(len(d_list)))
		#max_bin_counts.append(l_count/float(len(d_list)))
	else:
		for maf_i in range(len(maf_thresholds)):  
			count = 0
			for d in maf_d_list[maf_i]:
				if d_bin_thresholds[-1]<=d:
					count += 1
			plt.bar(3,count/float(len(d_list)),color=maf_colors[maf_i])
			if maf_i ==0:
				max_bin_counts.append(count/float(len(d_list)))

	if dist_origins:
		c_count = 0
		l_count = 0
		for j, d in enumerate(filtered_d_list):
			if d<0:
				c_count += 1
				if filtered_origins[j]==2:
					l_count += 1
		plt.bar(4,c_count/float(len(d_list)),color=origin_colors[0],label=origin_names[0])
		plt.bar(4,l_count/float(len(d_list)),color=origin_colors[1],label=origin_names[1])
		max_bin_counts.append(c_count/float(len(d_list)))
		#max_bin_counts.append(l_count/float(len(d_list)))
	else:
		maf_thresholds.append(0)
		for maf_i in range(len(maf_thresholds)-1):
			count = 0
			for d in maf_d_list[maf_i]:
				if 0>d: #on a different chr.
					count += 1
			plt.bar(len(d_bin_thresholds),count/float(len(d_list)),color=maf_colors[maf_i],label=str(maf_thresholds[maf_i-1])+"$<$MAF$\leq$"+str(maf_thresholds[maf_i]))
			if maf_i ==0:
				max_bin_counts.append(count/float(len(d_list)))

	minVal = 0
	maxVal = len(d_bin_thresholds)+1-0.2
	x_range = maxVal - minVal
	y_max = max(max_bin_counts)	 
	print max_bin_counts,y_max
	#plt.axis([minVal-0.035*x_range,maxVal+0.035*x_range,-0.035*y_max,1.2*y_max])
	plt.axis([minVal-0.035*x_range,maxVal+0.035*x_range,-0.035*1.18,1.035*1.18])
	
	x_ticks_pos = []
	for i in range(len(d_bin_thresholds)+1):
		x_ticks_pos.append(i+0.4)
	plt.xticks(x_ticks_pos,d_bin_str,fontsize="medium",rotation=45,horizontalalignment='right')
	if show_y_label:
		plt.ylabel("Frequency",fontsize='large')
	plt.text(0.2,1.06,fig_label,fontsize=15)
	if show_x_label:
		plt.xlabel("Distance to smallest p-value (kb)",fontsize='medium')
	fontProp = matplotlib.font_manager.FontProperties(size=11)
	if show_legend:
		if dist_origins:
			plt.legend(loc = 1, numpoints = 1, handlelen = 0.04, markerscale = 1,prop=fontProp)
		else:
			plt.legend(loc = 9, numpoints = 1, handlelen = 0.04, markerscale = 1,prop=fontProp)			
	pdfFile = filename+".dist.pdf"
	pngFile = filename+".dist.png"
	print "saving to figure:",pdfFile 
	plt.savefig(pdfFile, format = "pdf")
	print "saving to figure:",pngFile 
	plt.savefig(pngFile, format = "png")
	plt.clf()	
	



def plot_rank_3D_hist(r_list1,r_list2,min_pvals,filename,
		show_y_label=True,show_x_label=True):	
	"""
	Plot 3D histogram..
	"""
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	
	#Filter non-significant results out.
	sign_count = 0
	filtered_r_list1 = []
	filtered_r_list2 = []
	for r1,r2,mp in zip(r_list1,r_list2,min_pvals):
		if math.log10(mp)<=-6.6354837468149119: 
			filtered_r_list1.append(r1)
			filtered_r_list2.append(r2)
	
	print "Significant fractions were:",len(filtered_r_list1)/float(len(r_list1)),len(filtered_r_list2)/float(len(r_list2))

	
	plt.figure(figsize=(6,5))
	#plt.subplot(2,1,1)
	plt.axes([0.22,0.2,0.72,0.78])
	#plt.axes([0.14,0.13,0.81,0.83])
	maf_thresholds = [0.1,0.2,0.3,0.4,0.5]

	#Binning (manually)
	import numpy as np
	r_bin_str =["r$=$1","1$<$r$\leq$10", "10$<$r$\leq$10$^{2}$", "10$^{2}<$r$\leq$10$^{3}$", "10$^{3}<$r"]
	counts = np.zeros((len(r_bin_str),len(r_bin_str)))
	r_bin_thresholds = [1,2,11,101,1001]
	for r1,r2 in zip(r_list1,r_list2):
		ci1 = 4
		cur_thres = r_bin_thresholds[ci1]
		while ci1>0 and r1<cur_thres:
			ci1-=1
			cur_thres = r_bin_thresholds[ci1]
		ci2 = 4
		cur_thres = r_bin_thresholds[ci2]
		while ci2>0 and r2<cur_thres:
			ci2-=1
			cur_thres = r_bin_thresholds[ci2]
		counts[ci1,ci2]+=1
	tot_count = sum(sum(counts))
	for i in range(len(counts)):
		for j in range(len(counts)):
			counts[i,j]=counts[i,j]/tot_count
		
	plt.imshow(counts,interpolation='nearest')
	plt.colorbar()

	plt.axis([-0.5-len(r_bin_thresholds)*0.035,-0.5+len(r_bin_thresholds)*1.035,-0.5-len(r_bin_thresholds)*0.035,-0.5+len(r_bin_thresholds)*1.035])
	
	x_ticks_pos = []
	for i in range(len(r_bin_thresholds)):
		x_ticks_pos.append(i)
	plt.xticks(x_ticks_pos,r_bin_str,fontsize="medium",rotation=45)
	plt.yticks(x_ticks_pos,r_bin_str,fontsize="medium",rotation=45)
	if show_y_label:
		plt.ylabel("P-value rank of the first causative SNP",fontsize='medium')
	if show_x_label:
		plt.xlabel("P-value rank of the second causative SNP",fontsize='medium')
	pdfFile = filename+".3D_ranks.pdf"
	pngFile = filename+".3D_ranks.png"
	print "saving to figure:",pdfFile 
	plt.savefig(pdfFile, format = "pdf")
	print "saving to figure:",pngFile 
	plt.savefig(pngFile, format = "png")
	plt.clf()


def plot_dist_3D_hist(r_list1,r_list2,min_pvals,filename,
		show_y_label=True,show_x_label=True):	
	"""
	Plot 3D histogram..
	"""
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	
	#Filter non-significant results out.
	sign_count = 0
	filtered_r_list1 = []
	filtered_r_list2 = []
	for r1,r2,mp in zip(r_list1,r_list2,min_pvals):
		if math.log10(mp)<=-6.6354837468149119: 
			filtered_r_list1.append(r1)
			filtered_r_list2.append(r2)
	
	print "Significant fractions were:",len(filtered_r_list1)/float(len(r_list1)),len(filtered_r_list2)/float(len(r_list2))

	
	plt.figure(figsize=(6,5))
	#plt.subplot(2,1,1)
	plt.axes([0.22,0.2,0.72,0.78])
	#plt.axes([0.14,0.13,0.81,0.83])
	maf_thresholds = [0.1,0.2,0.3,0.4,0.5]

	#Binning (manually)
	import numpy as np
	r_bin_str =["d$=$0","0$<$d$\leq$5", "5$<$d$\leq$50", "50$<$d", "Other chr."]
	counts = np.zeros((len(r_bin_str),len(r_bin_str)))
	r_bin_thresholds = [0,1,5001,50001]
	for r1,r2 in zip(r_list1,r_list2):
		if r1<0:
			ci1=4
		else:
			ci1 = 3
			cur_thres = r_bin_thresholds[ci1]
			while ci1>0 and r1<cur_thres:
				ci1-=1
				cur_thres = r_bin_thresholds[ci1]
		if r2<0:
			ci2=4
		else:
			ci2 = 3
			cur_thres = r_bin_thresholds[ci2]
			while ci2>0 and r2<cur_thres:
				ci2-=1
				cur_thres = r_bin_thresholds[ci2]
		counts[ci1,ci2]+=1
	tot_count = sum(sum(counts))
	for i in range(len(counts)):
		for j in range(len(counts)):
			counts[i,j]=counts[i,j]/tot_count
		
	plt.imshow(counts,interpolation='nearest')
	plt.colorbar()

	plt.axis([-0.5-len(r_bin_thresholds)*0.035,-0.5+len(r_bin_thresholds)*1.035,-0.5-len(r_bin_thresholds)*0.035,-0.5+len(r_bin_thresholds)*1.035])
	
	x_ticks_pos = []
	for i in range(len(r_bin_thresholds)):
		x_ticks_pos.append(i)
	plt.xticks(x_ticks_pos,r_bin_str,fontsize="medium",rotation=45)
	plt.yticks(x_ticks_pos,r_bin_str,fontsize="medium",rotation=45)
	if show_y_label:
		plt.xlabel("Distance to smallest p-value from second causative SNP (kb)",fontsize='medium')
	if show_x_label:
		plt.xlabel("Distance to smallest p-value from first causative SNP (kb)",fontsize='medium')
	pdfFile = filename+".3D_ranks.pdf"
	pngFile = filename+".3D_ranks.png"
	print "saving to figure:",pdfFile 
	plt.savefig(pdfFile, format = "pdf")
	print "saving to figure:",pngFile 
	plt.savefig(pngFile, format = "png")
	plt.clf()




		
if __name__=='__main__':
	_run_()
	print "Done!"

