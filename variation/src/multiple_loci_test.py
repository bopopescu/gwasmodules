"""
Tests involving stepwise regressions, and model selection.

Option:

	-i ...			The run phenotype index/indices.  
	-o ...			The run_id, used as unique identifier for this run.  
	-s			Collect results and plot things (Does not generate pvalue files...) 

	-t ...			What data set is used. (75 is default)

	-n ...			Number of SNPs (phenotypes) per node, default is 1
	-d ...			Debug filter, random fraction of phenotypes or snps will be used.

	-l ...			Type of latent variable: random_snps (default), random, pc_split, etc..
	-h ...			Heritability in percentages, possible values are 1,10,25,50,75,90,99
	
	-m ...			How to generate the phenotypes: plus, xor, or

	--save_plots           	Plot Manhattan plots
	--phen_file=...		File where the phenotypes will be saved, and loaded from.
	--sim_phen		Simulate phenotype, write to phenotype file.
	--parallel		Run parallel on the cluster
	--num_steps=...		Number of steps in the regression. (Default is 10)
	--herit_plots=...	Plot heritabilities
	--var_plots		Plot variation analysis


Examples:
python multiple_loci_test.py run_index test -i 1,5 -a kw,emmax -b most_normal -r ~/Projects/Data/phenotypes/phen_raw_092910.tsv 

"""
import matplotlib
matplotlib.use('Agg')

import cPickle
import scipy as sp
import linear_models as lm
import scipy.linalg as linalg
import phenotypeData
import snpsdata
import sys
import os
import env
import random
import dataParsers as dp
import util
import gwaResults as gr
import analyze_gwas_results as agr
import traceback
import getopt
import time
import pdb


mapping_methods = ['LM', 'KW', 'EX', 'Stepw_LM', 'Stepw_EX'] #5 in total

def parse_parameters():
	'Parse the parameters into a dict, etc.'
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)

	long_options_list = ['save_plots', 'phen_file=', 'sim_phen', 'num_steps=', 'parallel', 'herit_plots=', 'var_plots']
	try:
		opts, args = getopt.getopt(sys.argv[1:], "i:o:t:k:n:d:l:m:h:s", long_options_list)

	except:
		traceback.print_exc()
		print __doc__
		sys.exit(2)

	p_dict = {'number_per_run':20, 'debug_filter':1.0, 'summarize':False,
		'latent_variable':'random_snp', 'phenotype_model':'plus', 'run_id':'mlt',
		'mapping_method':'emmax', 'heritability':50, 'save_plots':False, 'call_method_id':75,
		'phen_file':env.env['data_dir'] + 'multi_locus_phen.pickled', 'num_steps':10,
		'phen_index':None, 'sim_phen':False, 'parallel':False, 'herit_plots':None,
		'var_plots':False}


	for opt, arg in opts:
		if opt in ('-i'): p_dict['phen_index'] = util.parse_ids(arg)
		elif opt in ('-o'): p_dict['run_id'] = arg
		elif opt in ('-t'): p_dict['call_method_id'] = int(arg)
		elif opt in ('-n'): p_dict['number_per_run'] = int(arg)
		elif opt in ('-m'): p_dict['phenotype_model'] = arg
		elif opt in ('-d'): p_dict['debug_filter'] = float(arg)
		elif opt in ('-l'): p_dict['latent_variable'] = arg
		elif opt in ("-s"): p_dict['summarize'] = True
		elif opt in ('-h'): p_dict['heritability'] = int(arg)
		elif opt in ("--phen_file"): p_dict['phen_file'] = arg
		elif opt in ("--save_plots"): p_dict['save_plots'] = True
		elif opt in ("--sim_phen"): p_dict['sim_phen'] = True
		elif opt in ("--num_steps"): p_dict['num_steps'] = int(arg)
		elif opt in ("--parallel"): p_dict['parallel'] = True
		elif opt in ("--herit_plots"): p_dict['herit_plots'] = util.parse_ids(arg)
		elif opt in ("--var_plots"): p_dict['var_plots'] = True
		else:
			print "Unkown option!!\n"
			print __doc__
			sys.exit(2)

	print p_dict, args
	return p_dict, args




def run_parallel(run_id, start_i, stop_i, latent_var, heritability, phen_model,
		phen_file, summary_run, call_method_id, num_steps, cluster='gmi'):
	"""
	If no mapping_method, then analysis run is set up.
	"""
	phen_ids_str = '%d-%d' % (start_i, stop_i)
	job_id = '%s_%s_lv%s_h%d_m%s_ns%d_t%d' % (run_id, phen_ids_str , latent_var,
						heritability, phen_model, num_steps, call_method_id)
	file_prefix = env.env['results_dir'] + job_id

	#Cluster specific parameters	
	if cluster == 'gmi': #GMI cluster.  #Sumit to the memory node
		shstr = '#!/bin/bash\n'
		shstr += '#$ -S /bin/bash\n'
		shstr += '#$ -N %s\n' % job_id
		#shstr += "#$ -q q.norm@blade*\n"
		shstr += '#$ -o %s_job_$JOB_ID.out\n' % file_prefix
		shstr += '#$ -e %s_job_$JOB_ID.err\n' % file_prefix
		shstr += 'source /etc/modules-env.sh\n'
		shstr += 'module load scipy/GotoBLAS2/0.9.0\n'
		shstr += 'module load matplotlib/1.0.0\n'
		shstr += 'module load mysqldb/1.2.3\n'
		shstr += 'export GOTO_NUM_THREADS=1\n'
		#shstr += '#$ -cwd /home/GMI/$HOME\n'
		#shstr += '#$ -M bjarni.vilhjalmsson@gmi.oeaw.ac.at\n\n'

	elif cluster == 'usc':  #USC cluster.
		shstr = "#!/bin/csh\n"
		shstr += "#PBS -l walltime=%s \n" % '72:00:00'
		shstr += "#PBS -l mem=%s \n" % '2950mb'
		shstr += "#PBS -q cmb\n"
		shstr += "#PBS -N p%s \n" % job_id

	shstr += "python %smultiple_loci_test.py -i %s -l %s -h %d -m %s --phen_file=%s -t %d --num_steps=%d " % \
			(env.env['script_dir'], phen_ids_str, latent_var, heritability, phen_model, phen_file, call_method_id, num_steps)
	if summary_run:
		shstr += '-s '

	#shstr += "> " + file_prefix + "_job.out) >& " + file_prefix + "_job.err\n"
	print '\n', shstr, '\n'
	script_file_name = run_id + ".sh"
	f = open(script_file_name, 'w')
	f.write(shstr)
	f.close()

	#Execute qsub script
	if cluster == 'gmi':
		os.system("qsub " + script_file_name)
	elif cluster == 'usc':
		os.system("qsub " + script_file_name)




def get_snps_heritabilities(snps, phenotype):
	Y = sp.mat(phenotype).T
	rss_0 = sp.var(Y) * len(phenotype)
	X_0 = sp.mat(sp.ones((len(phenotype), 1)))

	h_expl = []
	for snp in snps:
		rss = linalg.lstsq(sp.hstack([X_0, sp.mat(snp).T]), Y)[1]
		h_expl.append(1 - (rss / rss_0))
	return h_expl


def __get_latent_snps__(ets):
	ecotype_info_dict = phenotypeData.get_ecotype_id_info_dict()
	north_south_split_snp = []
	lats = [ecotype_info_dict[int(et)][2] for et in ets]
	m = sp.median(lats)
	for et in ets:
		latitude = ecotype_info_dict[int(et)][2]
		north_south_split_snp.append(1) if latitude > m else north_south_split_snp.append(0)
	pc_snp = []
	K = dp.load_kinship() #All ecotypes
	(evals, evecs) = linalg.eigh(K)
	pc = (sp.mat(evecs).T[-1]).tolist()[0]
	m = sp.median(pc)
	for v in pc:
		pc_snp.append(1) if v > m else pc_snp.append(0)
	return sp.array(north_south_split_snp, dtype='int8'), sp.array(pc_snp, dtype='int8')


def simulate_phenotypes(phen_file, sd, mac_threshold=0, debug_filter=1.0, num_phens=1000):
	"""
	Simulate the phenotypes
	"""
	print 'Generating the phenotypes'
	latent_var_keys = ['random_snp', 'random', 'north_south_split', 'pc_split']
	phenotype_models = ['xor', 'or', 'plus', 'xor2']
	heritabilities = [1, 2, 5, 10, 15, 20, 25, 50] #in percentages

	if mac_threshold > 0:
		sd.filter_mac_snps(mac_threshold)
	num_lines = len(sd.accessions)  #Number of lines
	mafs = sd.get_mafs()["marfs"]
	if debug_filter:
		sd.sample_snps(debug_filter)
	snp_chr_pos_maf_list = sd.get_all_snp_w_info()
	all_indices = range(len(snp_chr_pos_maf_list))
	snp_indices = random.sample(all_indices, num_phens)
	map(all_indices.remove, snp_indices)

	#The first locus..
	snp_chr_pos_maf_list = [snp_chr_pos_maf_list[i] for i in snp_indices]

	#Invert every other SNP (randomize the SNP decoding)
	all_indices = range(len(snp_chr_pos_maf_list))
	invert_indices = random.sample(all_indices, num_phens / 2)
	for i in invert_indices:
		snp, chr, pos, maf = snp_chr_pos_maf_list[i]
		snp_chr_pos_maf_list[i] = (lm.get_anti_snp(snp), chr, pos, maf)

	north_south_split_snp, pc_snp = __get_latent_snps__(sd.accessions)

	phen_dict = {'snp_chr_pos_maf_list': snp_chr_pos_maf_list, 'snp_indices':snp_indices,
			'north_south_split_snp':north_south_split_snp, 'pc_snp':pc_snp}
	for latent_var in latent_var_keys:
		d = {}
		if latent_var == 'random_snp':
			l_snp_indices = random.sample(all_indices, num_phens)
			latent_snps = [snp_chr_pos_maf_list[i][0] for i in l_snp_indices]
			d['latent_chr_pos_maf_list'] = \
				[(snp_chr_pos_maf_list[i][1], snp_chr_pos_maf_list[i][2], \
				snp_chr_pos_maf_list[i][3]) for i in l_snp_indices]
			d['latent_snps'] = latent_snps

		elif latent_var == 'random':
			latent_snps = []
			for i in range(num_phens):
				num_ones = random.randint(1, num_lines - 1)
				l_snp = [0] * num_lines
				one_indices = random.sample(range(num_lines), num_ones)
				for i in one_indices:
					l_snp[i] = 1
				latent_snps.append(sp.array(l_snp, dtype='int8'))
			d['latent_snps'] = latent_snps

		elif latent_var == 'north_south_split':
			latent_snp = north_south_split_snp
			d['latent_snp'] = latent_snp

		elif latent_var == 'pc_snp':
			latent_snp = pc_snp
			d['latent_snp'] = latent_snp

		for h in heritabilities:
			her = h / 100.0
			d2 = {}
			for phen_model in phenotype_models:  #Simulate all phenotype models.
				d3 = {'phenotypes': [], 'h_estimates': [], 'h_loci_est_list': []}
				for i in range(num_phens):
					if latent_var in ['random_snp', 'random']:
						latent_snp = latent_snps[i]
					snp = snp_chr_pos_maf_list[i][0]
					if phen_model == 'xor':
						phenotype = snp ^ latent_snp
					elif phen_model == 'or':
						phenotype = snp | latent_snp
					elif phen_model == 'plus':
						phenotype = snp + latent_snp
					elif phen_model == 'xor2':
						phenotype = (snp ^ latent_snp) + 0.5 * (snp & latent_snp)
					if len(sp.unique(phenotype)) > 1:
						phen_var = sp.var(phenotype, ddof=1)
						error_vector = sp.random.normal(0, 1, size=num_lines)
						error_var = sp.var(error_vector, ddof=1)
						scalar = sp.sqrt((phen_var / error_var) * ((1 - her) / her))
						phenotype = phenotype + error_vector * scalar
						h_est = phen_var / sp.var(phenotype, ddof=1)
						h_est_snp1 = sp.corrcoef(snp, phenotype)[0, 1]
						h_est_snp2 = sp.corrcoef(latent_snp, phenotype)[0, 1]
						#print phen_model, latent_var, her, h_est, h_est_snp1 ** 2, h_est_snp2 ** 2
						d3['h_loci_est_list'].append(h_est)
						d3['h_estimates'].append((h_est_snp1 ** 2, h_est_snp2 ** 2))
					else:
						print 'encountered invalid phenotype for phen_model: %s' % phen_model
						phenotype = None
					d3['phenotypes'].append(phenotype)
				d2[phen_model] = d3
			d[h] = d2
		phen_dict[latent_var] = d


	#phenotype_models for loop ends.
	f = open(phen_file, "wb")
	print "dumping phenotypes to file:", f
	cPickle.dump(phen_dict, f, protocol=2)
	f.close()


def load_phenotypes(phen_file):
	print 'Loading phenotypes and related data'
        f = open(phen_file, "rb")
	phen_dict = cPickle.load(f)
	f.close()
        print 'Loading done..'
	return phen_dict


def summarize_runs(file_prefix, latent_var, heritability, phen_model, phen_d, index_list=None):
	"""
	Summarize runs.. duh
	"""
	pd = phen_d[latent_var][heritability][phen_model]
	if not index_list:
		index_list = range(len(pd['phenotypes']))

	num_pthres = len(pval_thresholds)
	num_winsizes = len(window_sizes)
	summary_dict = {'p_her':[]}
	analysis_methods = ['LM', 'KW', 'EX', 'Stepw_LM_Bonf', 'Stepw_LM_EBIC', 'Stepw_LM_MBIC',
				'Stepw_EX_Bonf', 'Stepw_EX_EBIC', 'Stepw_EX_MBIC']

	for am in analysis_methods:
		if am in ['LM', 'EX', 'KW']:
			d = {'fdrs':sp.zeros((num_pthres, num_winsizes), dtype='double'),
				'tprs':sp.zeros((num_pthres, num_winsizes), dtype='double')}
			d['ks'] = []
			d['medp'] = []
		else:
			if am in ['Stepw_EX_EBIC', 'Stepw_EX_MBIC', 'Stepw_LM_EBIC', 'Stepw_LM_MBIC']:
				d = {'fdrs':sp.zeros((num_winsizes), dtype='double'),
					'tprs':sp.zeros((num_winsizes), dtype='double')}
			else:
				d = {'fdrs':sp.zeros((num_pthres, num_winsizes), dtype='double'),
					'tprs':sp.zeros((num_pthres, num_winsizes), dtype='double')}
			d['perc_var_expl'] = []
			if 'EX' in am:
				d['rem_p_her'] = []
				d['perc_phen_var_expl'] = []
				d['perc_err_var_expl'] = []

		summary_dict[am] = d

	criteria_map = {'mbonf':('Stepw_EX_Bonf', 'Stepw_LM_Bonf'), 'ebics': ('Stepw_EX_EBIC', 'Stepw_LM_EBIC'),
			'mbics': ('Stepw_EX_MBIC', 'Stepw_LM_MBIC')}

	num_files_found = 0
	print '%s_%d_%s_%s' % (file_prefix, heritability, latent_var, phen_model)
	for i in index_list:
		pickled_file = '%s_%d_%s_%s_%dresults.pickled' % (file_prefix, heritability, latent_var, phen_model, i)
		if os.path.isfile(pickled_file):
			num_files_found += 1
			with open(pickled_file) as f:
				r = cPickle.load(f)
			p_her = r['p_her']
			summary_dict['p_her'].append(p_her)
			for am in mapping_methods:
				if am == 'Stepw_LM':
					for criteria in ['mbonf', 'ebics', 'mbics']:
						for rs in ['fdrs', 'tprs']:
							rs_array = sp.array(r[am][criteria][rs])
							rs_array[rs_array == -1.0] = 0
							summary_dict[criteria_map[criteria][1]][rs] += rs_array

						summary_dict[criteria_map[criteria][1]]['perc_var_expl'].append(r[am][criteria]['perc_var_expl'])
				elif am == 'Stepw_EX':
					for criteria in ['mbonf', 'ebics', 'mbics']:
						for rs in ['fdrs', 'tprs']:
							rs_array = sp.array(r[am][criteria][rs])
							rs_array[rs_array == -1.0] = 0
							summary_dict[criteria_map[criteria][0]][rs] += rs_array
						if 'pseudo_heritability' in r[am][criteria]:
							rem_p_her = r[am][criteria]['pseudo_heritability']
						else:
							rem_p_her = r[am][criteria]['remaining_p_her']
						perc_var_expl = r[am][criteria]['perc_var_expl']
						summary_dict[criteria_map[criteria][0]]['rem_p_her'].append(rem_p_her)
						summary_dict[criteria_map[criteria][0]]['perc_var_expl'].append(perc_var_expl)
						rss_0 = r[am]['step_info_list'][0]['rss']
						rss = r[am]['step_info_list'][r[am][criteria]['opt_i']]['rss']

						if p_her > 0.01:
							summary_dict[criteria_map[criteria][0]]['perc_phen_var_expl'].append(1.0 - ((rss * rem_p_her) / (rss_0 * p_her)))
							summary_dict[criteria_map[criteria][0]]['perc_err_var_expl'].append(1.0 - ((rss * (1 - rem_p_her)) / (rss_0 * (1 - p_her))))
				elif am in ['LM', 'KW', 'EX']:
					summary_dict[am]['fdrs'] += sp.array(r[am]['fdrs'])
					summary_dict[am]['tprs'] += sp.array(r[am]['tprs'])
					summary_dict[am]['ks'].append(r[am]['ks_stat']['D'])
					summary_dict[am]['medp'].append(r[am]['med_pval'])
	print 'Found %d results' % num_files_found
	for am in analysis_methods:
		for k in ['fdrs', 'tprs']:
			summary_dict[am][k] = summary_dict[am][k] / float(num_files_found)

	return summary_dict


def plot_tprs_fdrs(file_prefix, summary_dict):
	"""
	Plot various things relating to run summaries
	"""
	import pylab
	import matplotlib.font_manager
	prop = matplotlib.font_manager.FontProperties(size=10)

	#Heritabilities..
	# - histogram of each category
	# - pseudoheritabilities vs. ks and med pval. of KW and LM

	# TPRs vs. FDRs
	am_list = ['LM', 'KW', 'EX', 'Stepw_LM_Bonf', 'Stepw_EX_Bonf']
	am_colors = ['r', 'g', 'b', 'r', 'b']
	am_ls = ['--', '--', '--', '-', '-']
	am_dot_list = ['Stepw_EX_EBIC', 'Stepw_EX_MBIC', 'Stepw_LM_EBIC', 'Stepw_LM_MBIC']
	am_dot_colors = ['#22DD66', '#22DD66', '#DD2266', '#DD2266']
	am_dot_marker = ['s', '^', 's', '^']

	for w_i, ws in enumerate(window_sizes):
		pylab.figure(figsize=(10, 8))
		pylab.axes([0.06, 0.05, 0.9, 0.94])
		for am, amc, amls in zip(am_list, am_colors, am_ls):
			xs = sp.zeros(len(pval_thresholds))
			ys = sp.zeros(len(pval_thresholds))
			for pt_i, pt in enumerate(pval_thresholds):
				ys[pt_i] = summary_dict[am]['tprs'][pt_i][w_i]
				xs[pt_i] = summary_dict[am]['fdrs'][pt_i][w_i]
			pylab.plot(xs, ys, label=am, color=amc, ls=amls, alpha=0.8, marker='.')
		for am, amc, amm in zip(am_dot_list, am_dot_colors, am_dot_marker):
			pylab.plot(summary_dict[am]['fdrs'][w_i], summary_dict[am]['tprs'][w_i], label=am, marker=amm,
				ls='', color=amc)
		png_file = '%s_w%d.png' % (file_prefix, ws)
		pylab.ylabel('Power')
		pylab.xlabel('FDR')
		pylab.legend(loc=4, prop=prop, numpoints=1, scatterpoints=1)
		x_min, x_max = pylab.xlim()
		x_range = x_max - x_min
		y_min, y_max = pylab.ylim()
		y_range = y_max - y_min
		pylab.axis([x_min - 0.025 * x_range, x_max + 0.025 * x_range,
				y_min - 0.025 * y_range, y_max + 0.025 * y_range])

		pylab.savefig(png_file)
		pylab.clf()



def plot_herit_hist(file_prefix, her_dict, latent_var, phen_model):
	import pylab
	import matplotlib.font_manager
	prop = matplotlib.font_manager.FontProperties(size=10)
	file_prefix += '_%s_%s' % (latent_var, phen_model)
	png_file_name = file_prefix + '_h%s_hist.png' % ('_'.join(map(str, her_dict.keys())))
	max_bin_count = 0
	for h in sorted(her_dict):
		bin_counts, bins, patch_list = pylab.hist(her_dict[h]['p_her'], range=(0, 0.8), bins=25, alpha=0.6)
		max_bin_count = max(max_bin_count, max(bin_counts))
		pylab.axvline((h / 100.0), color='k', alpha=0.8, ls='-.')
		pylab.axvline(sp.median(her_dict[h]['p_her']), color='#DD3311', alpha=0.8, ls='-.')
	y_range = max_bin_count - 0
	pylab.axis([-0.8 * 0.025, 0.8 * 1.025, -0.025 * max_bin_count, max_bin_count * 1.025])
	pylab.xlabel('heritability')
	pylab.ylabel('Counts')
	pylab.savefig(png_file_name)

	pylab.figure()
	png_file_name = file_prefix + '_h%s_ks_her_scatter.png' % ('_'.join(map(str, her_dict.keys())))
	for h in sorted(her_dict):
		pylab.plot(her_dict[h]['p_her'], her_dict[h]['LM']['ks'], ls='', marker='.', alpha=0.5, label='herit.=%0.2f' % (h / 100.0))
	pylab.xlabel('pseudo-heritability')
	pylab.ylabel('Kolmogorov-Smirnov statistic')
	pylab.legend(loc=2, prop=prop, numpoints=1, scatterpoints=1)
	pylab.savefig(png_file_name)

	pylab.figure()
	png_file_name = file_prefix + '_h%s_pmed_her_scatter.png' % ('_'.join(map(str, her_dict.keys())))
	for h in sorted(her_dict):
		pylab.plot(her_dict[h]['p_her'], her_dict[h]['LM']['medp'], ls='', marker='.', alpha=0.5, label='herit.=%0.2f' % (h / 100.0))
	pylab.xlabel('pseudo-heritability')
	pylab.ylabel('Median pvalue bias')
	pylab.legend(loc=2, prop=prop, numpoints=1, scatterpoints=1)
	pylab.savefig(png_file_name)


	x = []
	for h in sorted(her_dict):
		pylab.figure()
		png_file_name = file_prefix + '_h%d_hist.png' % (h)
		p_her_bias = (sp.array(her_dict[h]['p_her']) - h / 100.0)
		bin_counts, bins, patch_list = pylab.hist(p_her_bias, range=(-0.25, 0.25), bins=25, alpha=0.6)
		max_bin_count = max(bin_counts)
		pylab.axvline(0.0, color='k', alpha=0.6, ls='-.')
		y_range = max_bin_count - 0
		pylab.axis([-0.25 - 0.025 * 0.5, 0.25 + 0.5 * 0.025, -0.025 * max_bin_count, max_bin_count * 1.025])
		pylab.xlabel('pseudo-heritability bias')
		pylab.ylabel('Counts')
		x.append(p_her_bias)
		pylab.savefig(png_file_name)

	#Box-plot version
	png_file_name = file_prefix + '_h%s_boxplot.png' % ('_'.join(map(str, her_dict.keys())))
	pylab.figure()
	pylab.boxplot(x)
	pylab.ylabel('psuedo-heritability bias')
	pylab.axhline(0.0, color='k', alpha=0.6, ls='-.')
	pylab.xticks(range(1, len(x) + 1), map(str, sorted(her_dict.keys())))
	pylab.xlabel('Heritability')
	pylab.savefig(png_file_name)


def plot_var(file_prefix, d, latent_variable, heritability, phen_models):
	"""
	Plots variance explained by the model for plus,  xor, and or.
	"""
	import pylab

	#Plot remaining heritability
	file_prefix += '_%s_%d' % (latent_variable, heritability)
	for am in ['Stepw_EX_Bonf', 'Stepw_EX_EBIC', 'Stepw_EX_MBIC']:
		png_file_name = file_prefix + '_%s_pher_boxplot.png' % (am)
		pylab.figure()
		rem_pher_list = []
		for m in phen_models:
			rem_pher_list.append(d[m][am]['rem_p_her'])
		pylab.boxplot(rem_pher_list)
		pylab.axhline(0.0, color='k', alpha=0.6, ls='-.')
		pylab.xticks(range(1, len(phen_models) + 1), phen_models)
		pylab.ylabel('Remaining pseudo-heritability in the model')
		pylab.savefig(png_file_name)

	#Plot % explained variance.
	for am in ['Stepw_EX_Bonf', 'Stepw_EX_EBIC', 'Stepw_EX_MBIC', 'Stepw_LM_Bonf', 'Stepw_LM_EBIC', 'Stepw_LM_MBIC']:
		png_file_name = file_prefix + '_%s_var_boxplot.png' % (am)
		pylab.figure()
		perc_var_list = []
		for m in phen_models:
			perc_var_list.append(d[m][am]['perc_var_expl'])
		pylab.boxplot(perc_var_list)
		pylab.axhline(heritability / 100.0, color='k', alpha=0.6, ls='-.')
		pylab.xticks(range(1, len(phen_models) + 1), phen_models)
		pylab.ylabel('Percentage of variance explained in the model')
		pylab.savefig(png_file_name)

	#Plot % explained phenot. variance
	for am in ['Stepw_EX_Bonf', 'Stepw_EX_EBIC', 'Stepw_EX_MBIC']:
		png_file_name = file_prefix + '_%s_phen_var_boxplot.png' % (am)
		pylab.figure()
		perc_var_list = []
		for m in phen_models:
			perc_var_list.append(d[m][am]['perc_phen_var_expl'])
		pylab.boxplot(perc_var_list)
		pylab.axhline(1.0, color='k', alpha=0.6, ls='-.')
		pylab.xticks(range(1, len(phen_models) + 1), phen_models)
		pylab.ylabel('Percentage of phenotypic variance explained in the model')
		pylab.savefig(png_file_name)
	#Plot % explained error variance
	for am in ['Stepw_EX_Bonf', 'Stepw_EX_EBIC', 'Stepw_EX_MBIC']:
		png_file_name = file_prefix + '_%s_err_var_boxplot.png' % (am)
		pylab.figure()
		perc_var_list = []
		for m in phen_models:
			perc_var_list.append(d[m][am]['perc_err_var_expl'])
		pylab.boxplot(perc_var_list)
		pylab.axhline(0.0, color='k', alpha=0.6, ls='-.')
		pylab.xticks(range(1, len(phen_models) + 1), phen_models)
		pylab.ylabel('Percentage of error variance explained in the model')
		pylab.savefig(png_file_name)




def __get_thresholds(min_thres=10, max_thres=1, num_thres=18):
	thres_step = (min_thres - max_thres) / float(num_thres)
	pval_thresholds = []
	for i in range(num_thres):
		pval_thresholds.append(max_thres + i * thres_step)
	return pval_thresholds

pval_thresholds = __get_thresholds()

window_sizes = [0, 1000, 5000, 10000, 25000, 50000, 100000]

def _update_stats_(gwa_res, c_chr, c_pos, l_chr=None, l_pos=None, significance_threshold=None, sign_res=None):
	"""
	Update result dictionary.
	"""
	res_dict = {}
	cpl = [(c_chr, c_pos)]#Causal chr_pos_list
	if l_chr != None:
		cpl.append((l_chr, l_pos))
	caus_indices = gwa_res.get_indices(cpl)
	gwa_res._rank_scores_()

	#Calculate KS and P-med..
	pvals = gwa_res.snp_results['scores'][:]
	res_dict['ks_stat'] = agr.calc_ks_stats(pvals)
	res_dict['med_pval'] = agr.calc_median(pvals)

	#Get causal p-values, and ranks
	res_dict['causal_pvals'] = [gwa_res.snp_results['scores'][i] for i in caus_indices]
	res_dict['causal_ranks'] = [gwa_res.ranks[i] for i in caus_indices]

	#Get significant chrom_pos_pvals..
	if (not sign_res) and significance_threshold :
		sign_res = gwa_res.filter_attr('scores', significance_threshold, reversed=True, return_clone=True)
		res_dict['sign_chr_pos'] = sign_res.get_chr_pos_score_list()
		res_dict['causal_dist_matrix'] = sign_res.get_distances(cpl)
	elif sign_res:
		res_dict['sign_chr_pos'] = sign_res.get_chr_pos_score_list()
		res_dict['causal_dist_matrix'] = sign_res.get_distances(cpl)


	#Of all SNPs ranked higher than the second causative... which is farthest from a nearest causative.
	dist = gwa_res.get_farthest_w_stronger_association(cpl)
	res_dict['dist_f_w_s_a'] = -1 if dist[0] > 0 else dist[1]

	#Perform power (sensitivity, TPR), FDR, FPR calculations..
	gwa_res.neg_log_trans()
	tprs_list = []
	fdrs_list = []
	for pval_thres in pval_thresholds:
		#Filter data
		gwa_res.filter_attr('scores', pval_thres)
		tprs, fdrs = gwa_res.get_power_analysis(cpl, window_sizes)
		tprs_list.append(tprs)
		fdrs_list.append(fdrs)
	res_dict['tprs'] = tprs_list #[p_valthreshold][window_size]
	res_dict['fdrs'] = fdrs_list #[p_valthreshold][window_size]
	return res_dict


def _update_sw_stats_(res_dict, step_info_list, opt_dict, c_chr, c_pos, l_chr=None, l_pos=None,
			significance_threshold=None, type='LM'):
	"""
	Update result dictionary for a stepwise result.
	"""
	res_dict['step_info_list'] = step_info_list
	cpl = [(c_chr, c_pos)]#Causal chr_pos_list
	if l_chr != None:
		cpl.append((l_chr, l_pos))

	for criteria in ['mbonf', 'mbics', 'ebics']:
		opt_i = opt_dict[criteria]
		d = {'opt_i':opt_i}
		si = step_info_list[opt_i]
		if criteria == 'mbonf':
			tprs_list = []
			fdrs_list = []
			t_opt_i_list = []
			num_steps = len(step_info_list) / 2
			max_cof_pvals = -sp.log10([step_info_list[i]['mbonf'] for i in range(2 * num_steps)])
			for pval_thres in pval_thresholds:
				t_opt_i = 0
				for i in range(num_steps + 1):
					if max_cof_pvals[i] >= pval_thres:
						t_opt_i = i
				for j in range(1, num_steps):
					i = 2 * num_steps - j
					if max_cof_pvals[i] >= pval_thres:
						if j > t_opt_i:
							t_opt_i = i
				if t_opt_i == 0:
					tprs = [-1 for ws in window_sizes]
					fdrs = [-1 for ws in window_sizes]
				else:
					t_si = step_info_list[t_opt_i]
					cpst = map(list, zip(*t_si['cofactors']))
					sign_res = gr.Result(scores=cpst[2], chromosomes=cpst[0], positions=cpst[1])
					tprs, fdrs = sign_res.get_power_analysis(cpl, window_sizes)
				tprs_list.append(tprs)
				fdrs_list.append(fdrs)
				t_opt_i_list.append(t_opt_i)
			d['tprs'] = tprs_list #[p_valthreshold][window_size]
			d['fdrs'] = fdrs_list #[p_valthreshold][window_size]
			d['t_opt_i_list'] = t_opt_i_list


		if opt_i == 0:
			#Set default values (Are these appropriate?)
			d['tprs'] = [-1 for ws in window_sizes]
			d['fdrs'] = [-1 for ws in window_sizes]
			d['sign_chr_pos'] = []
			d['causal_dist_matrix'] = []
		else:
			cpst = map(list, zip(*si['cofactors']))
			#Create a result object..
			sign_res = gr.Result(scores=cpst[2], chromosomes=cpst[0], positions=cpst[1])
			d['sign_chr_pos'] = sign_res.get_chr_pos_score_list()
			d['causal_dist_matrix'] = sign_res.get_distances(cpl)
			if criteria == 'mbonf':
				d['mbonf_tprs'], d['mbonf_fdrs'] = sign_res.get_power_analysis(cpl, window_sizes)
			else:
				d['tprs'], d['fdrs'] = sign_res.get_power_analysis(cpl, window_sizes)
		d['kolmogorov_smirnov'] = si['kolmogorov_smirnov']
		d['pval_median'] = si['pval_median']
		d['perc_var_expl'] = 1.0 - si['rss'] / step_info_list[0]['rss']
		d['num_cofactors'] = len(si['cofactors'])
		if type == 'EX':
			d['remaining_p_her'] = si['pseudo_heritability']
			d['perc_her_var_expl'] = (si['pseudo_heritability'] / step_info_list[0]['pseudo_heritability']) * d['perc_var_expl']
			d['perc_err_var_expl'] = ((1 - si['pseudo_heritability']) / (1 - step_info_list[0]['pseudo_heritability'])) * d['perc_var_expl']
		res_dict[criteria] = d




def run_analysis(sd, K, file_prefix, latent_var, heritability, phen_model, phen_index, phen_d,
		call_method_id, num_steps=10, pickle_results=True, save_plots=False):
	"""
	Perform the GWA mapping..
	using the different methods..
	
	Linear model, 
	Kruskal-Wallis
	EMMA
	
	Stepwise Linear Model (bonf. and ext. BIC)
	Stepwise EMMA (bonf. and ext. BIC)
	"""
	file_prefix += '_%d_%s_%s_%d' % (heritability, latent_var, phen_model, phen_index)

	pd = phen_d[latent_var][heritability][phen_model]


	result_dict = {}
	for mm in mapping_methods:
		result_dict[mm] = {}

	print "Loading SNPS dataset (again)"
	bonferroni_threshold = 1.0 / (20.0 * sd.num_snps())

	snps_list = sd.getSnps()
	phen_vals = pd['phenotypes'][phen_index]
	(c_snp, c_chr, c_pos, c_maf) = phen_d['snp_chr_pos_maf_list'][phen_index] #Causal SNP
	highlight_loci = [(c_chr, c_pos)]
	if latent_var == 'random_snp':
		(l_chr, l_pos, l_maf) = phen_d[latent_var]['latent_chr_pos_maf_list'][phen_index]
		highlight_loci.append((l_chr, l_pos))
	else:
		l_chr, l_pos = None, None

	print "Running Analysis"
	print 'Running KW'
	p_vals = util.kruskal_wallis(snps_list, phen_vals)['ps'].tolist()
	print len(p_vals)
	kw_res = gr.Result(snps_data=sd, scores=p_vals)
	if save_plots:
		kw_file_prefix = file_prefix + '_kw'
		kw_res.plot_manhattan(png_file=kw_file_prefix + '.png', highlight_loci=highlight_loci, neg_log_transform=True,
					plot_bonferroni=True)
		agr.plot_simple_qqplots(kw_file_prefix, [kw_res], result_labels=['Kruskal-Wallis'])


	print 'Updating stats for KW'
	result_dict['KW'] = _update_stats_(kw_res, c_chr, c_pos, l_chr, l_pos,
					significance_threshold=bonferroni_threshold)
	print 'Running SW LM'
	if save_plots:
		lm_file_prefix = file_prefix + '_lm'
	else:
		lm_file_prefix = None
	ret_dict = lm.lm_step_wise(phen_vals, sd, num_steps=num_steps, file_prefix=lm_file_prefix,
					with_qq_plots=save_plots, highlight_loci=highlight_loci)
	lm_step_info = ret_dict['step_info_list']
	lm_pvals = ret_dict['first_lm_res']['ps'].tolist()
	lm_opt_dict = ret_dict['opt_dict']
	lm_res = gr.Result(scores=lm_pvals, snps_data=sd)
	print 'Updating stats for LM'
	result_dict['LM'] = _update_stats_(lm_res, c_chr, c_pos, l_chr, l_pos,
					significance_threshold=bonferroni_threshold)
	print 'Updating stats for SW LM'
	_update_sw_stats_(result_dict['Stepw_LM'], lm_step_info, lm_opt_dict, c_chr, c_pos, l_chr, l_pos,
					significance_threshold=bonferroni_threshold)

	print 'Running SW EX'
	if save_plots:
		emmax_file_prefix = file_prefix + '_emmax'
	else:
		emmax_file_prefix = None
	ret_dict = lm.emmax_step_wise(phen_vals, K, sd, num_steps=num_steps, file_prefix=emmax_file_prefix,
					with_qq_plots=True, highlight_loci=highlight_loci)
	emmax_step_info = ret_dict['step_info_list']
	emmax_pvals = ret_dict['first_emmax_res']['ps'].tolist()
	emmax_opt_dict = ret_dict['opt_dict']
	emmax_res = gr.Result(scores=emmax_pvals, snps_data=sd)
	print 'Updating stats for EX'
	result_dict['EX'] = _update_stats_(emmax_res, c_chr, c_pos, l_chr, l_pos,
					significance_threshold=bonferroni_threshold)
	print 'Updating stats for SW EX'
	_update_sw_stats_(result_dict['Stepw_EX'], emmax_step_info, emmax_opt_dict, c_chr, c_pos, l_chr, l_pos,
					significance_threshold=bonferroni_threshold, type='EX')

	#Record trait pseudo-heritability:
	result_dict['p_her'] = emmax_step_info[0]['pseudo_heritability']

	if pickle_results == True:
		pickled_results_file = file_prefix + 'results.pickled'
		print 'Pickling result dict in file: %s' % pickled_results_file
		with open(pickled_results_file, 'wb') as f:
			cPickle.dump(result_dict, f, protocol=2)

	return result_dict




def _run_():
	p_dict, args = parse_parameters()
	print args
	file_prefix = env.env['results_dir'] + p_dict['run_id']

	if p_dict['sim_phen'] and p_dict['phen_file']: 	#Simulate phenotypes
		print 'Setting up phenotype simulations'
		sd = dp.load_250K_snps(p_dict['call_method_id'])
		simulate_phenotypes(p_dict['phen_file'], sd, debug_filter=p_dict['debug_filter'])

	elif p_dict['parallel']:
		#set up parallel runs
		if p_dict['phen_index'] == None:
			phen_d = load_phenotypes(p_dict['phen_file'])
			phenotypes = phen_d[p_dict['latent_variable']][p_dict['heritability']][p_dict['phenotype_model']]['phenotypes']
			start_i = 0
			end_i = len(phenotypes)
		else:
			start_i = p_dict['phen_index'][0]
			end_i = p_dict['phen_index'][-1]
		num_per_run = p_dict['number_per_run']
		for i in range(start_i, end_i, num_per_run):
			run_parallel(p_dict['run_id'], i, i + num_per_run - 1, p_dict['latent_variable'],
					p_dict['heritability'], p_dict['phenotype_model'],
					p_dict['phen_file'], p_dict['summarize'], p_dict['call_method_id'], p_dict['num_steps'],
					cluster='gmi')

	elif p_dict['phen_index']: #Run things..
		sd = dp.load_250K_snps(p_dict['call_method_id'], debug_filter=p_dict['debug_filter'])
		K = dp.load_kinship(p_dict['call_method_id'])
		results_list = []
		for pid in p_dict['phen_index']:
			result_dict = run_analysis(sd, K, file_prefix, p_dict['latent_variable'], p_dict['heritability'],
						p_dict['phenotype_model'], pid, load_phenotypes(p_dict['phen_file']),
						p_dict['call_method_id'], num_steps=p_dict['num_steps'])
			results_list.append(result_dict)
		#Save as pickled
	else:
		if p_dict['summarize']:
			results_list = []
			file_prefix = '/storage/mlt_results/' + p_dict['run_id']
			summary_dict = summarize_runs(file_prefix, p_dict['latent_variable'], p_dict['heritability'],
							p_dict['phenotype_model'], load_phenotypes(p_dict['phen_file']),
							index_list=p_dict['phen_index'])
			plot_file_prefix = '%s_%d_%s_%s' % (file_prefix, p_dict['heritability'], p_dict['latent_variable'],
								p_dict['phenotype_model'])
			plot_tprs_fdrs(plot_file_prefix, summary_dict)

		if p_dict['herit_plots'] != None:
			d = {}
			file_prefix = '/storage/mlt_results/' + p_dict['run_id']
			pd = load_phenotypes(p_dict['phen_file'])
			for her in p_dict['herit_plots']:
				d[her] = summarize_runs(file_prefix, p_dict['latent_variable'], her,
							p_dict['phenotype_model'], pd, index_list=p_dict['phen_index'])
			plot_herit_hist(file_prefix, d, p_dict['latent_variable'], p_dict['phenotype_model'])


		if p_dict['var_plots']:
			d = {}
			file_prefix = '/storage/mlt_results/' + p_dict['run_id']
			pd = load_phenotypes(p_dict['phen_file'])
			for mod in ['plus', 'or', 'xor']:
				d[mod] = summarize_runs(file_prefix, p_dict['latent_variable'], p_dict['heritability'], mod, pd, index_list=p_dict['phen_index'])
			plot_var(file_prefix, d, p_dict['latent_variable'], p_dict['heritability'], ['plus', 'or', 'xor'])



#def _run_vincent_scripts_():
#	type = sys.argv[1]
#	start = int(sys.argv[2])
#	end = int(sys.argv[3])
#	if type == 'add_emmax':
#		exec_str = 'sim2loci_add_fwdbwdemmax.sh'
#	elif type == 'add_lm':
#		exec_str = 'sim2loci_add_fwdbwdlm.sh'
#	elif type == 'or_emmax':
#		exec_str = 'sim2loci_or_fwdbwdemmax.sh'
#	elif type == 'or_lm':
#		exec_str = 'sim2loci_or_fwdbwdlm.sh'
#	elif type == 'xor_emmax':
#		exec_str = 'sim2loci_xor_fwdbwdemmax.sh'
#	elif type == 'xor_lm':
#		exec_str = 'sim2loci_xor_fwdbwdlm.sh'
#	for i in range(start, end + 1):
#		exec_st = 'qsub -q cmb -l walltime=6:00:00 -l mem=2950mb /home/cmbpanfs-01/bvilhjal/vincent/jobs/' + exec_str + ' -v VARIABLE=' + str(i)
#		print exec_st
#		os.system(exec_st)


if __name__ == '__main__':
	_run_()
#	sd = dp.load_250K_snps()
#	simulate_phenotypes(env.env['tmp_dir'] + 'simulated_phenotypes.pickled', sd)
	print "Done!!\n"
