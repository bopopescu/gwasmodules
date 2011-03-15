"""
Tests involving stepwise regressions, and model selection.

Option:

	-i ...			The run_id.  Note that if this option is used then no output file should be specified.
	-a ...			Mapping method, e.g emmax, lm
	-s			Collect results and plot things (Does not generate pvalue files...) 

	-t ...			What data set is used. (54 is default)
	-k ...			For running EMMAX.

	-n ...			Number of SNPs (phenotypes) per node, default is 100
	-f ...			Random fraction (given as parameters) of phenotypes will be used.

	-l ...			Type of latent variable: random_snps (default), etc..
	-e ...			Percentage of variance, due to error.
	-g ...			Fraction of error, du to kinship term.
	
	-m ...			How to generate the phenotypes: additive, xor, or

	--maf_filter=...	Generate the phenotypes on the fly, using SNPs with MAF greater than the given value.	
	--plot_pvals           	Plot Manhattan plots
	--phen_file=...		File where the phenotypes will be saved, and loaded from.
	--sim_phen		Simulate phenotype, write to phenotype file.
	--num_steps=...		Number of steps in the regression.

	-h			show this help

Additional parameters are ...

"""

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


#Parse Vincent's phenotype file...?  Why?
#Or create my phenotypes..
#3 models. 5 heritabilites.. 1000 phenotypes.

def parse_parameters():
	'Parse the parameters into a dict, etc.'
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)

	long_options_list = ["maf_filter=", 'plot_pvals', 'phen_file=', 'sim_phen', 'num_steps=']
	try:
		opts, args = getopt.getopt(sys.argv[1:], "a:i:t:k:n:f:l:e:g:m:hs", long_options_list)

	except:
		traceback.print_exc()
		print __doc__
		sys.exit(2)

	p_dict = {'number_per_run':100, 'filter':None, 'summarize':False, 'maf_filter':0,
		'latent_variable':'random_snps', 'phenotype_model':None, 'run_id':None,
		'mapping_method':'emmax', 'kinship_file':None, 'phenotype_error':0,
		'kinship_error':0, 'plot_pvals':False, 'call_method_id':54, 'phen_file':None,
		'num_steps':5}


	for opt, arg in opts:
		if opt in ("-h"):
			print __doc__
			return
		elif opt in ('-i'): p_dict['run_id'] = arg
		elif opt in ('-t'): p_dict['call_method_id'] = int(arg)
		elif opt in ('-n'): p_dict['number_per_run'] = int(arg)
		elif opt in ('-m'): p_dict['phenotype_model'] = int(arg)
		elif opt in ('-f'): p_dict['filter'] = float(arg)
		elif opt in ('-l'): p_dict['latent_variable'] = arg
		elif opt in ("-s"): p_dict['summarize'] = arg
		elif opt in ('-a'):p_dict['mapping_method'] = arg
		elif opt in ('-k'):p_dict['kinship_file'] = arg
		elif opt in ('-e'):p_dict['phenotype_error'] = float(arg)
		elif opt in ('-g'):p_dict['kinship_error'] = float(arg)
		elif opt in ("--phen_file"): p_dict['phen_file'] = arg
		elif opt in ("--plot_pvals"): p_dict['plot_pvals'] = True
		elif opt in ("--maf_filter"): p_dict['maf_filter'] = float(arg)
		elif opt in ("--sim_phen"): p_dict['sim_phen'] = True
		elif opt in ("--num_steps"): p_dict['num_steps'] = int(arg)
		else:
			print "Unkown option!!\n"
			print __doc__
			sys.exit(2)

	print p_dict, args
	return p_dict, args


def run_parallel(p_dict, phen_index=None, summary_run=False):
	"""
	Set up a parallel run on the cluster.
	"""
	raise NotImplementedError


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
	K = lm.load_kinship_from_file(env.env['data_dir'] + 'kinship_matrix_cm72.pickled') #All ecotypes
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
	heritabilities = [0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99]

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
						scalar = sp.sqrt((phen_var / error_var) * ((1 - h) / h))
						phenotype = phenotype + error_vector * scalar
						h_est = phen_var / sp.var(phenotype, ddof=1)
						h_est_snp1 = sp.corrcoef(snp, phenotype)[0, 1]
						h_est_snp2 = sp.corrcoef(latent_snp, phenotype)[0, 1]
						#print phen_model, latent_var, h, h_est, h_est_snp1 ** 2, h_est_snp2 ** 2
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


def summarize_runs(p_dict):
	raise NotImplementedError


def __get_thresholds(min_thres=10, max_thres=4, num_thres=12):
	thres_step = (min_thres - max_thres) / float(num_thres)
	pval_thresholds = []
	for i in range(num_thres):
		pval_thresholds.append(max_thres + i * thres_step)
	return pval_thresholds

__pval_thresholds = __get_thresholds()

__window_sizes = [0, 1000, 5000, 10000, 20000, 50000, 100000]

def _update_stats_(res_dict, gwa_res, c_chr, c_pos, l_chr=None, l_pos=None, significance_threshold=None, sign_res=None):
	"""
	Update result dictionary.
	"""
	cpl = [(c_chr, c_pos)]
	if l_chr != None:
		cpl.append((l_chr, l_pos))
	caus_indices = gwa_res.get_indices(cpl)
	gwa_res._rank_scores_()

	#Get causal p-values, and ranks
	caus_pvals = [gwa_res.snp_results['scores'][i] in caus_indices]
	caus_ranks = [gwa_res.ranks[i] in caus_indices]

	#Get significant chrom_pos_pvals..
	if (not sign_res) and significance_threshold :
		sign_res = gwa_res.filter_attr('scores', significance_threshold, reversed=True, return_clone=True)
		sign_chr_pos_score = sign_res.get_chr_pos_score_list()
		sign_dist = sign_res.get_distances(cpl)

	#Of all SNPs ranked higher than the second causative... which is farthest from a nearest causative.
	dist = gwa_res.get_farthest_w_stronger_association(cpl)
	dist = -1 if dist[0] > 0 else dist[1]

	#Perform power (sensitivity, TPR), FDR, FPR calculations..
	gwa_res.neg_log_trans()
	for pval_thres in __pval_thresholds:
		#Filter data
		gwa_res.filter_attr('scores', pval_thres)
		for window_size in __window_sizes:
			#calculate Power (sensitivity, TPR), FDR, FPR
			min_dist = gwa_res.get_min_distances(cpl)



def run_analysis(latent_var, heritability, phen_model, phen_index, phen_d):
	"""
	Perform the GWA mapping..
	using the different methods..
	
	Linear model, 
	Kruskal-Wallis
	EMMA
	
	Stepwise Linear Model (bonf. and ext. BIC)
	Stepwise EMMA (bonf. and ext. BIC)
	"""
	pd = phen_d[latent_var][heritability][phen_model]
	mapping_methods = ['LM', 'KW', 'EX', 'Stepw_LM_Bonf', 'Stepw_LM_extBIC', 'Stepw_EX_Bonf', 'Stepw_EX_extBIC'] #7 in total

	#What to save...
	#Distance of significant SNPs to all causative SNPs
	#Ranks and p-values of significant SNPs
	#Chromosome and position of all significant SNPs
	#Ranks and p-values of causal SNPs
	#Summarize the following:
	#	- Power (sensitivity, TPR), FDR, FPR, for different bonf. thresholds.. at 0, 5, 10, 20, and 100 kb window
	#	- Of all SNPs ranked higher than the second causative... which is farthest from a nearest causative.

	result_dict = {}
	for mm in mapping_methods:
		d = {'stats':{}, 'sign_chr_pos': [], 'sign_pvals':[], 'causal_dist_matrix':None}
		if latent_var in ['random_snp', 'random']:
			d['causal_pvals'] = []
			d['causal_ranks'] = []
		result_dict[mm] = d

	print "Loading SNPS dataset (again)"
	sd = dp.load_250K_snps(call_method_id=72)
	bonferroni_threshold = 1.0 / (20.0 * sd.num_snps())

	snps_list = sd.getSnps()
	phen_vals = pd['phenotypes'][phen_index]
	(c_snp, c_chr, c_pos, c_maf) = phen_d['snp_chr_pos_maf_list'][phen_index] #Causal SNP

	print "Running Analysis"
	print 'Running KW'
	p_vals = util.kruskal_wallis(snps_list, phen_vals)['ps']
	kw_res = gr.Result(snps_data=sd, scores=p_values)
	if latent_var == 'random_snp':
		(l_chr, l_pos, l_maf) = phen_d[latent_var]['latent_chr_pos_maf_list']
		_update_stats_(result_dict['KW'], kw_res, c_chr, c_pos, l_chr, l_pos,
				significance_threshold=bonferroni_threshold)
	else:
		_update_stats_(result_dict['KW'], kw_res, c_chr, c_pos, significance_threshold=bonferroni_threshold)


	print 'Running LM'
	#First step..
	print 'Running SW LM'
	#The other steps.. 10 steps..

	kinship_file = snpsdata.get_call_method_kinship_file(72)
	K = lm.load_kinship_from_file(kinship_file)
	print 'Running EX'
	#First step..
	print 'Running SW EX'
	#The other steps.. until pseudo-heritability is 0 or at most 10 steps.. 





        print 'Setting up global runs'
	for i, phenotype in enumerate(phenotypes):
		highlight_loci = [phen_chr_pos[i]]
		if latent_loci_snp_chr_pos_mafs:
			latent_snp = latent_loci_snp_chr_pos_mafs[i][0]
			latent_chr_pos = (latent_loci_snp_chr_pos_mafs[i][1], latent_loci_snp_chr_pos_mafs[i][2])
			highlight_loci.append(latent_chr_pos)
			print highlight_loci
                print "\nThe %d'th phenotype, variance=%0.3f :" % (i + phen_index, phenotype.var(ddof=1))
                #print phenotype
		sys.stdout.flush()

		if mapping_method == 'emmax':
			print 'Estimated heritability: %0.4f' % h_estimates[i]
 			file_prefix = '%s_%s_%d' % (outputFile, mapping_method, i + phen_index)
			results.append(emmax_step_wise(phenotypes, K, snps=snps_list, positions=snps_positions,
						chromosomes=chromosomes, file_prefix=file_prefix,
						num_steps=p_dict['num_steps']))
			print 'Performing the KS test'
			snps_sample = random.sample(snps_list, len(snps_list) / 5)
		       	kw_pvals = util.kruskal_wallis(snps_sample, phenotype, verbose=False)["ps"]
			perm_kw_pvals = gwa.get_perm_pvals(snps_sample, phenotype, snps_filter=0.25)
			ks_res = gwaResults._calcKS_(kw_pvals, perm_kw_pvals)
			ks_pval_statistic.append(ks_res['D'])
			print 'KS test finished, D=%0.4f' % (ks_res['D'])
			print 'Calculating the k_correlation.'
			k_corr = util.snp_kinship_correlation([causative_snps[i], latent_snp], k)['corr']
			k_correlation.append(k_corr)
			print 'Finished: corr=%0.4f' % k_corr
		elif mapping_method == 'lm':
			raise NotImplementedError
	print "len(results):", len(results)


	# ------------------------------------- SUMMARIZING RESULTS!!!!!! -------------------------------------
#	#Calculating various statistics.

	#Needs to be expanded for 
	min_dists_to_closest_causative = [[] for i in range(p_dict['num_steps'] + 1)]
	min_dists_to_second_causative = [[] for i in range(p_dict['num_steps'])]
	for i, phenotype in enumerate(phenotypes):
		step_infos = results[i]
		causatives = [phen_chr_pos[i], (latent_loci_snp_chr_pos_mafs[i][1], latent_loci_snp_chr_pos_mafs[i][2])]
		for j, si in enumerate(step_infos):
			distances = map(tuple,
				list(sp.absolute(sp.array(si['min_pval_chr_pos']) - sp.array(causatives))))
			min_dist = min(distances)
			min_dists_to_closest_causative[j].append(min_dist)
			distances.remove(min_dist)
			min_dists_to_second_causative[j].append(min(distances))







#	#Distances to min p-val.
#	distances_to_min_pval = []
#	#Ranks of true SNPs
#	true_snps_ranks = []
#	min_pvals = []
#	obs_pvals = []
#	rank_statistics = []  #For the max ranks histogram!
#	sign_statistics = []  #For the max significant histogram!
#	sign_fractions = []
#	if latent_corr:
#		latent_distances_to_min_pval = []
#		latent_snps_ranks = []
#		latent_pvals = []
#		second_dist_to_min_pval = []
#	print len(results), len(phen_mafs)
#	for i in range(len(results)):
#		sys.stdout.write("Summarizing results for the " + str(i) + "'th phenotype.\n")
#		sys.stdout.flush()
#		pvals = results[i]
#		sign_count = 0
#		tot_count = 0
#		pval_chr_pos = []
#		sign_pval_chr_pos = []
#		for j, pval in enumerate(pvals):
#			tot_count += 1
#			if pval < (0.05 / len(snps_list)):
#				sign_count += 1
#				sign_pval_chr_pos.append((pval, chr_pos_list[j][0], chr_pos_list[j][1]))
#			pval_chr_pos.append((pval, chr_pos_list[j][0], chr_pos_list[j][1]))
#		sign_fractions.append(sign_count / float(tot_count))
#		pval_chr_pos.sort()
#
#		min_pval = min(pvals)
#		min_pvals.append(min_pval)
#		min_pvals_indices = list(*sp.where(pvals == min_pval))
#		if len(min_pvals_indices) > 1:
#			print 'Found two minimum p-values:'
#			for j in min_pvals_indices:
#				print phen_chr_pos[i]
#			print 'Using the first one.'
#		j = min_pvals_indices[0]
#
#		#chromosome and position of the causative SNP
#		(phen_chr, phen_pos) = phen_chr_pos[i]
#
#		#chromosome and position of the most significant SNP
#		(chr, pos) = chr_pos_list[j]
#
#		if chr == phen_chr:
#			distances_to_min_pval.append(abs(int(phen_pos) - int(pos)))
#		else:
#			distances_to_min_pval.append(-abs(int(phen_chr) - int(chr)))
#		try:
#        		snp_i = chr_pos_list.index((phen_chr, phen_pos))
#        		obs_pval = pvals[snp_i]
#        		obs_pvals.append(obs_pval)
#        		r_list = util.getRanks(pvals)
#        		snp_rank = int(r_list[snp_i])
#        		true_snps_ranks.append(snp_rank)
#
#		except Exception, err_str:
#        		print "Didn't find causative SNP:", phen_pos
#        		print err_str
#        		true_snps_ranks.append(-1)
#
#		if latent_corr:
#			(latent_snp, latent_chr, latent_pos, latent_maf) = latent_loci_snp_chr_pos_mafs[i]
#			if chr == latent_chr:
#				latent_distances_to_min_pval.append(abs(int(latent_pos) - int(pos)))
#			else:
#				latent_distances_to_min_pval.append(-abs(int(latent_chr) - int(chr)))
#			try:
#				latent_snp_i = chr_pos_list.index((latent_chr, latent_pos))
#				obs_pval = pvals[latent_snp_i]
#				latent_pvals.append(obs_pval)
#				latent_rank = int(r_list[latent_snp_i])
#				latent_snps_ranks.append(latent_rank)
#			except Exception, err_str:
#				print "Didn't find latent causative SNP:", phen_pos
#				print err_str
#				true_snps_ranks.append(-1)
#
#			# Which is the closest one to the min. pval?
#			if latent_distances_to_min_pval[-1] < 0 and distances_to_min_pval[-1] < 0:
#				second_dist_to_min_pval.append(-3)  #Both causative loci on diff. chromosomes from the min pval.
#			elif latent_distances_to_min_pval[-1] < 0:
#				second_chr = latent_chr
#				second_pos = latent_pos
#			elif distances_to_min_pval[-1] < 0:
#				second_chr = phen_chr
#				second_pos = phen_pos
#			else:
#				second_dist_to_min_pval.append(-2)  #Both causative loci on same chromosome.
#
#			# New min pval!!
#			elif chr != first_chr:
#				if chr != second_chr:
#					second_dist_to_min_pval.append(-1)
#				else:
#					second_dist_to_min_pval.append(abs(int(second_pos)-pos))
#			else:
#				
#
#
#
#		#Calculating the dist to the farthest SNP with rank greater or equal to the second ranked causative SNP.
#		if latent_corr:
#			print latent_corr
#			max_rank = int(math.ceil(max(latent_rank, snp_rank)) + 0.001)
#		else:
#			max_rank = int(math.ceil(snp_rank) + 0.001)
#
#		im_pval_chr_pos = pval_chr_pos[:max_rank]
#		max_dist = 0
#		diff_chr = False
#		for (im_pval, im_ch, im_pos) in im_pval_chr_pos:
#			if latent_corr and im_ch == phen_chr and im_ch == latent_chr:
#				max_dist = max(max_dist, min(abs(im_pos - phen_pos), abs(im_pos - latent_pos)))
#			elif im_ch == phen_chr:
#				max_dist = max(max_dist, abs(im_pos - phen_pos))
#			elif latent_corr and im_ch == latent_chr:
#				max_dist = max(max_dist, abs(im_pos - latent_pos))
#			else:
#				diff_chr = True
#				break
#		if diff_chr:
#			rank_statistics.append(-1)
#		else:
#			rank_statistics.append(max_dist)
#
#		#Calculating the max distance to a significant p-val.
#		if sign_pval_chr_pos:
#			max_dist = 0
#			diff_chr = False
#			for (s_pval, s_ch, s_pos) in sign_pval_chr_pos:
#				if s_ch == phen_chr:
#					if latent_corr and phen_chr == latent_chr:
#						max_dist = max(max_dist, min(abs(s_pos - phen_pos), abs(phen_pos - latent_pos)))
#					else:
#						max_dist = max(max_dist, abs(s_pos - phen_pos))
#				elif latent_corr and s_ch == latent_chr:
#					max_dist = max(max_dist, abs(phen_pos - latent_pos))
#				else:
#					diff_chr = True
#					break
#			if diff_chr:
#				sign_statistics.append(-1)
#			else:
#				sign_statistics.append(max_dist)
##					if max_dist == 0:
##						print "sign_pval_chr_pos:",sign_pval_chr_pos
##						print "pvals[snp_i]:",pvals[snp_i]
##						print "pvals[latent_snp_i]:",pvals[latent_snp_i]
##						print "phen_pos,phen_chr:",phen_pos,phen_chr
##						print "latent_pos,latent_chr:",latent_pos,latent_chr
#		else:
#			sign_statistics.append(-2)
#
#
#
#
#
#
#
#		gc.collect()  #Calling garbage collector, in an attempt to clean up memory..
#	print "Distances:", distances_to_min_pval
#	print "Ranks:", true_snps_ranks
#	print "Obs. pvals:", obs_pvals
#	print "Min pvals:", min_pvals
#	print "Significant fractions:", sign_fractions
#	print "Significant pvals statistics:", sign_statistics
#	if latent_corr:
#		print "Distances:", latent_distances_to_min_pval
#		print "Latent ranks:", latent_snps_ranks
#		print "Latent obs. pvals:", latent_pvals
#		print "Intermittent rank statistics:", rank_statistics
#
#
#	#Discarding the missing p-values due to local association mapping. 
#	#results = map(list,zip(*results))
#	#print "len(results):",len(results), len(chr_pos)
#
#	#Write results to file
#
#	if not noPvals and local:
#		print "Writing p-values to file:"
#		new_results = []
#		for pvals in results: #For all phenotypes (results)
#			chr_pos_pvals = []
#			for i in range(0, len(chr_pos_list)): #For all SNPs
#				pval = pvals[i]
#				if pval != 2 and pval <= pvalueThreshold:
#					(chr, pos) = chr_pos_list[i]
#					chr_pos_pvals.append((chr, pos, pval))
#			new_results.append(chr_pos_pvals)
#		results = new_results
#		l = [phen_positions, results]
#		pvalFile = outputFile + ".pvals"
#		print "Writing p-values to file:", pvalFile
#		f = open(pvalFile, "w")
#		cPickle.dump(l, f)
#		f.close()
#
#	p_d_r_dict = {}
#	if latent_corr:
#		p_d_r_dict['latent_pvals'] = latent_pvals
#		p_d_r_dict['latent_snps_ranks'] = latent_snps_ranks
#		p_d_r_dict['latent_distances_to_min_pval'] = latent_distances_to_min_pval
#		p_d_r_dict['latent_loci_snp_chr_pos_mafs'] = latent_loci_snp_chr_pos_mafs
#
#	p_d_r_dict['phen_positions'] = phen_positions
#	p_d_r_dict['distances_to_min_pval'] = distances_to_min_pval
#	p_d_r_dict['true_snps_ranks'] = true_snps_ranks
#	p_d_r_dict['phen_mafs'] = phen_mafs
#	p_d_r_dict['obs_pvals'] = obs_pvals
#	p_d_r_dict['min_pvals'] = min_pvals
#	p_d_r_dict['sign_fractions'] = sign_fractions
#	p_d_r_dict['rank_statistics'] = rank_statistics
#	p_d_r_dict['sign_statistics'] = sign_statistics
#	if mapping_method == 'emmax':
#		p_d_r_dict['pseudo_heritabilities'] = pseudo_heritabilities
#		p_d_r_dict['h_estimates'] = h_estimates
#		p_d_r_dict['k_correlation'] = k_correlation
#		p_d_r_dict['ks_statistic'] = ks_pval_statistic
#
#	filename = outputFile + ".stats"
#	print "Writing results to file:", filename
#	f = open(filename, "w")
#	cPickle.dump(p_d_r_dict, f)
#	f.close



def _run_():
	p_dict, args = parse_parameters()

	phenotype_models = range(1, 5)
	p_dict['snps_dataset'] = None

	if p_dict['summarize']: #Set up a summary run..

		if not p_dict['phenotype_model']: #Then run on cluster, summarize for all phenotype models
			for pm in phenotype_models:
				p_dict['phenotype_model'] = pm
				run_parallel(p_dict, summary_run=True)
		else:
			summarise_runs(p_dict)
		return  #Exit...

	else: #Then load SNPs data
		snps_data_file = snpsdata.get_call_method_dataset_file(p_dict['call_method_id'], binary_format=True)
		p_dict['snps_dataset'] = dataParsers.parse_binary_snp_data(snpsDataFile)

	if p_dict['phen_file']:
		phen_file = p_dict['phen_file']
	else:
		phen_file = env.env['tmp_dir'] + p_dict['run_id'] + ".phen"

	if p_dict['phenotype_model']:  #If phenotype model is specified, then use that, otherwise run all..
		phenotype_models = [p_dict['phenotype_model']]


	if p_dict['sim_phen']:
		phen_d = simulate_phenotype(phen_file, p_dict)
	else:
		phen_d = load_phenotypes(phen_file)

	if args:
		phen_index = int(args[0])
	else:
		#Running on the cluster..
		phen_indices = range(0, num_phenotypes, numberPerRun)
		for pm in phenotype_models:
			print "Submitting jobs for phenotype model:", pm
			p_dict['phenotype_model'] = pm
			for phen_index in phen_indices:
				run_parallel(p_dict, phen_index=phen_index)
		return #Exiting

	print "phen_index:", phen_index
	print "\nStarting analysis now!\n"
	print "run_id:", run_id
	run_analysis(phen_index, phen_d, p_dict, phenotype_models)



def _run_vincent_scripts_():
	type = sys.argv[1]
	start = int(sys.argv[2])
	end = int(sys.argv[3])
	if type == 'add_emmax':
		exec_str = 'sim2loci_add_fwdbwdemmax.sh'
	elif type == 'add_lm':
		exec_str = 'sim2loci_add_fwdbwdlm.sh'
	elif type == 'or_emmax':
		exec_str = 'sim2loci_or_fwdbwdemmax.sh'
	elif type == 'or_lm':
		exec_str = 'sim2loci_or_fwdbwdlm.sh'
	elif type == 'xor_emmax':
		exec_str = 'sim2loci_xor_fwdbwdemmax.sh'
	elif type == 'xor_lm':
		exec_str = 'sim2loci_xor_fwdbwdlm.sh'
	for i in range(start, end + 1):
		exec_st = 'qsub -q cmb -l walltime=6:00:00 -l mem=2950mb /home/cmbpanfs-01/bvilhjal/vincent/jobs/' + exec_str + ' -v VARIABLE=' + str(i)
		print exec_st
		os.system(exec_st)


if __name__ == '__main__':
	#_run_()
	sd = dp.load_250K_snps()
	simulate_phenotypes(env.env['tmp_dir'] + 'simulated_phenotypes.pickled', sd)
	print "Done!!\n"
