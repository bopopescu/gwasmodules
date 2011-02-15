"""
Simulations for joint analysis of correlated traits.
"""
import scipy as sp
from scipy import stats
import random
import pdb
import sys
import dataParsers as dp
import phenotypeData as pd
import analyze_gwas_results as agr
import env
import os
import cPickle
import linear_models as lm
import math
import pylab
import time
import gwaResults as gr

#Simulate phenotypes..
def simulate_traits(sd, num_traits=200, heritabilities=[0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99], num_effects=100):
	"""
	Return a dictionary object of traits
	
	100 causal SNPs, with exponentially distributed effect sizes.
	"""

	snps = sd.getSnps()
	cpl = sd.getChrPosList()
	num_snps = len(snps)
	num_accessions = len(sd.accessions)

	sim_traits_dict = {}
	for h in heritabilities:
		trait_pairs = []
		sample_indices_list = []
		snp_effects_list = []
		num_non_overlapping_list = []
		rho_est_list = []
		trait1_perc_var_list = []
		trait2_perc_var_list = []
		for i in range(num_traits):
			if (i + 1) % int(num_traits / 100) == 0:
				sys.stdout.write('.')
				sys.stdout.flush()
				#print (i + 1) / int(num_traits / 100)
			#Simulate trait pair..

			num_non_overlapping = random.randint(0, num_effects)
			num_non_overlapping_list.append(num_non_overlapping)
			num_causatives = num_effects + num_non_overlapping
			sample_indices = random.sample(range(num_snps), num_causatives)
			chosen_snps = sp.array([snps[i] for i in sample_indices])
			c = sp.random.random_integers(0, 1, (num_causatives, 1))
			chosen_snps = sp.absolute(c - chosen_snps)
			exp_effects = stats.expon.rvs(scale=1, size=(num_causatives, 1))
			#exp_effects = stats.norm.rvs(scale=1, size=(num_causatives, 1))
			snp_effects = chosen_snps * exp_effects
			snp_effects1 = snp_effects[:num_effects]
			snp_effects2 = snp_effects[-num_effects:]
			trait1 = sp.sum(snp_effects1, 0)
			trait2 = sp.sum(snp_effects2, 0)

			gv = sp.var(trait1, ddof=1)
			error = stats.norm.rvs(0, 1, size=num_accessions)
			ev = sp.var(error, ddof=1)
			n_trait1 = trait1 + error * sp.sqrt(((1.0 - h) / h) * (gv / ev))
			trait1_perc_var_list.append(sp.var(snp_effects1, 1) / sp.var(n_trait1))
			n_trait1 = (n_trait1 - sp.mean(n_trait1)) / sp.std(n_trait1)
			gv = sp.var(trait2, ddof=1)
			error = stats.norm.rvs(0, 1, size=num_accessions)
			ev = sp.var(error, ddof=1)
			n_trait2 = trait2 + error * sp.sqrt(((1.0 - h) / h) * (gv / ev))
			trait2_perc_var_list.append(sp.var(snp_effects2, 1) / sp.var(n_trait2))
			n_trait2 = (n_trait2 - sp.mean(n_trait2)) / sp.std(n_trait2)
			trait_pairs.append((n_trait1, n_trait2))
			sample_indices_list.append(sample_indices)
			snp_effects_list.append(exp_effects)
			rho_est = sp.corrcoef(trait1, trait2)
			rho_est_list.append(rho_est)

			#Variance contributions.


		sim_traits_dict[h] = {'trait_pairs':trait_pairs, 'sample_indices_list':sample_indices_list,
				'num_non_overlapping_list':num_non_overlapping_list, 'snp_effects_list':snp_effects_list,
				'rho_est_list':rho_est_list, 'trait1_perc_var_list':trait1_perc_var_list,
				'trait2_perc_var_list':trait2_perc_var_list}
	return sim_traits_dict



def _load_sim_data_(use_pickle=True):
	sim_data_file = env.env['data_dir'] + 'corr_trait_sim_data.pickled'
	if use_pickle and os.path.isfile(sim_data_file):
		with open(sim_data_file) as f:
			sim_traits_dict = cPickle.load(f)
	else:
		import dataParsers as dp
		sd = dp.load_250K_snps()
		sim_traits_dict = simulate_traits(sd)
		if use_pickle:
			with open(sim_data_file, 'wb') as f:
				cPickle.dump(sim_traits_dict, f)
	return sim_traits_dict




def run_joint_analysis(start_i, stop_i, heritability, mac_threshold=15, debug_filter=1.0):
	sim_traits_dict = _load_sim_data_()[heritability]

	sd = dp.load_250K_snps(debug_filter=debug_filter)
	cpl = sd.getChrPosList() #For simulated SNPs lookup.
	sd.filter_mac_snps(mac_threshold)
	snps = sd.getSnps()
 	positions = sd.getPositions()
	chromosomes = sd.get_chr_list()
	K = lm.load_kinship_from_file(env.env['data_dir'] + 'kinship_matrix_cm72.pickled', sd.accessions)
	num_ecotypes = len(sd.accessions)

	file_prefix = env.env['results_dir'] + 'corr_sim_h%0.2f' % heritability
	res_dict = {'ks_stats_her':[], 'ks_stats_pher':[], 'ks_stats_indep':[],
			'pval_medians_her':[], 'pval_medians_pher':[], 'pval_medians_indep':[],
			'rho':[], 'rho_her':[], 'rho_pher':[],
			'p_hers':[], 'vgs':[], 'trait_corr':[]}
			#Add some summary statistics on the results.
	for i in range(start_i, stop_i):
		print i - start_i
		trait_pair = sim_traits_dict['trait_pairs'][i]
		sample_indices = sim_traits_dict['sample_indices_list'][i]
		num_non_overlapping = sim_traits_dict['num_non_overlapping_list'][i]
		snp_effects = sim_traits_dict['snp_effects_list'][i]
		rho = sim_traits_dict['rho_est_list'][i][0, 1]


		joint_phen_vals = trait_pair[0].tolist() + trait_pair[1].tolist()
		joint_ecotypes = sd.accessions + sd.accessions

		trait1_res = lm.emmax(snps, trait_pair[0], K)
		trait2_res = lm.emmax(snps, trait_pair[1], K)
		f_prefix = file_prefix + '_i%d' % i
		qq_file_name = f_prefix + '_qq_plot.png'
		log_qq_file_name = f_prefix + '_log_qq_plot.png'

		res1 = gr.Result(scores=trait1_res['ps'].tolist(), positions=positions, chromosomes=chromosomes)
		res2 = gr.Result(scores=trait2_res['ps'].tolist(), positions=positions, chromosomes=chromosomes)
		(areas, medians_h01) = agr.qq_plot({'EMMAX_trait1':res1, 'EMMAX_trait2':res2},
				1000, method_types=['emma', 'emma'], mapping_labels=['EMMAX_trait1', 'EMMAX_trait2'],
				phenName='marginal EMMAX', pngFile=qq_file_name)
		(ds_h01, areas, slopes) = agr.log_qq_plot({'EMMAX_trait1':res1, 'EMMAX_trait2':res2},
				1000, 7, method_types=['emma', 'emma'], mapping_labels=['EMMAX_trait1', 'EMMAX_trait2'],
				phenName='marginal EMMAX', pngFile=log_qq_file_name)

		png_file_name = f_prefix + '_manhattan_res1.png'
		res1.neg_log_trans()
		res1.plot_manhattan(png_file=png_file_name, plot_bonferroni=True)
		png_file_name = f_prefix + '_manhattan_res2.png'
		res2.neg_log_trans()
		res2.plot_manhattan(png_file=png_file_name, plot_bonferroni=True)

		gen_var_list = [trait1_res['vg'], trait2_res['vg']]
		err_var_list = [trait1_res['ve'], trait2_res['ve']]
		her_list = [trait1_res['pseudo_heritability'], trait2_res['pseudo_heritability']]
		res_dict['p_hers'].append(her_list)
		res_dict['vgs'].append(gen_var_list)

		min_ps = sp.minimum(trait1_res['ps'], trait2_res['ps'])




		#FINISH EMMA stuff..  from here on.. for comparison..

		#E in the model.
		joint_X = [[1] * (2 * num_ecotypes), [0] * num_ecotypes + [1] * num_ecotypes]
		E = sp.array([0] * num_ecotypes + [1] * num_ecotypes)

		#Get correlations between traits..
		corr_mat = sp.corrcoef(trait_pair[0], trait_pair[1])
		print her_list, gen_var_list, err_var_list, corr_mat[0, 1]
		res_dict['trait_corr'].append(corr_mat[0, 1])



		#Construct the full variance matrix
		V = sp.zeros((2 * num_ecotypes, 2 * num_ecotypes))
		V_2 = sp.zeros((2 * num_ecotypes, 2 * num_ecotypes))
		V_3 = sp.zeros((2 * num_ecotypes, 2 * num_ecotypes))
		V[0:num_ecotypes, 0:num_ecotypes] = gen_var_list[0] * K + err_var_list[0] * sp.eye(num_ecotypes)
		V[num_ecotypes:2 * num_ecotypes, num_ecotypes:2 * num_ecotypes] = \
				gen_var_list[1] * K + err_var_list[1] * sp.eye(num_ecotypes)
		V_2[0:num_ecotypes, 0:num_ecotypes] = V[0:num_ecotypes, 0:num_ecotypes]
		V_2[num_ecotypes:2 * num_ecotypes, num_ecotypes:2 * num_ecotypes] = \
				V[num_ecotypes:2 * num_ecotypes, num_ecotypes:2 * num_ecotypes]
		V_3[0:num_ecotypes, 0:num_ecotypes] = V[0:num_ecotypes, 0:num_ecotypes]
		V_3[num_ecotypes:2 * num_ecotypes, num_ecotypes:2 * num_ecotypes] = \
				V[num_ecotypes:2 * num_ecotypes, num_ecotypes:2 * num_ecotypes]

		rho_est = corr_mat[0, 1] / math.sqrt(her_list[0] * her_list[1])
		if rho_est > 1: rho_est = 1.0
		if rho_est < -1: rho_est = -1.0
		rho_est_2 = corr_mat[0, 1] / heritability
		if rho_est_2 > 1: rho_est_2 = 1.0
		if rho_est_2 < -1: rho_est_2 = -1.0

		print rho, rho_est, rho_est_2
		res_dict['rho'].append(rho)
		res_dict['rho_pher'].append(rho_est)
		res_dict['rho_her'].append(rho_est_2)
		V[0:num_ecotypes, num_ecotypes:2 * num_ecotypes] = \
				rho_est * math.sqrt(gen_var_list[0] * gen_var_list[1]) * K
		V[num_ecotypes:2 * num_ecotypes, 0:num_ecotypes] = V[0:num_ecotypes, num_ecotypes:2 * num_ecotypes]
		V_2[0:num_ecotypes, num_ecotypes:2 * num_ecotypes] = \
				rho_est_2 * math.sqrt(gen_var_list[0] * gen_var_list[1]) * K
		V_2[num_ecotypes:2 * num_ecotypes, 0:num_ecotypes] = V_2[0:num_ecotypes, num_ecotypes:2 * num_ecotypes]

		print 'Performing Cholesky decomposition'
		H_sqrt = lm.cholesky(V)
		H_sqrt_inv = sp.mat(H_sqrt.T).I
		H_sqrt_2 = lm.cholesky(V_2)
		H_sqrt_inv_2 = sp.mat(H_sqrt_2.T).I

		H_sqrt_3 = lm.cholesky(V_3) #Independence..
		H_sqrt_inv_3 = sp.mat(H_sqrt_3.T).I

		#pdb.set_trace()

		#Set up analysis..
		lmm = lm.LinearMixedModel(joint_phen_vals)
		lmm.set_factors(joint_X, False)

		#Doubling the SNPs!!!
		Z = sp.int16(sp.mat(joint_ecotypes).T == sp.mat(sd.accessions))

		#Running EMMAX(s)
		t0 = time.time()
		print 'Running EMMAX full model'
		res1 = lmm.emmax_full_model_gxe(snps, E, H_sqrt_inv, Z)
		t = time.time() - t0
		print 'EMMAX took %d minutes and %0.2f seconds.' % (int(t / 60), t % 60)
		t0 = time.time()
		print 'Running EMMAX full model'
		res2 = lmm.emmax_full_model_gxe(snps, E, H_sqrt_inv_2, Z)
		t = time.time() - t0
		print 'Second EMMAX took %d minutes and %0.2f seconds.' % (int(t / 60), t % 60)
		t0 = time.time()
		print 'Running EMMAX full model'
		res3 = lmm.emmax_full_model_gxe(snps, E, H_sqrt_inv_3, Z)
		t = time.time() - t0
		print 'Indep EMMAX took %d minutes and %0.2f seconds.' % (int(t / 60), t % 60)

		res1_h01 = gr.Result(scores=res1['ps_h01'].tolist(), positions=positions, chromosomes=chromosomes)
		res1_h02 = gr.Result(scores=res1['ps_h02'].tolist(), positions=positions, chromosomes=chromosomes)
		res1_h12 = gr.Result(scores=res1['ps_h12'].tolist(), positions=positions, chromosomes=chromosomes)

		res2_h01 = gr.Result(scores=res2['ps_h01'].tolist(), positions=positions, chromosomes=chromosomes)
		res2_h02 = gr.Result(scores=res2['ps_h02'].tolist(), positions=positions, chromosomes=chromosomes)
		res2_h12 = gr.Result(scores=res2['ps_h12'].tolist(), positions=positions, chromosomes=chromosomes)

		res3_h01 = gr.Result(scores=res3['ps_h01'].tolist(), positions=positions, chromosomes=chromosomes)
		res3_h02 = gr.Result(scores=res3['ps_h02'].tolist(), positions=positions, chromosomes=chromosomes)
		res3_h12 = gr.Result(scores=res3['ps_h12'].tolist(), positions=positions, chromosomes=chromosomes)

		qq_file_name = f_prefix + '_qq_plot_h01.png'
		log_qq_file_name = f_prefix + '_log_qq_plot_h01.png'
		(areas, medians_h01) = agr.qq_plot({'EMMAX_joint_her':res2_h01, 'EMMAX_joint_pher':res1_h01, 'EMMAX_joint_indep':res3_h01},
				1000, method_types=['emma', 'emma', 'emma'], mapping_labels=['EMMAX_joint_her', 'EMMAX_joint_pher', 'EMMAX_joint_indep'],
				phenName='G', pngFile=qq_file_name)
		(ds_h01, areas, slopes) = agr.log_qq_plot({'EMMAX_joint_her':res2_h01, 'EMMAX_joint_pher':res1_h01, 'EMMAX_joint_indep':res3_h01},
				1000, 7, method_types=['emma', 'emma', 'emma'], mapping_labels=['EMMAX_joint_her', 'EMMAX_joint_pher', 'EMMAX_joint_indep'],
				phenName='G', pngFile=log_qq_file_name)
		qq_file_name = f_prefix + '_qq_plot_h02.png'
		log_qq_file_name = f_prefix + '_log_qq_plot_h02.png'
		(areas, medians_h02) = agr.qq_plot({'EMMAX_joint_her':res2_h02, 'EMMAX_joint_pher':res1_h02, 'EMMAX_joint_indep':res3_h02},
				1000, method_types=['emma', 'emma', 'emma'], mapping_labels=['EMMAX_joint_her', 'EMMAX_joint_pher', 'EMMAX_joint_indep'],
				phenName='G+GxE', pngFile=qq_file_name)
		(ds_h02, areas, slopes) = agr.log_qq_plot({'EMMAX_joint_her':res2_h02, 'EMMAX_joint_pher':res1_h02, 'EMMAX_joint_indep':res3_h02},
				1000, 7, method_types=['emma', 'emma', 'emma'], mapping_labels=['EMMAX_joint_her', 'EMMAX_joint_pher', 'EMMAX_joint_indep'],
				phenName='G+GxE', pngFile=log_qq_file_name)
		qq_file_name = f_prefix + '_qq_plot_h12.png'
		log_qq_file_name = f_prefix + '_log_qq_plot_h12.png'
		(areas, medians_h12) = agr.qq_plot({'EMMAX_joint_her':res2_h12, 'EMMAX_joint_pher':res1_h12, 'EMMAX_joint_indep':res3_h12},
				1000, method_types=['emma', 'emma', 'emma'], mapping_labels=['EMMAX_joint_her', 'EMMAX_joint_pher', 'EMMAX_joint_indep'],
				phenName='GxE', pngFile=qq_file_name)
		(ds_h12, areas, slopes) = agr.log_qq_plot({'EMMAX_joint_her':res2_h12, 'EMMAX_joint_pher':res1_h12, 'EMMAX_joint_indep':res3_h12},
				1000, 7, method_types=['emma', 'emma', 'emma'], mapping_labels=['EMMAX_joint_her', 'EMMAX_joint_pher', 'EMMAX_joint_indep'],
				phenName='GxE', pngFile=log_qq_file_name)
#		pdb.set_trace()
#		qq_file_name = f_prefix + '_qq_full_model_plot.png'
#		log_qq_file_name = f_prefix + '_log_qq_full_model_plot.png'
#		(areas, medians) = agr.qq_plot({'EMMAX_joint_her_full_model':res2_h02, 'EMMAX_joint_pher_full_model':res1_h02}, 1000,
#				method_types=['emma', 'emma'],
#				mapping_labels=['EMMAX_joint_her_full_model', 'EMMAX_joint_pher_full_model'],
#				phenName='joint', pngFile=qq_file_name)
#		(ds, areas, slopes) = agr.log_qq_plot({'EMMAX_joint_her_full_model':res2_h02, 'EMMAX_joint_pher_full_model':res1_h02}, 1000,
#				7, method_types=['emma', 'emma'],
#				mapping_labels=['EMMAX_joint_her_full_model', 'EMMAX_joint_pher_full_model'],
#				phenName='joint', pngFile=log_qq_file_name)
		res_dict['ks_stats_her'].append([ds_h01[0], ds_h02[0], ds_h12[0]])
		res_dict['ks_stats_pher'].append([ds_h01[1], ds_h02[1], ds_h12[1]])
		res_dict['ks_stats_indep'].append([ds_h01[2], ds_h02[2], ds_h12[2]])
		res_dict['pval_medians_her'].append([medians_h01[0], medians_h02[0], medians_h12[0]])
		res_dict['pval_medians_pher'].append([medians_h01[1], medians_h02[1], medians_h12[1]])
		res_dict['pval_medians_indep'].append([medians_h01[2], medians_h02[2], medians_h12[2]])

		f_prefix = file_prefix + '_pher_i%d' % i
		png_file_name = f_prefix + '_manhattan_h01.png'
		res1_h01.neg_log_trans()
		res1_h01.plot_manhattan(png_file=png_file_name, plot_bonferroni=True)
		png_file_name = f_prefix + '_manhattan_h02.png'
		res1_h02.neg_log_trans()
		res1_h02.plot_manhattan(png_file=png_file_name, plot_bonferroni=True)
		png_file_name = f_prefix + '_manhattan_h12.png'
		res1_h12.neg_log_trans()
		res1_h12.plot_manhattan(png_file=png_file_name, plot_bonferroni=True)

		f_prefix_2 = file_prefix + '_her_i%d' % i
		png_file_name_2 = f_prefix_2 + '_manhattan_h01.png'
		res2_h01.neg_log_trans()
		res2_h01.plot_manhattan(png_file=png_file_name_2, plot_bonferroni=True)
		png_file_name_2 = f_prefix_2 + '_manhattan_h02.png'
		res2_h02.neg_log_trans()
		res2_h02.plot_manhattan(png_file=png_file_name_2, plot_bonferroni=True)
		png_file_name_2 = f_prefix_2 + '_manhattan_h12.png'
		res2_h12.neg_log_trans()
		res2_h12.plot_manhattan(png_file=png_file_name_2, plot_bonferroni=True)

		f_prefix_3 = file_prefix + '_indep_i%d' % i
		png_file_name_3 = f_prefix_3 + '_manhattan_h01.png'
		res3_h01.neg_log_trans()
		res3_h01.plot_manhattan(png_file=png_file_name_3, plot_bonferroni=True)
		png_file_name_3 = f_prefix_3 + '_manhattan_h02.png'
		res3_h02.neg_log_trans()
		res3_h02.plot_manhattan(png_file=png_file_name_3, plot_bonferroni=True)
		png_file_name_3 = f_prefix_3 + '_manhattan_h12.png'
		res3_h12.neg_log_trans()
		res3_h12.plot_manhattan(png_file=png_file_name_3, plot_bonferroni=True)

		#Process results.
		#Do we detect the causal loci?
#		trait1_only_causal =
#		common_causal =
#		trait1_perc_var = sim_traits_dict['trait1_perc_var_list'][i]
#		trait2_perc_var = sim_traits_dict['trait2_perc_var_list'][i]


		#Where are the peaks?
		#Etc..



	#Basic plots
	pylab.figure(figsize=(8, 6))
	pylab.plot(res_dict['rho'], res_dict['rho_pher'], 'g.', label='pseudo-herit. est.', alpha=0.6)
	pylab.plot(res_dict['rho'], res_dict['rho_her'], 'b.', label='herit. est.', alpha=0.6)
	min_rho = min(res_dict['rho'])
	max_rho = max(res_dict['rho'])
	range_rho = max_rho - min_rho
	max_y = max(max(res_dict['rho_pher']), max(res_dict['rho_her']))
	min_y = min(min(res_dict['rho_pher']), min(res_dict['rho_her']))
	range_y = max_y - min_y
	pylab.xlabel(r'$\rho$ (direct estimate)')
	pylab.ylabel(r'$\rho$ estimated using herit.')
	pylab.rcParams['legend.fontsize'] = 'small'
	pylab.legend(loc=2, numpoints=1)
	pylab.axis([min_rho - 0.05 * range_rho, max_rho + 0.05 * range_rho,
			min_y - 0.05 * range_y, max_y + 0.05 * range_y])
	lim = [min(pylab.xlim()[0], pylab.ylim()[0]), max(pylab.xlim()[1], pylab.ylim()[1])]
	pylab.plot(lim, lim, alpha=0.6, color='k', ls='--')
	pylab.savefig(env.env['tmp_dir'] + 'rho_est_plot_h%0.2f.png' % heritability)
	pylab.clf()

	rho_dist_pher = sp.absolute(sp.array(res_dict['rho']) - sp.array(res_dict['rho_pher']))
	rho_dist_her = sp.absolute(sp.array(res_dict['rho']) - sp.array(res_dict['rho_her']))
	min_rho = min(min(rho_dist_pher), min(rho_dist_her))
	max_rho = max(max(rho_dist_pher), max(rho_dist_her))
	range_rho = max_rho - min_rho
	pylab.plot(rho_dist_pher, res_dict['pval_medians_pher'], 'g.', label='pseudo-herit. est.', alpha=0.6)
	pylab.plot(rho_dist_her, res_dict['pval_medians_her'], 'b.', label='true herit. est.', alpha=0.6)
	pylab.plot(rho_dist_her, res_dict['pval_medians_indep'], 'b.', label='indep', alpha=0.6)
	max_y = max(max(res_dict['pval_medians_pher']), max(res_dict['pval_medians_her']), max(res_dict['pval_medians_indep']))
	min_y = min(min(res_dict['pval_medians_pher']), min(res_dict['pval_medians_her']), min(res_dict['pval_medians_indep']))
	range_y = max_y - min_y
	pylab.xlabel(r'$\rho - \rho_{\mathrm{est}}$')
	pylab.ylabel('Median p-value deviation.')
	pylab.rcParams['legend.fontsize'] = 'small'
	pylab.legend(loc=2, numpoints=1)
	pylab.axis([min_rho - 0.05 * range_rho, max_rho + 0.05 * range_rho,
			min_y - 0.05 * range_y, max_y + 0.05 * range_y])
	pylab.plot(pylab.xlim(), [0, 0], alpha=0.6, color='k', ls='--')
	pylab.savefig(env.env['tmp_dir'] + 'rho_median_pval_plot_h%0.2f.png' % heritability)
	pylab.clf()

	pylab.plot(rho_dist_pher, res_dict['ks_stats_pher'], 'g.', label='pseudo-herit. est.', alpha=0.6)
	pylab.plot(rho_dist_her, res_dict['ks_stats_her'], 'b.', label='herit. est.', alpha=0.6)
	pylab.plot(rho_dist_her, res_dict['ks_stats_indep'], 'b.', label='indep', alpha=0.6)
	max_y = max(max(res_dict['ks_stats_pher']), max(res_dict['ks_stats_her']), max(res_dict['ks_stats_indep']))
	min_y = min(min(res_dict['ks_stats_pher']), min(res_dict['ks_stats_her']), min(res_dict['ks_stats_indep']))
	range_y = max_y - min_y
	pylab.xlabel(r'$\rho - \rho_{\mathrm{est}}$')
	pylab.ylabel('Kolmogorov-Smirnov statistic.')
	pylab.rcParams['legend.fontsize'] = 'small'
	pylab.legend(loc=2, numpoints=1)
	pylab.axis([min_rho - 0.05 * range_rho, max_rho + 0.05 * range_rho,
			min_y - 0.05 * range_y, max_y + 0.05 * range_y])
	pylab.savefig(env.env['tmp_dir'] + 'rho_ks_stat_plot_h%0.2f.png' % heritability)



def _test_():
	#sim_traits_dict = _load_sim_data_()
	run_joint_analysis(0, 10, 0.99, debug_filter=0.05)

if __name__ == '__main__':
	_test_()


