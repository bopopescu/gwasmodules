"""
Basic analysis of RNA data.

Bjarni J. Vilhjalmsson (bjarni.vilhjalmsson@gmi.oeaw.ac.at)
"""

#TODO
# - Run calculations on cluster (using full sequence data):
# - Draw cool plots..
# - Using pylons, create interactive web based plots.

import phenotypeData as pd
import dataParsers as dp
import gwaResults as gr
import linear_models as lm
import scipy as sp
from env import *
import sys
import cPickle
import util
import gwaResults
import pylab
import pdb
import itertools as it
import math
import analyze_gwas_results as agr
#import ipdb

#Annoyingly ad hoc constants
near_const_filter = 20
phen_file_prefix = env['phen_dir'] + 'rna_seq_081411'
#phen_file_prefix = env['phen_dir'] + 'rna_seq_061611'


def run_parallel(x_start_i, x_stop_i, temperature, call_method_id, cluster='gmi', run_id='rna_seq'):
	"""
	If no mapping_method, then analysis run is set up.
	"""
	job_id = '%s_%d_%d_%s' % (run_id, x_start_i, x_stop_i, temperature)
	file_prefix = env['results_dir'] + job_id

	#Cluster specific parameters	
	if cluster == 'gmi': #GMI cluster.  
		shstr = '#!/bin/bash\n'
		shstr += '#$ -S /bin/bash\n'
		shstr += '#$ -N %s\n' % job_id
		#shstr += '#$ -o %s_job_$JOB_ID.out\n' % file_prefix
		#shstr += '#$ -e %s_job_$JOB_ID.err\n' % file_prefix
		shstr += '#$ -o %s_job.out\n' % file_prefix
		shstr += '#$ -e %s_job.err\n' % file_prefix
		shstr += 'source /etc/modules-env.sh\n'
		shstr += 'module load scipy/GotoBLAS2/0.9.0\n'
		shstr += 'module load matplotlib/1.0.0\n'
		shstr += 'module load mysqldb/1.2.3\n'
		shstr += 'module load h5py/2.0.0\n'
		shstr += 'export GOTO_NUM_THREADS=1\n'

	elif cluster == 'usc':  #USC cluster.
		shstr = "#!/bin/csh\n"
		shstr += "#PBS -l walltime=%s \n" % '72:00:00'
		shstr += "#PBS -l mem=%s \n" % '2950mb'
		shstr += "#PBS -q cmb\n"
		shstr += "#PBS -N p%s \n" % job_id

	shstr += "python %sanalyze_rna_seq.py %d %d %s %d %s" % \
			(env['script_dir'], x_start_i, x_stop_i, temperature, call_method_id, run_id)

	#shstr += "> " + file_prefix + "_job.out) >& " + file_prefix + "_job.err\n"
	print '\n', shstr, '\n'
	script_file_name = run_id + ".sh"
	f = open(script_file_name, 'w')
	f.write(shstr)
	f.close()

	#Execute qsub script
	os.system("qsub " + script_file_name)






def summarize_stepwise(summary_dict, gene, step_info_list, opt_dict):
	#Store results for PPAs, MBONF, EBIC
	sw_d = {}
	for criteria in ['ebics', 'mbonf', 'min_cof_ppa']:
		i_opt = opt_dict[criteria]
		step_info = step_info_list[i_opt]
		cof_list = step_info['cofactors']
		ppa_cof_list = step_info['ppa_cofactors']
		cofactors = [(chrom, pos) for chrom, pos, pval in cof_list]
		cof_pvals = [pval for chrom, pos, pval in cof_list]
		cof_ppas = [ppa for chrom, pos, ppa in ppa_cof_list]

		# How close to the gene are the selected loci.
		da = []
		for chrom, pos in cofactors:
			if gene.chromosome == chrom:
				if gene.startPos < pos < gene.endPos:
					da.append(0)
				else:
					da.append(min(abs(pos - gene.startPos), abs(pos - gene.endPos)))
			else:
				da.append(-1)

		# Number of selected loci near the gene (bins)
		bin_distances = [100000, 50000, 25000, 10000, 5000, 1000, 0]
		bin_counts = [da.count(-1)]
		bin_count = 0
		for d in da:
			if d > bin_distances[0]: bin_count += 1
		bin_counts.append(bin_count)
		for bin_dist in bin_distances:
			bin_count = 0
			for d in da:
				if d == -1: continue
				elif d <= bin_dist: bin_count += 1
			bin_counts.append(bin_count)

		# Percentage of variance (error and genetic) explained by the selected loci (bin it as well)
		pass

		d = {'cofactors':cofactors, 'cof_pvals':cof_pvals, 'cof_ppas':cof_ppas, 'cof_gene_dist':da,
		'bin_counts':bin_counts, 'i_opt':i_opt}
		#print d
		sw_d[criteria] = d
		sw_d['step_info_list'] = step_info_list
	summary_dict['SW'] = sw_d




def run_gwas(file_prefix, phen_file, start_i, stop_i, temperature, mac_threshold=15, filter_threshold=0.02,
		call_method_id=79, data_format='diploid_int', debug_filter=1.0, near_const_filter=20):
	"""
	GWAS
	"""
	phed = pd.parse_phenotype_file(phen_file, with_db_ids=False)  #load phenotype file
	phed.filter_near_const_phens(near_const_filter)
	phed.convert_to_averages()
	num_traits = phed.num_traits()
	pids = phed.phen_ids[start_i :stop_i]
	sd = dp.load_snps_call_method(call_method_id=call_method_id, data_format=data_format, debug_filter=debug_filter)
	indices_to_keep = sd.coordinate_w_phenotype_data(phed, 1, coord_phen=False)  #All phenotypes are ordered the same way, so we pick the first one.
	phed.filter_ecotypes(indices_to_keep, pids=pids)
	print len(sd.accessions)
	K = sd.get_ibs_kinship_matrix()
	#K = dp.load_kinship(call_method_id=call_method_id, data_format=data_format, sd=sd, method='ibs')

	sd.filter_mac_snps(mac_threshold)
	snps = sd.getSnps()
	positions = sd.getPositions()
	chromosomes = sd.get_chr_list()
	r = sd.get_mafs()
	macs = r['mafs']
	mafs = r['marfs']

	print 'In total there are %d SNPs to be mapped.' % len(snps)
	gene_dict = dp.parse_tair_gff_file()#_load_genes_list_('rna_seq_031311_%sC' % temperature)
	for i, pid in enumerate(pids):
		if not pid in phed.phen_ids: continue
		gene_tair_id = phed.get_name(pid)
#		exons = []
#		for isoform in d:
#			for exon in isoform['exons']:
#				exons.append((d['chromosome'], exon['start_pos'], exon['end_pos']))


		d = gene_dict[gene_tair_id]
		gene_strand = d['strand']
		try:
			chrom = int(d['chromosome'])
		except Exception:
			raise
		gene = gwaResults.Gene(chromosome=int(d['chromosome']), startPos=d['start_pos'],
				endPos=d['end_pos'], name=gene_tair_id, description=None, dbRef=gene_tair_id,
				tairID=gene_tair_id)
		print i, pid, gene
		curr_file_prefix = '%s_mac%d_pid%d_%s' % (file_prefix, mac_threshold, pid, gene_tair_id)

		trans_type, shapiro_pval = phed.most_normal_transformation(pid)
		print 'Most normal transformation was: %s' % trans_type
		#trans_type = 'None'
		summary_dict = {'transformation_type':trans_type, 'transformation_shapiro_pval':shapiro_pval}
		#summary_dict = {'transformation_type':trans_type, 'transformation_shapiro_pval':0}


		print'Applying Kruskal-Wallis'
		phen_vals = phed.get_values(pid)
		res = util.kruskal_wallis(snps, phen_vals)
		pvals = res['ps'].tolist()
		kw_res = gr.Result(scores=pvals, macs=macs, mafs=mafs, positions=positions, chromosomes=chromosomes)
		print 'Summarizing KW'
		summary_dict['KW'] = kw_res.get_gene_analysis(gene)
		summary_dict['KW']['kolmogorov_smirnov'] = agr.calc_ks_stats(res['ps'])
		summary_dict['KW']['pval_median'] = agr.calc_median(res['ps'])


		print 'Applying LM'
		res = lm.linear_model(snps, phen_vals)
		pvals = res['ps'].tolist()
		perc_var_expl = res['var_perc'].tolist()
		lm_res = gr.Result(scores=pvals, macs=macs, mafs=mafs, positions=positions, chromosomes=chromosomes,
				perc_var_expl=perc_var_expl)
		print 'Summarizing LM'
		summary_dict['LM'] = lm_res.get_gene_analysis(gene)
		summary_dict['LM']['kolmogorov_smirnov'] = agr.calc_ks_stats(res['ps'])
		summary_dict['LM']['pval_median'] = agr.calc_median(res['ps'])


		print 'Applying EX Stepwise'
		snp_priors = sd.get_cand_genes_snp_priors([gene])
		ex_sw_res = lm.emmax_step_wise(phen_vals, K, macs=macs, mafs=mafs, positions=positions,
					chromosomes=chromosomes, snps=snps, num_steps=5, cand_gene_list=[gene],
					with_qq_plots=False, log_qq_max_val=6.0, save_pvals=True, snp_priors=snp_priors)
		print 'Summarizing the step-wise mixed model'
		pvals = ex_sw_res['first_emmax_res']['ps'].tolist()
		perc_var_expl = ex_sw_res['first_emmax_res']['var_perc'].tolist()
		ex_res = gr.Result(scores=pvals, macs=macs, mafs=mafs, positions=positions, chromosomes=chromosomes,
				perc_var_expl=perc_var_expl)
		summary_dict['EX'] = ex_res.get_gene_analysis(gene)
		summary_dict['pseudo_heritability'] = ex_sw_res['step_info_list'][0]['pseudo_heritability']
		summary_dict['EX']['kolmogorov_smirnov'] = agr.calc_ks_stats(ex_sw_res['first_emmax_res']['ps'])
		summary_dict['EX']['pval_median'] = agr.calc_median(ex_sw_res['first_emmax_res']['ps'])

		#Does the linear mixed model fit the data better?
		summary_dict['MM_LRT'] = lm.mm_lrt_test(phen_vals, K)

		#FINISH summarizing the stepwise!!!
		summarize_stepwise(summary_dict, gene, ex_sw_res['step_info_list'], ex_sw_res['opt_dict'])

		cvt_dict = {'radius':{}, 'tss_upstream':{}}
		print 'Comparing cis vs. trans kinship'
		#Check 1 mb, 200kb, 100kb, 50kb, 20kb, 10kb, 2kb, 0kb
		for radius in [500000, 100000, 50000, 25000, 10000, 5000, 1000, 0]:
			print radius
			r_start_pos = max(gene.startPos - radius, 0)
			r_end_pos = gene.endPos + radius
			d = sd.get_region_split_kinships([(chrom, r_start_pos, r_end_pos)],
							kinship_method='ibs', global_kinship=K)
			reg_k = d['regional_k']
			glob_k = d['global_k']
			if reg_k != None:
				cvt_dict['radius'][radius] = lm.local_vs_global_mm(phen_vals, reg_k, glob_k, K)
			else:
				cvt_dict['radius'][radius] = None
			print cvt_dict['radius'][radius]

		#Check TSS, 100kb, 50kb,25kb, 10kb,5kb,0kb, (all upstream)
		for dist in [200000, 100000, 50000, 25000, 10000, 5000, 1000]:
			print dist, gene_strand
			if gene_strand == '+':
				r_start_pos = max(gene.startPos - dist, 0)
				r_end_pos = gene.startPos
			else:
				r_start_pos = gene.endPos
				r_end_pos = gene.endPos + dist
			d = sd.get_region_split_kinships([(chrom, r_start_pos, r_end_pos)],
							kinship_method='ibs', global_kinship=K)
			reg_k = d['regional_k']
			glob_k = d['global_k']
			if reg_k != None:
				cvt_dict['tss_upstream'][dist] = lm.local_vs_global_mm(phen_vals, reg_k, glob_k, K)
			else:
				cvt_dict['tss_upstream'][dist] = None
			print cvt_dict['tss_upstream'][dist]

		summary_dict['CVT'] = cvt_dict

		#Write info to file..
		cPickle.dump(summary_dict, open(curr_file_prefix + '_info.pickled', 'w'), protocol=2)

		f_prefix = curr_file_prefix + '_hist'
		phed.plot_histogram(pid, title='Gene expressions for %s' % gene_tair_id,
				png_file=f_prefix + '.png', p_her=summary_dict['pseudo_heritability'],
				x_label='RNA seq expression levels (%s transformed)' % trans_type)
		#Plot GWAs...
		for res, method_name in [(kw_res, 'KW'), (lm_res, 'LM'), (ex_res, 'EX')]:
			res.filter_percentile(filter_threshold, reversed=True)
			res.write_to_file('%s_%s_.pvals' % (curr_file_prefix, method_name), only_pickled=True)
			if ex_res.min_score() < 10e-10:
				#print [cg.tairID for cg in cgs]
				f_prefix = '%s_%s_manhattan' % (curr_file_prefix, method_name)
				res.plot_manhattan(png_file=f_prefix + '.png', percentile=0, cand_genes=[gene],
						plot_bonferroni=True, neg_log_transform=True)




def plot(temperature=10, call_method_id=75, mapping_method='EX', mac_threshold=15, min_score=5,
		near_const_filter=20, data_format='binary', plot_data=True):
	#Load in chromosome dict..

	#file_prefix = '/srv/lab/data/rna_seq_062911/%dC/cm_%d/' % (temperature, call_method_id)
	file_prefix = '/srv/lab/data/rna_seq_083011/%dC/cm_%d/' % (temperature, call_method_id)

	results_dict_file = '%sresults_%s_mac%d.pickled' % (file_prefix, mapping_method, mac_threshold)
	res_dict = cPickle.load(open(results_dict_file))

	phen_file = '%s_%dC.csv' % (phen_file_prefix, temperature)
	phed = pd.parse_phenotype_file(phen_file, with_db_ids=False)  #load phenotype file
	phed.filter_near_const_phens(near_const_filter)
	phed.convert_to_averages()
	num_traits = phed.num_traits()
	pids = phed.phen_ids
	sd = dp.load_snps_call_method(call_method_id=call_method_id, data_format=data_format, debug_filter=0.01)
	indices_to_keep = sd.coordinate_w_phenotype_data(phed, 1, coord_phen=False)  #All phenotypes are ordered the same way, so we pick the first one.
	phed.filter_ecotypes(indices_to_keep, pids=pids)

	chrom_dict = {}
	for x_chrom in [1, 2, 3, 4, 5]:
		for y_chrom in [1, 2, 3, 4, 5]:
			chrom_dict[(x_chrom, y_chrom)] = {'scores':[], 'x_positions':[], 'y_positions':[],
								'tair_ids':[], 'r2':[], 'mac':[]}
	scores = []
	for x_chrom, x_pos in res_dict:
		d = res_dict[(x_chrom, x_pos)]
		tair_id = d['tair_id']
		for y_chrom in [1, 2, 3, 4, 5]:
			cps_d = d['chrom_pos_score'][y_chrom]
			for i in range(len(cps_d['scores'])):
				s = cps_d['scores'][i]
				if s > min_score:
					if s > 25:
						s = 25
					scores.append(s)
					chrom_dict[(x_chrom, y_chrom)]['scores'].append(s)
					chrom_dict[(x_chrom, y_chrom)]['tair_ids'].append(tair_id)
					chrom_dict[(x_chrom, y_chrom)]['x_positions'].append(x_pos)
					chrom_dict[(x_chrom, y_chrom)]['y_positions'].append(cps_d['positions'][i])

	#Write chrom_dict to file..
	if not plot_data:
		for x_chrom in [1, 2, 3, 4, 5]:
			for y_chrom in [1, 2, 3, 4, 5]:
				file_name = file_prefix + 'result_plots/pvalues_chrom%d_chrom%d_%s_min%d.txt' % (x_chrom, y_chrom, mapping_method, min_score)
				print 'Writing to file:', file_name
				with open(file_name, 'w') as f:
					d = chrom_dict[(x_chrom, y_chrom)]
					f.write('x_position, y_position, score, tair_id\n')
					l = zip(d['x_positions'], d['y_positions'], d['scores'], d['tair_ids'])
					l.sort()
					for t in l:
						f.write('%d,%d,%f,%s\n' % t)



	chrom_sizes = [30425061, 19694800, 23456476, 18578714, 26974904]
	cum_chrom_sizes = [sum(chrom_sizes[:i]) for i in range(5)]
	tot_num_bases = float(sum(chrom_sizes))
	rel_chrom_sizes = map(lambda x: 0.925 * (x / tot_num_bases), chrom_sizes)
	rel_cum_chrom_sizes = map(lambda x: 0.925 * (x / tot_num_bases), cum_chrom_sizes)
	for i in range(5):
		rel_cum_chrom_sizes[i] = rel_cum_chrom_sizes[i] + 0.02 + 0.01 * i

	chromosome_ends = {1:30.425061, 2:19.694800, 3:23.456476, 4:18.578714, 5:26.974904}
	print rel_chrom_sizes, rel_cum_chrom_sizes

	#Filter data..

	#Now plot data!!
	if plot_data:
		alpha = 0.8
		linewidths = 0
		vmin = min_score
		f = pylab.figure(figsize=(40, 35))
		chromosomes = [1, 2, 3, 4, 5]
		plot_file_name = file_prefix + 'result_plots/pvalues_%s_min%d.png' % (mapping_method, min_score)
		label = '$-log_{10}$(p-value)'
		vmax = max(scores)

		for yi, chr2 in enumerate(chromosomes):
			for xi, chr1 in enumerate(chromosomes):

				l = chrom_dict[(chr1, chr2)]['scores']
				if len(l) == 0:
					continue
				ax = f.add_axes([0.96 * (rel_cum_chrom_sizes[xi] + 0.01), rel_cum_chrom_sizes[yi] - 0.02,
						0.96 * (rel_chrom_sizes[xi]), rel_chrom_sizes[yi] ])
				ax.spines['right'].set_visible(False)
				ax.spines['bottom'].set_visible(False)
				#ax.tick_params(fontsize='x-large')
				if xi > 0:
					ax.spines['left'].set_visible(False)
					ax.yaxis.set_visible(False)
				else:
					ax.yaxis.set_ticks_position('left')
					ax.set_ylabel('Chromosome %d (Mb)' % chr2, fontsize='x-large')
				if yi < 4:
					ax.spines['top'].set_visible(False)
					ax.xaxis.set_visible(False)
				else:
					ax.xaxis.set_ticks_position('top')
					ax.xaxis.set_label_position('top')
					ax.set_xlabel('Chromosome %d (Mb)' % chr1, fontsize='x-large')
					#ax.set_xlabel('Chromosome %d' % chr1)

				#l = -sp.log10(l)
				#l = l.tolist()
				l_zxy = zip(l, chrom_dict[(chr1, chr2)]['x_positions'],
					chrom_dict[(chr1, chr2)]['y_positions'])
				l_zxy.sort()
				l = map(list, zip(*l_zxy))
				zs = l[0]
				xs = map(lambda x: x / 1000000.0, l[1])
				ys = map(lambda x: x / 1000000.0, l[2])

				scatter_plot = ax.scatter(xs, ys, c=zs, alpha=alpha, linewidths=linewidths, vmin=vmin,
							vmax=vmax)
				ax.axis([-0.025 * chromosome_ends[chr1], 1.025 * chromosome_ends[chr1],
					- 0.025 * chromosome_ends[chr2], 1.025 * chromosome_ends[chr2]])

		cax = f.add_axes([0.965, 0.7, 0.01, 0.2])
		cb = pylab.colorbar(scatter_plot, cax=cax)
		cb.set_label(label, fontsize='xx-large')
		#cb.set_tick_params(fontsize='x-large')
		f.text(0.005, 0.47, 'Associated SNP position', size='xx-large', rotation='vertical')
		f.text(0.47, 0.988, 'Expressed gene position', size='xx-large')
		print 'Saving figure:', plot_file_name
		f.savefig(plot_file_name, format='png')








def load_and_plot_info_files(call_method_id=75, temperature=10, mac_threshold=15, debug_filter=1,
			near_const_filter=20, data_format='binary'):
	import random

	phen_file = '%s_%dC.csv' % (phen_file_prefix, temperature)
	phed = pd.parse_phenotype_file(phen_file, with_db_ids=False)  #load phenotype file
	phed.filter_near_const_phens(near_const_filter)
	phed.convert_to_averages()
	num_traits = phed.num_traits()
	pids = phed.phen_ids
	sd = dp.load_snps_call_method(call_method_id=call_method_id, data_format=data_format, debug_filter=0.01)
	indices_to_keep = sd.coordinate_w_phenotype_data(phed, 1, coord_phen=False)  #All phenotypes are ordered the same way, so we pick the first one.
	phed.filter_ecotypes(indices_to_keep, pids=pids)

	print 'Loading the gene annotation dictionary'
	gene_dict = dp.parse_tair_gff_file()
	run_id = 'rna_seq'
	#run_id = 'rs_%d' % call_method_id


	file_prefix = '/srv/lab/data/rna_seq_083011/%dC/cm_%d/' % (temperature, call_method_id)


	num_genes = 0

	radii = [500000, 100000, 50000, 25000, 10000, 5000, 1000, 0]
	tss_dists = [200000, 100000, 50000, 25000, 10000, 5000, 1000]
	cvt_summary_dict = {'radius':{'avg_cis_trans_var_ratio':[0.0 for r in radii],
					'avg_cis_herit':[0.0 for r in radii],
					'avg_trans_herit':[0.0 for r in radii],
					'counts':[0.0 for td in radii]},
				'radius_herit':{'avg_cis_trans_var_ratio':[0.0 for r in radii],
					'avg_cis_herit':[0.0 for r in radii],
					'avg_trans_herit':[0.0 for r in radii],
					'counts':[0.0 for td in radii]},
				'tss_dist':{'avg_cis_trans_var_ratio':[0.0 for td in tss_dists],
					'avg_cis_herit':[0.0 for td in tss_dists],
					'avg_trans_herit':[0.0 for td in tss_dists],
					'counts':[0.0 for td in tss_dists]}}
	heritabilities = []
	transformations = []
	shapiro_wilk_pvals = []
	tair_ids = []
	pval_infl_dict = {}
	dist_min_pval_dict = {}
	distance_bins = [(0, 5000), (0, 10000), (0, 25000), (0, 50000), (0, 100000), (1, -1), (6, -1)]
	radius_bins = [0, 1000, 5000, 10000, 25000, 50000, 100000]
	bonf_sign_bin_dict = {}
	res_dict = {}
	sign_count = {}
	for mm in ['EX', 'LM', 'KW']:
		pval_infl_dict[mm] = {'kolmogorov_smirnov':[], 'median_pvals':[]}
		dist_min_pval_dict[mm] = {}
		for bin in distance_bins:
			dist_min_pval_dict[mm][bin] = 0
		bonf_sign_bin_dict[mm] = {}
		for bin in radius_bins:
			bonf_sign_bin_dict[mm][bin] = {'count':0.0, 'total':0.0}
		sign_count[mm] = 0

	cofactor_count_dict = {}
	for criteria in ['ebics', 'mbonf', 'min_cof_ppa']:
		cofactor_count_dict[criteria] = {'num_cofactor_list':[], 'bin_counts':sp.zeros(9),
						'num_cis_cofactor_list':[], 'num_found':0}

	pickle_file_dict = {}
	for mm in ['EX', 'LM', 'KW']:
		pickle_file_dict[mm] = {}
		pickle_file_dict[mm]['file_name'] = '%sresults_%s_mac%d.pickled' % (file_prefix, mm, mac_threshold)
		pickle_file_dict[mm]['res_dict'] = {}

	pids = phed.get_pids()
	for i, pid in enumerate(pids):
		tair_id = phed.get_name(pid)
		chrom = int(tair_id[2])
		curr_file_prefix = '%schr_%d/rna_seq_%s_%dC_mac%d_pid%d_%s' % \
					(file_prefix, chrom, run_id, temperature, mac_threshold, pid, tair_id)
		info_file_name = '%s_info.pickled' % curr_file_prefix
		for mm in ['EX', 'LM', 'KW']:
			res_dict[mm] = '%s_%s_.pvals' % (curr_file_prefix, mm)
		if random.random() > debug_filter:
			continue
		if os.path.isfile(info_file_name) and os.path.isfile(res_dict['EX'] + ".pickled") \
				and os.path.isfile(res_dict['LM'] + ".pickled") and os.path.isfile(res_dict['KW'] + ".pickled"):
			print 'Loading info file: %s' % info_file_name
			num_genes += 1
			info_dict = cPickle.load(open(info_file_name)) #Loading the info dict
			for mm in ['EX', 'LM', 'KW']:
				res_dict[mm] = gr.Result(res_dict[mm]) #Loading the result

			#Saving some basic statistics
			transformations.append(info_dict['transformation_type'])
			shapiro_wilk_pvals.append(info_dict['transformation_shapiro_pval'])
			heritabilities.append(info_dict['pseudo_heritability'])

			#cis vs. trans stuff
			cvt_dict = info_dict['CVT']
			for r_i, r in enumerate(radii):
				if cvt_dict['radius'][r] != None:
					pvg = cvt_dict['radius'][r]['perc_var1']
					pvl = cvt_dict['radius'][r]['perc_var2']
					herit = cvt_dict['radius'][r]['pseudo_heritability1']
					cvt_summary_dict['radius']['avg_cis_trans_var_ratio'][r_i] += pvl / (pvl + pvg)
					cvt_summary_dict['radius']['avg_cis_herit'][r_i] += pvl * herit
					cvt_summary_dict['radius']['avg_trans_herit'][r_i] += pvg * herit
					cvt_summary_dict['radius']['counts'][r_i] += 1.0

			for r_i, r in enumerate(radii):
				if cvt_dict['radius'][r] != None:
					herit = cvt_dict['radius'][r]['pseudo_heritability1']
					if herit > 0.05:
						pvg = cvt_dict['radius'][r]['perc_var1']
						pvl = cvt_dict['radius'][r]['perc_var2']
						cvt_summary_dict['radius_herit']['avg_cis_trans_var_ratio'][r_i] += pvl / (pvl + pvg)
						cvt_summary_dict['radius_herit']['avg_cis_herit'][r_i] += pvl * herit
						cvt_summary_dict['radius_herit']['avg_trans_herit'][r_i] += pvg * herit
						cvt_summary_dict['radius_herit']['counts'][r_i] += 1.0

			for td_i, td in enumerate(tss_dists):
				if cvt_dict['tss_upstream'][td] != None:
					pvg = cvt_dict['tss_upstream'][td]['perc_var1']
					pvl = cvt_dict['tss_upstream'][td]['perc_var2']
					herit = cvt_dict['tss_upstream'][td]['pseudo_heritability1']
					cvt_summary_dict['tss_dist']['avg_cis_trans_var_ratio'][td_i] += pvl / (pvl + pvg)
					cvt_summary_dict['tss_dist']['avg_cis_herit'][td_i] += pvl * herit
					cvt_summary_dict['tss_dist']['avg_trans_herit'][td_i] += pvg * herit
					cvt_summary_dict['tss_dist']['counts'][td_i] += 1.0



			tair_ids.append(tair_id)
			for mm in ['EX', 'LM', 'KW']:
				pval_infl_dict[mm]['kolmogorov_smirnov'].append(info_dict[mm]['kolmogorov_smirnov']['D'])
				pval_infl_dict[mm]['median_pvals'].append(info_dict[mm]['pval_median'])
				dist_min_pval = tuple(info_dict[mm]['dist_to_min_pval'])
				if res_dict[mm].min_score() < 1 / (20.0 * res_dict[mm].num_scores()):
					sign_count[mm] += 1
					for bin in distance_bins:
						if dist_min_pval <= bin:
							dist_min_pval_dict[mm][bin] += 1
							break

				for bin in radius_bins:
					pval = info_dict[mm]['bin_dict'][bin]['min_pval']
					num_snps = info_dict[mm]['bin_dict'][bin]['num_snps']
					if num_snps > 0:
						bonf_sign_bin_dict[mm][bin]['total'] += 1
						if pval < 1.0 / (20 * num_snps):
							bonf_sign_bin_dict[mm][bin]['count'] += 1

			#Stepwise stuff 
			for criteria in ['ebics', 'mbonf', 'min_cof_ppa']:
				num_cofactors = len(info_dict['SW'][criteria]['cofactors'])
				cofactor_count_dict[criteria]['num_cofactor_list'].append(num_cofactors)
				if num_cofactors > 0:
					cofactor_count_dict[criteria]['num_found'] += 1
					cofactor_count_dict[criteria]['bin_counts'] += sp.array(info_dict['SW'][criteria]['bin_counts'])
					cofactor_count_dict[criteria]['num_cis_cofactor_list'].append(info_dict['SW'][criteria]['bin_counts'][2])


			#Pre-process the results..
			for mm in ['EX', 'LM', 'KW']:
				res = res_dict[mm]
				#Trim results
				res.neg_log_trans()
				if mm == 'EX':
					res.filter_attr('scores', 2.5) #Filter everything below 10^-2.5
				else:
					res.filter_attr('scores', 4) #Filter everything below 10^-4
				if res.num_scores() == 0:
					print "Skipping file since nothing is below 10^-5"
					continue

				gene_d = gene_dict[tair_id]
				avg_g_pos = (gene_d['start_pos'] + gene_d['end_pos']) / 2.0
				chrom = int(gene_d['chromosome']) #Current gene chromosome


				#Prepare for plotting results.. x,y style, where gene is x, and y is p-values
				chrom_pos_score_dict = res.get_chrom_score_pos_dict()

				dist_dict = {}
				for score_threshold in [5, 6, 7]: #negative log10 thresholds.
					if len(res.snp_results['scores']) == 0:
						dist_dict[score_threshold] = -2 #No results
					else:
						res.filter_attr('scores', score_threshold)
						if len(res.snp_results['scores']) == 0:
							dist_dict[score_threshold] = -2 #No results
						else:
							cps_dict = res.get_chrom_score_pos_dict()
							pos_list = cps_dict[chrom]['positions']
							if len(pos_list) > 0:
								distances = sp.absolute(sp.array(pos_list) - avg_g_pos)
								d_i = sp.argmin(distances)
								dist_dict[score_threshold] = distances[d_i] #Min distance.
							else:
								dist_dict[score_threshold] = -1 #Different chromosome

				pickle_file_dict[mm]['res_dict'][(chrom, avg_g_pos)] = {'tair_id':tair_id,
							'chrom_pos_score':chrom_pos_score_dict, 'dist_dict':dist_dict,
							'pid':pid}
				print dist_dict
		else:
			print "Didn't find file: %s or %s" % (info_file_name, res_dict['EX'] + ".pickled")

	for mm in ['EX', 'LM', 'KW']:
		cPickle.dump(pickle_file_dict[mm]['res_dict'], open(pickle_file_dict[mm]['file_name'], 'wb'), protocol=2)


	for r_i, r in enumerate(radii):
		r_counts = cvt_summary_dict['radius']['counts'][r_i]
		cvt_summary_dict['radius']['avg_cis_trans_var_ratio'][r_i] = \
			cvt_summary_dict['radius']['avg_cis_trans_var_ratio'][r_i] / r_counts
		cvt_summary_dict['radius']['avg_cis_herit'][r_i] = \
			cvt_summary_dict['radius']['avg_cis_herit'][r_i] / r_counts
		cvt_summary_dict['radius']['avg_trans_herit'][r_i] = \
			cvt_summary_dict['radius']['avg_trans_herit'][r_i] / r_counts


	for r_i, r in enumerate(radii):
		r_counts = cvt_summary_dict['radius_herit']['counts'][r_i]
		cvt_summary_dict['radius_herit']['avg_cis_trans_var_ratio'][r_i] = \
			cvt_summary_dict['radius_herit']['avg_cis_trans_var_ratio'][r_i] / r_counts
		cvt_summary_dict['radius_herit']['avg_cis_herit'][r_i] = \
			cvt_summary_dict['radius_herit']['avg_cis_herit'][r_i] / r_counts
		cvt_summary_dict['radius_herit']['avg_trans_herit'][r_i] = \
			cvt_summary_dict['radius_herit']['avg_trans_herit'][r_i] / r_counts


	for td_i, td in enumerate(tss_dists):
		td_counts = cvt_summary_dict['tss_dist']['counts'][td_i]
		cvt_summary_dict['tss_dist']['avg_cis_trans_var_ratio'][td_i] = \
			cvt_summary_dict['tss_dist']['avg_cis_trans_var_ratio'][td_i] / td_counts
		cvt_summary_dict['tss_dist']['avg_cis_herit'][td_i] = \
			cvt_summary_dict['tss_dist']['avg_cis_herit'][td_i] / td_counts
		cvt_summary_dict['tss_dist']['avg_trans_herit'][td_i] = \
			cvt_summary_dict['tss_dist']['avg_trans_herit'][td_i] / td_counts



	results_prefix = env['results_dir'] + 'RNAseq_summary_%dC_cm%d' % (temperature, call_method_id)

	pylab.figure()
	pylab.plot(cvt_summary_dict['radius']['avg_cis_trans_var_ratio'])
	pylab.ylabel('Avg. perc. of cis genetic var.')
	pylab.xlabel('Dist. from gene (kb)')
	pylab.xticks(range(8), [500, 100, 50, 25, 10, 5, 1, 0])
	pylab.savefig(results_prefix + '_avg_perc_cis_gen_var_rad.png')
	pylab.clf()

	pylab.figure()
	pylab.plot(cvt_summary_dict['tss_dist']['avg_cis_trans_var_ratio'])
	pylab.ylabel('Avg. perc. of cis genetic var.')
	pylab.xlabel('Dist. upstream from gene TSS (kb)')
	pylab.xticks(range(7), [200, 100, 50, 25, 10, 5, 1])
	pylab.savefig(results_prefix + '_avg_perc_cis_gen_var_td.png')
	pylab.clf()

#	pylab.figure()
#	pylab.plot(cvt_summary_dict['tss_dist']['avg_cis_herit'])
#	pylab.ylabel('Avg. cis heritability')
#	pylab.xlabel('Dist. upstream from gene TSS (kb)')
#	pylab.xticks(range(7), [200, 100, 50, 25, 10, 5, 1])
#	pylab.savefig(results_prefix + 'avg_cis_herit_td.png')
#	pylab.clf()
#
#
#	pylab.figure()
#	pylab.plot(cvt_summary_dict['tss_dist']['avg_trans_herit'])
#	pylab.ylabel('Avg. remaining heritability')
#	pylab.xlabel('Dist. upstream from gene TSS (kb)')
#	pylab.xticks(range(7), [200, 100, 50, 25, 10, 5, 1])
#	pylab.savefig(results_prefix + 'avg_trans_herit_td.png')
#	pylab.clf()


#	pylab.figure()
#	pylab.plot(cvt_summary_dict['radius']['avg_trans_herit'])
#	pylab.ylabel('Avg. remaining heritability')
#	pylab.xlabel('Dist. from gene (kb)')
#	pylab.xticks(range(8), [500, 100, 50, 25, 10, 5, 1, 0])
#	pylab.savefig(results_prefix + 'avg_trans_herit_rad.png')
#	pylab.clf()
#
#	pylab.figure()
#	pylab.plot(cvt_summary_dict['radius']['avg_cis_herit'])
#	pylab.ylabel('Avg. cis heritability')
#	pylab.xlabel('Dist. from gene (kb)')
#	pylab.xticks(range(8), [500, 100, 50, 25, 10, 5, 1, 0])
#	pylab.savefig(results_prefix + 'avg_cis_herit_rad.png')
#	pylab.clf()

	tot_herit = sp.array(cvt_summary_dict['radius']['avg_cis_herit']) + \
		sp.array(cvt_summary_dict['radius']['avg_trans_herit'])
	cis_herit = sp.array(cvt_summary_dict['radius']['avg_cis_herit'])
	pylab.figure(figsize=(10, 6))
	pylab.axes([0.06, 0.08, 0.92, 0.90])
	pylab.fill_between([0, 7], 0, 1, color='#DD3333', alpha=0.8, label='Error')
	pylab.fill_between(sp.arange(8), 0, tot_herit, color='#22CC44', alpha=0.8, label='Heritable variance')
	pylab.fill_between(sp.arange(8), 0, cis_herit, color='#2255AA', \
				alpha=0.8, label='Heritable variance (cis)')
	pylab.ylabel('Average partition of variance')
	pylab.xlabel('Dist. from gene (kb)')
	pylab.xticks(range(8), [500, 100, 50, 25, 10, 5, 1, 0])
	pylab.legend(loc=1, ncol=3, shadow=True)
	pylab.axis([0, 7, 0, 1])
	pylab.savefig(results_prefix + 'avg_herit_rad.png')

	tot_herit = sp.array(cvt_summary_dict['radius_herit']['avg_cis_herit']) + \
		sp.array(cvt_summary_dict['radius_herit']['avg_trans_herit'])
	cis_herit = sp.array(cvt_summary_dict['radius_herit']['avg_cis_herit'])
	pylab.figure(figsize=(10, 6))
	pylab.axes([0.06, 0.08, 0.92, 0.90])
	pylab.fill_between([0, 7], 0, 1, color='#DD3333', alpha=0.8, label='Error')
	pylab.fill_between(sp.arange(8), 0, tot_herit, color='#22CC44', alpha=0.8, label='Heritable variance')
	pylab.fill_between(sp.arange(8), 0, cis_herit, color='#2255AA', \
				alpha=0.8, label='Heritable variance (cis)')
	pylab.ylabel('Average partition of variance')
	pylab.xlabel('Dist. from gene (kb)')
	pylab.xticks(range(8), [500, 100, 50, 25, 10, 5, 1, 0])
	pylab.legend(loc=1, ncol=3, shadow=True)
	pylab.axis([0, 7, 0, 1])
	pylab.savefig(results_prefix + 'avg_herit_2_rad.png')



	tot_herit = sp.array(cvt_summary_dict['tss_dist']['avg_cis_herit']) + \
		sp.array(cvt_summary_dict['tss_dist']['avg_trans_herit'])
	cis_herit = sp.array(cvt_summary_dict['tss_dist']['avg_cis_herit'])
	pylab.figure(figsize=(10, 6))
	pylab.axes([0.06, 0.08, 0.92, 0.90])
	pylab.fill_between([0, 6], 0, 1, color='#DD3333', alpha=0.8, label='Error')
	pylab.fill_between(sp.arange(7), 0, tot_herit, color='#22CC44', alpha=0.8, label='Heritable variance')
	pylab.fill_between(sp.arange(7), 0, cis_herit, color='#2255AA', \
				alpha=0.8, label='Heritable variance (cis)')
	pylab.ylabel('Average partition of variance')
	pylab.xlabel('Dist. upstream from gene TSS (kb)')
	pylab.xticks(range(7), [200, 100, 50, 25, 10, 5, 1])
	pylab.legend(loc=1, ncol=3, shadow=True)
	pylab.axis([0, 6, 0, 1])
	pylab.savefig(results_prefix + 'avg_herit_td.png')



	pylab.figure()
	pylab.hist(heritabilities, bins=20, alpha=0.7)
	pylab.xlabel('Pseudo-heritability')
	pylab.xlim((-0.025, 1.025))
	pylab.savefig(results_prefix + '_herits_hist.png')
	pylab.clf()

	ks_list = []
	pm_list = []
	for mm in ['EX', 'LM', 'KW']:
		ks_list.append(pval_infl_dict[mm]['kolmogorov_smirnov'])
		pm_list.append(pval_infl_dict[mm]['median_pvals'])

	png_file_name = results_prefix + '_kolmogorov_smirnov_boxplot.png'
	pylab.figure()
	pylab.boxplot(ks_list)
	pylab.axhline(0, color='k', alpha=0.6, ls='-.')
	pylab.xticks(range(1, 4), ['EX', 'LM', 'KW'])
	pylab.ylabel('Kolmogorov-Smirnov statistic D.')
	pylab.savefig(png_file_name)
	pylab.clf()


	png_file_name = results_prefix + '_median_pvals_boxplot.png'
	pylab.figure()
	pylab.boxplot(pm_list)
	pylab.axhline(0, color='k', alpha=0.6, ls='-.')
	pylab.xticks(range(1, 4), ['EX', 'LM', 'KW'])
	pylab.ylabel('Median p-value bias')
	pylab.savefig(png_file_name)
	pylab.clf()


	x_positions = sp.arange(len(distance_bins), dtype='d64')
	width = 0.25
	png_file_name = results_prefix + '_dist_min_pval_hist.png'
	pylab.axes([0.08, 0.2, 0.91, 0.75])
	for mm, color in zip(['EX', 'LM', 'KW'], ['b', 'c', 'g']):
		l = [dist_min_pval_dict[mm][bin] for bin in distance_bins]
		tot_sum = sum(l)
		l = map(lambda x: x / float(tot_sum), l)
		pylab.bar(x_positions, l, width, color=color, alpha=0.7, label=mm)
		x_positions += width


	pylab.ylabel('Frequency')
	pylab.xticks(x_positions - 3 * width / 2.0, (r'$d \leq 5$', r'$5< d \leq 10$', r'$10< d \leq 25$', \
						r'$25< d \leq 50$', r'$50< d \leq 100$', r'$d>100$', \
						'Other chrom.'), rotation='45')
	pylab.xlabel('Distance $d$ (kb) to the smallest p-value from the gene.')
	pylab.xlim((-0.25, len(distance_bins)))
	pylab.legend(loc=2)
	pylab.savefig(png_file_name)
	pylab.clf()


	x_positions = sp.arange(len(radius_bins) + 1, dtype='d64')
	width = 0.25
	png_file_name = results_prefix + 'bonf_sign_bin_hist.png'
	pylab.axes([0.08, 0.22, 0.91, 0.73])
	for mm, color in zip(['EX', 'LM', 'KW'], ['b', 'c', 'g']):
		l = [bonf_sign_bin_dict[mm][bin]['count'] / bonf_sign_bin_dict[mm][bin]['total'] for bin in radius_bins]
		l.append(sign_count[mm] / float(num_genes))
		pylab.bar(x_positions, l, width, color=color, alpha=0.7, label=mm)
		x_positions += width


	pylab.ylabel('Fraction of sign. results')
	pylab.xticks(x_positions - 3 * width / 2.0, ('Within gene', r'$d \leq 1$', r'$d \leq 5$', \
						r'$d \leq 10$', r'$d \leq 25$', r'$d \leq 50$', \
						r'$d \leq 100$', 'Whole genome'), rotation='45')
	pylab.xlabel(r'Among SNPs with distance $d$ (kb) from gene.')
	pylab.xlim((-0.25, len(radius_bins) + 1))
	pylab.legend(loc=2)
	pylab.savefig(png_file_name)
	pylab.clf()


	png_file_name = results_prefix + 'cofactor_count_hist.png'
	x_positions = sp.arange(6, dtype='d64')
	width = 0.25
	for criteria, color in zip(['ebics', 'mbonf', 'min_cof_ppa'], ['b', 'c', 'g']):
		bin_counts = list(sp.bincount(cofactor_count_dict[criteria]['num_cofactor_list']))
		while len(bin_counts) < 6:
			bin_counts.append(0)
		pylab.bar(x_positions, bin_counts, width, color=color, alpha=0.7, label=criteria)
		x_positions += width
	pylab.xlabel('Number of cofactor SNPs')
	pylab.ylabel('Number of genes')
	pylab.xticks(x_positions - 3 * width / 2.0, ('0', '1', '2', '3', '4', '5'))
	pylab.legend(loc=1)
	pylab.xlim((-0.2, 6))
	pylab.savefig(png_file_name)
	pylab.clf()


	png_file_name = results_prefix + 'cis_cofactor_count_hist.png'
	x_positions = sp.arange(6, dtype='d64')
	for criteria, color in zip(['ebics', 'mbonf', 'min_cof_ppa'], ['b', 'c', 'g']):
		bin_counts = list(sp.bincount(cofactor_count_dict[criteria]['num_cis_cofactor_list']))
		while len(bin_counts) < 6:
			bin_counts.append(0)
		pylab.bar(x_positions, bin_counts, width, color=color, alpha=0.7, label=criteria)
		x_positions += width
	pylab.xlabel('Number of cis cofactor SNPs')
	pylab.ylabel('Number of genes')
	pylab.xticks(x_positions - 3 * width / 2.0, ('0', '1', '2', '3', '4', '5'))
	pylab.legend(loc=1)
	pylab.xlim((-0.2, 6))
	pylab.savefig(png_file_name)
	pylab.clf()


	png_file_name = results_prefix + 'cofactor_bin_count_hist.png'
	x_positions = sp.arange(9, dtype='d64')
	width = 0.25
	pylab.axes([0.08, 0.2, 0.91, 0.75])
	for criteria, color in zip(['ebics', 'mbonf', 'min_cof_ppa'], ['b', 'c', 'g']):
		cofactor_count_dict[criteria]['bin_counts'] = \
			cofactor_count_dict[criteria]['bin_counts'] / cofactor_count_dict[criteria]['num_found']
		l = list(cofactor_count_dict[criteria]['bin_counts'])
		l.reverse()
		pylab.bar(x_positions, l, width, color=color, alpha=0.7, label=criteria)
		x_positions += width
	pylab.ylabel('Fraction all genes with cofactors.')
	pylab.xlabel(r'Distance $d$ (kb) to cofactor from gene.')
	pylab.xticks(x_positions - 3 * width / 2.0, ('Within gene', r'$1\geq d$', r'$5\geq d$', r'$10\geq d$', \
						r'$25\geq d$', r'$50\geq d$', r'$100\geq d$', \
						r'$d>100$', 'Other chrom.'), rotation='45')
	pylab.xlim((-0.2, 9))
	pylab.legend(loc=2)
	pylab.savefig(png_file_name)
	pylab.clf()








def run_parallel_rna_seq_gwas():
	if len(sys.argv) > 4:
		run_id = sys.argv[5]
		call_method_id = int(sys.argv[4])
		temperature = int(sys.argv[3])
		phen_file = '%s_%dC.csv' % (phen_file_prefix, temperature)
		file_prefix = env['results_dir'] + 'rna_seq_%s_%dC' % (run_id, temperature)
		run_gwas(file_prefix, phen_file, int(sys.argv[1]), int(sys.argv[2]), temperature,
			data_format='binary', call_method_id=call_method_id, near_const_filter=near_const_filter)
	else:
		call_method_id = int(sys.argv[3])
		temperature = sys.argv[2]
		phen_file = '%s_%sC.csv' % (phen_file_prefix, temperature)
		phed = pd.parse_phenotype_file(phen_file, with_db_ids=False)  #load phenotype file
		phed.filter_near_const_phens(near_const_filter)
		num_traits = phed.num_traits()
		print 'Found %d traits' % num_traits
		chunck_size = int(sys.argv[1])
		for i in range(0, num_traits, chunck_size):
			run_parallel(i, i + chunck_size, temperature, call_method_id)


def _test_parallel_():
	run_parallel(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3])


if __name__ == '__main__':
#	print sys.argv
	if  len(sys.argv) > 3:
		run_parallel_rna_seq_gwas()
#	else:
#		_test_parallel_()
#	sys.exit(0)
#	_test_()
	#load_and_plot_info_files(temperature=10, call_method_id=79, debug_filter=1)
#	plot(min_score=1, temperature=16, mapping_method='KW', call_method_id=79, plot_data=False)
#	plot(min_score=7, temperature=10, mapping_method='KW')
#	plot(min_score=10, temperature=16, mapping_method='KW')
	#plot(min_score=11, temperature=16, mapping_method='KW')
#	plot(min_score=1, temperature=16, mapping_method='LM', call_method_id=79, plot_data=False)
#	plot(min_score=7, temperature=16, mapping_method='LM')
#	plot(min_score=10, temperature=16, mapping_method='LM')
#	plot(min_score=11, temperature=16, mapping_method='LM')
#	plot(min_score=3, temperature=10, mapping_method='EX', plot_data=False)
#	plot(min_score=1, temperature=16, mapping_method='EX', call_method_id=79, plot_data=False)
	print  'Done'

