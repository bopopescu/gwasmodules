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
phen_file_prefix = env['phen_dir'] + 'rna_seq_061611'


def run_parallel(x_start_i, x_stop_i, temperature, cluster='gmi', run_id='rs'):
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
		shstr += 'export GOTO_NUM_THREADS=1\n'

	elif cluster == 'usc':  #USC cluster.
		shstr = "#!/bin/csh\n"
		shstr += "#PBS -l walltime=%s \n" % '72:00:00'
		shstr += "#PBS -l mem=%s \n" % '2950mb'
		shstr += "#PBS -q cmb\n"
		shstr += "#PBS -N p%s \n" % job_id

	shstr += "python %sanalyze_rna_seq.py %d %d %s %s" % \
			(env['script_dir'], x_start_i, x_stop_i, temperature, run_id)

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





def run_gwas(file_prefix, phen_file, start_i, stop_i, temperature, mac_threshold=15, filter_threshold=0.05,
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
		d = gene_dict[gene_tair_id]
		gene = gwaResults.Gene(chromosome=int(d['chromosome']), startPos=d['start_pos'],
				endPos=d['end_pos'], name=gene_tair_id, description=None, dbRef=gene_tair_id,
				tairID=gene_tair_id)
		print i, pid, gene
		curr_file_prefix = '%s_mac%d_pid%d_%s' % (file_prefix, mac_threshold, pid, gene_tair_id)

		trans_type, shapiro_pval = phed.most_normal_transformation(pid)
		print 'Most normal transformation was: %s' % trans_type
		summary_dict = {'transformation_type':trans_type, 'transformation_shapiro_pval':shapiro_pval}

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
		lm_res = gr.Result(scores=pvals, macs=macs, mafs=mafs, positions=positions, chromosomes=chromosomes)
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
		ex_res = gr.Result(scores=pvals, macs=macs, mafs=mafs, positions=positions, chromosomes=chromosomes)
		summary_dict['EX'] = ex_res.get_gene_analysis(gene)
		summary_dict['pseudo_heritability'] = ex_sw_res['step_info_list'][0]['pseudo_heritability']
		summary_dict['EX']['kolmogorov_smirnov'] = agr.calc_ks_stats(ex_sw_res['first_emmax_res']['ps'])
		summary_dict['EX']['pval_median'] = agr.calc_median(ex_sw_res['first_emmax_res']['ps'])

		#FINISH summarizing the stepwise!!!
		summarize_stepwise(summary_dict, gene, ex_sw_res['step_info_list'], ex_sw_res['opt_dict'])

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





def load_results(mapping_method, temperature, file_prefix='', data_type='quan_seq_data', mac_threshold=15, debug_filter=1.0):
	pickle_file = '%s_%s_%s_mac%d_res.pickled' % (file_prefix, temperature, mapping_method, mac_threshold)
	if os.path.isfile(pickle_file):
		with open(pickle_file) as f:
			d = cPickle.load(f)
	else:
		phen_file = env['phen_dir'] + 'rna_seq_031311_%s.csv' % temperature
		phen_pickle_file = phen_file + 'sd_overlap.pickled'
		if os.path.isfile(phen_pickle_file):
			with file(phen_pickle_file) as f:
				phed = cPickle.load(f)
		else:
			phed = pd.parse_phenotype_file(phen_file, with_db_ids=False)  #load phenotype file
			phed.filter_near_const_phens(15)
			phed.convert_to_averages()
			if data_type == 'use_1001_data':
				sd = dp.load_1001_full_snps(debug_filter=debug_filter)
			elif data_type == 'quan_seq_data':
				sd = dp.load_quan_data(debug_filter=debug_filter)
			else:
				sd = dp.load_250K_snps(debug_filter=debug_filter)
			indices_to_keep = sd.coordinate_w_phenotype_data(phed, 1, coord_phen=False)  #All phenotypes are ordered the same way, so we pick the first one.
			phed.filter_ecotypes(indices_to_keep)
			with file(phen_pickle_file, 'wb') as f:
				cPickle.dump(phed, f)
		#sd.filter_mac_snps(mac_threshold)

		gene_dict = _load_genes_list_('rna_seq_031311_%s' % temperature)
		res_dict = {}
		pids = phed.get_pids()
		for pid in pids:
			print 'pid: %d' % pid
			if not pid in phed.phen_ids: continue
			gene_tair_id = phed.get_name(pid)
			curr_file_prefix = '%s_%s_%s_mac%d_pid%d_%s' % \
				(file_prefix, temperature, mapping_method, mac_threshold, pid, gene_tair_id)
			if phed.is_constant(pid):
				print "Skipping RNA expressions for %s since it's constant." % gene_tair_id
				continue
			if phed.is_near_constant(pid, 15):
				print "Skipping RNA expressions for %s since it's almost constant." % gene_tair_id
				continue
			print 'Loading file'
			if os.path.isfile(curr_file_prefix + '.pvals') or os.path.isfile(curr_file_prefix + '.pvals.pickled'):
				res = gwaResults.Result(curr_file_prefix + '.pvals')
				#Trim results..
				res.neg_log_trans()

				res.filter_attr('scores', 5) #Filter everything below 10^-5
				if res.num_scores() == 0:
					print "Skipping file since nothing is below 10^-5"
					continue
				cgs = gene_dict[pid]
				avg_g_pos = sp.mean([(cg.startPos + cg.endPos) / 2.0 for cg in cgs])
				chrom = cgs[0].chromosome #Current gene chromosome
				print 'Working on gene: %s, chrom=%d, pos=%0.1f' % (gene_tair_id, chrom, avg_g_pos)


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

				res_dict[(chrom, avg_g_pos)] = {'tair_id':gene_tair_id, 'genes':cgs,
							'chrom_pos_score':chrom_pos_score_dict, 'dist_dict':dist_dict}
				print dist_dict
				if 0 < dist_dict[5] < 50000:
					res = gwaResults.Result(curr_file_prefix + '.pvals')
					#Trim results..
					res.neg_log_trans()
					res.plot_manhattan(png_file=curr_file_prefix + '.png', percentile=50,
							cand_genes=cgs, plot_bonferroni=True)

				#print res_dict
			else:
				print 'File not found:', curr_file_prefix + '.pvals'

		#Now pickle everything...
		d = {'res_dict':res_dict, 'phed':phed}
		with open(pickle_file, 'wb') as f:
			cPickle.dump(d, f)
	return d



def plot(file_prefix, results_file_prefix, temperature, mapping_method, min_score=5):
	#Load in chromosome dict..
	chrom_dict = {}
	for x_chrom in [1, 2, 3, 4, 5]:
		for y_chrom in [1, 2, 3, 4, 5]:
			chrom_dict[(x_chrom, y_chrom)] = {'scores':[], 'x_positions':[], 'y_positions':[], 'tair_ids':[]}
	d = _load_results_(mapping_method, temperature, file_prefix=results_file_prefix)
	res_dict = d['res_dict']
	phed = d['phed']
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
	for x_chrom in [1, 2, 3, 4, 5]:
		for y_chrom in [1, 2, 3, 4, 5]:
			file_name = file_prefix + '_chrom%d_chrom%d_%s_%d.txt' % (x_chrom, y_chrom, mapping_method, min_score)
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
	alpha = 0.8
	linewidths = 0
	vmin = min_score
	f = pylab.figure(figsize=(40, 35))
	chromosomes = [1, 2, 3, 4, 5]
	plot_file_name = file_prefix + '_%s_%d.png' % (mapping_method, min_score)
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








def load_and_plot_info_files(file_prefix='', call_method=75, temperature=10, mac_threshold=15, debug_filter=1.0, near_const_filter=20):
	phen_file = env['phen_dir'] + 'rna_seq_031311_%s.csv' % temperature
	phen_pickle_file = phen_file + 'sd_overlap.pickled'
	if os.path.isfile(phen_pickle_file):
		with file(phen_pickle_file) as f:
			phed = cPickle.load(f)
	else:
		phed = pd.parse_phenotype_file(phen_file, with_db_ids=False)  #load phenotype file
		phed.filter_near_const_phens(near_const_filter)
		phed.convert_to_averages()
		if data_type == 'use_1001_data':
			sd = dp.load_1001_full_snps(debug_filter=debug_filter)
		elif data_type == 'quan_seq_data':
			sd = dp.load_quan_data(debug_filter=debug_filter)
		else:
			sd = dp.load_250K_snps(debug_filter=debug_filter)
		indices_to_keep = sd.coordinate_w_phenotype_data(phed, 1, coord_phen=False)  #All phenotypes are ordered the same way, so we pick the first one.
		phed.filter_ecotypes(indices_to_keep)
		with file(phen_pickle_file, 'wb') as f:
			cPickle.dump(phed, f)

	gene_dict = _load_genes_list_('rna_seq_031311_%s' % temperature)
	res_dict = {}
	pids = phed.get_pids()
	heritabilities = []
	transformations = []
	shapiro_wilk_pvals = []
	tair_ids = []
	for i, pid in enumerate(pids):
		if not pid in phed.phen_ids: continue

		gene_tair_id = phed.get_name(pid)
		tair_ids.append(gene_tair_id)
		curr_file_prefix = '%s_%s_%s_mac%d_pid%d_%s' % \
			(file_prefix, temperature, mapping_method, mac_threshold, pid, gene_tair_id)
#		if phed.is_constant(pid):
#			print "Skipping RNA expressions for %s since it's constant." % gene_tair_id
#			continue
		if phed.is_near_constant(pid, 15):
#			print "Skipping RNA expressions for %s since it's almost constant." % gene_tair_id
			continue
		if (i + 1) % (len(pids) / 100) == 0:
			sys.stdout.write('.')
			sys.stdout.flush()
		info_filename = curr_file_prefix + '_info.txt'
		if os.path.isfile(curr_file_prefix + '_info.txt'):
			with open(info_filename) as f:
				line = f.next()
				l = map(str.strip, line.split())
				heritabilities.append(float(l[1]))
				line = f.next()
				l = map(str.strip, line.split())
				transformations.append(l[1])
				line = f.next()
				l = map(str.strip, line.split())
				shapiro_wilk_pvals.append(float(l[1]))
	print ''
	print len(heritabilities), sp.mean(heritabilities)
	pylab.hist(heritabilities, alpha=0.6, bins=20)
	ymin, ymax = pylab.ylim()
	y_range = ymax - 0
	pylab.axis([-.025, 1.025, -.025 * y_range, y_range])
	fig_file = env['tmp_dir'] + 'rna_seq_031311_%s_herits.png' % temperature
	pylab.xlabel('Expression heritability')
	pylab.ylabel('Frequency')
	pylab.savefig(fig_file)

	info_file_name = env['tmp_dir'] + 'rna_seq_info%s.csv' % temperature
	with open(info_file_name, 'w') as f:
		f.write('tair_id, heritability, shapiro_wilks_pval, transformation')
		for tair_id, h, swp, trans in it.izip(tair_ids, heritabilities, shapiro_wilk_pvals, transformations):
			f.write('%s,%f,%f,%s\n' % (tair_id, h, swp, trans))





def run_parallel_rna_seq_gwas():
	if len(sys.argv) > 4:
		run_id = sys.argv[4]
		temperature = int(sys.argv[3])
		phen_file = '%s_%dC.csv' % (phen_file_prefix, temperature)
		file_prefix = env['results_dir'] + 'rna_seq_%s_%dC' % (run_id, temperature)
		run_gwas(file_prefix, phen_file, int(sys.argv[1]), int(sys.argv[2]), temperature,
			data_format='binary', call_method_id=79, near_const_filter=near_const_filter)
	else:
		run_id = sys.argv[3]
		temperature = sys.argv[2]
		phen_file = '%s_%sC.csv' % (phen_file_prefix, temperature)
		phed = pd.parse_phenotype_file(phen_file, with_db_ids=False)  #load phenotype file
		phed.filter_near_const_phens(near_const_filter)
		num_traits = phed.num_traits()
		print 'Found %d traits' % num_traits
		chunck_size = int(sys.argv[1])
		for i in range(0, num_traits, chunck_size):
			run_parallel(i, i + chunck_size, temperature, run_id=run_idr)


def _test_parallel_():
	run_parallel(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3])


if __name__ == '__main__':
	#_load_results_('lm', '16C', file_prefix='/storage/rna_seq_results_032011/rna_seq')
	#load_and_plot_info_files('lm', '16C', file_prefix='/storage/rna_seq_results_032011/rna_seq')
	#plot('/tmp/rna_seq_10C', '/storage/rna_seq_results_032011/rna_seq', '10C', 'emmax', 5)
	#plot('/tmp/rna_seq_10C', '/storage/rna_seq_results_032011/rna_seq', '10C', 'emmax', 7)
	#plot('/tmp/rna_seq_10C', '/storage/rna_seq_results_032011/rna_seq', '10C', 'emmax', 10)
	#plot('/tmp/rna_seq_10C', '/storage/rna_seq_results_032011/rna_seq', '10C', 'emmax', 11)
#	plot('/tmp/rna_seq_16C', '/storage/rna_seq_results_032011/rna_seq', '16C', 'lm', 7)
#	plot('/tmp/rna_seq_16C', '/storage/rna_seq_results_032011/rna_seq', '16C', 'emmax', 7)
#	plot('/tmp/rna_seq_16C', '/storage/rna_seq_results_032011/rna_seq', '16C', 'lm', 10)
#	plot('/tmp/rna_seq_16C', '/storage/rna_seq_results_032011/rna_seq', '16C', 'emmax', 10)
#	plot('/tmp/rna_seq_16C', '/storage/rna_seq_results_032011/rna_seq', '16C', 'lm', 12)
#	plot('/tmp/rna_seq_16C', '/storage/rna_seq_results_032011/rna_seq', '16C', 'emmax', 12)
#	plot('/tmp/rna_seq_10C', '/storage/rna_seq_results_032011/rna_seq', '10C', 'lm', 12)
	#_load_genes_list_()
	#_test_()
	#print sys.argv
	run_parallel_rna_seq_gwas()
#	if  len(sys.argv) > 3:
#		run_parallel_rna_seq_gwas()
#	else:
#		_test_parallel_()
#	sys.exit(0)
	#_test_()



