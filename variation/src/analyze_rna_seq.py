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

def run_parallel(mapping_method, x_start_i, x_stop_i, temperature, cluster='usc'):
	"""
	If no mapping_method, then analysis run is set up.
	"""
	run_id = 'rs'
	job_id = '%s_%d_%d_%s' % (run_id, x_start_i, x_stop_i, temperature)
	file_prefix = env['results_dir'] + job_id

	#Cluster specific parameters	
	if cluster == 'gmi': #GMI cluster.
		shstr = '#!/bin/sh\n'
		shstr += '#$ -N %s\n' % job_id
		shstr += "#$ -q q.norm@blade*\n"
		shstr += '#$ -o %s.log\n' % job_id
		#shstr += '#$ -cwd /home/GMI/$HOME\n'
		#shstr += '#$ -M bjarni.vilhjalmsson@gmi.oeaw.ac.at\n\n'

	elif cluster == 'usc':  #USC cluster.
		shstr = "#!/bin/csh\n"
		shstr += "#PBS -l walltime=%s \n" % '72:00:00'
		shstr += "#PBS -l mem=%s \n" % '4950mb'
		shstr += "#PBS -q cmb\n"
		shstr += "#PBS -N p%s \n" % job_id

	shstr += "(python %sanalyze_rna_seq.py %s %d %d %s" % \
			(env['script_dir'], mapping_method, x_start_i, x_stop_i, temperature)

	shstr += "> " + file_prefix + "_job.out) >& " + file_prefix + "_job.err\n"
	print '\n', shstr, '\n'
	script_file_name = run_id + ".sh"
	f = open(script_file_name, 'w')
	f.write(shstr)
	f.close()

	#Execute qsub script
	os.system("qsub " + script_file_name)




def run_gwas(file_prefix, mapping_method, start_i, stop_i, temperature, mac_threshold=15, filter_threshold=0.05,
		debug_filter=1.0, use_1001_data=True):
	if mapping_method not in ['emmax', 'kw']:
		raise Exception('Mapping method unknown')
	phen_file = env['phen_dir'] + 'rna_seq_031311_%s.csv' % temperature
	phed = pd.parse_phenotype_file(phen_file, with_db_ids=False)  #load phenotype file
	phed.filter_near_const_phens(15)
	phed.convert_to_averages()
	num_traits = phed.num_traits()
	pids = phed.phen_ids[start_i :stop_i]
	if use_1001_data:
		sd = dp.load_1001_full_snps(debug_filter=debug_filter)
	else:
		sd = dp.load_250K_snps(debug_filter=debug_filter)
	indices_to_keep = sd.coordinate_w_phenotype_data(phed, 1, coord_phen=False)  #All phenotypes are ordered the same way, so we pick the first one.
	phed.filter_ecotypes(indices_to_keep, pids=pids)
	print len(indices_to_keep)
	if mapping_method == 'emmax':
		k_file = env['data_1001_dir'] + 'kinship_matrix.pickled'
		K = lm.load_kinship_from_file(k_file, sd.accessions)
	sd.filter_mac_snps(mac_threshold)
	snps = sd.getSnps()
	positions = sd.getPositions()
	chromosomes = sd.get_chr_list()
	r = sd.get_mafs()
	macs = r['mafs']
	mafs = r['marfs']
	print 'In total there are %d SNPs to be mapped.' % len(snps)
	gene_dict = _load_genes_list_('rna_seq_031311_%s' % temperature)
	for i, pid in enumerate(pids):
		if not pid in phed.phen_ids: continue
		cgs = gene_dict[pid]
		print i
		gene_tair_id = phed.get_name(pid)
		curr_file_prefix = '%s_%s_mac%d_pid%d_%s' % (file_prefix, mapping_method, mac_threshold, pid, gene_tair_id)

		if mapping_method == 'kw':
			phen_vals = phed.get_values(pid)
			kw_res = util.kruskal_wallis(snps, phen_vals)
			pvals = kw_res['ps'].tolist()

		elif mapping_method == 'emmax':
			#Identify the right transformation
			trans_type, shapiro_pval = phed.most_normal_transformation(pid)
			phen_vals = phed.get_values(pid)
			#Get pseudo-heritabilities
			res = lm.get_emma_reml_estimates(phen_vals, K)

			f_prefix = curr_file_prefix + '_hist'
			phed.plot_histogram(pid, title='Gene expressions for %s' % gene_tair_id,
					png_file=f_prefix + '.png', p_her=res['pseudo_heritability'],
					x_label='RNA seq expression levels (%s transformed)' % trans_type)
			res = lm.emmax(snps, phen_vals, K)
			pvals = res['ps'].tolist()
			p_her = res['pseudo_heritability']


		else:
			raise Exception(mapping_method)
			continue

		res = gr.Result(scores=pvals, macs=macs, mafs=mafs, positions=positions, chromosomes=chromosomes)
		num_scores = res.get_num_scores()

		#Record distance to the most significant SNP from the gene.
		min_i = res.arg_min_attr()
		min_chr = chromosomes[min_i]
		min_pos = positions[min_i]
		min_dist = -1
		min_pos = 100000000
		max_pos = 0
		cur_chrom = cgs[0].chromosome
		for g in cgs:
			if min_dist != 0 and g.chromosome == min_chr:
				if g.startPos <= min_pos <= g.endPos:
					min_dist = 0
				else:
					dist = min(abs(g.startPos - min_pos), abs(g.endPos - min_pos))
					if min_dist == -1 or min_dist > dist:
						min_dist = dist
			min_pos = min(g.startPos, min_pos)
			max_pos = max(g.endPos, max_pos)

		#Record most significant p-value within 0kb, 5kb, 10kb, 25kb, 50kb, and 100kb window, and SNP count.
		window_sizes = [0, 5000, 10000, 25000, 50000, 100000]
		window_dict = {}
		for window_size in window_sizes:
			reg_res = res.get_region_result(cur_chrom, min_pos, max_pos, buffer=window_size)
			num_snps = reg_res.num_scores()
			if num_scores:
				min_pval = reg_res.min_score()
			else:
				min_pval = -1
			window_dict[window_size] = {'min_pval':min_pval, 'num_snps':num_snps}


		#Write info to file..
		with open(curr_file_prefix + '_info.txt', 'w') as f:
			if mapping_method == 'emmax':
				f.write('pseudo_heritability: %f \n' % p_her)
				f.write('transformation_type: %s \n' % trans_type)
				f.write('transformation_shapiro_pval: %f \n' % shapiro_pval)
			f.write('dist_to_cand_gene: %f \n' % min_dist)
			for window_size in window_sizes:
				min_pval = window_dict[window_size]['min_pval']
				num_snps = window_dict[window_size]['num_snps']
				f.write('window_size: %d, min_pval: %f \n' % (window_size, min_pval))
				f.write('window_size: %d, num_snps: %d \n' % (window_size, num_snps))



		#filter, top 10%...
		res.filter_percentile(filter_threshold, reversed=True)
		res.write_to_file(curr_file_prefix + '.pvals')



		#Plot GWAs...
		if res.min_score() < 10e-9:
			#print [cg.tairID for cg in cgs]
			f_prefix = curr_file_prefix + '_manhattan'
			res.neg_log_trans()
			res.plot_manhattan(png_file=f_prefix + '.png', percentile=50, cand_genes=cgs, plot_bonferroni=True,
					b_threshold= -math.log10(1.0 / (num_traits * num_scores * 20.0)))





def _load_results_(mapping_method, temperature, file_prefix='', use_1001_data=True, mac_threshold=15, debug_filter=1.0):
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
			if use_1001_data:
				sd = dp.load_1001_full_snps(debug_filter=debug_filter)
			else:
				sd = dp.load_250K_snps(debug_filter=debug_filter)
			indices_to_keep = sd.coordinate_w_phenotype_data(phed, 1, coord_phen=False)  #All phenotypes are ordered the same way, so we pick the first one.
			phed.filter_ecotypes(indices_to_keep)
			with file(phen_pickle_file, 'wb') as f:
				cPickle.dump(phed, f)
		#sd.filter_mac_snps(mac_threshold)

		gene_dict = _load_genes_list_(temperature=temperature)
		res_dict = {}
		pids = phed.get_pids()
		for pid in pids:
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
			if os.path.isfile(curr_file_prefix + '.pvals'):
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
			file_name = file_prefix + '_chrom%d_chrom%d.txt' % (x_chrom, y_chrom)
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
	f = pylab.figure(figsize=(50, 45))
	chromosomes = [1, 2, 3, 4, 5]
	plot_file_name = file_prefix + '_emmax_%d.png' % (min_score)
	label = '$-log_{10}$(p-value)'
	vmax = max(scores)

	for yi, chr2 in enumerate(chromosomes):
		for xi, chr1 in enumerate(chromosomes):

			l = chrom_dict[(chr1, chr2)]['scores']
			if len(l) == 0:
				continue
			ax = f.add_axes([0.97 * (rel_cum_chrom_sizes[xi] + 0.01), rel_cum_chrom_sizes[yi] - 0.02,
					0.97 * (rel_chrom_sizes[xi]), rel_chrom_sizes[yi] ])
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

	cax = f.add_axes([0.975, 0.7, 0.01, 0.2])
	cb = pylab.colorbar(scatter_plot, cax=cax)
	cb.set_label(label, fontsize='xx-large')
	#cb.set_tick_params(fontsize='x-large')
	f.text(0.005, 0.47, 'Expressed gene position', size='xx-large', rotation='vertical')
	f.text(0.47, 0.99, 'Associated SNP position', size='xx-large')
	print 'Saving figure:', plot_file_name
	f.savefig(plot_file_name, format='png')




def load_info_files(mapping_method, temperature, file_prefix='', use_1001_data=True, debug_filter=1.0):
	pickle_file = '%s_%s_%s_mac%d_res.pickled' % (file_prefix, temperature, mapping_method)
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
			if use_1001_data:
				sd = dp.load_1001_full_snps(debug_filter=debug_filter)
			else:
				sd = dp.load_250K_snps(debug_filter=debug_filter)
			indices_to_keep = sd.coordinate_w_phenotype_data(phed, 1, coord_phen=False)  #All phenotypes are ordered the same way, so we pick the first one.
			phed.filter_ecotypes(indices_to_keep)
			with file(phen_pickle_file, 'wb') as f:
				cPickle.dump(phed, f)

		gene_dict = _load_genes_list_(temperature=temperature)
		res_dict = {}
		pids = phed.get_pids()
		heritabilities = []
		for pid in pids:
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
			info_filename = curr_file_prefix + '_info.txt'
			if os.path.isfile(curr_file_prefix + '_info.txt'):
				with open(info_filename) as f:
					f.next()




def run_parallel_rna_seq_gwas():
	if len(sys.argv) > 4:
		run_gwas(env['results_dir'] + 'rna_seq_10C', sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
	else:
		temperature = sys.argv[3]
		phen_file = env['phen_dir'] + 'rna_seq_031311_%s.csv' % temperature
		phed = pd.parse_phenotype_file(phen_file, with_db_ids=False)  #load phenotype file
		phed.filter_near_const_phens(15)
		num_traits = phed.num_traits()
		print 'Found %d traits' % num_traits
		chunck_size = int(sys.argv[2])
		for i in range(0, num_traits, chunck_size):
			run_parallel(sys.argv[1], i, i + chunck_size, temperature)


def _gene_list_to_file_(file_prefix='rna_seq_031311_10C'):
	import bisect
	phen_file = env['phen_dir'] + file_prefix + '.csv'
	phed = pd.parse_phenotype_file(phen_file, with_db_ids=False)  #load phenotype file
	phed.filter_near_const_phens(15)
	tair_gene_versions = phed.get_names()
	tair_ids = [s.split('.')[0] for s in tair_gene_versions]
	pids = phed.phen_ids
	gl = gr.get_gene_list(include_intron_exons=False)
	g_tair_ids = [g.tairID for g in gl]
	l = zip(g_tair_ids, gl)
	l.sort()
	(g_tair_ids, gl) = map(list, zip(*l))
	gene_dict = {}
	g_chrom_pos = []
	for pid, tid in zip(pids, tair_ids):
		g_i = bisect.bisect_left(g_tair_ids, tid)
		if tid != g_tair_ids[g_i]:
			print tid, g_tair_ids[g_i - 1], g_tair_ids[g_i]
			gene_dict[pid] = [gl[g_i - 1], gl[g_i]]
			pos_l = [gl[g_i - 1].startPos, gl[g_i - 1].endPos, gl[g_i].startPos, gl[g_i].endPos]
		else:
			pos_l = [gl[g_i].startPos, gl[g_i].endPos]
			gene_dict[pid] = [gl[g_i]]
		chr_pos_t = (gl[g_i].chromosome, min(pos_l), max(pos_l))
		if chr_pos_t[2] - chr_pos_t[1] > 50000 or chr_pos_t[2] - chr_pos_t[1] < 0:
			print chr_pos_t
		g_chrom_pos.append(chr_pos_t)
	print len(gene_dict), len(tair_ids), len(g_chrom_pos)
	rna_gene_pickle_file = env['phen_dir'] + file_prefix + '.genes'
	with open(rna_gene_pickle_file, 'wb') as f:
		cPickle.dump(gene_dict, f)


def _load_genes_list_(file_prefix='rna_seq_031311', temperature='10C'):
	file_prefix = file_prefix + '_' + temperature
	rna_gene_pickle_file = env['phen_dir'] + file_prefix + '.genes'
	if not os.path.isfile(rna_gene_pickle_file):
		_gene_list_to_file_(file_prefix)
	with open(rna_gene_pickle_file) as f:
		gene_dict = cPickle.load(f)
	return gene_dict






def _test_():
	run_gwas(env['results_dir'] + 'rna_seq', sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), debug_filter=0.01)

def _test_parallel_():
	run_parallel(sys.argv[1], int(sys.argv[2]), int(sys.argv[2]) + 1)

if __name__ == '__main__':
	#_load_results_('emmax', '16C', file_prefix='/storage/rna_seq_results/rna_seq')
	plot('/tmp/rna_seq_16C', '/storage/rna_seq_results/rna_seq', '16C', 'emmax', 11)
	#_load_genes_list_()
	#_test_()
	#print sys.argv
	#run_parallel_rna_seq_gwas()
#	if  len(sys.argv) > 3:
#		run_parallel_rna_seq_gwas()
#	else:
#		_test_parallel_()
#	sys.exit(0)




