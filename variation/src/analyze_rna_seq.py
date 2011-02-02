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
from env import *
import sys
import cPickle
import util
import gwaResults

def run_parallel(mapping_method, x_start_i, x_stop_i, cluster='usc'):
	"""
	If no mapping_method, then analysis run is set up.
	"""
	run_id = 'rs'
	job_id = '%s_%d_%d' % (run_id, x_start_i, x_stop_i)
	file_prefix = env['results_dir'] + run_id + '_' + str(x_start_i) + '_' + str(x_stop_i)

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
		shstr += "#PBS -l mem=%s \n" % '2950mb'
		shstr += "#PBS -q cmb\n"
		shstr += "#PBS -N p%s \n" % job_id

	shstr += "(python %sanalyze_rna_seq.py %s %d %d " % (env['script_dir'], mapping_method, x_start_i, x_stop_i)

	shstr += "> " + file_prefix + "_job.out) >& " + file_prefix + "_job.err\n"
	print '\n', shstr, '\n'
	script_file_name = run_id + ".sh"
	f = open(script_file_name, 'w')
	f.write(shstr)
	f.close()

	#Execute qsub script
	os.system("qsub " + script_file_name)




def run_gwas(file_prefix, mapping_method, start_pid, stop_pid, mac_threshold=15, filter_threshold=0.1,
		debug_filter=1.0, use_1001_data=True):
	if mapping_method not in ['emmax', 'kw']:
		raise Exception('Mapping method unknown')
	phen_file = env['phen_dir'] + 'rna_seq.csv'
	phed = pd.parse_phenotype_file(phen_file, with_db_ids=False)  #load phenotype file
	phed.convert_to_averages()
	if use_1001_data:
		sd = dp.load_1001_full_snps(debug_filter=debug_filter)
	else:
		sd = dp.load_250K_snps(debug_filter=debug_filter)
	indices_to_keep = sd.coordinate_w_phenotype_data(phed, 1, coord_phen=False)  #All phenotypes are ordered the same way, so we pick the first one.
	phed.filter_ecotypes(indices_to_keep, pids=range(start_pid, stop_pid))
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
	g_list = _load_genes_list_()
	for pid in range(start_pid, stop_pid):
		gene_tair_id = phed.get_name(pid)
		curr_file_prefix = '%s_%s_mac%d_pid%d_%s' % (file_prefix, mapping_method, mac_threshold, pid, gene_tair_id)
		if phed.is_constant(pid):
			print "Skipping expressions for %s since it's constant." % gene_tair_id
			continue

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

		else:
			raise Exception(mapping_method)
			continue
		res = gr.Result(scores=pvals, macs=macs, mafs=mafs, positions=positions,
				chromosomes=chromosomes)

		#filter, top 10%...
		res.filter_percentile(filter_threshold, reversed=True)
		res.write_to_file(curr_file_prefix + '.pvals')

		#Plot GWAs...
		cgs = g_list[pid - 1]
		print [cg.tairID for cg in cgs]
		f_prefix = curr_file_prefix + '_manhattan'
		res.neg_log_trans()
		res.plot_manhattan(png_file=f_prefix + '.png', percentile=50, cand_genes=cgs, plot_bonferroni=True)

		#Process top regions, where are they...



def _load_results_(mapping_method, file_prefix='', use_1001_data=True, mac_threshold=15):
	pickle_file = '%s_%s_mac%d_res.pickled' % (file_prefix, mapping_method, mac_threshold)
	if os.path.isfile(pickle_file):
		with open(pickle_file) as f:
			d = cPickle.load(f)
	else:
		phen_file = env['phen_dir'] + 'rna_seq.csv'
		phed = pd.parse_phenotype_file(phen_file, with_db_ids=False)  #load phenotype file
		phed.convert_to_averages()
		if use_1001_data:
			sd = dp.load_1001_full_snps(debug_filter=debug_filter)
		else:
			sd = dp.load_250K_snps(debug_filter=debug_filter)
		indices_to_keep = sd.coordinate_w_phenotype_data(phed, 1, coord_phen=False)  #All phenotypes are ordered the same way, so we pick the first one.
		phed.filter_ecotypes(indices_to_keep, pids=range(start_pid, stop_pid))
		print len(indices_to_keep)
		sd.filter_mac_snps(mac_threshold)


		g_list = _load_genes_list_()
		res_dict = {}
		for pid in range(start_pid, stop_pid):
			gene_tair_id = phed.get_name(pid)
			curr_file_prefix = '%s_%s_mac%d_pid%d_%s' % (file_prefix, mapping_method, mac_threshold, pid, gene_tair_id)
			if phed.is_constant(pid):
				print "Skipping expressions for %s since it's constant." % gene_tair_id
				continue
			if phed.is_near_constant(10):
				print "Skipping expressions for %s since it's almost constant." % gene_tair_id
				continue
			print 'Loading file'
			res = gwaResults.Result(curr_file_prefix + '.pvals')
			#Trim results..
			res.neg_log_trans()
			res.filter_attr('scores', 4) #Filter everything below 10^-4
			cgs = g_list[pid - 1]
			avg_g_pos = sp.mean([(cg.startPos + cg.endPos) / 2.0 for cg in cgs])
			chrom = cgs[0].chromosome #Current gene chromosome

			#Prepare for plotting results.. x,y style, where gene is x, and y is p-values
			chrom_pos_score_dict = res.get_chrom_score_pos_dict()
			dist_dict = {}
			for score_threshold in [4, 5, 6, 7]: #negative log10 thresholds.
				res.filter_attr('scores', score_threshold)
				cps_dict = res.get_chrom_score_pos_dict()
				pos_list = cps_dict[chrom]['positions']
				if len(pos_list):
					distances = sp.absolute(sp.array(pos_list) - avg_g_pos)
					d_i = sp.argmin(distances)
					dist_dict[score_threshold] = distances[d_i] #Min distance.
				else:
					dist_dict[score_threshold] = -1 #Different chromosome

			res_dict[(chrom, avg_g_pos)] = {'tair_id':gene_tair_id, 'genes':cgs,
						'chrom_pos_score':chrom_pos_score_dict, 'dist_dict':dist_dict}

		#Now pickle everything...
		d = {'res_dict':res_dict, 'phed':phed}
		with open(pickle_file, 'wb') as f:
			cPickle.dump(d, f)

	return d



#def plot():
#	res_dict = d['res_dict']
#	phed = d['phed']




def run_parallel_rna_seq_gwas():
	if len(sys.argv) > 3:
		run_gwas(env['results_dir'] + 'rna_seq', sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
	else:
		phen_file = env['phen_dir'] + 'rna_seq.csv'
		phed = pd.parse_phenotype_file(phen_file, with_db_ids=False)  #load phenotype file
		num_traits = phed.num_traits()
		print 'Found %d traits' % num_traits
		chunck_size = int(sys.argv[2])
		for i in range(0, num_traits, chunck_size):
			run_parallel(sys.argv[1], i + 1, i + chunck_size + 1)


def _gene_list_to_file_():
	import bisect
	phen_file = env['phen_dir'] + 'rna_seq.csv'
	phed = pd.parse_phenotype_file(phen_file, with_db_ids=False)  #load phenotype file
	tair_gene_versions = phed.get_names()
	tair_ids = [s.split('.')[0] for s in tair_gene_versions]
	gl = gr.get_gene_list(include_intron_exons=False)
	g_tair_ids = [g.tairID for g in gl]
	l = zip(g_tair_ids, gl)
	l.sort()
	(g_tair_ids, gl) = map(list, zip(*l))
	genes = []
	g_chrom_pos = []
	for tid in tair_ids:
		g_i = bisect.bisect_left(g_tair_ids, tid)
		if tid != g_tair_ids[g_i]:
			print tid, g_tair_ids[g_i - 1], g_tair_ids[g_i]
			genes.append([gl[g_i - 1], gl[g_i]])
			pos_l = [gl[g_i - 1].startPos, gl[g_i - 1].endPos, gl[g_i].startPos, gl[g_i].endPos]
		else:
			pos_l = [gl[g_i].startPos, gl[g_i].endPos]
			genes.append([gl[g_i]])
		chr_pos_t = (gl[g_i].chromosome, min(pos_l), max(pos_l))
		if chr_pos_t[2] - chr_pos_t[1] > 50000 or chr_pos_t[2] - chr_pos_t[1] < 0:
			print chr_pos_t
		g_chrom_pos.append(chr_pos_t)
	print len(genes), len(tair_ids), len(g_chrom_pos)
	rna_gene_pickle_file = env['phen_dir'] + 'rna_seq.genes'
	with open(rna_gene_pickle_file, 'wb') as f:
		cPickle.dump(genes, f)


def _load_genes_list_():

	rna_gene_pickle_file = env['phen_dir'] + 'rna_seq.genes'
	with open(rna_gene_pickle_file) as f:
		genes = cPickle.load(f)
	return genes






def _test_():
	run_gwas(env['results_dir'] + 'rna_seq', sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), debug_filter=0.01)

def _test_parallel_():
	run_parallel(sys.argv[1], int(sys.argv[2]), int(sys.argv[2]) + 1)

if __name__ == '__main__':
	_load_results_()
	#_gene_list_to_file_()
	#_test_()
	#run_parallel_rna_seq_gwas()
#	print len(sys.argv)
#	if  len(sys.argv) > 3:
#		run_parallel_rna_seq_gwas()
#	else:
#		_test_parallel_()
#	sys.exit(0)




