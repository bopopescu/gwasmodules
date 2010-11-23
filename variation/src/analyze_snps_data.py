"""
Contains basic functions to analyze SNP data, e.g. LD levels, polymorphism rates, Fst, recombination, etc.
"""
import sys
from env import *
import dataParsers as dp
import scipy as sp

def calc_r2_levels(file_prefix, x_start_i, x_stop_i, mac_filter=10, emma_r2_threshold=0.2, min_emma_dist=100000,
		save_threshold=0.2, debug_filter=1):
	"""
	Returns statistics on LD levels, and plot them.
	
	Calculates emma_
	"""
	import scipy as sp
	import scipy.stats as st
	from scipy import linalg
	import linear_models as lm

	sd = dp.parse_numerical_snp_data(env['data_dir'] + '250K_t72.csv.binary',
					filter=debug_filter)
	sd.filter_mac_snps(mac_filter)
	K = lm.load_kinship_from_file(env['data_dir'] + 'kinship_matrix_cm72.pickled',
				sd.accessions)
	cps_list = sd.getChrPosSNPList()
	x_cps = cps_list[x_start_i:x_stop_i]
	y_cps = cps_list
	result_list = []
	i = 0
	j = 0
	q = 1  # Single SNP is being tested
	p = 2
	n = len(sd.accessions)
	n_p = n - p
	for (x_c, x_p, x_snp) in x_cps:
		#get the EMMA REML
		res = lm.get_emma_reml_estimates(x_snp, K)
		h_sqrt_inv = res['H_sqrt_inv']
		h0_X = res['X_t']
		Y = res['Y_t']	#The transformed outputs.
		(h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y)
		print i
		i += 1
		#j = 0
		for (y_c, y_p, y_snp) in y_cps:
			if (x_c, x_p) < (y_c, y_p):
				(r, pearson_pval) = st.pearsonr(x_snp, y_snp)
				r2 = r * r
				if r2 > save_threshold:
					emma_pval = 2
					if (x_c != y_c or abs(x_p - y_p) > min_emma_dist) and r2 > emma_r2_threshold:
						#Do EMMAX
						yt = y_snp * h_sqrt_inv
						(b, rss, p, s) = linalg.lstsq(sp.hstack([h0_X, sp.matrix(yt).T]), Y)
						if rss:
							f_stats = (h0_rss / rss[0] - 1) * n_p / float(q)
							emma_pval = st.f.sf(f_stats, q, n_p)[0]
					if emma_pval < 0.01:
						print [x_c, x_p, y_c, y_p, r2, pearson_pval, emma_pval]
					result_list.append([x_c, x_p, y_c, y_p, r2, pearson_pval, emma_pval])
	file_name = file_prefix + '_x_' + str(x_start_i) + '_' + str(x_stop_i) + ".csv"
	f = open(file_name, 'w')
	for r in result_list:
		st = ','.join(map(str, r))
		f.write(st + '\n')
	f.close()
	return result_list




def run_parallel(x_start_i, x_stop_i):
	"""
	If no mapping_method, then analysis run is set up.
	"""

	#Cluster specific parameters
	run_id = 'r2_250k'


	shstr = "#!/bin/csh\n"
	shstr += "#PBS -l walltime=%s \n" % '72:00:00'
	shstr += "#PBS -l mem=%s \n" % '1200mb'
	shstr += "#PBS -q cmb\n"

	job_id = '%s_%d_%d' % (run_id, x_start_i, x_stop_i)
	shstr += "#PBS -N p%s \n" % job_id
	shstr += "(python %sanalyze_snps_data.py %d %d" % (env['script_dir'], x_start_i, x_stop_i)

	file_prefix = env['results_dir'] + run_id + '_' + str(x_start_i) + '_' + str(x_stop_i)
	shstr += "> " + file_prefix + "_job.out) >& " + file_prefix + "_job.err\n"
	print '\n', shstr, '\n'
	script_file_name = run_id + ".sh"
	f = open(script_file_name, 'w')
	f.write(shstr)
	f.close()

	#Execute qsub script
	os.system("qsub " + script_file_name)


def run_r2_calc():
	if len(sys.argv) > 2:
		calc_r2_levels(env['results_dir'] + '250K_r2_min02', int(sys.argv[1]), int(sys.argv[2]))
	else:
		sd = dp.parse_numerical_snp_data(env['data_dir'] + '250K_t72.csv.binary')
		num_snps = len(sd.getSnps())
		chunck_size = int(sys.argv[1])
		for i in range(0, num_snps, chunck_size):
			run_parallel(i, i + chunck_size)


def load_and_plot_r2_results():
	import pylab
	import os
	import cPickle
	import math
	min_float = 5e-324
	sd = dp.parse_numerical_snp_data(env['data_dir'] + '250K_t72.csv.binary')
	num_snps = len(sd.getSnps())
	file_prefix = '/Users/bjarni.vilhjalmsson/Projects/250K_r2/results/250K_r2_min02'
	headers = ['x_chr', 'x_pos', 'y_chr', 'y_pos', 'r2', 'pval', 'emmax_pval']
	if os.path.isfile(file_prefix + '.pickled'):
		f = open(file_prefix + '.pickled', 'rb')
		res_dict = cPickle.load(f)
		f.close()
	else:
		chunck_size = int(sys.argv[1])
		res_dict = {}
		for h in headers:
			res_dict[h] = []
		delim = ','
		for i in range(0, num_snps, chunck_size):
			file_name = file_prefix + '_x_' + str(i) + '_' + str(i + chunck_size) + ".csv"
			print i
			f = open(file_name)
			for line in f:
				l = map(str.strip, line.split(delim))
				for j, st in enumerate(l):
					h = headers[j]
					if h in ['x_chr', 'x_pos', 'y_chr', 'y_pos']:
						res_dict[h].append(int(st))
					elif h in ['pval', 'emmax_pval']:
						v = float(st)
						res_dict[h].append(v if v != 0.0 else min_float)
					elif h == 'r2':
						res_dict[h].append(float(st))
					else:
						raise Exception()
		f = open(file_prefix + '.pickled', 'wb')
		cPickle.dump(res_dict, f, 2)
		f.close()

	max_pval = -math.log10(min_float)
	max_emmax_pval = -math.log10(min(res_dict['emmax_pval']))
	#Filter data..
	num_res = len(res_dict['x_chr'])
	chromosomes = [1, 2, 3, 4, 5]
	chr_res_dict = {}
	for chr2 in chromosomes:
		for chr1 in chromosomes[:chr2]:
			d = {}
			for h in headers:
				d[h] = []
			chr_res_dict[(chr1, chr2)] = d
	num_retained = 0
	chr_pos_set = set()
	for i in range(num_res):
		x_chr = res_dict['x_chr'][i]
		y_chr = res_dict['y_chr'][i]
		x_pos = res_dict['x_pos'][i]
		y_pos = res_dict['y_pos'][i]
		x_chr_pos = (x_chr, x_pos)
		y_chr_pos = (y_chr, y_pos)
		if x_chr == y_chr and  y_pos - x_pos < 25000:
			if res_dict['r2'][i] > 0.6:
				for h in headers:
					chr_res_dict[(x_chr, y_chr)][h].append(res_dict[h][i])
				num_retained += 1
				chr_pos_set.add((x_chr, x_pos))
				chr_pos_set.add((y_chr, y_pos))
		if x_chr == y_chr and  y_pos - x_pos < 50000:
			if res_dict['r2'][i] > 0.5:
				for h in headers:
					chr_res_dict[(x_chr, y_chr)][h].append(res_dict[h][i])
				num_retained += 1
				chr_pos_set.add((x_chr, x_pos))
				chr_pos_set.add((y_chr, y_pos))
		elif x_chr == y_chr and  y_pos - x_pos < 100000:
			if res_dict['r2'][i] > 0.4:
				for h in headers:
					chr_res_dict[(x_chr, y_chr)][h].append(res_dict[h][i])
				num_retained += 1
				chr_pos_set.add((x_chr, x_pos))
				chr_pos_set.add((y_chr, y_pos))

		elif res_dict['r2'][i] > 0.2:
			for h in headers:
					chr_res_dict[(x_chr, y_chr)][h].append(res_dict[h][i])
			num_retained += 1
			chr_pos_set.add((x_chr, x_pos))
			chr_pos_set.add((y_chr, y_pos))

	print num_retained
	print len(chr_pos_set)
	#Now plot data!!
	alpha = 0.8
	linewidths = 0
	vmin = 0.0
	for chr2 in chromosomes:
		for chr1 in chromosomes[:chr2]:
			x_min = min(chr_res_dict[(chr1, chr2)]['x_pos'])
			x_max = max(chr_res_dict[(chr1, chr2)]['x_pos'])
			y_min = min(chr_res_dict[(chr1, chr2)]['y_pos'])
			y_max = max(chr_res_dict[(chr1, chr2)]['y_pos'])
			x_range = x_max - x_min
			y_range = y_max - y_min
			x_lab = 'Chromosome %d' % chr1
			y_lab = 'Chromosome %d' % chr2
			left = 0.05
			bottom = 0.04
			width = 0.94
			height = 0.94
			r2_plot_file_name = file_prefix + 'c_' + str(chr1) + 'x' + str(chr2) + '_r2s.png'
			pval_file_name = file_prefix + 'c_' + str(chr1) + 'x' + str(chr2) + '_pvals.png'
			emma_pval_file_name = file_prefix + 'c_' + str(chr1) + 'x' + str(chr2) + '_emmax_pvals.png'
#			print r2_plot_file_name, pval_file_name , emma_pval_file_name
			pylab.figure(figsize=(18, 16))
			pylab.axes([left, bottom, width, height])

			l_zxy = zip(chr_res_dict[(chr1, chr2)]['r2'], chr_res_dict[(chr1, chr2)]['x_pos'],
				chr_res_dict[(chr1, chr2)]['y_pos'])
			l_zxy.sort()
			l = map(list, zip(*l_zxy))
			zs = l[0]
			xs = l[1]
			ys = l[2]
#			print len(chr_res_dict[(chr1, chr2)]['x_pos']), len(chr_res_dict[(chr1, chr2)]['y_pos']), \
#				len(chr_res_dict[(chr1, chr2)]['r2'])
			pylab.scatter(xs, ys, c=zs, alpha=alpha, linewidths=linewidths, vmin=vmin, vmax=1.0)
			pylab.axis([x_min - 0.025 * x_range, x_max + 0.025 * x_range, y_min - 0.025 * y_range,
				y_max + 0.025 * y_range])
			pylab.colorbar()
			pylab.xlabel(x_lab)
			pylab.ylabel(y_lab)
			pylab.savefig(r2_plot_file_name, format='png')
			pylab.clf()
			pylab.figure(figsize=(18, 16))
			pylab.axes([left, bottom, width, height])
			log_pvals = map(lambda x:-math.log10(x), chr_res_dict[(chr1, chr2)]['pval'])
			l_zxy = zip(log_pvals, chr_res_dict[(chr1, chr2)]['x_pos'], chr_res_dict[(chr1, chr2)]['y_pos'])
			l_zxy.sort()
			l = map(list, zip(*l_zxy))
			zs = l[0]
			xs = l[1]
			ys = l[2]
			pylab.scatter(xs, ys, c=zs, alpha=alpha, linewidths=linewidths, vmin=vmin, vmax=max_pval)
			pylab.axis([x_min - 0.025 * x_range, x_max + 0.025 * x_range, y_min - 0.025 * y_range,
				y_max + 0.025 * y_range])
			pylab.colorbar()
			pylab.xlabel(x_lab)
			pylab.ylabel(y_lab)
			pylab.savefig(pval_file_name, format='png')
			pylab.clf()
			pylab.figure(figsize=(18, 16))
			pylab.axes([left, bottom, width, height])
			log_pvals = map(lambda x:-math.log10(x), chr_res_dict[(chr1, chr2)]['emmax_pval'])
			l_zxy = zip(log_pvals, chr_res_dict[(chr1, chr2)]['x_pos'], chr_res_dict[(chr1, chr2)]['y_pos'])
			l_zxy.sort()
			i = 0
			while l_zxy[i] < 0:
				i += 1
			l_zxy = l_zxy[i:]
			l = map(list, zip(*l_zxy))
			zs = l[0]
			xs = l[1]
			ys = l[2]
			pylab.scatter(xs, ys, c=zs, alpha=alpha, linewidths=linewidths, vmin=vmin, vmax=max_emmax_pval)
			pylab.axis([x_min - 0.025 * x_range, x_max + 0.025 * x_range, y_min - 0.025 * y_range,
				y_max + 0.025 * y_range])
			pylab.colorbar()
			pylab.xlabel(x_lab)
			pylab.ylabel(y_lab)
			pylab.savefig(emma_pval_file_name, format='png')

#plt.subplot(121)
#plt.scatter(xyc[:13], xyc[:13], c=xyc[:13], s=35, vmin=0, vmax=20)
#plt.colorbar()
#plt.xlim(0, 20)
#plt.ylim(0, 20)
#
#plt.subplot(122)
#plt.scatter(xyc[8:20], xyc[8:20], c=xyc[8:20], s=35, vmin=0, vmax=20)   
#plt.colorbar()
#plt.xlim(0, 20)
#plt.ylim(0, 20)









if __name__ == "__main__":
	load_and_plot_r2_results()
