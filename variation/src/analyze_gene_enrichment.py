import gwaResults as gr
import phenotypeData as pd
import random
#import pylab as plt
#import matplotlib
import scipy as sp
import random, math
import rpy
from scipy import stats as st
import env
import sys
import dataParsers

#def num_enrich_cand_genes(cand_gene_indices, genes, result, pvalue_threshold):
#	num_enriched_cand_genes = 0
#	last_start = 0
#	for cgi in cand_gene_indices:
#		gene = genes[cgi]
#		i = last_start
#		while i < len(result.positions) and (result.chromosomes[i] < gene.chromosome or result.positions[i] < gene.startPos):
#			i += 1
#		last_start = i
#		mean_gene_pval = 0.0
#		num_pvals = 0
#		while i < len(result.positions) and result.chromosomes[i] == gene.chromosome and result.positions[i] <= gene.endPos:
#			mean_gene_pval += result.scores[i]
#			num_pvals += 1
#			i += 1
#		if mean_gene_pval != 0.0:
#			mean_gene_pval = mean_gene_pval / num_pvals
#		if mean_gene_pval > pvalue_threshold:
#			num_enriched_cand_genes += 1
#	#print num_enriched_cand_genes
#	return num_enriched_cand_genes
#
#
#
#def cand_genes_enrichment(phen_id, cgl_id, phenotypeFile="/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phenotypes_all_raw_081009.tsv"):
#	#Get p-values, and their pos.
#	phed = pd.readPhenotypeFile(phenotypeFile, delimiter='\t')
#	phenName = phed.getPhenotypeName(phen_id)
#	res_path = "/Users/bjarnivilhjalmsson/Projects/Data/gwas_results/"
#	resultType = gr.ResultType("KW", ".pvals", "raw", res_path + "kw_results/")
#	result_file = resultType.getFileName(phed, phen_id)
#	kw_result = gr.Result(result_file, name=str(resultType) + "_" + phenName, resultType=resultType, phenotypeID=phen_id)
#	kw_result.filterMARF(0.1)
#	kw_result.negLogTransform()
#
#	resultType = gr.ResultType("Emma", ".pvals", "logTransform", res_path + "emma_results/")
#	result_file = resultType.getFileName(phed, phen_id)
#	emma_result = gr.Result(result_file, name=str(resultType) + "_" + phenName, resultType=resultType, phenotypeID=phen_id)
#	emma_result.filterMARF(0.1)
#	emma_result.negLogTransform()
#
#
#	#Get genes, (using the gene structure)
#	genes = gr.get_gene_list(passwd='bjazz32')
#	#Mark Cand. genes.
#	cand_genes = gr.getCandidateGeneList(cgl_id, passwd='bjazz32')
#
#	cand_gene_indices = []
#	for i, g in enumerate(genes):
#		for cg in cand_genes:
#			if g.dbRef == cg.dbRef:
#				#print g.dbRef, cg.dbRef
#				cand_gene_indices.append(i)
#				break
#	num_cand_genes = len(cand_genes)
#	print num_cand_genes, cand_gene_indices
#
#	#Calc. statistics.
#	kw_mean_pval = sum(kw_result.scores) / float(len(kw_result.scores))
#	kw_s = num_enrich_cand_genes(cand_gene_indices, genes, kw_result, 2)
#	print "Observed:", kw_s
#
#	num_perm = 1000
#	stat_list = []
#	for i in range(num_perm):
#		cand_gene_indices = random.sample(range(len(genes)), num_cand_genes)
#		cand_gene_indices.sort()
#		#print cand_gene_indices
#		kw_s = num_enrich_cand_genes(cand_gene_indices, genes, kw_result, 2)
#		stat_list.append(kw_s)
#		#print kw_s		
#	plt.hist(stat_list, 50)
#	plt.savefig("/Users/bjarnivilhjalmsson/tmp/kw_or.pdf")
#
#	emma_mean_pval = sum(emma_result.scores) / float(len(emma_result.scores))
#	emma_s = num_enrich_cand_genes(cand_gene_indices, genes, emma_result, 1)
#	print "Observed:", emma_s
#	stat_list = []
#	for i in range(num_perm):
#		cand_gene_indices = random.sample(range(len(genes)), num_cand_genes)
#		cand_gene_indices.sort()
#		emma_s = num_enrich_cand_genes(cand_gene_indices, genes, emma_result, 0.5)
#		stat_list.append(emma_s)
#	plt.hist(stat_list, 50)
#	plt.savefig("/Users/bjarnivilhjalmsson/tmp/emma_or.pdf")
#
#
#
#def mean_pvalues(result, genes, gene_window_size):
#	"""
#	Genes are assumed to be sorted by starting position.
#	"""
#	chr_pos_list = result.getChrPos()
#	mean_log_pvals = []
#	last_k = 0
#	for j, g in enumerate(genes):
#		g_start = g.startPos - gene_window_size
#		g_end = g.endPos + gene_window_size
#		k = last_k
#		(chr, pos) = chr_pos_list[k]
#		while k < len(chr_pos_list) - 1 and (chr < g.chromosome or (chr == g.chromosome and pos < g_start)):
#			k += 1
#			(chr, pos) = chr_pos_list[k]
#		last_k = k
#		pval_sum = 0
#		while k < len(chr_pos_list) and (chr == g.chromosome and g_start <= pos < g_end):
#			pval_sum += result.scores[k]
#			k += 1
#			if k < len(chr_pos_list):
#				(chr, pos) = chr_pos_list[k]
#		if k - last_k > 0:
#			pval = pval_sum / float(k - last_k)
#			mean_log_pvals.append(pval)
#	return mean_log_pvals
#
#
#
#
#
#def plot_enrichment_stat(phenotype_ids, cgl_id, phenotypeFile="/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phenotypes_all_raw_081009.tsv",
#				gene_window_size=5000):
#	"""
#	Plotting enrichment statistics for different phenotypes
#	"""
#
#	#Loading the phenotypes
#	phed = pd.readPhenotypeFile(phenotypeFile, delimiter='\t')
#	#Results settings
#	res_path = "/Users/bjarnivilhjalmsson/Projects/Data/gwas_results/"
#	result_types = [gr.ResultType("KW", ".pvals", "raw", res_path + "kw_results/"), gr.ResultType("Emma", ".pvals", "logTransform", res_path + "emma_results/")]
#	result_labels = ["KW", "EMMA"]
#	min_marf = 0.1
#
#	#Loading one result, to get SNPs..	
#	result_file = result_types[0].getFileName(phed, 1)
#	r = gr.Result(result_file, name=str(result_types[0]) + "_" + phed.getPhenotypeName(1), resultType=result_types[0], phenotypeID=1)
#
#	#Loading and filtering cand. genes
#	cand_genes = gr.getCandidateGeneList(cgl_id, passwd='bjazz32')
#	new_cand_genes = []
#	i = 0
#	cur_pos = r.positions[i]
#	cur_chr = r.chromosomes[i]
#	for g in cand_genes:
#		while r.chromosomes[i] < g.chromosome or (r.chromosomes[i] == g.chromosome and r.positions[i] < g.startPos - gene_window_size):
#			i += 1
#		if r.chromosomes[i] == g.chromosome and r.positions[i] < g.endPos + gene_window_size:
#			new_cand_genes.append(g)
#	print "len(new_cand_genes):", len(new_cand_genes), ", len(cand_genes):", len(cand_genes)
#	cand_genes = new_cand_genes
#
#	#Loading and filtering the full gene list
#	all_genes = gr.get_gene_list(passwd='bjazz32')
#	new_all_genes = []
#	i = 0
#	cur_pos = r.positions[i]
#	cur_chr = r.chromosomes[i]
#	for g in all_genes:
#		while r.chromosomes[i] < g.chromosome or (r.chromosomes[i] == g.chromosome and r.positions[i] < g.startPos - gene_window_size):
#			i += 1
#		if r.chromosomes[i] == g.chromosome and r.positions[i] < g.endPos + gene_window_size:
#			new_all_genes.append(g)
#	print "len(new_all_genes):", len(new_all_genes), ", len(all_genes):", len(all_genes)
#	all_genes = new_all_genes
#
#
#
#
#	#For each phenotype
#	for pi in phenotype_ids:
#		phenName = phed.getPhenotypeName(pi)
#		ers = []
#		cutoffs_list = []
#		pvals_list = []
#		exp_obs_diff_list = []
#		#For each result:
#		min_x = 100
#		max_x = 0
#		for rt_i, rt in enumerate(result_types):
#			result_file = rt.getFileName(phed, pi)
#			#Load result
#			result = gr.Result(result_file, name=str(rt) + "_" + phenName, resultType=rt, phenotypeID=pi)
#
#			result.filterMARF(min_marf)
#			result.negLogTransform()
#
#			#Calculate the mean pvalues for the two gene groups.
#			cg_mean_pvals = mean_pvalues(result, cand_genes, gene_window_size)
#			g_mean_pvals = mean_pvalues(result, all_genes, gene_window_size)
#
#			mean_pval_diff = sum(cg_mean_pvals) / len(cg_mean_pvals) - sum(g_mean_pvals) / len(g_mean_pvals)
#			print "mean_pval_diff:", mean_pval_diff
#
#			plt.figure(figsize=(8, 5))
#			plt.hist(g_mean_pvals, 1000)
#			plt.hist(cg_mean_pvals, 20)
#			plt.savefig("/Users/bjarnivilhjalmsson/tmp/mean_pval_p" + str(pi) + "_rt" + str(rt_i) + ".pdf")
#			plt.clf()
#
#
#
#
#def plot_enrichment_histogram(phenotype_ids, cgl_id, phenotypeFile="/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phenotypes_all_raw_081009.tsv",
#				gene_window_size=5000, top_quantiles=[0.075, 0.09, 0.1, 0.12]):
#	"""
#	Plotting enrichment ratios for different phenotypes
#	"""
#	fontProp = matplotlib.font_manager.FontProperties(size=10)
#
#
#	phed = pd.readPhenotypeFile(phenotypeFile, delimiter='\t')
#	res_path = "/Users/bjarnivilhjalmsson/Projects/Data/gwas_results/"
#	result_types = [gr.ResultType("KW", ".pvals", "raw", res_path + "kw_results/"), gr.ResultType("Emma", ".pvals", "logTransform", res_path + "emma_results/")]
#	result_labels = ["KW", "EMMA"]
#	min_marf = 0.1
#
#	result_file = result_types[0].getFileName(phed, 1)
#	r = gr.Result(result_file, name=str(result_types[0]) + "_" + phed.getPhenotypeName(1), resultType=result_types[0], phenotypeID=1)
#	r.filterMARF(min_marf)
#
#	cand_genes = gr.getCandidateGeneList(cgl_id, passwd='bjazz32')
#	new_cand_genes = []
#	i = 0
#	cur_pos = r.positions[i]
#	cur_chr = r.chromosomes[i]
#	for g in cand_genes:
#		while r.chromosomes[i] < g.chromosome or (r.chromosomes[i] == g.chromosome and r.positions[i] < g.startPos - gene_window_size):
#			i += 1
#		if r.chromosomes[i] == g.chromosome and r.positions[i] < g.endPos + gene_window_size:
#			new_cand_genes.append(g)
#	print "len(new_cand_genes):", len(new_cand_genes), ", len(cand_genes):", len(cand_genes)
#	cand_genes = new_cand_genes
#
#	all_genes = gr.get_gene_list(passwd='bjazz32')
#	new_all_genes = []
#	i = 0
#	cur_pos = r.positions[i]
#	cur_chr = r.chromosomes[i]
#	for g in all_genes:
#		while r.chromosomes[i] < g.chromosome or (r.chromosomes[i] == g.chromosome and r.positions[i] < g.startPos - gene_window_size):
#			i += 1
#		if r.chromosomes[i] == g.chromosome and r.positions[i] < g.endPos + gene_window_size:
#			new_all_genes.append(g)
#	print "len(new_all_genes):", len(new_all_genes), ", len(all_genes):", len(all_genes)
#	all_genes = new_all_genes
#	print "len(new_all_genes):", len(new_all_genes), ", len(all_genes):", len(all_genes)
#
#
#	cand_gene_indices = []
#	for i, g in enumerate(all_genes):
#		for cg in cand_genes:
#			if g.dbRef == cg.dbRef:
#				#print g.dbRef, cg.dbRef
#				cand_gene_indices.append(i)
#				break
#	num_genes = len(all_genes)
#	num_cand_genes = len(cand_genes)
#	num_non_cand_genes = num_genes - num_cand_genes
#	cand_gene_frac = num_cand_genes / float(num_genes)
#
#	num_steps = len(top_quantiles)
#
#	oeds = []
#
#	#For each phenotype
#	for pi in phenotype_ids:
#		phenName = phed.getPhenotypeName(pi)
#		ers = []
#		cutoffs_list = []
#		pvals_list = []
#		exp_obs_diff_list = []
#		#For each result:
#		min_x = 100
#		max_x = 0
#		for rt_i, rt in enumerate(result_types):
#			result_file = rt.getFileName(phed, pi)
#			#Load result
#			result = gr.Result(result_file, name=str(rt) + "_" + phenName, resultType=rt, phenotypeID=pi)
#
#			result.filterMARF(min_marf)
#			result.negLogTransform()
#
#			cutoffs = []
#			result._rankScores_()
#			for q in top_quantiles:
#				cutoffs.append(result.scores[result.orders[int(q * len(result.scores))]])
#
#
#			#print max_log_pval, min_log_pval, log_pval_step
#			result.filterScoreCutoff(cutoffs[-1])
#			print "Number of p-vals:", len(result.scores)
#			min_x = min(min_x, cutoffs[-1])
#			max_x = max(max_x, cutoffs[0])
#			#print max_x, min_x
#
#			#Pre process the p-values
#			pval_dict = [[] for i in range(num_steps)]
#			chr_pos_list = result.getChrPos()
#			for i, log_pval in enumerate(result.scores):
#				cat = num_steps - 1
#				while cat > 0 and log_pval > cutoffs[cat - 1]:
#					cat = cat - 1
#				pval_dict[cat].append(chr_pos_list[i])
#			for i, c in enumerate(cutoffs):
#				print "cutoff", c, ", num pvals:", len(pval_dict[i])
#
#			#calculate enrichment ratios
#			def calc_enrichement_ratios(cand_gene_indices, default_er= -1.0):
#				enrichm_ratios = [0.0] * num_steps
#				num_cand_found = [0] * num_steps
#				gene_indices_found = [set() for i in range(num_steps)]
#				cand_indices_found = [set() for i in range(num_steps)]
#				num_genes_found = [0] * num_steps
#				pvals = []
#				exp_obs_cand_gene_diff = []
#				#pval_counts = np.zeros(len(cutoffs)*len(all_genes)).reshape(len(cutoffs),len(all_genes))
#				for i, chr_pos_list in enumerate(pval_dict):
#					k = 0
#					if len(chr_pos_list) > 0:
#						for j, g in enumerate(all_genes):
#							g_start = g.startPos - gene_window_size
#							g_end = g.endPos + gene_window_size
#							(chr, pos) = chr_pos_list[k]
#							while k < len(chr_pos_list) - 1 and (chr < g.chromosome or (chr == g.chromosome and pos < g_start)):
#								k += 1
#								(chr, pos) = chr_pos_list[k]
#
#							if chr == g.chromosome and g_start <= pos < g_end:
#								gene_indices_found[i].add(j)
#								if j in cand_gene_indices:
#									cand_indices_found[i].add(j)
#					if i > 0:
#						gene_indices_found[i] = gene_indices_found[i].union(gene_indices_found[i - 1])
#						cand_indices_found[i] = cand_indices_found[i].union(cand_indices_found[i - 1])
#
#					num_genes_found[i] = len(gene_indices_found[i])
#					num_cand_found[i] = len(cand_indices_found[i])
#					if num_genes_found[i] >= 5 and num_cand_found[i] >= 1:
#						enrichm_ratios[i] = (num_cand_found[i] / float(num_cand_genes)) / (num_genes_found[i] / float(num_genes))
#						pval = 1.0 - rpy.r.phyper(num_cand_found[i], num_cand_genes, num_non_cand_genes, num_genes_found[i])
#						eod = num_cand_found[i] - cand_gene_frac * num_genes_found[i]
#						print num_cand_found[i] - 1, num_cand_genes, num_non_cand_genes, num_genes_found[i], pval, eod
#						if pval == 0.0:
#							pvals.append(1.0)
#						else:
#							pvals.append(pval)
#						exp_obs_cand_gene_diff.append(eod)
#					else:
#						enrichm_ratios[i] = default_er
#						exp_obs_cand_gene_diff.append(0.0)
#						pvals.append(1.0)
#
#
#				results = {"enrichment_ratios": enrichm_ratios, "num_genes_found": num_genes_found, "num_cand_found": num_cand_found, "pvals":pvals, "exp_obs_cand_gene_diff":exp_obs_cand_gene_diff}
#				return results
#
#			enrichment_res = calc_enrichement_ratios(cand_gene_indices)
#			pvals = enrichment_res["pvals"]
#			pvals_list.append(pvals)
#
#			eod = enrichment_res["exp_obs_cand_gene_diff"]
#			oer = enrichment_res["enrichment_ratios"]
#			print oer
#			print pvals
#			print eod
#
#			ers.append(oer)
#			cutoffs_list.append(cutoffs)
#			exp_obs_diff_list.append(eod)
#		oeds.append(exp_obs_diff_list)
#
#	print oeds #[pi][rt][tq]
#	oeds_list = [] #[tq][rt][pi]
#	for k in range(len(top_quantiles)):
#		x_s = []
#		for j in range(len(result_types)):
#			x = []
#			for i in range(len(phenotype_ids)):
#				x.append(oeds[i][j][k])
#			x_s.append(x)
#		oeds_list.append(x_s)
#
#	oeds = oeds_list
#
#	x_s = []
#	for i in range(len(result_types)):
#		x_s.append(range(i, len(phenotype_ids) * (len(result_types) + 1), len(result_types) + 1))
#	print x_s
#
#	colors = ["blue", "red", "green"]
#	labels = ["KW", "EMMA"]
#	xticks_pos = range(len(result_types) / 2, len(phenotype_ids) * (len(result_types) + 1), len(result_types) + 1)
#	for i in range(len(top_quantiles)):
#		plt.figure(figsize=(8, 5))
#		for j in range(len(result_types)):
#			plt.bar(x_s[j], oeds[i][j], color=colors[j], label=labels[j])
#		plt.legend()
#		plt.xticks(xticks_pos, phenotype_ids)
#		plt.savefig("/Users/bjarnivilhjalmsson/tmp/EOD_hist_q" + str(top_quantiles[i]) + "_w" + str(gene_window_size) + ".pdf")
#		plt.clf()
#
#
#
#def plot_enrichment_ratios(phenotype_ids, cgl_id, phenotypeFile="/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phenotypes_all_raw_081009.tsv",
#				gene_window_size=5000, use_fixed_pvalue_limits=True, max_log_pval=10.0, min_log_pval=1.0, num_steps=30,
#				min_num_pvals=10, max_num_pvals=50000):
#	"""
#	Plotting enrichment ratios for different phenotypes
#	"""
#	fontProp = matplotlib.font_manager.FontProperties(size=10)
#
#
#	phed = pd.readPhenotypeFile(phenotypeFile, delimiter='\t')
#	res_path = "/Users/bjarnivilhjalmsson/Projects/Data/gwas_results/"
#	result_types = [gr.ResultType("KW", ".pvals", "raw", res_path + "kw_results/"), gr.ResultType("Emma", ".pvals", "logTransform", res_path + "emma_results/")]
#	result_labels = ["KW", "EMMA"]
#	min_marf = 0.1
#	if use_fixed_pvalue_limits:
#		log_pval_step = (max_log_pval - min_log_pval) / (num_steps - 1)
#		cutoffs = np.arange(max_log_pval, min_log_pval - log_pval_step, -log_pval_step).tolist()
#
#
#	result_file = result_types[0].getFileName(phed, 1)
#	r = gr.Result(result_file, name=str(result_types[0]) + "_" + phed.getPhenotypeName(1), resultType=result_types[0], phenotypeID=1)
#	r.filterMARF(min_marf)
#
#	cand_genes = gr.getCandidateGeneList(cgl_id, passwd='bjazz32')
#	new_cand_genes = []
#	i = 0
#	cur_pos = r.positions[i]
#	cur_chr = r.chromosomes[i]
#	for g in cand_genes:
#		while r.chromosomes[i] < g.chromosome or (r.chromosomes[i] == g.chromosome and r.positions[i] < g.startPos - gene_window_size):
#			i += 1
#		if r.chromosomes[i] == g.chromosome and r.positions[i] < g.endPos + gene_window_size:
#			new_cand_genes.append(g)
#	print "len(new_cand_genes):", len(new_cand_genes), ", len(cand_genes):", len(cand_genes)
#	cand_genes = new_cand_genes
#
#	all_genes = gr.get_gene_list(passwd='bjazz32')
#	new_all_genes = []
#	i = 0
#	cur_pos = r.positions[i]
#	cur_chr = r.chromosomes[i]
#	for g in all_genes:
#		while r.chromosomes[i] < g.chromosome or (r.chromosomes[i] == g.chromosome and r.positions[i] < g.startPos - gene_window_size):
#			i += 1
#		if r.chromosomes[i] == g.chromosome and r.positions[i] < g.endPos + gene_window_size:
#			new_all_genes.append(g)
#	print "len(new_all_genes):", len(new_all_genes), ", len(all_genes):", len(all_genes)
#	all_genes = new_all_genes
#	print "len(new_all_genes):", len(new_all_genes), ", len(all_genes):", len(all_genes)
#
#
#	cand_gene_indices = []
#	for i, g in enumerate(all_genes):
#		for cg in cand_genes:
#			if g.dbRef == cg.dbRef:
#				#print g.dbRef, cg.dbRef
#				cand_gene_indices.append(i)
#				break
#	num_genes = len(all_genes)
#	num_cand_genes = len(cand_genes)
#	num_non_cand_genes = num_genes - num_cand_genes
#	cand_gene_frac = num_cand_genes / float(num_genes)
#
#	#For each phenotype
#	for pi in phenotype_ids:
#		phenName = phed.getPhenotypeName(pi)
#		ers = []
#		cutoffs_list = []
#		pvals_list = []
#		exp_obs_diff_list = []
#		#For each result:
#		min_x = 100
#		max_x = 0
#		for rt in result_types:
#			result_file = rt.getFileName(phed, pi)
#			#Load result
#			result = gr.Result(result_file, name=str(rt) + "_" + phenName, resultType=rt, phenotypeID=pi)
#
#			result.filterMARF(min_marf)
#			result.negLogTransform()
#
#			if not use_fixed_pvalue_limits:
#
#				result._rankScores_()
#				min_log_pval = result.scores[result.orders[max_num_pvals]]
#				max_log_pval = result.scores[result.orders[min_num_pvals]]
#				log_pval_step = (max_log_pval - min_log_pval) / (num_steps - 1.0)
#				cutoffs = np.arange(max_log_pval, min_log_pval - log_pval_step, -log_pval_step).tolist()
#
#				if len(cutoffs) > num_steps:#Strange python bug?
#					cutoffs = cutoffs[:num_steps]
#				#print cutoffs, len(cutoffs), num_steps
#
#			#print max_log_pval, min_log_pval, log_pval_step
#			result.filterScoreCutoff(min_log_pval)
#			print "Number of p-vals:", len(result.scores)
#			min_x = min(min_x, cutoffs[-1])
#			max_x = max(max_x, cutoffs[0])
#			#print max_x, min_x
#
#			#Pre process the p-values
#			pval_dict = [[] for i in range(num_steps)]
#			chr_pos_list = result.getChrPos()
#			for i, log_pval in enumerate(result.scores):
#				cat = num_steps - 1
#				while cat > 0 and log_pval > cutoffs[cat - 1]:
#					cat = cat - 1
#				pval_dict[cat].append(chr_pos_list[i])
#			for i, c in enumerate(cutoffs):
#				print "cutoff", c, ", num pvals:", len(pval_dict[i])
#
#			#calculate enrichment ratios
#			def calc_enrichement_ratios(cand_gene_indices, default_er= -1.0):
#				enrichm_ratios = [0.0] * num_steps
#				num_cand_found = [0] * num_steps
#				gene_indices_found = [set() for i in range(num_steps)]
#				cand_indices_found = [set() for i in range(num_steps)]
#				num_genes_found = [0] * num_steps
#				pvals = []
#				exp_obs_cand_gene_diff = []
#				#pval_counts = np.zeros(len(cutoffs)*len(all_genes)).reshape(len(cutoffs),len(all_genes))
#				for i, chr_pos_list in enumerate(pval_dict):
#					k = 0
#					if len(chr_pos_list) > 0:
#						for j, g in enumerate(all_genes):
#							g_start = g.startPos - gene_window_size
#							g_end = g.endPos + gene_window_size
#							(chr, pos) = chr_pos_list[k]
#							while k < len(chr_pos_list) - 1 and (chr < g.chromosome or (chr == g.chromosome and pos < g_start)):
#								k += 1
#								(chr, pos) = chr_pos_list[k]
#
#							if chr == g.chromosome and g_start <= pos < g_end:
#								gene_indices_found[i].add(j)
#								if j in cand_gene_indices:
#									cand_indices_found[i].add(j)
#					if i > 0:
#						gene_indices_found[i] = gene_indices_found[i].union(gene_indices_found[i - 1])
#						cand_indices_found[i] = cand_indices_found[i].union(cand_indices_found[i - 1])
#
#					num_genes_found[i] = len(gene_indices_found[i])
#					num_cand_found[i] = len(cand_indices_found[i])
#					if num_genes_found[i] >= 5 and num_cand_found[i] >= 1:
#						enrichm_ratios[i] = (num_cand_found[i] / float(num_cand_genes)) / (num_genes_found[i] / float(num_genes))
#						pval = 1.0 - rpy.r.phyper(num_cand_found[i], num_cand_genes, num_non_cand_genes, num_genes_found[i])
#						eod = num_cand_found[i] - cand_gene_frac * num_genes_found[i]
#						print num_cand_found[i] - 1, num_cand_genes, num_non_cand_genes, num_genes_found[i], pval, eod
#						if pval == 0.0:
#							pvals.append(1.0)
#						else:
#							pvals.append(pval)
#						exp_obs_cand_gene_diff.append(eod)
#					else:
#						enrichm_ratios[i] = default_er
#						exp_obs_cand_gene_diff.append(0.0)
#						pvals.append(1.0)
#
#
#				results = {"enrichment_ratios": enrichm_ratios, "num_genes_found": num_genes_found, "num_cand_found": num_cand_found, "pvals":pvals, "exp_obs_cand_gene_diff":exp_obs_cand_gene_diff}
#				return results
#
#			enrichment_res = calc_enrichement_ratios(cand_gene_indices)
#			pvals = enrichment_res["pvals"]
#			pvals_list.append(pvals)
#
#			eod = enrichment_res["exp_obs_cand_gene_diff"]
#			oer = enrichment_res["enrichment_ratios"]
#			print oer
#			print pvals
#			print eod
#
#			ers.append(oer)
#			cutoffs_list.append(cutoffs)
#			exp_obs_diff_list.append(eod)
#
#
#		plt.figure(figsize=(8, 5))
#		max_y = 0
#		for i, er in enumerate(ers):
#			x_s = []
#			y_s = []
#			cutoffs = cutoffs_list[i]
#			for j in range(num_steps):
#				if er[j] != -1:
#					x_s.append(cutoffs[j])
#					y_s.append(er[j])
#			plt.plot(x_s, y_s, '-', label=result_labels[i])
#			max_y = max(max_y, max(y_s))
#		y_range = max_y
#		x_range = max_x - min_x
#		plt.axis([min_x - 0.025 * x_range, max_x + 0.1 * x_range, -0.025 * y_range, 1.15 * y_range])
#		plt.axhline(y=1.0, ls=":", color='black')
#		plt.xlabel("Association mapping pvalue ($-log_{10}(p)$)")
#		plt.ylabel("Enrichment ratio")
#		plt.legend(numpoints=1, handlelen=0.005, prop=fontProp)
#		plt.savefig("/Users/bjarnivilhjalmsson/tmp/ER_cgl" + str(cgl_id) + "_w" + str(gene_window_size) + "_pi" + str(pi) + ".pdf")
#
#		#Plotting p-values 
#		plt.clf()
#		max_y = 0
#		for i, pvals in enumerate(pvals_list):
#			er = ers[i]
#			x_s = []
#			y_s = []
#			cutoffs = cutoffs_list[i]
#			for j in range(len(cutoffs)):
#				if er[j] != -1:
#					x_s.append(cutoffs[j])
#					y_s.append(-(math.log10(pvals[j])))
#			plt.plot(x_s, y_s, '-', label=result_labels[i])
#			max_y = max(max_y, max(y_s))
#		y_range = max_y
#		plt.axis([min_x - 0.025 * x_range, max_x + 0.1 * x_range, -0.025 * y_range, 1.15 * y_range])
#		plt.axhline(y= -math.log10(0.05), ls=":", color='black')
#		plt.xlabel("Association mapping pvalue ($-log_{10}(p)$)")
#		plt.ylabel("Enrichment ratio pvalue ($-log_{10}(p)$)")
#		plt.legend(numpoints=1, handlelen=0.005, prop=fontProp)
#		plt.savefig("/Users/bjarnivilhjalmsson/tmp/ER_cgl" + str(cgl_id) + "_w" + str(gene_window_size) + "_pi" + str(pi) + "_pvals.pdf")
#		plt.clf()
#
#		#Plotting expected-observed cand. gene number differences 
#		plt.clf()
#		max_y = 0
#		min_y = 0
#		for i, eods in enumerate(exp_obs_diff_list):
#			er = ers[i]
#			x_s = []
#			y_s = []
#			cutoffs = cutoffs_list[i]
#			for j in range(len(cutoffs)):
#				if er[j] != -1:
#					x_s.append(cutoffs[j])
#					y_s.append(eods[j])
#			plt.plot(x_s, y_s, '-', label=result_labels[i])
#			max_y = max(max_y, max(y_s))
#			min_y = min(min_y, min(y_s))
#		y_range = max_y - min_y
#		plt.axis([min_x - 0.025 * x_range, max_x + 0.1 * x_range, min_y - 0.025 * y_range, max_y + 0.15 * y_range])
#		plt.axhline(y=0.0, ls=":", color='black')
#		plt.xlabel("Association mapping pvalue ($-log_{10}(p)$)")
#		plt.ylabel("Expected-observed cand. gene number difference")
#		plt.legend(numpoints=1, handlelen=0.005, prop=fontProp)
#		plt.savefig("/Users/bjarnivilhjalmsson/tmp/EOD_cgl" + str(cgl_id) + "_w" + str(gene_window_size) + "_pi" + str(pi) + ".pdf")
#		plt.clf()

def plot_gene_distances(cgl_id=145):
	cand_genes = gr.getCandidateGeneList(cgl_id, passwd='bjazz32')
	cand_gene_distances = []
	for i in range(len(cand_genes)):
		g1 = cand_genes[i]
		for j in range(i + 1, len(cand_genes)):
			g2 = cand_genes[j]
			if g1.chromosome == g2.chromosome:
				dist = max(max(g1.startPos - g2.endPos, 0), max(g2.startPos - g1.endPos, 0))
				if dist < 500000:
					cand_gene_distances.append(dist)
	print cand_gene_distances
	plt.figure(figsize=(8, 5))
	plt.hist(cand_gene_distances, bins=100)
	plt.savefig("/Users/bjarnivilhjalmsson/tmp/Gene_distances_cgl" + str(cgl_id) + ".pdf")
	plt.clf()
	all_genes = gr.get_gene_list(passwd='bjazz32')
	all_gene_distances = []
	for i in range(len(all_genes)):
		g1 = all_genes[i]
		for j in range(i + 1, len(all_genes)):
			g2 = all_genes[j]
			if g1.chromosome == g2.chromosome:
				dist = max(max(g1.startPos - g2.endPos, 0), max(g2.startPos - g1.endPos, 0))
				if dist < 500000:
					all_gene_distances.append(dist)
		print i
	plt.hist(all_gene_distances, bins=100)
	plt.savefig("/Users/bjarnivilhjalmsson/tmp/Gene_distances.pdf")
	(D_stat, pval) = stats.ks_2samp(all_gene_distances, cand_gene_distances)
	print pval




#def plot_gene_pvalue_count(cgl_id=145, window_size=20000):
#	chr_pos_list = dataParsers.parse_chr_pos_list("/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_f13_012609.csv")
#	cand_genes = gr.getCandidateGeneList(cgl_id, passwd='bjazz32')
#	cand_gene_snp_counts = []
#	next_start = 0
#	for i in range(len(cand_genes)):
#		g = cand_genes[i]
#		j = next_start
#		chr_pos = chr_pos_list[j]
#		while  chr_pos < (g.chromosome, g.startPos - window_size):
#			j += 1
#			if j == len(chr_pos_list):
#				break
#			else:
#				chr_pos = chr_pos_list[j]
#		next_start = j
#		snp_count = 0
#		while chr_pos <= (g.chromosome, g.endPos + window_size):
#			snp_count += 1
#			j += 1
#			if j == len(chr_pos_list):
#				break
#			else:
#				chr_pos = chr_pos_list[j]
#		cand_gene_snp_counts.append(snp_count)
#
#	plt.figure(figsize=(8, 5))
#	plt.hist(cand_gene_snp_counts, bins=100)
#	plt.savefig("/Users/bjarnivilhjalmsson/tmp/SNP_counts_cgl" + str(cgl_id) + "_w" + str(window_size) + ".pdf")
#	plt.clf()
#
#	all_genes = gr.get_gene_list(passwd='bjazz32')
#	all_gene_snp_counts = []
#	next_start = 0
#	for i in range(len(all_genes)):
#		g = all_genes[i]
#		j = next_start
#		chr_pos = chr_pos_list[j]
#		while  chr_pos < (g.chromosome, g.startPos - window_size):
#			j += 1
#			if j == len(chr_pos_list):
#				break
#			else:
#				chr_pos = chr_pos_list[j]
#		next_start = j
#		snp_count = 0
#		while chr_pos <= (g.chromosome, g.endPos + window_size):
#			snp_count += 1
#			j += 1
#			if j == len(chr_pos_list):
#				break
#			else:
#				chr_pos = chr_pos_list[j]
#		all_gene_snp_counts.append(snp_count)
#
#	plt.hist(all_gene_snp_counts, bins=100)
#	plt.savefig("/Users/bjarnivilhjalmsson/tmp/SNP_counts_w" + str(window_size) + ".pdf")
#	(D_stat, pval) = stats.ks_2samp(all_gene_snp_counts, cand_gene_snp_counts)
#	print pval



def calc_enrichment(all_genes, cg_indices, regions):
	"""
	Calculate the enrichment (efficiently iterating over all genes.)
	"""
	num_genes = len(all_genes)
	num_cand_genes = len(cg_indices)
	g_iter = enumerate(all_genes)
	g_i, g = g_iter.next()
	g_end_chr_pos = sp.array([g.chromosome, g.endPos])
	g_start_chr_pos = sp.array([g.chromosome, g.startPos])
	r_iter = enumerate(regions)
	r_i, r = r_iter.next()
	r_start_chr_pos = r[0]
	r_end_chr_pos = r[1]
	num_close_genes = 0
	num_close_cand_genes = 0
	cg_iter = enumerate(cg_indices)
	cg_ii, cg_i = cg_iter.next()
	obs_cg_indices = []
	while g_i < num_genes - 1 and r_i < len(regions) - 1:
		#count until overlap
		#if g_i % 100 == 0: print g_i
		while g_i < num_genes - 1 and tuple(g_end_chr_pos) < tuple(r_start_chr_pos):
			g_i, g = g_iter.next()
			g_end_chr_pos = sp.array([g.chromosome, g.endPos])
			g_start_chr_pos = sp.array([g.chromosome, g.startPos])

		while r_i < len(regions) - 1 and  tuple(r_end_chr_pos) < tuple(g_start_chr_pos):
			r_i, r = r_iter.next()
			r_start_chr_pos = r[0]
			r_end_chr_pos = r[1]

		while g_i < num_genes - 1 and tuple(g_end_chr_pos) >= tuple(r_start_chr_pos) \
				and tuple(r_end_chr_pos) >= tuple(g_start_chr_pos):
			#there is an overlap
			num_close_genes += 1
			while cg_ii < num_cand_genes - 1 and cg_i < g_i:
				cg_ii, cg_i = cg_iter.next()
			if g_i == cg_i:
				num_close_cand_genes += 1
				obs_cg_indices.append(g_i)
			g_i, g = g_iter.next()
			g_end_chr_pos = sp.array([g.chromosome, g.endPos])
			g_start_chr_pos = sp.array([g.chromosome, g.startPos])
		#print g_i, r_i
		#print num_close_genes, num_close_cand_genes
		#pdb.set_trace()

	#r1 = (num_close_cand_genes / float(num_close_genes))
	#r2 = (num_cand_genes / float(num_genes))
	return (num_close_cand_genes, num_close_genes, obs_cg_indices)


def get_chi_square_pval(num_close_cand_genes, num_close_genes, num_cand_genes, num_genes):
	num_exp_ccg = (num_close_genes) * num_cand_genes / float(num_genes)
	num_exp_cg = (num_close_genes) * (1 - num_cand_genes / float(num_genes))
	num_exp_dcg = (float(num_genes) - num_close_genes) * (num_cand_genes / float(num_genes))
	num_exp_dg = (float(num_genes) - num_close_genes) * (1 - num_cand_genes / float(num_genes))
	num_obs_ccg = num_close_cand_genes
	num_obs_cg = num_close_genes - num_close_cand_genes
	num_obs_dcg = num_cand_genes - num_close_cand_genes
	num_obs_dg = num_genes - num_close_genes
	f_obs = sp.array([num_obs_ccg, num_obs_cg, num_obs_dcg, num_obs_dg])
	f_exp = sp.array([num_exp_ccg, num_exp_cg, num_exp_dcg, num_exp_dg])
	chi_sq_stat, chi_sq_pval = st.chisquare(f_obs, f_exp, 2)
	if num_obs_ccg >= num_exp_ccg:
		chi_sq_pval = chi_sq_pval / 2
	else:
		chi_sq_pval = 1 - chi_sq_pval / 2
	print 'Cand. gene enrichment p-value from Chi-square test:', chi_sq_pval
	return chi_sq_pval, chi_sq_stat


def get_gene_perm_pval(obs_stat, regions, all_genes, cand_gene_indices, num_perm=500, early_stop_threshold=25):
	print "Doing gene permutations"
	num_genes = len(all_genes)
	num_cand_genes = len(cand_gene_indices)
	r2 = num_cand_genes / float(num_genes)
	perm_stats = []
	if obs_stat != float('-inf'):
		sign_count = 0
		for perm_i in range(num_perm):
			perm_cand_gene_indices = random.sample(range(num_genes), num_cand_genes)
			perm_cand_gene_indices.sort()
			perm_enrichments = calc_enrichment(all_genes, perm_cand_gene_indices, regions)
			r1 = perm_enrichments[0] / float(perm_enrichments[1])
			if r1 == 0.0:
				perm_stat = float('-inf')
			else:
				perm_stat = sp.log(r1 / r2)
			perm_stats.append(perm_stat)
			if perm_stat >= obs_stat:
				sign_count += 1
				if sign_count == early_stop_threshold:
					break
			sys.stdout.write('.')
			sys.stdout.flush()
		sys.stdout.write('\n')
		sys.stdout.flush()
		perm_stats.sort()
		p_val = sign_count / float(perm_i + 1)
	else:
		p_val = 1.0
	print 'Permutation p-value estimate: % f' % p_val
	return p_val, perm_stats



def get_snps_perm_pval(obs_stat, regions, all_genes, cand_gene_indices, chrom_ends, num_perm=500,
		early_stop_threshold=25):
	print "Doing SNPs permutations"
	num_genes = len(all_genes)
	num_cand_genes = len(cand_gene_indices)
	r2 = num_cand_genes / float(num_genes)
	def _get_perm_regions_(region_dict):
		new_regions = []
		for chrom in region_dict:
			chrom_end = chrom_ends[chrom - 1]
			regions = region_dict[chrom]
			while True:
				shift = int(random.random()*chrom_end)
				for region in regions:
					if (region[0][1] + shift) < chrom_end \
							and (region[1][1] + shift) > chrom_end:
						break #for loop
				else:
					break #while loop
			start_pos_list = []
			for i, region in enumerate(regions):
				start_pos = (region[0][1] + shift) % chrom_end
				region[0][1] = start_pos
				region[1][1] = (region[1][1] + shift) % chrom_end
				start_pos_list.append((start_pos, i))
			start_pos_list.sort()
			for sp, i in start_pos_list:
				new_regions.append(regions[i])
		return new_regions


	region_dict = {1:[], 2:[], 3:[], 4:[], 5:[]}
	for region in regions:
		region_dict[region[0][0]].append(region)

	perm_stats = []
	if obs_stat != float('-inf'):
		sign_count = 0
		for perm_i in range(num_perm):
			perm_regions = _get_perm_regions_(region_dict)
			perm_enrichments = calc_enrichment(all_genes, cand_gene_indices, perm_regions)
			r1 = perm_enrichments[0] / float(perm_enrichments[1])
			if r1 == 0.0:
				perm_stat = float('-inf')
			else:
				perm_stat = sp.log(r1 / r2)
			perm_stats.append(perm_stat)
			if perm_stat >= obs_stat:
				sign_count += 1
				if sign_count == early_stop_threshold:
					break
			sys.stdout.write('.')
			sys.stdout.flush()
		sys.stdout.write('\n')
		sys.stdout.flush()
		perm_stats.sort()
		p_val = sign_count / float(perm_i + 1)
	else:
		p_val = 1.0
	print 'Permutation p-value estimate: % f' % p_val
	return p_val, perm_stats


def _test_cand_gene_enrichment_():
	#load result
	res_files = []
	#res_files.append((env.env['tmp_dir'] + 'seeds_pid1_seed_size_2n_emma_none.pvals', 'seed_size_2n EMMAX'))
	#res_files.append((env.env['tmp_dir'] + 'seeds_pid2_seed_size_3n_2x4_emmax_none.pvals', 'seed_size_3n_2x4 EMMAX'))
	#res_files.append((env.env['tmp_dir'] + 'seeds_pid3_seed_size_3n_4x2_emmax_none.pvals', 'seed_size_3n_4x2 EMMAX'))
	#res_files.append((env.env['tmp_dir'] + 'seeds_pid4_seed_size_spss_emmax_none.pvals', 'seed_size_spss EMMAX'))
	#res_files.append((env.env['tmp_dir'] + 'seeds_pid1_seed_size_2n_kw_none.pvals', 'seed_size_2n KW'))
#	res_files.append((env.env['tmp_dir'] + 'seeds_pid2_seed_size_3n_2x4_kw_none.pvals', 'seed_size_3n_2x4 KW'))
#	res_files.append((env.env['tmp_dir'] + 'seeds_pid3_seed_size_3n_4x2_kw_none.pvals', 'seed_size_3n_4x2 KW'))
#	res_files.append((env.env['tmp_dir'] + 'seeds_pid4_seed_size_spss_kw_none.pvals', 'seed_size_spss KW'))
	res_files.append((env.env['results_dir'] + 'donald_duck_pid0_25-DHBAG_emmax_log_trans.pvals', '25-DHBAG'))
	res_files.append((env.env['results_dir'] + 'dilkes_pid1_23-DHBAG_emmax_log_trans.pvals', '23-DHBAG'))
	res_files.append((env.env['results_dir'] + 'dilkes_pid2_25-DHBAP_emmax_log_trans.pvals', '25-DHBAP'))
	res_files.append((env.env['results_dir'] + 'dilkes_pid3_23-DHBAP_emmax_log_trans.pvals', '23-DHBAP'))
	res_files.append((env.env['results_dir'] + 'dilkes_pid4_DHBAG_emmax_log_trans.pvals', 'DHBAG'))
	res_files.append((env.env['results_dir'] + 'dilkes_pid5_DHBAP_emmax_none.pvals', 'DHBAP'))
	res_files.append((env.env['results_dir'] + 'dilkes_pid6_total_emmax_log_trans.pvals', 'total'))
	res_files.append((env.env['results_dir'] + 'dilkes_pid7_25-DHBAg_emmax_sqrt_trans.pvals', '25-DHBA_g'))
	res_files.append((env.env['results_dir'] + 'dilkes_pid8_23-DHBAg_emmax_log_trans.pvals', '23-DHBA_g'))
	res_files.append((env.env['results_dir'] + 'dilkes_pid9_Gpct_emmax_none.pvals', 'Gpct'))
	res_files.append((env.env['results_dir'] + 'dilkes_pid10_25pct_emmax_none.pvals', '25pct'))
	res_files.append((env.env['results_dir'] + 'dilkes_pid12_SM_emmax_none.pvals', 'SM'))
	res_files.append((env.env['results_dir'] + 'dilkes_pid14_unkown1_emmax_log_trans.pvals', 'unknown1'))
#	res_files.append((env.env['results_dir'] + 'donald_duck_pid1_telomere_length_emmax_none.pvals', 'telomere_emmax'))
#	res_files.append((env.env['results_dir'] + 'donald_duck_pid1_telomere_length_emmax_sqrt_trans.pvals', 'telomere_emmax_sqrt_trans'))
#	res_files.append((env.env['results_dir'] + 'donald_duck_pid1_telomere_length_kw_none.pvals', 'telomere_kw'))
	#cgl_file = env.env['phen_dir'] + 'seeds_cgl.csv'
	#cgl_file = env.env['phen_dir'] + 'teleomere_cgl.csv'
	cgl_file = env.env['phen_dir'] + 'Dilkes_data_candidate_gene_list.csv'
	for res_file, name in reversed(res_files):

		r = gr.Result(res_file)
		r.neg_log_trans()
		r.filter_attr('mafs', 10)
		res = r.candidate_gene_enrichments(cgl_file=cgl_file, methods=['snps_perm', 'gene_perm'],
					pval_thresholds=[ 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001],
					num_perm=500, obs_genes_file=env.env['tmp_dir'] + 'enrichments_' + name + '_obs_genes.csv',
					file_prefix=env.env['tmp_dir'] + 'enrichments_' + name, early_stop_threshold=25)



if __name__ == '__main__':
	pass
	#cand_genes_enrichment(5,43)
	#plot_gene_distances()
	#plot_gene_pvalue_count()
	#plot_enrichment_stat([1,2,3,4,5,6,7],145,gene_window_size=20000)
	#plot_enrichment_histogram([1,2,3,4,5,6,7],145)
	#plot_enrichment_ratios([1,2,3,4,5,6,7],145,use_fixed_pvalue_limits=False,gene_window_size=10000,max_num_pvals=160000,num_steps=160)
	_test_cand_gene_enrichment_()


