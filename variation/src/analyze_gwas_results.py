"""
A container for functions which aim to analyze or process gwas results, for some aim.
"""

import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
import warnings
import itertools as it
import env
import random
import phenotypeData as pd



def qq_plots(results, num_dots, max_log_val, file_prefix, method_types=['kw', 'emma'], mapping_labels=None,
	 	phen_name=None, perm_pvalues=None, is_binary=False, **kwargs):
	"""
	Plots both log and normal qq plots.
	"""
	log_pdf_file = None
	log_png_file = None
	pdf_file = None
	png_file = None
	if file_prefix:
		log_pdf_file = file_prefix + "_qq_log.pdf"
		log_png_file = file_prefix + "_qq_log.png"
		pdf_file = file_prefix + "_qq.pdf"
		png_file = file_prefix + "_qq.png"
	qq_plot(results, num_dots, method_types=method_types, mapping_labels=mapping_labels, phenName=phen_name,
		pdfFile=pdf_file, pngFile=png_file, perm_pvalues=perm_pvalues, isBinary=is_binary, kwargs=kwargs)
	log_qq_plot(results, num_dots, max_log_val, method_types=method_types, mapping_labels=mapping_labels,
		phenName=phen_name, pdfFile=log_pdf_file, pngFile=log_png_file, perm_pvalues=perm_pvalues,
		isBinary=is_binary, kwargs=kwargs)


def _getQuantiles_(scores, numQuantiles):
	scores.sort()
	quantiles = []
	for i in range(1, numQuantiles + 1):
		j = int(len(scores) * i / (numQuantiles + 2))
		quantiles.append(scores[j])
	return quantiles


def _calcMedian_(scores, exp_median=0.5):
        scores.sort()
	median = scores[len(scores) / 2]
	return (exp_median - median)

def _estAreaBetweenCurves_(quantiles, expQuantiles):
	area = 0
	for i in range(0, len(quantiles) - 1):
		area += (expQuantiles[i + 1] - expQuantiles[i]) * (abs(quantiles[i + 1] - expQuantiles[i + 1] + quantiles[i] - expQuantiles[i])) / 2.0
	return area

def _calcKS_(scores, exp_scores=None):
	ret = {}
	ret["D"] = -1
	try:
		from rpy import r
		if exp_scores:
			res = r.ks_test(scores, exp_scores)
		else:
			res = r.ks_test(scores, "punif")
		ret = res["statistic"]
		ret["p.value"] = res["p.value"]
	except Exception, message:
		print "Calculating KS failed??", message
	return ret


def _getExpectedPvalueQuantiles_(numQuantiles):
	quantiles = []
	for i in range(1, numQuantiles + 1):
		quantiles.append(float(i) / (numQuantiles + 2))
	return quantiles



def qq_plot(results, numQuantiles, method_types=["kw", "emma"], mapping_labels=None, phenName=None, pdfFile=None, pngFile=None,
	    perm_pvalues=None, **kwargs):
	"""
	QQ-plot for the given results.
	"""

	if not mapping_labels:
		mapping_labels = method_types

	plt.figure(figsize=(5, 4))
	#plt.figure(figsize=(10,8))
	#plt.figure(figsize=(4,3.5))
	plt.axes([0.15, 0.14, 0.82, 0.79])
	plt.plot([0, 1], [0, 1], "k", label="Expected")
	areas = []
	medians = []
	for method_type, label in zip(method_types, mapping_labels):
		result = results[label]
		newScores = result.scores[:]
		quantiles = _getQuantiles_(newScores, numQuantiles)
		if perm_pvalues and method_type in ['kw', 'ft']:
			print "Getting exp. quantiles for permuted p-values"
			expQuantiles = _getQuantiles_(perm_pvalues, numQuantiles)
			q_i = numQuantiles / 2
			if numQuantiles % 2 == 0: #even
				exp_median = (expQuantiles[q_i - 1] + expQuantiles[q_i]) / 2.0
			else: #odd
				exp_median = expQuantiles[q_i]

		else:
			exp_median = 0.5
			expQuantiles = _getExpectedPvalueQuantiles_(numQuantiles)
		area = _estAreaBetweenCurves_(quantiles, expQuantiles)
		median = _calcMedian_(newScores, exp_median)
		plt.plot(expQuantiles, quantiles, label=label + ", A=" + str(round(area, 3)) + ", M=" + str(round(median, 3)))
		areas.append(area)
		medians.append(median)

	if phenName:
		plt.title(phenName)
	fontProp = matplotlib.font_manager.FontProperties(size=8)
	plt.legend(loc=2, numpoints=4, handlelen=0.05, markerscale=1, prop=fontProp, pad=0.018)
	plt.axis([-0.01, 1.01, -0.01, 1.01])
	plt.xlabel("Expected $p$-value")
	plt.ylabel("Observed $p$-value")
	if pdfFile:
		plt.savefig(pdfFile, format="pdf")
	if pngFile:
		plt.savefig(pngFile, format="png", dpi=300)
	elif not pdfFile:
		plt.show()
	plt.clf()
	return (areas, medians)




def _getLogQuantilesMaxVal_(scores, maxScore=None):
	scores.sort()
	i = 0
	new_score = -math.log(scores[i], 10)
	score = new_score
	#print score, maxScore
	while i < len(scores) - 1 and new_score > maxScore:
		score = new_score
		i += 1
		new_score = -math.log(scores[i], 10)

	maxVal = math.log((len(scores)) / float(i + 1), 10)
	#print maxVal,i, score, maxScore
	return(maxVal)

def _getLogQuantiles_(scores, numDots, maxVal=None):
	scores.sort()
	quantiles = []
	for i in range(0, numDots):
		j = max(int(round(math.pow(10, -(float(i) / (numDots - 1)) * maxVal) * len(scores))) - 1, 0) #A bug fixed to make sure j was not less than 0
		val = -math.log10(scores[j])
		quantiles.append(val)
	return quantiles


def _estLogSlope_(ys, xs=None):
	if xs:
		q1 = _getQuantiles_(xs, 1000)
	else:
		q1 = _getExpectedPvalueQuantiles_(1000)
	q2 = _getQuantiles_(ys, 1000)
	b_sum = 0.0
	num_valid = 0
	for (x, y) in zip(q1, q2):
		if x < 1.0:
			b_sum += math.log(y, 10) / math.log(x, 10)
			num_valid += 1
	return(b_sum / num_valid)


def log_qq_plot(results, numDots, maxVal, method_types=['kw', 'emma'], mapping_labels=None, phenName=None, pdfFile=None,
	        pngFile=None, perm_pvalues=None, **kwargs):
	"""
	log-transformed QQ-plot for the given results.
	"""
	def _getExpectedLogQuantiles_():
		quantiles = []
		for i in range(1, numDots + 1):
			quantiles.append((float(i) / (numDots + 2.0)) * maxVal)
		return quantiles

	if not mapping_labels:
		mapping_labels = method_types
	plt.figure(figsize=(5, 4))
	plt.axes([0.15, 0.14, 0.82, 0.79])
	maxVal = min(math.log10(len(results[mapping_labels[0]].scores)), maxVal)
	minVal = (1.0 / numDots) * maxVal
	valRange = maxVal - minVal
	plt.plot([minVal, maxVal], [minVal, maxVal], "k", label="Expected")
	maxObsVals = []
	areas = []
	ds = []
	slopes = []
	for method_type, label in zip(method_types, mapping_labels):
		result = results[label]
		if perm_pvalues and method_type in ['kw', 'ft']:
			exp_maxVal = _getLogQuantilesMaxVal_(perm_pvalues[:], maxVal)
			expQuantiles = _getLogQuantiles_(perm_pvalues[:], numDots, exp_maxVal)
			ks_res = _calcKS_(result.scores, perm_pvalues)
			quantiles = _getLogQuantiles_(result.scores[:], numDots, exp_maxVal)
			slope = _estLogSlope_(result.scores[:], perm_pvalues)
		else:
			quantiles = _getLogQuantiles_(result.scores[:], numDots, maxVal)
			expQuantiles = _getExpectedLogQuantiles_()
			ks_res = _calcKS_(result.scores)
			slope = _estLogSlope_(result.scores[:])

		area = _estAreaBetweenCurves_(quantiles, expQuantiles)
		areas.append(area)
		slopes.append(slope)
		ds.append(ks_res["D"])
		#plt.plot(expQuantiles, quantiles, label = label+", A="+str(round(area,2))+", D="+str(round(ks_res["D"],3))+", S="+str(round(slope,3)))
		plt.plot(expQuantiles, quantiles, label=label + ", D=" + str(round(ks_res["D"], 3)) + ", S=" + str(round(slope, 3)))
		maxObsVals.append(max(quantiles))

	maxObsVal = max(maxObsVals)
	obsValRange = maxObsVal - minVal
	plt.axis([minVal - 0.025 * valRange, maxVal + 0.025 * valRange, minVal - 0.025 * obsValRange, maxObsVal + 0.025 * obsValRange])
	plt.ylabel("Observed $-log_{10}(p$-value$)$")
	plt.xlabel("Expected $-log_{10}(p$-value$)$")
	if phenName:
		plt.title(phenName)
	fontProp = matplotlib.font_manager.FontProperties(size=8)
	plt.legend(loc=2, numpoints=4, handlelen=0.05, markerscale=1, prop=fontProp, pad=0.018)
	if pdfFile:
		plt.savefig(pdfFile, format="pdf")
	if pngFile:
		plt.savefig(pngFile, format="png", dpi=300)
	elif not pdfFile:
		plt.show()
	plt.clf()
	return (ds, areas, slopes)




def identify_interesting_accessions(sd, snps, snp_chromosomes, snp_positions, snp_ecotypes, num_picked=100):
	"""
	Identifies accessions which share rare haplotype combinations of the given SNPs. 
	"""
	import bisect

	if len(snps) > 10:
		warnings.warn('Number of possible haplotypes is greater than 2^10.')
	snps_array = sp.array(snps, dtype='single')
	fs = sp.sum(snps_array, 1) / len(snp_ecotypes)  #Frequencies of 1's

	haplotype_map = {}
	for i in range(2 ** len(snps)):
		hl = map(int, list(bin(i)[2:]))
		l = [0] * (len(snps) - len(hl))
		l.extend(hl)
		h = tuple(l)
		f = 1
		for i, nt in enumerate(h):
			f *= fs[i] if nt == 1.0 else 1 - fs[i]
		haplotype_map[h] = {'f':f, 'c':0, 'et_occurrences':[]}

	haplotypes = zip(*snps) #list of haplotype tuples (hashable)
	for h in haplotypes:
		haplotype_map[h]['c'] += 1

	l = []
	for h in haplotype_map:
		hm = haplotype_map[h]
		l.append((hm['f'], hm['c'], h))
	l.sort()

	#Now locate the interesting SNPs in the snps data
	chr_pos_list = sd.getChrPosList()
	snps_indices = []
	for chr_pos in zip(snp_chromosomes, snp_positions):
		i = bisect.bisect(chr_pos_list, chr_pos) - 1
		if chr_pos_list[i] != chr_pos:
			raise Exception('The SNP at chr=%d, pos=%d, was not found in the snps data.' % chr_pos)
		snps_indices.append(i)
	f_snps = sd.getSnps() #full SNPs
	f_snps = [f_snps[i] for i in snps_indices]
	for et, h in it.izip(sd.accessions, zip(*f_snps)):
		if et in snp_ecotypes: continue
		haplotype_map[h]['et_occurrences'].append(et)

	et_dict = pd.get_ecotype_id_info_dict()

#	print 'expected_frequency, num_phenotyped, num_not_phenotyped, non_phenotyped_ecotypes..'
#	for f, c, h in l:
#		ets = map(int, haplotype_map[h]['et_occurrences'])
#		if len(haplotype_map[h]['et_occurrences']):
#			print '%f, %d, %d, %s' % (f, c, len(haplotype_map[h]['et_occurrences']),
#					','.join(map(str, zip(ets, [et_dict[et][0] for et in ets]))))


	snps = map(list, snps)
	haplotype_list = []
	sd_accessions = sd.accessions
	num_ecotypes = 1
	while len(snp_ecotypes) < len(sd_accessions):
		for i, t in enumerate(l):
			f, c, h = l[i]
			if len(haplotype_map[h]['et_occurrences']):
				break
		else:
			break

		f, c, h = l[i]
		ets = [haplotype_map[h]['et_occurrences'][0]]
		print 'Iteration %d: %f, %d, %d, %s' % (num_ecotypes, f, c, len(haplotype_map[h]['et_occurrences']),
						str((int(ets[0]), et_dict[int(ets[0])])))
		haplotype_list.append((f, c, ets))
		remove_ids = [sd_accessions.index(et) for et in ets]
		for snp, f_snp in zip(snps, f_snps):
			for i, nt in enumerate(f_snp):
				if i in remove_ids:
					snp.append(nt)
		for i in remove_ids:
			snp_ecotypes.append(sd_accessions[i])


		snps_array = sp.array(snps, dtype='single')
		fs = sp.sum(snps_array, 1) / len(snp_ecotypes)  #Frequencies of 1's

		haplotype_map = {}
		for i in range(2 ** len(snps)):
			hl = map(int, list(bin(i)[2:]))
			l = [0] * (len(snps) - len(hl))
			l.extend(hl)
			h = tuple(l)
			f = 1
			for i, nt in enumerate(h):
				f *= fs[i] if nt == 1.0 else 1 - fs[i]
			haplotype_map[h] = {'f':f, 'c':0, 'et_occurrences':[]}

		haplotypes = zip(*snps) #list of haplotype tuples (hashable)
		for h in haplotypes:
			haplotype_map[h]['c'] += 1

		l = []
		for h in haplotype_map:
			hm = haplotype_map[h]
			l.append((hm['f'], hm['c'], h))
		l.sort()
		for et, h in it.izip(sd.accessions, zip(*f_snps)):
			if et in snp_ecotypes: continue
			haplotype_map[h]['et_occurrences'].append(et)

#		print 'expected_frequency, num_phenotyped, num_not_phenotyped, non_phenotyped_ecotypes..'
#		for f, c, h in l:
#			ets = map(int, haplotype_map[h]['et_occurrences'])
#			if len(haplotype_map[h]['et_occurrences']):
#				print '%f, %d, %d, %s' % (f, c, len(haplotype_map[h]['et_occurrences']),
#						','.join(map(str, zip(ets, [et_dict[et][0] for et in ets]))))
		num_ecotypes += 1






def identify_interesting_haplotypes(chrom_pos_list, phenotype_file, pid):
	import dataParsers as dp
	import bisect
	sd = dp.load_250K_snps()
	#sd = dp.load_1001_full_snps()
	phed = pd.parse_phenotype_file(phenotype_file)
	phed.convert_to_averages()
	sd.coordinate_w_phenotype_data(phed, pid)
	cpl = sd.getChrPosList()
	all_snps = sd.getSnps()
	snps = []
	snp_chromosomes = []
	snp_positions = []
	for chrom_pos in chrom_pos_list:
		i = bisect.bisect(cpl, chrom_pos) - 1
		if cpl[i] != chrom_pos:
			raise Exception('SNP not found')
		snps.append(all_snps[i])
		snp_chromosomes.append(chrom_pos[0])
		snp_positions.append(chrom_pos[1])
	sd = dp.load_250K_snps()
	#sd = dp.load_1001_full_snps()
	identify_interesting_accessions(sd, snps, snp_chromosomes, snp_positions, phed.get_ecotypes(pid))



if __name__ == '__main__':
	#chrom_pos_list = [(1, 22349990), (1, 25296405), (4, 453759), (5, 24053984), (5, 25458236)]#KW
	#chrom_pos_list = [(1,5133207),(3,5423868),(4,5795054),(4,10864907),(5,894530)]#EMMAX
	chrom_pos_list = [(1, 5133207), (3, 5423868), (4, 5795054), (4, 10864907), (5, 894530), (1, 22349990),
			(1, 25296405), (4, 453759), (5, 24053984), (5, 25458236)]#Both
	#chrom_pos_list = [(1, 5133217), (3, 20581778), (4, 5458010), (4, 5727758), (5, 1001970), (1, 1427391), (1, 22360158), (1, 29293376), (5, 21020181), (5, 24054819)]
	pid = 1
	identify_interesting_haplotypes(chrom_pos_list, env.env['phen_dir'] + 'telomere_lengths_all.csv', pid)

