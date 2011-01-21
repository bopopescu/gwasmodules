"""
Contains classes to handle results of GWAS

020210
- TODO:
	- Result object should be coupled to a marker data, in a simple manner.  
	  E.g. with a pointer, and/or a list marking what markers are in the result data. 

"""


import pdb, gc
import dbutils
import csv
import math
import itertools as it
try:
	import scipy as sp
	import scipy.stats as st
except Exception, err_str:
	print 'scipy is missing:', err_str
import random
import sys
import env
import tair_converter as tc
#A dictionary for loaded results.. to avoid reloading. 
#Use carefully to avoid memory leaks!
loaded_results = {}



class ResultType(object):
	def __init__(self, resultType=None, fileType=None, datasetName=None, resultDir=None, mafCutoff=0, logTransform=None, name=None):
		self.resultType = resultType
		self.fileType = fileType
		self.datasetName = datasetName
		self.resultDir = resultDir
		self.mafCutoff = mafCutoff
		self.logTransform = logTransform
		if self.logTransform == None:
			self.logTransform = self.fileType == ".pvals"
		self.name = name
		if not name and datasetName:
			self.name = resultType + "_" + datasetName
		if resultDir and datasetName:
			self.filenamePrefix = resultDir + resultType + "_" + datasetName + "_"


	def getFileName(self, phed, phenotypeID, secondRun=False):
		print phenotypeID
		if secondRun:
			filename = self.filenamePrefix + str(phed.getPhenotypeName(phenotypeID)) + ".sr" + self.fileType
		else:
			filename = self.filenamePrefix + str(phed.getPhenotypeName(phenotypeID)) + self.fileType
		return filename

	def __str__(self):
		return self.resultType + "_" + self.datasetName


class Region(object):

	def __init__(self, chromosome, startPos, endPos, snps=None, snpsd_indices=None, snpsInfo=None, name=""):
		self.name = name
		self.chromosome = chromosome
		self.startPos = startPos
		self.endPos = endPos
		self.size = endPos - startPos
		self.snps = []
		if snps:
			self.snps = snps
		self.snpsd_indices = []  #Indicies in other data structures.
		if snpsd_indices:
			self.snpsd_indices = snpsd_indices
		self.snpsInfo = {}       #E.g. a dictionary of information about the snps.  
		if snpsInfo:
			self.snpsInfo = snpsInfo

	def __cmp__(self, other):
		return cmp((self.chromosome, self.startPos), (other.chromosome, other.startPos))

	def __str__(self):
		return "Chr.:" + str(self.chromosome) + ", start pos.:" + str(self.startPos) + ", end pos.:" + str(self.endPos)

	def get_chr_pos_str(self):
		return str(self.chromosome) + "_" + str(self.startPos) + "_" + str(self.endPos)

	def overlapping(self, region):
		if self.chromosome != region.chromosome:
			return False
		else:
			if self.startPos < region.startPos:
				return self.endPos >= region.startPos
			else:
				return region.endPos >= self.startPos

	def merge(self, region):
		if not self.overlapping(region):
			raise Exception
		new_start = min(self.startPos, region.startPos)
		new_stop = max(self.endPos, region.endPos)
		self.startPos = new_start
		self.endPos = new_stop
		self.size = self.endPos - self.startPos


fri_region_small = Region(4, 220000, 320000, name="FRI_small")
flc_region_small = Region(5, 3130000, 3230000, name="FLC_small")

fri_region = Region(4, 120000, 550000, name="FRI_large")
flc_region = Region(5, 2100000, 3300000, name="FLC_large")

fri_flc_regions = [fri_region_small, fri_region, flc_region_small, flc_region]

def getRegions(regionSet, window=[25000, 25000]):
	"""
	Converts a set of crh_pos into a list of regions objects.
	"""

	res_ls = list(regionSet)
	res_ls.sort()

	oldPos = 0
	countRegions = 0
	chr = -1
	regions = []
	curPosList = []
	for i in range(0, len(res_ls)):
		pos = res_ls[i][1]
		if chr != res_ls[i][0]:
			if len(curPosList):
				regions.append(curPosList)
			curPosList = []
			countRegions += 1
			chr = res_ls[i][0]
		elif pos - oldPos > sum(window):
			if len(curPosList):
				regions.append(curPosList)
			curPosList = []
			countRegions += 1

		curPosList.append((chr, pos))
		oldPos = pos

	print countRegions, len(regions)

	regionList = []
	for region in regions:
		chr = region[0][0]
		positions = []
		for (chr, pos) in region:
			positions.append(pos)
		regionList.append(Region(chr, min(positions) - window[0], max(positions) + window[1]))
	regionList.sort()
	return regionList

# Old result class 01/18/11
#
#class Result(object):
#	"""
#	Contains information on the result.  (The old object is renamed Result_old, for now..)  Poised to cause problems?
#	"""
#	def __init__(self, result_file=None, scores=None, snps_data=None, accessions=None, name=None,
#		     result_type=None, phen_id=None, positions=None, chromosomes=None, marfs=None,
#		     mafs=None, snps=None, **snp_results_info):
#		"""
#		032610: A new init function to clean up the previous mess...  
#
#		A simple init function.
#		
#		snps_data is assumed to match the results, (in term of positions, etc)
#		
#		accessions (when used) should fit the results!
#		"""
#
#		#Contain various information for every snp and position
#		#FIXME: To replace older lists.. such as positions, chromosomes, scores, etc. 
#		self.snp_results = {}
#
#		self.phen_id = phen_id
#		self.result_type = result_type
#		self.name = name
#		self.positions = positions
#		self.chromosomes = chromosomes
#		self.marfs = marfs #Minor allele relative frequencies.
#		self.mafs = mafs
#		self.snps = snps
#		self.accessions = accessions
#		if scores != None:
#			self.scores = list(scores) #Scores or p-values
#		else:
#			self.scores = []
#		if not self.positions:
#			self.positions = []
#		if not self.chromosomes:
#			self.chromosomes = []
#		if not self.marfs:
#			self.marfs = []
#		if not self.mafs:
#			self.mafs = []
#		if not self.snps:
#			self.snps = []
#		if not self.accessions:
#			self.accessions = []
#
#		self.orders = None
#		self.ranks = None
#		self.chromosome_ends = []
#
#		if result_file:
#			self._load_result_(result_file)
#		elif snps_data:
#			self._load_snps_data_(snps_data)
#
#		self.snp_results['chromosomes']	 = self.chromosomes
#		self.snp_results['scores'] = self.scores
#		self.snp_results['positions'] = self.positions
#		self.snp_results['marfs'] = self.marfs
#		self.snp_results['mafs'] = self.mafs
#		self.snp_results['snps'] = self.snps
#		if snp_results_info:
#			for info in snp_results_info:
#				self.snp_results[info] = snp_results_info[info]
#
#
#
#	def _load_result_(self, resultFile):
#		f = open(resultFile, "r")
#		lines = f.readlines()
#		print "Loading", len(lines), "lines of results."
#		try_delims = [",", "\t", " ", ]
#		i = 0
#		delim = try_delims[i]
#		while len(lines[0].split(delim)) == 1:
#			i += 1
#			delim = try_delims[i]
#		if len(lines[0].split(delim)) == 1:
#			raise Exception("Appropriate delimiter wasn't found.")
#		#print "Delimiter used:",str(delim)	
#		line = lines[0].split(delim)
#		withMAF = len(line) > 3
#		#print "withMAF:", withMAF
#		start = 0
#		currChrom = -1
#		lastPos = 0
#		if not withMAF:
#			for i in range(start, len(lines)):
#				line = lines[i].split(delim)
#				newChrom = int(line[0].strip())
#				if newChrom != currChrom:
#					currChrom = newChrom
#					if lastPos:
#						self.chromosome_ends.append(lastPos)
#				self.chromosomes.append(newChrom)
#				lastPos = int(line[1].strip())
#				self.positions.append(lastPos)
#				self.scores.append(float(line[2].strip()))
#		else:
#			extra_columns = []
#			if not (line[0].strip()).isdigit():
#				del lines[0]
#				#print "Header detected"
#				columns = map(str.strip, line)
#				#print columns
#				extra_columns = columns[5:]
#			self.snp_results = {}
#			for col in extra_columns:
#				self.snp_results[col] = []
#			for line in lines:
#				line_list = map(str.strip, line.split(delim))
#				newChrom = int(line_list[0])
#				if newChrom != currChrom:
#					currChrom = newChrom
#					if lastPos:
#						self.chromosome_ends.append(lastPos)
#				self.chromosomes.append(newChrom)
#				lastPos = int(line_list[1])
#				self.positions.append(lastPos)
#				self.scores.append(float(line_list[2]))
#				self.marfs.append(float(line_list[3]))
#				self.mafs.append(float(line_list[4]))
#				for i, col in enumerate(extra_columns):
#					self.snp_results[col].append(line_list[5 + i])
#
#			assert len(self.positions) == len(self.mafs) == len(self.scores), \
#			     "length of scores, positions, mafs, isn't equal"
#
#		self.chromosome_ends.append(lastPos)
#		self._rank_scores_()
#
#
#	def _load_snps_data_(self, snps_data):
#		self.chromosome_ends = []
#		for i, snpsd in enumerate(snps_data.snpsDataList):
#			self.chromosomes.extend([snps_data.chromosomes[i]] * len(snpsd.positions))
#			self.positions.extend(snpsd.positions)
#			self.snps.extend(snpsd.snps)
#			self.chromosome_ends.append(snpsd.positions[-1])
#		if snps_data.data_format == 'float':
#			self.mafs = [1] * len(self.snps)
#			self.marfs = [1] * len(self.snps)
#		else:
#			maf_d = snps_data.get_mafs()
#			self.mafs = maf_d['mafs']
#			self.marfs = maf_d['marfs']
#
#
#	def _rank_scores_(self):
#		"""
#		Generates two data structures:
#		self.orders (SNP indices ordered by rank)
#		self.ranks (ranks of the SNPs)
#		"""
#		rank_ls = zip(self.scores, range(len(self.scores)))
#		rank_ls.sort(reverse=True)
#		self.orders = []
#		for j in range(len(rank_ls)):
#			(s, i) = rank_ls[j]
#			self.orders.append(i)
#			rank_ls[j] = (i, j)
#
#		rank_ls.sort()
#		self.ranks = []
#		for (i, j) in rank_ls:
#			self.ranks.append(j + 1)
#
#
#	def _sort_by_chr_pos_(self):
#		if self.snps:
#			res_ls = zip(self.chromosomes, self.positions, self.scores, self.mafs, self.marfs, self.snps)
#		elif self.mafs:
#			res_ls = zip(self.chromosomes, self.positions, self.scores, self.mafs, self.marfs)
#		else:
#			res_ls = zip(self.chromosomes, self.positions, self.scores)
#		res_ls.sort()
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		new_snps = []
#		for res in res_ls:
#			newScores.append(res[2])
#			newPositions.append(res[1])
#			newChromosomes.append(res[0])
#			if self.mafs:
#				newMafs.append(res[3])
#				newMarfs.append(res[4])
#			if self.snps:
#				new_snps.append(res[5])
#
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#		self.snps = new_snps
#
#
#
#	def candidate_gene_enrichments(self, cgl=None, cgl_file=None, pval_thresholds=[0.01], gene_radius=20000,
#				methods=['chi_square'], num_perm=500, file_prefix=None,
#				obs_genes_file=None, early_stop_threshold=25, all_genes=None, cand_gene_indices=None):
#		"""
#		Performs CGR analysis on this results object.
#		
#		cgl is a list of genes.
#		cgl_file is a file with a list of genes.
#		
#		method: 
#			chi_square - statistics
#			multinomial - statistics
#			gene_perm - permute genes.
#			snps_perm - permute SNPs
#		"""
#		import pylab
#		import analyze_gene_enrichment as genr
#
#		if not 'chi_square' in methods:
#			methods.append('chi_square')
#
#		chrom_ends = self.get_chromosome_ends()
#		print 'chrom_ends', chrom_ends
#
#		#Parse cgl file
#		if cgl_file:
#			cgl, cg_tair_ids = load_cand_genes_file(cgl_file)
#		for cg in cgl:
#			print str(cg)
#
#
#		if not all_genes:
#			#Load genes from DB.
#			print 'Fetching all genes'
#			all_genes = get_gene_list(include_intron_exons=False, verbose=False)
#			print 'Fetched %d genes.' % len(all_genes)
#			num_genes = len(all_genes)
#
#		if not cand_gene_indices:
#			#Pre-process cgl.
#			cand_gene_indices = []
#			for i, g in enumerate(all_genes):
#				for cg in cgl:
#					if g.dbRef == cg.dbRef:
#						#print g.dbRef, cg.dbRef
#						cand_gene_indices.append(i)
#						break
#			num_cand_genes = len(cand_gene_indices)
#			#print num_cand_genes, cand_gene_indices
#
#
#		method_res_dict = {}
#		for m in methods:
#			method_res_dict[m] = {'statistics':[], 'pvals':[]}
#
#		if obs_genes_file:
#			obs_gene_str = ''
#
#		pval_thresholds.sort()
#		last_thres = 1.0
#		log_enrichments = []
#		for pval_threshold in reversed(pval_thresholds):
#			print 'Using p-value threshold % f' % pval_threshold
#			thres = pval_threshold / last_thres
#			print 'Using corrected threshold % f' % thres
#			last_thres = pval_threshold
#
#			#Filter pvalue file
#			self.filter_percentile(1 - thres)
#
#			#pre-process pvalues
#			regions = self.get_regions(gene_radius=gene_radius)
#
#			#Calculate observed candidate gene enrichment. 
#			obs_enrichments = genr.calc_enrichment(all_genes, cand_gene_indices, regions)
#			r1 = obs_enrichments[0] / float(obs_enrichments[1])
#			r2 = (num_cand_genes / float(num_genes))
#			obs_stat = sp.log(r1 / r2)
#			print 'Observed statistics % f' % obs_stat
#			log_enrichments.append(obs_stat)
#
#			#What cand. genes overlap with regions?
#			obs_cg_indices = obs_enrichments[2]
#			if obs_cg_indices:
#				obs_gene_str += str(pval_threshold) + ','
#				tair_ids = [all_genes[cgi].tairID for cgi in obs_cg_indices]
#				obs_gene_str += ','.join(tair_ids)
#				obs_gene_str += '\n'
#
#
#
#			for method in methods:
#				if method == 'chi_square':
#					chi_sq_pval, chi_sq_stat = genr.get_chi_square_pval(obs_enrichments[0],
#											obs_enrichments[1],
#											num_cand_genes, num_genes)
#					method_res_dict[method]['statistics'].append(chi_sq_stat)
#					method_res_dict[method]['pvals'].append(chi_sq_pval)
#
#				if method == 'multinomial':
#					pass
#				if method == 'gene_perm':
#					p_val, perm_stats = genr.get_gene_perm_pval(obs_stat, regions, all_genes,
#									cand_gene_indices, num_perm=num_perm,
#									early_stop_threshold=early_stop_threshold)
#					method_res_dict[method]['statistics'].append(perm_stats)
#					method_res_dict[method]['pvals'].append(p_val)
#
#				if method == 'snps_perm':
#					p_val, perm_stats = genr.get_snps_perm_pval(obs_stat, regions, all_genes,
#									cand_gene_indices, chrom_ends,
#									num_perm=num_perm,
#									early_stop_threshold=early_stop_threshold)
#
#					method_res_dict[method]['statistics'].append(perm_stats)
#					method_res_dict[method]['pvals'].append(p_val)
#
##						h_res = pylab.hist(perm_stats)
##						pylab.vlines(obs_stat, 0, max(h_res[0]), colors='r')
##						pylab.savefig(env.env['tmp_dir'] + 'test.pdf', format='pdf')
#
#
#		if obs_genes_file:
#			with open(obs_genes_file, 'w') as f:
#				f.write(obs_gene_str)
#
#		#Now the plotting of the results.
#		method_name_dict = {'chi_square':'Chi-square test', 'gene_perm':'Candidate gene permutations',
#					'snps_perm':'SNP positions permutation (chromosome rotation)' }
#
#
#		pval_thresholds.reverse()
#		pos_list = range(len(pval_thresholds))
#		for m in methods:
#			pylab.figure()
#			pvals = []
#			for p in method_res_dict[m]['pvals']:
#				if p != 0:
#					pvals.append(p)
#				else:
#					pvals.append(0.5 / num_perm)
#			neg_log_pvals = map(lambda x:-math.log10(x), pvals)
#			pylab.barh(pos_list, neg_log_pvals, align='center', color='g', alpha=0.6)
#			pylab.ylim((-1, len(pval_thresholds)))
#			pylab.yticks(pos_list, map(str, pval_thresholds))
#			pylab.xlabel('Enrichment -log(p-value)')
#			pylab.ylabel('p-value percentile threshold')
#			ymin, ymax = pylab.ylim()
#			xmin, xmax = pylab.xlim()
#			pylab.axvline(-math.log10(0.05), ymin=ymin, ymax=ymax, color='r')
#			pylab.xlim((0, max(-math.log10(0.05), max(neg_log_pvals)) * 1.05))
#			pylab.title(method_name_dict[m])
#			pylab.savefig(file_prefix + '_' + m + '.png', format='png')
#
#
#		pylab.figure()
#		pylab.barh(pos_list, log_enrichments, align='center', color='g', alpha=0.6)
#		pylab.ylim((-1, len(pval_thresholds)))
#		pylab.yticks(pos_list, map(str, pval_thresholds))
#		pylab.xlabel('Enrichment (log[ratio])')
#		pylab.ylabel('p-value percentile threshold')
#		pylab.savefig(file_prefix + '_enrichment_ratio.png', format='png')
#
#		return {'enr_stats':log_enrichments, 'method_res_dict':method_res_dict}
#
#
#
#
#	def _plot_small_manhattan_(self, pdf_file=None, png_file=None, min_score=None, max_score=None,
#				type="pvals", ylab="$-$log$_{10}(p-$value$)$", plot_bonferroni=False,
#				cand_genes=None, threshold=0, highlight_markers=None, chromosome=None):
#		import matplotlib
#		#matplotlib.use('Agg')
#		import matplotlib.pyplot as plt
#
#		displayed_unit = 1000.0 #kbs
#		scoreRange = max_score - min_score
#		plt.figure(figsize=(6, 3.5))
#		plt.axes([0.045, 0.12, 0.95, 0.7])
#		starPoints = [[], [], []]
#		color = 'b'
#		color_map = {1:'b', 2:'g', 3:'r', 4:'c', 5:'m'}
#		if chromosome:
#			color = color_map[chromosome]
#
#		chrom = self.snp_results['chromosomes'][0]
#		positions = map(lambda x: x / displayed_unit, self.snp_results['positions'])
#		scores = self.snp_results['scores']
#		for s_i, (score, pos) in enumerate(it.izip(scores, positions)):
#			if score > max_score:
#				starPoints[0].append(pos)
#				starPoints[1].append(max_score)
#				starPoints[2].append(score)
#				score = max_score
#			scores[s_i] = score
#
#		plt.plot(positions, scores, ".", markersize=3, alpha=0.7, color=color)
#
#		if cand_genes:
#			for cg in cand_genes:
#				plt.axvspan(cg.startPos / displayed_unit, cg.endPos / displayed_unit, facecolor='k', alpha=0.5, linewidth=0)
#
#
#		if highlight_markers:
#			ys = []
#			xs = []
#			for c, p, score in highlight_markers:
#				xs.append(p / displayed_unit)
#				if score > max_score:
#					plt.text(x, max_score * 1.1, str(round(score, 2)), rotation=45, size="small")
#					ys.append(max_score)
#				else:
#					ys.append(score)
#			plt.plot(xs, ys, ".", color="#ff9944", markersize=6, alpha=0.7)
#
#		if len(starPoints[0]) > 0:
#			plt.plot(starPoints[0], starPoints[1], ".", color="#ee9922", markersize=4)
#			i = 0
#			while i < len(starPoints[0]):
#				max_point = i
#				cur_pos = starPoints[0][i]
#				while i < len(starPoints[0]) and abs(starPoints[0][i] - cur_pos) < 1000000:
#					if starPoints[2][i] > starPoints[2][max_point]:
#						max_point = i
#					i += 1
#				plt.text(starPoints[0][max_point] - 200000, (starPoints[1][max_point] - 1) * 1.15, str(round(starPoints[2][max_point], 2)), rotation=45, size="small")
#
#
#
#
#		if plot_bonferroni:
#			b_threshold = -math.log10(1.0 / (len(scores) * 20.0))
#			if threshold :
#				plt.plot([0, max(positions)], [b_threshold, b_threshold], ":")
#				threshold = -math.log10(threshold)
#				plt.plot([0, max(positions)], [threshold, threshold], color='#6495ed', linestyle='-.')
#			#Bonferroni threshold
#			else:
#				plt.plot([0, max(positions)], [b_threshold, b_threshold], color='#000000', linestyle="-.")
#
#		plt.axis([min(positions), max(positions), min_score - 0.05 * scoreRange, max_score + 0.05 * scoreRange])
#		if not ylab:
#			if type == "pvals":
#				plt.ylabel('$ - log(p - $value$)$')
#
#			else:
#				plt.ylabel('score')
#		else:
#			plt.ylabel(ylab)
#		plt.xlabel("kilobases")
#		plt.title('Chromsome %d' % chrom)
#
#		if pdf_file:
#			plt.savefig(pdf_file, format="pdf")
#		if png_file:
#			plt.savefig(png_file, format="png", dpi=300, bbox_inches='tight')
#		if not (pdf_file or png_file):
#			plt.show()
#
#		plt.clf()
#		plt.close()
#
#
#	def plot_manhattan(self, pdf_file=None, png_file=None, min_score=None, max_score=None,
#		       percentile=90, type="pvals", ylab="$-$log$_{10}(p-$value$)$",
#		       plot_bonferroni=False, cand_genes=None, threshold=0, highlight_markers=None):
#
#		"""
#		Plots a 'Manhattan' style GWAs plot.
#		"""
#
#		import matplotlib
#		#matplotlib.use('Agg')
#		import matplotlib.pyplot as plt
#
#		"Plotting a Manhattan-style plot with %i markers." % len(self.scores)
#
#		num_scores = len(self.scores)
#		chromosome_ends = self.get_chromosome_ends()
#		result = self.simple_clone()
#
#		chrom_set = set(result.snp_results['chromosomes'])
#		chromosomes = list(chrom_set)
#		chromosomes.sort()
#		if len(chrom_set) == 1:
#			percentile = 0.0
#		if percentile != 0.0:
#			result.filter_percentile(percentile / 100.0)
#
#		if highlight_markers:
#			new_h_markers = []
#			for c, p, pval in highlight_markers:
#				new_h_markers.append((c, p, -math.log10(pval)))
#			highlight_markers = new_h_markers
#
#		if not max_score:
#			max_score = max(result.scores)
#			if highlight_markers:
#				h_scores = [s for c, p, s in highlight_markers]
#				max_score = max(max_score, max(h_scores))
#		if not min_score:
#			if type == "pvals":
#				min_score = 0
#			else:
#				min_score = min(result.scores)
#
#
#		if cand_genes:
#			#processing candidate genes by chromosome
#			chr_cand_genes = {}
#			for chrom in chromosomes:
#				chr_cand_genes[chrom] = []
#			for cg in cand_genes:
#				chr_cand_genes[cg.chromosome].append(cg)
#
#		if len(chrom_set) == 1:
#			chrom = chrom_set.pop()
#			if cand_genes:
#				cand_genes = chr_cand_genes[chrom]
#			return result._plot_small_manhattan_(pdf_file=pdf_file, png_file=png_file, min_score=min_score,
#						max_score=max_score, ylab=ylab, plot_bonferroni=plot_bonferroni,
#						cand_genes=cand_genes, threshold=threshold,
#						highlight_markers=highlight_markers, chromosome=chrom)
#
#
#		scoreRange = max_score - min_score
#		offset = 0
#		chromosomeSplits = result.get_chromosome_splits()
#
#		ticksList1 = []
#		ticksList2 = []
#		textPos = []
#		plt.figure(figsize=(12, 2.8))
#		plt.axes([0.045, 0.15, 0.95, 0.71])
#		starPoints = [[], [], []]
#		chr_offsets = []
#		for i, chromosome_end in enumerate(chromosome_ends):
#			chr_offsets.append(offset)
#			index1 = chromosomeSplits[i][0]
#			index2 = chromosomeSplits[i + 1][0]
#			scoreList = result.scores[index1:index2]
#			posList = result.positions[index1:index2]
#			chrom = chromosomes[i]
#
#			if cand_genes:
#				for cg in chr_cand_genes[chrom]:
#					plt.axvspan(offset + cg.startPos, offset + cg.endPos, facecolor='k', alpha=0.5)
#
#			newPosList = [offset + pos for pos in posList]
#
#			for s_i, (score, pos) in enumerate(it.izip(scoreList, newPosList)):
#				if score > max_score:
#					starPoints[0].append(pos)
#					starPoints[1].append(max_score)
#					starPoints[2].append(score)
#					score = max_score
#				scoreList[s_i] = score
#
#			plt.plot(newPosList, scoreList, ".", markersize=2, alpha=0.7)
#			oldOffset = offset
#			textPos.append(offset + chromosome_end / 2 - 2000000)
#			offset += chromosome_end
#			for j in range(oldOffset, offset, 2000000):
#				ticksList1.append(j)
#			#pdb.set_trace()
#			for j in range(0, chromosome_end, 2000000):
#				if j % 4000000 == 0 and j < chromosome_end - 2000000 :
#					ticksList2.append(j / 1000000)
#				else:
#					ticksList2.append("")
#			#pdb.set_trace()
#
#
#
#		plt.plot(starPoints[0], starPoints[1], ".", color="#ee9922", markersize=4)
#		if len(starPoints[0]) > 0:
#			i = 0
#			while i < len(starPoints[0]):
#				max_point = i
#				cur_pos = starPoints[0][i]
#				while i < len(starPoints[0]) and abs(starPoints[0][i] - cur_pos) < 3000000:
#					if starPoints[2][i] > starPoints[2][max_point]:
#						max_point = i
#					i += 1
#				plt.text(starPoints[0][max_point] - 1000000, (starPoints[1][max_point] - 1) * 1.15, str(round(starPoints[2][max_point], 2)), rotation=45, size="small")
#
#
#		if highlight_markers:
#			ys = []
#			xs = []
#			for c, p, score in highlight_markers:
#				x = chr_offsets[c - 1] + p
#				xs.append(x)
#				if score > max_score:
#					plt.text(x, max_score * 1.1, str(round(score, 2)), rotation=45, size="small")
#					ys.append(max_score)
#				else:
#					ys.append(score)
#			plt.plot(xs, ys, ".", color="#ff9944", markersize=6, alpha=0.8)
#
#
#		if plot_bonferroni:
#			b_threshold = -math.log10(1.0 / (num_scores * 20.0))
#			if threshold :
#				plt.plot([0, sum(result.chromosome_ends)], [b_threshold, b_threshold], ":")
#				threshold = -math.log10(threshold)
#				plt.plot([0, sum(result.chromosome_ends)], [threshold, threshold], color='#6495ed', linestyle='-.')
#			#Bonferroni threshold
#			else:
#				plt.plot([0, sum(result.chromosome_ends)], [b_threshold, b_threshold], color='#000000', linestyle="-.")
#
#		plt.axis([0, sum(result.chromosome_ends), min_score - 0.05 * scoreRange, max_score + 0.05 * scoreRange])
#		plt.xticks(ticksList1, ticksList2, fontsize='x-small')
#		#pdb.set_trace()
#		if not ylab:
#			if type == "pvals":
#				plt.ylabel('$ - log(p - $value$)$')
#
#			else:
#				plt.ylabel('score')
#		else:
#			plt.ylabel(ylab)
#		plt.xlabel("Mb")
#
#		if pdf_file:
#			plt.savefig(pdf_file, format="pdf")
#		if png_file:
#			plt.savefig(png_file, format="png", dpi=300, bbox_inches='tight')
#		if not (pdf_file or png_file):
#			plt.show()
#
#		plt.clf()
#		plt.close()
#
#
#
#	def get_chromosome_splits(self):
#		"""
#		Returns list of indices (and prev chromosome), for the when the chromosomes
#		change in the scores, positions indices.
#		USE WITH CARE
#		"""
#		oldChrom = 0
#		chromosome_splits = []
#		for i in range(0, len(self.scores)):
#			newChrom = self.chromosomes[i]
#			if oldChrom != newChrom:
#				while oldChrom < newChrom:
#					oldChrom += 1
#					chromosome_splits.append((i, oldChrom))
#		chromosome_splits.append((i, -1))
#		return chromosome_splits
#
#
#	def neg_log_trans(self):
#		"""
#		apply - log(x) to the pvalues (scores)
#		"""
#		import math, warnings
#		for i, score in enumerate(self.scores):
#			if score != 0.0:
#				self.scores[i] = -math.log(score, 10)
#			else:
#				warnings.warn("P-value is 0, log transformation is invalid.  Score is arbitrary set to 50.")
#				self.scores[i] = 50
#
#
#
#	def filter_percentile(self, percentile):
#		"""
#		Filter top percentile.
#		"""
#		new_scores = []
#		for score in self.scores:
#			new_scores.append(score)
#		new_scores.sort()
#		score_cutoff = new_scores[int(len(new_scores) * percentile)]
#		self.filter_attr("scores", score_cutoff)
#
#
#
#	def filter_common_top_scores(self, result, top_fraction=0.1):
#		import bisect
#		self._rank_scores_()
#		result._rank_scores_()
#		keep_indices_1 = set()
#		keep_indices_2 = set()
#		chr_pos_list_1 = self.get_chr_pos_list()
#		chr_pos_list_2 = result.get_chr_pos_list()
#		for i in self.orders[:int(top_fraction * len(self.orders))]:
#			j = bisect.bisect(chr_pos_list_2, chr_pos_list_1[i]) - 1
#			if chr_pos_list_2[j] == chr_pos_list_1[i]:
#				keep_indices_1.add(i)
#
#		for i in result.orders[:int(top_fraction * len(result.orders))]:
#			j = bisect.bisect(chr_pos_list_1, chr_pos_list_2[i]) - 1
#			if chr_pos_list_1[j] == chr_pos_list_2[i]:
#				keep_indices_2.add(i)
#
#		keep_indices_list_1 = list(keep_indices_1)
#		for i in keep_indices_2:
#			keep_indices_1.add(bisect.bisect(chr_pos_list_1, chr_pos_list_2[i]) - 1)
#
#		for i in keep_indices_list_1:
#			keep_indices_2.add(bisect.bisect(chr_pos_list_2, chr_pos_list_1[i]) - 1)
#
#		print len(keep_indices_1), len(keep_indices_2)
#		self.filter_indices(list(keep_indices_1))
#		result.filter_indices(list(keep_indices_2))
#		return keep_indices_1, keep_indices_2
#
#
#	def filter_indices(self, indices_to_keep):
#		indices_to_keep.sort()
#		new_snp_results = {}
#		for info in self.snp_results:
#			new_snp_results[info] = []
#		count = len(self.scores)
#		for i in indices_to_keep:
#			for info in self.snp_results:
#				if self.snp_results[info]:
#					new_snp_results[info].append(self.snp_results[info][i])
#
#		self.snp_results = new_snp_results
#		self.scores = self.snp_results['scores']
#		self.positions = self.snp_results['positions']
#		self.chromosomes = self.snp_results['chromosomes']
#		self.mafs = self.snp_results['mafs']
#		self.marfs = self.snp_results['marfs']
#		self.snps = self.snp_results['snps']
#		print "%i scores were removed." % (count - len(self.scores))
#
#
#	def filter_attr(self, attr_name, attr_threshold, reversed=False):
#		"""
#		Filter out scores / pvalues etc. which have attr < attr_threshold.
#
#		attr are e.g.
#		'mafs', 'marfs', 'scores', etc.
#		"""
#		print "Filtering for attribute ' % s' with threshold: %g" % (attr_name, attr_threshold)
#		attr = getattr(self, attr_name)
#		new_snp_results = {}
#		for info in self.snp_results:
#			new_snp_results[info] = []
#		count = len(self.scores)
#		for i in range(len(self.scores)):
#			if reversed:
#				if attr[i] <= attr_threshold:
#					for info in self.snp_results:
#						if len(self.snp_results[info]) > 0:
#							new_snp_results[info].append(self.snp_results[info][i])
#			else:
#				if attr[i] >= attr_threshold:
#					for info in self.snp_results:
#						if len(self.snp_results[info]) > 0:
#							new_snp_results[info].append(self.snp_results[info][i])
#
#		self.snp_results = new_snp_results
#		self.scores = self.snp_results['scores']
#		self.positions = self.snp_results['positions']
#		self.chromosomes = self.snp_results['chromosomes']
#		if self.mafs:
#			self.mafs = self.snp_results['mafs']
#			self.marfs = self.snp_results['marfs']
#		if self.snps:
#			self.snps = self.snp_results['snps']
#		print "%i scores were removed out of %i." % (count - len(self.scores), count)
#		return len(self.scores)
#
#
#
#	def filter_non_segregating_snps(self, ecotype1, ecotype2, accessions=None):
#		"""
#		Filter out all SNPs which are not segregating in the two accessions.
#
#		Assumes the accessions map the results objects SNPs. (and that they are defined)
#		"""
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		new_snps = []
#
#		if accessions:
#			ecotypes = accessions
#		else:
#			ecotypes = self.accessions
#
#		e_i1 = ecotypes.index(ecotype1)
#		e_i2 = ecotypes.index(ecotype2)
#
#		for i in range(len(self.snps)):
#			snp = self.snps[i]
#			if snp[e_i1] != snp[e_i2]:
#				newScores.append(self.scores[i])
#				newPositions.append(self.positions[i])
#				newChromosomes.append(self.chromosomes[i])
#				newMafs.append(self.mafs[i])
#				newMarfs.append(self.marfs[i])
#				new_snps.append(snp)
#
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#		self.snps = new_snps
#
#
#
##	def filter_score_cutoff(self, score_cutoff):
##
##		newScores = []
##		newPositions = []
##		newChromosomes = []
##		newMafs = []
##		newMarfs = []
##		new_snps = []
##
##		for i in range(0,len(self.scores)):
##			if self.scores[i]>score_cutoff:
##				newScores.append(self.scores[i])
##				newPositions.append(self.positions[i])
##				newChromosomes.append(self.chromosomes[i])
##				newMafs.append(self.mafs[i])
##				newMarfs.append(self.marfs[i])
##				if self.snps:
##					new_snps.append(self.snps[i])
##
##		self.scores = newScores
##		self.positions = newPositions
##		self.chromosomes = newChromosomes
##		self.mafs = newMafs
##		self.marfs = newMarfs
##		self.snps = new_snps
#
#
#
#	def get_top_snps_result(self, n):
#		"""
#		returns top n SNPs
#		"""
#		import copy
#		result = copy.deepcopy(self) #Cloning
#		result.filter_top_snps(n)
#		return result
#
#	def get_top_snps(self, n):
#		"""
#		returns top n SNPs
#		"""
#		self._rank_scores_() #Making sure the ranks are updated
#		chr_pos_list = []
#		for i in self.orders[:n]:
#			chr_pos_list.append((self.chromosomes[i], self.positions[i]))
#		return chr_pos_list
#
#
#	def get_top_genes(self, n, window_size=5000, conn=None):
#		"""
#		Returns a set of (chromosome, start_pos, end_pos), for genes found.
#		"""
#		self._rank_scores_() #Making sure the ranks are updated
#		genes = set()
#		snp_ix = []
#		i = 0
#		while len(genes) < n:
#			snp_i = self.orders[i]
#			p = self.positions[snp_i]
#			c = self.chromosomes[snp_i]
#			c_genes = get_gene_list(p - window_size, p + window_size, c, include_intron_exons=False, conn=conn)
#			for g in c_genes:
#				genes.add((g.chromosome, g.startPos, g.endPos, g.tairID))
#			#print len(genes)
#			i += 1
#		print 'found % d genes' % len(genes)
#		return genes
#
#
#	def get_chr_pos_score_list(self):
#		return zip(self.chromosomes, self.positions, self.scores)
#
#
#	def get_chr_pos_list(self):
#		return zip(self.chromosomes, self.positions)
#
#
#	def get_top_regions(self, n, distance_threshold=25000):
#		"""
#		Returns a list of regions, defined by (chromosome, start_pos, end_pos).
#		"""
#		self._rank_scores_()
#		chromosome_ends = self.get_chromosome_ends()
#
#		i = 0 #how many SNPs are needed
#		region_count = 0
#		regions = [[], [], [], [], []]  #a list of (chromosome,start,stop)
#		while sum(map(len, regions)) < n:
#			snp_i = self.orders[i]
#			snp_pos = self.positions[snp_i]
#			chromosome = self.chromosomes[snp_i]
#			new_snp_reg = (chromosome, max(snp_pos - distance_threshold, 0), min(snp_pos + distance_threshold, chromosome_ends[chromosome - 1]))
#			covered = False
#			for reg_i, reg in enumerate(regions[chromosome - 1]):
#				if (new_snp_reg[1] < reg[2] and new_snp_reg[2] > reg[2]) or (new_snp_reg[1] < reg[1] and new_snp_reg[2] > reg[1]): #then merge					
#					regions[chromosome - 1].pop(reg_i) #Removing the region which we're merging with.
#					new_snp_reg = (reg[0], min(reg[1], new_snp_reg[1]), max(reg[2], new_snp_reg[2]))
#				elif (new_snp_reg[1] >= reg[1] and new_snp_reg[2] <= reg[2]):
#					covered = True  #It's already covered
#					break
#			if not covered:
#				regions[chromosome - 1].append(new_snp_reg)
#			i += 1
#		print "It took %i SNPS to find %i regions." % (i, n)
#		regions_list = regions[0] + regions[1] + regions[2] + regions[3] + regions[4]
#		return regions_list
#
#
#	def get_regions(self, gene_radius=20000):
#		"""
#		Returns a list of [[chr,start],[chr,end]] elements, 
#		"""
#		regions = []
#		num_scores = len(self.snp_results['scores'])
#		iter = enumerate(it.izip(self.snp_results['chromosomes'], self.snp_results['positions']))
#		th = sp.array([0, gene_radius])
#		(pos_i, t) = iter.next()
#		chr_pos = sp.array(t)
#		while pos_i < num_scores - 1:
#			start_chr_pos = chr_pos
#			end_chr_pos = start_chr_pos
#			while pos_i < num_scores - 1 and sp.all((chr_pos - end_chr_pos) <= th):
#				end_chr_pos = chr_pos
#				(pos_i, t) = iter.next()
#				chr_pos = sp.array(t)
#			if tuple(end_chr_pos) < tuple(start_chr_pos):
#				pdb.set_trace()
#			if pos_i == num_scores - 1: #Last one..
#				if sp.all((chr_pos - end_chr_pos) <= th):
#					regions.append([start_chr_pos - th, chr_pos + th])
#				else:
#					regions.append([start_chr_pos - th, end_chr_pos + th])
#					regions.append([chr_pos - th, chr_pos + th])
#
#			else:
#				regions.append([start_chr_pos - th, end_chr_pos + th])
#		start_chr_pos = chr_pos
#		end_chr_pos = start_chr_pos
#		regions.append([start_chr_pos - th, end_chr_pos + th])
#		print 'Found % d regions' % len(regions)
#		return regions
#
#
#	def get_region_result(self, chromosome, start_pos, end_pos, buffer=0):
#		"""
#		returns a result object with only the SNPs, etc. within the given boundary.
#		"""
#		positions = []
#		scores = []
#		snps = []
#		mafs = []
#		marfs = []
#		chromosomes = []
#		i = 0
#		start_pos -= buffer
#		end_pos += buffer
#		while i < len(self.chromosomes) and self.chromosomes[i] != chromosome:
#			i += 1
#		while i < len(self.chromosomes) and self.positions[i] < start_pos:
#			i += 1
#		if i == len(self.chromosomes):
#			raise Exception("region results never found!")
#
#		while i < len(self.chromosomes) and self.positions[i] <= end_pos and self.chromosomes[i] == chromosome:
#			chromosomes.append(chromosome)
#			positions.append(self.positions[i])
#			scores.append(self.scores[i])
#			if self.snps:
#				snps.append(self.snps[i])
#			mafs.append(self.mafs[i])
#			marfs.append(self.marfs[i])
#			i += 1
#
#		return Result(result_type=self.result_type, phen_id=self.phen_id,
#			      scores=scores, chromosomes=chromosomes, positions=positions,
#			      snps=snps, accessions=self.accessions, marfs=marfs, mafs=mafs)
#
#
#
#	def get_top_region_results(self, n, distance_threshold=25000, buffer=25000):
#		reg_results = []
#		i = 0
#		for reg in self.get_top_regions(n, distance_threshold):
#			reg_results.append(self.get_region_result(*reg, buffer=buffer))
#		print "Regions were retrieved."
#		return reg_results
#
#
#	def filter_top_snps(self, n):
#		"""
#		Filter all but top n SNPs, sorted by score.
#		"""
#		self._rank_scores_() #Making sure the ranks are updated
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		new_snps = []
#		l = self.orders[0:n]
#		l.sort()
#
#		for i in l:
#			newScores.append(self.scores[i])
#			newPositions.append(self.positions[i])
#			newChromosomes.append(self.chromosomes[i])
#			newMafs.append(self.mafs[i])
#			newMarfs.append(self.marfs[i])
#			if self.snps:
#				new_snps.append(self.snps[i])
#
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#		self.snps = new_snps
#
#
#
#	def get_chromosome_ends(self):
#		if not self.chromosome_ends:
#			self._sort_by_chr_pos_()
#			i = 0
#			chromosome_ends = []
#			while i < len(self.chromosomes):
#				curr_chr = self.chromosomes[i]
#				while i < len(self.chromosomes) and self.chromosomes[i] == curr_chr:
#					i += 1
#				chromosome_ends.append(self.positions[i - 1])
#			self.chromosome_ends = chromosome_ends
#		print self.chromosome_ends
#		return self.chromosome_ends
#
#
#
#	def clone(self):
#		import copy
#		result = copy.deepcopy(self) #Cloning
#		return result
#
#
#	def simple_clone(self):
#		result = Result(scores=self.scores[:], positions=self.positions[:], chromosomes=self.chromosomes[:],
#				marfs=self.marfs[:], mafs=self.mafs[:], accessions=self.accessions[:])
#		result.chromosome_ends = self.chromosome_ends[:]
#		return result
#
#
#
#
##	def insert_into_db(self,result_id,host='gmi-ara-devel-be'):
##		"""
##		Insert pvalues into DB.. (TEST) 
##		"""
##		columns = ['result_id', 'chromosome', 'position', 'pValue']
##		if self.marfs:
##			columns.append('MAF')
##		if self.mafs:
##			columns.append('MAC')
##		import dbutils
##		sql_statement_prefix = "INSERT INTO stock_snp_test.result_pvalues ("+",".join(columns)+")"
##		c_p_pv_list = self.get_chr_pos_score_list()
##		print len(c_p_pv_list)
##		conn = dbutils.connect_to_db(host,'stock_snp_test')
##		cursor = conn.cursor()
##		for i, (c,p,pv) in enumerate(c_p_pv_list):
##			sql_statement = sql_statement_prefix+"VALUES (%d, %d, %d, %f"%(result_id,c,p,pv)
##			if self.marfs and self.mafs:
##				sql_statement +=", %f, %d)"%(self.marfs[i],self.mafs[i])
##			else:
##				sql_statement +=")"
##		
##			#print sql_statement
##			numRows = int(cursor.execute(sql_statement))
##			row = cursor.fetchone()
##			if row:
##				print row
##		conn.commit()
##		cursor.close ()
##		conn.close ()
#
#
#
#	def na_mafs(self, min_maf=10):
#		"""
#		NA scores/pvalues which have maf<minMaf.		
#		"""
#		for i in range(0, len(self.scores)):
#			if self.mafs[i] < min_maf:
#				self.scores[i] = "NA"
#
#
#	def get_max_snp(self):
#		max_val = max(self.scores)
#		mi = self.scores.index(max_val)
#		return (self.snps[mi], self.scores[mi], self.chromosomes[mi], self.positions[mi])
#
#
#
#	def write_to_file_old(self, filename, with_extra_info=False):
#		header = ['chromosomes', 'positions', 'scores', 'marfs', 'mafs']
#		if with_extra_info:
#			for info in self.snp_results:
#				if not info in header:
#					header.append(info)
#		f = open(filename, "w")
#
#		f.write("Chromosome,Position,Score,MARF,MAF \n")
#
#		for (ch, pos, score, marf, maf) in zip(self.chromosomes, self.positions, self.scores, self.marfs, self.mafs):
#			l = map(str, [ch, pos, score, marf, maf])
#			f.write(",".join(l) + "\n")
#		f.close()
#
#
#	def write_to_file(self, filename, additional_columns=None):
#		columns = ['chromosomes', 'positions', 'scores', 'marfs', 'mafs']
#		if additional_columns:
#			for info in additional_columns:
#				if info in self.snp_results:
#					columns.append(info)
#		try:
#			f = open(filename, "w")
#
#			f.write(','.join(columns) + "\n")
#
#			for i in range(len(self.snp_results[columns[0]])):
#				l = []
#				for c in columns:
#					l.append(self.snp_results[c][i])
#				l = map(str, l)
#				f.write(",".join(l) + "\n")
#			f.close()
#		except Exception, err_str:
#			print 'Failed writing the resultfile:', err_str
#			print 'Make sure the given path is correct, and you have write rights.'


class Result(object):
	"""
	Contains information on the result.  (The old object is renamed Result_old, for now..)  Poised to cause problems?
	"""
	def __init__(self, result_file=None, snp_results=None, scores=None, snps_data=None, accessions=None, name=None,
		     result_type=None, phen_id=None, positions=None, chromosomes=None, mafs=None, macs=None,
		     snps=None, **snp_results_info):
		"""
		A simple init function.
		snps_data is assumed to match the results, (in term of positions, etc)
		
		accessions (when used) should fit the results!
		"""

		if snp_results:
			self.snp_results = snp_results
		else:
			self.snp_results = {}  #Main data container
		self.phen_id = phen_id
		self.result_type = result_type
		self.name = name
		self.keys = []
		self.accessions = accessions
		if result_file:
			self._load_result_(result_file)
		else:

			if chromosomes:
				self.keys.append('chromosomes')
				self.snp_results['chromosomes']	 = chromosomes
			if positions:
				self.keys.append('positions')
				self.snp_results['positions'] = positions
			if scores:
				self.keys.append('scores')
				self.snp_results['scores'] = scores
			if mafs:
				self.keys.append('mafs')
				self.snp_results['mafs'] = mafs
			if macs:
				self.keys.append('macs')
				self.snp_results['macs'] = macs

		self.orders = None
		self.ranks = None
		self.chromosome_ends = []

		if snps_data:
			self._load_snps_(snps_data)
			if 'mafs' not in self.snp_results:
				self._load_mafs_(snps_data)
		else:
			if snps:
				self.keys.append('snps')
				self.snp_results['snps'] = snps


		#Adding various info...
		if snp_results_info:
			for info in snp_results_info:
				self.keys.append(info)
				self.snp_results[info] = snp_results_info[info]



	def _load_result_(self, result_file, try_delims=[",", "\t", " ", ]):
		with open(result_file, "r") as f:
			print 'Loading results..'
			line = f.next()
			for delim in try_delims:
				if len(line.split(delim)) > 1:
					break
			else:
				raise Exception("Appropriate delimiter wasn't found.")

			header = line.split(delim)
			for key in header:  #Handling old data formats..
				key = key.lower()
				if key == 'chromosome':
					key = 'chromosomes'
				if key == 'position':
					key = 'positions'
				if key == 'score':
					key = 'scores'
				if key == 'marf':
					key = 'mafs'
				if key == 'maf':
					key = 'macs'
				self.keys.append(key)

			for key in self.keys:
				self.snp_results[key] = []

			chrom = 0
			chrom_splits = []
			for i, line in enumerate(f):
				l = map(float, map(str.strip, line.split(delim)))
				for v, key in it.izip(l, keys):
					self.snp_results[key].append(v)
				if int(l[0]) != chrom:
					chrom_splits.append(i)
					chrom = int(l[0])

			for si in chrom_splits[1:]:
				self.chromosome_ends.append(self.snp_results['positions'][si - 1])
			self.chromosome_ends.append(self.snp_results['positions'][-1])

			self._rank_scores_()


	def _load_snps_(self, snps_data):
		"""
		This only loads SNPs...
		"""
		snps = []
		positions = []
		chromosomes = []
		for i, snpsd in enumerate(snps_data.snpsDataList):
			snps.extend(snpsd.snps)
			pl = snpsd.positions
			positions.extend(pl)
			chromosomes.extend([snps_data.chromosomes[i]] * len(pl))
			self.chromosome_ends.append(snpsd.positions[-1])
		self.snp_results['snps'] = snps
		self.snp_results['positions'] = positions
		self.snp_results['chromosomes'] = chromosomes
		self.keys.append('chromosomes')
		self.keys.append('positions')
		self.keys.append('snps')


	def _load_mafs_(self, snps_data):
		snps = self.snp_results['snps']
		if snps_data.data_format == 'float':
			self.snp_results['mafs'] = [1] * len(snps)
			self.snp_results['macs'] = [1] * len(snps)
		else:
			maf_d = snps_data.get_mafs()
			self.snp_results['mafs'] = maf_d['marfs']
			self.snp_results['macs'] = maf_d['mafs']
		self.keys.append('mafs')
		self.keys.append('macs')

	def _rank_scores_(self):
		"""
		Generates two data structures:
		self.orders (SNP indices ordered by rank)
		self.ranks (ranks of the SNPs)
		"""
		rank_ls = zip(self.scores, range(len(self.scores)))
		rank_ls.sort(reverse=True)
		self.orders = []
		for j in range(len(rank_ls)):
			(s, i) = rank_ls[j]
			self.orders.append(i)
			rank_ls[j] = (i, j)

		rank_ls.sort()
		self.ranks = []
		for (i, j) in rank_ls:
			self.ranks.append(j + 1)


	def _sort_by_chr_pos_(self):
		res_ls = []
		for key in self.keys:
			res_ls.append(self.snp_results[key])
		res_ls.sort()

		snp_results = {}
		for ls, key in it.izip(res_ls, self.keys):
			snp_results[key] = ls


	def candidate_gene_enrichments(self, cgl=None, cgl_file=None, pval_thresholds=[0.01], gene_radius=20000,
				methods=['chi_square'], num_perm=500, file_prefix=None,
				obs_genes_file=None, early_stop_threshold=25, all_genes=None, cand_gene_indices=None):
		"""
		Performs CGR analysis on this results object.
		
		cgl is a list of genes.
		cgl_file is a file with a list of genes.
		
		method: 
			chi_square - statistics
			multinomial - statistics
			gene_perm - permute genes.
			snps_perm - permute SNPs
		"""
		import pylab
		import analyze_gene_enrichment as genr

		if not 'chi_square' in methods:
			methods.append('chi_square')

		chrom_ends = self.get_chromosome_ends()
		print 'chrom_ends', chrom_ends

		#Parse cgl file
		if cgl_file:
			cgl, cg_tair_ids = load_cand_genes_file(cgl_file)
		for cg in cgl:
			print str(cg)


		if not all_genes:
			#Load genes from DB.
			print 'Fetching all genes'
			all_genes = get_gene_list(include_intron_exons=False, verbose=False)
			print 'Fetched %d genes.' % len(all_genes)
			num_genes = len(all_genes)

		if not cand_gene_indices:
			#Pre-process cgl.
			cand_gene_indices = []
			for i, g in enumerate(all_genes):
				for cg in cgl:
					if g.dbRef == cg.dbRef:
						#print g.dbRef, cg.dbRef
						cand_gene_indices.append(i)
						break
			num_cand_genes = len(cand_gene_indices)
			#print num_cand_genes, cand_gene_indices


		method_res_dict = {}
		for m in methods:
			method_res_dict[m] = {'statistics':[], 'pvals':[]}

		if obs_genes_file:
			obs_gene_str = ''

		pval_thresholds.sort()
		last_thres = 1.0
		log_enrichments = []
		for pval_threshold in reversed(pval_thresholds):
			print 'Using p-value threshold % f' % pval_threshold
			thres = pval_threshold / last_thres
			print 'Using corrected threshold % f' % thres
			last_thres = pval_threshold

			#Filter pvalue file
			self.filter_percentile(1 - thres)

			#pre-process pvalues
			regions = self.get_regions(gene_radius=gene_radius)

			#Calculate observed candidate gene enrichment. 
			obs_enrichments = genr.calc_enrichment(all_genes, cand_gene_indices, regions)
			r1 = obs_enrichments[0] / float(obs_enrichments[1])
			r2 = (num_cand_genes / float(num_genes))
			obs_stat = sp.log(r1 / r2)
			print 'Observed statistics % f' % obs_stat
			log_enrichments.append(obs_stat)

			#What cand. genes overlap with regions?
			obs_cg_indices = obs_enrichments[2]
			if obs_cg_indices:
				obs_gene_str += str(pval_threshold) + ','
				tair_ids = [all_genes[cgi].tairID for cgi in obs_cg_indices]
				obs_gene_str += ','.join(tair_ids)
				obs_gene_str += '\n'



			for method in methods:
				if method == 'chi_square':
					chi_sq_pval, chi_sq_stat = genr.get_chi_square_pval(obs_enrichments[0],
											obs_enrichments[1],
											num_cand_genes, num_genes)
					method_res_dict[method]['statistics'].append(chi_sq_stat)
					method_res_dict[method]['pvals'].append(chi_sq_pval)

				if method == 'multinomial':
					pass
				if method == 'gene_perm':
					p_val, perm_stats = genr.get_gene_perm_pval(obs_stat, regions, all_genes,
									cand_gene_indices, num_perm=num_perm,
									early_stop_threshold=early_stop_threshold)
					method_res_dict[method]['statistics'].append(perm_stats)
					method_res_dict[method]['pvals'].append(p_val)

				if method == 'snps_perm':
					p_val, perm_stats = genr.get_snps_perm_pval(obs_stat, regions, all_genes,
									cand_gene_indices, chrom_ends,
									num_perm=num_perm,
									early_stop_threshold=early_stop_threshold)

					method_res_dict[method]['statistics'].append(perm_stats)
					method_res_dict[method]['pvals'].append(p_val)

#						h_res = pylab.hist(perm_stats)
#						pylab.vlines(obs_stat, 0, max(h_res[0]), colors='r')
#						pylab.savefig(env.env['tmp_dir'] + 'test.pdf', format='pdf')


		if obs_genes_file:
			with open(obs_genes_file, 'w') as f:
				f.write(obs_gene_str)

		#Now the plotting of the results.
		method_name_dict = {'chi_square':'Chi-square test', 'gene_perm':'Candidate gene permutations',
					'snps_perm':'SNP positions permutation (chromosome rotation)' }


		pval_thresholds.reverse()
		pos_list = range(len(pval_thresholds))
		for m in methods:
			pylab.figure()
			pvals = []
			for p in method_res_dict[m]['pvals']:
				if p != 0:
					pvals.append(p)
				else:
					pvals.append(0.5 / num_perm)
			neg_log_pvals = map(lambda x:-math.log10(x), pvals)
			pylab.barh(pos_list, neg_log_pvals, align='center', color='g', alpha=0.6)
			pylab.ylim((-1, len(pval_thresholds)))
			pylab.yticks(pos_list, map(str, pval_thresholds))
			pylab.xlabel('Enrichment -log(p-value)')
			pylab.ylabel('p-value percentile threshold')
			ymin, ymax = pylab.ylim()
			xmin, xmax = pylab.xlim()
			pylab.axvline(-math.log10(0.05), ymin=ymin, ymax=ymax, color='r')
			pylab.xlim((0, max(-math.log10(0.05), max(neg_log_pvals)) * 1.05))
			pylab.title(method_name_dict[m])
			pylab.savefig(file_prefix + '_' + m + '.png', format='png')


		pylab.figure()
		pylab.barh(pos_list, log_enrichments, align='center', color='g', alpha=0.6)
		pylab.ylim((-1, len(pval_thresholds)))
		pylab.yticks(pos_list, map(str, pval_thresholds))
		pylab.xlabel('Enrichment (log[ratio])')
		pylab.ylabel('p-value percentile threshold')
		pylab.savefig(file_prefix + '_enrichment_ratio.png', format='png')

		return {'enr_stats':log_enrichments, 'method_res_dict':method_res_dict}




	def _plot_small_manhattan_(self, pdf_file=None, png_file=None, min_score=None, max_score=None,
				type="pvals", ylab="$-$log$_{10}(p-$value$)$", plot_bonferroni=False,
				cand_genes=None, threshold=0, highlight_markers=None, chromosome=None):
		import matplotlib
		#matplotlib.use('Agg')
		import matplotlib.pyplot as plt

		displayed_unit = 1000.0 #kbs
		scoreRange = max_score - min_score
		plt.figure(figsize=(6, 3.5))
		plt.axes([0.045, 0.12, 0.95, 0.7])
		starPoints = [[], [], []]
		color = 'b'
		color_map = {1:'b', 2:'g', 3:'r', 4:'c', 5:'m'}
		if chromosome:
			color = color_map[chromosome]

		chrom = self.snp_results['chromosomes'][0]
		positions = map(lambda x: x / displayed_unit, self.snp_results['positions'])
		scores = self.snp_results['scores']
		for s_i, (score, pos) in enumerate(it.izip(scores, positions)):
			if score > max_score:
				starPoints[0].append(pos)
				starPoints[1].append(max_score)
				starPoints[2].append(score)
				score = max_score
			scores[s_i] = score

		plt.plot(positions, scores, ".", markersize=3, alpha=0.7, color=color)

		if cand_genes:
			for cg in cand_genes:
				plt.axvspan(cg.startPos / displayed_unit, cg.endPos / displayed_unit, facecolor='k', alpha=0.5, linewidth=0)


		if highlight_markers:
			ys = []
			xs = []
			for c, p, score in highlight_markers:
				xs.append(p / displayed_unit)
				if score > max_score:
					plt.text(x, max_score * 1.1, str(round(score, 2)), rotation=45, size="small")
					ys.append(max_score)
				else:
					ys.append(score)
			plt.plot(xs, ys, ".", color="#ff9944", markersize=6, alpha=0.7)

		if len(starPoints[0]) > 0:
			plt.plot(starPoints[0], starPoints[1], ".", color="#ee9922", markersize=4)
			i = 0
			while i < len(starPoints[0]):
				max_point = i
				cur_pos = starPoints[0][i]
				while i < len(starPoints[0]) and abs(starPoints[0][i] - cur_pos) < 1000000:
					if starPoints[2][i] > starPoints[2][max_point]:
						max_point = i
					i += 1
				plt.text(starPoints[0][max_point] - 200000, (starPoints[1][max_point] - 1) * 1.15, str(round(starPoints[2][max_point], 2)), rotation=45, size="small")




		if plot_bonferroni:
			b_threshold = -math.log10(1.0 / (len(scores) * 20.0))
			if threshold :
				plt.plot([0, max(positions)], [b_threshold, b_threshold], ":")
				threshold = -math.log10(threshold)
				plt.plot([0, max(positions)], [threshold, threshold], color='#6495ed', linestyle='-.')
			#Bonferroni threshold
			else:
				plt.plot([0, max(positions)], [b_threshold, b_threshold], color='#000000', linestyle="-.")

		plt.axis([min(positions), max(positions), min_score - 0.05 * scoreRange, max_score + 0.05 * scoreRange])
		if not ylab:
			if type == "pvals":
				plt.ylabel('$ - log(p - $value$)$')

			else:
				plt.ylabel('score')
		else:
			plt.ylabel(ylab)
		plt.xlabel("kilobases")
		plt.title('Chromsome %d' % chrom)

		if pdf_file:
			plt.savefig(pdf_file, format="pdf")
		if png_file:
			plt.savefig(png_file, format="png", dpi=300, bbox_inches='tight')
		if not (pdf_file or png_file):
			plt.show()

		plt.clf()
		plt.close()


	def plot_manhattan(self, pdf_file=None, png_file=None, min_score=None, max_score=None,
		       percentile=90, type="pvals", ylab="$-$log$_{10}(p-$value$)$",
		       plot_bonferroni=False, cand_genes=None, threshold=0, highlight_markers=None):

		"""
		Plots a 'Manhattan' style GWAs plot.
		"""

		import matplotlib
		#matplotlib.use('Agg')
		import matplotlib.pyplot as plt
		num_scores = len(self.snp_results['scores'])

		"Plotting a Manhattan-style plot with %i markers." % num_scores

		chromosome_ends = self.get_chromosome_ends()
		result = self.simple_clone()

		chrom_set = set(result.snp_results['chromosomes'])
		chromosomes = list(chrom_set)
		chromosomes.sort()
		if len(chrom_set) == 1:
			percentile = 0.0
		if percentile != 0.0:
			result.filter_percentile(percentile / 100.0)

		if highlight_markers:
			new_h_markers = []
			for c, p, pval in highlight_markers:
				new_h_markers.append((c, p, -math.log10(pval)))
			highlight_markers = new_h_markers

		if not max_score:
			max_score = max(result.snp_results['scores'])
			if highlight_markers:
				h_scores = [s for c, p, s in highlight_markers]
				max_score = max(max_score, max(h_scores))
		if not min_score:
			if type == "pvals":
				min_score = 0
			else:
				min_score = min(result.snp_results['scores'])


		if cand_genes:
			#processing candidate genes by chromosome
			chr_cand_genes = {}
			for chrom in chromosomes:
				chr_cand_genes[chrom] = []
			for cg in cand_genes:
				chr_cand_genes[cg.chromosome].append(cg)

		if len(chrom_set) == 1:
			chrom = chrom_set.pop()
			if cand_genes:
				cand_genes = chr_cand_genes[chrom]
			return result._plot_small_manhattan_(pdf_file=pdf_file, png_file=png_file, min_score=min_score,
						max_score=max_score, ylab=ylab, plot_bonferroni=plot_bonferroni,
						cand_genes=cand_genes, threshold=threshold,
						highlight_markers=highlight_markers, chromosome=chrom)


		scoreRange = max_score - min_score
		offset = 0
		chromosome_splits = result.get_chromosome_splits()

		ticksList1 = []
		ticksList2 = []
		textPos = []
		plt.figure(figsize=(12, 2.8))
		plt.axes([0.045, 0.15, 0.95, 0.71])
		starPoints = [[], [], []]
		chr_offsets = []
		for i, chromosome_end in enumerate(chromosome_ends):
			chr_offsets.append(offset)
			index1 = chromosome_splits[i]
			index2 = chromosome_splits[i + 1]
			scoreList = result.snp_results['scores'][index1:index2]
			posList = result.snp_results['positions'][index1:index2]
			chrom = chromosomes[i]
			newPosList = [offset + pos for pos in posList]

			for s_i, (score, pos) in enumerate(it.izip(scoreList, newPosList)):
				if score > max_score:
					starPoints[0].append(pos)
					starPoints[1].append(max_score)
					starPoints[2].append(score)
					score = max_score
				scoreList[s_i] = score

			plt.plot(newPosList, scoreList, ".", markersize=2, alpha=0.7)

			if cand_genes:
				for cg in chr_cand_genes[chrom]:
					plt.axvspan(offset + cg.startPos, offset + cg.endPos, facecolor='#FF9900', alpha=0.6)

			oldOffset = offset
			textPos.append(offset + chromosome_end / 2 - 2000000)
			offset += chromosome_end
			for j in range(oldOffset, offset, 2000000):
				ticksList1.append(j)
			#pdb.set_trace()
			for j in range(0, chromosome_end, 2000000):
				if j % 4000000 == 0 and j < chromosome_end - 2000000 :
					ticksList2.append(j / 1000000)
				else:
					ticksList2.append("")
			#pdb.set_trace()



		plt.plot(starPoints[0], starPoints[1], ".", color="#ee9922", markersize=4)
		if len(starPoints[0]) > 0:
			i = 0
			while i < len(starPoints[0]):
				max_point = i
				cur_pos = starPoints[0][i]
				while i < len(starPoints[0]) and abs(starPoints[0][i] - cur_pos) < 3000000:
					if starPoints[2][i] > starPoints[2][max_point]:
						max_point = i
					i += 1
				plt.text(starPoints[0][max_point] - 1000000, (starPoints[1][max_point] - 1) * 1.15, str(round(starPoints[2][max_point], 2)), rotation=45, size="small")


		if highlight_markers:
			ys = []
			xs = []
			for c, p, score in highlight_markers:
				x = chr_offsets[c - 1] + p
				xs.append(x)
				if score > max_score:
					plt.text(x, max_score * 1.1, str(round(score, 2)), rotation=45, size="small")
					ys.append(max_score)
				else:
					ys.append(score)
			plt.plot(xs, ys, ".", color="#ff9944", markersize=6, alpha=0.8)


		if plot_bonferroni:
			b_threshold = -math.log10(1.0 / (num_scores * 20.0))
			if threshold :
				plt.plot([0, sum(result.chromosome_ends)], [b_threshold, b_threshold], ":")
				threshold = -math.log10(threshold)
				plt.plot([0, sum(result.chromosome_ends)], [threshold, threshold], color='#6495ed', linestyle='-.')
			#Bonferroni threshold
			else:
				plt.plot([0, sum(result.chromosome_ends)], [b_threshold, b_threshold], color='#000000', linestyle="-.")

		plt.axis([0, sum(result.chromosome_ends), min_score - 0.05 * scoreRange, max_score + 0.05 * scoreRange])
		plt.xticks(ticksList1, ticksList2, fontsize='x-small')
		#pdb.set_trace()
		if not ylab:
			if type == "pvals":
				plt.ylabel('$ - log(p - $value$)$')

			else:
				plt.ylabel('score')
		else:
			plt.ylabel(ylab)
		plt.xlabel("Mb")

		if pdf_file:
			plt.savefig(pdf_file, format="pdf")
		if png_file:
			plt.savefig(png_file, format="png", dpi=300, bbox_inches='tight')
		if not (pdf_file or png_file):
			plt.show()

		plt.clf()
		plt.close()



	def get_chromosome_splits(self):
		"""
		Returns list of indices (and prev chromosome), for the when the chromosomes
		change in the scores, positions indices.
		"""
		last_chrom = 0
		chromosome_splits = []
		for i, chrom in enumerate(self.snp_results['chromosomes']):
			if last_chrom != chrom:
				while last_chrom < chrom:
					last_chrom += 1
					chromosome_splits.append(i)
		chromosome_splits.append(i)
		return chromosome_splits


	def neg_log_trans(self):
		"""
		Apply - log(x) to the pvalues (scores)
		"""
		import math
		f = lambda x: 323.0 if x == 0.0 else - math.log10(x)
		self.snp_results['scores'] = map(f, self.snp_results['scores'])


	def filter_percentile(self, percentile, reversed=False):
		"""
		Filter above (or below) percentile.
		"""
		sl = self.snp_results['scores'][:]
		sl.sort()
		score_cutoff = sl[int(len(sl) * percentile)]
		self.filter_attr("scores", score_cutoff, reversed=reversed)


	def filter_common_top_scores(self, result, top_fraction=0.1):
		import bisect
		self._rank_scores_()
		result._rank_scores_()
		keep_indices_1 = set()
		keep_indices_2 = set()
		chr_pos_list_1 = self.get_chr_pos_list()
		chr_pos_list_2 = result.get_chr_pos_list()
		for i in self.orders[:int(top_fraction * len(self.orders))]:
			j = bisect.bisect(chr_pos_list_2, chr_pos_list_1[i]) - 1
			if chr_pos_list_2[j] == chr_pos_list_1[i]:
				keep_indices_1.add(i)

		for i in result.orders[:int(top_fraction * len(result.orders))]:
			j = bisect.bisect(chr_pos_list_1, chr_pos_list_2[i]) - 1
			if chr_pos_list_1[j] == chr_pos_list_2[i]:
				keep_indices_2.add(i)

		keep_indices_list_1 = list(keep_indices_1)
		for i in keep_indices_2:
			keep_indices_1.add(bisect.bisect(chr_pos_list_1, chr_pos_list_2[i]) - 1)

		for i in keep_indices_list_1:
			keep_indices_2.add(bisect.bisect(chr_pos_list_2, chr_pos_list_1[i]) - 1)

		print len(keep_indices_1), len(keep_indices_2)
		self.filter_indices(list(keep_indices_1))
		result.filter_indices(list(keep_indices_2))
		return keep_indices_1, keep_indices_2


	def filter_indices(self, indices_to_keep):
		indices_to_keep.sort()
		snp_results = {}
		for info in self.snp_results:
			snp_results[info] = []
		count = len(self.scores)
		for i in indices_to_keep:
			for info in self.snp_results:
				if self.snp_results[info]:
					snp_results[info].append(self.snp_results[info][i])

		self.snp_results = snp_results
		print "%i scores were removed." % (count - len(self.scores))


	def filter_attr(self, attr_name, attr_threshold, reversed=False):
		"""
		Filter out scores / pvalues etc. which have attr < attr_threshold.

		attr are e.g.
		'mafs', 'macs', 'scores', etc.
		"""
		print "Filtering for attribute ' % s' with threshold: %g" % (attr_name, attr_threshold)
		attr = self.snp_results[attr_name]
		snp_results = {}
		for info in self.snp_results:
			snp_results[info] = []
		count = len(self.snp_results['scores'])
		for i in range(count):
			if reversed:
				if attr[i] <= attr_threshold:
					for info in self.snp_results:
						if len(self.snp_results[info]) > 0:
							snp_results[info].append(self.snp_results[info][i])
			else:
				if attr[i] >= attr_threshold:
					for info in self.snp_results:
						if len(self.snp_results[info]) > 0:
							snp_results[info].append(self.snp_results[info][i])

		self.snp_results = snp_results
		num_scores = len(self.snp_results['scores'])
		print "%i scores were removed out of %i." % (count - num_scores, count)
		return num_scores



#	def filter_non_segregating_snps(self, ecotype1, ecotype2, accessions=None):
#		"""
#		Filter out all SNPs which are not segregating in the two accessions.
#
#		Assumes the accessions map the results objects SNPs. (and that they are defined)
#		"""
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		new_snps = []
#
#		if accessions:
#			ecotypes = accessions
#		else:
#			ecotypes = self.accessions
#
#		e_i1 = ecotypes.index(ecotype1)
#		e_i2 = ecotypes.index(ecotype2)
#
#		for i in range(len(self.snps)):
#			snp = self.snps[i]
#			if snp[e_i1] != snp[e_i2]:
#				newScores.append(self.scores[i])
#				newPositions.append(self.positions[i])
#				newChromosomes.append(self.chromosomes[i])
#				newMafs.append(self.mafs[i])
#				newMarfs.append(self.marfs[i])
#				new_snps.append(snp)
#
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#		self.snps = new_snps


	def get_top_snps_result(self, n):
		"""
		returns top n SNPs
		"""
		import copy
		result = copy.deepcopy(self) #Cloning
		result.filter_top_snps(n)
		return result


	def get_top_snps(self, n):
		"""
		returns top n SNPs
		"""
		self._rank_scores_() #Making sure the ranks are updated
		chr_pos_list = []
		for i in self.orders[:n]:
			chr_pos_list.append((self.snp_results['chromosomes'][i], self.snp_results['positions'][i]))
		return chr_pos_list


	def get_top_genes(self, n, window_size=5000, conn=None):
		"""
		Returns a set of (chromosome, start_pos, end_pos), for genes found.
		"""
		self._rank_scores_() #Making sure the ranks are updated
		genes = set()
		snp_ix = []
		i = 0
		while len(genes) < n:
			snp_i = self.orders[i]
			p = self.snp_results['positions'][snp_i]
			c = self.snp_results['chromosomes'][snp_i]
			c_genes = get_gene_list(p - window_size, p + window_size, c, include_intron_exons=False, conn=conn)
			for g in c_genes:
				genes.add((g.chromosome, g.startPos, g.endPos, g.tairID))
			#print len(genes)
			i += 1
		print 'found % d genes' % len(genes)
		return genes


	def get_chr_pos_score_list(self):
		return zip(self.snp_results['chromosomes'], self.snp_results['positions'], self.snp_results['scores'])


	def get_chr_pos_list(self):
		return zip(self.snp_results['chromosomes'], self.snp_results['positions'])


	def get_top_regions(self, n, distance_threshold=25000):
		"""
		Returns a list of regions, defined by (chromosome, start_pos, end_pos).
		"""
		self._rank_scores_()
		chromosome_ends = self.get_chromosome_ends()

		i = 0 #how many SNPs are needed
		region_count = 0
		regions = [[], [], [], [], []]  #a list of (chromosome,start,stop)
		while sum(map(len, regions)) < n:
			snp_i = self.orders[i]
			snp_pos = self.snp_results['positions'][snp_i]
			chromosome = self.snp_results['chromosomes'][snp_i]
			new_snp_reg = (chromosome, max(snp_pos - distance_threshold, 0), min(snp_pos + distance_threshold, chromosome_ends[chromosome - 1]))
			covered = False
			for reg_i, reg in enumerate(regions[chromosome - 1]):
				if (new_snp_reg[1] < reg[2] and new_snp_reg[2] > reg[2]) or (new_snp_reg[1] < reg[1] and new_snp_reg[2] > reg[1]): #then merge					
					regions[chromosome - 1].pop(reg_i) #Removing the region which we're merging with.
					new_snp_reg = (reg[0], min(reg[1], new_snp_reg[1]), max(reg[2], new_snp_reg[2]))
				elif (new_snp_reg[1] >= reg[1] and new_snp_reg[2] <= reg[2]):
					covered = True  #It's already covered
					break
			if not covered:
				regions[chromosome - 1].append(new_snp_reg)
			i += 1
		print "It took %i SNPS to find %i regions." % (i, n)
		regions_list = regions[0] + regions[1] + regions[2] + regions[3] + regions[4]
		return regions_list


	def get_regions(self, gene_radius=20000):
		"""
		Returns a list of [[chr,start],[chr,end]] elements, 
		"""
		regions = []
		num_scores = len(self.snp_results['scores'])
		iter = enumerate(it.izip(self.snp_results['chromosomes'], self.snp_results['positions']))
		th = sp.array([0, gene_radius])
		(pos_i, t) = iter.next()
		chr_pos = sp.array(t)
		while pos_i < num_scores - 1:
			start_chr_pos = chr_pos
			end_chr_pos = start_chr_pos
			while pos_i < num_scores - 1 and sp.all((chr_pos - end_chr_pos) <= th):
				end_chr_pos = chr_pos
				(pos_i, t) = iter.next()
				chr_pos = sp.array(t)
			if tuple(end_chr_pos) < tuple(start_chr_pos):
				pdb.set_trace()
			if pos_i == num_scores - 1: #Last one..
				if sp.all((chr_pos - end_chr_pos) <= th):
					regions.append([start_chr_pos - th, chr_pos + th])
				else:
					regions.append([start_chr_pos - th, end_chr_pos + th])
					regions.append([chr_pos - th, chr_pos + th])

			else:
				regions.append([start_chr_pos - th, end_chr_pos + th])
		start_chr_pos = chr_pos
		end_chr_pos = start_chr_pos
		regions.append([start_chr_pos - th, end_chr_pos + th])
		print 'Found % d regions' % len(regions)
		return regions


	def get_region_result(self, chromosome, start_pos, end_pos, buffer=0):
		"""
		returns a result object with only the SNPs, etc. within the given boundary.
		"""
		snp_results = {}
		for k in keys:
			snp_results[k] = []
		i = 0
		start_pos -= buffer
		end_pos += buffer
		chromosomes = self.snp_results['chromosomes']
		positions = self.snp_results['positions']
		while i < len(chromosomes) and chromosomes[i] != chromosome:
			i += 1
		while i < len(chromosomes) and positions[i] < start_pos:
			i += 1
		if i == len(chromosomes):
			raise Exception("region results never found!")

		while i < len(chromosomes) and positions[i] <= end_pos and chromosomes[i] == chromosome:
			for k in keys:
				snp_results[k].append(self.snp_results[k][i])
			i += 1

		return Result(result_type=self.result_type, snp_results=snp_results, phen_id=self.phen_id,
			accessions=self.accessions)



	def get_top_region_results(self, n, distance_threshold=25000, buffer=25000):
		reg_results = []
		i = 0
		for reg in self.get_top_regions(n, distance_threshold):
			reg_results.append(self.get_region_result(*reg, buffer=buffer))
		print "Regions were retrieved."
		return reg_results




	def get_chromosome_ends(self):
		if not self.chromosome_ends:
			self._sort_by_chr_pos_()
			i = 0
			chromosome_ends = []
			chromosomes = self.snp_results['chromosomes']
			while i < len(chromosomes):
				curr_chr = chromosomes[i]
				while i < len(chromosomes) and chromosomes[i] == curr_chr:
					i += 1
				chromosome_ends.append(self.snp_results['positions'][i - 1])
			self.chromosome_ends = chromosome_ends
		return self.chromosome_ends



	def clone(self):
		import copy
		result = copy.deepcopy(self) #Cloning
		return result


	def simple_clone(self):
		snp_results = {}
		for k in self.keys:
			snp_results[k] = self.snp_results[k][:]
		if self.accessions:
			accessions = self.accessions[:]
		else:
			accessions = None
		result = Result(snp_results=snp_results, accessions=accessions)
		result.chromosome_ends = self.chromosome_ends[:]
		return result


#
#
#	def na_mafs(self, min_maf=10):
#		"""
#		NA scores/pvalues which have maf<minMaf.		
#		"""
#		for i in range(0, len(self.scores)):
#			if self.mafs[i] < min_maf:
#				self.scores[i] = "NA"


	def get_max_snp(self):
		max_val = max(self.snp_results['scores'])
		mi = self.scores.index(max_val)
		return (self.snp_results['snps'][mi], self.snp_results['scores'][mi],
			self.snp_results['chromosomes'][mi], self.snp_results['positions'][mi])




	def write_to_file(self, filename, additional_columns=None):
		columns = ['chromosomes', 'positions', 'scores', 'mafs', 'macs']
		if additional_columns:
			for info in additional_columns:
				if info in self.snp_results:
					columns.append(info)
		try:
			f = open(filename, "w")

			f.write(','.join(columns) + "\n")

			for i in range(len(self.snp_results[columns[0]])):
				l = []
				for c in columns:
					l.append(self.snp_results[c][i])
				l = map(str, l)
				f.write(",".join(l) + "\n")
			f.close()
		except Exception, err_str:
			print 'Failed writing the resultfile:', err_str
			print 'Make sure the given path is correct, and you have write rights.'



#
#class SNPResult(Result):
#	"""
#	Contains information on the result.
#	"""
#
#	def _loadSnpsData_(self, snpsds):
#		"""
#		Loads the SNP data.
#		"""
#		self.snps = []
#		self.accessions = snpsds[0].accessions
#		i = 0 #result index
#		chr = -1 # chromosome (index)
#		while i < len(self.scores):
#			if chr != self.chromosomes[i] - 1:
#				chr = self.chromosomes[i] - 1
#				j = 0 #snpsdata index
#			pos = self.positions[i]
#			#print i,chr, j,len(snpsds[chr].positions)
#			while j < len(snpsds[chr].positions) and pos > snpsds[chr].positions[j]:
#				j += 1
#			if j < len(snpsds[chr].positions) and pos == snpsds[chr].positions[j]:
#				self.snps.append(snpsds[chr].snps[j])
#			i += 1
#
#		if pos > snpsds[chr].positions[j]:
#			while j < len(snpsds[chr].positions) and pos > snpsds[chr].positions[j]:
#				j += 1
#			if j < len(snpsds[chr].positions) and pos == snpsds[chr].positions[j]:
#				self.snps.append(snpsds[chr].snps[j])
#
#		if i != len(self.snps):
#			print "Problems with loading SNPs", i, len(self.snps)
#
#		print "Loaded", len(self.snps), " SNPs and", len(self.accessions), "accessions."
#
#
#
#	def filterNicePeaks(self, scoreThreshold, singletonScoreThreshold, window=[20000, 20000], method=1):
#		currScoreWeight = 0.2
#		currChrom = -1
#		newScores = []
#		newSnps = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		singletonCount = 0
#		lastSecondaryEnd = 0
#		if method == 1:
#			for i in range(0, len(self.positions)):
#				if currChrom != self.chromosomes[i]: #Restart
#					currChrom = self.chromosomes[i]
#					startIndex = i #windowIndices 
#					stopIndex = i #windowIndices 
#					curScoreSum = currScoreWeight * self.scores[i]
#					oldScore = 0
#					numSNPs = 1
#				currPos = self.positions[i]
#
#				while currPos - self.positions[startIndex] > window[0]:
#					curScoreSum -= self.scores[startIndex]
#					startIndex += 1
#					numSNPs -= 1
#
#				while stopIndex + 1 < len(self.positions) and self.positions[stopIndex + 1] - currPos < window[1]:
#					stopIndex += 1
#					curScoreSum += self.scores[stopIndex]
#					numSNPs += 1
#
#				curScoreSum -= oldScore
#				oldScore = currScoreWeight * numSNPs * self.scores[i]
#				curScoreSum += oldScore
#
#				if (numSNPs < 5 and self.scores[i] > singletonScoreThreshold) or  (numSNPs > 5 and curScoreSum / float(numSNPs) > scoreThreshold):
#					if (numSNPs < 5 and self.scores[i] > singletonScoreThreshold):
#						singletonCount += 1
#					newScores.append(self.scores[i])
#					newPositions.append(self.positions[i])
#					newChromosomes.append(self.chromosomes[i])
#					newMafs.append(self.mafs[i])
#					newMarfs.append(self.marfs[i])
#					newSnps.append(self.snps[i])
#
#		elif method == 2:
#			for i in range(0, len(self.scores)):
#				if self.scores[i] >= singletonScoreThreshold:
#					newScores.append(self.scores[i])
#					newPositions.append(self.positions[i])
#					newChromosomes.append(self.chromosomes[i])
#					newMafs.append(self.mafs[i])
#					newMarfs.append(self.marfs[i])
#					newSnps.append(self.snps[i])
#
#		# The following code locates the regions before and after the "nice" SNPs.
#		j = 0
#		for i in range(0, len(self.positions)):
#			pos = self.positions[i]
#			chr = self.chromosomes[i]
#			if j < len(newPositions) and pos == newPositions[j] and chr == newChromosomes[j]:
#				k = 0
#				while  i + k > lastSecondaryEnd and pos - self.positions[i + k - 1] < window[0] and self.chromosomes[i + k - 1] == chr:
#					k -= 1
#				while i + k < len(self.positions) - 1 and self.positions[i + k] - pos < window[1] and self.chromosomes[i + k] == chr:
#					if i + k > lastSecondaryEnd:
#						self.secondaryScores.append(self.scores[i + k])
#						self.secondaryPositions.append(self.positions[i + k])
#						self.secondaryChromosomes.append(self.chromosomes[i + k])
#					k += 1
#				lastSecondaryEnd = i + k - 1
#				j += 1
#
#
#		self.snps = newSnps
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#
#		print singletonCount, "singletons were added"
#
#
#
#	def filterMAF(self, minMaf=15):
#		"""
#		Filter out scores/pvalues which have maf<minMaf.		
#		"""
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		newSnps = []
#		for i in range(0, len(self.scores)):
#			if self.mafs[i] >= minMaf:
#				newScores.append(self.scores[i])
#				newPositions.append(self.positions[i])
#				newChromosomes.append(self.chromosomes[i])
#				newMafs.append(self.mafs[i])
#				newMarfs.append(self.marfs[i])
#				newSnps.append(self.snps[i])
#		del self.snps
#		del self.scores
#		del self.positions
#		del self.chromosomes
#		del self.mafs
#		del self.marfs
#
#		self.snps = newSnps
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#
#
#
#	def filterMARF(self, minMarf=0.1):
#		"""
#		Filter out scores/pvalues which have maf<minMaf.		
#		"""
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		newSnps = []
#		for i in range(0, len(self.scores)):
#			if self.marfs[i] >= minMarf:
#				newScores.append(self.scores[i])
#				newPositions.append(self.positions[i])
#				newChromosomes.append(self.chromosomes[i])
#				newMafs.append(self.mafs[i])
#				newMarfs.append(self.marfs[i])
#				newSnps.append(self.snps[i])
#		del self.snps
#		del self.scores
#		del self.positions
#		del self.chromosomes
#		del self.mafs
#		del self.marfs
#
#		self.snps = newSnps
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#
#
#	def filterScoreCutoff(self, scoreCutoff):
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		newSnps = []
#		for i in range(0, len(self.scores)):
#			if self.scores[i] > scoreCutoff:
#				newScores.append(self.scores[i])
#				newPositions.append(self.positions[i])
#				newChromosomes.append(self.chromosomes[i])
#				newMafs.append(self.mafs[i])
#				newMarfs.append(self.marfs[i])
#				newSnps.append(self.snps[i])
#		del self.snps
#		del self.scores
#		del self.positions
#		del self.chromosomes
#		del self.mafs
#		del self.marfs
#
#		self.snps = newSnps
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#
#
#	def getRegionSNPs(self, window=[25000, 25000], snpsDataFile=None):
#		if self.snps:
#			regionSNPs = []
#			self.getRegions(window=window)
#			for (snpIndex, startPos, endPos, chr, size, maxScore, maxPos, snpRank, regionRank) in self.regions:
#				snp = SNP(self.positions[snpIndex], self.chromosomes[snpIndex], alleles=self.snps[snpIndex], accessions=self.accessions, score=maxScore)
#				regionSNPs.append(snp)
#			return regionSNPs
#
#
#	def filterTopSNPs(self, n):
#		self._rankScores_() #Making sure the ranks are updated
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		newSnps = []
#		l = self.orders[0:n]
#		l.sort()
#
#		for i in l:
#			newScores.append(self.scores[i])
#			newPositions.append(self.positions[i])
#			newChromosomes.append(self.chromosomes[i])
#			newMafs.append(self.mafs[i])
#			newMarfs.append(self.marfs[i])
#			newSnps.append(self.snps[i])
#		del self.snps
#		del self.scores
#		del self.positions
#		del self.chromosomes
#		del self.mafs
#		del self.marfs
#
#		self.snps = newSnps
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#
#
#	def _sortByChrPos_(self):
#		res_ls = zip(self.chromosomes, self.positions, self.scores, self.mafs, self.marfs, self.snps)
#		res_ls.sort()
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		newSnps = []
#		for (chr, pos, score, maf, marf, snp) in res_ls:
#			newScores.append(score)
#			newPositions.append(pos)
#			newChromosomes.append(chr)
#			newMafs.append(maf)
#			newMarfs.append(marf)
#			newSnps.append(snp)
#
#		del self.snps
#		del self.scores
#		del self.positions
#		del self.chromosomes
#		del self.mafs
#		del self.marfs
#
#		self.snps = newSnps
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#
#
#	def filterTopRegions(self, n, window=[25000, 25000], minScore=None):
#		self._rankScores_()
#		oldScores = self.scores
#		oldPositions = self.positions
#		oldChromosomes = self.chromosomes
#		oldMafs = self.mafs
#		oldMarfs = self.marfs
#		oldSnps = self.snps
#		oldGlobalRanks = self.globalRanks
#		self.snps = []
#		self.scores = []
#		self.positions = []
#		self.chromosomes = []
#		self.mafs = []
#		self.marfs = []
#		self.globalRanks = []
#
#		regionCount = 0
#		i = 0
#
#		res_ls = []
#		if minScore:
#			while regionCount < n:
#				res_ls.append((oldChromosomes[self.orders[i]], oldPositions[self.orders[i]]))
#				res_ls.sort()
#				if oldScores[self.orders[i]] < minScore:
#					break
#				self.scores.append(oldScores[self.orders[i]])
#				self.positions.append(oldPositions[self.orders[i]])
#				self.chromosomes.append(oldChromosomes[self.orders[i]])
#				self.mafs.append(oldMafs[self.orders[i]])
#				self.marfs.append(oldMarfs[self.orders[i]])
#				self.snps.append(oldSnps[self.orders[i]])
#				self.globalRanks.append(oldGlobalRanks[self.orders[i]])
#				regionCount = self._countRegions_(res_ls, window=window)
#				i += 1
#		else:
#			while regionCount < n:
#				res_ls.append((oldChromosomes[self.orders[i]], oldPositions[self.orders[i]]))
#				res_ls.sort()
#				self.scores.append(oldScores[self.orders[i]])
#				self.positions.append(oldPositions[self.orders[i]])
#				self.chromosomes.append(oldChromosomes[self.orders[i]])
#				self.mafs.append(oldMafs[self.orders[i]])
#				self.marfs.append(oldMarfs[self.orders[i]])
#				self.snps.append(oldSnps[self.orders[i]])
#				self.globalRanks.append(oldGlobalRanks[self.orders[i]])
#				regionCount = self._countRegions_(res_ls, window=window)
#				i += 1
#
#		del oldScores
#		del oldPositions
#		del oldChromosomes
#		del oldMafs
#		del oldMarfs
#		del oldSnps
#		del oldGlobalRanks
#
#		self._sortByChrPos_()
#
##	def get_snps(self):
##		"""
##		Returns this result's SNPs as a list of SNP objects.
##		"""
##		raise  NotImplementedError
##		for i, snp in enumerate(self.snps):
##			SNP(position,chromosome,accessions=None,alleles=None,snpsds=None,score=None,rank=None,regionRank=None)
#
#
#	def mergeWith(self, snpResult):
#		#pdb.set_trace()
#
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		newSnps = []
#
#		if len(snpResult.scores) == 0:
#			return
#		elif len(self.scores) == 0:
#			newPositions = snpResult.positions
#			newChromosomes = snpResult.chromosomes
#			newScores = snpResult.scores
#			newMafs = snpResult.mafs
#			newMarfs = snpResult.marfs
#			newSnps = snpResult.snps
#		else:
#
#			i1 = 0
#
#			pos1 = self.positions[i1]
#			chr1 = self.chromosomes[i1]
#
#			for i2 in range(0, len(snpResult.positions)):
#				pos2 = snpResult.positions[i2]
#				chr2 = snpResult.chromosomes[i2]
#				while i1 < len(self.positions) - 1 and (chr1, pos1) < (chr2, pos2):
#					newPositions.append(self.positions[i1])
#					newChromosomes.append(self.chromosomes[i1])
#					newScores.append(self.scores[i1])
#					newMafs.append(self.mafs[i1])
#					newMarfs.append(self.marfs[i1])
#					newSnps.append(self.snps[i1])
#					i1 += 1
#					pos1 = self.positions[i1]
#					chr1 = self.chromosomes[i1]
#
#				if i1 < len(self.positions) - 1 and (chr1, pos1) == (chr2, pos2):
#					i1 += 1
#					pos1 = self.positions[i1]
#					chr1 = self.chromosomes[i1]
#
#				newPositions.append(snpResult.positions[i2])
#				newChromosomes.append(snpResult.chromosomes[i2])
#				newScores.append(snpResult.scores[i2])
#				newMafs.append(snpResult.mafs[i2])
#				newMarfs.append(snpResult.marfs[i2])
#				newSnps.append(snpResult.snps[i2])
#				#pdb.set_trace()
#
#		del self.snps
#		del self.scores
#		del self.positions
#		del self.chromosomes
#		del self.mafs
#		del self.marfs
#		del snpResult
#
#		self.snps = newSnps
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs


#
#class RegionsTable(object):
#	"""
#	A table or regions X methods/phenotypes.
#	"""
#	def __init__(self, result_ls, window=[25000, 25000]):
#		merged_result = result_ls[0].clone()
#		for i in range(1, len(result_ls)):
#			result = result_ls[i]
#			merged_result.mergeWith(result)
#		merged_result.getRegions(window=window)
#		totalRegSize = 0
#		if len(merged_result.positions):
#			for reg in merged_result.regions:
#				totalRegSize += reg[4]
#		#print "The union of the results: Number of sign. SNPs: "+str(len(merged_result.scores))+", number of sign. regions: "+str(len(merged_result.regions))+", ave. region size: "+str(totalRegSize/float(len(merged_result.regions)))+".\n"
#
#		self.regions = []
#		for (snpIndex, startPos, endPos, chr, size, maxScore, maxPos, snpRank, regionRank) in merged_result.regions:
#			self.regions.append((chr, startPos, endPos, size))
#		del merged_result
#
#		self.region_by_methods_table = [] # list[region_index][method/phenotype_index][snp_index]
#		for (chr, startPos, endPos, size) in self.regions:  #For all regions
#			methods_snps_ls = []
#			for m_i in range(0, len(result_ls)):      #for all results
#				result = result_ls[m_i]
#				snps_ls = []
#
#				if result.snps:
#					##Finding all SNPs in region of interest
#					regionRank = None
#					snpRank = None
#					result.getRegions(window=window)
#					for (si, startPos2, endPos2, chr2, size2, mScore, mPos, sRank, rRank) in result.regions:
#						if chr2 == chr and startPos <= startPos2 and endPos2 <= endPos:
#							#print startPos,startPos2,endPos2,endPos
#							regionRank = rRank
#							snpRank = sRank
#					#if not regionRank:
#						#pdb.set_trace()
#
#				##Finding all SNPs in region of interest
#				for i in range(0, len(result.positions)):
#					if result.chromosomes[i] == chr and startPos < result.positions[i] < endPos:
#						if result.snps:
#							#print result.snps[i]
#							snp = SNP(result.positions[i], chr, alleles=result.snps[i], score=result.scores[i], rank=snpRank, regionRank=regionRank)
#						else:
#							snp = (chr, result.positions[i], result.scores[i], i)  #(chr,pos,score,index)
#						snps_ls.append(snp)
#				methods_snps_ls.append(snps_ls)
#			self.region_by_methods_table.append(methods_snps_ls)
#
#		self.resultNames = []
#		for result in result_ls:
#			self.resultNames.append(result.name)



class SNP(object):
	"""
	A class to represent a SNP. 

	It's used only when analysing a SNP.
	"""

	def __init__(self, position, chromosome, accessions=None, alleles=None, snpsds=None, score=None, rank=None, regionRank=None):
		self.position = position
		self.chromosome = chromosome

		self.alleles = None
		self.accessions = None
		self.score = None
		self.rank = None
		self.regionRank = None

		if not alleles and snpsds:
			self._getAllele_(snpsds)
		else:
			self.alleles = alleles
		if accessions:
			self.accessions = accessions
		if score:
			self.score = score
		if rank:
			self.rank = rank
		if regionRank:
			self.regionRank = regionRank

	def _getAllele_(self, snpsds):
		chr = self.chromosome - 1
		snpsd = snpsds[chr]
		self.snpsdIndex = -1
		self.alleles = None
		for i in range(0, len(snpsd.positions)):
			if snpsd.position == snpsd.positions[i]:
				self.alleles = snpsd.snps[i]
				self.snpsdIndex = i
				break
		if not self.alleles:
			print "The corresponding allele was not found in the data."




class Gene(object):
	"""
	A class which encompasses basic information about a gene.
	"""
	def __init__(self, chromosome=None, startPos=None, endPos=None, name="", description=None, dbRef="", tairID=""):
		self.chromosome = chromosome
		self.startPos = startPos
		self.endPos = endPos
		self.exons = []
		self.introns = []
		self.tairID = tairID
		self.dbRef = dbRef
		self.name = name
		self.description = description
		self.functionDescriptions = []
		self.shortDescriptions = []
		self.direction = None
		self.highlight = False


	def __str__(self):
		if not self.description:
			return "Chromosome=" + str(self.chromosome) + ", position=(" + str(self.startPos) + "," + str(self.endPos) + "), tair ID=" + self.tairID + ", short descriptions=" + str(self.shortDescriptions) + ", function descriptions=" + str(self.functionDescriptions) + "."
		else:
			return "Chromosome=" + str(self.chromosome) + ", position=(" + str(self.startPos) + "," + str(self.endPos) + "), tdbRef=" + self.dbRef + ", name=" + str(self.name) + ", description=" + str(self.description) + "."


	def _update_introns_(self):
		"""
		updates the introns, given the exons.  It uses the Region object..
		"""
		introns = []
		for i in range(len(self.exons) - 1):
			e1 = self.exons[i]
			e2 = self.exons[i + 1]
			intron = Region(self.chromosome, e1.endPos, e2.startPos)
			introns.append(intron)
		self.introns = introns


def getCandidateGeneList(cgl_id, host="papaya.usc.edu", user="bvilhjal", passwd="bjazz32", db="stock_250k", tair9=True):
	import MySQLdb
	#Load cand. gene list.	
	print "Connecting to db, host=" + host
	if not user:
		import sys
		sys.stdout.write("Username: ")
		user = sys.stdin.readline().rstrip()
	if not passwd:
		import getpass
		passwd = getpass.getpass()
	try:
		conn = MySQLdb.connect (host=host, user=user, passwd=passwd, db=db)
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()
	#Retrieve the filenames
	print "Fetching data"

	#select c.locustag, b.start, b.stop, a.comment from genome.gene_commentary a, genome.entrezgene_mapping b, genome.gene c where b.start > 25000 and b.stop < 75000 and b.chromosome=1 and b.gene_id = c.gene_id and c.gene_id = a.gene_id and a.gene_commentary_type_id = 8
	#select distinct t8_fd.tair_id, t8.chromosome, t8.start, t8.end, t8_fd.type, t8_fd.short_description from T8_annotation_TH.t8_063008 t8, T8_annotation_TH.t8_func_desc t8_fd, stock_250k.candidate_gene_list cgl where t8.pub_locus+'.1' = t8_fd.tair_id and cgl.list_type_id=129  and cgl.original_name=t8.pub_locus and t8.chromosome =1 order by t8.chromosome, t8.start
	#select distinct gm.chromosome, gm.start, gm.stop, g.locustag from genome.entrezgene_mapping gm, genome.gene g, stock_250k.candidate_gene_list cgl where cgl.list_type_id=129 and gm.gene_id = g.gene_id and cgl.gene_id=g.gene_id order by gm.chromosome, gm.start, gm.stop

	numRows = int(cursor.execute("select distinct gm.chromosome, gm.start, gm.stop, g.locustag, g.gene_symbol, g.description, g.dbxrefs from genome.entrezgene_mapping gm, genome.gene g, stock_250k.candidate_gene_list cgl where cgl.list_type_id=" + str(cgl_id) + " and gm.gene_id = g.gene_id and cgl.gene_id=g.gene_id order by gm.chromosome, gm.start, gm.stop"))
	if tair9:
		t_map = tc.tair8_to_tair9_map()
	candGenes = []
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		if tair9:
			start_pos = t_map.get_tair9_pos(int(row[1]))
			end_pos = t_map.get_tair9_pos(int(row[2]))
		else:
			start_pos = int(row[1])
			end_pos = int(row[2])
		gene = Gene(int(row[0]), start_pos, end_pos, name=row[4], description=row[5], dbRef=row[6])
		candGenes.append(gene)
	cursor.close ()
	conn.close ()
	print "Candiate gene-lists fetched"
	return candGenes


def get_gene_list(start_pos=None, end_pos=None, chr=None, include_intron_exons=True, \
		verbose=True, conn=None, tair9=True):
	"""
	Fetch genes within a region or all genes from DB.
	"""
	import dbutils
	if not conn:
		new_conn = dbutils.connect_to_default_lookup('genome')
		cursor = new_conn.cursor()
	else:
		cursor = conn.cursor()
	#Retrieve the filenames
	#if verbose:
	#	print "Fetching data"  
	#print "Fetching data"  

	#select c.locustag, b.start, b.stop, a.comment from genome.gene_commentary a, genome.entrezgene_mapping b, genome.gene c where b.start > 25000 and b.stop < 75000 and b.chromosome=1 and b.gene_id = c.gene_id and c.gene_id = a.gene_id and a.gene_commentary_type_id = 8
	#select distinct t8_fd.tair_id, t8.chromosome, t8.start, t8.end, t8_fd.type, t8_fd.short_description from T8_annotation_TH.t8_063008 t8, T8_annotation_TH.t8_func_desc t8_fd, stock_250k.candidate_gene_list cgl where t8.pub_locus+'.1' = t8_fd.tair_id and cgl.list_type_id=129  and cgl.original_name=t8.pub_locus and t8.chromosome =1 order by t8.chromosome, t8.start
	#select distinct gm.chromosome, gm.start, gm.stop, g.locustag from genome.entrezgene_mapping gm, genome.gene g, stock_250k.candidate_gene_list cgl where cgl.list_type_id=129 and gm.gene_id = g.gene_id and cgl.gene_id=g.gene_id order by gm.chromosome, gm.start, gm.stop
	if chr and start_pos and end_pos:
		sql_statement = "SELECT DISTINCT gm.chromosome, gm.start, gm.stop, g.locustag, \
		g.gene_symbol, g.description, g.dbxrefs FROM genome.entrezgene_mapping gm, genome.gene g WHERE \
		gm.gene_id = g.gene_id AND gm.chromosome=" + str(chr) + " AND gm.stop>" + str(start_pos) + " AND \
		gm.start<" + str(end_pos) + " ORDER BY gm.chromosome, gm.start, gm.stop"
		numRows = int(cursor.execute(sql_statement))
	else:
		sql_statement = "select distinct gm.chromosome, gm.start, gm.stop, g.locustag, g.gene_symbol, \
			g.description, g.dbxrefs from genome.entrezgene_mapping gm, genome.gene g \
			where gm.gene_id = g.gene_id order by gm.chromosome, gm.start, gm.stop"
		numRows = int(cursor.execute(sql_statement))
	if numRows == 0:
		pass
		#print sql_statment

	if tair9:
		t_map = tc.tair8_to_tair9_map()
	genes = []
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		#try:
			#chr, start, stop, gene_symbol, description, dbref,  
		if row[1] and  row[2] and row[0] in ['1', '2', '3', '4', '5']:
			chrom = int(row[0])
			if tair9:
				start_pos = t_map.get_tair9_pos(chrom, int(row[1]))
				end_pos = t_map.get_tair9_pos(chrom, int(row[2]))
			else:
				start_pos = int(row[1])
				end_pos = int(row[2])

			gene = Gene(int(row[0]), start_pos, end_pos, name=row[4], description=row[5], dbRef=row[6], tairID=row[3])
			gene.tairID = row[6][5:]
			genes.append(gene)
#		except Exception, err_str:
#			#pass
#			if verbose:
#				print err_str, ':'
#				print row

	if include_intron_exons:
		for g in genes:
			sql_stat = "select distinct gs.start, gs.stop, g.dbxrefs from genome.entrezgene_mapping gm, genome.gene g, genome.gene_segment gs, genome.gene_commentary gc where g.dbxrefs='" + str(g.dbRef) + "' and gm.gene_id = g.gene_id and gs.gene_commentary_id=gc.id and gc.gene_id=gm.gene_id order by gs.start, gs.stop"
			numRows = int(cursor.execute(sql_stat))
			segments = []
			while(1):
				row = cursor.fetchone()
				if not row:
					break;
				try:
					if tair9:
						start_pos = t_map.get_tair9_pos(g.chromosome, int(row[0]))
						end_pos = t_map.get_tair9_pos(g.chromosome, int(row[1]))
					else:
						start_pos = int(row[0])
						end_pos = int(row[1])
					segments.append(Region(g.chromosome, start_pos, end_pos))
				except Exception, err_str:
					print err_str, ':'
					print row
			exons = []
			i = 1
			#print "len(segments):",len(segments)
			while i < len(segments):
				curr_exon = segments[i - 1]
				while i < len(segments) and curr_exon.overlapping(segments[i]):
					curr_exon.merge(segments[i])
					i += 1
				exons.append(curr_exon)
				i += 1
			g.exons = exons
			g._update_introns_()
			#print "len(exons):",len(exons)
			#for e in g.exons:
			#	print e.startPos, e.endPos
			#print "len(g.introns):",len(g.introns)
			#for e in g.introns:
			#	print e.startPos, e.endPos
	cursor.close()
	if not conn:
		new_conn.close ()
	if verbose:
		print "Gene-lists fetched"
	return genes


#"""
#select distinct gm.chromosome, gm.start, gm.stop, g.locustag, g.gene_symbol, g.description, g.dbxrefs, gs.start, gs.stop
#from genome.entrezgene_mapping gm, genome.gene g, gene_segment gs, gene_commentary gc
#where gm.gene_id = g.gene_id and gm.chromosome=5 and gm.stop>3173000 and gm.start<3181000 and gs.gene_commentary_id=gc.id and gc.gene_id=gm.gene_id 
#order by gm.chromosome, gm.start, gm.stop, gs.start, gs.stop
#"""


def load_cand_genes_file(file_name, format=1):
	"""
	Loads a candidate gene list from a csv file...
	"""
	f = open(file_name, "r")
	#print(f.readline())
	reader = csv.reader(f)
	tair_ids = []
	gene_names = []
	if format == 1:
		for row in reader:
		      tair_ids.append(row[0].upper())
		      if len(row) > 1:
		      	gene_names.append(row[1])
	f.close()
	return get_genes_w_tair_id(tair_ids), tair_ids




def get_genes_w_tair_id(tair_ids, tair9=True):
	conn = dbutils.connect_to_default_lookup("genome")
	cursor = conn.cursor()
	if tair9:
		t_map = tc.tair8_to_tair9_map()
	genes = []
	#print tair_ids
	tair_ids.sort()
	for tair_id in tair_ids:
		sql_statment = "select distinct gm.chromosome, gm.start, gm.stop, g.locustag, g.gene_symbol, g.description, g.dbxrefs from genome.entrezgene_mapping gm, genome.gene g where g.dbxrefs='TAIR:" + tair_id.upper() + "' and gm.gene_id = g.gene_id order by gm.chromosome, gm.start, gm.stop"
		#print sql_statment
		numRows = int(cursor.execute(sql_statment))
		if numRows > 1:
			print "Found 2 copies:", sql_statment
		while(1):
			row = cursor.fetchone()
			if not row:
				break;
			try:
				if row[1] and  row[2]:
					chrom = int(row[0])
					if tair9:
						start_pos = t_map.get_tair9_pos(chrom, int(row[1]))
						end_pos = t_map.get_tair9_pos(chrom, int(row[2]))
					else:
						start_pos = int(row[1])
						end_pos = int(row[2])
					#chr, start, stop, gene_symbol, description, dbref,  
					gene = Gene(chrom, start_pos, end_pos, name=row[4], description=row[5], dbRef=row[6], tairID=row[3])
					gene.tairID = row[6][5:]
					genes.append(gene)
			except Exception, err_str:
				pass
				print err_str, ':'
				print row

	cursor.close()
	conn.close()
	return genes


def getResultsFilename(host, user, passwd, callMethodID, phenotypeMethodID, analysisMethodID):
	"""
	Retrieve the filename with the results.
	"""
	import MySQLdb
	print "Connecting to db, host=" + host
	if not user:
		import sys
		sys.stdout.write("Username: ")
		user = sys.stdin.readline().rstrip()
	if not passwd:
		import getpass
		passwd = getpass.getpass()
	try:
		conn = MySQLdb.connect (host=host, user=user, passwd=passwd, db="stock_250k")
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()

	#Retrieve the filenames
	print "Fetching data"
	numRows = int(cursor.execute("select rm.filename from stock_250k.results_method rm where rm.call_method_id=" + str(callMethodID) + " and rm.phenotype_method_id=" + str(phenotypeMethodID) + " and analysis_method_id=" + str(analysisMethodID) + " "))
	filenames = []
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		filenames.append(row[0])
	cursor.close ()
	conn.close ()
	return filenames




def _getStandardResultTypes_():
	res_path = "/Network/Data/250k/tmp-bvilhjal/"
	resultTypes = []
	resultTypes.append(ResultType("KW", ".pvals", "newDataset", res_path + "kw_results/"))
	resultTypes.append(ResultType("Emma", ".pvals", "newDataset", res_path + "emma_results/", mafCutoff=15))
	resultTypes.append(ResultType("Marg", ".score", "newDataset", res_path + "marg_results/"))
	resultTypes.append(ResultType("RF", ".imp", "newDataset", res_path + "rf_results/", mafCutoff=15))
	return resultTypes

def _getStandardResultTypes2_():
	res_path = "/Network/Data/250k/tmp-bvilhjal/"
	resultTypes = []
	resultTypes.append(ResultType("KW", ".pvals", "raw", res_path + "kw_results/"))
	resultTypes.append(ResultType("Emma", ".pvals", "newDataset", res_path + "emma_results/", mafCutoff=15))
	resultTypes.append(ResultType("Marg", ".score", "newDataset", res_path + "marg_results/"))
	resultTypes.append(ResultType("RF", ".imp", "newDataset", res_path + "rf_results/", mafCutoff=15))
	return resultTypes

def _getStandardResultTypes3_():
	res_path = "/Network/Data/250k/tmp-bvilhjal/"
	resultTypes = []
	resultTypes.append(ResultType("KW", ".pvals", "raw", res_path + "kw_results/"))
	resultTypes.append(ResultType("Emma", ".pvals", "newDataset", res_path + "emma_results/", mafCutoff=15))
	return resultTypes

def _getStandardResultTypes4_():
	res_path = "/Users/bjarnivilhjalmsson/Projects/Data/gwas_results/"
	resultTypes = []
	resultTypes.append(ResultType("KW", ".pvals", "raw", res_path + "kw_results/"))
	resultTypes.append(ResultType("Emma", ".pvals", "trans", res_path + "emma_results/", mafCutoff=0.1))
	return resultTypes


def _getStandardBinaryResultTypes_():
	res_path = "/Network/Data/250k/tmp-bvilhjal/"
	resultTypes = []
	resultTypes.append(ResultType("KW", ".pvals", "newDataset", res_path + "kw_results/"))
	resultTypes.append(ResultType("Marg", ".score", "newDataset", res_path + "marg_results/"))
	resultTypes.append(ResultType("RF", ".imp", "newDataset", res_path + "rf_results/", mafCutoff=15))
	return resultTypes

def _getStandardSecondRunResultTypes_():
	res_path = "/Network/Data/250k/tmp-bvilhjal/"
	resultTypes = []
	resultTypes.append(ResultType("KW", ".pvals", "raw", res_path + "kw_results/"))
	resultTypes.append(ResultType("KW", ".pvals", "raw", res_path + "kw_results/"))
#	resultTypes.append(ResultType("Emma",".pvals","logTransform",res_path+"emma_results/",mafCutoff=20))
#	resultTypes.append(ResultType("Emma",".pvals","logTransform",res_path+"emma_results/",mafCutoff=20))
#	resultTypes.append(ResultType("Emma",".pvals","raw",res_path+"emma_results/",mafCutoff=15))
#	resultTypes.append(ResultType("Emma",".pvals","raw",res_path+"emma_results/",mafCutoff=15))
	return resultTypes


def load_result_from_db(pid, aid, cmid=54, host='gmi-ara-devel-be', conn=None):
	"""
	Imports a result object from the filesystem/DB.
	"""
	import dbutils
	if conn:
		cursor = conn.cursor()
	else:
		new_conn = dbutils.connect_to_db(host, 'stock_250k')
		cursor = new_conn.cursor()
	sql_statement = "SELECT short_name, filename, id FROM stock_250k.results_method \
			 WHERE call_method_id=%d and phenotype_method_id=%d and analysis_method_id=%d"\
			 % (cmid, pid, aid)
	#print sql_statement
	numRows = int(cursor.execute(sql_statement))
	row = cursor.fetchone()
	r, r_id = None, None
	if row:
		fname = row[1]
		r_id = int(row[2])
		print "File for %s, found at:%s" % (row[0], fname)
		r = Result(fname)

	else:
		print "Result not found with pid=%d, aid=%d, cmid=%d" % (pid, aid, cmid)

	cursor.close ()
	if not conn:
		new_conn.close ()

	return r, r_id


def loadResults(phenotypeIndices, resultTypes=None, phed=None, snpsds=None, filterPercentile=None, filterCutoffs=None, phenotypeFile="/Network/Data/250k/dataFreeze_080608/phenotypes_all_raw_111008.tsv", secondRun=False):

	if not phed:
		phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')

	if not resultTypes:
		if secondRun:
			resultTypes = _getStandardSecondRunResultTypes_()
		else:
			resultTypes = _getStandardResultTypes4_()

	results_map = {}
	for i in phenotypeIndices:

		results = []
		for j in range(0, len(resultTypes)):
			resultType = resultTypes[j]
			phenName = phed.getPhenotypeName(i)
			if phenName:
				resultFile = resultType.getFileName(phed, i, secondRun=(secondRun and j % 2 == 1))  #Modify back to get old results 120708
				try:
					print "Loading result file", resultFile
					if snpsds:
						result = SNPResult(resultFile, snpsds=snpsds, name=str(resultType) + "_" + phenName, resultType=resultType, phenotypeID=i)
					else:
						result = Result(resultFile, name=str(resultType) + "_" + phenName, resultType=resultType, phenotypeID=i)
					if resultType.logTransform:
						print "Log transformed the p-values"
						result.negLogTransform()

					result.filterMARF(minMaf=resultType.mafCutoff)
					#result.filterMAF(minMaf=resultType.mafCutoff)
					if filterPercentile:
						result.filterPercentile(filterPercentile)
					elif filterCutoffs:
						result.filterScoreCutoff(filterCutoffs[j])

					results.append(result)
				except Exception, e:
					print e.message
					print "Couldn't load", resultFile

		results_map[i] = results
		gc.collect()  #Calling garbage collector, in an attempt to clean up memory..
	return results_map



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
	import matplotlib
	import matplotlib.pyplot as plt

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
	import matplotlib
	import matplotlib.pyplot as plt
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


def _add_results_to_db_():
	"""
	TEST
	"""
	host = 'gmi-ara-devel-be'
	for pid in range(1, 2):
		for aid in range(20):
			for cmid in range(60):
				try:
					r, r_id = load_result_from_db(pid, aid, cmid, host=host)
					if not r:
						raise Exception
					r.insert_into_db(r_id)
		                except Exception, err_str:
		                	print "File not found, pid=%d, aid=%d, cmid=%d, err_str=%s" % (pid, aid, cmid, err_str)






if __name__ == "__main__":
	#load_cand_genes_file("/Users/bjarnivilhjalmsson/Projects/Ales_Pecinka/UV_cand_genes_021710.csv")
	#load_result_from_db(513,1)
	_add_results_to_db_()
