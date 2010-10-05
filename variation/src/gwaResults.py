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


class Result(object):
	"""
	Contains information on the result.  (The old object is renamed Result_old, for now..)  Poised to cause problems?
	"""
	def __init__(self, result_file=None, scores=None, snps_data=None, accessions=None, name=None,
		     result_type=None, phen_id=None, positions=None, chromosomes=None, marfs=None,
		     mafs=None, snps=None, recall_snps=False, **snp_results_info):
		"""
		032610: A new init function to clean up the previous mess...  

		A simple init function.
		
		snps_data is assumed to match the results, (in term of positions, etc)
		
		accessions (when used) should fit the results!
		
		recall_snps: if this is flagged, then the result object stores all the SNPs, 
			otherwise SNP data is just uesd for MAF calculations. 
		"""

		#Contain various information for every snp and position
		#FIXME: To replace older lists.. such as positions, chromosomes, scores, etc. 
		self.snp_results = {}

		self.phen_id = phen_id
		self.result_type = result_type
		self.name = name
		self.positions = positions
		self.chromosomes = chromosomes
		self.marfs = marfs #Minor allele relative frequencies.
		self.mafs = mafs
		self.snps = snps
		self.accessions = accessions
		if scores != None:
			self.scores = list(scores) #Scores or p-values
		else:
			self.scores = []
		if not self.positions:
			self.positions = []
		if not self.chromosomes:
			self.chromosomes = []
		if not self.marfs:
			self.marfs = []
		if not self.mafs:
			self.mafs = []
		if not self.snps:
			self.snps = []
		if not self.accessions:
			self.accessions = []

		self.orders = None
		self.ranks = None
		self.chromosome_ends = []

		if result_file:
			self._load_result_(result_file)
		elif snps_data:
			self._load_snps_data_(snps_data)

		self.snp_results['chromosomes']	 = self.chromosomes
		self.snp_results['scores'] = self.scores
		self.snp_results['positions'] = self.positions
		self.snp_results['marfs'] = self.marfs
		self.snp_results['mafs'] = self.mafs
		self.snp_results['snps'] = self.snps
		if snp_results_info:
			for info in snp_results_info:
				self.snp_results[info] = snp_results_info[info]



	def _load_result_(self, resultFile):
		f = open(resultFile, "r")
		lines = f.readlines()
		print "Loading", len(lines), "lines of results."
		try_delims = [",", "\t", " ", ]
		i = 0
		delim = try_delims[i]
		while len(lines[0].split(delim)) == 1:
			i += 1
			delim = try_delims[i]
		if len(lines[0].split(delim)) == 1:
			raise Exception("Appropriate delimiter wasn't found.")
		#print "Delimiter used:",str(delim)	
		line = lines[0].split(delim)
		withMAF = len(line) > 3
		#print "withMAF:", withMAF
		start = 0
		currChrom = -1
		lastPos = 0
		extra_columns = []
		if not (line[0].strip()).isdigit():
			start = 1
			print "Header detected"
			columns = map(str.strip, line)
			print columns
			#extra_columns
		if not withMAF:
			for i in range(start, len(lines)):
				line = lines[i].split(delim)
				newChrom = int(line[0].strip())
				if newChrom != currChrom:
					currChrom = newChrom
					if lastPos:
						self.chromosome_ends.append(lastPos)
				self.chromosomes.append(newChrom)
				lastPos = int(line[1].strip())
				self.positions.append(lastPos)
				self.scores.append(float(line[2].strip()))
		else:
			for i in range(start, len(lines)):
				line = lines[i].split(delim)
				newChrom = int(line[0].strip())
				if newChrom != currChrom:
					currChrom = newChrom
					if lastPos:
						self.chromosome_ends.append(lastPos)
				self.chromosomes.append(newChrom)
				lastPos = int(line[1].strip())
				self.positions.append(lastPos)
				self.scores.append(float(line[2].strip()))
				self.marfs.append(float(line[3].strip()))
				self.mafs.append(float(line[4].strip()))

			assert len(self.positions) == len(self.mafs) == len(self.scores), \
			     "length of scores, positions, mafs, isn't equal"

		self.chromosome_ends.append(lastPos)
		self._rank_scores_()


	def _load_snps_data_(self, snps_data):
		self.chromosome_ends = []
		for i, snpsd in enumerate(snps_data.snpsDataList):
			self.chromosomes.extend([snps_data.chromosomes[i]] * len(snpsd.positions))
			self.positions.extend(snpsd.positions)
			self.snps.extend(snpsd.snps)
			self.chromosome_ends.append(snpsd.positions[-1])
		maf_d = snps_data.get_mafs()
		self.mafs = maf_d['mafs']
		self.marfs = maf_d['marfs']


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
		if self.snps:
			res_ls = zip(self.chromosomes, self.positions, self.scores, self.mafs, self.marfs, self.snps)
		else:
			res_ls = zip(self.chromosomes, self.positions, self.scores, self.mafs, self.marfs)
		res_ls.sort()
		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		new_snps = []
		for res in res_ls:
			newScores.append(res[0])
			newPositions.append(res[1])
			newChromosomes.append(res[2])
			newMafs.append(res[3])
			newMarfs.append(res[4])
			if self.snps:
				new_snps.append(res[5])

		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs
		self.snps = new_snps




	def plot_manhattan(self, pdf_file=None, png_file=None, min_score=None, max_score=None,
		       percentile=98, type="pvals", ylab="$-$log$_{10}(p-$value$)$",
		       plot_bonferroni=False, cand_genes=None):

		import matplotlib
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt
		"""
		Plots a 'Manhattan' style GWAs plot.
		"""

		"Plotting a Manhattan-style plot with %i markers." % len(self.scores)
		if cand_genes:
			#processing candidate genes by chromosome
			chr_cand_genes = []
			for ch in [1, 2, 3, 4, 5]:
				cgs = []
				for cg in cand_genes:
					if cg.chromosome == ch:
						cgs.append(cg)
				chr_cand_genes.append(cgs)

		num_scores = len(self.scores)
		result = self.simple_clone()

		result.filter_percentile(percentile / 100.0)
		if percentile < 100 and len(result.scores) == len(self.scores):
			raise Exception()

		if not max_score:
			max_score = max(result.scores)
		if not min_score:
			if type == "pvals":
				min_score = 0
			else:
				min_score = min(result.scores)

		scoreRange = max_score - min_score
		offset = 0
		chromosomeSplits = result.get_chromosome_splits()
		ticksList1 = []
		ticksList2 = []
		textPos = []
		plt.figure(figsize=(12, 2.8))
		plt.axes([0.045, 0.15, 0.95, 0.71])
		starPoints = [[], [], []]
		for i in range(0, len(result.chromosome_ends)):
			index1 = chromosomeSplits[i][0]
			index2 = chromosomeSplits[i + 1][0]
			scoreList = result.scores[index1:index2]
			posList = result.positions[index1:index2]

			if cand_genes:
				for cg in chr_cand_genes[i]:
					plt.axvspan(offset + cg.startPos, offset + cg.endPos, facecolor='k', alpha=0.5)


			newPosList = []
			for pos in posList:
				newPosList.append(offset + pos)

			for s_i, (score, pos) in enumerate(zip(scoreList, newPosList)):
				if score > max_score:
					starPoints[0].append(pos)
					starPoints[1].append(max_score)
					starPoints[2].append(score)
					score = max_score
				scoreList[s_i] = score

			plt.plot(newPosList, scoreList, ".", markersize=5, alpha=0.7)
			oldOffset = offset
			textPos.append(offset + result.chromosome_ends[i] / 2 - 2000000)
			offset += result.chromosome_ends[i]
			for j in range(oldOffset, offset, 5000000):
				ticksList1.append(j)
			for j in range(0, result.chromosome_ends[i], 5000000):
				if j % 10000000 == 0 and j < result.chromosome_ends[i] - 2500000 :
					ticksList2.append(j / 1000000)
				else:
					ticksList2.append("")



		plt.plot(starPoints[0], starPoints[1], ".", color="#ee9922", markersize=6)
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




		if plot_bonferroni:
			import math
			bonferroni_threshold = -math.log10(1.0 / (num_scores * 20.0))
			plt.plot([0, sum(result.chromosome_ends)], [bonferroni_threshold, bonferroni_threshold], "k-.")

		plt.axis([0, sum(result.chromosome_ends), min_score - 0.05 * scoreRange, max_score + 0.05 * scoreRange])
		plt.xticks(ticksList1, ticksList2)
		if not ylab:
			if type == "pvals":
				plt.ylabel('$-log(p-$value$)$', size="large")

			else:
				plt.ylabel('score')
		else:
			plt.ylabel(ylab)
		plt.xlabel("Mb", size="large")

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
		USE WITH CARE
		"""
		oldChrom = 0
		chromosome_splits = []
		for i in range(0, len(self.scores)):
			newChrom = self.chromosomes[i]
			if oldChrom != newChrom:
				while oldChrom < newChrom:
					oldChrom += 1
					chromosome_splits.append((i, oldChrom))
		chromosome_splits.append((i, -1))
		return chromosome_splits


	def neg_log_trans(self):
		"""
		apply - log(x) to the pvalues (scores)
		"""
		import math, warnings
		for i, score in enumerate(self.scores):
			if score != 0.0:
				self.scores[i] = -math.log(score, 10)
			else:
				warnings.warn("P-value is 0, log transformation is invalid.  Score is arbitrary set to 50.")
				self.scores[i] = 50



	def filter_percentile(self, percentile):
		new_scores = []
		for score in self.scores:
			new_scores.append(score)
		new_scores.sort()
		score_cutoff = new_scores[int(len(new_scores) * percentile)]
		self.filter_attr("scores", score_cutoff)



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
		new_snp_results = {}
		for info in self.snp_results:
			new_snp_results[info] = []
		count = len(self.scores)
		for i in indices_to_keep:
			for info in self.snp_results:
				if self.snp_results[info]:
					new_snp_results[info].append(self.snp_results[info][i])

		self.snp_results = new_snp_results
		self.scores = self.snp_results['scores']
		self.positions = self.snp_results['positions']
		self.chromosomes = self.snp_results['chromosomes']
		self.mafs = self.snp_results['mafs']
		self.marfs = self.snp_results['marfs']
		self.snps = self.snp_results['snps']
		print "%i scores were removed." % (count - len(self.scores))


	def filter_attr(self, attr_name, attr_threshold):
		"""
		Filter out scores / pvalues etc. which have attr < attr_threshold.

		attr are e.g.
		'mafs', 'marfs', 'scores', etc.
		"""
		print "Filtering for attribute '%s' with threshold: %g" % (attr_name, attr_threshold)
		attr = getattr(self, attr_name)
		print len(attr)
		new_snp_results = {}
		for info in self.snp_results:
			new_snp_results[info] = []
		count = len(self.scores)
		for i in range(0, len(self.scores)):
			if attr[i] >= attr_threshold:
				for info in self.snp_results:
					if len(self.snp_results[info]) > 0:
						new_snp_results[info].append(self.snp_results[info][i])

		self.snp_results = new_snp_results
		self.scores = self.snp_results['scores']
		self.positions = self.snp_results['positions']
		self.chromosomes = self.snp_results['chromosomes']
		self.mafs = self.snp_results['mafs']
		self.marfs = self.snp_results['marfs']
		self.snps = self.snp_results['snps']
		print "%i scores were removed." % (count - len(self.scores))
		return len(self.scores)



	def filter_non_segregating_snps(self, ecotype1, ecotype2, accessions=None):
		"""
		Filter out all SNPs which are not segregating in the two accessions.

		Assumes the accessions map the results objects SNPs. (and that they are defined)
		"""
		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		new_snps = []

		if accessions:
			ecotypes = accessions
		else:
			ecotypes = self.accessions

		e_i1 = ecotypes.index(ecotype1)
		e_i2 = ecotypes.index(ecotype2)

		for i in range(len(self.snps)):
			snp = self.snps[i]
			if snp[e_i1] != snp[e_i2]:
				newScores.append(self.scores[i])
				newPositions.append(self.positions[i])
				newChromosomes.append(self.chromosomes[i])
				newMafs.append(self.mafs[i])
				newMarfs.append(self.marfs[i])
				new_snps.append(snp)

		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs
		self.snps = new_snps



#	def filter_score_cutoff(self, score_cutoff):
#
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		new_snps = []
#
#		for i in range(0,len(self.scores)):
#			if self.scores[i]>score_cutoff:
#				newScores.append(self.scores[i])
#				newPositions.append(self.positions[i])
#				newChromosomes.append(self.chromosomes[i])
#				newMafs.append(self.mafs[i])
#				newMarfs.append(self.marfs[i])
#				if self.snps:
#					new_snps.append(self.snps[i])
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
			chr_pos_list.append((self.chromosomes[i], self.positions[i]))
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
			p = self.positions[snp_i]
			c = self.chromosomes[snp_i]
			c_genes = get_gene_list(p - window_size, p + window_size, c, include_intron_exons=False, conn=conn)
			for g in c_genes:
				genes.add((g.chromosome, g.startPos, g.endPos, g.tairID))
			#print len(genes)
			i += 1
		print 'found %d genes' % len(genes)
		return genes


	def get_chr_pos_score_list(self):
		return zip(self.chromosomes, self.positions, self.scores)


	def get_chr_pos_list(self):
		return zip(self.chromosomes, self.positions)


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
			snp_pos = self.positions[snp_i]
			chromosome = self.chromosomes[snp_i]
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


	def get_region_result(self, chromosome, start_pos, end_pos, buffer=0):
		"""
		returns a result object with only the SNPs, etc. within the given boundary.
		"""
		positions = []
		scores = []
		snps = []
		mafs = []
		marfs = []
		chromosomes = []
		i = 0
		start_pos -= buffer
		end_pos += buffer
		while i < len(self.chromosomes) and self.chromosomes[i] != chromosome:
			i += 1
		while i < len(self.chromosomes) and self.positions[i] < start_pos:
			i += 1
		if i == len(self.chromosomes):
			raise Exception("region results never found!")

		while i < len(self.chromosomes) and self.positions[i] <= end_pos and self.chromosomes[i] == chromosome:
			chromosomes.append(chromosome)
			positions.append(self.positions[i])
			scores.append(self.scores[i])
			if self.snps:
				snps.append(self.snps[i])
			mafs.append(self.mafs[i])
			marfs.append(self.marfs[i])
			i += 1

		return Result(result_type=self.result_type, phen_id=self.phen_id,
			      scores=scores, chromosomes=chromosomes, positions=positions,
			      snps=snps, accessions=self.accessions, marfs=marfs, mafs=mafs)



	def get_top_region_results(self, n, distance_threshold=25000, buffer=25000):
		reg_results = []
		i = 0
		for reg in self.get_top_regions(n, distance_threshold):
			reg_results.append(self.get_region_result(*reg, buffer=buffer))
		print "Regions were retrieved."
		return reg_results


	def filter_top_snps(self, n):
		self._rank_scores_() #Making sure the ranks are updated
		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		new_snps = []
		l = self.orders[0:n]
		l.sort()

		for i in l:
			newScores.append(self.scores[i])
			newPositions.append(self.positions[i])
			newChromosomes.append(self.chromosomes[i])
			newMafs.append(self.mafs[i])
			newMarfs.append(self.marfs[i])
			if self.snps:
				new_snps.append(self.snps[i])

		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs
		self.snps = new_snps



	def get_chromosome_ends(self):
		if not self.chromosome_ends:
			self._sort_by_chr_pos_()
			i = 0
			ch_i = 0
			curr_chr = self.chromosomes[i]
			chromosome_ends = []
			while ch_i < 5:
				while i < len(self.chromosomes) and self.chromosomes[i] == curr_chr:
					i += 1
				chromosome_ends.append(self.positions[i - 1])
			self.chromosome_ends = chromosome_ends
		return self.chromosome_ends



	def clone(self):
		import copy
		result = copy.deepcopy(self) #Cloning
		return result


	def simple_clone(self):
		result = Result(scores=self.scores[:], positions=self.positions[:], chromosomes=self.chromosomes[:],
				marfs=self.marfs[:], mafs=self.mafs[:], accessions=self.accessions[:])
		result.chromosome_ends = self.chromosome_ends[:]
		return result




#	def insert_into_db(self,result_id,host='gmi-ara-devel-be'):
#		"""
#		Insert pvalues into DB.. (TEST) 
#		"""
#		columns = ['result_id', 'chromosome', 'position', 'pValue']
#		if self.marfs:
#			columns.append('MAF')
#		if self.mafs:
#			columns.append('MAC')
#		import dbutils
#		sql_statement_prefix = "INSERT INTO stock_snp_test.result_pvalues ("+",".join(columns)+")"
#		c_p_pv_list = self.get_chr_pos_score_list()
#		print len(c_p_pv_list)
#		conn = dbutils.connect_to_db(host,'stock_snp_test')
#		cursor = conn.cursor()
#		for i, (c,p,pv) in enumerate(c_p_pv_list):
#			sql_statement = sql_statement_prefix+"VALUES (%d, %d, %d, %f"%(result_id,c,p,pv)
#			if self.marfs and self.mafs:
#				sql_statement +=", %f, %d)"%(self.marfs[i],self.mafs[i])
#			else:
#				sql_statement +=")"
#		
#			#print sql_statement
#			numRows = int(cursor.execute(sql_statement))
#			row = cursor.fetchone()
#			if row:
#				print row
#		conn.commit()
#		cursor.close ()
#		conn.close ()



	def na_mafs(self, min_maf=10):
		"""
		NA scores/pvalues which have maf<minMaf.		
		"""
		for i in range(0, len(self.scores)):
			if self.mafs[i] < min_maf:
				self.scores[i] = "NA"


	def get_max_snp(self):
		max_val = max(self.scores)
		mi = self.scores.index(max_val)
		return (self.snps[mi], self.scores[mi], self.chromosomes[mi], self.positions[mi])



	def write_to_file_old(self, filename, with_extra_info=False):
		header = ['chromosomes', 'positions', 'scores', 'marfs', 'mafs']
		if with_extra_info:
			for info in self.snp_results:
				if not info in header:
					header.append(info)
		f = open(filename, "w")

		f.write("Chromosome,Position,Score,MARF,MAF \n")

		for (ch, pos, score, marf, maf) in zip(self.chromosomes, self.positions, self.scores, self.marfs, self.mafs):
			l = map(str, [ch, pos, score, marf, maf])
			f.write(",".join(l) + "\n")
		f.close()


	def write_to_file(self, filename, additional_columns=None):
		columns = ['chromosomes', 'positions', 'scores', 'marfs', 'mafs']
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
#class Result_old(object):
#	"""
#	Contains information on the result.
#	"""
#
#	def __init__(self,resultFile=None,snpsds=None,chromosome_list=None,name=None,
#		     resultType=None,phenotypeID=None,interactionResult=False,load_snps=True,
#		     scores=[],snps=[],chromosomes=[],positions=[],accessions=None,
#		     marfs=[],mafs=[]):
#		"""
#		if snps is given as argument, then they must match the other data...
#		"""
#
#		self.phenotypeID =phenotypeID
#		self.resultType=resultType
#		self.name = name
#		self.scores = scores #Scores or p-values
#		self.positions = positions
#		self.chromosomes = chromosomes		
#		self.marfs = marfs #Minor allele relative frequencies.
#		self.mafs = mafs #Minor allele frequencies.
#		self.accessions = accessions
#
#		self.secondaryScores = []
#		self.secondaryPositions = []
#		self.secondaryChromosomes = []
#		self.orders = None
#		self.ranks = None
#		self.chromosomeEnds = []
#		self.snps = snps
#
#		self.interactionResult=interactionResult
#		if not chromosomes and snpsds:
#			chromosomes = range(1,len(snpsds)+1)
#		if interactionResult:
#			self.interactionPositions = []			
#		if resultFile:
#			self._loadResult_(resultFile)
#		elif scores and snpsds:
#			print "Loading p-values directly!"
#			self._load_pvals_(scores,snpsds,chromosome_list)
#		if snpsds and load_snps:
#			self._loadSnpsData_(snpsds,chromosome_list)
#
#
#
#	def _loadResult_(self,resultFile):
#		f = open(resultFile,"r")
#		lines = f.readlines()
#		try_delims = [",","\t"," ",]
#		i = 0
#		delim = try_delims[i]
#		while len(lines[0].split(delim))==1:
#			i += 1
#			delim = try_delims[i]
#		if len(lines[0].split(delim))==1:
#			raise Exception("Appropriate delimiter wasn't found.")
#		print "Delimiter used:",str(delim)	
#		line = lines[0].split(delim)
#		withMAF = len(line)>3
#		print "withMAF:", withMAF
#		start = 0
#		currChrom = -1
#		lastPos = 0
#		if not (line[0].strip()).isdigit():
#			start = 1
#			#print "Header detected"
#			#print line[0].strip()
#		if not withMAF:
#			for i in range(start,len(lines)):
#				line = lines[i].split(delim)
#				newChrom = int(line[0].strip())
#				if newChrom!=currChrom:
#					currChrom = newChrom
#					if lastPos:
#						self.chromosomeEnds.append(lastPos)
#				self.chromosomes.append(newChrom)
#				if self.interactionResult:
#					iPos = [int(line[-1].strip())]
#					self.interactionPositions.append(iPos)
#				lastPos = int(line[1].strip())
#				self.positions.append(lastPos)
#				self.scores.append(float(line[2].strip()))
#		else:
#			for i in range(start,len(lines)):
#				line = lines[i].split(delim)
#				newChrom = int(line[0].strip())
#				if newChrom!=currChrom:
#					currChrom = newChrom
#					if lastPos:
#						self.chromosomeEnds.append(lastPos)
#				self.chromosomes.append(newChrom)
#				if self.interactionResult:
#					iPos = [int(line[-1].strip())]
#					self.interactionPositions.append(iPos)
#				lastPos = int(line[1].strip())
#				self.positions.append(lastPos)
#				self.scores.append(float(line[2].strip()))
#				self.marfs.append(float(line[3].strip()))
#				self.mafs.append(float(line[4].strip()))
#								
#		self.chromosomeEnds.append(lastPos)
#		self._calcGlobalRanks_()
#		if not (len(self.positions)==len(self.mafs)==len(self.scores)):
#			print len(self.positions),len(self.mafs),len(self.scores)
#			raise Exception
#
#	
#
#	def _load_pvals_(self,result_pvals,snpsds,chromosomes):
#		print "Loading the pvalues, and snpsdata's to construct a Result object."
#		for i, snpsd in enumerate(snpsds):
#			r = snpsd.get_mafs()
#			self.mafs += r["mafs"]
#			self.marfs += r["marfs"]
#			self.positions += snpsd.positions
#			self.chromosomes += [chromosomes[i]]*len(snpsd.positions)
#			self.chromosomeEnds.append(snpsd.positions[-1])
#		self.scores = result_pvals  # Should I log transform? 
#
#	def _calcGlobalRanks_(self):
#		"""
#		This function should only be called right after loading the full results.
#		"""
#		self._rankScores_()
#		self.globalRanks = self.ranks
#		self.ranks = None
#
#
#	def _loadSnpsData_(self,snpsds,chromosomes=None):
#		"""
#		Loads the SNP data.  (only loads those which are in the results)
#		"""		
#		self.snps = []
#		self.accessions=snpsds[0].accessions
#		chr_i = 0 # chromosome index
#		chr_set = set(self.chromosomes)
#		print chr_set
#		print len(self.scores)
#		for chr_i in range(len(chromosomes)):
#			if chromosomes[chr_i] in chr_set:
#				i = 0 #result index
#				while i<len(self.scores):
#					while i<len(self.chromosomes) and chromosomes[chr_i] != self.chromosomes[i]:
#						i += 1
#					j = 0
#					pos = self.positions[i]
#					while j<len(snpsds[chr_i].positions) and pos > snpsds[chr_i].positions[j]:
#						j += 1
#					if j<len(snpsds[chr_i].positions) and pos == snpsds[chr_i].positions[j]:
#						self.snps.append(snpsds[chr_i].snps[j])
#					i += 1
#		
#				if pos > snpsds[chr_i].positions[j]:
#					while j<len(snpsds[chr_i].positions) and pos > snpsds[chr_i].positions[j]:
#						j += 1
#					if j<len(snpsds[chr_i].positions) and pos == snpsds[chr_i].positions[j]:
#						self.snps.append(snpsds[chr_i].snps[j])
#
#		if i!= len(self.snps):
#			print "Problems with loading SNPs",i,len( self.snps)
#		else:
#			print "Loading SNPs appears to have worked?"
#			
#		print "Loaded",len(self.snps)," SNPs and",len(self.accessions),"accessions."
#
#	def _rankScores_(self):
#		"""
#		Generates two data structures: 
#		self.orders (SNP indices ordered by rank)
#		self.ranks (ranks of the SNPs)
#		"""
#		rank_ls = zip(self.scores, range(0,len(self.scores)))
#		rank_ls.sort()
#		rank_ls.reverse()
#		self.orders = []
#		for j in range(0,len(rank_ls)):
#			(s,i) = rank_ls[j]
#			self.orders.append(i)
#			rank_ls[j] = (i,j)
#
#		rank_ls.sort()
#		self.ranks = []
#		for (i,j) in rank_ls:
#			self.ranks.append(j+1)
#
#
#	def getChromosomeSplit(self):
#		"""
#		"""
#		oldChrom = 0
#		chromosomeSplits = []
#		for i in range(0,len(self.scores)):
#			newChrom = self.chromosomes[i]
#			if oldChrom != newChrom:
#				while oldChrom < newChrom:
#					oldChrom += 1
#					chromosomeSplits.append((i,oldChrom))
#		chromosomeSplits.append((i,-1))
#		return chromosomeSplits
#
#			
#	def getSecondaryChromosomeSplit(self):
#		"""
#		"""
#		oldChrom = 0
#		chromosomeSplits = []
#		for i in range(0,len(self.secondaryScores)):
#			newChrom = self.secondaryChromosomes[i]
#			if oldChrom != newChrom:
#				while oldChrom < newChrom:
#					oldChrom += 1
#					chromosomeSplits.append((i,oldChrom))
#		chromosomeSplits.append((i,-1))
#		return chromosomeSplits
# 
#
#	def negLogTransform(self):
#		"""
#		apply -log(x) to the pvalues (scores)
#		"""
#		import math
#		newScores = []
#		for score in self.scores:
#			if score != 0.0:
#				newScore = -math.log(score,10)
#				newScores.append(newScore)
#			else:
#				newScores.append(50)  				
#		self.scores = newScores
#		
#
#	def filterNicePeaks(self,scoreThreshold,singletonScoreThreshold,window=[20000,20000], method=1):
#		currScoreWeight=0.2
#		currChrom = -1
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		singletonCount = 0
#		lastSecondaryEnd = 0
#		if method == 1:
#			for i in range(0,len(self.positions)):
#				if currChrom != self.chromosomes[i]: #Restart
#					currChrom = self.chromosomes[i]
#					startIndex = i #windowIndices 
#					stopIndex = i #windowIndices 
#					curScoreSum = currScoreWeight*self.scores[i]				
#					oldScore=0
#					numSNPs = 1
#				currPos = self.positions[i]
#
#				while currPos - self.positions[startIndex]>window[0]:
#					curScoreSum -= self.scores[startIndex]
#					startIndex += 1
#					numSNPs -= 1
#				
#				while stopIndex+1 < len(self.positions) and self.positions[stopIndex+1]-currPos<window[1]:
#					stopIndex += 1
#					curScoreSum += self.scores[stopIndex]
#					numSNPs += 1
#
#				curScoreSum -= oldScore				
#				oldScore = currScoreWeight*numSNPs*self.scores[i]
#				curScoreSum += oldScore
#
#				if (numSNPs < 5 and self.scores[i] > singletonScoreThreshold) or  (numSNPs > 5 and curScoreSum/float(numSNPs) > scoreThreshold):
#					if (numSNPs < 5 and self.scores[i] > singletonScoreThreshold):
#						singletonCount +=1 
#					newScores.append(self.scores[i])
#					newPositions.append(self.positions[i])
#					newChromosomes.append(self.chromosomes[i])
#					newMafs.append(self.mafs[i])
#					newMarfs.append(self.mafs[i])
#
#		elif method==2:
#			for i in range(0,len(self.scores)):
#				if self.scores[i]>=singletonScoreThreshold:
#					newScores.append(self.scores[i])
#					newPositions.append(self.positions[i])
#					newChromosomes.append(self.chromosomes[i])
#					newMafs.append(self.mafs[i])
#					newMarfs.append(self.mafs[i])
#
#		# The following code locates the regions before and after the "nice" SNPs.
#		j = 0
#		for i in range(0,len(self.positions)):
#			pos = self.positions[i]
#			chr = self.chromosomes[i]
#			if j < len(newPositions) and pos == newPositions[j] and chr==newChromosomes[j]:
#				k = 0				
#				while  i+k > lastSecondaryEnd and pos-self.positions[i+k-1]<window[0] and self.chromosomes[i+k-1]==chr:
#					k -= 1
#				while i+k < len(self.positions)-1 and self.positions[i+k]-pos < window[1] and self.chromosomes[i+k]==chr:
#					if i+k > lastSecondaryEnd:
#						self.secondaryScores.append(self.scores[i+k])
#						self.secondaryPositions.append(self.positions[i+k])
#						self.secondaryChromosomes.append(self.chromosomes[i+k])						
#					k+= 1 
#				lastSecondaryEnd = i+k-1
#				j += 1 
#				
#		
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#
#		print singletonCount,"singletons were added"
#
#
#	def filterPercentile(self,percentile):
#		newScores = []
#		for score in self.scores:
#			newScores.append(score)
#		newScores.sort() 
#		scoreCutoff = newScores[int(len(newScores)*percentile)]
#		self.filterScoreCutoff(scoreCutoff)
#
#
#	def filter(self, quantile=0.98, window=[25000,25000], singletonScoreCutoff=None, nicePeaks = False, method=1):
#		"""
#		Filter out scores/pvalues.
#				
#		"""
#		originalSize = len(self.scores)
#		newScores = []
#		for score in self.scores:
#			newScores.append(score)
#		newScores.sort() 
#		
#		if nicePeaks: 
#			basic_quantile = 0.90
#			top_quantile = 0.998
#			singleton_quantile=0.9998
#			
#			scoreCutoff = newScores[int(len(newScores)*basic_quantile)]
#			self._filterScoreCutoff_(scoreCutoff)
#			scoreCutoff = newScores[int(len(newScores)*top_quantile)]				
#			if not singletonScoreCutoff:
#				singletonScoreCutoff = newScores[int(len(newScores)*singleton_quantile)]
#			self.filterNicePeaks(scoreCutoff,singletonScoreCutoff,window,method)
#
#			
#		scoreCutoff = newScores[int(len(newScores)*quantile)]
#		self._filterScoreCutoff_(scoreCutoff)
#
#		finalSize = len(self.scores)
#
#		print "Original results size =",originalSize
#		print "results size after filtration =",finalSize
#
#
#	def filterMAF(self, minMaf=15):
#		"""
#		Filter out scores/pvalues which have maf<minMaf.		
#		"""
#		if self.interactionResult:
#			newIPos = []
#			for i in range(0,len(self.scores)):
#				if self.mafs[i]>= minMaf:
#					newIPos.append(self.interactionPositions[i])
#			del self.interactionPositions
#			self.interactionPositions = newIPos
#
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		for i in range(0,len(self.scores)):
#			if self.mafs[i]>= minMaf:
#				newScores.append(self.scores[i])
#				newPositions.append(self.positions[i])
#				newChromosomes.append(self.chromosomes[i])
#				newMafs.append(self.mafs[i])
#				newMarfs.append(self.marfs[i])
#
#		del self.scores
#		del self.positions
#		del self.chromosomes
#		del self.mafs
#		del self.marfs
#		
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#
#
#	def filterMARF(self, minMaf=0.1):
#		"""
#		Filter out scores/pvalues which have maf<minMaf.		
#		"""
#		if self.interactionResult:
#			newIPos = []
#			for i in range(0,len(self.scores)):
#				if self.mafs[i]>= minMaf:
#					newIPos.append(self.interactionPositions[i])
#			del self.interactionPositions
#			self.interactionPositions = newIPos
#
#		snps_indices_to_keep = []
#
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		for i in range(0,len(self.scores)):
#			if self.marfs[i]>= minMaf:
#				newScores.append(self.scores[i])
#				newPositions.append(self.positions[i])
#				newChromosomes.append(self.chromosomes[i])
#				newMafs.append(self.mafs[i])
#				newMarfs.append(self.marfs[i])
#				snps_indices_to_keep.append(i)
#
#		del self.scores
#		del self.positions
#		del self.chromosomes
#		del self.mafs
#		del self.marfs
#		
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#		return snps_indices_to_keep
#
#	
#	def filterNonSegregatingSnps(self, ecotype1, ecotype2, snpsd):
#		"""
#		Filter out all SNPs which are not segregating in the two accessions.		
#		"""
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#
#		ecotypes = snpsd.accessions
#		
#		e_i1 = ecotypes.index(ecotype1)
#		e_i2 = ecotypes.index(ecotype2)
#
#		res_chr_pos_list = zip(self.chromosomes,self.positions)
#		snpsd_chr_pos_snp_list = snpsd.getChrPosSNPList()
#
#		i = 0 #index in snpsd
#		j = 0 #index in result
#		
#		while i < len(snpsd_chr_pos_snp_list):
#			(chr,pos,snp) = snpsd_chr_pos_snp_list[i]
#			if j<len(res_chr_pos_list) and (chr,pos)==res_chr_pos_list[j]:
#				if snp[e_i1]!=snp[e_i2]:
#					newScores.append(self.scores[j])
#					newPositions.append(self.positions[j])
#					newChromosomes.append(self.chromosomes[j])
#					newMafs.append(self.mafs[j])
#					newMarfs.append(self.marfs[j])					
#				j += 1				
#			elif j<len(res_chr_pos_list) and (chr,pos)>res_chr_pos_list[j]:
#				import pdb;pdb.set_trace()
#				
#				print "ERROR!!!"
#			i += 1
#			
#		n = len(self.scores)
#		del self.scores
#		del self.positions
#		del self.chromosomes
#		del self.mafs
#		del self.marfs
#
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#		
#		print n-len(newScores),"SNPs were removed, with",len(newScores),"remaining."
#
#			
#	def filterScoreCutoff(self, scoreCutoff):
#
#		if self.interactionResult:
#			newIPos = []
#			for i in range(0,len(self.scores)):
#				if self.scores[i]>scoreCutoff:
#					newIPos.append(self.interactionPositions[i])
#			del self.interactionPositions
#			self.interactionPositions = newIPos
#
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		for i in range(0,len(self.scores)):
#			if self.scores[i]>scoreCutoff:
#				newScores.append(self.scores[i])
#				newPositions.append(self.positions[i])
#				newChromosomes.append(self.chromosomes[i])
#				newMafs.append(self.mafs[i])
#				newMarfs.append(self.marfs[i])
#		del self.scores
#		del self.positions
#		del self.chromosomes
#		del self.mafs
#		del self.marfs
#
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#		
#
#	def getRegions(self,window=[25000,25000]):
#		self._rankScores_()
#		snpIndex = 0
#		startPos = startPos = max(0,self.positions[0]-window[0])
#		endPos = self.positions[0]+window[1]
#		regions = []
#		snpRank = self.ranks[0]
#		chr = self.chromosomes[0]
#		maxScore = self.scores[0]
#		maxPos = self.positions[0]
#		oldPos = self.positions[0]
#		for i in range(1,len(self.positions)):
#			pos = self.positions[i]
#			if chr!=self.chromosomes[i]:
#				size = endPos-startPos
#				regions.append((snpRank,snpIndex,startPos,endPos,chr,size,maxScore,maxPos))
#				snpIndex = -1
#				maxScore = 0
#				startPos = max(0,pos-window[0])
#				chr = self.chromosomes[i]
#			elif pos-oldPos>sum(window):
#				size = endPos-startPos
#				regions.append((snpRank,snpIndex,startPos,endPos,chr,size,maxScore,maxPos))
#				maxScore = 0
#				startPos = pos-window[0]
#			if self.scores[i]>maxScore:
#				snpIndex = i 
#				snpRank = self.ranks[snpIndex]
#				maxScore = self.scores[i]
#				maxPos = self.positions[i]
#			
#			endPos = pos+window[1]
#			oldPos = pos
#		
#		regions.append((snpRank,snpIndex,startPos,endPos,chr,size,maxScore,maxPos))
#		
#		regions.sort()
#		
#		self.regions = []
#		for i in range(0,len(regions)):
#			(snpRank,snpIndex,startPos,endPos,chr,size,maxScore,maxPos) = regions[i]
#			#self.regions.append((snpIndex,startPos,endPos,chr,size,maxScore,maxPos,snpRank,i+1))
#			#def __init__(self,chromosome,startPos,endPos,snps=None,snpsd_indices=None,maxScores=None,maxPositions=None,maxSnpIndices=None,ranks=None):
#			self.regions.append(Region(chr,startPos,endPos))#,maxScores=[maxScore],maxPositions=[maxPos],maxSnpIndices=[snpIndex],ranks=[snpRank]))
#		
#		#pdb.set_trace()
#		self.regions.sort()
#		#print self.regions		
#		return self.regions
#		
#
#	def getTopSnps(self,n):
#		"""
#		returns top n SNPs
#		"""
#		import copy
#		result = copy.deepcopy(self) #Cloning
#		result.filterTopSNPs(n) 
#		return result
#	
#	def get_snps(self):
#		"""
#		Returns this result's SNPs as a list of SNP objects.
#		"""
#		raise  NotImplementedError
#		for i, snp in enumerate(self.snps):
#			SNP(position,chromosome,accessions=None,alleles=None,snpsds=None,score=None,rank=None,regionRank=None)
#			
#		
#			
#
#	def _countRegions_(self,res_ls,window=[25000,25000]):
#		oldPos = 0
#		countRegions = 0
#		chr = -1
#		for i in range(0,len(res_ls)):
#			pos = res_ls[i][1]
#			if chr!=res_ls[i][0]:
#				countRegions += 1
#				chr = res_ls[i][0]
#			elif pos-oldPos>sum(window):
#				countRegions += 1
#			oldPos = pos
#		
#		return countRegions
#
#
#	def getChrScorePos(self,chromosome):
#		"""
#		returns a list of (score,pos) tuples.
#		"""
#		i = 0
#		while i<len(self.chromosomes) and self.chromosomes[i]<chromosome:
#			i += 1
#		
#		posList = []
#		scoreList = []
#		while i<len(self.chromosomes) and self.chromosomes[i]==chromosome:
#			posList.append(self.positions[i])
#			scoreList.append(self.scores[i])
#			i += 1
#
#		return zip(scoreList,posList)
#
#	def getChrPos(self):
#		"""
#		returns a list of (chr,pos) tuples
#		"""
#		posList = []
#		chrList = []
#		for i in range(0,len(self.positions)):
#			posList.append(self.positions[i])
#			chrList.append(self.chromosomes[i])
#		return zip(chrList,posList)
#
#
#
#	def get_top_regions(self,n,distance_threshold=25000):
#		"""
#		returns a regionSet top n SNPs
#		"""
#		self._rankScores_() 
#		chromosome_ends = self.get_chromosome_ends()
#			
#		i = 0 #how many SNPs are needed
#		region_count = 0
#		regions = {1:[],2:[],3:[],4:[],5:[]}  #a list of (chromosome,start,stop)
#		while region_count < n:
#			snp_i = self.orders[i]
#			snp_pos = self.positions[snp_i]
#			chromosome = self.chromosomes[snp_i]
#			new_snp_reg = (chromosome,max(snp_pos-distance_threshold,0),min(snp_pos+distance_threshold,chromosome_ends[chromosome-1]))
#			merged=False
#			for reg_i,reg in enumerate(regions[chromosome]):
#				if (new_snp_reg[1]<reg[2] and new_snp_reg[2]>reg[2]) or (new_snp_reg[1]<reg[1] and new_snp_reg[2]>reg[1]): #then merge					
#					merged_reg = (reg[0],min(reg[1],new_snp_reg[1]),max(reg[2],new_snp_reg[2]))
#					print "Merged region:",merged_reg
#					regions[chromosome][reg_i] = merged_reg     
#					merged=True
#					break
#				elif (new_snp_reg[1]>reg[2] and new_snp_reg[2]<reg[2]):
#					merged=True  #It's already covered
#					break
#			if not merged:
#				regions[chromosome].append(new_snp_reg)
#				region_count+=1
#			i += 1 
#		
#		regions_list = regions[1]+regions[2]+regions[3]+regions[4]+regions[5]
#		return regions_list 
#	
#	def get_region_result(self,chromosome,start_pos,end_pos):
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
#		while i<len(self.chromosomes) and self.chromosomes[i]!=chromosome: 
#			i += 1
#		while i<len(self.chromosomes) and self.positions[i]<start_pos:
#			i += 1
#		if i == len(self.chromosomes):
#			raise Exception("region results never found!")
#		
#		while i<len(self.chromosomes) and self.positions[i]<=end_pos and self.chromosomes[i]==chromosome:
#			chromosomes.append(chromosome)
#			positions.append(self.positions[i])
#			scores.append(self.scores[i])
#			snps.append(self.snps[i])
#			mafs.append(self.mafs[i])
#			marfs.append(self.marfs[i])	
#		
#		return Result(resultType=self.resultType,phenotypeID=self.phenotypeID,
#			      scores=scores,chromosomes=chromosomes,positions=positions,
#			      snps=snps,accessions=self.accessions,marfs=marfs,mafs=mafs)
#		
#		
#		
#	def get_top_region_results(self,n,distance_threshold=25000):
#		reg_results = []
#		i = 0
#		for reg in self.get_top_regions(n, distance_threshold):
#			reg_results.append(self.get_region_result(*reg))
#		return reg_results
#	
#
#	def filterTopSNPs(self,n):
#		self._rankScores_() #Making sure the ranks are updated
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		l = self.orders[0:n]
#		l.sort()
#	
#		for i in l:
#			newScores.append(self.scores[i])
#			newPositions.append(self.positions[i])
#			newChromosomes.append(self.chromosomes[i])
#			newMafs.append(self.mafs[i])
#			newMarfs.append(self.marfs[i])
#
#		del self.scores
#		del self.positions
#		del self.chromosomes
#		del self.mafs
#		del self.marfs
#
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#
#
#	def _sortByChrPos_(self):
#		res_ls = zip(self.chromosomes,self.positions,self.scores,self.mafs,self.marfs)
#		res_ls.sort()
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		for (chr,pos,score,maf,marf) in res_ls:
#			newScores.append(score)
#			newPositions.append(pos)
#			newChromosomes.append(chr)
#			newMafs.append(maf)
#			newMarfs.append(marf)
#
#		del self.scores
#		del self.positions
#		del self.chromosomes
#		del self.mafs
#		del self.marfs
#
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#
#
#	def get_chromosome_ends(self):
#		if not self.chromosomeEnds:
#			self._sortByChrPos_()
#			i = 0
#			ch_i = 0
#			curr_chr = self.chromosomes[i]
#			chromosome_ends = [] 
#			while ch_i<5:
#				while i<len(self.chromosomes) and self.chromosomes[i]==curr_chr:
#					i += 1
#				chromosome_ends.append(self.positions[i-1])
#			self.chromosomeEnds = chromosome_ends
#		return self.chromosomeEnds
#		
#
#	def filterTopRegions(self,n,window=[25000,25000],minScore=None):
#		self._rankScores_() 
#		oldScores = self.scores
#		oldPositions = self.positions
#		oldChromosomes = self.chromosomes
#		oldMafs = self.mafs
#		oldMarfs = self.marfs
#		oldGlobalRanks = self.globalRanks
#
#		self.scores = []
#		self.positions = []
#		self.chromosomes = []
#		self.mafs = []
#		self.marfs = []
#		elf.globalRanks = []
#		regionCount = 0 
#		i = 0 
#		
#		res_ls = []
#		while regionCount < n:
#			res_ls.append((oldChromosomes[self.orders[i]],oldPositions[self.orders[i]]))
#			res_ls.sort()
#			if minScore and oldScores[self.orders[i]]<minScore:
#				break
#			self.scores.append(oldScores[self.orders[i]])
#			self.positions.append(oldPositions[self.orders[i]])
#			self.chromosomes.append(oldChromosomes[self.orders[i]])
#			self.mafs.append(oldMafs[self.orders[i]])
#			self.marfs.append(oldMarfs[self.orders[i]])
#			self.globalRanks.append(oldGlobalRanks[self.orders[i]])
#			regionCount = self._countRegions_(res_ls,window=window)
#			i += 1 
#		
#		del oldScores
#		del oldPositions
#		del oldChromosomes
#		del oldMafs
#		del oldMarfs
#		del oldGlobalRanks
#
#		self._sortByChrPos_()
#
#
#		
#			
#
#	def mergeWith(self,snpResult):
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		
#		i1 = 0
#		
#		pos1 = self.positions[i1] 
#		chr1 = self.chromosomes[i1]
#
#		for i2 in range(0,len(snpResult.positions)):
#			pos2 = snpResult.positions[i2]
#			chr2 = snpResult.chromosomes[i2]
#			while i1 <len(self.positions)-1 and (chr1,pos1)<(chr2,pos2):
#				newPositions.append(self.positions[i1])
#				newChromosomes.append(self.chromosomes[i1])
#				newScores.append(self.scores[i1])
#				newMafs.append(self.mafs[i1])
#				newMarfs.append(self.marfs[i1])
#				i1 += 1
#				pos1 = self.positions[i1]
#				chr1 = self.chromosomes[i1]
#
#			if i1 <len(self.positions)-1 and (chr1,pos1)==(chr2,pos2):
#				i1 += 1
#				pos1 = self.positions[i1]
#				chr1 = self.chromosomes[i1]
#
#			newPositions.append(snpResult.positions[i2])
#			newChromosomes.append(snpResult.chromosomes[i2])
#			newScores.append(snpResult.scores[i2])
#			newMafs.append(snpResult.mafs[i2])
#			newMarfs.append(snpResult.marfs[i2])
#
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#
#	def clone(self):
#		if self.snps:
#			result = SNPResult(name=self.name)
#			result.snps = self.snps + []
#		else:
#			result = Result(name=self.name)
#		result.resultType = self.resultType
#		result.name = self.name
#		result.scores = self.scores + []
#		result.positions = self.positions + []
#		result.chromosomes = self.chromosomes + []
#		result.chromosomeEnds = self.chromosomeEnds + []
#		result.marfs = self.marfs + [] #Minor allele relative frequencies.
#		result.mafs = self.mafs + [] #Minor allele frequencies.
#		if self.accessions:
#			result.accessions = self.accessions + []
#		if self.orders:
#		   result.orders = self.orders + []
#		if self.ranks:
#		   result.ranks = self.ranks + []
#		return result
#
#	
#	def naMAF(self, minMaf=10):
#		"""
#		NA scores/pvalues which have maf<minMaf.		
#		"""
#		for i in range(0,len(self.scores)):
#			if self.mafs[i]< minMaf:
#				self.scores[i]= "NA"
#
#
#	def alexFiltering(self,emmaScores,cutoff=6,window=[50000,50000]): 
#		"""
#		
#		"""
#		scoreList = zip(emmaScores,self.scores)
#		
#		scoreList.sort()
#		scoreList.reverse()
#		
#		cScores = []
#		i=0
#		emmaScore = scoreList[0][0]
#		while emmaScore>cutoff: 
#			cScores.append(scoreList[i][1])
#			i += 1
#			emmaScore = scoreList[i][0]
#
#
#		#Always top 5 regions, and at most 20 regions,  continue until cutoff is reached.
#
#
#		if len(cScores):
#			first5Regions = self.clone()
#			first5Regions.filterTopRegions(5,window=window)
#			cScoreCutoff = min(cScores)
#			self.filterTopRegions(20,window=window)
#			self.filterScoreCutoff(cScoreCutoff)
#			self.mergeWith(first5Regions)
#		else:
#			self.filterTopRegions(5,window=window)
#			
#	def updateRegions(self,regionList):
#		self._rankScores_()
#		i=0
#		rl_i=0 #region list index
#		while rl_i<len(regionList):
#			region = regionList[rl_i]
#			cp1=(self.chromosomes[i],self.positions[i])
#			cp_start=(region.chromosome,region.startPos)
#			while cp1 < cp_start:
#				i += 1
#				cp1=(self.chromosomes[i],self.positions[i])
#			cp_end=(region.chromosome,region.endPos)
#
#			maxScore = 0
#			maxRank = 0
#			maxPos = 0 
#			while cp1<=cp_end:
#				"""Update current region!"""
#				if maxScore < self.scores[i]:
#					maxScore = self.scores[i]
#					maxRank = self.ranks[i]
#					maxPos = self.positions[i]
#				i += 1
#				cp1=(self.chromosomes[i],self.positions[i])
#						
#			region.snpsInfo[self.resultType.name] = {"maxScore":maxScore,"maxRank":maxRank,"maxPos":maxPos}
#			rl_i += 1
#			
#	
#	def writeToFile(self,filename,format="simple"):
#		f = open(filename,"w")
#		f.write("Chromosome,Position,Score,MARF,MAF \n")
#		for (ch,pos,score,marf,maf) in zip(self.chromosomes,self.positions,self.scores,self.marfs,self.mafs):
#			l = map(str,[ch,pos,score,marf,maf])
#			f.write(",".join(l)+"\n")
#		f.close()
#		
#		
#		
#class EmmaResult(Result):
#	"""
#	Loads information on variance estimates and maximum likelihood as well...
#	"""			
#	def __init__(self,resultFile=None,snpsds=None,name=None,resultType=None,phenotypeID=None,interactionResult=False):
#
#		self.phenotypeID =phenotypeID
#		self.resultType=resultType
#		self.name = name
#		self.scores = [] #Scores or p-values
#		self.positions = []
#		self.chromosomes = []
#		self.chromosomeEnds = []
#		self.marfs = [] #Minor allele relative frequencies.
#		self.mafs = [] #Minor allele frequencies.
#		self.var_expl = [] #Variance explained
#		self.ml = [] #Maximum likelihood
#		self.vg = [] #Genetic variance scalar
#		self.ve = [] #Error variance scalar
#		self.accessions = None
#
#		self.secondaryScores = []
#		self.secondaryPositions = []
#		self.secondaryChromosomes = []
#		self.orders = None
#		self.ranks = None
#		self.snps = None
#		self.interactionResult=interactionResult
#		if interactionResult:
#			self.interactionPositions = []			
#		if resultFile:
#			self._loadResult_(resultFile)
#		if snpsds:
#			self._loadSnpsData_(snpsds)
#
#	
#	def _loadResult_(self,resultFile):
#		f = open(resultFile,"r")
#		lines = f.readlines()
#		line = lines[0].split(",")
#		withMAF = len(line)>3
#		start = 0
#		currChrom = -1
#		lastPos = 0
#		if not (line[0].strip()).isdigit():
#			start = 1
#
#		for i in range(start,len(lines)):
#			line = lines[i].split(",")
#			newChrom = int(line[0].strip())
#			if newChrom!=currChrom:
#				currChrom = newChrom
#				if lastPos:
#					self.chromosomeEnds.append(lastPos)
#			self.chromosomes.append(newChrom)
#			if self.interactionResult:
#				iPos = [int(line[-1].strip())]
#				self.interactionPositions.append(iPos)
#			lastPos = int(line[1].strip())
#			self.positions.append(lastPos)
#			self.scores.append(float(line[2].strip()))
#			self.marfs.append(float(line[3].strip()))
#			self.mafs.append(float(line[4].strip()))
#			self.var_expl.append(float(line[5].strip()))
#			self.ml.append(float(line[6].strip()))
#			self.vg.append(float(line[7].strip()))
#			self.ve.append(float(line[8].strip()))
#								
#		self.chromosomeEnds.append(lastPos)
#		self._calcGlobalRanks_()
#
#	def filterMARF(self, minMaf=0.1):
#		"""
#		Filter out scores/pvalues which have maf<minMaf.		
#		"""
#		if self.interactionResult:
#			newIPos = []
#			for i in range(0,len(self.scores)):
#				if self.mafs[i]>= minMaf:
#					newIPos.append(self.interactionPositions[i])
#			del self.interactionPositions
#			self.interactionPositions = newIPos
#
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		newVarExpl = []
#		newML = []
#		newVG = []
#		newVE = []
#		for i in range(0,len(self.scores)):
#			if self.marfs[i]>= minMaf:
#				newScores.append(self.scores[i])
#				newPositions.append(self.positions[i])
#				newChromosomes.append(self.chromosomes[i])
#				newMafs.append(self.mafs[i])
#				newMarfs.append(self.marfs[i])
#				newVarExpl.append(self.var_expl[i])
#				newML.append(self.ml[i])
#				newVG.append(self.vg[i])
#				newVE.append(self.ve[i])
#
#		del self.scores
#		del self.positions
#		del self.chromosomes
#		del self.mafs
#		del self.marfs
#		del self.var_expl
#		del self.ml
#		del self.vg
#		del self.ve
#		
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#		self.var_expl = newVarExpl
#		self.ml = newML
#		self.vg = newVG
#		self.ve = newVE


class SNPResult(Result):
	"""
	Contains information on the result.
	"""

	def _loadSnpsData_(self, snpsds):
		"""
		Loads the SNP data.
		"""
		self.snps = []
		self.accessions = snpsds[0].accessions
		i = 0 #result index
		chr = -1 # chromosome (index)
		while i < len(self.scores):
			if chr != self.chromosomes[i] - 1:
				chr = self.chromosomes[i] - 1
				j = 0 #snpsdata index
			pos = self.positions[i]
			#print i,chr, j,len(snpsds[chr].positions)
			while j < len(snpsds[chr].positions) and pos > snpsds[chr].positions[j]:
				j += 1
			if j < len(snpsds[chr].positions) and pos == snpsds[chr].positions[j]:
				self.snps.append(snpsds[chr].snps[j])
			i += 1

		if pos > snpsds[chr].positions[j]:
			while j < len(snpsds[chr].positions) and pos > snpsds[chr].positions[j]:
				j += 1
			if j < len(snpsds[chr].positions) and pos == snpsds[chr].positions[j]:
				self.snps.append(snpsds[chr].snps[j])

		if i != len(self.snps):
			print "Problems with loading SNPs", i, len(self.snps)

		print "Loaded", len(self.snps), " SNPs and", len(self.accessions), "accessions."



	def filterNicePeaks(self, scoreThreshold, singletonScoreThreshold, window=[20000, 20000], method=1):
		currScoreWeight = 0.2
		currChrom = -1
		newScores = []
		newSnps = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		singletonCount = 0
		lastSecondaryEnd = 0
		if method == 1:
			for i in range(0, len(self.positions)):
				if currChrom != self.chromosomes[i]: #Restart
					currChrom = self.chromosomes[i]
					startIndex = i #windowIndices 
					stopIndex = i #windowIndices 
					curScoreSum = currScoreWeight * self.scores[i]
					oldScore = 0
					numSNPs = 1
				currPos = self.positions[i]

				while currPos - self.positions[startIndex] > window[0]:
					curScoreSum -= self.scores[startIndex]
					startIndex += 1
					numSNPs -= 1

				while stopIndex + 1 < len(self.positions) and self.positions[stopIndex + 1] - currPos < window[1]:
					stopIndex += 1
					curScoreSum += self.scores[stopIndex]
					numSNPs += 1

				curScoreSum -= oldScore
				oldScore = currScoreWeight * numSNPs * self.scores[i]
				curScoreSum += oldScore

				if (numSNPs < 5 and self.scores[i] > singletonScoreThreshold) or  (numSNPs > 5 and curScoreSum / float(numSNPs) > scoreThreshold):
					if (numSNPs < 5 and self.scores[i] > singletonScoreThreshold):
						singletonCount += 1
					newScores.append(self.scores[i])
					newPositions.append(self.positions[i])
					newChromosomes.append(self.chromosomes[i])
					newMafs.append(self.mafs[i])
					newMarfs.append(self.marfs[i])
					newSnps.append(self.snps[i])

		elif method == 2:
			for i in range(0, len(self.scores)):
				if self.scores[i] >= singletonScoreThreshold:
					newScores.append(self.scores[i])
					newPositions.append(self.positions[i])
					newChromosomes.append(self.chromosomes[i])
					newMafs.append(self.mafs[i])
					newMarfs.append(self.marfs[i])
					newSnps.append(self.snps[i])

		# The following code locates the regions before and after the "nice" SNPs.
		j = 0
		for i in range(0, len(self.positions)):
			pos = self.positions[i]
			chr = self.chromosomes[i]
			if j < len(newPositions) and pos == newPositions[j] and chr == newChromosomes[j]:
				k = 0
				while  i + k > lastSecondaryEnd and pos - self.positions[i + k - 1] < window[0] and self.chromosomes[i + k - 1] == chr:
					k -= 1
				while i + k < len(self.positions) - 1 and self.positions[i + k] - pos < window[1] and self.chromosomes[i + k] == chr:
					if i + k > lastSecondaryEnd:
						self.secondaryScores.append(self.scores[i + k])
						self.secondaryPositions.append(self.positions[i + k])
						self.secondaryChromosomes.append(self.chromosomes[i + k])
					k += 1
				lastSecondaryEnd = i + k - 1
				j += 1


		self.snps = newSnps
		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs

		print singletonCount, "singletons were added"



	def filterMAF(self, minMaf=15):
		"""
		Filter out scores/pvalues which have maf<minMaf.		
		"""
		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		newSnps = []
		for i in range(0, len(self.scores)):
			if self.mafs[i] >= minMaf:
				newScores.append(self.scores[i])
				newPositions.append(self.positions[i])
				newChromosomes.append(self.chromosomes[i])
				newMafs.append(self.mafs[i])
				newMarfs.append(self.marfs[i])
				newSnps.append(self.snps[i])
		del self.snps
		del self.scores
		del self.positions
		del self.chromosomes
		del self.mafs
		del self.marfs

		self.snps = newSnps
		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs



	def filterMARF(self, minMarf=0.1):
		"""
		Filter out scores/pvalues which have maf<minMaf.		
		"""
		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		newSnps = []
		for i in range(0, len(self.scores)):
			if self.marfs[i] >= minMarf:
				newScores.append(self.scores[i])
				newPositions.append(self.positions[i])
				newChromosomes.append(self.chromosomes[i])
				newMafs.append(self.mafs[i])
				newMarfs.append(self.marfs[i])
				newSnps.append(self.snps[i])
		del self.snps
		del self.scores
		del self.positions
		del self.chromosomes
		del self.mafs
		del self.marfs

		self.snps = newSnps
		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs


	def filterScoreCutoff(self, scoreCutoff):
		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		newSnps = []
		for i in range(0, len(self.scores)):
			if self.scores[i] > scoreCutoff:
				newScores.append(self.scores[i])
				newPositions.append(self.positions[i])
				newChromosomes.append(self.chromosomes[i])
				newMafs.append(self.mafs[i])
				newMarfs.append(self.marfs[i])
				newSnps.append(self.snps[i])
		del self.snps
		del self.scores
		del self.positions
		del self.chromosomes
		del self.mafs
		del self.marfs

		self.snps = newSnps
		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs


	def getRegionSNPs(self, window=[25000, 25000], snpsDataFile=None):
		if self.snps:
			regionSNPs = []
			self.getRegions(window=window)
			for (snpIndex, startPos, endPos, chr, size, maxScore, maxPos, snpRank, regionRank) in self.regions:
				snp = SNP(self.positions[snpIndex], self.chromosomes[snpIndex], alleles=self.snps[snpIndex], accessions=self.accessions, score=maxScore)
				regionSNPs.append(snp)
			return regionSNPs


	def filterTopSNPs(self, n):
		self._rankScores_() #Making sure the ranks are updated
		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		newSnps = []
		l = self.orders[0:n]
		l.sort()

		for i in l:
			newScores.append(self.scores[i])
			newPositions.append(self.positions[i])
			newChromosomes.append(self.chromosomes[i])
			newMafs.append(self.mafs[i])
			newMarfs.append(self.marfs[i])
			newSnps.append(self.snps[i])
		del self.snps
		del self.scores
		del self.positions
		del self.chromosomes
		del self.mafs
		del self.marfs

		self.snps = newSnps
		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs


	def _sortByChrPos_(self):
		res_ls = zip(self.chromosomes, self.positions, self.scores, self.mafs, self.marfs, self.snps)
		res_ls.sort()
		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		newSnps = []
		for (chr, pos, score, maf, marf, snp) in res_ls:
			newScores.append(score)
			newPositions.append(pos)
			newChromosomes.append(chr)
			newMafs.append(maf)
			newMarfs.append(marf)
			newSnps.append(snp)

		del self.snps
		del self.scores
		del self.positions
		del self.chromosomes
		del self.mafs
		del self.marfs

		self.snps = newSnps
		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs


	def filterTopRegions(self, n, window=[25000, 25000], minScore=None):
		self._rankScores_()
		oldScores = self.scores
		oldPositions = self.positions
		oldChromosomes = self.chromosomes
		oldMafs = self.mafs
		oldMarfs = self.marfs
		oldSnps = self.snps
		oldGlobalRanks = self.globalRanks
		self.snps = []
		self.scores = []
		self.positions = []
		self.chromosomes = []
		self.mafs = []
		self.marfs = []
		self.globalRanks = []

		regionCount = 0
		i = 0

		res_ls = []
		if minScore:
			while regionCount < n:
				res_ls.append((oldChromosomes[self.orders[i]], oldPositions[self.orders[i]]))
				res_ls.sort()
				if oldScores[self.orders[i]] < minScore:
					break
				self.scores.append(oldScores[self.orders[i]])
				self.positions.append(oldPositions[self.orders[i]])
				self.chromosomes.append(oldChromosomes[self.orders[i]])
				self.mafs.append(oldMafs[self.orders[i]])
				self.marfs.append(oldMarfs[self.orders[i]])
				self.snps.append(oldSnps[self.orders[i]])
				self.globalRanks.append(oldGlobalRanks[self.orders[i]])
				regionCount = self._countRegions_(res_ls, window=window)
				i += 1
		else:
			while regionCount < n:
				res_ls.append((oldChromosomes[self.orders[i]], oldPositions[self.orders[i]]))
				res_ls.sort()
				self.scores.append(oldScores[self.orders[i]])
				self.positions.append(oldPositions[self.orders[i]])
				self.chromosomes.append(oldChromosomes[self.orders[i]])
				self.mafs.append(oldMafs[self.orders[i]])
				self.marfs.append(oldMarfs[self.orders[i]])
				self.snps.append(oldSnps[self.orders[i]])
				self.globalRanks.append(oldGlobalRanks[self.orders[i]])
				regionCount = self._countRegions_(res_ls, window=window)
				i += 1

		del oldScores
		del oldPositions
		del oldChromosomes
		del oldMafs
		del oldMarfs
		del oldSnps
		del oldGlobalRanks

		self._sortByChrPos_()

#	def get_snps(self):
#		"""
#		Returns this result's SNPs as a list of SNP objects.
#		"""
#		raise  NotImplementedError
#		for i, snp in enumerate(self.snps):
#			SNP(position,chromosome,accessions=None,alleles=None,snpsds=None,score=None,rank=None,regionRank=None)


	def mergeWith(self, snpResult):
		#pdb.set_trace()

		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		newSnps = []

		if len(snpResult.scores) == 0:
			return
		elif len(self.scores) == 0:
			newPositions = snpResult.positions
			newChromosomes = snpResult.chromosomes
			newScores = snpResult.scores
			newMafs = snpResult.mafs
			newMarfs = snpResult.marfs
			newSnps = snpResult.snps
		else:

			i1 = 0

			pos1 = self.positions[i1]
			chr1 = self.chromosomes[i1]

			for i2 in range(0, len(snpResult.positions)):
				pos2 = snpResult.positions[i2]
				chr2 = snpResult.chromosomes[i2]
				while i1 < len(self.positions) - 1 and (chr1, pos1) < (chr2, pos2):
					newPositions.append(self.positions[i1])
					newChromosomes.append(self.chromosomes[i1])
					newScores.append(self.scores[i1])
					newMafs.append(self.mafs[i1])
					newMarfs.append(self.marfs[i1])
					newSnps.append(self.snps[i1])
					i1 += 1
					pos1 = self.positions[i1]
					chr1 = self.chromosomes[i1]

				if i1 < len(self.positions) - 1 and (chr1, pos1) == (chr2, pos2):
					i1 += 1
					pos1 = self.positions[i1]
					chr1 = self.chromosomes[i1]

				newPositions.append(snpResult.positions[i2])
				newChromosomes.append(snpResult.chromosomes[i2])
				newScores.append(snpResult.scores[i2])
				newMafs.append(snpResult.mafs[i2])
				newMarfs.append(snpResult.marfs[i2])
				newSnps.append(snpResult.snps[i2])
				#pdb.set_trace()

		del self.snps
		del self.scores
		del self.positions
		del self.chromosomes
		del self.mafs
		del self.marfs
		del snpResult

		self.snps = newSnps
		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs




class RegionsTable(object):
	"""
	A table or regions X methods/phenotypes.
	"""
	def __init__(self, result_ls, window=[25000, 25000]):
		merged_result = result_ls[0].clone()
		for i in range(1, len(result_ls)):
			result = result_ls[i]
			merged_result.mergeWith(result)
		merged_result.getRegions(window=window)
		totalRegSize = 0
		if len(merged_result.positions):
			for reg in merged_result.regions:
				totalRegSize += reg[4]
		#print "The union of the results: Number of sign. SNPs: "+str(len(merged_result.scores))+", number of sign. regions: "+str(len(merged_result.regions))+", ave. region size: "+str(totalRegSize/float(len(merged_result.regions)))+".\n"

		self.regions = []
		for (snpIndex, startPos, endPos, chr, size, maxScore, maxPos, snpRank, regionRank) in merged_result.regions:
			self.regions.append((chr, startPos, endPos, size))
		del merged_result

		self.region_by_methods_table = [] # list[region_index][method/phenotype_index][snp_index]
		for (chr, startPos, endPos, size) in self.regions:  #For all regions
			methods_snps_ls = []
			for m_i in range(0, len(result_ls)):      #for all results
				result = result_ls[m_i]
				snps_ls = []

				if result.snps:
					##Finding all SNPs in region of interest
					regionRank = None
					snpRank = None
					result.getRegions(window=window)
					for (si, startPos2, endPos2, chr2, size2, mScore, mPos, sRank, rRank) in result.regions:
						if chr2 == chr and startPos <= startPos2 and endPos2 <= endPos:
							#print startPos,startPos2,endPos2,endPos
							regionRank = rRank
							snpRank = sRank
					#if not regionRank:
						#pdb.set_trace()

				##Finding all SNPs in region of interest
				for i in range(0, len(result.positions)):
					if result.chromosomes[i] == chr and startPos < result.positions[i] < endPos:
						if result.snps:
							#print result.snps[i]
							snp = SNP(result.positions[i], chr, alleles=result.snps[i], score=result.scores[i], rank=snpRank, regionRank=regionRank)
						else:
							snp = (chr, result.positions[i], result.scores[i], i)  #(chr,pos,score,index)
						snps_ls.append(snp)
				methods_snps_ls.append(snps_ls)
			self.region_by_methods_table.append(methods_snps_ls)

		self.resultNames = []
		for result in result_ls:
			self.resultNames.append(result.name)




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


def getCandidateGeneList(cgl_id, host="papaya.usc.edu", user="bvilhjal", passwd="bjazz32", db="stock_250k"):
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
	candGenes = []
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		gene = Gene(int(row[0]), int(row[1]), int(row[2]), name=row[4], description=row[5], dbRef=row[6])
		candGenes.append(gene)
	cursor.close ()
	conn.close ()
	print "Candiate gene-lists fetched"
	return candGenes


def get_gene_list(start_pos=None, end_pos=None, chr=None, host="gmi-ara-devel-be", include_intron_exons=True, \
		verbose=False, conn=None):
	import dbutils
	if not conn:
		new_conn = dbutils.connect_to_db(host, 'genome')
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
		sql_statement = "select distinct gm.chromosome, gm.start, gm.stop, g.locustag, \
		g.gene_symbol, g.description, g.dbxrefs from genome.entrezgene_mapping gm, genome.gene g where \
		gm.gene_id = g.gene_id and gm.chromosome=" + str(chr) + " and gm.stop>" + str(start_pos) + " and \
		gm.start<" + str(end_pos) + " order by gm.chromosome, gm.start, gm.stop"
		numRows = int(cursor.execute(sql_statement))
	else:
		sql_statement = "select distinct gm.chromosome, gm.start, gm.stop, g.locustag, g.gene_symbol, \
			g.description, g.dbxrefs from genome.entrezgene_mapping gm, genome.gene g \
			where gm.gene_id = g.gene_id order by gm.chromosome, gm.start, gm.stop"
		numRows = int(cursor.execute(sql_statement))
	if numRows == 0:
		pass
		#print sql_statment
		#print cursor.fetchone()
	        #import pdb;pdb.set_trace()

	genes = []
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		try:
			#chr, start, stop, gene_symbol, description, dbref,  
			gene = Gene(int(row[0]), int(row[1]), int(row[2]), name=row[4], description=row[5], dbRef=row[6], tairID=row[3])
			gene.tairID = row[6][5:]
			genes.append(gene)
		except Exception, err_str:
			#pass
			print err_str, ':'
			print row

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
					segments.append(Region(g.chromosome, row[0], row[1]))
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
	print(f.readline())
	reader = csv.reader(f)
	tair_ids = []
	gene_names = []
	if format == 1:
		for row in reader:
		      tair_ids.append(row[0].upper())
		      gene_names.append(row[1])
	f.close()
	return get_genes_w_tair_id(tair_ids), tair_ids




def get_genes_w_tair_id(tair_ids):
	conn = dbutils.connect_to_papaya("genome")
	cursor = conn.cursor()
	genes = []
	for tair_id in tair_ids:
		sql_statment = "select distinct gm.chromosome, gm.start, gm.stop, g.locustag, g.gene_symbol, g.description, g.dbxrefs from genome.entrezgene_mapping gm, genome.gene g where g.dbxrefs='TAIR:" + tair_id.upper() + "' and gm.gene_id = g.gene_id order by gm.chromosome, gm.start, gm.stop"
		numRows = int(cursor.execute(sql_statment))
		if numRows > 1:
			print "Found 2 copies:", sql_statment
		while(1):
			row = cursor.fetchone()
			if not row:
				break;
			try:
				#chr, start, stop, gene_symbol, description, dbref,  
				gene = Gene(int(row[0]), int(row[1]), int(row[2]), name=row[4], description=row[5], dbRef=row[6], tairID=row[3])
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



def qq_plots(results, num_dots, max_log_val, file_prefix, method_types=['kw', 'emma'], phen_name=None,
	      perm_pvalues=None, is_binary=False, **kwargs):
	"""
	Plots both log and normal qq plots.
	"""
	if file_prefix:
		log_pdf_file = file_prefix + "_qq_log.pdf"
		log_png_file = file_prefix + "_qq_log.png"
		pdf_file = file_prefix + "_qq.pdf"
		png_file = file_prefix + "_qq.png"
	else:
		log_pdf_file = None
		log_png_file = None
		pdf_file = None
		png_file = None
	qq_plot(results, num_dots, method_types=method_types, phenName=phen_name, pdfFile=pdf_file, pngFile=png_file,
	        perm_pvalues=perm_pvalues, isBinary=is_binary, kwargs=kwargs)
	log_qq_plot(results, num_dots, max_log_val, method_types=method_types, phenName=phen_name, pdfFile=log_pdf_file,
		    pngFile=log_png_file, perm_pvalues=perm_pvalues, isBinary=is_binary, kwargs=kwargs)


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

def qq_plot(results, numQuantiles, method_types=["kw", "emma"], phenName=None, pdfFile=None, pngFile=None,
	    perm_pvalues=None, **kwargs):
	"""
	QQ-plot for the given results.
	"""
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt

	plt.figure(figsize=(5, 4))
	#plt.figure(figsize=(10,8))
	#plt.figure(figsize=(4,3.5))
	plt.axes([0.15, 0.14, 0.82, 0.79])
	plt.plot([0, 1], [0, 1], "k", label="Expected")
	areas = []
	medians = []
	for method_type in method_types:
		result = results[method_type]
		label = method_type
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


def log_qq_plot(results, numDots, maxVal, method_types=['kw', 'emma'], phenName=None, pdfFile=None,
	        pngFile=None, perm_pvalues=None, **kwargs):
	"""
	log-transformed QQ-plot for the given results.
	"""
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	def _getExpectedLogQuantiles_():
		quantiles = []
		for i in range(1, numDots + 1):
			quantiles.append((float(i) / (numDots + 2.0)) * maxVal)
		return quantiles

	plt.figure(figsize=(5, 4))
	#plt.figure(figsize=(10,8))
	#plt.figure(figsize=(4,3.5))
	plt.axes([0.15, 0.14, 0.82, 0.79])
	maxVal = min(math.log10(len(results[method_types[0]].scores)), maxVal)
	minVal = (1.0 / numDots) * maxVal
	valRange = maxVal - minVal
	plt.plot([minVal, maxVal], [minVal, maxVal], "k", label="Expected")
	maxObsVals = []
	areas = []
	ds = []
	slopes = []
	for method_type in method_types:
		result = results[method_type]
		label = method_type
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
