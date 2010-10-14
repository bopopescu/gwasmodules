"""
Functions for ploting GWA results in pylab.
"""

import gwaResults
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math

def _get_chromosome_splits_(chromosomes, chr_set):
	"""
	Assumes chromosome list is sorted.
	"""
	import bisect
	chromosome_splits = [0]
	for c in chr_set:
		chromosome_splits.append(bisect.bisect(chromosomes, c))
	return chromosome_splits


def plotFilteredResult(result, pdfFile, minScore=0, maxScore=10, plotBonferroni=False, usePylab=True):

	newScores = []
	for score in result.scores:
		if score > maxScore:
			score = maxScore
		newScores.append(score)
	newSecondaryScores = []
	for score in result.secondaryScores:
		if score > maxScore:
			score = maxScore
		newSecondaryScores.append(score)
	minScore = min([min(newScores), min(newSecondaryScores), minScore])
	scoreRange = maxScore - minScore
	print result.chromosomeEnds
	offset = 0
	chromosomeSplits = result.getChromosomeSplit()
	secondaryChrSplits = result.getSecondaryChromosomeSplit()
	print chromosomeSplits
	print secondaryChrSplits
	ticksList1 = []
	ticksList2 = []
	textPos = []
	plt.figure(1, figsize=(20, 4))
	for i in range(0, len(result.chromosomeEnds)):
		index1 = chromosomeSplits[i][0]
		index2 = chromosomeSplits[i + 1][0]
		secondaryIndex1 = secondaryChrSplits[i][0]
		secondaryIndex2 = secondaryChrSplits[i + 1][0]
		scoreList = newScores[index1:index2]
		posList = result.positions[index1:index2]
		secScoreList = newSecondaryScores[secondaryIndex1:secondaryIndex2]
		secPosList = result.secondaryPositions[secondaryIndex1:secondaryIndex2]

		newPosList = []
		for pos in posList:
			newPosList.append(offset + pos)

		newSecPosList = []
		for pos in secPosList:
			newSecPosList.append(offset + pos)

		plt.plot(newSecPosList, secScoreList, ".", color=(0.75, 0.75, 0.75, 0.2))
		plt.plot(newPosList, scoreList, ".")
		oldOffset = offset
		textPos.append(offset + result.chromosomeEnds[i] / 2 - 2000000)
		offset += result.chromosomeEnds[i]
		if i < 4:
			plt.plot([offset, offset], [minScore - 0.05 * scoreRange, maxScore + 0.05 * scoreRange], "k--")
		for j in range(oldOffset, offset, 1000000):
			ticksList1.append(j)
		for j in range(0, result.chromosomeEnds[i], 1000000):
			if j % 5000000 == 0 and j < result.chromosomeEnds[i] - 1500000 :
				ticksList2.append(j / 1000000)
			else:
				ticksList2.append("")


	if plotBonferroni:
		plt.plot([0, sum(result.chromosomeEnds)], [6.68, 6.68], "k--")


	plt.axis([0, sum(result.chromosomeEnds), minScore - 0.05 * scoreRange, maxScore + 0.05 * scoreRange])
	plt.xticks(ticksList1, ticksList2)
	for i in range(0, len(textPos)):
		plt.text(textPos[i], minScore - scoreRange * 0.2, "Chr " + str(i + 1))
	plt.subplots_adjust(right=0.98)
	plt.subplots_adjust(left=0.03)
	plt.subplots_adjust(bottom=0.15)
	plt.subplots_adjust(top=0.9)
	plt.text(offset / 2, maxScore + scoreRange * 0.1, 'Position')
	plt.ylabel('-log(p-value)')
	plt.savefig(pdfFile, format="pdf")

	#   result.chromosomeEnds[i]

	"""
	plt.title(statistic.name)
	plt.xlabel('Observed theta')
	plt.ylabel('Estimated theta')
	plt.clf()
	"""

def plotResult(result, pdfFile=None, pngFile=None, minScore=None, maxScore=8.4,
	       percentile=98, type="pvals", ylab="$-$log$_{10}(p-$value$)$",
	       plotBonferroni=False, usePylab=False, cand_genes=None):

	"""
	type is either 'pvals' or 'scores'.
	
	"""

	if cand_genes:
		#processing candidate genes by chromosome
		chr_cand_genes = []
		for ch in [1, 2, 3, 4, 5]:
			cgs = []
			for cg in cand_genes:
				if cg.chromosome == ch:
					cgs.append(cg)
			chr_cand_genes.append(cgs)

	result = result.clone() #Very slow and memory inefficient

	result.filterPercentile(percentile / 100.0)

	newScores = result.scores
	if not maxScore:
		maxScore = max(newScores)
	if minScore:
	   pass
	   #minScore = min([min(newScores),minScore])
	else:
		if type == "pvals":
			minScore = 0
		else:
			minScore = min(newScores)
	scoreRange = maxScore - minScore
	print result.chromosomeEnds
	offset = 0
	chromosomeSplits = result.getChromosomeSplit()
	print chromosomeSplits
	ticksList1 = []
	ticksList2 = []
	textPos = []
	#plt.figure(figsize=(20,4))
	plt.figure(figsize=(12, 2.8))
	#plt.figure(figsize=(18,4))
	#plt.axes([0.03,0.105,0.967,0.79])
	plt.axes([0.045, 0.15, 0.95, 0.71])
	starPoints = [[], [], []]
	for i in range(0, len(result.chromosomeEnds)):
		index1 = chromosomeSplits[i][0]
		index2 = chromosomeSplits[i + 1][0]
		scoreList = newScores[index1:index2]
		posList = result.positions[index1:index2]

		if cand_genes:
			for cg in chr_cand_genes[i]:
				plt.axvspan(offset + cg.startPos, offset + cg.endPos, facecolor='k', alpha=0.5)


		newPosList = []
		for pos in posList:
			newPosList.append(offset + pos)

		newScoreList = []
		for (score, pos) in zip(scoreList, newPosList):
			if score > maxScore:
				starPoints[0].append(pos)
				starPoints[1].append(maxScore)
				starPoints[2].append(score)
				score = maxScore
			newScoreList.append(score)
		scoreList = newScoreList

		#plt.plot(newPosList,scoreList,".",markersize=5)
		plt.plot(newPosList, scoreList, ".", markersize=5, alpha=0.7)
		oldOffset = offset
		textPos.append(offset + result.chromosomeEnds[i] / 2 - 2000000)
		offset += result.chromosomeEnds[i]
		#if i<4:
		#	plt.plot([offset,offset],[minScore-0.05*scoreRange,maxScore+0.05*scoreRange],"k--")
		for j in range(oldOffset, offset, 5000000):
			ticksList1.append(j)
		for j in range(0, result.chromosomeEnds[i], 5000000):
			if j % 10000000 == 0 and j < result.chromosomeEnds[i] - 2500000 :
				ticksList2.append(j / 1000000)
			else:
				ticksList2.append("")


#				
#		for j in range(oldOffset,offset,3000000):
#			ticksList1.append(j)
#		for j in range(0,result.chromosomeEnds[i],3000000):
#			if j%6000000 == 0 and j < result.chromosomeEnds[i]-2500000 :
#				ticksList2.append(j/1000000)
#			else:
#				ticksList2.append("")

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




	if plotBonferroni:
		plt.plot([0, sum(result.chromosomeEnds)], [6.63, 6.63], "k-.")

	plt.axis([0, sum(result.chromosomeEnds), minScore - 0.05 * scoreRange, maxScore + 0.05 * scoreRange])
	plt.xticks(ticksList1, ticksList2)
	#for i in range(0,len(textPos)):
	#	plt.text(textPos[i],maxScore+scoreRange*0.08,"Chr "+str(i+1))
	#plt.subplots_adjust(right=0.99)
	#plt.subplots_adjust(left=0.05)
	#plt.subplots_adjust(bottom=0.15)
	#plt.subplots_adjust(top=0.9)
	#plt.text(offset/2,maxScore+scoreRange*0.1,'Position')
	if not ylab:
		if type == "pvals":
			plt.ylabel('$-log(p-$value$)$', size="large")

		else:
			plt.ylabel('score')
	else:
		plt.ylabel(ylab)
	plt.xlabel("Mb", size="large")

	if pdfFile:
		plt.savefig(pdfFile, format="pdf")
	if pngFile:
		plt.savefig(pngFile, format="png", dpi=300, bbox_inches='tight')
	if not (pdfFile or pngFile):
		plt.show()

	plt.clf()




def plotResultWithSecondRun(result, secondRunResult, pdfFile=None, pngFile=None, minScore=None,
						maxScore=None, percentile=90, srPercentile=90, type="pvals", ylab=None,
						plotBonferroni=False, bonferroniCutoffs=(6.64, 9.02), usePylab=False):

	result.filterPercentile(percentile / 100.0)
	secondRunResult.filterPercentile(srPercentile / 100.0)
	#print len(secondRunResult.scores)

	if maxScore:
		newScores = []
		newSecondRunScores = []
		for score in result.scores:
			if score > maxScore:
				score = maxScore
			newScores.append(score)
		for score in secondRunResult.scores:
			if score > maxScore:
				score = maxScore
			newSecondRunScores.append(score)
	else:
		newScores = result.scores
		newSecondRunScores = secondRunResult.scores
		maxScore = max(max(newScores), max(newSecondRunScores))
	if minScore:
	   minScore = min([min(newScores), minScore, min(newSecondRunScores)])
	else:
		if type == "pvals":
			minScore = 0
		else:
			minScore = min(newScores)
	scoreRange = maxScore - minScore
	print result.chromosomeEnds
	offset = 0
	chromosomeSplits = result.getChromosomeSplit()
	secondRunChromosomeSplits = secondRunResult.getChromosomeSplit()
	print chromosomeSplits
	print secondRunChromosomeSplits
	ticksList1 = []
	ticksList2 = []
	textPos = []
	plt.figure(figsize=(20, 4))
	for i in range(0, len(result.chromosomeEnds)):
		index1 = chromosomeSplits[i][0]
		index2 = chromosomeSplits[i + 1][0]
		srIndex1 = secondRunChromosomeSplits[i][0]
		srIndex2 = secondRunChromosomeSplits[i + 1][0]
		scoreList = newScores[index1:index2]
		srScoreList = newSecondRunScores[srIndex1:srIndex2]
		posList = result.positions[index1:index2]
		srPosList = secondRunResult.positions[srIndex1:srIndex2]

		newPosList = []
		for pos in posList:
			newPosList.append(offset + pos)

		newSRPosList = []
		for pos in srPosList:
			newSRPosList.append(offset + pos)

		if i % 2 == 0:
			basicColor = 'b'
			signColor = 'r'
			srBasicColor = 'y'
			srSignColor = 'm'
		else:
			basicColor = 'g'
			signColor = 'r'
			srBasicColor = 'c'
			srSignColor = 'm'

		for (pos, score) in zip(newSRPosList, srScoreList):
			if score >= bonferroniCutoffs[1]:
				plt.plot([pos], [score], "+", color=srSignColor)
			else:
				plt.plot([pos], [score], "+", color=srBasicColor)

		for (pos, score) in zip(newPosList, scoreList):
			if score >= bonferroniCutoffs[0]:
				plt.plot([pos], [score], ".", color=signColor)
			else:
				plt.plot([pos], [score], ".", color=basicColor)

		oldOffset = offset
		textPos.append(offset + result.chromosomeEnds[i] / 2 - 2000000)
		offset += result.chromosomeEnds[i]
		if i < 4:
			plt.plot([offset, offset], [minScore - 0.05 * scoreRange, maxScore + 0.05 * scoreRange], "k--")
		for j in range(oldOffset, offset, 1000000):
			ticksList1.append(j)
		for j in range(0, result.chromosomeEnds[i], 1000000):
			if j % 5000000 == 0 and j < result.chromosomeEnds[i] - 1500000 :
				ticksList2.append(j / 1000000)
			else:
				ticksList2.append("")


	if plotBonferroni:
		plt.plot([0, sum(result.chromosomeEnds)], [bonferroniCutoffs[0], bonferroniCutoffs[0]], "k-.")
		plt.plot([0, sum(result.chromosomeEnds)], [bonferroniCutoffs[1], bonferroniCutoffs[1]], "k-.")

	plt.axis([0, sum(result.chromosomeEnds), minScore - 0.05 * scoreRange, maxScore + 0.05 * scoreRange])
	plt.xticks(ticksList1, ticksList2)
	for i in range(0, len(textPos)):
		plt.text(textPos[i], minScore - scoreRange * 0.2, "Chr " + str(i + 1))
	plt.subplots_adjust(right=0.98)
	plt.subplots_adjust(left=0.05)
	plt.subplots_adjust(bottom=0.15)
	plt.subplots_adjust(top=0.9)
	plt.text(offset / 2, maxScore + scoreRange * 0.1, 'Position')
	if not ylab:
		if type == "pvals":
			plt.ylabel('-log(p-value)')
		else:
			plt.ylabel('score')
	else:
		plt.ylabel(ylab)

	if pdfFile:
		plt.savefig(pdfFile, format="pdf")
	if pngFile:
		plt.savefig(pngFile, format="png")
	if not (pdfFile or pngFile):
		plt.show()

	plt.clf()

	"""
	plt.title(statistic.name)
	plt.xlabel('Observed theta')
	plt.ylabel('Estimated theta')
	"""


def plot_raw_result(p_vals, chromosomes, positions, pdf_file=None, png_file=None, p_value_filter=0.02,
		min_score=None, max_score=None, plot_bonferroni=True, highlight_markers=None, highlight_loci=None,
		ylab="$-$log$_{10}(p-$value$)$"):
	""" 
	Plots a 'Manhattan' style GWAs plot.
	"""
	assert len(p_vals) == len(chromosomes) == len(positions), \
	'p_vals, chromosomes, positions, are not of equal length!'

	"Plotting a Manhattan-style plot with %i markers." % len(p_vals)
	num_scores = len(p_vals)

	chromosome_ends = []
	for i in range(len(positions) - 1):
		if chromosomes[i] != chromosomes[i + 1]:
			chromosome_ends.append(positions[i])
	chromosome_ends.append(positions[-1])
	scores = []
	new_chromosomes = []
	new_positions = []
	for i, p in enumerate(p_vals):
		if p <= p_value_filter:
			try:
				s = -math.log10(p)
			except Exception, err_str:
				print p, err_str
				s = 324
			scores.append(s)
			new_positions.append(positions[i])
			new_chromosomes.append(chromosomes[i])
	positions = new_positions
	chromosomes = new_chromosomes

	chromosome_splits = [0]
	for i in range(len(positions) - 1):
		if chromosomes[i] != chromosomes[i + 1]:
			chromosome_splits.append(i)
	chromosome_splits.append(len(chromosomes) - 1)

	if highlight_markers:
		new_h_markers = []
		for c, p, pval in highlight_markers:
			new_h_markers.append((c, p, -math.log10(pval)))
		highlight_markers = new_h_markers

	if not max_score:
		max_score = max(scores)
		if highlight_markers:
			h_scores = [s for c, p, s in highlight_markers]
			max_score = max(max_score, max(h_scores))
	if not min_score:
		min_score = 0


	score_range = max_score - min_score
	min_y = min_score - 0.05 * score_range
	max_y = max_score + 0.05 * score_range
	offset = 0
	chr_offsets = []
	ticksList1 = []
	ticksList2 = []
	textPos = []
	plt.figure(figsize=(12, 2.8))
	plt.axes([0.045, 0.15, 0.95, 0.71])
	starPoints = [[], [], []]
	for i in range(len(chromosome_ends)):
		chr_offsets.append(offset)
		index1 = chromosome_splits[i]
		index2 = chromosome_splits[i + 1]
		scoreList = scores[index1:index2]
		posList = positions[index1:index2]
		newPosList = [offset + pos for pos in posList]

		for s_i, (score, pos) in enumerate(zip(scoreList, newPosList)):
			if score > max_score:
				starPoints[0].append(pos)
				starPoints[1].append(max_score)
				starPoints[2].append(score)
				score = max_score
			scoreList[s_i] = score

		plt.plot(newPosList, scoreList, ".", markersize=5, alpha=0.7)
		oldOffset = offset
		textPos.append(offset + chromosome_ends[i] / 2 - 2000000)
		offset += chromosome_ends[i]
		for j in range(oldOffset, offset, 5000000):
			ticksList1.append(j)
		for j in range(0, chromosome_ends[i], 5000000):
			if j % 10000000 == 0 and j < chromosome_ends[i] - 2500000 :
				ticksList2.append(j / 1000000)
			else:
				ticksList2.append("")


	if highlight_loci:
		for c, p in highlight_loci:
			x_pos = chr_offsets[c - 1] + p
			plt.axvline(x=x_pos, ymin=min_y, ymax=max_y, color='#33ee66')


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
		plt.plot(xs, ys, ".", color="#ff9944", markersize=9, alpha=0.8)



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
			plt.text(starPoints[0][max_point] - 1000000, (starPoints[1][max_point] - 1) * 1.15, \
				 str(round(starPoints[2][max_point], 2)), rotation=45, size="small")




	if plot_bonferroni:
		bonferroni_threshold = -math.log10(1.0 / (num_scores * 20.0))
		plt.plot([0, sum(chromosome_ends)], [bonferroni_threshold, bonferroni_threshold], "k-.")

	plt.axis([0, sum(chromosome_ends), min_y, max_y])
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


def _plot_snp_pair_scatter_(plot_file_name, x_pos_list, y_pos_list, score_array, top_filter=0.05, fig_format='png'):
	print 'Plotting scatter plot.'
	plt.figure(figsize=(14, 12))
	x_min = min(x_pos_list)
	x_max = max(x_pos_list)
	x_range = x_max - x_min
	y_min = min(y_pos_list)
	y_max = max(y_pos_list)
	y_range = y_max - y_min
	plt.axis([x_min - 0.025 * x_range, x_max + 0.025 * x_range, y_min - 0.025 * y_range, y_max + 0.025 * y_range])

	sxy_list = []
	for i, x_pos, in enumerate(x_pos_list):
		for j, y_pos in enumerate(y_pos_list):
			sxy_list.append((score_array[i, j], x_pos, y_pos))

	sxy_list.sort()
	sxy_list = sxy_list[int((1.0 - top_filter) * len(sxy_list)):] #Filtering the top SNPs.
	l = map(list, zip(*sxy_list))
	cs = l[0]
	xs = l[1]
	ys = l[2]

	plt.plot([x_min, x_max], [y_min, y_max], 'k--')
	plt.scatter(xs, ys, c=cs, alpha=0.8, linewidths=0)
	plt.colorbar()
	plt.savefig(plot_file_name + '.' + fig_format, format=fig_format)
	plt.clf()



def plot_snp_pair_result(chromosomes, positions, score_array, file_prefix, all_in_one=False):
	"""
	Plots 5x5 2D scatter plots
	
	Or, all-chromosomes-in-one plot.
	"""
	chr_set = set(chromosomes)
	assert not all_in_one, "all_in_one option isn't implemented yet"

	chromosome_splits = _get_chromosome_splits_(chromosomes, chr_set)
	for i, chr_1 in enumerate(chr_set):
		for j, chr_2 in enumerate(chr_set):
			ci1, ci2 = chromosome_splits[i], chromosome_splits[i + 1]
			cj1, cj2 = chromosome_splits[j], chromosome_splits[j + 1]
			plot_file_name = file_prefix + '_chr%d_chr%d' % (chr_1, chr_2)
			x_pos_list = positions[ci1:ci2]
			y_pos_list = positions[cj1:cj2]
			score_array_slice = score_array[ci1:ci2, cj1:cj2]
			_plot_snp_pair_scatter_(plot_file_name, x_pos_list, y_pos_list, score_array_slice)


