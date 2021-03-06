#!/usr/bin/env python
"""
Examples:
	# draw map (-w), draw snp pca (-q), need_convert_alleles2binary (-x)
	PlotGroupOfSNPs.py -i /Network/Data/250k/tmp-yh/call_method_17.tsv -n /tmp/phenotype.tsv
		-o /tmp/phenotype_1_top_200_snps_a$j -s /Network/Data/A_lyrata/250K.ancestral -y 1 -a $j
		-w -x -q
	
	#2010-4-21 For BorevitzLab's seasoning data. phenotype is outputted by JBDataGWA.outputPhenotypePrediction() called within
	# JBDataGWA.checkEpistasisInJBLabData() of misc.py
	# include all the SNPs (-j 1_1-5_30000000) from the input file, no map, no binary-allele conversion.
	# no result_method from db due to a inexisting call method 100 (-l 100)
	dir=~/script/variation/data/JBLabSeasonFlowering/
	
	~/script/variation/src/PlotGroupOfSNPs.py
		-i $dir/DaysToFlower16replicates_full_model_9_logPhenotype_genotype_from_call_method_49_core482_with_FRI_del_chr_order_one_time_impute_yu_format.tsv
		-n $dir/DaysToFlower16replicates_full_model_9_logPhenotype_phenotype_predicted.tsv
		-o ~/script/variation/data/JBLabSeasonFlowering/full_model_9 -y 3 -a 4 -l 100 -j 1_1-5_30000000 -q -A -L
	
	# 2010-4-21 ditto. draw SNP id labels (-L) and sort the SNPs according PC1 (-S) but neither strain_PCA plot nor snp_pca plot
	# SNP label font size =5 (-F 5)
	~/script/variation/src/PlotGroupOfSNPs.py
		-i $dir/DaysToFlower16replicates_full_model_9_logPhenotype_genotype_from_call_method_49_core482_with_FRI_del_chr_order_one_time_impute_yu_format.tsv
		-n $dir/DaysToFlower16replicates_full_model_9_logPhenotype_phenotype_predicted.tsv
		-o ~/script/variation/data/JBLabSeasonFlowering/full_model_9 -y 3 -a 4 -l 100 -j 1_1-5_30000000 -L -S -F 5
	
	# all the swedish phenotypes
	~/script/variation/src/PlotGroupOfSNPs.py
		-i $dir/DaysToFlower16replicates_full_model_9_logPhenotype_genotype_from_call_method_49_core482_with_FRI_del_chr_order_one_time_impute_yu_format.tsv
		-n $dir/DaysToFlower16replicates_full_model_9_logPhenotype_phenotype_predicted.tsv -o ~/script/variation/data/JBLabSeasonFlowering/full_model_9_Swe
		-y 1,3,4,6,7,9,10,12,13,15,16,18,19,21,22,24 -a 4 -l 100 -j 1_1-5_30000000 -L -S -F 5

Description:
	Program to plot groups of SNPs for one phenotype/one analysis method. both Strains and SNPs go under PCA clustering.
"""

import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib as mpl; mpl.use("Agg")	#to avoid popup and collapse in X11-disabled environment
import time, csv, cPickle
import warnings, traceback
import pylab
import ImageColor
import numpy
from sets import Set
from matplotlib.patches import Polygon, CirclePolygon, Ellipse, Wedge
from pymodule import PassingData, figureOutDelimiter, getColName2IndexFromHeader, getListOutOfStr, SNPData, read_data,\
	assignMatPlotlibHueColorToLs, drawName2FCLegend
from pymodule import number2complement
from pymodule.algorithm	import pca_module
from variation.src import Stock_250kDB, StockDB
from variation.src.enrichment.GeneListRankTest import GeneListRankTest	#GeneListRankTest.getGeneList()
from variation.src.association.Association import Association
from variation.src.popgen.PhenotypeOfAncestralDerivedAllele import PhenotypeOfAncestralDerivedAllele
from variation.src.common import get_chr_id2size, get_chr_id2cumu_size, getEcotypeInfo
from variation.src.db.output.OutputPhenotype import OutputPhenotype
from variation.src import AbstractVariationMapper

class PlotGroupOfSNPs(GeneListRankTest):
	__doc__ = __doc__
	option_default_dict = AbstractVariationMapper.option_default_dict.copy()
	option_default_dict.update({
							('phenotype_fname', 1, ): [None, 'n', 1, 'phenotype file', ],\
							('ancestral_allele_fname', 0, ): [None, 's', 1, 'File containing the ancestral alleles to polarize SNPs', ],\
							('minus_log_pvalue', 0, ): [0, '', 0, 'toggle -log(pvalue)', ],\
							('phenotype_method_id_ls', 1, ): [None, 'y', 1, 'a comma-dash-separated list of phenotype_method ids in the phenotype file. i.e. 1,4-7',],\
							('call_method_id', 0, int):[17, 'l', 1, 'Restrict results based on this call_method. Default is no such restriction.'],\
							('analysis_method_id', 0, int):[7, 'a', 1, 'Restrict results based on this analysis_method. Default is no such restriction.'],\
							('no_of_top_hits', 1, int): [200, 'f', 1, 'how many number of top hits based on score or -log(pvalue).'],\
							('snp_region_tup', 0, ):[None, 'j', 1, 'designating a region on a chromosome where SNPs should be drawn. like 4_200000-5_700000. if this is given, no_of_top_hits will not be used.'],\
							('country_order_type', 1, int): [1, 'g', 1, 'How to order countries from where strains are from. 1: order by latitude. 2: by longitude.'],\
							('draw_map', 0, int):[0, 'w', 0, 'toggle to draw a geographic map of strains on the right side'],\
							('need_convert_alleles2binary', 0, int):[0, 'x', 0, 'toggle to binary-convert SNPs alleles'],\
							('draw_snp_pca', 0, int):[0, 'q', 0, 'toggle to draw a 2D PCA clustering of SNPs on the top. this would toggle sortSNPByPCA_value.'],\
							('draw_strain_pca', 0, int):[0, 'A', 0, 'toggle to draw a 2D PCA clustering of strains on the left'],\
							('draw_snp_id_label', 0, int):[0, 'L', 0, 'toggle to label each SNP column with SNP ID when drawing the SNP matrix '],\
							('sortSNPByPCA_value', 0, int):[0, 'S', 0, 'toggle to sort the SNPs by PC1 from PCA when drawing the SNP matrix '],\
							('snp_label_font_size', 1, int): [1, 'F', 1, 'the size of the font used to label SNPs when drawing SNP matrix'],\
							('drawGrid', 0, int):[0, '', 0, 'toggle to draw a grid for the SNP-matrix/haplotype'],\
							})
	country_order_type2name = {1: 'Country ordered by latitude',
							2: 'Country ordered by longitude'}
	def __init__(self,  **keywords):
		"""
		2008-10-14
			add option country_order_type, draw_map
		2008-10-07
		"""
		AbstractVariationMapper.__init__(self, inputFnameLs=None, **keywords)
		if self.phenotype_method_id_ls:
			self.phenotype_method_id_ls = getListOutOfStr(self.phenotype_method_id_ls, data_type=int)
	
	snp_value2five_color = {-1:'k', 0:'w', 1:'b', 2:'g', 3:'r'}	#deletion, NA, allele1, allele2, heterozygous
	
	@classmethod
	def drawSNPMtrix(cls, axe_snp_matrix, subSNPData, top_snp_data, StrainID2PCAPosInfo, SNPID2PCAPosInfo, \
					ecotype_info, strain_id_label_x_offset=0.9, snp_id_label_y_offset=0.9, strain_snp_label_font_size=1, \
					draw_snp_id_label=True, strain_id_label_x_offset_extra=None, snpData_before_impute=None, \
					useAlleleToDetermineAlpha=False,\
					snp_value2color=None, snp_label_font_size=None, drawGrid=False):
		"""
		2012.3.27 default index_type for each SNP (if not present in snp_id2index_type) is 2 (based on MAF).
		2010-4-23
			add argument drawGrid: whether to draw a grid (horizontal, vertical lines across the whole matrix)
		2010-4-21
			add argument snp_label_font_size. default is None, which means it'll be equal to strain_snp_label_font_size.
		2010-2-2
			becomes classmethod
		2009-4-30
			add argument snp_value2color
		2009-3-23
			add argument useAlleleToDetermineAlpha, which decides whether to use allele in subSNPData or association score to decide alpha
			add a grid to the whole matrix
			change the color coding for major/minor from white/red to blue/yellow
		2008-12-09
			the SNP matrix in subSNPData, if it's binary, it now uses 0,1 to represent alleles, rather than 1,2.
		2008-12-01
			add option snpData_before_impute
		2008-11-25
			allow the cell_x_len to vary by allowing values of SNPID2PCAPosInfo.snp_id2img_x_pos to become 2-element tuple or list
				which stores the start and stop x-position of the cell
			add option strain_id_label_x_offset_extra
				if strain_id_label_x_offset_extra is not None, the odd/even-row (strain) labels stagger/zigzag with an extra x offset.
		2008-10-07
			1. allele=1 no color (white), 2(red)
			2. snpData.snp_id2index_type.
				index_type=1(polarized), shape=rectangle.
				index_type=2(MAF), circle.
				index_type=3(weirdo), triangle
			3. alpha of each cell determined by the score in top_snp_data.score_ls
			4. label strain, snp plainlys
		"""
		sys.stderr.write("Drawing SNP matrix ...")
		no_of_snps = len(top_snp_data.snp_id_ls)
		max_score = max(top_snp_data.score_ls)
		min_score = min(top_snp_data.score_ls)
		if min_score==max_score:
			min_score=max_score/2.0
			max_score = numpy.max(subSNPData.data_matrix)
			min_score = numpy.min(subSNPData.data_matrix)
		cell_x_len = getattr(SNPID2PCAPosInfo, 'step', None)
		cell_y_len = StrainID2PCAPosInfo.step
		if cell_x_len is not None and cell_y_len is not None:
			radius = min(cell_x_len, cell_y_len)/2.	#half of the smaller one
		if not snp_label_font_size:	#2010-4-21
			snp_label_font_size = strain_snp_label_font_size
		#2009-3-23 2 variables below is used to draw horizontal lines as part of the grid
		min_img_x_pos = None
		max_img_x_pos = None
		for i in range(no_of_snps):
			snp_id = top_snp_data.snp_id_ls[i]
			#if snp_id not in SNPID2PCAPosInfo.snp_id2img_x_pos:	#skip it if it's not in the dictionary
			# 2012.3.26 if snp_id not in snp_id2img_x_pos, something is wrong (gwas result file has a locus id not in snp_id2img_x_pos).
			#	continue
			#draw snp label
			snp_img_x_pos = SNPID2PCAPosInfo.snp_id2img_x_pos[snp_id]
			if type(snp_img_x_pos)==tuple or type(snp_img_x_pos)==list:	#2008-11-25 it's tuple used in DrawSNPRegion.py
				cell_x_len = abs(snp_img_x_pos[1]-snp_img_x_pos[0])
				snp_img_x_pos = min(snp_img_x_pos)
				#reset radius
				radius = min(cell_x_len, cell_y_len)/2.	#half of the smaller one
			
			if drawGrid:
				#2009-3-23 vertical line for the grid
				axe_snp_matrix.vlines([snp_img_x_pos], ymin=0, ymax=1, linewidth=0.15)
			#2009-3-23 for horizontal lines
			if min_img_x_pos is None or snp_img_x_pos<min_img_x_pos:
				min_img_x_pos=snp_img_x_pos
			if max_img_x_pos is None or snp_img_x_pos>max_img_x_pos:
				max_img_x_pos=snp_img_x_pos
			
			if draw_snp_id_label:
				axe_snp_matrix.text(snp_img_x_pos, snp_id_label_y_offset, snp_id, rotation='vertical', \
								horizontalalignment ='left', verticalalignment='bottom', size=snp_label_font_size)
			
			col_index = subSNPData.col_id2col_index[snp_id]
			index_type = subSNPData.snp_id2index_type[snp_id]
			#index_type = subSNPData.snp_id2index_type.get(snp_id, 2)	#2012.3.27 default is 2 (based on MAF)
			score = top_snp_data.score_ls[i]
			alpha = (score-min_score)/(max_score-min_score)*(1-0.2)+0.2	#alpha from 0.2 to 1. can't be low.
			for j in range(len(StrainID2PCAPosInfo.strain_id_ls)):
				strain_id = StrainID2PCAPosInfo.strain_id_ls[j]
				strain_img_y_pos = StrainID2PCAPosInfo.strain_id2img_y_pos[strain_id]
				row_index = subSNPData.row_id2row_index[strain_id]
				allele = subSNPData.data_matrix[row_index][col_index]
				_linewidth=0
				if snpData_before_impute:	#2008-12-01	snpData_before_impute is available. highlight missing data.
					before_impute_row_index = snpData_before_impute.row_id2row_index.get(strain_id)
					before_impute_col_index = snpData_before_impute.col_id2col_index.get(snp_id)

					if before_impute_row_index is not None and before_impute_col_index is not None:
						before_impute_allele = snpData_before_impute.data_matrix[before_impute_row_index][before_impute_col_index]
						if before_impute_allele==-2 or before_impute_allele==0 or numpy.isnan(before_impute_allele):
							#allele = -1
							_linewidth = 0.3	#mark the imputed with non-zero linewidth
					else:	#NA
						#allele = -1
						_linewidth = 0.3
					"""
					#2008-12-05 temporary for CNV output
					if before_impute_row_index is not None and before_impute_col_index is not None:
						before_impute_allele = snpData_before_impute.data_matrix[before_impute_row_index][before_impute_col_index]
						if before_impute_allele==-2 or numpy.isnan(before_impute_allele):
							allele=1
						elif before_impute_allele==0:	#white color
							allele=1
						elif before_impute_allele==-1:	#blue
							allele=3
						elif before_impute_allele==1:	#red
							allele=2
						else:
							allele=1	#white default
					else:	#NA
						allele = 1
					"""
				if allele==0:
					facecolor = 'b'
				elif allele==1:
					facecolor = 'y'
				else:
					facecolor = 'g'
				if snp_value2color is not None:	#2009-3-27
					facecolor = snp_value2color.get(allele)
					if facecolor is None:
						facecolor = 'g'
				
				if useAlleleToDetermineAlpha:	#2009-3-23
					alpha = (allele-min_score)/(max_score-min_score)	#2008-12-05 use transparency to show CNV data intensity
				
				if i ==0:	#draw strain label, first SNP
					ecotype_id = int(strain_id)
					ecotype_obj = ecotype_info.ecotype_id2ecotype_obj.get(ecotype_id)
					if ecotype_obj:
						nativename = ecotype_obj.nativename
						country = ecotype_obj.country
					else:
						nativename = ''
						country = 'UNK'
					#nativename = ecotype_info.ecotype_id2ecotype_obj[ecotype_id].nativename
					strain_label = '%s_%s_%s'%(strain_id, nativename, country)	#add country into strain label

					if j%2==1 and strain_id_label_x_offset_extra is not None:	#even rows, push labels a bit to right
						x_pos = strain_id_label_x_offset + strain_id_label_x_offset_extra
					else:
						x_pos = strain_id_label_x_offset
					axe_snp_matrix.text(x_pos, strain_img_y_pos, strain_label, \
							horizontalalignment ='left', verticalalignment='bottom', size=strain_snp_label_font_size)
				center = (snp_img_x_pos+cell_x_len/2., strain_img_y_pos+cell_y_len/2.)
				if index_type==1:	#(polarized, ancestral vs. derived), shape=rectangle.
					xs = [snp_img_x_pos, snp_img_x_pos+cell_x_len, snp_img_x_pos+cell_x_len, snp_img_x_pos]
					ys = [strain_img_y_pos, strain_img_y_pos, strain_img_y_pos+cell_y_len, strain_img_y_pos+cell_y_len]
					patch = Polygon(zip(xs,ys), facecolor=facecolor, linewidth=_linewidth, alpha=alpha)	#
				elif index_type==2:	#based on MAF
					patch = Ellipse(center, cell_x_len, cell_y_len, facecolor=facecolor, linewidth=_linewidth, alpha=alpha)
				else:	#others
					xs = [snp_img_x_pos, snp_img_x_pos+cell_x_len, snp_img_x_pos+cell_x_len/2.]
					ys = [strain_img_y_pos, strain_img_y_pos, strain_img_y_pos+cell_y_len]
					patch = Polygon(zip(xs,ys), facecolor=facecolor, linewidth=_linewidth, alpha=alpha)
				axe_snp_matrix.add_patch(patch)
		
		if drawGrid:
			#2009-3-23 draw a grid, horizontal lines first. vertical lines are tricky
			# cuz value of SNPID2PCAPosInfo.snp_id2img_x_pos might be a tuple, rather than a simple integer
			axe_snp_matrix.hlines(StrainID2PCAPosInfo.strain_id2img_y_pos.values(), xmin=min_img_x_pos, \
							xmax=max_img_x_pos, linewidth=0.15)
		sys.stderr.write("Done.\n")
	
	@classmethod
	def drawStrainPCA(cls, axe_strain_pca, axe_strain_map, axe_strain_map_pca_cover, axe_strain_pca_legend, \
					StrainID2PCAPosInfo, ecotype_info, phenData, \
					phenotype_col_index, phenotype_cmap, phenotype_norm, rightmost_x_value=1.05,\
					country_order_name='', alpha=0.8, strain_color_type=2, pca2map_line_color='g', ecotype_width_on_map=4,\
					draw_lines_to_axe_snp_matrix=True, strain_size_on_axe_strain_pca=14, draw_axe_strain_map=True, \
					pic_area=[-15,30,38,66], map_pca_line_alpha=0.2, map_pca_linewidth=0.2):
		"""
		2010-4-26
			reduce the ecotype_width_on_map default value from 9 to 4.
		2010-2-2
			becomes classmethod
		2008-11-14
			add option pca2map_line_color, color specification for the lines that connect objects from the axe_strain_pca
				to axe_strain_map. default is 'g'
			add options: ecotype_width_on_map=9, draw_lines_to_axe_snp_matrix=True, strain_size_on_axe_strain_pca=14, draw_axe_strain_map=True
		2008-11-08
			add axe_strain_map right above axe_strain_pca to show where the strains are located geographically.
			
			axe_strain_map_pca_cover has lines drawn to link the same strain between axe_strain_pca and axe_strain_map
		2008-10-14
			add phenData, phenotype_col_index, phenotype_cmap, phenotype_norm
			add strain_color_type: 1=country(order by either latitude or longitude); 2=phenotype
			add code to color strain according to phenotype value
		2008-10-08
			order strain country according to ecotype_info.country2order
		2008-10-07
			1. find all distinct countries
			2. assign colors to the countries and draw a country legend
			3. plot each strain on the axe by pca x/y value, color by its country
			4. draw line linking each strain from axe_snp_matrix to the point in axe_strain_pca
		"""
		sys.stderr.write("Drawing strain PCA ...")
		if draw_axe_strain_map:
			map_data = cls.justDrawMap(axe_strain_map, pic_area=pic_area)
		
		if strain_color_type==1:
			country_abbr_set = Set()
			for strain_id in StrainID2PCAPosInfo.strain_id_ls:
				ecotype_id = int(strain_id)
				ecotype_obj = ecotype_info.ecotype_id2ecotype_obj.get(ecotype_id)
				if ecotype_obj:
					country = ecotype_obj.country
				else:
					country = 'UNK'	#unknown , same abbr as in database for unknown country origin
				if country in ecotype_info.country2order:
					order_of_country = ecotype_info.country2order[country]
				else:
					order_of_country = -1	#unknown order, put before any other country
				country_abbr_set.add((order_of_country, country))
			order_country_ls = list(country_abbr_set)
			order_country_ls.sort()
			country_ls = [row[1] for row in order_country_ls]
			country2fc = assignMatPlotlibHueColorToLs(country_ls)
			drawName2FCLegend(axe_strain_pca_legend, country_ls, country2fc, no_edge_color=True, title=country_order_name, font_size=2, alpha=alpha)
		
		#freeze transformations of all 3 axes
		for ax in [axe_strain_pca, axe_strain_map, axe_strain_map_pca_cover]:
			if ax:
				#ax.transData.freeze()  # eval the lazy objects
				#ax.transAxes.freeze()
				pass	# 2010-2-27
		
		for strain_id in StrainID2PCAPosInfo.strain_id_ls:
			x_value = StrainID2PCAPosInfo.strain_id2pca_x[strain_id]
			y_value = StrainID2PCAPosInfo.strain_id2pca_y[strain_id]
			
			img_y_pos = StrainID2PCAPosInfo.strain_id2img_y_pos[strain_id]
			if draw_lines_to_axe_snp_matrix:
				#connect from axe_snp_matrix to axe_strain_pca
				axe_strain_pca.plot([rightmost_x_value, x_value], [img_y_pos, y_value], c='b', linestyle='--', alpha=0.2, \
								linewidth=0.2)
			
			ecotype_id = int(strain_id)
			ecotype_obj = ecotype_info.ecotype_id2ecotype_obj.get(ecotype_id)
			if ecotype_obj:
				lat, lon = ecotype_obj.latitude, ecotype_obj.longitude
			else:
				sys.stderr.write("Warning: Ecotype %s not in ecotype_info (fetched from stock db).\n"%ecotype_id)
				continue
			#strain color according to phenotype
			phenotype_row_index = phenData.row_id2row_index[strain_id]
			phenotype = phenData.data_matrix[phenotype_row_index][phenotype_col_index]
			strain_fc = phenotype_cmap(phenotype_norm(phenotype))
			if numpy.isnan(phenotype):
				linewidth=0.2
				strain_fc = 'w'
				edgecolor = 'k'
				_alpha = 0	#facecolor gets very transparent
			else:
				linewidth=0
				strain_fc = strain_fc
				edgecolor = 'k'
				_alpha = alpha
			if draw_axe_strain_map:
				if lat and lon:
					x, y = map_data.m(lon, lat)
					
					#mark the strain on the map
					axe_strain_map.scatter([x],[y], s=ecotype_width_on_map, linewidth=linewidth, facecolor=strain_fc, alpha=_alpha, zorder=10)
					canvas_x, canvas_y = axe_strain_map.transData.xy_tup((x,y))
					strain_map_cover_x,  strain_map_cover_y = axe_strain_map_pca_cover.transData.inverse_xy_tup((canvas_x,canvas_y))
					
					canvas_x, canvas_y = axe_strain_pca.transData.xy_tup((x_value, y_value))
					strain_pca_cover_x,  strain_pca_cover_y = axe_strain_map_pca_cover.transData.inverse_xy_tup((canvas_x,canvas_y))
					#connect the strain on axe_strain_pca to axe_strain_map
					if pca2map_line_color:
						axe_strain_map_pca_cover.plot([strain_pca_cover_x, strain_map_cover_x], [strain_pca_cover_y, strain_map_cover_y], c=pca2map_line_color, \
												linestyle='--', alpha=map_pca_line_alpha, linewidth=map_pca_linewidth)
					else:
						axe_strain_map_pca_cover.plot([strain_pca_cover_x, strain_map_cover_x], [strain_pca_cover_y, strain_map_cover_y],\
												linestyle='--', alpha=map_pca_line_alpha, linewidth=map_pca_linewidth)
			
			if strain_color_type==1:	#change strain color
				if ecotype_obj:
					country = ecotype_obj.country
				else:
					country = 'UNK'
				strain_fc = country2fc[country]
			
			#plot strain dot at last to avoid being shadowed by lines above
			axe_strain_pca.scatter([x_value],[y_value], s=strain_size_on_axe_strain_pca, linewidth=linewidth, facecolor=strain_fc, alpha=_alpha, zorder=10)
		
		#free transformations of all 3 axes
		for ax in [axe_strain_pca, axe_strain_map, axe_strain_map_pca_cover]:
			if ax:
				#ax.transData.thaw()  # eval the lazy objects
				#ax.transAxes.thaw()
				pass	# 2010-2-27
		if draw_axe_strain_map:
			axe_strain_map.set_xlim(map_data.map_xlim)
			axe_strain_map.set_ylim(map_data.map_ylim)
		
		#axe_strain_pca.set_title('Strain by PC2(x-axis) vs PC1(y-axis)')
		axe_strain_pca.set_xlabel('Strain PC2 %.2f%s'%(StrainID2PCAPosInfo.x_var*100, '%'))
		axe_strain_pca.set_ylabel('Strain PC1 %.2f%s'%(StrainID2PCAPosInfo.y_var*100, '%'))
		
		sys.stderr.write("Done.\n")
	
	def drawSNPPCA(self, axe_snp_pca, axe_snp_pca_legend, SNPID2PCAPosInfo, chr_pos_info, lowest_y_value=-0.05, \
				highest_y_value=1.02, lowest_x_value=0, alpha=0.8):
		"""
		2008-10-08
			connect SNP dots in axe_snp_pca to 5 chromosomes in the same axe
		2008-10-07
		"""
		sys.stderr.write("Drawing SNP PCA ...")
		chr_set = Set()
		for snp_id in SNPID2PCAPosInfo.snp_id_ls:
			chr, pos = snp_id.split('_')
			chr_set.add(chr)
		chr_ls = list(chr_set)
		chr_ls.sort()
		chr2fc = assignMatPlotlibHueColorToLs(chr_ls)
		drawName2FCLegend(axe_snp_pca_legend, chr_ls, chr2fc, no_edge_color=True, title='SNP chromosome legend', font_size=4, alpha=alpha)
		
		#draw chromosomes
		start_pos = lowest_x_value
		self.drawChromosome(axe_snp_pca, chr_pos_info.chr_id_ls, start_pos, chr_pos_info.scale_chr_id2cumu_size, highest_y_value, \
						chr_pos_info.scale_chr_gap, verticalalignment='bottom', font_size=4)
		
		for snp_id in SNPID2PCAPosInfo.snp_id_ls:
			chr, pos = snp_id.split('_')
			chr_fc = chr2fc[chr]
			x_value = SNPID2PCAPosInfo.snp_id2pca_x[snp_id]
			y_value = SNPID2PCAPosInfo.snp_id2pca_y[snp_id]
			
			img_x_pos = SNPID2PCAPosInfo.snp_id2img_x_pos[snp_id]
			axe_snp_pca.plot([img_x_pos, x_value], [lowest_y_value, y_value], c='b', linestyle='--', alpha=0.2, linewidth=0.2)
			
			#connect SNP dots in axe_snp_pca to 5 chromosomes in the same axe
			chr = int(chr)	#watch this conversion
			pos = int(pos)
			chr_size = float(chr_pos_info.chr_id2size[chr])
			cumu_pos = chr_pos_info.scale_chr_id2cumu_size[chr]-chr_pos_info.scale_chr_id2size[chr]+pos/chr_pos_info.total_cumu_size
			axe_snp_pca.plot([cumu_pos, x_value], [highest_y_value, y_value], c='g', linestyle='--', alpha=0.2, linewidth=0.2)
			#plot the SNP dot at last to avoid being shadowed by those lines
			axe_snp_pca.scatter([x_value],[y_value], s=10, linewidth=0, facecolor=chr_fc, alpha=alpha)
		axe_snp_pca.set_xlabel('PC1 %.2f%s'%(SNPID2PCAPosInfo.x_var*100, '%'), size=3)
		axe_snp_pca.set_ylabel('PC2 %.2f%s'%(SNPID2PCAPosInfo.y_var*100, '%'))
		sys.stderr.write("Done.\n")
	
	@classmethod
	def drawPhenotype(cls, axe, StrainID2PCAPosInfo, phenData, phenotype_col_index, phenotype_method_id=None, ecotype_info=None,
					phenotype_cmap=None, phenotype_norm=None, phenotypeTitle=None):
		"""
		2010-4-21
			add phenotypeTitle, if not None, add a title
			argument phenotype_method_id is not used.
		2010-2-5
			add arguments phenotype_cmap & phenotype_norm offer an option to draw a phenotype bars in color
		2010-2-2
			becomes classmethod
		2009-3-23
			draw blue line for NA phenotype
		2008-11-30
			fix a bug.
				axe.hlines draws a gradient-diminishing line when phenotype is numpy.nan.
				now use axe.plot to draw. NA-phenotype simply won't show up with anything.
		2008-10-09
			split out of drawMap()
			1. draw each phenotype with a bar (temporary solution)
		"""
		sys.stderr.write("Drawing phenotype ...")
		for strain_id in StrainID2PCAPosInfo.strain_id_ls:
			img_y_pos = StrainID2PCAPosInfo.strain_id2img_y_pos[strain_id]
			phenotype_row_index = phenData.row_id2row_index[strain_id]
			
			phenotype = phenData.data_matrix[phenotype_row_index][phenotype_col_index]
			if phenotype_cmap and phenotype_norm:	# 2010-2-5 if phenotype_cmap and phenotype_norm are given, use a colored bar with same length
				color = phenotype_cmap(phenotype_norm(phenotype))
				phenotype = phenotype_norm.vmax	# 
			elif numpy.isnan(phenotype):	#draw blue line for NA phenotype
				phenotype = numpy.nanmax(phenData.data_matrix[:, phenotype_col_index])
				color = 'b'
			else:
				color = 'k'
			axe.plot([0, phenotype], [img_y_pos, img_y_pos], linewidth=0.2, c=color)
		if phenotypeTitle:
			axe.set_title(phenotypeTitle, fontsize=3)
		sys.stderr.write("Done.\n")
	
	@classmethod
	def justDrawMap(cls, axe, pic_area=[-15,30,38,66]):
		"""
		2010-2-2
			becomes classmethod
		"""
		from matplotlib.toolkits.basemap import Basemap
		
		#pic_area = [left longitude, bottom latitude, right longitude, upper latitude]
		#pic_area=[-15,30,38,66]	#[-125,-38,180,70] covers everything, [-125,10,90,70] covers US and europe.
		m = Basemap(llcrnrlon=pic_area[0],llcrnrlat=pic_area[1],urcrnrlon=pic_area[2],urcrnrlat=pic_area[3],\
				resolution='c',projection='mill', ax=axe)	#resolution: c=crude, l=low, i=intermediate, h=high
		
		#m.drawcoastlines()
		m.drawparallels(pylab.arange(-90,90,30), labels=[1,1,0,1], size=2, linewidth=0.05)
		m.drawmeridians(pylab.arange(-180,180,30), labels=[1,1,0,1], size=2, linewidth=0.05)
		m.fillcontinents()
		m.drawcountries(linewidth=0.05)
		#m.drawstates()
		map_xlim = axe.get_xlim()
		map_ylim = axe.get_ylim()
		return_data = PassingData(m=m, map_xlim=map_xlim, map_ylim=map_ylim)
		return return_data
		
	@classmethod
	def drawMap(cls, axe_map_frame, axe_map, StrainID2PCAPosInfo, phenData, phenotype_col_index, phenotype_method_id, \
			ecotype_info, phenotype_cmap, phenotype_norm):
		"""
		2010-2-2
			becomes classmethod
		2009-3-23
			draw empty black circles for NA phenotypes
		2008-10-07
			1. draw map, locate ecotype and color according to phentoype
			2. because the axe_map would be contracted to the center(default) of the original axe_map due to Basemap,
				x-lim, y-lim are also totally different (physical distance on the earth in meters),
				have to connect strain y-pos to each strain on the map in another axe, axe_map_frame through coordinates transformation
		"""
		sys.stderr.write("Drawing map ...")
		map_data = cls.justDrawMap(axe_map)		
		
		#fix the two transformations before doing cross-axe drawings
		axe_map.transData.freeze()  # eval the lazy objects
		axe_map.transAxes.freeze()
		axe_map_frame.transData.freeze()  # eval the lazy objects
		axe_map_frame.transAxes.freeze()
		for strain_id in StrainID2PCAPosInfo.strain_id_ls:
			img_y_pos = StrainID2PCAPosInfo.strain_id2img_y_pos[strain_id]
			phenotype_row_index = phenData.row_id2row_index[strain_id]
			phenotype = phenData.data_matrix[phenotype_row_index][phenotype_col_index]
			ecotype_id = int(strain_id)
			ecotype_obj = ecotype_info.ecotype_id2ecotype_obj.get(ecotype_id)
			if ecotype_obj:
				lat, lon = ecotype_obj.latitude, ecotype_obj.longitude
			else:
				sys.stderr.write("Warning: Ecotype %s not in ecotype_info (fetched from stock db).\n"%ecotype_id)
				continue
			if lat and lon:
				x, y = map_data.m(lon, lat)
				if numpy.isnan(phenotype):	#2009-3-23
					linewidth=0.2
					color = 'w'
					edgecolor = 'k'
					_alpha = 0	#facecolor gets very transparent
				else:
					linewidth=0
					color = phenotype_cmap(phenotype_norm(phenotype))
					edgecolor = 'k'
					_alpha = 1
				
				#axe_map.plot([0, x], [img_y_pos, y], linestyle='--', alpha=0.2, linewidth=0.2)
				axe_map.scatter([x],[y], s=1, linewidth=linewidth, facecolor=color, zorder=10, alpha=_alpha, edgecolor=edgecolor)
				canvas_x, canvas_y = axe_map.transData.xy_tup((x,y))
				map_frame_x,  map_frame_y = axe_map_frame.transData.inverse_xy_tup((canvas_x,canvas_y))
				axe_map_frame.plot([0, map_frame_x], [img_y_pos, map_frame_y], linestyle='--', alpha=0.2, linewidth=0.2)
		#release two transformations
		axe_map.transData.thaw()  # eval the lazy objects
		axe_map.transAxes.thaw()
		axe_map_frame.transData.thaw()  # eval the lazy objects
		axe_map_frame.transAxes.thaw()	
		#set back to the original x/y lim
		axe_map.set_xlim(map_data.map_xlim)
		axe_map.set_ylim(map_data.map_ylim)
		axe_map_frame.set_xlim([0,1])
		axe_map_frame.set_ylim([0,1])
		phenotype_method = Stock_250kDB.PhenotypeMethod.get(phenotype_method_id)
		axe_map_frame.set_title('Phenotype %s %s'%(phenotype_method.id, phenotype_method.short_name))
		sys.stderr.write("Done.\n")
	
	def scale_chr_size(cls, chr_id2size, total_cumu_size, scale_range=[0,1]):
		new_chr_id2size = {}
		for chr_id, size in chr_id2size.iteritems():
			scale_min, scale_max = scale_range
			new_size = size/float(total_cumu_size)*(scale_max-scale_min) + scale_min
			new_chr_id2size[chr_id] = new_size
		return new_chr_id2size
	scale_chr_size = classmethod(scale_chr_size)
	
	def drawChromosome(cls, axe, chr_id_ls, start_pos, scale_chr_id2cumu_size, chr_y_value, scale_chr_gap, \
					verticalalignment='top', font_size=4):
		"""
		2008-10-08
			refatored out of connectSNPPCA2Chromosome()
			draw horizontal lines as chromosomes with gaps between them
		"""
		#chr_index = 1	#chr_index=0 is not interesting. it's the fake '0' chromosome.
		#axe.hlines(chr_y_value, start_pos, scale_chr_id2cumu_size[chr_id_ls[chr_index]])
		
		for chr_index in range(1, len(chr_id_ls)):
			prev_chr_id = chr_id_ls[chr_index-1]
			if chr_index==1:
				chr_start_pos = start_pos
			else:
				chr_start_pos = start_pos + scale_chr_id2cumu_size[prev_chr_id] + scale_chr_gap
			curr_chr_id = chr_id_ls[chr_index]
			chr_stop_pos = start_pos + scale_chr_id2cumu_size[curr_chr_id]
			axe.hlines(chr_y_value, chr_start_pos, chr_stop_pos, linewidth=0.5)
			axe.text((chr_start_pos+chr_stop_pos)/2., chr_y_value, 'chr %s'%curr_chr_id, \
					horizontalalignment ='center', verticalalignment=verticalalignment, size=font_size)
	
	drawChromosome = classmethod(drawChromosome)
	
	def connectSNPmatrixColumn2Chromosome(self, axe_chromosome, SNPID2PCAPosInfo, chr_id2size, chr_y_value = 0, scale_range=[0,1]):
		"""
		2008-10-08
			renamed from drawChromosome(), which is a split-out function
			return all chromosome info for connectSNPPCA2Chromosome() to use
		2008-10-07
			lay down five chromosomes
		"""
		sys.stderr.write("Drawing chromosomes ...")
		"""
		chr_size_ls = chr_id2size.values()
		min_chr_size = min(chr_size_ls)
		chr_gap = min_chr_size/3.
		"""
		chr_id2cumu_size, chr_gap, chr_id_ls = get_chr_id2cumu_size(chr_id2size)
		last_chr_id = chr_id_ls[-1]
		
		#rescale everything within [0,1]
		total_cumu_size = float(chr_id2cumu_size[last_chr_id])
		scale_min, scale_max = scale_range
		scale_chr_id2size = self.scale_chr_size(chr_id2size, total_cumu_size, scale_range)
		scale_chr_id2cumu_size = self.scale_chr_size(chr_id2cumu_size, total_cumu_size, scale_range)
		scale_chr_gap = chr_gap/total_cumu_size*(scale_max-scale_min) + scale_min
		
		start_pos = 0
		self.drawChromosome(axe_chromosome, chr_id_ls, start_pos, scale_chr_id2cumu_size, chr_y_value, \
						scale_chr_gap, verticalalignment='top', font_size=4)
		
		for snp_id in SNPID2PCAPosInfo.snp_id_ls:
			chr, pos = snp_id.split('_')
			chr = int(chr)
			pos = int(pos)
			img_x_pos = SNPID2PCAPosInfo.snp_id2img_x_pos[snp_id]
			chr_size = float(chr_id2size[chr])
			
			cumu_pos = scale_chr_id2cumu_size[chr]-scale_chr_id2size[chr]+pos/total_cumu_size
			
			axe_chromosome.plot([img_x_pos, cumu_pos], [1, chr_y_value], '--', alpha=0.2, linewidth=0.3)
		
		chr_pos_info = PassingData()
		chr_pos_info.scale_chr_id2size = scale_chr_id2size
		chr_pos_info.scale_chr_id2cumu_size = scale_chr_id2cumu_size
		chr_pos_info.scale_chr_gap = scale_chr_gap
		chr_pos_info.chr_id2size = chr_id2size
		chr_pos_info.chr_id2cumu_size = chr_id2cumu_size
		chr_pos_info.chr_id_ls = chr_id_ls
		chr_pos_info.total_cumu_size = total_cumu_size
		sys.stderr.write("Done.\n")
		return chr_pos_info
	
	@classmethod
	def findOutWhichPhenotypeColumn(cls, phenData, phenotype_method_id_set):
		"""
		2010-4-21
			main code split into SNPData.getColIndexLsGivenQuerySet() and OutputPhenotype.extractPhenotypeIDFromMethodIDName()
			
		2010-2-2
			becomes classmethod
		2009-4-30
			phenData from HaplotypeView (pylons web framework use integer as id)
		2008-11-25
			change option phenotype_method_id to phenotype_method_id_set
			return a list
			become classmethod
		2008-10-07
		"""
		col_index_to_return_ls = phenData.getColIndexLsGivenQuerySet(phenotype_method_id_set, \
																	colIDHashFunction=OutputPhenotype.extractPhenotypeIDFromMethodIDName)
		return col_index_to_return_ls
	
	@classmethod
	def assignAncestralAlleleSmallerNumber(cls, snpData, chr_pos2ancestral_allele):
		"""
		2010-2-2
			becomes classmethod
		2008-10-07
			index_type = 1: Ancestral-allele=1 if it's in chr_pos2ancestral_allele and snpData has that allele.
			index_type = 2: no ancestral in chr_pos2ancestral_allele
			index_type = 3: ancestral exists in chr_pos2ancestral_allele but doesn't appear in snpData
		"""
		snp_id2index_type = {}
		ancestral_allele_used = 1
		for snp_id, snp_allele2index in snpData.snp_id2snp_allele2index.iteritems():
			col_index = snpData.col_id2col_index[snp_id]
			snp_index2allele = snpData.snp_id2snp_allele2index[snp_id]
			if chr_pos2ancestral_allele and snp_id in chr_pos2ancestral_allele:
				ancestral_allele = chr_pos2ancestral_allele[snp_id]
				ancestral_allele_complement=number2complement[ancestral_allele]	#get the complement
				if ancestral_allele in snp_allele2index:
					if snp_allele2index[ancestral_allele]!=1:	#ancestral allele is not assigned the smallest number. swap. assume only 2 alleles
						#max_index_in_this_col = numpy.max(snpData.data_matrix[:,col_index])
						snpData.data_matrix[:,col_index] = numpy.abs(snpData.data_matrix[:,col_index]-3)	#reverse the index
					index_type = 1
				elif ancestral_allele_complement in snp_allele2index:
					if snp_allele2index[ancestral_allele_complement]!=1:
						snpData.data_matrix[:,col_index] = numpy.abs(snpData.data_matrix[:,col_index]-3)
					index_type = 1
				else:
					index_type = 3
				
			else:
				index_type = 2
			snp_id2index_type[snp_id] = index_type
		snpData.snp_id2index_type = snp_id2index_type
		return snpData
	
	@classmethod
	def getSubStrainSNPMatrix(cls, snpData, phenData, phenotype_method_id=None, phenotype_col_index=None, snp_id_ls=None, \
							chr_pos2ancestral_allele=None, need_convert_alleles2binary=True, skip_strains_with_NA_phenotype=False):
		"""
		2010-4-21 getRowIndexGivenRowID is more versatile than row_id2row_index.
			row_index = snpData.getRowIndexGivenRowID(row_id)
		2010-2-2
			becomes classmethod
		2009-3-23
			add argument skip_strains_with_NA_phenotype
		2008-12-09
			add option need_convert_alleles2binary and replace 1,2 with 0,1 to represent binary alleles.
		2008-11-30
			fix a bug
				previously, NA-phenotype strains were not thrown out due to numpy.nan!=numpy.nan.
				now use numpy.isnan() to judge whether a data point is NA.
		2008-10-07
			Ancestral-allele=1 if it exists
			Major allele=1 if no ancestral allele info
			remove all NA to do PCA (assume snpData is imputed and has no NA, so use phenData to kick out strain with NA phenotype.)
		"""
		sys.stderr.write("Getting a sub strain-snp matrix ...")
		if phenotype_col_index is None:
			sys.stderr.write("Phenotype %s not in phenotype file.\n"%phenotype_method_id)
			return None
		sub_row_id_ls = []
		
		for row_id, row_index in phenData.row_id2row_index.iteritems():
			#if phenData.data_matrix[row_index][phenotype_col_index]!=numpy.nan:	#WATCH: numpy.nan!=numpy.nan, weird!
			if skip_strains_with_NA_phenotype and numpy.isnan(phenData.data_matrix[row_index][phenotype_col_index]):	#skip rows with NA.
				continue
			sub_row_id_ls.append(row_id)
		no_of_sub_rows = len(sub_row_id_ls)
		no_of_sub_snps = len(snp_id_ls)
		sub_matrix = numpy.zeros([no_of_sub_rows, no_of_sub_snps], snpData.data_matrix.dtype)
		snp_id2snp_index2allele = {}
		snp_id2snp_allele2index = {}
		snpDataColIDType = type(snpData.col_id2col_index.keys()[0])		#2012.3.8
		for j in range(no_of_sub_snps):
			snp_id = snp_id_ls[j]
			snp_id = snpDataColIDType(snp_id)		#2012.3.8
			col_index = snpData.col_id2col_index.get(snp_id)
			if col_index is None:
				sys.stderr.write("Warning: col_index for snp: %s is None and so is skipped.\n"%snp_id)
				continue
			snp_allele2index = {}
			snp_allele2count = {}
			snp_index2allele = {}
			for i in range(no_of_sub_rows):
				row_id = sub_row_id_ls[i]
				row_index = snpData.getRowIndexGivenRowID(row_id)	#2010-4-21 getRowIndexGivenRowID is more versatile than row_id2row_index.
				if row_index is None:
					sys.stderr.write("Warning: row_index for row_id: %s is None and so is skipped.\n"%row_id)
					continue
				allele = snpData.data_matrix[row_index][col_index]
				if need_convert_alleles2binary:
					try:
						if allele not in snp_allele2index:	#if not given how to map allele1 and allele2
							snp_allele2index[allele] = len(snp_allele2index)
							snp_index = snp_allele2index[allele]
							snp_allele2count[allele] = 0
							snp_index2allele[snp_index] = allele
					except:
						sys.stderr.write('Except: %s\n'%repr(sys.exc_info()))
						traceback.print_exc()
						import pdb
						pdb.set_trace()
					#increase the counter
					snp_index = snp_allele2index[allele]
					snp_allele2count[allele] += 1
					allele = snp_index
				sub_matrix[i][j] = allele	#2008-12-05
			if len(snp_allele2count)>2:
				sys.stderr.write("Warnings: This SNP %s has %s (>2) alleles.\n"%(snp_id, len(snp_allele2count)))
			
			#make sure minor allele gets bigger number (assume newly-derived allele)
			if len(snp_allele2count)==2:
				allele1 = snp_index2allele[0]
				allele2 = snp_index2allele[1]
				if snp_allele2count[allele1]<snp_allele2count[allele2]:	#minor allele got assigned the smaller number, reverse it
					sub_matrix[:,j] = numpy.abs(sub_matrix[:,j]-1)	#reverse the index, 0,1 to 1,0
					snp_allele2index[allele1] = 1
					snp_allele2index[allele2] = 0
					snp_index2allele[1] = allele1
					snp_index2allele[0] = allele2
			
			snp_id2snp_index2allele[snp_id] = snp_index2allele
			snp_id2snp_allele2index[snp_id] = snp_allele2index
		subSNPData = SNPData(col_id_ls=snp_id_ls, row_id_ls=sub_row_id_ls, data_matrix=sub_matrix)
		subSNPData.snp_id2snp_index2allele = snp_id2snp_index2allele
		subSNPData.snp_id2snp_allele2index = snp_id2snp_allele2index
		subSNPData = cls.assignAncestralAlleleSmallerNumber(subSNPData, chr_pos2ancestral_allele)
		
		sys.stderr.write("%s X %s.\n"%(subSNPData.data_matrix.shape[0], subSNPData.data_matrix.shape[1]))
		return subSNPData
	
	@classmethod
	def getTopSNPData(cls, genome_wide_result=None, no_of_top_hits=None, snp_region_tup=[], chr_pos_ls=[],\
					sortChrPosLs=True, snpIDType=str):
		"""
		2012.3.27
			give more stdout/err message regarding the input.
		2012.3.7
			return the db_id as snp_id, instead of chr_pos
				cast it to str (argument snpIDType) format because col_id in SNPData.col_id2col_index is of str type.
		2010-4-6
			add argument sortChrPosLs. default to True which means sort genome_chr_pos_ls.
			check if chr_pos is of type str before proceedng.
				chr_pos_ls could ['1_1342','2_434'] or [(1,1342), (2,434)] or [('1',1342), ('2',434)]
		2010-2-2
			becomes classmethod
		2009-4-30
			can handle chr,pos,offset in chr_pos_ls by ignoring offset
		2008-11-14
			add option snp_region_tup.
			snp_region_tup = [start_snp_chr=None, start_snp_pos=None, stop_snp_chr=None, stop_snp_pos=None]
			if snp_region_tup is given, no_of_top_hits is defunct.
		2008-10-07
		"""
		lengthFunctor = getattr(chr_pos_ls, '__len__', None)	#2012.3.27
		if lengthFunctor:
			no_of_given_chr_pos = lengthFunctor()
		else:
			no_of_given_chr_pos = None
		sys.stderr.write("Getting Top %s SNPs within %s or from %s (chr_pos)s ..."%(no_of_top_hits, snp_region_tup, no_of_given_chr_pos))
		top_snp_data = PassingData()
		snp_id_ls = []
		score_ls = []
		if snp_region_tup:
			start_snp_chr, start_snp_pos, stop_snp_chr, stop_snp_pos = snp_region_tup
			if chr_pos_ls:	#2012.3.8 this route should not be activated. because snp is identified by chr_pos, rather than db-id
				genome_chr_pos_ls = chr_pos_ls
				if sortChrPosLs:
					genome_chr_pos_ls.sort()
				for chr_pos in genome_chr_pos_ls:
					if type(chr_pos)==str:	# 2010-4-6
						chr_pos = chr_pos.split('_')
						chr_pos = map(int, chr_pos)
					chr, pos, stopPosition = chr_pos[:3]
					chr = str(chr)
					include_this_SNP = False
					#if chr>start_snp_chr and chr<stop_snp_chr:	# no need to check position, between two chromosomes
					#	include_this_SNP = True
					if (chr==start_snp_chr and stopPosition>=start_snp_pos) or (chr==stop_snp_chr and pos<=stop_snp_pos):
						include_this_SNP = True
					
					if include_this_SNP:
						chr_pos = map(str, chr_pos)
						snp_id_ls.append('_'.join(chr_pos))	#2012.3.27 this type of SNP id only fits old-style.
						score_ls.append(1.0)
			else:
				genome_chr_pos_ls = genome_wide_result.chr_pos2index.keys()
				if sortChrPosLs:
					genome_chr_pos_ls.sort()
				for chr_pos in genome_chr_pos_ls:
					chr, pos, stopPosition = chr_pos[:3]
					if chr==start_snp_chr and chr==stop_snp_chr and stopPosition>=start_snp_pos and pos<=stop_snp_pos:
						data_obj = genome_wide_result.get_data_obj_by_chr_pos(chr, pos, stopPosition)
						snp_id_ls.append(data_obj.db_id)	#2012.3.7 use db_id instead of chr_start_stop
						#snp_id_ls.append('%s_%s_%s'%(data_obj.chromosome, data_obj.position, stopPosition))
						score_ls.append(data_obj.value)
		else:
			for i in range(no_of_top_hits):
				data_obj = genome_wide_result.get_data_obj_at_given_rank(i+1)	#rank starts from 1
				snp_id_ls.append(data_obj.db_id)	#2012.3.7 use db_id instead of chr_start_stop
				#snp_id_ls.append('%s_%s_%s'%(data_obj.chromosome, data_obj.position, data_obj.stopPosition))
				score_ls.append(data_obj.value)
		snp_id_ls = map(snpIDType, snp_id_ls)
		top_snp_data.snp_id_ls = snp_id_ls
		top_snp_data.score_ls = score_ls
		sys.stderr.write("%s SNPs.\n"%len(snp_id_ls))
		return top_snp_data
	
	def sortObjByPCA_value_ls(cls, obj_id2index, pca_value_ls, pca_range):
		"""
		2008-10-07
			position on axe_strain_pca or axe_snp_pca, the constraint axis (adjacent to axe_snp_matrix)
		"""
		tuple_to_sort_ls = []
		for obj_id, index in obj_id2index.iteritems():
			pca_value = pca_value_ls[index]
			tuple_to_sort_ls.append((pca_value, obj_id))
		tuple_to_sort_ls.sort()
		obj_id2pos = {}
		obj_id_ls = []
		pca_value_min = tuple_to_sort_ls[0][0]
		pca_value_max = tuple_to_sort_ls[-1][0]
		for pca_value, obj_id in tuple_to_sort_ls:
			if pca_range:
				pca_range_min, pca_range_max = pca_range
				pos = ((pca_value-pca_value_min)/(pca_value_max-pca_value_min))*(pca_range_max-pca_range_min)+pca_range_min
			else:
				pos = pca_value
			obj_id_ls.append(obj_id)
			obj_id2pos[obj_id] = pos
		return obj_id_ls, obj_id2pos
	sortObjByPCA_value_ls=classmethod(sortObjByPCA_value_ls)
	
	@classmethod
	def getObj2posFromPCA_value_ls(cls, obj_id2index, pca_value_ls, pca_range):
		"""
		2008-10-07
			position on axe_strain_pca or axe_snp_pca, the non-constraint axis (not adjacent to axe_snp_matrix)
		"""
		obj_id2pos = {}
		tuple_to_sort_ls = []
		for obj_id, index in obj_id2index.iteritems():
			pca_value = pca_value_ls[index]
			tuple_to_sort_ls.append((pca_value, obj_id))
		tuple_to_sort_ls.sort()
		pca_value_min = tuple_to_sort_ls[0][0]
		pca_value_max = tuple_to_sort_ls[-1][0]
		
		for obj_id, index in obj_id2index.iteritems():
			pca_value = pca_value_ls[index]
			if pca_range:
				pca_range_min, pca_range_max = pca_range
				pos = ((pca_value-pca_value_min)/(pca_value_max-pca_value_min))*(pca_range_max-pca_range_min)+pca_range_min
			else:
				pos = pca_value
			obj_id2pos[obj_id] = pos
		return 	obj_id2pos
	
	def getObj2ImgPos(cls, obj_id_ls, img_pos_range=[0., 0.9]):
		"""
		2008-10-07
			position on axe_snp_matrix (either x or y axis)
		"""
		obj_id2img_pos = {}
		no_of_objs = len(obj_id_ls)
		step = (img_pos_range[1]-img_pos_range[0])/no_of_objs
		for i in range(no_of_objs):
			obj_id = obj_id_ls[i]
			obj_id2img_pos[obj_id] = i*step
		return obj_id2img_pos, step
	getObj2ImgPos = classmethod(getObj2ImgPos)
	
	
	def getStrainID2PCAPosInfo(cls, subSNPData, pca_range=[0,1], snp_id_label_y_offset=0.95, explained_var=None, T=None):
		"""
		2012.3.27
			in the case that subSNPData has only 1 SNP/locus in it. PC matrix T will only have one column.
		2008-10-08
			fix a bug that using PC2 and PC3 to get coordinates.
		2008-10-07
			figure out strain positions in axe_strain_pca and axe_snp_matrix
		"""
		sys.stderr.write("Getting  StrainID2PCAPosInfo ... ")
		StrainID2PCAPosInfo = PassingData()
		if explained_var is None or T is None:
			T, P, explained_var = pca_module.PCA_svd(subSNPData.data_matrix, standardize=False)
		#import PCA
		#T, P, explained_var = PCA.eig(subSNPData.data_matrix)	#2008-11-30 try the PCA from pymodule
		
		#T[:,0] and T[:,1] are new positions corresponding to PC1 and PC2
		strain_id_ls, strain_id2pca_y = cls.sortObjByPCA_value_ls(subSNPData.row_id2row_index, T[:,0], pca_range)
		if T.shape[1]>1:	#2012.3.27 subSNPData has >1 SNPs, and then >1 PCs in T.
			strain_id2pca_x = cls.getObj2posFromPCA_value_ls(subSNPData.row_id2row_index, T[:,1], pca_range)
		else:	#2012.3.27 subSNPData has only one column, and thus only one PC in T.
			strain_id2pca_x = cls.getObj2posFromPCA_value_ls(subSNPData.row_id2row_index, T[:,0], pca_range)
			
		StrainID2PCAPosInfo.strain_id_ls = strain_id_ls
		StrainID2PCAPosInfo.strain_id2pca_y = strain_id2pca_y
		StrainID2PCAPosInfo.strain_id2pca_x = strain_id2pca_x
		if len(explained_var)>1:
			StrainID2PCAPosInfo.x_var = explained_var[1]
		else:
			StrainID2PCAPosInfo.x_var = explained_var[0]
		StrainID2PCAPosInfo.y_var = explained_var[0]
		StrainID2PCAPosInfo.PC_matrix = T	#2008-12-08
		StrainID2PCAPosInfo.explained_var = explained_var	#2008-12-08
		obj_id2img_pos, step = cls.getObj2ImgPos(strain_id_ls, img_pos_range=[0., snp_id_label_y_offset])
		StrainID2PCAPosInfo.strain_id2img_y_pos = obj_id2img_pos
		StrainID2PCAPosInfo.step = step
		sys.stderr.write("Done.\n")
		return StrainID2PCAPosInfo
	
	getStrainID2PCAPosInfo = classmethod(getStrainID2PCAPosInfo)
	
	def getSNPID2PCAPosInfo(self, subSNPData, pca_range=[0,1], strain_id_label_x_offset=0.95, sortObjByPCA_value=True):
		"""
		2012.3.27
			in the case that subSNPData has only 1 individual/row in it. PC matrix T will only have one column.
		2008-11-19
			add option sortObjByPCA_value to allow user decide whether to sort SNPs according PCA or not.
		2008-10-08
			fix a bug that using PC2 and PC3 to get coordinates.
		2008-10-07
			figure out SNP positions in axe_snp_pca and axe_snp_matrix
		"""
		sys.stderr.write("Getting  SNPID2PCAPosInfo ... ")
		SNPID2PCAPosInfo = PassingData()
		T, P, explained_var = pca_module.PCA_svd(numpy.transpose(subSNPData.data_matrix), standardize=False)
		#T[:,0] and T[:,1] are new positions corresponding to PC1 and PC2
		if sortObjByPCA_value:
			snp_id_ls, snp_id2pca_x = self.sortObjByPCA_value_ls(subSNPData.col_id2col_index, T[:,0], pca_range)
		else:
			snp_id_ls = subSNPData.col_id_ls
			snp_id2pca_x  = {}
			pca_range_length = pca_range[1] - pca_range[0]
			pca_range_gap = pca_range_length/float(len(snp_id_ls))
			for i in range(len(snp_id_ls)):
				snp_id = snp_id_ls[i]
				snp_id2pca_x[snp_id] = pca_range[0] + i*pca_range_gap
		if T.shape[1]>1:	#2012.3.27
			snp_id2pca_y = self.getObj2posFromPCA_value_ls(subSNPData.col_id2col_index, T[:,1], pca_range)
		else:
			snp_id2pca_y = self.getObj2posFromPCA_value_ls(subSNPData.col_id2col_index, T[:,0], pca_range)
		SNPID2PCAPosInfo.snp_id_ls = snp_id_ls
		SNPID2PCAPosInfo.snp_id2pca_x = snp_id2pca_x
		SNPID2PCAPosInfo.snp_id2pca_y = snp_id2pca_y
		SNPID2PCAPosInfo.x_var = explained_var[0]
		SNPID2PCAPosInfo.y_var = explained_var[1]
		obj_id2img_pos, step = self.getObj2ImgPos(snp_id_ls, img_pos_range=[0., strain_id_label_x_offset])
		SNPID2PCAPosInfo.snp_id2img_x_pos = obj_id2img_pos
		SNPID2PCAPosInfo.step = step
		sys.stderr.write("Done.\n")
		return SNPID2PCAPosInfo
	
	def run(self):
		"""
		2010-4-21
		2008-10-07
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		db = self.db_250k
		session = db.session
		rm = Stock_250kDB.ResultsMethod.query.filter_by(phenotype_method_id=self.phenotype_method_id_ls[0], \
													call_method_id=self.call_method_id, analysis_method_id=self.analysis_method_id).first()
		
		if rm:
			genome_wide_result = self.getResultMethodContent(rm, self.data_dir, min_MAF=0, construct_chr_pos2index=True)
		elif self.snp_region_tup:	#2010-4-21 if this region is provided, no need to get genome_wide_result
			genome_wide_result = None
		else:
			sys.stderr.write("No genome association results for phenotype_method_id=%s, call_method_id=%s, analysis_method_id=%s.\n"%\
							(self.phenotype_method_id_ls[0], self.call_method_id, self.analysis_method_id))
			sys.exit(3)
		
		initData = Association.readInData(self.phenotype_fname, self.inputFname, eigen_vector_fname=None,\
										phenotype_method_id_ls=self.phenotype_method_id_ls, test_type=1, report=0,\
										snpAlleleOrdinalConversion=False, removeUnPhenotypedSNPData=False,\
										ignore_2nd_column=True)
		snpData = initData.snpData
		phenData = initData.phenData
		which_phenotype_index_ls = initData.which_phenotype_ls
		if self.ancestral_allele_fname:
			chr_pos2ancestral_allele = PhenotypeOfAncestralDerivedAllele.get_chr_pos2ancestral_allele(self.ancestral_allele_fname)
		else:
			chr_pos2ancestral_allele = {}
		
		#phenotype_col_index = self.findOutWhichPhenotypeColumn(phenData, Set([self.phenotype_method_id]))[0]
		phenotype_col_index = which_phenotype_index_ls[0]
		
		if self.draw_snp_pca:
			self.sortSNPByPCA_value = True
			axe_y_offset3 = 0.72
			axe_y_offset4 = 0.98
		else:
			#squeeze axe_snp_pca out as it won't be drawn
			axe_y_offset3 = 0.96
			axe_y_offset4 = 0.98
		
		if self.snp_region_tup:
			#2008-11-14 region around FRI
			#snp_region_tup = [4, 200000, 4, 700000]
			start_snp, stop_snp = self.snp_region_tup.split('-') #2010-4-21
			start_snp = start_snp.split('_')
			start_snp = map(int, start_snp)
			stop_snp = stop_snp.split('_') 
			stop_snp = map(int, stop_snp)
			"""
			chr, snp_start_stop = self.snp_region_tup.split(':')
			chr = int(chr)
			snp_start_stop = snp_start_stop.split('-')
			snp_start_stop = map(int, snp_start_stop)
			snp_region_tup = [chr, snp_start_stop[0], chr, snp_start_stop[1]]
			"""
			snp_region_tup = [start_snp[0], start_snp[1], stop_snp[0], stop_snp[1]]
		else:
			snp_region_tup = []
		
		axe_y_offset1 = 0.01
		axe_y_offset2 = 0.09
		axe_height1 = axe_y_offset2 - axe_y_offset1	#height of axe_chromosome, twice height of axe_map_phenotype_legend
		axe_height2 = axe_y_offset3 - axe_y_offset2	#height of axe_strain_pca, axe_snp_matrix, axe_map
		axe_height3 = axe_y_offset4 - axe_y_offset3	#height of axe_snp_pca
		
		axe_x_offset1 = 0.02
		if self.draw_strain_pca:
			axe_x_offset2 = 0.27
			axe_x_offset3 = 0.76	#width of axe_chromosome, axe_snp_matrix, axe_snp_pca
			strain_id_label_x_offset=0.9	# space for the stain labels for the SNP matrix
			
		else:
			axe_x_offset2 = axe_x_offset1	# squeeze out axe_strain_pca
			axe_x_offset3 = 0.28
			strain_id_label_x_offset=0.8	# more space for the stain labels
		axe_snp_matrix_margin = (1-strain_id_label_x_offset)/2.0	#offset to push strain labels on even rows further right
		#the two offsets below decides where the label of strains/snps should start in axe_snp_matrix
		snp_id_label_y_offset=0.95
		
		axe_width1 = axe_x_offset2 - axe_x_offset1	#width of axe_strain_pca
		axe_width2 = axe_x_offset3 - axe_x_offset2	#width of axe_chromosome, axe_snp_matrix, axe_snp_pca
		if self.draw_map:
			axe_x_offset4 = 0.78
			axe_x_offset5 = 0.99
		else:
			axe_x_offset4 = 0.99
			axe_x_offset5 = axe_x_offset4
		axe_width3 = axe_x_offset4 - axe_x_offset3	#width of axe_phenotype
		axe_width4 = axe_x_offset5 - axe_x_offset4	#width of axe_map, axe_map_frame
		
		if genome_wide_result is None:
			chr_pos_ls = snpData.col_id_ls
		else:
			chr_pos_ls = []
		top_snp_data = self.getTopSNPData(genome_wide_result, self.no_of_top_hits, snp_region_tup=snp_region_tup, \
										chr_pos_ls=chr_pos_ls)
		subSNPData = self.getSubStrainSNPMatrix(snpData, phenData, self.phenotype_method_id_ls[0], \
											phenotype_col_index, top_snp_data.snp_id_ls, chr_pos2ancestral_allele,\
											need_convert_alleles2binary=self.need_convert_alleles2binary)
		
		ecotype_info = getEcotypeInfo(db, self.country_order_type)
		chr_id2size = get_chr_id2size(db.metadata.bind)
		
		StrainID2PCAPosInfo = self.getStrainID2PCAPosInfo(subSNPData, pca_range=[0,1], snp_id_label_y_offset=snp_id_label_y_offset)
		SNPID2PCAPosInfo = self.getSNPID2PCAPosInfo(subSNPData, pca_range=[0,1], strain_id_label_x_offset=strain_id_label_x_offset, \
														sortObjByPCA_value=self.sortSNPByPCA_value)
		
		no_of_axes_drawn = 0
		
		axe_chromosome = pylab.axes([axe_x_offset2, axe_y_offset1, axe_width2, axe_height1], frameon=False)
		axe_chromosome.set_xticks([])
		axe_chromosome.set_yticks([])
		chr_pos_info = self.connectSNPmatrixColumn2Chromosome(axe_chromosome, SNPID2PCAPosInfo, chr_id2size, \
															chr_y_value = 0, scale_range=[0,1])
		#axe_chromosome.set_xlim([0,1])
		#axe_chromosome.set_ylim([0,1])
		no_of_axes_drawn += 1
		#pylab.savefig('%s_%s.png'%(self.outputFnamePrefix, no_of_axes_drawn), dpi=400)
		

		if self.draw_strain_pca or self.draw_map:	# either needs this phenotype color map
			#legend for the ecotype color (by phenotype value) on the axe_map/axe_strain_pca
			phenotype_cmap = mpl.cm.jet
			max_phenotype = numpy.nanmax(phenData.data_matrix[:,phenotype_col_index])
			min_phenotype = numpy.nanmin(phenData.data_matrix[:,phenotype_col_index])
			phenotype_gap = max_phenotype - min_phenotype
			phenotype_jitter = phenotype_gap/10.
			phenotype_norm = mpl.colors.Normalize(vmin=min_phenotype-phenotype_jitter, vmax=max_phenotype+phenotype_jitter)
			axe_map_phenotype_legend = pylab.axes([axe_x_offset4, axe_y_offset2-axe_height1/2., axe_width4, axe_height1/2.], frameon=False)
			cb = mpl.colorbar.ColorbarBase(axe_map_phenotype_legend, cmap=phenotype_cmap,
									norm=phenotype_norm,
									orientation='horizontal')
			cb.set_label('Phenotype Legend On the Map')
		
		axe_snp_matrix = pylab.axes([axe_x_offset2, axe_y_offset2, axe_width2, axe_height2], frameon=False, sharex=axe_chromosome)
		axe_snp_matrix.set_xticks([])
		axe_snp_matrix.set_yticks([])
		#self.drawSNPMtrix(axe_snp_matrix, subSNPData, top_snp_data, StrainID2PCAPosInfo, SNPID2PCAPosInfo, \
		#				ecotype_info, strain_id_label_x_offset, snp_id_label_y_offset, draw_snp_id_label=draw_snp_id_label)
		self.drawSNPMtrix(axe_snp_matrix, subSNPData, top_snp_data, StrainID2PCAPosInfo, SNPID2PCAPosInfo, \
							ecotype_info, strain_id_label_x_offset, snp_id_label_y_offset, \
							strain_id_label_x_offset_extra=axe_snp_matrix_margin, \
							draw_snp_id_label=self.draw_snp_id_label, \
							strain_snp_label_font_size=1, \
							useAlleleToDetermineAlpha=False, \
							snp_value2color=None, snp_label_font_size=self.snp_label_font_size,\
							drawGrid=self.drawGrid)
		axe_snp_matrix.set_xlim([0,1])
		#axe_snp_matrix.set_ylim([0,1])
		no_of_axes_drawn += 1
		#pylab.savefig('%s_%s.png'%(self.outputFnamePrefix, no_of_axes_drawn), dpi=400)
		
		if self.draw_strain_pca:
			axe_strain_map = pylab.axes([axe_x_offset1, axe_y_offset3, axe_width1, axe_height3], frameon=False)
			axe_strain_pca = pylab.axes([axe_x_offset1, axe_y_offset2, axe_width1, axe_height2], frameon=False, sharey=axe_snp_matrix)
			axe_strain_map_pca_cover = pylab.axes([axe_x_offset1, axe_y_offset2, axe_width1, axe_height2+axe_height3], frameon=False, \
												sharex=axe_strain_pca)	#cover both axe_strain_map and axe_strain_pca
			axe_strain_map_pca_cover.set_yticks([])
			axe_strain_pca_xlim = [-0.05,1.05]
			axe_strain_pca_ylim = [0, 1]
			axe_strain_pca.set_xlim(axe_strain_pca_xlim)
			axe_strain_pca.set_ylim(axe_strain_pca_ylim)
			axe_strain_map_pca_cover_ylim = [0, (axe_height2+axe_height3)/axe_height2]	#set it accordingly
			axe_strain_map_pca_cover.set_ylim(axe_strain_map_pca_cover_ylim)
			
			axe_strain_map.set_xticks([])
			axe_strain_map.set_yticks([])
					
			axe_strain_pca.grid(True, alpha=0.3)
			axe_strain_pca.set_xticks([])
			axe_strain_pca.set_yticks([])
			axe_strain_pca_legend = pylab.axes([axe_x_offset1, axe_y_offset3+0.03, 0.1, axe_height3-0.03], frameon=False)
			axe_strain_pca_legend.set_xticks([])
			axe_strain_pca_legend.set_yticks([])
			country_order_name = self.country_order_type2name[self.country_order_type]
			
			self.drawStrainPCA(axe_strain_pca, axe_strain_map, axe_strain_map_pca_cover, axe_strain_pca_legend, StrainID2PCAPosInfo, \
							ecotype_info, phenData, \
						phenotype_col_index, phenotype_cmap, phenotype_norm, rightmost_x_value=axe_strain_pca_xlim[1],\
						country_order_name=country_order_name, strain_color_type=2, draw_axe_strain_map=False)
			
			axe_strain_pca.set_xlim(axe_strain_pca_xlim)
			axe_strain_pca.set_ylim(axe_strain_pca_ylim)
			axe_strain_map_pca_cover.set_ylim(axe_strain_map_pca_cover_ylim)
			no_of_axes_drawn += 1
			#pylab.savefig('%s_%s.png'%(self.outputFnamePrefix, no_of_axes_drawn), dpi=400)
		
		if self.draw_map:
			#mark ecotypes on the map colored according to phenotype
			axe_map = pylab.axes([axe_x_offset4, axe_y_offset2, axe_width4, axe_height2], frameon=False)
			#axe_map_frame is used to connect strains from axe_phenotype to dot on the axe_map (another axe due to reasons stated in drawMap())
			axe_map_frame = pylab.axes([axe_x_offset4, axe_y_offset2, axe_width4, axe_height2], frameon=False, sharey=axe_strain_pca)
			axe_map_frame.set_xticks([])
			axe_map_frame.set_yticks([])
			self.drawMap(axe_map_frame, axe_map, StrainID2PCAPosInfo, phenData, phenotype_col_index, self.phenotype_method_id_ls[0], \
						ecotype_info, phenotype_cmap, phenotype_norm)
			#axe_map.set_ylim([0,1])
			no_of_axes_drawn += 1
			#pylab.savefig('%s_%s.png'%(self.outputFnamePrefix, no_of_axes_drawn), dpi=400)
		
		no_of_phenotypes = len(self.phenotype_method_id_ls)
		gap_between_each_phenotype_axe = 0.005
		width_for_each_phenotype_axe = float(axe_width3)/no_of_phenotypes
		if width_for_each_phenotype_axe<gap_between_each_phenotype_axe:
			gap_between_each_phenotype_axe = width_for_each_phenotype_axe/3.0
		width_for_each_phenotype_axe -= gap_between_each_phenotype_axe
		for i in range(no_of_phenotypes):
			# 2010-4-21 do it one by one , rather than just take which_phenotype_index_ls
			# because which_phenotype_index_ls is not necessarily in the order of self.phenotype_method_id_ls.
			# The order is not guaranteed by phenData.getColIndexLsGivenQuerySet().
			phenotype_id = self.phenotype_method_id_ls[i]
			which_phenotype_ls = phenData.getColIndexLsGivenQuerySet(Set([phenotype_id]), \
												colIDHashFunction=OutputPhenotype.extractPhenotypeIDFromMethodIDName)
			phenotype_col_index = which_phenotype_ls[0]
			axe_x_offset = axe_x_offset3 + i*gap_between_each_phenotype_axe + i*width_for_each_phenotype_axe
			axe_phenotype = pylab.axes([axe_x_offset, axe_y_offset2, width_for_each_phenotype_axe, axe_height2], \
									frameon=True, sharey=axe_snp_matrix)
			axe_phenotype.set_yticks([])
			#if self.draw_strain_pca or self.draw_map:	# if either is toggled, it means this phenotype color map is drawn
			max_phenotype = numpy.nanmax(phenData.data_matrix[:,phenotype_col_index])
			min_phenotype = numpy.nanmin(phenData.data_matrix[:,phenotype_col_index])
			#axe_phenotype.set_xticks([real_min_phenotype, real_max_phenotype])
			#axe_phenotype.set_xticklabels([str(real_max_phenotype), str(real_max_phenotype)], fontsize=2)
			axe_phenotype.set_xticklabels([])	#no tick labels (axe_map_phenotype_legend has it already)
			axe_phenotype.set_xlabel('%.2f-%.2f'%(min_phenotype, max_phenotype), fontsize=3)
			self.drawPhenotype(axe_phenotype, StrainID2PCAPosInfo, phenData, phenotype_col_index, ecotype_info=ecotype_info,\
							phenotypeTitle=phenData.col_id_ls[phenotype_col_index].split('_')[1])
			no_of_axes_drawn += 1
			axe_phenotype.set_ylim([0,1])
		axe_snp_matrix.set_ylim([0,1])	#without this, ylim of all 3 axes are set to [0,0.9] because axe_map automatically adjust to 0-0.9
		#pylab.savefig('%s_%s.png'%(self.outputFnamePrefix, no_of_axes_drawn), dpi=400)
		
		
		if self.draw_snp_pca:
			axe_snp_pca = pylab.axes([axe_x_offset2, axe_y_offset3, axe_width2, axe_height3], frameon=False, sharex=axe_snp_matrix)
			axe_snp_pca.grid(True, alpha=0.3)
			axe_snp_pca.set_xticks([])
			axe_snp_pca.set_yticks([])
			axe_snp_pca_legend = pylab.axes([axe_x_offset3, axe_y_offset3+0.05, 0.1, 0.15], frameon=False)
			axe_snp_pca_legend.set_xticks([])
			axe_snp_pca_legend.set_yticks([])
			axe_snp_pca_ylim = [-0.2,1.3]
			self.drawSNPPCA(axe_snp_pca, axe_snp_pca_legend, SNPID2PCAPosInfo, chr_pos_info, lowest_y_value=axe_snp_pca_ylim[0], highest_y_value=axe_snp_pca_ylim[1]-0.03)
			axe_snp_pca.set_ylim(axe_snp_pca_ylim)
		
			analysis_method = Stock_250kDB.AnalysisMethod.get(self.analysis_method_id)
			axe_snp_pca.title.set_text('Top %s SNPs by %s. Layout by PC1(x-axis) vs PC2(y-axis)'%(self.no_of_top_hits, analysis_method.short_name))
		else:
			analysis_method = Stock_250kDB.AnalysisMethod.get(self.analysis_method_id)
			#axe_snp_matrix.title.set_text('SNPs from chr %s.'%(self.snp_region_tup))
		png_output_fname = '%s.png'%self.outputFnamePrefix
		pylab.savefig(png_output_fname, dpi=600)
		#pylab.savefig('%s.svg'%self.outputFnamePrefix)
		
if __name__ == '__main__':
	from matplotlib import rcParams
	rcParams['font.size'] = 6
	rcParams['legend.fontsize'] = 6
	#rcParams['text.fontsize'] = 6	#deprecated. use font.size instead
	rcParams['axes.labelsize'] = 4
	rcParams['axes.titlesize'] = 8
	rcParams['xtick.labelsize'] = 4
	rcParams['ytick.labelsize'] = 4
	
	from pymodule import ProcessOptions
	main_class = PlotGroupOfSNPs
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
