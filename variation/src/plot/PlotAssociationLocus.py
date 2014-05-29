#!/usr/bin/env python
"""
Examples:
	#2012.11.14
	locusID=16919;
	%s --association_locus_id $locusID --outputFname /tmp/association_locus_$locusID\.png
		-u yh --list_type_id_list 129 
		--db_passwd SECRET --hostname banyan  --gene_annotation_pickleFname /Network/Data/250k/tmp-yh/at_gene_model_pickelf
		--locusExtensionDistance 15000
	
	%s 

Description:
	2012.11.14
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])


bit_number = math.log(sys.maxint)/math.log(2)
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import copy
import pylab
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter
from pymodule import yh_matplotlib, GenomeDB, utils, SNPData, read_data
from pymodule.plot.AbstractPlot import AbstractPlot
from variation.src.mapper.AbstractVariationMapper import AbstractVariationMapper
from variation.src import Stock_250kDB
from variation.src.plot.DrawSNPRegion import DrawSNPRegion, SNPPassingData
from variation.src.common import getEcotypeInfo
from variation.src.association.Kruskal_Wallis import Kruskal_Wallis

class PlotAssociationLocus(AbstractPlot, AbstractVariationMapper, DrawSNPRegion):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractPlot.option_default_dict)
	option_default_dict[('whichColumnPlotLabel', 0, )][0] = '#SNPs in 100kb window'
	option_default_dict[('whichColumn', 0, int)][0] = 3
	option_default_dict[('xColumnPlotLabel', 0, )][0] = 'position'
	option_default_dict[('xColumnHeader', 1, )][0] = 'BIN_START'
	option_default_dict.update(AbstractVariationMapper.option_default_dict)
	option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.update({
						('association_locus_id', 1, int): [None, '', 1, 'AssociationLocus.id'],\
						('locusExtensionDistance', 1, int): [5000, '', 1, 'in drawing association landscape of one locus, how far to extend on each end'],\
						('association_landscape_type_id', 1, int): [1, '', 1, 'what kind of gwas landscape'],\
						
						('snpInfoPickleFname', 0, ): ['', '', 1, 'The file to contain pickled SNPInfo.'],\
						('snp_matrix_fname', 0, ): ['', 'I', 1, 'genotype matrix. Strain X SNP format.', ],\
						('snp_matrix_data_type', 1, int):[1, 'V', 1, 'what type of data in snp_matrix_fname? \n\
	1: SNP data (snps are subset of what db has). \n \
	2: CNV discrete call data (1 or -1); \n \
	3: CNV amplitude data; \n\
	4: arbitrary non-diallelic SNP data (SNPs are not required to be included in db). \n\
	input file has to contain start, stop position for 2 & 3. Providing central position in input works only for SNP data.'],\
						('phenotype_fname', 0, ): [None, 'N', 1, 'phenotype file, if snp_matrix_fname is given, this is needed as well.', ],\
						("list_type_id_list", 1, ): [None, '', 1, 'comma/dash separated list of Gene list type. must be in table gene_list_type.'],\
						("gene_annotation_pickleFname", 0, ): [None, '', 1, 'If the file does not exist, store a pickled gene_annotation into it.\n\
	If the file exists, load gene_annotation out of it.'],\
						
						})
	
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractVariationMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		AbstractPlot.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		
		if self.list_type_id_list:
			self.list_type_id_list = getListOutOfStr(self.list_type_id_list, data_type=int)
		else:
			self.list_type_id_list = []
	
	def getXY(self, genome_wide_result=None, **keywords):
		"""
		2012.11.14 override DrawSNPRegion.getXY
		"""
		x_ls = []
		y_ls = []
		data_obj = None
		for data_obj in genome_wide_result.data_obj_ls:
			x_ls.append(data_obj.position)
			value = data_obj.value
			y_ls.append(value)
			#append the stop of last data_obj (for CNVs)
			stop = getattr(data_obj, 'stop_position', None)
			if stop and stop!=x_ls[-1]:	#stop not null and different from start
				x_ls.append(stop)
				y_ls.append(data_obj.value)
		return x_ls, y_ls
	
	analysis_method_id_locus_type_id2color = {(1,1):'k', (32,1):'r', (1,2):'b',(32,2):'g', 
							0:'#87CEFA'}	#lightskyblue: #87CEFA	R=135 G=206	B=250 ACCESS=16436871
	
	def drawPvalue(self, axe_pvalue=None, axe_to_put_pvalue_legend=None, landscape_gwr_ls=None, \
				legend_loc='upper right', **keywords):
		"""
		2012.11.14 override DrawSNPRegion.drawPvalue
		"""
		sys.stderr.write("\t Drawing pvalues  ...")
		legendKey2Data = {}
		
		for genome_wide_result in landscape_gwr_ls:
			x_ls, y_ls = self.getXY(genome_wide_result=genome_wide_result)
			if x_ls and y_ls:
				color_key = (genome_wide_result.db_entry.analysis_method.id, genome_wide_result.db_entry.locus_type.id)
				plotObject = axe_pvalue.plot(x_ls, y_ls, '-', linewidth=0.6,\
									color=self.analysis_method_id_locus_type_id2color.get(color_key, '#87CEFA'), alpha=0.4)
				plotObject = plotObject[0]
				#legend_ls.append(self.analysis_method_id2name[analysis_method_id])
				legendKey = '%s,%s'%(genome_wide_result.db_entry.analysis_method.short_name, genome_wide_result.db_entry.locus_type.short_name)
				if legendKey not in legendKey2Data:
					legendKey2Data[legendKey] = plotObject
		axe_pvalue.set_ylabel(r'-log(P-value)')
		
		plotObject_ls = []
		legend_ls = []
		for legendKey, data in legendKey2Data.iteritems():
			legend_ls.append(legendKey)
			plotObject_ls.append(data)
			
		legend = axe_to_put_pvalue_legend.legend(plotObject_ls, legend_ls, loc=legend_loc, fancybox=True, handlelength=0.02)	#cut the legend length to 0.02, default 0.05 (5% of x-axis).
		frame  = legend.get_frame()
		frame.set_alpha(0.60)	#2010-4-15 set the legend box semi-transparent
		
		sys.stderr.write("Done.\n")
		#return legend_ls
	
	def draw(self, centralLocus=None, candidate_gene_set=None, gene_annotation=None,
			snp_info=None, \
			landscape_gwr_ls=None, \
			LD_info=None, min_distance=40000, \
			label_gene=False,\
			snpData=None, phenData=None, ecotype_info=None,\
			snpData_before_impute=None, snp_matrix_data_type=1, \
			drawMap=False, drawStrainPCA=False, markSNPInGenes=True, drawGrid=False, phenotypeDrawType=1,\
			phenotype_method_id=None,\
			**keywords):
		"""
		2012.11.14
			# draw the gwas-landscape in the gwas-axe
			# draw the gene model
			# draw the SNP/deletion matrix
		"""
		
		if snp_info:
			if getattr(centralLocus, 'stop', None):
				snps_within_this_region = self.findSNPsInRegion(snp_info, chromosome=centralLocus.chromosome, start=centralLocus.position, \
													stop=centralLocus.stop)
				#snps_within_this_region_snpData = self.findSNPsInRegion(snpData, chromosome=centralLocus.chromosome, start=centralLocus.position, \
				#									stop=centralLocus.stop)
			else:
				snps_within_this_region = self.findSNPsInRegion(snp_info, chromosome=centralLocus.chromosome, start=centralLocus.position-min_distance, \
															stop=centralLocus.stop+min_distance, center_snp_position=centralLocus.position)
			if len(snps_within_this_region.chr_pos_ls)==0:
				return None
		else:
			snps_within_this_region = centralLocus
		snps_within_this_region_snpData = snps_within_this_region
		
		pylab.clf()
		#fig = pylab.figure()
		axe_y_offset1 = 0.05	#y_offset for axe_LD, axe_strain_pca, axe_phenotype, axe_map
		if (snpData and phenData and phenotype_method_id) or LD_info:
			axe_y_offset2 = 0.6
			axe_height1 = axe_y_offset2 - axe_y_offset1	#height of axe_LD or axe_snp_matrix
			axe_y_offset3 = 0.7
		else:
			axe_y_offset2 = 0.06
			axe_height1 = axe_y_offset2 - axe_y_offset1	#height of axe_LD or axe_snp_matrix (squeezed to almost nothing)
			axe_y_offset3 = 0.15
		axe_height2 = axe_y_offset3 - axe_y_offset2	#height of axe_gene_model
		axe_y_offset4 = 0.95
		axe_height3 = axe_y_offset4 - axe_y_offset3	#height of axe_pvalue
		
		axe_x_offset1 = 0.06	#
		if drawStrainPCA:	#2010-4-19
			axe_x_offset2 = axe_x_offset1 + 0.2
		else:
			axe_x_offset2 = axe_x_offset1
		axe_width1 = axe_x_offset2 - axe_x_offset1	#width of axe_strain_pca
		if drawMap:
			axe_x_offset3 = 0.77
			axe_width2 = axe_x_offset3 - axe_x_offset2	#width of axe_pvalue, axe_LD, or axe_snp_matrix
			axe_x_offset4 = 0.79 
			axe_width3 = axe_x_offset4 - axe_x_offset3	#width of axe_phenotype
		else:	#squeeze out the map panel
			axe_x_offset3 = 0.95
			axe_width2 = axe_x_offset3 - axe_x_offset2	#width of axe_pvalue, axe_LD, or axe_snp_matrix
			axe_x_offset4 = 0.97
			axe_width3 = axe_x_offset4 - axe_x_offset3	#width of axe_phenotype
		axe_x_offset5 = 0.99
		axe_width4 = axe_x_offset5 - axe_x_offset4	#width of axe_map, axe_map_frame
		no_of_axes_drawn = 0
		
		axe_pvalue = pylab.axes([axe_x_offset2, axe_y_offset3, axe_width2, axe_height3], frameon=False)	#left gap, bottom gap, width, height, axes for pvalue, gene models
		axe_pvalue.grid(True, alpha=0.3)
		axe_pvalue.set_xticklabels([])	#remove xtick labels on axe_pvalue because axe_LD's xtick labels cover this.
		axe_LD_center_SNP = pylab.twinx()	#axes for LD with center SNP, copy axe_pvalue's
		axe_LD_center_SNP.set_xticklabels([])
		axe_gene_model = pylab.axes([axe_x_offset2, axe_y_offset2, axe_width2, axe_height2], frameon=False, sharex=axe_pvalue)
		#axe_gene_model.set_xticks([])	#this will set axe_pvalue's xticks off as well because the x-axis is shared.
		axe_gene_model.set_yticks([])
		
		# use actual loci data to restrict the region to be shown
		#snp_region_tup = [snps_within_this_region_snpData.chr_pos_ls[0][0], snps_within_this_region_snpData.chr_pos_ls[0][1],\
		#				snps_within_this_region_snpData.chr_pos_ls[-1][0], snps_within_this_region_snpData.chr_pos_ls[-1][1]]
		# 2012.3.8 region decided by snps_within_this_region_snpData
		snp_region_tup = [snps_within_this_region_snpData.chromosome, snps_within_this_region_snpData.start,\
						snps_within_this_region_snpData.chromosome, snps_within_this_region_snpData.stop]
		axe_snp_matrix_margin = abs(snp_region_tup[3]-snp_region_tup[1])/15.	#offset to push strain labels on even rows further right
		if LD_info:
			axe_LD = pylab.axes([axe_x_offset2, axe_y_offset1, axe_width2, axe_height1], frameon=False)	#axes for LD
			axe_LD_legend = pylab.axes([axe_x_offset3-0.1, axe_y_offset1+0.03, 0.1, 0.13], frameon=False)	#axes for the legend of LD
			axe_LD_legend.set_xticks([])
			axe_LD_legend.set_yticks([])
			axe_to_put_pvalue_legend = axe_LD
			legend_loc = 'lower left'
			axe_pvalue_xlim = [snp_region_tup[1]-axe_snp_matrix_margin, snp_region_tup[3]+axe_snp_matrix_margin]
		elif snpData and phenData and phenotype_method_id:
			phenotype_col_index = self.findOutWhichPhenotypeColumn(phenData, set([phenotype_method_id]))[0]
			genome_wide_result = landscape_gwr_ls[0]
			if snp_matrix_data_type==1:
				chr_pos_ls = []
			elif snp_matrix_data_type==3:
				chr_pos_ls = snpData.chr_pos2index.keys()
				#2008-12-08 for CNV probes. use snpData.chr_pos2index.keys() to locate top_snp_data because here snpData doesn't match genome_wide_result.
			else:	#2012.3.27 for other types, genome_wide_result matches snpData. 
				chr_pos_ls = []
			top_snp_data = self.getTopSNPData(genome_wide_result, None, snp_region_tup, chr_pos_ls=chr_pos_ls)
			if len(top_snp_data.snp_id_ls)==0:	#2012.3.23 no plot
				return None
			if snp_matrix_data_type==3:	#2009-3-23 for CNV amplitude data, don't convert alleles into binary 0/1=major/minor form and use allele/amplitude to determine alpha
				need_convert_alleles2binary = False
				useAlleleToDetermineAlpha = True
			elif snp_matrix_data_type==4:
				#2009-3-27, for arbitrary non-diallelic SNP matrix
				need_convert_alleles2binary = False
				useAlleleToDetermineAlpha = False
			else:
				need_convert_alleles2binary = True
				useAlleleToDetermineAlpha = False
			
			subSNPData = self.getSubStrainSNPMatrix(snpData, phenData, phenotype_method_id=phenotype_method_id, \
												phenotype_col_index=phenotype_col_index, \
												snp_id_ls=top_snp_data.snp_id_ls, \
												need_convert_alleles2binary=need_convert_alleles2binary)	#2009-3-23 last argument is for CNV intensity matrix
			snp_value2color = None
			if snp_matrix_data_type==4:
				##2009-3-27 it's for SNP matrix inferred from raw sequences, might have >2 alleles, heterozygous calls, deletions etc.
				from DrawSNPMatrix import DrawSNPMatrix
				subSNPData.data_matrix = DrawSNPMatrix.transformMatrixIntoTwoAllelesAndHetero(subSNPData.data_matrix)
				snp_value2color = self.snp_value2five_color
			
			#the two offsets below decides where the label of strains/snps should start in axe_snp_matrix
			last_chr_pos = snps_within_this_region_snpData.chr_pos_ls[-1]
			strain_id_label_x_offset=snps_within_this_region_snpData.chr_pos2adjacent_window[last_chr_pos][1]	#right next to the rightmost SNP
			snp_id_label_y_offset=0.95
			
			StrainID2PCAPosInfo = self.getStrainID2PCAPosInfo(subSNPData, pca_range=[0,1], snp_id_label_y_offset=snp_id_label_y_offset)
			
			#fake one SNPID2PCAPosInfo only for drawSNPMtrix()
			SNPID2PCAPosInfo = PassingData(step=None, snp_id2img_x_pos={})
			for locus_id, adjacent_window in snps_within_this_region_snpData.locus_id2adjacent_window.iteritems():
				#2012.3.8 use locus_id, rather than chr_pos,
				#chr_pos = map(str, chr_pos)
				#snp_id = '_'.join(chr_pos)
				SNPID2PCAPosInfo.snp_id2img_x_pos[str(locus_id)] = adjacent_window
			
			axe_snp_matrix = pylab.axes([axe_x_offset2, axe_y_offset1, axe_width2, axe_height1], frameon=False)
			#axe_snp_matrix.set_xticks([])
			axe_snp_matrix.set_yticks([])
			self.drawSNPMtrix(axe_snp_matrix, subSNPData, top_snp_data, StrainID2PCAPosInfo, SNPID2PCAPosInfo, \
							ecotype_info, strain_id_label_x_offset, snp_id_label_y_offset, \
							strain_id_label_x_offset_extra=axe_snp_matrix_margin, \
							draw_snp_id_label=False, snpData_before_impute=snpData_before_impute, \
							strain_snp_label_font_size=1, \
							useAlleleToDetermineAlpha=useAlleleToDetermineAlpha, \
							snp_value2color=snp_value2color, drawGrid=drawGrid)	#2008-11-14 turn draw_snp_id_label off
			#axe_snp_matrix.set_xlim([0,1])
			#axe_snp_matrix.set_ylim([0,1])
			no_of_axes_drawn += 1
			#pylab.savefig('%s_%s.png'%(self.output_fname_prefix, no_of_axes_drawn), dpi=400)
			
			axe_snp_matrix.set_ylim([0,1])	#without this, ylim of all 3 axes are set to [0,0.9] because axe_map automatically adjust to 0-0.9
			#pylab.savefig('%s_%s.png'%(self.output_fname_prefix, no_of_axes_drawn), dpi=400)
			
			axe_to_put_pvalue_legend = axe_pvalue	#axe_LD is gone. put legend into axe_pvalue itself.
			legend_loc = 'upper right'
			axe_LD = None
			axe_LD_legend = None
			
			axe_pvalue_xlim = [snp_region_tup[1]-axe_snp_matrix_margin, snp_region_tup[3]+axe_snp_matrix_margin*2]
		else:
			axe_to_put_pvalue_legend = axe_pvalue	#axe_LD is gone. put legend into axe_pvalue itself.
			legend_loc = 'upper right'
			axe_LD = None
			axe_LD_legend = None
			axe_pvalue_xlim = [snp_region_tup[1]-axe_snp_matrix_margin, snp_region_tup[3]+axe_snp_matrix_margin*2]
			
		fig_title = 'chr %s: %s-%s'%(centralLocus.chromosome, centralLocus.start, centralLocus.stop)
		axe_pvalue.title.set_text(fig_title)	#main title using this snp.
		
		#2009-5-29 temporarily snps_within_this_region is replaced by snps_within_this_region_snpData below
		self.drawPvalue(axe_pvalue=axe_pvalue, axe_to_put_pvalue_legend=axe_to_put_pvalue_legend, \
					landscape_gwr_ls=landscape_gwr_ls, legend_loc=legend_loc)
		gene_position_cycle = 5
		base_y_value = 1
		gene_width=0.8
		gene_box_text_gap = min_distance*2*0.005
		
		skip_gene_model = False
		"""
		if len(snps_within_this_region.chr_pos_ls)>0:
			_snps_within_this_region = snps_within_this_region
		elif len(snps_within_this_region_snpData.chr_pos_ls)>0:
			_snps_within_this_region = snps_within_this_region_snpData
		else:
			skip_gene_model = True
		"""
		if not skip_gene_model:
			return_data = self.drawGeneModel(axe_gene_model, region=centralLocus, gene_annotation=gene_annotation, \
									candidate_gene_set=candidate_gene_set, gene_width=gene_width, gene_position_cycle=gene_position_cycle, \
									base_y_value=base_y_value, gene_box_text_gap=gene_box_text_gap,\
									label_gene=label_gene)
		matrix_of_gene_descriptions = return_data.matrix_of_gene_descriptions
		gene_model_min_y = base_y_value-gene_width
		gene_model_max_y = gene_position_cycle + base_y_value -1 + gene_width	#"-1" because genes never sit on y=gene_position_cycle + base_y_value
		
		#adjust x, y limits and etc
		axe_pvalue.set_xlim(axe_pvalue_xlim)
		axe_pvalue_ylim = axe_pvalue.get_ylim()
		axe_pvalue.set_ylim((0, axe_pvalue_ylim[1]))	#set axe_pvalue to 0 to sit right above axe_gene_model
		axe_gene_model.set_ylim((gene_model_min_y, gene_model_max_y))	#LD panel right under gene models
		
		if snpData:
			axe_snp_matrix.set_xlim(axe_pvalue.get_xlim())
		
		if self.debug:
			pylab.show()
		sys.stderr.write("Done.\n")
	
	def loadDataStructure(self, db_250k=None, association_locus_id=None, association_landscape_type_id=None, \
						locusExtensionDistance=5000,\
						data_dir=None, list_type_id_list=None, gene_annotation_pickleFname=None, \
						snpInfoPickleFname=None, locus_type_id=1, snp_matrix_fname=None, snp_matrix_data_type=None, \
						phenotype_fname=None):
		"""
		2012.11.14
		"""
		sys.stderr.write("Fetching GWAS landscape for association-locus %s, landscape type %s ..."%(association_locus_id, association_landscape_type_id))
		# fetch the associationLocus
		associationLocus = Stock_250kDB.AssociationLocus.get(association_locus_id)
		associationLandscapeType = Stock_250kDB.AssociationLandscapeType.get(association_landscape_type_id)
		
		# fetch all result-peaks
		landscape_gwr_ls = []
		# fetch landscape within this interval
		start = max(1, associationLocus.start-locusExtensionDistance)
		stop = associationLocus.stop + locusExtensionDistance
		pd = PassingData(min_MAF=associationLandscapeType.min_MAF, data_dir=data_dir, \
						need_chr_pos_ls=0, chromosome=associationLocus.chromosome, \
						start=start, stop=stop, report=False)	#report controls whether getResultMethodContent() will report progress.
		association_landscape_id_set = set()
		
		for association_peak in associationLocus.association_peak_ls:
			association_landscape = db_250k.getAssociationLandscape(result_id=association_peak.result_id, association_landscape_type_id=associationLandscapeType.id)
			if association_landscape and association_landscape.id not in association_landscape_id_set:
				association_landscape_id_set.add(association_landscape.id)
				genome_wide_result = db_250k.getResultMethodContent(association_landscape=association_landscape, data_dir=data_dir, \
												construct_chr_pos2index=True, pdata=pd)
				landscape_gwr_ls.append(genome_wide_result)
				sys.stderr.write(" %s%s "%('\x08'*80, len(landscape_gwr_ls)))
		sys.stderr.write("%s landscapes.\n"%(len(landscape_gwr_ls)))
		
		centralLocus = SNPPassingData(chromosome=associationLocus.chromosome, position=start, \
						snps_id=associationLocus.id, start=start, stop=stop,
						fileNamePrefix="")
		
		LD_info = None
		gene_annotation = DrawSNPRegion.dealWithGeneAnnotation(gene_annotation_pickleFname)
		if snpInfoPickleFname:
			snp_info = db_250k.dealWithSNPInfo(snpInfoPickleFname, locus_type_id=locus_type_id)	#2012.3.8
		else:
			snp_info = None
		
		candidate_gene_set = set()
		if list_type_id_list:
			for list_type_id in list_type_id_list:
				candidate_gene_list = db_250k.getGeneList(list_type_id)
				candidate_gene_set |= set(candidate_gene_list)
		
		if snp_matrix_fname and phenotype_fname:
			if snp_matrix_data_type==3:
				matrix_data_type=float		#2009-3-23 for CNV amplitude file
			else:
				matrix_data_type=int
			snpData = SNPData(input_fname=snp_matrix_fname, turn_into_integer=1, turn_into_array=1, ignore_2nd_column=1,\
							matrix_data_type=matrix_data_type)
			if snpData.data_matrix is None:
				sys.stderr.write("Error. snpData.data_matrix is None.\n")
				sys.exit(3)
			header_phen, strain_acc_list_phen, category_list_phen, data_matrix_phen = read_data(phenotype_fname, turn_into_integer=0)
			
			phenData = SNPData(header=header_phen, strain_acc_list=snpData.strain_acc_list, data_matrix=data_matrix_phen)
			#row label is that of the SNP matrix, because the phenotype matrix is gonna be re-ordered in that way
			
			phenData.data_matrix = Kruskal_Wallis.get_phenotype_matrix_in_data_matrix_order(snpData.row_id_ls, \
																		strain_acc_list_phen, phenData.data_matrix)
			#tricky, using strain_acc_list_phen
			
			#2008-12-05 fake a snp_info for findSNPsInRegion
			DrawSNPRegion.construct_chr_pos2index_forSNPData(snpData, snp_info=snp_info)
			ecotype_info = getEcotypeInfo(db_250k)
		else:
			snpData = None
			phenData = None
			ecotype_info = None
		
		return_data = PassingData(associationLocus=associationLocus, associationLandscapeType=associationLandscapeType, \
								landscape_gwr_ls=landscape_gwr_ls, \
								gene_annotation=gene_annotation, snp_info=snp_info, LD_info=LD_info, \
								candidate_gene_set=candidate_gene_set, snpData=snpData, phenData=phenData,\
								ecotype_info=ecotype_info, centralLocus=centralLocus)
		return return_data
		
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		self.setup()
		
		pdata = self.loadDataStructure(db_250k=self.db_250k, association_locus_id=self.association_locus_id, \
						association_landscape_type_id=self.association_landscape_type_id, \
						locusExtensionDistance=self.locusExtensionDistance,\
						data_dir=self.data_dir, list_type_id_list=self.list_type_id_list, \
						gene_annotation_pickleFname=self.gene_annotation_pickleFname, \
						snpInfoPickleFname=self.snpInfoPickleFname, locus_type_id=1, \
						snp_matrix_fname=self.snp_matrix_fname, snp_matrix_data_type=self.snp_matrix_data_type, \
						phenotype_fname=self.phenotype_fname)
		
		self.draw(centralLocus=pdata.centralLocus, candidate_gene_set=pdata.candidate_gene_set, \
				gene_annotation=pdata.gene_annotation, snp_info=pdata.snp_info, landscape_gwr_ls=pdata.landscape_gwr_ls, \
				LD_info=pdata.LD_info, \
				min_distance=10000, label_gene=False, snpData=pdata.snpData, \
				phenData=pdata.phenData, ecotype_info=pdata.ecotype_info, snpData_before_impute=None, \
				snp_matrix_data_type=self.snp_matrix_data_type, \
				drawMap=False, drawStrainPCA=False, markSNPInGenes=False)
		
		#self.handleTitle()
		#self.handleXLabel()
		#self.handleYLabel()
		
		self.reduce()
		
if __name__ == '__main__':
	main_class = PlotAssociationLocus
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()