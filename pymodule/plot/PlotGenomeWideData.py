#!/usr/bin/env python
"""
Examples:
	# 2012 input is hdf5 association table
	secret=...
	%s --xColumnHeader start --whichColumnHeader score
		-i /Network/Data/250k/db/results/type_1/56_call6_pheno38_ana1.h5
		-o /Network/Data/250k/db/results/type_1/56_call6_pheno38_ana1.h5.png
		--drivername postgresql --hostname localhost --dbname vervetdb --schema stock_250k --db_user yh --db_passwd $secret
		--genome_drivername=postgresql --genome_hostname=localhost --genome_dbname=vervetdb --genome_schema=genome --genome_db_user=yh
		--genome_db_passwd=$secret --tax_id=3702 --minNoOfTotal 1 --figureDPI 200 --h5TableName association
		--logY 2
		--drawCentromere
		#--chromosomeHeader chromosome
	
	#2013.1.7 input is plain tsv file (assuming it has columns, chromosome, start, score.)
	%s   --xColumnHeader start --whichColumnHeader score
		-i /Network/Data/250k/db/results/type_1/7314_call75_pheno2_ana32.tsv
		-o /Network/Data/250k/db/results/type_1/7314_call75_pheno2_ana32.tsv.png
		--drivername postgresql --hostname localhost --dbname vervetdb --schema stock_250k --db_user yh --db_passwd $secret
		--genome_drivername=postgresql --genome_hostname=localhost --genome_dbname=vervetdb --genome_schema=genome
		--genome_db_user=yh --genome_db_passwd=$secret --tax_id=3702 --minNoOfTotal 1 --figureDPI 200  --inputFileFormat 1
		--logY 2
	

Description:
	2012.12.11 a program that makes generic GWAS-like plots
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])


#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:	   #64bit
#	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
#	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
#else:   #32bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import copy
import pylab
from pymodule import ProcessOptions, PassingData
from pymodule import utils
from pymodule.db import GenomeDB
from AbstractPlot import AbstractPlot

class PlotGenomeWideData(AbstractPlot):
	__doc__ = __doc__
	#						
	option_default_dict = copy.deepcopy(AbstractPlot.option_default_dict)
	option_default_dict.update(AbstractPlot.db_option_dict.copy())
	option_default_dict.update(AbstractPlot.genome_db_option_dict.copy())
	
	# change
	option_default_dict[('xColumnPlotLabel', 0, )][0] = 'genome position'
	# change default of file-format
	option_default_dict[('inputFileFormat', 0, int)][0] = 2
	option_default_dict[('minNoOfTotal', 1, int)][0] = 1
	option_default_dict[('defaultFigureWidth', 1, float)][0] = 40
	option_default_dict[('defaultFigureHeight', 1, float)][0] = 5
	option_default_dict[('xmargin', 0, float)][0] = 0.02
	option_default_dict[('ymargin', 0, float)][0] = 0.1
	option_default_dict[('plotLeft', 0, float)][0] = 0.02
	option_default_dict[('plotRight', 0, float)][0] = 0.98
	option_default_dict[('defaultFontLabelSize', 1, int)][0] = 16
	
	option_default_dict.update({
					('chromosomeHeader', 1, ): ['chromosome', '', 1, 'header of the column that designates the chromosome' ],\
					('tax_id', 0, int): [3702, '', 1, 'taxonomy ID of the organism from which to retrieve the chromosome info', ],\
					('xtickInterval', 0, int): [1000000, '', 1, 'add a tick on the x-axis every this interval within each chromosome', ],\
					('drawCentromere', 0, int): [0, '', 0, 'toggle to plot centromere as semi-transparent band', ],\
					})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractPlot.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def getNumberOfData(self, pdata):
		"""
		2012.12.9
		"""
		return len(pdata.chr2xy_ls)
	
	def preFileFunction(self, **keywords):
		"""
		2012.12.7 setup chr2xy_ls
		"""
		pdata = AbstractPlot.preFileFunction(self, **keywords)
		pdata.chr2xy_ls = {}
		return pdata
	
	def processRow(self, row=None, pdata=None):
		"""
		2012.12.7 add data
		"""
		col_name2index = getattr(pdata, 'col_name2index', None)
		
		x_ls = getattr(pdata, 'x_ls', None)
		y_ls = getattr(pdata, 'y_ls', None)
		
		chromosomeIndex = col_name2index.get(self.chromosomeHeader)
		xColumnIndex = col_name2index.get(self.xColumnHeader)
		yColumnIndex = col_name2index.get(self.whichColumnHeader, self.whichColumn)
		
		chromosome = row[chromosomeIndex]
		if chromosome not in pdata.chr2xy_ls:
			pdata.chr2xy_ls[chromosome] = [[],[]]
		
		xValue = row[xColumnIndex]
		
		#stopPosition = row.stop
		#width = abs(stopPosition - xValue)
		
		#add cumu start to xValue so that it'll be part of the whole genome
		cumuStart = self.chr_id2cumu_start.get(chromosome)
		xValue += cumuStart
		xValue = self.processValue(value=xValue, processType=self.logX)
		pdata.chr2xy_ls[chromosome][0].append(xValue)
		
		yValue = row[yColumnIndex]
		yValue = self.processValue(value=yValue, processType=self.logY)
		pdata.chr2xy_ls[chromosome][1].append(yValue)
		
		return 1
	
	def plot(self, x_ls=None, y_ls=None, pdata=None):
		"""
		2012.12.7
		
		"""
		ax = pylab.gca()
		chr_ls = pdata.chr2xy_ls.keys()
		chr_ls.sort()
		max_y = None
		min_y = None
		for chromosome in chr_ls:
			# draw chromosome boundaries
			ax.axvline(self.chr_id2cumu_start[chromosome], linestyle='--', color='k', linewidth=0.8)
			#2013.2.17 draw the centromere as transparent band
			if self.drawCentromere:
				cumuStart = self.chr_id2cumu_start.get(chromosome)
				centromere = self.chr_id2centromere.get(chromosome)
				if centromere:	#centromere is an instance of GenomeAnnotation of GenomeDB.py 
					centromereCumuStart = cumuStart + centromere.start
					centromereCumuStop = cumuStart + centromere.stop
					
					pylab.axvspan(centromereCumuStart, centromereCumuStop, facecolor='0.5', alpha=0.3)	#gray color
			
			
			x_ls, y_ls = pdata.chr2xy_ls[chromosome][:2]
			
			self.setGlobalMinVariable(extremeVariableName='xMin', givenExtremeValue=min(x_ls))
			self.setGlobalMaxVariable(extremeVariableName='xMax', givenExtremeValue=max(x_ls))
			self.setGlobalMinVariable(extremeVariableName='yMin', givenExtremeValue=min(y_ls))
			self.setGlobalMaxVariable(extremeVariableName='yMax', givenExtremeValue=max(y_ls))
			
			if x_ls and y_ls:
				ax.plot(x_ls, y_ls, '.', markeredgewidth=0, alpha=0.6)	#markersize=5, 
		
		"""
		if drawBonferroni:
			#draw the bonferroni line
			bonferroni_value = -math.log10(0.01/len(genome_wide_result.data_obj_ls))
			ax.axhline(bonferroni_value, linestyle='--', color='k', linewidth=0.8)
		"""
	def setup(self, **keywords):
		"""
		2012.10.15
			run before anything is run
		"""
		#without commenting out db_vervet connection code. schema "genome" wont' be default path.
		#db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
		#				password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
		db_genome = GenomeDB.GenomeDatabase(drivername=self.genome_drivername, username=self.genome_db_user,
						db_passwd=self.genome_db_passwd, hostname=self.genome_hostname, dbname=self.genome_dbname, \
						schema=self.genome_schema)
		db_genome.setup(create_tables=False)
		#chrOrder=1 is to order chromosomes alphabetically
		self.oneGenomeData = db_genome.getOneGenomeData(tax_id=self.tax_id, chr_gap=0, chrOrder=1, sequence_type_id=1)
		self.chr_id2cumu_start = self.oneGenomeData.chr_id2cumu_start
		#2013.2.18 get centromere locations
		self.chr_id2centromere = self.oneGenomeData.chr_id2centromere
		
		#for marking the X-axis
		self.xtick_locs = []
		self.xtick_labels = []
		self.chrID2labelXPosition = {}	#the x-axis coordinate for each chromosome's label below x-axis
		
		for chromosome in self.oneGenomeData.chr_id_ls:
			cumuStart = self.chr_id2cumu_start.get(chromosome)
			chr_size = self.oneGenomeData.chr_id2size.get(chromosome)
			self.chrID2labelXPosition[chromosome] = cumuStart + chr_size/2	#chromosome label sits in the middle

			for j in range(cumuStart, cumuStart+ chr_size, self.xtickInterval):	#tick at each interval
				self.xtick_locs.append(j)
			for j in range(0, chr_size, self.xtickInterval):
				#label only at 5 X xtickInterval
				if j % (5*self.xtickInterval) == 0 and j < (chr_size - 1.5*self.xtickInterval):
					self.xtick_labels.append(j / self.xtickInterval)
				else:
					self.xtick_labels.append("")
		
		AbstractPlot.setup(self, **keywords)
	
	def handleXLabel(self, **keywords):
		"""
		2012.8.6
		"""
		#add proper ticks
		if self.xtick_locs and self.xtick_labels:
			pylab.xticks(self.xtick_locs, self.xtick_labels)
		#mark the chromosome
		if self.yMin is not None:
			for chromosome, labelXPosition in self.chrID2labelXPosition.iteritems():
				pylab.text(labelXPosition, self.yMin, "Chr %s"%(chromosome),\
					horizontalalignment='center',
				 	verticalalignment='center', )	#transform = ax.transAxes
		
		return AbstractPlot.handleXLabel(self, **keywords)
		"""
		if getattr(self, 'xColumnPlotLabel', None):
			xlabel = self.xColumnPlotLabel
		else:
			xlabel = getattr(self, "xColumnHeader", "")
		pylab.xlabel(xlabel)
	
		plt.text(offset / 2, maxScore + scoreRange * 0.1, 'Position')
		return xlabel
		"""


if __name__ == '__main__':
	main_class = PlotGenomeWideData
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()