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
import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
import pylab
from pymodule import ProcessOptions, PassingData
from pymodule import utils
from pymodule.db import GenomeDB
from AbstractPlot import AbstractPlot
from pymodule.yhio.AbstractGenomeFileWalker import AbstractGenomeFileWalker

parentClass = AbstractPlot
class PlotGenomeWideData(parentClass, AbstractGenomeFileWalker):
	__doc__ = __doc__
	#
	option_default_dict = copy.deepcopy(parentClass.option_default_dict)
	#option_default_dict.update(parentClass.db_option_dict.copy())
	option_default_dict.update(parentClass.genome_db_option_dict.copy())
	
	# 2013.07.31 use AbstractGenomeFileWalker's genome_option_dict
	option_default_dict.update(AbstractGenomeFileWalker.genome_option_dict.copy())
	
	option_default_dict.update({
					('xtickInterval', 0, int): [None, '', 1, 'add a tick on the x-axis every interval within each chromosome. Default is 25 ticks throughout the genome.', ],\
					('drawCentromere', 0, int): [0, '', 0, 'toggle to plot centromere as semi-transparent band', ],\
					})
	# 2013.07.31 no need for this
	option_default_dict.pop(('positionHeader', 1, ))
	
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
	
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		parentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def getNumberOfData(self, pdata):
		"""
		2012.12.9
		"""
		return len(pdata.chr2xy_ls)
	
	def preFileFunction(self, **keywords):
		"""
		2012.12.7 setup chr2xy_ls
		"""
		pdata = parentClass.preFileFunction(self, **keywords)
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
		
		xValue = float(row[xColumnIndex])
		
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
		2013.07.30 plot chromosomes in the order of self.oneGenomeData.chr_id_ls
		2012.12.7
		
		"""
		ax = pylab.gca()
		chr_ls = pdata.chr2xy_ls.keys()
		chr_ls.sort()
		max_y = None
		min_y = None
		for chromosome in self.oneGenomeData.chr_id_ls:
			if chromosome not in pdata.chr2xy_ls:	#2013.07.30 skip this chromosome if it doesnot have data
				continue
			
			#2013.2.17 draw the centromere as transparent band
			if self.drawCentromere:
				cumuStart = self.chr_id2cumu_start.get(chromosome)
				centromere = self.chr_id2centromere.get(chromosome)
				if centromere:	#centromere is an instance of GenomeAnnotation of GenomeDB.py 
					centromereCumuStart = cumuStart + centromere.start
					centromereCumuStop = cumuStart + centromere.stop
					
					pylab.axvspan(centromereCumuStart, centromereCumuStop, facecolor='0.5', alpha=0.3)	#gray color
			"""
			else:	#2013.08.01 only draw boundaries when drawCentromere is not on
				# draw chromosome boundaries
				ax.axvline(self.chr_id2cumu_start[chromosome], linestyle='--', color='k', linewidth=0.8)
			"""
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
		AbstractGenomeFileWalker._loadGenomeStructureFromDB(self, **keywords)
		
		if hasattr(self, 'xtickInterval') and self.xtickInterval is None:
			self.xtickInterval = int(sum(self.oneGenomeData.chr_id2size.values())/25.0)	#25 ticks on average		
		
		#for marking the X-axis
		self.xtick_locs = []
		self.xtick_labels = []
		self.chrID2labelXPosition = {}	#the x-axis coordinate for each chromosome's label below x-axis
		
		for chromosome in self.oneGenomeData.chr_id_ls:
			cumuStart = self.chr_id2cumu_start.get(chromosome)
			chr_size = self.oneGenomeData.chr_id2size.get(chromosome)
			self.chrID2labelXPosition[chromosome] = cumuStart + chr_size/2	#chromosome label sits in the middle
			if self.xtickInterval:
				for j in range(cumuStart, cumuStart+ chr_size, self.xtickInterval):	#tick at each interval
					self.xtick_locs.append(j)
				for j in range(0, chr_size, self.xtickInterval):
					#label only at 5 X xtickInterval
					#if j % (5*self.xtickInterval) == 0 and j < (chr_size - 1.5*self.xtickInterval):
					#	self.xtick_labels.append(j / self.xtickInterval)
					#else:
					self.xtick_labels.append("")
		
		parentClass.setup(self, **keywords)
	
	def handleXLabel(self, **keywords):
		"""
		2012.8.6
		"""
		#add proper ticks
		if self.xtick_locs and self.xtick_labels:
			pylab.xticks(self.xtick_locs, self.xtick_labels)
		#mark the chromosome
		if self.yMin is not None:
			if (self.yMin==0 and self.yScaleLog>0):
				"""
			###2013.2.24 this is causing some error when yScaleLog =1 or 2 and yMin=0: 
			File &quot;/u/local/python/2.6/lib64/python2.6/site-packages/matplotlib/backends/backend_agg.py&quot;, line 154, in draw_text
				self._renderer.draw_text_image(font.get_image(), int(x), int(y) + 1, angle, gc)
			File &quot;/u/home/eeskin/polyacti/lib/python/numpy-1.6.2-py2.6-linux-x86_64.egg/numpy/ma/core.py&quot;, line 3795, in __int__
				raise MaskError, &apos;Cannot convert masked element to a Python int.&apos;
			numpy.ma.core.MaskError: Cannot convert masked element to a Python int.
				"""
				pass
			else:
				for chromosome, labelXPosition in self.chrID2labelXPosition.iteritems():
					pylab.text(labelXPosition, self.yMin, "%s"%(chromosome),\
						horizontalalignment='center',
						verticalalignment='top', )	#transform = ax.transAxes
		
		return parentClass.handleXLabel(self, **keywords)
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
