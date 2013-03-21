#!/usr/bin/env python
"""
Examples:
	%s
	
	%s  -s 1  -o locus_frequency_vs_min_score.png /tmp/count_association_locus.h5
		 --logY --positiveLog
	

Description:
	2012/11/25
		input is output of AssociationPeak2AssociationLocus.py.
		Plot the number of results in which one association locus is significant on Y-axis as a vertical bar.
			genomic position on X-axis. It's like a GWAS manhattan plot with bars replacing the dots.
		Multiple input will be drawn in different colors.
		Preferably, this program should be called with --logY and --positiveLog. 
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
import csv
import numpy, random, pylab
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter
from pymodule import yh_matplotlib, utils, GenomeDB
from pymodule.plot.PlotGenomeWideData import PlotGenomeWideData
from PlotAssociationLocusFrequencyVsAssociationThreshold import PlotAssociationLocusFrequencyVsAssociationThreshold
class PlotAssociationLocusFrequencyOnGenome(PlotGenomeWideData, PlotAssociationLocusFrequencyVsAssociationThreshold):
	__doc__ = __doc__
	#						
	option_default_dict = PlotAssociationLocusFrequencyVsAssociationThreshold.option_default_dict.copy()
	option_default_dict.update(PlotAssociationLocusFrequencyVsAssociationThreshold.db_option_dict.copy())
	option_default_dict.update(PlotGenomeWideData.option_default_dict.copy())
	
	# change
	option_default_dict[('xColumnHeader', 1, )][0] = 'start'
	option_default_dict[('xColumnPlotLabel', 0, )][0] = 'genome position'
	option_default_dict[('whichColumnHeader', 0, )][0] = 'no_of_results'
	option_default_dict[('whichColumnPlotLabel', 0, )][0] = 'numberOfResults'
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		PlotAssociationLocusFrequencyVsAssociationThreshold.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		self.min_score_set = set()
	
	def preFileFunction(self, **keywords):
		"""
		2012.11.23
		"""
		pdata = PlotAssociationLocusFrequencyVsAssociationThreshold.preFileFunction(self, **keywords)
		pdata.chr2xy_ls = {}		
		pdata.widthList = []
		#reset these attributes
		pdata.min_score_set = set()
		return pdata
	
	def processRow(self, row=None, pdata=None):
		"""
		2012.8.2
			handles each row in each file
		"""
		#if not row['chromosome']:	#2012.12.3 ignore rows with empty chromosomes (stopgap bugfix, )
		#	return
		PlotGenomeWideData.processRow(self, row=row, pdata=pdata)
		
		reader = getattr(pdata, 'reader', None)
		tableObject =reader.getTableObject()
		
		"""
		xValue = getattr(row, self.xColumnHeader, None)	#start position
		stopPosition = row.stop
		width = abs(stopPosition - xValue)
		#add cumu start to xValue so that it'll be part of the whole genome
		cumuStart = self.chr_id2cumu_start.get(row.chromosome)
		xValue += cumuStart
		xValue = self.processValue(value=xValue, processType=self.logX)
		
		yValue = getattr(row, self.whichColumnHeader, None)	#no_of_results
		yValue = self.processValue(value=yValue, processType=self.logY)
		
		pdata.x_ls.append(xValue)
		pdata.y_ls.append(yValue)
		pdata.widthList.append(width)
		"""
		
		#add call/phenotype/analysis/min_overlap_ratio/total_no_of_results into self set.
		tableObject.addObjectListAttributeToSet(attributeName='call_method_id_ls', \
									setVariable=self.call_method_id_set)
		tableObject.addObjectListAttributeToSet(attributeName='phenotype_method_id_ls', \
									setVariable=self.phenotype_method_id_set)
		tableObject.addObjectListAttributeToSet(attributeName='analysis_method_id_ls', \
									setVariable=self.analysis_method_id_set)
		tableObject.addObjectAttributeToSet(attributeName='min_overlap_ratio', \
									setVariable=self.min_overlap_ratio_set)
		tableObject.addObjectAttributeToSet(attributeName='min_score', \
									setVariable=self.min_score_set)
		tableObject.addObjectAttributeToSet(attributeName='total_no_of_results', \
									setVariable=self.total_no_of_results_set)
		
		tableObject.addObjectListAttributeToSet(attributeName='call_method_id_ls', \
									setVariable=pdata.call_method_id_set)
		tableObject.addObjectListAttributeToSet(attributeName='phenotype_method_id_ls', \
									setVariable=pdata.phenotype_method_id_set)
		tableObject.addObjectListAttributeToSet(attributeName='analysis_method_id_ls', \
									setVariable=pdata.analysis_method_id_set)
		tableObject.addObjectAttributeToSet(attributeName='min_overlap_ratio', \
									setVariable=pdata.min_overlap_ratio_set)
		tableObject.addObjectAttributeToSet(attributeName='min_score', \
									setVariable=pdata.min_score_set)
		tableObject.addObjectAttributeToSet(attributeName='total_no_of_results', \
									setVariable=pdata.total_no_of_results_set)
		return 1
	
	def handleTitle(self,):
		"""
		2012.8.16
			add min_overlap_ratio, no of phenotypes, total_no_of_results into the title.
		"""
		if self.title:
			title = self.title
		else:
			phenotype_method_id_ls = list(self.phenotype_method_id_set)
			no_of_phenotypes = len(phenotype_method_id_ls)
			title = 'call %s, analysis %s, min_overlap %s, min_score %s, #results %s, #phenotypes %s'%\
					(utils.getStrOutOfList(list(self.call_method_id_set)),\
					utils.getStrOutOfList(list(self.analysis_method_id_set)),\
					utils.getStrOutOfList(list(self.min_overlap_ratio_set)),\
					utils.getStrOutOfList(list(self.min_score_set)),\
					utils.getStrOutOfList(list(self.total_no_of_results_set)), no_of_phenotypes)
			#title = yh_matplotlib.constructTitleFromTwoDataSummaryStat(self.invariantPData.x_ls, self.invariantPData.y_ls)
		pylab.title(title)
		return title
	
	def plot(self, x_ls=None, y_ls=None, pdata=None):
		"""
		2011-9-30
			get called by the end of fileWalker() for each inputFname.
		"""
		PlotGenomeWideData.plot(self, x_ls=x_ls, y_ls=y_ls, pdata=pdata)
		"""
		#handle the color properly
		association_group_key = (list(pdata.call_method_id_set)[0], list(pdata.analysis_method_id_set)[0])
		propertyData = self.association_group_key2propertyData.get(association_group_key, \
									self.association_group_key2propertyData.get(0))
		color = propertyData.color
		label = propertyData.label
		#plotObject = pylab.plot(x_ls, y_ls, self.formatString, c=color, alpha=0.8)[0]
		ind = numpy.array(x_ls)
		plotObject = pylab.plot(ind, y_ls, '.', color=color, alpha=0.6)[0]	#width=pdata.widthList, 
		#plotObject = pylab.bar(ind, y_ls, bottom=0, color=color, linewidth=0)	#width=pdata.widthList, 
		#linewidth=0 makes the bars without edges.
		
		self.setGlobalMinVariable(extremeVariableName='xMin', givenExtremeValue=min(x_ls))
		self.setGlobalMaxVariable(extremeVariableName='xMax', givenExtremeValue=max(x_ls))
		self.setGlobalMinVariable(extremeVariableName='yMin', givenExtremeValue=min(y_ls))
		self.setGlobalMaxVariable(extremeVariableName='yMax', givenExtremeValue=max(y_ls))
		
		#add the label, plot object of each curve into the global legend
		self.plotObjectLs.append(plotObject)
		#add call_method & analysis method into legend
		self.plotObjectLegendLs.append(label)
		"""
	
	def setup(self, **keywords):
		"""
		2012.10.15
			run before anything is run
		"""
		PlotGenomeWideData.setup(self, **keywords)
		
		PlotAssociationLocusFrequencyVsAssociationThreshold.setup(self, **keywords)
	

	def handleYLim(self,):
		"""
		2012.12.4 set the Y-axis in a more readable format.
		"""
		#self._handleYLim()	#set ylim according to ymin and ymax
		originalYMin, originalYMax = pylab.ylim()
		pylab.ylim(ymin=0.9, ymax=originalYMax)	#set ymin to 0.9 (= 9x10^-1).	setting it to 0 will eliminate the Y-axis frame (0 is no good for logarithm scale). 
		#pass
		#pylab.gca().margins(x=0.1, y=0.1)	#this is within the plotting area, set the padding on the x and y-direction.
		#no the padding outside the plot, which is set through "pylab.axes([-0.1, -0.1, 1.2, 1.2], frameon=False)"
		
if __name__ == '__main__':
	main_class = PlotAssociationLocusFrequencyOnGenome
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()