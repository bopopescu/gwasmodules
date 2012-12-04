#!/usr/bin/env python
"""
Examples:
	%s
	
	%s  -s 1  -o locus_frequency_vs_min_score.png /tmp/count_association_locus.h5
		 --logY --positiveLog
	

Description:
	2012/11/22
		This program draws multiple curves of no_of_association_loci vs min_score (association threshold).
		Each curve is for one association group (call-method, analysis-method).
		Input is output of CountAssociationLocus.py.
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
from pymodule import yh_matplotlib, utils
from pymodule.plot.AbstractPlot import AbstractPlot
#from variation.src.plot.PlotAssociationLocus import PlotAssociationLocus

class PlotAssociationLocusFrequencyVsAssociationThreshold(AbstractPlot):
	__doc__ = __doc__
#						
	option_default_dict = AbstractPlot.option_default_dict.copy()
	# change
	option_default_dict[('whichColumnHeader', 0, )][0] = 'no_of_association_loci'
	option_default_dict[('whichColumnPlotLabel', 0, )][0] = 'numberOfPeaks'
	option_default_dict[('xColumnHeader', 1, )][0] = 'min_score'
	option_default_dict[('xColumnPlotLabel', 0, )][0] = 'minimum -log(Pvalue)'
	# change default of file-format
	option_default_dict[('inputFileFormat', 0, int)][0] = 2
	option_default_dict[('formatString', 1, )][0] = '.-'
	option_default_dict[('minNoOfTotal', 1, int)][0] = 1
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractPlot.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		
		self.call_method_id_set = set()
		self.phenotype_method_id_set = set()
		self.analysis_method_id_set = set()
		self.min_overlap_ratio_set = set()
		self.total_no_of_results_set = set()
	
	#this is from plot/PlotAssociationLocus.py
	analysis_method_id_locus_type_id2color = {(1,1):'k', (32,1):'r', (1,2):'b',(32,2):'g', 
							0:'#87CEFA'}	#lightskyblue: #87CEFA	R=135 G=206	B=250 ACCESS=16436871
	#association_group_key is (call_method_id, analysis_method_id)
	association_group_key2propertyData = {(75,1):PassingData(color='k', label='SNP-KW'), \
						(75,32): PassingData(color='r',label='SNP-EMMA'), \
						(57,1): PassingData(color='b', label='del-KW'),\
						(57,32): PassingData(color='g', label='del-EMMA'),\
						0:PassingData(color='#87CEFA', label='Unknown')}
	#lightskyblue: #87CEFA	R=135 G=206	B=250 ACCESS=16436871
	
	def preFileFunction(self, **keywords):
		"""
		2012.11.23
		"""
		pdata = AbstractPlot.preFileFunction(self, **keywords)
		#reset these attributes
		pdata.call_method_id_set = set()
		pdata.phenotype_method_id_set = set()
		pdata.analysis_method_id_set = set()
		pdata.min_overlap_ratio_set = set()
		pdata.total_no_of_results_set = set()
		return pdata
	
	def processRow(self, row=None, pdata=None):
		"""
		2012.8.2
			handles each row in each file
		"""
		col_name2index = getattr(pdata, 'col_name2index', None)
		x_ls = getattr(pdata, 'x_ls', None)
		y_ls = getattr(pdata, 'y_ls', None)
		
		reader = getattr(pdata, 'reader', None)
		groupObject =reader.getGroupObject()
		
		xValue = getattr(row, self.xColumnHeader, None)
		yValue = getattr(row, self.whichColumnHeader, None)
		xValue = self.processValue(value=xValue, processType=self.logX)
		yValue = self.processValue(value=yValue, processType=self.logY)
		
		x_ls.append(xValue)
		y_ls.append(yValue)
		
		#add call/phenotype/analysis/min_overlap_ratio/total_no_of_results into self set.
		utils.addObjectListAttributeToSet(objectVariable=row, attributeName='call_method_id_ls', \
									setVariable=self.call_method_id_set, data_type=int)
		utils.addObjectListAttributeToSet(objectVariable=row, attributeName='phenotype_method_id_ls', \
									setVariable=self.phenotype_method_id_set, data_type=int)
		utils.addObjectListAttributeToSet(objectVariable=row, attributeName='analysis_method_id_ls', \
									setVariable=self.analysis_method_id_set, data_type=int)
		utils.addObjectAttributeToSet(objectVariable=row, attributeName='min_overlap_ratio', \
									setVariable=self.min_overlap_ratio_set)
		utils.addObjectAttributeToSet(objectVariable=row, attributeName='total_no_of_results', \
									setVariable=self.total_no_of_results_set)
		
		utils.addObjectListAttributeToSet(objectVariable=row, attributeName='call_method_id_ls', \
									setVariable=pdata.call_method_id_set, data_type=int)
		utils.addObjectListAttributeToSet(objectVariable=row, attributeName='phenotype_method_id_ls', \
									setVariable=pdata.phenotype_method_id_set, data_type=int)
		utils.addObjectListAttributeToSet(objectVariable=row, attributeName='analysis_method_id_ls', \
									setVariable=pdata.analysis_method_id_set, data_type=int)
		utils.addObjectAttributeToSet(objectVariable=row, attributeName='min_overlap_ratio', \
									setVariable=pdata.min_overlap_ratio_set)
		utils.addObjectAttributeToSet(objectVariable=row, attributeName='total_no_of_results', \
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
			title = 'min_overlap %s, #results %s, #phenotypes %s'%(utils.getStrOutOfList(list(self.min_overlap_ratio_set)),\
						utils.getStrOutOfList(list(self.total_no_of_results_set)), no_of_phenotypes)
			#title = yh_matplotlib.constructTitleFromTwoDataSummaryStat(self.invariantPData.x_ls, self.invariantPData.y_ls)
		pylab.title(title)
		return title
	
	def plot(self, x_ls=None, y_ls=None, pdata=None):
		"""
		2011-9-30
			get called by the end of fileWalker() for each inputFname.
		"""
		#handle the color properly
		association_group_key = (list(pdata.call_method_id_set)[0], list(pdata.analysis_method_id_set)[0])
		propertyData = self.association_group_key2propertyData.get(association_group_key, \
									self.association_group_key2propertyData.get(0))
		color = propertyData.color
		label = propertyData.label
		plotObject = pylab.plot(x_ls, y_ls, self.formatString, c=color, alpha=0.8)[0]
		self.setGlobalMinVariable(extremeVariableName='xMin', givenExtremeValue=min(x_ls))
		self.setGlobalMaxVariable(extremeVariableName='xMax', givenExtremeValue=max(x_ls))
		self.setGlobalMinVariable(extremeVariableName='yMin', givenExtremeValue=min(y_ls))
		self.setGlobalMaxVariable(extremeVariableName='yMax', givenExtremeValue=max(y_ls))
		
		#add the label, plot object of each curve into the global legend
		#add call_method & analysis method into legend
		self.addPlotLegend(plotObject=plotObject, legend=label, pdata=pdata)
		
		#AbstractPlot.plot(self, x_ls=x_ls, y_ls=y_ls, pdata=pdata)
	

if __name__ == '__main__':
	main_class = PlotAssociationLocusFrequencyVsAssociationThreshold
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()