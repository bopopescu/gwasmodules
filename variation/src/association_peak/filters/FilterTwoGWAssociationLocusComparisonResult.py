#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i association_min_score_3.0/compare_association_call_75_to_call_57_ana_32_loci_min_score_3.0_min_overlap_0.8.h5
		--associationLocusCategory 1  --phenotypeIDList 1-500
		-o association_min_score_3.0/compare_association_call_75_to_call_57_ana_32_loci_min_score_3.0_min_overlap_0.8_associationLocusCategory1.h5

Description:
	2013.2.24 This program filters output of CompareTwoGWAssociationLocusByPhenotypeVector.py via phenotypeIDList or associationLocusCategory.
	
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import copy
import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, getColName2IndexFromHeader, figureOutDelimiter,\
	YHFile, AbstractMatrixFileWalker, castPyTablesRowIntoPassingData
from pymodule import RBDict, CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
from variation.src.association_peak.CompareTwoGWAssociationLocusByPhenotypeVector import TwoGenomeWideAssociationLocusMapTable


class FilterTwoGWAssociationLocusComparisonResult(AbstractMatrixFileWalker):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractMatrixFileWalker.option_default_dict)
	# change default of file-format
	option_default_dict[('inputFileFormat', 0, int)][0] = 2
	option_default_dict[('outputFileFormat', 0, int)][0] = 2
	option_default_dict.update({
					('phenotypeIDList', 0, ): [None, '', 1, 'phenotypes not from this list will be tossed from total_phenotype_ls_in_str, overlap_phenotype_ls_in_str, input1_only_phenotype_ls_in_str, input2_only_phenotype_ls_in_str' ],\
					('associationLocusCategory', 0, int): [0, '', 1, '0: all; 1: confirmatory; 2: del-gwas-novel; 3:SNP-gwas-novel; 4:specific-in-both; 5: others' ],\
				})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractMatrixFileWalker.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		if hasattr(self, 'phenotypeIDList'):
			self.phenotypeIDList = utils.getListOutOfStr(self.phenotypeIDList)
		else:
			self.phenotypeIDList = None
	
	def setup(self, **keywords):
		"""
		2012.11.25 setup the output
		"""
		writer = None
		if self.outputFileFormat==1:
			suffix = os.path.splitext(self.outputFname)[1]
			if self.outputFname and suffix!='.png':
				writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		else:	#HDF5MatrixFile
			writer = YHFile(self.outputFname, openMode='w', rowDefinition=TwoGenomeWideAssociationLocusMapTable)
		self.writer = writer
	
	def calculateNormalizedElipseDistanceToCenter(self, tableObject=None, centerAndRotationScalarAndRadius=None):
		"""
		2013.2.10
			centerAndRotationScalarAndRadius defines the eclipse to which one category of peaks belong.
			It's in the format of (center_x, center_y, rotationScalar, majorRadius, minorRadius).
				i.e. 0,0,0, 0.3, 0.3
		"""
		center_x, center_y = centerAndRotationScalarAndRadius[:2]
		x = tableObject.fraction_of_input1_only_phenotypes
		y = tableObject.fraction_of_input2_only_phenotypes
		
		rotationScalar = centerAndRotationScalarAndRadius[2]
		majorRadius, minorRadius = centerAndRotationScalarAndRadius[3:5]
		if rotationScalar!=0:	#there's a rotation in terms of the elipse
			x_prime = x*math.cos(rotationScalar) + y*math.sin(rotationScalar)
			y_prime = -x*math.sin(rotationScalar) + y*math.cos(rotationScalar)
			x = x_prime
			y = y_prime
		
		normalizedElipseDistanceToCenter = math.pow((x-center_x)/majorRadius,2) + math.pow((y-center_y)/minorRadius, 2)
		return normalizedElipseDistanceToCenter
		
	
	def filterAssociationLocusByPhenotypeOverlapFraction(self, tableObject=None, confirmatoryLocusCenter=(0,0, 0, 0.3, 0.3), \
										novelDelGWASLocusCenter=(0, 1.0, 0, 0.2, 0.7),\
										novelSNPGWASLocusCenter=(1.0, 0, 0, 0.7, 0.2), \
										bothSpecificGWASLocusCenter=(0.5,0.5, -math.pi/4.0, 0.4, 0.1)):
		"""
		2013.2.19
			every peak category occupies an eclipse on a 2D (fraction_of_input1_only_phenotypes X fraction_of_input2_only_phenotypes) plane define by 
				centerAndRotationScalarAndRadius, which is in the format of
				(center_x, center_y, rotationScalar, majorRadius, minorRadius).
			The eclipse for bothSpecificGWASLocus is a clock-wise 45-degree rotated eclipse.
		"""
		categoryID2centerAndRotationScalarAndRadius = {1: confirmatoryLocusCenter,
							2: novelDelGWASLocusCenter,
							3: novelSNPGWASLocusCenter,
							4: bothSpecificGWASLocusCenter}
		categoryID2AssociationLocusList = {}
		for categoryID, centerAndRotationScalarAndRadius in categoryID2centerAndRotationScalarAndRadius.iteritems():
			normalizedElipseDistanceToCenter = self.calculateNormalizedElipseDistanceToCenter(tableObject, centerAndRotationScalarAndRadius)
			if normalizedElipseDistanceToCenter<=1.0:
				if categoryID not in  categoryID2AssociationLocusList:
					categoryID2AssociationLocusList[categoryID] = []
				categoryID2AssociationLocusList[categoryID].append(tableObject)
		if len(categoryID2AssociationLocusList)==0:
			#5 is for all other un-classified association loci 
			categoryID2AssociationLocusList[5] = [tableObject]
		return categoryID2AssociationLocusList
		
	
	def filterTableRowByPhenotypeVector(self, tableObject=None, phenotypeIDList=None):
		"""
		2013.2.19
			modify the tableObject and gets rid of phenotypes in (overlap_phenotype_ls_in_str,
					input1_only_phenotype_ls_in_str, input2_only_phenotype_ls_in_str)
				that are not in phenotypeIDList
		"""
		if phenotypeIDList:
			phenotypeIDSet = set(phenotypeIDList)
		else:	#return it without modifying
			return tableObject
		phenotypeRelatedAttributeNameList = [['total_phenotype_ls_in_str', 'no_of_total_phenotypes', 'fraction_of_total_phenotypes'],\
						['overlap_phenotype_ls_in_str', 'no_of_overlap_phenotypes', 'fraction_of_overlap_phenotypes'],\
						['input1_only_phenotype_ls_in_str', 'no_of_input1_only_phenotypes', 'fraction_of_input1_only_phenotypes'], \
						['input2_only_phenotype_ls_in_str', 'no_of_input2_only_phenotypes', 'fraction_of_input2_only_phenotypes']]
		no_of_total_phenotypes = None	#no of total significant phenotypes
		for phenotypeRelatedAttributeName in phenotypeRelatedAttributeNameList:
			phenotypeListAttributeName = phenotypeRelatedAttributeName[0]
			dependentAttributeNameList = phenotypeRelatedAttributeName[1:]
			phenotypeListAttributeValue = getattr(tableObject, phenotypeListAttributeName, None)
			givenPhenotypeIDList = getListOutOfStr(phenotypeListAttributeValue)
			modified_phenotypeIDList = []
			for phenotype_id in givenPhenotypeIDList:
				if phenotype_id in phenotypeIDSet:
					modified_phenotypeIDList.append(phenotype_id)
			#modify the list attribute
			modified_phenotypeIDListInStr = utils.getSuccinctStrOutOfList(modified_phenotypeIDList)
			setattr(tableObject, phenotypeListAttributeName, modified_phenotypeIDListInStr)
			
			#modify the #phenotypes & phenotypeFraction attributes
			noOfPhenotypeAttributeName, phenotypeFractionAttributeName = dependentAttributeNameList[:2]
			noOfPhenotypes = len(modified_phenotypeIDList)
			#calculate fraction of phenotypes
			if noOfPhenotypeAttributeName=='no_of_total_phenotypes':	# this fraction is #significant-phenotypes/#all-phenotypes-in-association
				no_of_total_phenotypes = noOfPhenotypes
				oldNoOfPhenotypes = getattr(tableObject, noOfPhenotypeAttributeName, None)
				oldPhenotypeFraction = getattr(tableObject, phenotypeFractionAttributeName, None)
				phenotypeFraction = noOfPhenotypes/float(oldNoOfPhenotypes/oldPhenotypeFraction)
			elif no_of_total_phenotypes is not None and no_of_total_phenotypes>0:
				phenotypeFraction = noOfPhenotypes/float(no_of_total_phenotypes)
			else:
				phenotypeFraction = None
			setattr(tableObject, noOfPhenotypeAttributeName,  noOfPhenotypes)
			setattr(tableObject, phenotypeFractionAttributeName,  phenotypeFraction)
		return tableObject
	
	def processRow(self, row=None, pdata=None):
		"""
		2012.10.15
			return returnValue: 0 if not included or nothing is done on it.
				1 if included or something is carried out on it.
		2012.8.2
			handles each row in each file, here it replaces the yValue
		"""
		returnValue = 0
		col_name2index = getattr(pdata, 'col_name2index', None)
		y_ls = getattr(pdata, 'y_ls', None)
		pdata = castPyTablesRowIntoPassingData(row)
		if self.associationLocusCategory>0:
			categoryID2AssociationLocusList = self.filterAssociationLocusByPhenotypeOverlapFraction(pdata)
			toOutputAssociationLocusList = categoryID2AssociationLocusList.get(self.associationLocusCategory)
		else:
			toOutputAssociationLocusList = [pdata]
		
		if toOutputAssociationLocusList:
			for associationLocus in toOutputAssociationLocusList:
				if self.phenotypeIDList:
					self.filterTableRowByPhenotypeVector(associationLocus, phenotypeIDList=self.phenotypeIDList)
				if self.writer:
					self.writer.writeOneCell(associationLocus, cellType=2)
					returnValue = 1
				
		return returnValue
	

if __name__ == '__main__':
	main_class = FilterTwoGWAssociationLocusComparisonResult
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
