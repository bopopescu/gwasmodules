#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i /tmp/VRCSamples.tsv --whichColumnHeader=sampleID
		-s 1.0  --sampleSize 5
		--plinkIBDCheckOutputFname PlinkIBDCheck/PlinkIBDCheck_Method38_W50Z20R0.9.2012.9.13T102614/ibdCheckIBDCheck/LDPrunedMerged_ibdCheck.tsv
		-o /tmp/VRCSamples_sampled.tsv

Description:
	2012.10.15
		sample rows from a matrix-like file
		
	If "-i ..." is given, it is regarded as one of the input files (plus the ones in trailing arguments). 
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

import csv
import random
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter
from pymodule import yh_matplotlib, SNP
from pymodule.yhio.AbstractMatrixFileWalker import AbstractMatrixFileWalker
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper

class SampleRows(AbstractMatrixFileWalker):
	__doc__ = __doc__
	option_default_dict = AbstractMatrixFileWalker.option_default_dict.copy()
	option_default_dict.update({
						('sampleSize', 1, int): [None, '', 1, 'number of samples (rows) to be sampled from input, \n\
	uniformly if plinkIBDCheckOutputFname is not given.'],\
						('plinkIBDCheckOutputFname', 0, ): [None, '', 1, 'file that contains IBD check result, PI_HAT=relateness.\n\
	at least 3-columns with header: IID1, IID2, PI_HAT. IID1 and IID2 should match the whichColumn (whichColumnHeader) of inputFname.\n\
	The sampling will try to avoid sampling close pairs, PI_HAT(i,j)<=maxIBDSharing'],\
						('maxIBDSharing', 1, float): [0.1, '', 1, 'This argument caps the maximum IBD sharing among any pair within the sampled.'],\
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractMatrixFileWalker.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
		
	def setup(self, **keywords):
		"""
		2012.10.15
			run before anything is run
		"""
		AbstractMatrixFileWalker.setup(self, **keywords)
		#read in the IBD check result
		if self.plinkIBDCheckOutputFname:
			ibdData = SNP.readAdjacencyListDataIntoMatrix(inputFname=self.plinkIBDCheckOutputFname, rowIDHeader="IID1", colIDHeader="IID2", \
								rowIDIndex=None, colIDIndex=None, \
								dataHeader="PI_HAT", dataIndex=None, hasHeader=True)
		else:
			ibdData = None
		self.ibdData = ibdData
		self.data_matrix = []	#data structure to store all rows during fileWalker()
	
	def processRow(self, row=None, pdata=None):
		"""
		2012.10.7
		"""
		returnValue = 1
		self.data_matrix.append(row)
		self.col_name2index = getattr(pdata, 'col_name2index', None)
		
		return returnValue
	
	def reduce(self, **keywords):
		"""
		2012.10.15
			run after all files have been walked through
		"""
		#sample the data
		noOfTotalRows = len(self.data_matrix)
		real_counter = 0
		if self.sampleSize<noOfTotalRows:
			if self.ibdData:
				col_name2index = getattr(self, 'col_name2index', None)
				sampledRowIndexList = []
				if col_name2index:
					sampledRowIndexSet = set()
					while len(sampledRowIndexSet)<self.sampleSize:
						randomIndex = random.randint(0,noOfTotalRows-1)
						if randomIndex in sampledRowIndexSet:	#must not be in there already
							continue
						row = self.data_matrix[randomIndex]
						if self.whichColumnHeader:
							whichColumn = col_name2index.get(self.whichColumnHeader, None)
						elif self.whichColumn:
							whichColumn = self.whichColumn
						else:
							whichColumn = None
						if whichColumn is not None:
							rowID = row[whichColumn]
							includeInTheSampling = True
							for sampledRowIndex in sampledRowIndexSet:	#not too close to anyone previously sampled
								previousSampledRowID = self.data_matrix[sampledRowIndex][whichColumn]
								#getting the relateness
								relateness = self.ibdData.getCellDataGivenRowColID(rowID, previousSampledRowID)
								if relateness>=self.maxIBDSharing:
									includeInTheSampling = False
							if includeInTheSampling:
								sampledRowIndexSet.add(randomIndex)
					sampledRowIndexList = list(sampledRowIndexSet)
			else:
				sampledRowIndexList = random.sample(xrange(noOfTotalRows), self.sampleSize)
		else:	#take all
			sampledRowIndexList = xrange(noOfTotalRows)
		for i in sampledRowIndexList:
			row = self.data_matrix[i]
			self.invariantPData.writer.writerow(row)
			real_counter +=1 
		
		fraction = float(real_counter)/float(noOfTotalRows)
		sys.stderr.write("%s/%s (%.3f) selected.\n"%(real_counter, noOfTotalRows, fraction))
		
		#close the self.invariantPData.writer
		AbstractMatrixFileWalker.reduce(self, **keywords)

if __name__ == '__main__':
	main_class = SampleRows
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
