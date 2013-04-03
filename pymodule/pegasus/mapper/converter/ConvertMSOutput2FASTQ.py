#!/usr/bin/env python
"""
Examples:
	# 2013.2.11 input is the msHOT-lite Heng-Li custom output
	%s   -i 1534_788_2009098_GA_vs_524.msHOT_lite_output.txt -o 1534_788_2009098_GA_vs_524.msHOT_lite.fq.gz  --inputFileFormat 2
	
	# 2013.2.11 input is the ms output
	%s -i 1534_788_2009098_GA_vs_524.msHOT_lite_traditional_output.txt.gz -o 1534_788_2009098_GA_vs_524.msHOT_lite_traditional_output.fa.gz
	
Description:
	2013.2.10 input is ms/msHOT/msHOT-lite simulation output.
			output is a fastq file. (then fq2psmcfa from PSMC package, Li, Durbin 2011 could use) 

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
import h5py
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule import SNP
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper


class ConvertMSOutput2FASTQ(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.update({
							('defaultBase', 1, ): ["A", '', 1, "corresponding to allele 0",],\
							('alternativeBase', 1, ): ["G", '', 1, "corresponding to allele 1",],\
							
							('defaultBaseQuality', 1, ): ["z", '', 1, "base quality to be assigned to each base.\n\
		Assuming Sanger format. z is ascii no. 122, corresponding to quality 87 ",],\
							('ploidy', 1, int): [2, '', 1, "1: haploid, one sample, one individual; \
					2: diploid, take every two consecutive samples as one individual. Other ploids are not supported yet.", ],\
							('inputFileFormat', 1, int): [1, '', 1, "1: input is ms/msHOT output (no good, non-polymorphic sites ignored); 2: in Heng Li output", ],\
							})
	def __init__(self,  **keywords):
		"""
		"""
		AbstractMapper.__init__(self, **keywords)
		
		self.convertFuncDict = {1: self.convertMSOutput, 2: self.convertMSHOTLiteOutput}
	
	def convertMSOutput(self, inf=None, outf=None):
		"""
		2013.2.11
			not right, as it ignored all the non-polymorphic loci
		"""
		isSampleBegun =False	#True when the sample data is up in the next line
		chromosomeNumber = 1
		for line in inf:
			content = line.strip()
			if content=="//":
				segsitesLineContent = inf.next().strip()
				no_of_segsites = int(segsitesLineContent.split()[-1])
				inf.next()	#an empty line
				positionList = inf.next().strip().split()[1:]
				isSampleBegun =True
			elif isSampleBegun:
				individualAlleleList = []
				if self.ploidy==2:
					nextSampleContent = inf.next().strip()
					for i in xrange(len(content)):
						individualAlleleList.append('%s%s'%(content[i], nextSampleContent[i]))
				else:	#haploid
					for i in xrange(len(content)):
						individualAlleleList.append(content[i])
				outputBaseList = []
				for i in xrange(len(individualAlleleList)):
					if individualAlleleList[i]=='00':
						outputBaseList.append(self.defaultBase)
					elif individualAlleleList[i]=='11':
						outputBaseList.append(self.alternativeBase)
					else:	#heterozygous
						het_in_nt = '%s%s'%(self.defaultBase, self.alternativeBase)
						het_in_number = SNP.nt2number.get(het_in_nt)
						het_in_single_char_nt = SNP.number2single_char_nt.get(het_in_number)
						outputBaseList.append(het_in_single_char_nt)
				
				outputLine = ''.join(outputBaseList)
				outf.write('@Chr%s\n'%(chromosomeNumber))
				outf.write("%s\n"%outputLine)
				
				outf.write("+\n")
				outf.write("%s\n"%(self.defaultBaseQuality*len(outputBaseList)))
				chromosomeNumber += 1
		
	def convertMSHOTLiteOutput(self, inf=None, outf=None):
		"""
		2013.2.11
		./msHOT-lite 2 1 -t 84989.8346003745 -r 34490.1412746802 30000000 -l -en 0.0013 1 0.0670 -en 0.0022 1 0.3866 -en 0.0032 1 0.3446 -en 0.0044 1 0.21
				79 -en 0.0059 1 0.1513 -en 0.0076 1 0.1144 -en 0.0096 1 0.0910 -en 0.0121 1 0.0757 -en 0.0150 1 0.0662 -en 0.0184 1 0.0609 -en 0.0226 1 0.0583 -en
				 0.0275 1 0.0572 -en 0.0333 1 0.0571 -en 0.0402 1 0.0577 -en 0.0485 1 0.0589 -en 0.0583 1 0.0603 -en 0.0700 1 0.0615 -en 0.0839 1 0.0624 -en 0.100
				5 1 0.0632 -en 0.1202 1 0.0641 -en 0.1437 1 0.0651 -en 0.1716 1 0.0663 -en 0.2048 1 0.0678 -en 0.2444 1 0.0696 -en 0.2914 1 0.0719 -en 0.3475 1 0.
				0752 -en 0.4935 1 0.0794 
				//
				@begin 6422
				30000000
				1100    01
				6074    10
				
				29966899        10
				29971027        01
				29973740        01
				29982767        01
				29985696        10
				@end
		"""
		isSampleBegun =False	#True when the sample data is up in the next line
		chromosomeNumber = 1
		previousPolymorhicSitePosition = 0
		individualByGenotypeMatrix = []
		no_of_segsites_encountered = 0
		for line in inf:
			content = line.strip()
			if content=="//":
				segsitesLineContent = inf.next().strip()
				no_of_segsites = int(segsitesLineContent.split()[-1])
				no_of_sites = int(inf.next().strip())
				isSampleBegun =True
			elif content =="@end":
				break
			elif isSampleBegun:
				polymorphicPosition, haplotypeAlleleList = content.split()
				polymorphicPosition = int(polymorphicPosition)
				for i in xrange(0, len(haplotypeAlleleList), self.ploidy):	#go through each indivdiual at a time
					if no_of_segsites_encountered==0:	#first time adding genotype
						individualByGenotypeMatrix.append([])
					individualGenotype = haplotypeAlleleList[i*self.ploidy:(i+1)*self.ploidy]
					#turn it into ACGT
					if individualGenotype=='0'*self.ploidy:
						individualGenotype_in_nt = self.defaultBase
					elif individualGenotype=='1'*self.ploidy:
						individualGenotype_in_nt = self.alternativeBase
					else:	#heterozygous
						het_in_nt = '%s%s'%(self.defaultBase, self.alternativeBase)
						het_in_number = SNP.nt2number.get(het_in_nt)
						individualGenotype_in_nt = SNP.number2single_char_nt.get(het_in_number)
					
					for j in xrange(previousPolymorhicSitePosition+1, polymorphicPosition):
						individualByGenotypeMatrix[i].append(self.defaultBase)
					individualByGenotypeMatrix[i].append(individualGenotype_in_nt)
				no_of_segsites_encountered  += 1
				previousPolymorhicSitePosition = polymorphicPosition
				
		
		for row in individualByGenotypeMatrix:
				outputLine = ''.join(row)
				outf.write('@Chr%s\n'%(chromosomeNumber))
				outf.write("%s\n"%outputLine)
				
				outf.write("+\n")
				outf.write("%s\n"%(self.defaultBaseQuality*len(row)))
				chromosomeNumber += 1
		
	def run(self):
		"""
		2013.2.11
			input looks like (inputFileFormat=1)
				msHOT-lite 2 1 -t 4781.50413187402 -r 790.4466018 ...
				//
				segsites: 40567
				
				positions: 0.0002 0.0003
				001001101011011001...
				101001010100101111...
				...
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		inf = utils.openGzipFile(self.inputFname, 'r')
		
		outf = utils.openGzipFile(self.outputFname, openMode='w')
		self.convertFuncDict[self.inputFileFormat](inf=inf, outf=outf)
		
		inf.close()
		outf.close()
		

if __name__ == '__main__':
	main_class = ConvertMSOutput2FASTQ
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()