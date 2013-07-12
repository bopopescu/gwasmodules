#!/usr/bin/env python
"""
Examples:
	%s -i folderReduceLiftOverVCF/CAEY.sorted.vcf.gz -o CAEY.sameSite.concordance.tsv
	
	%s 
	
	%s 
	
Description:
	2013.07.12 program that calculates genotype concordance of same-position SNPs (chromosome, position) inside a VCF.

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from Bio.Seq import Seq
from pymodule import ProcessOptions, MatrixFile, PassingData
from pymodule.yhio.VCFFile import VCFFile
from pymodule.pegasus.mapper.AbstractVCFMapper import AbstractVCFMapper
from pymodule import SNP
from pymodule.yhio.SNP import nt2number

parentClass = AbstractVCFMapper
class CalculateSameSiteConcordanceInVCF(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update({
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		parentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	
	def readInSNPID2GenotypeVectorLs(self, inputFname=None):
		"""
		2013.07.11
		"""
		sys.stderr.write("Reading in same-SNP genotype data from %s ..."%(inputFname))
		
		reader = VCFFile(inputFname=inputFname)
		counter = 0
		real_counter = 0
		snp_pos2genotypeVectorLs = {}
		for vcfRecord in reader:
			key = (vcfRecord.chromosome, vcfRecord.position)
			if key not in snp_pos2genotypeVectorLs:
				snp_pos2genotypeVectorLs[key] = []
			else:
				real_counter += 1
			snp_pos2genotypeVectorLs[key].append(vcfRecord.data_row[1:])	#[0] is reference
			
			counter += 1
		del reader
		sys.stderr.write("%s snp coordinates from %s vcf records. %s entries with same-positions.\n"%\
						(len(snp_pos2genotypeVectorLs), counter, real_counter))
		return snp_pos2genotypeVectorLs
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		
		outputDir = os.path.split(self.outputFname)[0]
		if outputDir and not os.path.isdir(outputDir):
			os.makedirs(outputDir)
		
		snp_pos2genotypeVectorLs =self.readInSNPID2GenotypeVectorLs(self.inputFname)
		
		
		
		writer = MatrixFile(self.outputFname, openMode='w', delimiter='\t')
		header = ['chromosome', 'position', 'noOfMatches', 'noOfTotal', 'concordance']
		writer.writeHeader(header)
		
		
		counter = 0
		real_counter = 0
		no_of_pairs = 0
		snp_pos_ls = snp_pos2genotypeVectorLs.keys()
		snp_pos_ls.sort()
		for i in xrange(len(snp_pos_ls)):
			counter += 1
			key = snp_pos_ls[i]
			chromosome, position = snp_pos_ls[i][:2]
			genotypeVectorLs = snp_pos2genotypeVectorLs.get(key)
			if len(genotypeVectorLs)>1:
				real_counter += 1
				for k in xrange(0, len(genotypeVectorLs)-1):
					for l in xrange(k+1, len(genotypeVectorLs)):
						no_of_pairs +=1
						noOfMatches = 0
						noOfTotal = 0
						genotypeVector0 = genotypeVectorLs[k]
						genotypeVector1 = genotypeVectorLs[l]
						for j in xrange(len(genotypeVector0)):
							call1 = genotypeVector0[j]['GT']
							call2 = genotypeVector1[j]['GT']
							if call1!='NA' and call2!='NA':
								noOfTotal += 1
								if SNP.nt2number[call1]==SNP.nt2number[call2]:
									noOfMatches += 1
						if noOfTotal>0:
							concordance = float(noOfMatches)/float(noOfTotal)
						else:
							concordance = -1
						data_row = [chromosome, position,noOfMatches, noOfTotal, concordance ]
						writer.writerow(data_row)
		writer.close()
		sys.stderr.write("%s (out of %s, %s) snps have >1 same-position entries. %s pairs.\n"%(real_counter, counter, \
												real_counter/float(counter), no_of_pairs))

if __name__ == '__main__':
	main_class = CalculateSameSiteConcordanceInVCF
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()