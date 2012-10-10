#!/usr/bin/env python
"""
Examples:
	#run on hoffman2's condor
	%s -d ~/NetworkData/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta
		-i ~/script/vervet/data/OphoffMethylation/HumanMethylation450_15017482_v.1.1.fasta.gz 
		-a 3 -l hcondor -j hcondor   -C 1 -o workflow/BlastHumanMethylation450_ProbeSeqAgainst524_3Mismatches.xml
		-s ~/bin/blast/bin/blastall -f ~/bin/blast/bin/formatdb
	
	%s -d /Network/Data/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta -i /tmp/input.fasta
		-a 2 -o /tmp/output.tsv

Description:
	2012.5.23
		a wrapper around blastall (blastn after 2012.8.19), parse and filter.
		accept multiple input files trailing all arguments.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from AbstractMapper import AbstractMapper
from Bio.Blast import NCBIXML, NCBIStandalone
from Bio.Blast.Applications import NcbiblastnCommandline
import cStringIO

class BlastWrapper(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	#option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
							('blastallPath', 1, ): ['/usr/bin/blastall', 'l', 1, 'path to blastall', ],\
							#('blastnPath', 1, ): ['%s/bin/ncbi-blast/bin/blastn', 'l', 1, 'path to blastn', ],\
							('databaseFname', 1, ): ['', 'd', 1, 'filename of the database to blast against, must be indexed', ],\
							('minNoOfIdentities', 0, int): [None, 'm', 1, 'minimum number of identities between a query and target', ],\
							('maxNoOfMismatches', 0, int): [None, 'a', 1, 'minimum number of mismatches between a query and target', ],\
							('minIdentityPercentage', 0, float): [None, 'n', 1, 'minimum percentage of identities between a query and target', ],\
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		2012.5.23
		"""
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		#self.blastnPath =  self.insertHomePath(self.blastnPath, self.home_path)
		self.blastallPath =  self.insertHomePath(self.blastallPath, self.home_path)
	
	def runBlast(self, inputFname=None, databaseFname=None, outputFname=None, outputFnamePrefix=None, \
				blastallPath=None, minNoOfIdentities=None, \
				maxNoOfMismatches=None,\
				minIdentityPercentage=None, maxNoOfHits=10):
		"""
		2012.8.19
			output xml dump if outputFnamePrefix is given.
		2012.5.23
		  -p  Program Name [String]
		  -d  Database [String]
		    default = nr
		  -i  Query File [File In]
		    default = stdin
		  -e  Expectation value (E) [Real]
		    default = 10.0

			blastall align_view option values:
				0 = pairwise,
				1 = query-anchored showing identities,
				2 = query-anchored no identities,
				3 = flat query-anchored, show identities,
				4 = flat query-anchored, no identities,
				5 = query-anchored no identities and blunt ends,
				6 = flat query-anchored, no identities and blunt ends,
				7 = XML Blast output,
				8 = tabular, 
				9 tabular with commresult_handleent lines
				10 ASN, text
				11 ASN, binary [Integer]
				    default = 0
				    range from 0 to 11
		
		"""
		
		result_handle, error_info = NCBIStandalone.blastall(blastallPath, "blastn", databaseFname, inputFname, align_view=7)
		
		
		#blastn_cline = NcbiblastnCommandline(cmd=self.blastnPath, query=inputFname, db=databaseFname, evalue=0.001,\
		#									outfmt=5, out="opuntia.xml")	#outfmt 5 is xml output.
		
		#error_info = error_info.read()	#2010-4-14 this read() causes program to hang out forever. ???
		#if error_info:
		#	sys.stderr.write("%s"%error_info)
		if outputFnamePrefix:
			outf = open('%s.xml'%(outputFnamePrefix), 'w')
			blastContent = result_handle.read()
			outf.write(blastContent)
			outf.close()
			result_handle = cStringIO.StringIO(blastContent)
		blast_records = NCBIXML.parse(result_handle)
		
		if self.report:
			sys.stderr.write("finished blasting.\n")
		
		counter = 0
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		header = ['queryID', "queryStart", "queryEnd", 'queryLength', 'targetChr', 'targetStart', 'targetStop', \
				'targetLength', 'noOfIdentities', \
				'noOfMismatches', 'identityPercentage']
		writer.writerow(header)
		for blast_record in blast_records:
			no_of_hits = min(maxNoOfHits, len(blast_record.alignments))	# top 50 or the number of available alignments
			# each alignment is one chromosome (=one fasta record).
			for i in range(no_of_hits):
				alignment_title = blast_record.alignments[i].title
				targetChr = blast_record.alignments[i].hit_def
				targetLength = blast_record.alignments[i].length
				for hsp in blast_record.alignments[i].hsps:
					hitIsGood = True
					noOfMismatches = blast_record.query_length - hsp.identities
					identityPercentage = float(hsp.identities)/float(blast_record.query_length)
					if minNoOfIdentities is not None and hsp.identities < minNoOfIdentities:
						hitIsGood = False
					if maxNoOfMismatches is not None and noOfMismatches>maxNoOfMismatches:
						hitIsGood = False
					if minIdentityPercentage is not None and identityPercentage<minIdentityPercentage:
						hitIsGood = False
					if hitIsGood:
						counter += 1
						result_entry = [blast_record.query, hsp.query_start, hsp.query_end, blast_record.query_length,\
									targetChr, hsp.sbjct_start, hsp.sbjct_end, targetLength, hsp.identities, noOfMismatches,\
									identityPercentage]
						#20104-25 hsp.strand is always (None, None), hsp.frame is either (1,1) or (1, -1) when the query's end < start
							#[query name (probe id and pos) , alignment title , number of matches, pos in contig ]
						writer.writerow(result_entry)
		if self.report:
			sys.stderr.write("%s blast records, %s pass the filter.\n"%\
							(len(blast_records), counter))
		del writer
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		self.runBlast(inputFname=self.inputFname, databaseFname=self.databaseFname, outputFname=self.outputFname, \
					outputFnamePrefix = self.outputFnamePrefix,\
					blastallPath=self.blastallPath, minNoOfIdentities=self.minNoOfIdentities, \
					maxNoOfMismatches=self.maxNoOfMismatches,\
					minIdentityPercentage=self.minIdentityPercentage)
		

if __name__ == '__main__':
	main_class = BlastWrapper
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()