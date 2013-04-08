#!/usr/bin/env python
"""
2012.1.17
	a common class for pegasus workflows that work on VCF variant files
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq, utils
from pymodule.yhio.MatrixFile import MatrixFile
from pymodule.yhio.VCFFile import VCFFile
from AbstractNGSWorkflow import AbstractNGSWorkflow

class AbstractVCFWorkflow(AbstractNGSWorkflow):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractNGSWorkflow.option_default_dict)
	option_default_dict.update({
						('inputDir', 0, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('minDepth', 0, float): [0, 'm', 1, 'minimum depth for a call to regarded as non-missing', ],\
						
						('intervalOverlapSize', 1, int): [300000, 'U', 1, 'overlap #bps/#loci between adjacent intervals from one contig/chromosome,\
				only used for TrioCaller, not for SAMtools/GATK', ],\
						('intervalSize', 1, int): [5000000, 'Z', 1, '#bps/#loci for adjacent intervals from one contig/chromosome (alignment or VCF)', ],\
						("ligateVcfPerlPath", 1, ): ["%s/bin/umake/scripts/ligateVcf.pl", '', 1, 'path to ligateVcf.pl'],\
						})
	def __init__(self,  **keywords):
		"""
		2012.1.17
		"""
		AbstractNGSWorkflow.__init__(self, **keywords)
	
	def extra__init__(self):
		"""
		2013.2.14
		"""
		AbstractNGSWorkflow.extra__init__(self)
		if getattr(self, "inputDir", None):
			self.inputDir = os.path.abspath(self.inputDir)
		
		if hasattr(self, 'ligateVcfPerlPath'):
			self.ligateVcfPerlPath =  self.insertHomePath(self.ligateVcfPerlPath, self.home_path)
		
	
	def addAddVCFFile2DBJob(self, executable=None, inputFile=None, genotypeMethodShortName=None,\
						logFile=None, format=None, data_dir=None, checkEmptyVCFByReading=None, commit=False, \
						parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.12.12 use extended argument name
		2012.10.6 use addGenericDBJob() instead of addGenericJob()
		2012.8.30 moved from vervet/src/AddVCFFolder2DBWorkflow.py
		2012.6.27
		"""
		extraArgumentList = ['--format', format]
		if logFile:
			extraArgumentList.extend(["--logFilename", logFile])
		if data_dir:
			extraArgumentList.extend(['--data_dir', data_dir])
		if checkEmptyVCFByReading:
			extraArgumentList.extend(['--checkEmptyVCFByReading'])
		if genotypeMethodShortName:
			extraArgumentList.extend(['--genotypeMethodShortName', genotypeMethodShortName, ])
		if commit:
			extraArgumentList.append('--commit')
		if extraArguments:
			extraArgumentList.append(extraArguments)
		
		job= self.addGenericDBJob(executable=executable, inputFile=inputFile, outputFile=None, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=[logFile],\
						transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, **keywords)
		return job
	
	def registerExecutables(self, workflow=None):
		"""
		"""
		AbstractNGSWorkflow.registerExecutables(self, workflow)
		
		if not workflow:
			workflow = self
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		clusters_size = self.clusters_size
		site_handler = self.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		
		PlotVCFtoolsStat = Executable(namespace=namespace, name="PlotVCFtoolsStat", version=version, os=operatingSystem, arch=architecture, installed=True)
		PlotVCFtoolsStat.addPFN(PFN("file://" +  os.path.join(self.vervetSrcPath, "plot/PlotVCFtoolsStat.py"), site_handler))
		executableClusterSizeMultiplierList.append((PlotVCFtoolsStat, 0))
		
		SplitVCFFile = Executable(namespace=namespace, name="SplitVCFFile", version=version, os=operatingSystem, arch=architecture, installed=True)
		SplitVCFFile.addPFN(PFN("file://" +  os.path.join(self.pymodulePath, "pegasus/mapper/splitter/SplitVCFFile.py"), site_handler))
		executableClusterSizeMultiplierList.append((SplitVCFFile, 1))
		
		#2012.8.30 moved from vervet/src/AddVCFFolder2DBWorkflow.py
		AddVCFFile2DB = Executable(namespace=namespace, name="AddVCFFile2DB", \
											version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		AddVCFFile2DB.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "db/input/AddVCFFile2DB.py"), site_handler))
		executableClusterSizeMultiplierList.append((AddVCFFile2DB, 1))
		
		FilterVCFSNPCluster = Executable(namespace=namespace, name="FilterVCFSNPCluster", version=version, os=operatingSystem, arch=architecture, installed=True)
		FilterVCFSNPCluster.addPFN(PFN("file://" +  os.path.join(self.pymodulePath, "pegasus/mapper/filter/FilterVCFSNPCluster.py"), site_handler))
		executableClusterSizeMultiplierList.append((FilterVCFSNPCluster, 1))
		
		ExtractSamplesFromVCF = Executable(namespace=namespace, name="ExtractSamplesFromVCF", \
											version=version, os=operatingSystem, arch=architecture, installed=True)
		ExtractSamplesFromVCF.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "mapper/ExtractSamplesFromVCF.py"), \
													site_handler))
		executableClusterSizeMultiplierList.append((ExtractSamplesFromVCF, 1))
		
		JuxtaposeAlleleFrequencyFromMultiVCFInput = Executable(namespace=namespace, name="JuxtaposeAlleleFrequencyFromMultiVCFInput", \
											version=version, os=operatingSystem, arch=architecture, installed=True)
		JuxtaposeAlleleFrequencyFromMultiVCFInput.addPFN(PFN("file://" + os.path.join(self.pymodulePath, \
											"pegasus/mapper/extractor/JuxtaposeAlleleFrequencyFromMultiVCFInput.py"), \
											site_handler))
		executableClusterSizeMultiplierList.append((JuxtaposeAlleleFrequencyFromMultiVCFInput, 1))
		
		ExtractInfoFromVCF = Executable(namespace=namespace, name="ExtractInfoFromVCF", version=version, os=operatingSystem,\
									arch=architecture, installed=True)
		ExtractInfoFromVCF.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/extractor/ExtractInfoFromVCF.py"), site_handler))
		executableClusterSizeMultiplierList.append((ExtractInfoFromVCF, 1))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
	def registerCommonExecutables(self, workflow=None):
		"""
		"""
		AbstractNGSWorkflow.registerCommonExecutables(self, workflow=workflow)
	
	def registerAllInputFiles(self, workflow=None, inputDir=None, input_site_handler=None, \
					checkEmptyVCFByReading=False, pegasusFolderName='',\
					maxContigID=None, minContigID=None, db_vervet=None, needToKnowNoOfLoci=False,
					minNoOfLoci=None, includeIndelVCF=True):
		"""
		2013.3.1 flip includeIndelVCF to true (now indel and SNP vcf files from AlignmentToCall workflows are in separate folders.
		2012.8.15 add argument db_vervet, needToKnowNoOfLoci, to get no_of_loci by parsing inputFname and find db-entry...
			argument minNoOfLoci, only used when it's not None and needToKnowNoOfLoci is True
		2012.8.10 add maxContigID and minContigID to restrict input
		2012.7.27 add attribute file to each object in returnData.jobDataLs
		2012.5.9
			register the tbi file if it exists
		2012.3.1
			moved from CalculateTrioInconsistencyPipeline.py
		2012.1.9
			the returning data structure is changed to conform to some standard used across several functions
		2011-9-29
			vcf files only
		"""
		sys.stderr.write("Registering input files from %s ..."%(inputDir))
		if workflow is None:
			workflow = self
		returnData = PassingData(jobDataLs = [])
		if inputDir and os.path.isdir(inputDir):	#2013.04.07
			fnameLs = os.listdir(inputDir)
			counter = 0
			real_counter = 0
			previous_reported_real_counter = ''
			for fname in fnameLs:
				counter += 1
				inputFname = os.path.realpath(os.path.join(inputDir, fname))
				if (maxContigID is not None and maxContigID!=0) or (minContigID is not None and minContigID!=0):
					try:
						contigID = int(self.getContigIDFromFname(os.path.basename(fname)))
						if (maxContigID is not None and maxContigID!=0) and contigID>maxContigID:
							continue
						if (minContigID is not None and minContigID!=0) and contigID<minContigID:
							continue
					except:
						sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
						import traceback
						traceback.print_exc()
				if NextGenSeq.isFileNameVCF(fname, includeIndelVCF=includeIndelVCF) and \
						not NextGenSeq.isVCFFileEmpty(inputFname, checkContent=checkEmptyVCFByReading):
					real_counter += 1
					inputBaseFname = os.path.basename(inputFname)
					inputF = File(os.path.join(pegasusFolderName, inputBaseFname))
					inputF.addPFN(PFN("file://" + inputFname, input_site_handler))
					inputF.absPath = inputFname
					inputF.abspath = inputFname
					no_of_loci = None
					no_of_individuals = None
					if needToKnowNoOfLoci:
						if db_vervet:
							genotype_file = db_vervet.parseGenotypeFileGivenDBAffiliatedFilename(filename=inputFname)
							if genotype_file and inputFname.find(genotype_file.path)>=0:	#2012.9.6 make sure same file
								no_of_loci = genotype_file.no_of_loci
								no_of_individuals = genotype_file.no_of_individuals
						if no_of_loci is None:
							#do file parsing
							vcfFile = VCFFile(inputFname=inputFname, report=False)
							counter = 0
							no_of_loci = vcfFile.getNoOfLoci()
							no_of_individuals = len(vcfFile.getSampleIDList())
							vcfFile.close()
					inputF.noOfLoci = no_of_loci
					inputF.no_of_loci = no_of_loci
					inputF.no_of_individuals = no_of_individuals
					inputF.noOfIndividuals = no_of_individuals
					
					if minNoOfLoci is None or inputF.no_of_loci is None or (minNoOfLoci and inputF.no_of_loci and  inputF.no_of_loci >=minNoOfLoci):
						workflow.addFile(inputF)
						
						tbi_F_absPath = "%s.tbi"%inputFname
						if os.path.isfile(tbi_F_absPath):	#it exists
							tbi_F = File(os.path.join(pegasusFolderName, "%s.tbi"%inputBaseFname))
							tbi_F.addPFN(PFN("file://" + tbi_F_absPath, input_site_handler))
							tbi_F.abspath = tbi_F_absPath
							workflow.addFile(tbi_F)
						else:
							tbi_F = None
						inputF.tbi_F = tbi_F
						returnData.jobDataLs.append(PassingData(job=None, jobLs=[], \
													vcfFile=inputF, tbi_F=tbi_F, file=inputF, fileLs=[inputF, tbi_F]))
					if real_counter%200==0:
						sys.stderr.write("%s%s"%('\x08'*len(previous_reported_real_counter), real_counter))
						previous_reported_real_counter = repr(real_counter)
		sys.stderr.write("  %s non-empty VCF out of %s files.\n"%(len(returnData.jobDataLs), counter))
		return returnData

	def addPlotVCFtoolsStatJob(self, workflow=None, executable=None, inputFileList=None, outputFnamePrefix=None, \
							whichColumn=None, whichColumnHeader=None, whichColumnPlotLabel=None, need_svg=False, \
							logWhichColumn=True, positiveLog=False, valueForNonPositiveYValue=-1, \
							xColumnPlotLabel=None, xColumnHeader=None, chrLengthColumnHeader=None, chrColumnHeader=None, \
							minChrLength=1000000, minNoOfTotal=100,\
							figureDPI=300, ylim_type=2, samplingRate=0.0001, logCount=False,\
							parentJobLs=None, \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=True, job_max_memory=2000, sshDBTunnel=False, **keywords):
		"""
		2012.10.6 use addGenericDBJob() instead of addGenericJob()
		2012.8.31 add argument positiveLog and valueForNonPositiveYValue
		# whichColumnPlotLabel and xColumnPlotLabel should not contain spaces or ( or ). because they will disrupt shell commandline
		
		2012.8.2 moved from vervet/src/CalculateVCFStatPipeline.py
		2012.8.1
			
			('whichColumn', 0, int): [3, 'w', 1, 'data from this column (index starting from 0) is plotted as y-axis value'],\
			('whichColumnHeader', 0, ): ["", 'W', 1, 'column label (in the header) for the data to be plotted as y-axis value, substitute whichColumn'],\
			('logWhichColumn', 0, int): [0, 'g', 0, 'whether to take log of whichColumn'],\
			('need_svg', 0, ): [0, 'n', 0, 'whether need svg output', ],\
			('whichColumnPlotLabel', 1, ): ['#SNPs in 100kb window', 'D', 1, 'plot label for data of the whichColumn', ],\
			('xColumnPlotLabel', 1, ): ['position', 'x', 1, 'x-axis label (posColumn) in manhattan plot', ],\
			('chrLengthColumnHeader', 1, ): ['chrLength', 'c', 1, 'label of the chromosome length column', ],\
			('chrColumnHeader', 1, ): ['CHR', 'C', 1, 'label of the chromosome column', ],\
			('minChrLength', 1, int): [1000000, 'm', 1, 'minimum chromosome length for one chromosome to be included', ],\
			('xColumnHeader', 1, ): ['BIN_START', 'l', 1, 'label of the position column, BIN_START for binned vcftools output. POS for others.', ],\
			('outputFnamePrefix', 0, ): [None, 'O', 1, 'output filename prefix (optional).'],\
			
				('minNoOfTotal', 1, int): [100, 'i', 1, 'minimum no of total variants (denominator of inconsistent rate)'],\
				('title', 1, ): [None, 't', 1, 'title for the figure.'],\
				('figureDPI', 1, int): [200, 'f', 1, 'dpi for the output figures (png)'],\
				('formatString', 1, ): ['-', '', 1, 'formatString passed to matplotlib plot'],\
				('ylim_type', 1, int): [1, 'y', 1, 'y-axis limit type, 1: 0 to max. 2: min to max'],\
				('samplingRate', 1, float): [0.001, 's', 1, 'how often you include the data'],\
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		if inputFileList:
			extraDependentInputLs.extend(inputFileList)
		extraArgumentList = ["--outputFnamePrefix %s"%outputFnamePrefix, '--minNoOfTotal %s'%(minNoOfTotal), \
							'--figureDPI %s'%(figureDPI), '--ylim_type %s'%(ylim_type), '--samplingRate %s'%(samplingRate), \
							'--xColumnHeader %s'%(xColumnHeader)]
		extraOutputLs = [File('%s.png'%(outputFnamePrefix)), File('%s_hist.png'%(outputFnamePrefix))]
		if need_svg:
			extraOutputLs.append(File('%s.svg'%(outputFnamePrefix)))
		key2ObjectForJob = {}
		if minChrLength is not None:
			extraArgumentList.append('--minChrLength %s'%(minChrLength))
		if whichColumnHeader:
			extraArgumentList.append("--whichColumnHeader %s"%(whichColumnHeader))
		if whichColumn:
			extraArgumentList.append("--whichColumn %s"%(whichColumn))
		if logWhichColumn:
			extraArgumentList.append('--logWhichColumn')
			if positiveLog:
				extraArgumentList.append('--positiveLog')
		if whichColumnPlotLabel:
			extraArgumentList.append("--whichColumnPlotLabel %s"%(whichColumnPlotLabel))
		if xColumnPlotLabel:
			extraArgumentList.append("--xColumnPlotLabel %s"%(xColumnPlotLabel))
		if chrLengthColumnHeader:
			extraArgumentList.append("--chrLengthColumnHeader %s"%(chrLengthColumnHeader))
		if chrColumnHeader:
			extraArgumentList.append("--chrColumnHeader %s"%(chrColumnHeader))
		if logCount:
			extraArgumentList.append("--logCount")
		if valueForNonPositiveYValue:
			extraArgumentList.append("--valueForNonPositiveYValue %s"%(valueForNonPositiveYValue))
		if extraArguments:
			extraArgumentList.append(extraArguments)
		job= self.addGenericDBJob(executable=executable, inputFile=None, outputFile=None, \
				inputFileList=inputFileList, \
				parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, \
				sshDBTunnel=sshDBTunnel, objectWithDBArguments=self, **keywords)
		return job
	
	
	def addSplitVCFFileJob(self, workflow=None, executable=None, inputFile=None, outputFnamePrefix=None, \
					noOfOverlappingSites=1000, noOfSitesPerUnit=5000, noOfTotalSites=10000, \
					parentJobLs=None, \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=True, job_max_memory=2000, sshDBTunnel=False, **keywords):
		"""
		2012.8.26
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		
		#turn them into nonnegative	
		noOfOverlappingSites = abs(noOfOverlappingSites)
		noOfSitesPerUnit = abs(noOfSitesPerUnit)
		noOfTotalSites = abs(noOfTotalSites)
		
		if noOfSitesPerUnit>noOfTotalSites:
			noOfSitesPerUnit = noOfTotalSites
		if noOfOverlappingSites>noOfSitesPerUnit:
			noOfOverlappingSites = noOfSitesPerUnit
		key2ObjectForJob = {}
		extraArgumentList = ["-O %s"%outputFnamePrefix,]
		extraOutputLs = []
		if noOfOverlappingSites is not None:
			extraArgumentList.append('--noOfOverlappingSites %s'%(noOfOverlappingSites))
		if noOfSitesPerUnit is not None:
			extraArgumentList.append('--noOfSitesPerUnit %s'%(noOfSitesPerUnit))
		if noOfTotalSites is not None:
			extraArgumentList.append('--noOfTotalSites %s'%(noOfTotalSites))			
		noOfUnits = max(1, utils.getNoOfUnitsNeededToCoverN(N=noOfTotalSites, s=noOfSitesPerUnit, o=noOfOverlappingSites)-1)
		
		suffixAndNameTupleList = []	# a list of tuples , in each tuple, 1st element is the suffix. 2nd element is the proper name of the suffix.
			#job.$nameFile will be the way to access the file.
			#if 2nd element (name) is missing, suffix[1:].replace('.', '_') is the name (dot replaced by _) 
		for i in xrange(1, noOfUnits+1):
			suffixAndNameTupleList.append(['_unit%s.vcf'%(i), 'unit%s'%(i)])
		if extraArguments:
			extraArgumentList.append(extraArguments)
		self.setupMoreOutputAccordingToSuffixAndNameTupleList(outputFnamePrefix=outputFnamePrefix, suffixAndNameTupleList=suffixAndNameTupleList, \
													extraOutputLs=extraOutputLs, key2ObjectForJob=key2ObjectForJob)
		
		job= self.addGenericJob(executable=executable, inputFile=inputFile, outputFile=None, \
				parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, \
				sshDBTunnel=sshDBTunnel, **keywords)
		return job
	
	def getChr2IntervalDataLsBySelectVCFFile(self, vcfFname=None, noOfSitesPerUnit=5000, noOfOverlappingSites=1000, \
											folderName=None, parentJobLs= None):
		"""
		2012.8.9 update it so that the interval encompassing all lines in one block/unit is known.
			good for mpileup to only work on that interval and then "bcftools view" select from sites from the block.
			TODO: offer partitioning by equal-chromosome span, rather than number of sites.
				Some sites could be in far from each other in one block, which could incur long-running mpileup. goal is to skip these deserts.
		2012.8.8 bugfix add -1 to the starting number below cuz otherwise it's included in the next block's start
				blockStopLineNumber = min(startLineNumber+(i+1)*noOfSitesPerUnit-1, stopLineNumber)	
		2012.8.14
			1.
			2. folderName is the relative path of the folder in the pegasus workflow, that holds vcfFname.
				it'll be created upon file stage-in. no mkdir job for it.
				
			get the number of lines in vcfFname.
			get chr2StartStopDataLsTuple
			for each chr, split its lines into units  that don't exceed noOfSitesPerUnit
				add the split job
				 
		"""
		sys.stderr.write("Splitting %s into blocks of %s lines ... "%(vcfFname, noOfSitesPerUnit))
		#from pymodule import utils
		#noOfLines = utils.getNoOfLinesInOneFileByWC(vcfFname)
		chr2StartStopDataLs = {}
		reader = MatrixFile(vcfFname)
		#csv.reader(inf, delimiter=figureOutDelimiter(vcfFname))
		lineNumber = 0
		previousChromosome = None
		previousLine = None
		chromosome = None
		for row in reader:
			lineNumber += 1
			chromosome, start, stop = row[:3]
			start = int(start)	#0-based, starting base
			stop = int(stop)	#0-based, stopping base but not inclusive, i.e. [start, stop)
			
			if previousLine is None or chromosome!=previousLine.chromosome:	#first line or different chromosome
				if previousLine is not None and previousLine.chromosome is not None:
					
					prevChrLastStartStopData = chr2StartStopDataLs[previousLine.chromosome][-1]
					if prevChrLastStartStopData.stopLineNumber is None:
						prevChrLastStartStopData.stopLineNumber = previousLine.lineNumber
						prevChrLastStartStopData.stopLineStart = previousLine.start
						prevChrLastStartStopData.stopLineStop = previousLine.stop
					
				if chromosome not in chr2StartStopDataLs:
					StartStopData = PassingData(startLineNumber=lineNumber, startLineStart=start, startLineStop=stop, \
											stopLineNumber=None, stopLineStart=None, stopLineStop=None)
					chr2StartStopDataLs[chromosome] = [StartStopData]
			else:	#same chromosome and not first line
				lastStartStopData = chr2StartStopDataLs[chromosome][-1]
				if lastStartStopData.stopLineNumber is None:	#last block hasn't been closed yet.
					noOfLinesInCurrentBlock = lineNumber - lastStartStopData.startLineNumber +1
					if noOfLinesInCurrentBlock>=noOfSitesPerUnit:	#time to close it
						lastStartStopData.stopLineNumber = lineNumber
						lastStartStopData.stopLineStart = start
						lastStartStopData.stopLineStop = stop
				else:	#generate a new block
					StartStopData = PassingData(startLineNumber=lineNumber, startLineStart=start, startLineStop=stop, \
											stopLineNumber=None, stopLineStart=None, stopLineStop=None)
					chr2StartStopDataLs[chromosome].append(StartStopData)
			previousLine = PassingData(chromosome = chromosome, start=start, stop=stop, lineNumber=lineNumber)
		#final closure
		if previousLine is not None:	#vcfFname is not empty
			lastStartStopData = chr2StartStopDataLs[previousLine.chromosome][-1]
			if lastStartStopData.stopLineNumber is None:	#last block hasn't been closed yet.
				#close it regardless of whether it has enough lines in it or not.
				lastStartStopData.stopLineNumber = previousLine.lineNumber
				lastStartStopData.stopLineStart = previousLine.start
				lastStartStopData.stopLineStop = previousLine.stop
		sys.stderr.write("%s chromosomes out of %s lines.\n"%(len(chr2StartStopDataLs), lineNumber))
		
		intervalFile = self.registerOneInputFile(inputFname=vcfFname, folderName=folderName)
		chr2IntervalDataLs = {}
		counter = 0
		for chr, startStopDataLs in chr2StartStopDataLs.iteritems():
			for startStopData in startStopDataLs:
				blockStartLineNumber = startStopData.startLineNumber
				blockStopLineNumber = startStopData.stopLineNumber
				# 2012.8.9 the large interval that encompasses all BED lines 
				interval = '%s:%s-%s'%(chr, startStopData.startLineStart, startStopData.stopLineStop)
				blockIntervalFile = File(os.path.join(folderName, '%s_line_%s_%s_bed.tsv'%(chr, \
												blockStartLineNumber, blockStopLineNumber)))
				blockIntervalJob = self.addSelectLineBlockFromFileJob(executable=self.SelectLineBlockFromFile, \
						inputFile=intervalFile, outputFile=blockIntervalFile,\
						startLineNumber=blockStartLineNumber, stopLineNumber=blockStopLineNumber, \
						parentJobLs=parentJobLs, extraDependentInputLs=None, \
						transferOutput=False, job_max_memory=500)
				intervalFnameSignature = '%s_%s_%s'%(chr, blockStartLineNumber, blockStopLineNumber)
				if chr not in chr2IntervalDataLs:
					chr2IntervalDataLs[chr] = []
				intervalData = PassingData(file=blockIntervalFile, intervalFnameSignature=intervalFnameSignature, interval=interval,\
										chr=chr, jobLs=[blockIntervalJob], job=blockIntervalJob)
				chr2IntervalDataLs[chr].append(intervalData)
				counter += 1
		sys.stderr.write("%s intervals and %s SelectLineBlockFromFile jobs.\n"%(counter, counter))
		return chr2IntervalDataLs
	
	
	def connectDB(self):
		"""
		2012.9.24
			place holder. AbstractVervetMapper.py will use it 
		"""
		self.registerReferenceData = None
		self.refFastaFList= None
	
	def preReduce(self, workflow=None, outputDirPrefix="", passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		
		mapDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, \
										outputDir="%sMap"%(outputDirPrefix))
		passingData.mapDirJob = mapDirJob
		returnData.mapDirJob = mapDirJob
		
		reduceOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, \
												outputDir="%sReduce"%(outputDirPrefix))
		passingData.reduceOutputDirJob = reduceOutputDirJob
		
		returnData.reduceOutputDirJob = reduceOutputDirJob
		
		return returnData
	
	def mapEachChromosome(self, workflow=None, chromosome=None,\
				VCFFile=None, passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def mapEachInterval(self, **keywords):
		"""
		2012.9.22
		"""
		
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def mapEachVCF(self, workflow=None, chromosome=None, passingData=None, transferOutput=False, **keywords):
		"""
		2012.9.22
			default is to split each VCF into intervals
		"""
		
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		
		topOutputDirJob = passingData.topOutputDirJob
		fnamePrefix = passingData.fnamePrefix
		VCFFile = passingData.VCFFile
		jobData = passingData.jobData
		intervalOverlapSize = passingData.intervalOverlapSize
		intervalSize = passingData.intervalSize
		
		outputFnamePrefix = os.path.join(topOutputDirJob.output, '%s_splitVCF'%fnamePrefix)
		
		splitVCFJob = self.addSplitVCFFileJob(executable=self.SplitVCFFile, inputFile=VCFFile, \
											outputFnamePrefix=outputFnamePrefix, \
				noOfOverlappingSites=intervalOverlapSize, noOfSitesPerUnit=intervalSize, noOfTotalSites=VCFFile.noOfLoci, \
				parentJobLs=jobData.jobLs+[topOutputDirJob], \
				extraDependentInputLs=[jobData.tbi_F], \
				extraArguments=None, transferOutput=transferOutput, job_max_memory=2000)
		self.no_of_jobs +=1
		returnData.jobDataLs.append(PassingData(jobLs=[splitVCFJob], file=splitVCFJob.output, \
											fileList=splitVCFJob.outputLs))
		returnData.splitVCFJob = splitVCFJob
		
		return returnData
	
	def linkMapToReduce(self, workflow=None, mapEachIntervalData=None, preReduceReturnData=None, passingData=None, transferOutput=True, **keywords):
		"""
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def reduceEachChromosome(self, workflow=None, chromosome=None, passingData=None, mapEachVCFDataLs=None, reduceEachVCFDataLs=None,\
						transferOutput=True, \
						**keywords):
		"""
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.mapEachVCFDataLs = mapEachVCFDataLs
		returnData.reduceEachVCFDataLs = reduceEachVCFDataLs
		return returnData

	def reduceEachVCF(self, workflow=None, chromosome=None, passingData=None, mapEachIntervalDataLs=None,\
					transferOutput=True, **keywords):
		"""
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.mapEachIntervalDataLs = mapEachIntervalDataLs
		return returnData
	
	def reduce(self, workflow=None, reduceEachChromosomeDataLs=None, \
			mapEachChromosomeDataLs=None, passingData=None, transferOutput=True, \
			**keywords):
		"""
		2012.9.17
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.mapEachChromosomeDataLs = mapEachChromosomeDataLs
		returnData.reduceEachChromosomeDataLs = reduceEachChromosomeDataLs
		return returnData
	
	def addAllJobs(self, workflow=None, inputVCFData=None, chr2IntervalDataLs=None, \
				GenomeAnalysisTKJar=None, samtools=None, \
				CreateSequenceDictionaryJava=None, CreateSequenceDictionaryJar=None, \
				BuildBamIndexFilesJava=None, BuildBamIndexJar=None,\
				mv=None, \
				registerReferenceData=None, \
				needFastaIndexJob=False, needFastaDictJob=False, \
				data_dir=None, no_of_gatk_threads = 1, \
				intervalSize=3000, intervalOverlapSize=0, \
				outputDirPrefix="", transferOutput=True, job_max_memory=2000, **keywords):
		"""
		2012.7.26
			architect of the whole map-reduce framework
		"""
		chr2jobDataLs = {}
		for jobData in inputVCFData.jobDataLs:
			inputF = jobData.file
			chr = self.getChrFromFname(os.path.basename(inputF.name))
			if chr not in chr2jobDataLs:
				chr2jobDataLs[chr] = []
			chr2jobDataLs[chr].append(jobData)
		
		sys.stderr.write("Adding jobs for %s chromosomes/contigs of %s VCF files..."%(len(chr2jobDataLs), len(inputVCFData.jobDataLs)))
		refFastaFList = registerReferenceData.refFastaFList
		if refFastaFList:
			refFastaF = refFastaFList[0]
		else:
			refFastaF = None
		
		topOutputDir = "%sMap"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		
		if needFastaDictJob or registerReferenceData.needPicardFastaDictJob:
			fastaDictJob = self.addRefFastaDictJob(workflow, CreateSequenceDictionaryJava=CreateSequenceDictionaryJava, \
												refFastaF=refFastaF)
			refFastaDictF = fastaDictJob.refFastaDictF
		else:
			fastaDictJob = None
			refFastaDictF = registerReferenceData.refPicardFastaDictF
		
		if needFastaIndexJob or registerReferenceData.needSAMtoolsFastaIndexJob:
			fastaIndexJob = self.addRefFastaFaiIndexJob(workflow, samtools=samtools, refFastaF=refFastaF)
			refFastaIndexF = fastaIndexJob.refFastaIndexF
		else:
			fastaIndexJob = None
			refFastaIndexF = registerReferenceData.refSAMtoolsFastaIndexF
		
		returnData = PassingData()
		returnData.jobDataLs = []
		
		#2012.9.22 
		# 	mapEachAlignmentDataLs is never reset.
		#	mapEachChromosomeDataLs is reset right after a new alignment is chosen.
		#	mapEachIntervalDataLs is reset right after each chromosome is chosen.
		#	all reduce dataLs never gets reset.
		passingData = PassingData(fnamePrefix=None, \
					chromosome=None, \
					
					topOutputDirJob=topOutputDirJob, \
					outputDirPrefix=outputDirPrefix, \
					intervalFnamePrefix=None,\
					
					refFastaFList=refFastaFList, \
					registerReferenceData=registerReferenceData, \
					refFastaF=refFastaFList[0],\
					
					fastaDictJob = fastaDictJob,\
					refFastaDictF = refFastaDictF,\
					fastaIndexJob = fastaIndexJob,\
					refFastaIndexF = refFastaIndexF,\
					
					intervalOverlapSize =intervalOverlapSize, intervalSize=intervalSize,\
					jobData=None,\
					VCFFile=None,\
					splitVCFFile=None,\
					preReduceReturnData=None,\
					
					mapEachIntervalData=None,\
					mapEachVCFData=None,\
					mapEachChromosomeData=None, \
					
					mapEachIntervalDataLs=None,\
					mapEachVCFDataLs=None,\
					
					mapEachIntervalDataLsLs=[],\
					mapEachVCFDataLsLs=[],\
					
					mapEachChromosomeDataLs=[], \
					
					reduceEachVCFData=None,\
					reduceEachChromosomeData=None,\
					
					reduceEachVCFDataLs=None,\
					
					reduceEachVCFDataLsLs=[],\
					
					reduceEachChromosomeDataLs=[],\
					
					chr2jobDataLs = chr2jobDataLs,\
					)
		# mapEachIntervalDataLsLs is list of mapEachIntervalDataLs by each VCF file.
		# mapEachVCFDataLsLs is list of mapEachVCFDataLs by each chromosome
		# reduceEachVCFDataLsLs is list of reduceEachVCFDataLs by each chromosome
		
		preReduceReturnData = self.preReduce(workflow=workflow, outputDirPrefix=outputDirPrefix, \
									passingData=passingData, transferOutput=False,\
									**keywords)
		passingData.preReduceReturnData = preReduceReturnData
		
		#gzip folder jobs (to avoid repeatedly creating the same folder
		gzipReduceEachVCFFolderJob = None
		gzipReduceEachChromosomeFolderJob = None
		gzipReduceFolderJob = None
		gzipPreReduceFolderJob = None
		for chr, jobDataLs in chr2jobDataLs.iteritems():
			passingData.chromosome = chr
			mapEachChromosomeData = self.mapEachChromosome(workflow=workflow, chromosome=chr, \
										passingData=passingData, \
										transferOutput=False, **keywords)
			passingData.mapEachChromosomeData = mapEachChromosomeData
			passingData.mapEachChromosomeDataLs.append(mapEachChromosomeData)
			
			passingData.mapEachVCFDataLsLs.append([])
			#the last one from the double list is the current one
			passingData.mapEachVCFDataLs = passingData.mapEachVCFDataLsLs[-1]
			
			passingData.reduceEachVCFDataLsLs.append([])
			passingData.reduceEachVCFDataLs = passingData.reduceEachVCFDataLsLs[-1]
			
			for i in xrange(len(jobDataLs)):
				jobData = jobDataLs[i]
				passingData.jobData = jobData
				
				VCFFile = jobData.vcfFile
				inputFBaseName = os.path.basename(VCFFile.name)
				commonPrefix = inputFBaseName.split('.')[0]
				
				passingData.fnamePrefix = commonPrefix
				passingData.VCFFile = VCFFile
				
				mapEachVCFData = self.mapEachVCF(workflow=workflow,\
												VCFFile=VCFFile, passingData=passingData, \
												transferOutput=False, **keywords)
				passingData.mapEachVCFData = mapEachVCFData
				passingData.mapEachVCFDataLs.append(mapEachVCFData)
				
				passingData.mapEachIntervalDataLsLs.append([])
				passingData.mapEachIntervalDataLs = passingData.mapEachIntervalDataLsLs[-1]
				
				noOfUnits = max(1, utils.getNoOfUnitsNeededToCoverN(N=inputF.noOfLoci, s=intervalSize, o=intervalOverlapSize)-1)
				for unitNumber in xrange(1, noOfUnits+1):
					splitVCFJob = mapEachVCFData.splitVCFJob
					splitVCFFile = getattr(splitVCFJob, 'unit%sFile'%(unitNumber), None)
					if splitVCFFile is not None:
						passingData.splitVCFFile = splitVCFFile
						passingData.unitNumber = unitNumber
						passingData.intervalFnamePrefix = '%s_%s_splitVCF_u%s'%(chr, commonPrefix, unitNumber)
						
						mapEachIntervalData = self.mapEachInterval(workflow=workflow, \
												VCFFile=splitVCFFile, passingData=passingData, transferOutput=False, \
												**keywords)
						passingData.mapEachIntervalData = mapEachIntervalData
						passingData.mapEachIntervalDataLs.append(mapEachIntervalData)
						
						linkMapToReduceData = self.linkMapToReduce(workflow=workflow, mapEachIntervalData=mapEachIntervalData, \
											preReduceReturnData=preReduceReturnData, \
											passingData=passingData, \
											**keywords)
				
				reduceEachVCFData = self.reduceEachVCF(workflow=workflow, chromosome=chr, passingData=passingData, \
								mapEachIntervalDataLs=passingData.mapEachIntervalDataLs,\
								transferOutput=False, data_dir=data_dir, \
								**keywords)
				passingData.reduceEachVCFData = reduceEachVCFData
				passingData.reduceEachVCFDataLs.append(reduceEachVCFData)
				
				gzipReduceEachVCFData = self.addGzipSubWorkflow(workflow=workflow, \
					inputData=reduceEachVCFData, transferOutput=transferOutput,\
					outputDirPrefix="%sReduceEachVCF"%(outputDirPrefix), topOutputDirJob=gzipReduceEachVCFFolderJob, \
					report=False)
				gzipReduceEachVCFFolderJob = gzipReduceEachVCFData.topOutputDirJob
			reduceEachChromosomeData = self.reduceEachChromosome(workflow=workflow, chromosome=chr, passingData=passingData, \
								mapEachVCFDataLs=passingData.mapEachVCFDataLs,\
								reduceEachVCFDataLs=passingData.reduceEachVCFDataLs,\
								transferOutput=False, data_dir=data_dir, \
								**keywords)
			passingData.reduceEachChromosomeData = reduceEachChromosomeData
			passingData.reduceEachChromosomeDataLs.append(reduceEachChromosomeData)
			
			gzipReduceEachChromosomeData = self.addGzipSubWorkflow(workflow=workflow, \
					inputData=reduceEachChromosomeData, transferOutput=transferOutput,\
					outputDirPrefix="%sReduceEachChromosome"%(outputDirPrefix), \
					topOutputDirJob=gzipReduceEachChromosomeFolderJob, report=False)
			gzipReduceEachChromosomeFolderJob = gzipReduceEachChromosomeData.topOutputDirJob
			
		reduceReturnData = self.reduce(workflow=workflow, passingData=passingData, transferOutput=False, \
							mapEachChromosomeDataLs=passingData.mapEachVCFDataLs,\
							reduceEachChromosomeDataLs=passingData.reduceEachChromosomeDataLs,\
							**keywords)
		passingData.reduceReturnData = reduceReturnData
		
		
		gzipPreReduceReturnData = self.addGzipSubWorkflow(workflow=workflow, inputData=preReduceReturnData, transferOutput=transferOutput,\
						outputDirPrefix="%sPreReduce"%(outputDirPrefix), \
						topOutputDirJob= gzipPreReduceFolderJob, report=False)
		gzipPreReduceFolderJob = gzipPreReduceReturnData.topOutputDirJob
		
		gzipReduceReturnData = self.addGzipSubWorkflow(workflow=workflow, inputData=reduceReturnData, transferOutput=transferOutput,\
						outputDirPrefix="%sReduce"%(outputDirPrefix), \
						topOutputDirJob=gzipReduceFolderJob, report=False)
		gzipReduceFolderJob = gzipReduceReturnData.topOutputDirJob
		sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
		return reduceReturnData
	
	def setup_run(self):
		"""
		2013.04.07 wrap all standard pre-run() related functions into this function.
			setting up for run(), called by run()
		"""
		pdata = AbstractNGSWorkflow.setup_run(self)
		workflow = pdata.workflow
		
		#self.chr2size = {}
		#self.chr2size = set(['Contig149'])	#temporary when testing Contig149
		#self.chr2size = set(['1MbBAC'])	#temporary when testing the 1Mb-BAC (formerly vervet_path2)
		chrLs = self.chr2size.keys()
		chr2IntervalDataLs = self.getChr2IntervalDataLsBySplitChrSize(chr2size=self.chr2size, \
													intervalSize=self.intervalSize, \
													intervalOverlapSize=self.intervalOverlapSize)
		
		inputData = self.registerAllInputFiles(workflow, self.inputDir, input_site_handler=self.input_site_handler, \
											checkEmptyVCFByReading=self.checkEmptyVCFByReading,\
											pegasusFolderName=self.pegasusFolderName,\
											maxContigID=self.maxContigID, \
											minContigID=self.minContigID,\
											db_vervet=getattr(self, 'db_vervet', None), \
											needToKnowNoOfLoci=True,\
											minNoOfLoci=getattr(self, 'minNoOfLoci', 10))
		
		registerReferenceData = self.getReferenceSequence()
		
		
		return PassingData(workflow=workflow, inputData=inputData,\
						chr2IntervalDataLs=chr2IntervalDataLs, registerReferenceData=registerReferenceData)
	
	def run(self):
		"""
		2012.9.24
		"""
		pdata = self.setup_run()
		workflow = pdata.workflow
		
		
		inputData=pdata.inputData
		
		if len(inputData.jobDataLs)<=0:
			sys.stderr.write("No VCF files in this folder , %s.\n"%self.inputDir)
			sys.exit(0)
				
		self.addAllJobs(workflow=workflow, inputVCFData=inputData, \
					chr2IntervalDataLs=None, samtools=workflow.samtools, \
				GenomeAnalysisTKJar=workflow.GenomeAnalysisTKJar, \
				CreateSequenceDictionaryJava=workflow.CreateSequenceDictionaryJava, CreateSequenceDictionaryJar=workflow.CreateSequenceDictionaryJar, \
				BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexJar=workflow.BuildBamIndexJar,\
				mv=workflow.mv, \
				registerReferenceData=pdata.registerReferenceData,\
				needFastaIndexJob=getattr(self, 'needFastaIndexJob',False), needFastaDictJob=getattr(self, 'needFastaDictJob', False), \
				data_dir=self.data_dir, no_of_gatk_threads = 1,\
				intervalSize=self.intervalSize, intervalOverlapSize=self.intervalOverlapSize, \
				outputDirPrefix=self.pegasusFolderName, transferOutput=True,)
		
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)