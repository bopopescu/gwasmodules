#!/usr/bin/env python
"""
Examples:
	# 2011-10-10 call method 80, analysis method 1, min_score 5
	%s  --call_method_id_ls 80 --analysis_method_id_ls 1 --min_score 5
		-o call80analysis1MinScore5.xml -u yh -z banyan -l condorpool
	
	# 2011-10-16 call method 57 or cnv method 20, analysis method 1, min score 4
	%s --call_method_id_ls 57 --cnv_method_id 20 --analysis_method_id_ls 1 --min_score 4
		-o call57cnv20analysis1MinScore4.xml -u yh -z banyan -l condorpool
	
	#ditto but analysis method 7 and min score 3 (make sure -c is added, otherwise nothing will be stored in db)
	%s --call_method_id_ls 57 --cnv_method_id 20 --analysis_method_id_ls 7 --min_score 3
		-o call57cnv20analysis7MinScore3.xml -u yh -z banyan -j condorpool -l condorpool
	
	#ditto but analysis method 1,32
	%s  --call_method_id_ls 57,75 --analysis_method_id_ls 1,32 --min_score_ls 3,4,5 --min_overlap_ratio_ls 0.05,0.20,0.50,0.80
		--phenotype_method_id_ls 1-35,39-48,57-82,158-159,161-179,182-186,272-283,314-351,362-380,418-589
		-o workflow/BootstrapAssociationPeak/Bootstrap_call57_75_ana_1_32_min_score_3_4_5_min_overlap_0.05_0.2_0.5_0.8.xml
		-u yh  --db_passwd secret -z banyan -j condorpool -l condorpool

Description:
	2012.11.21 a workflow that checks association landscape/peak/loci under different threshold settings.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])


#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:	   #64bit
#	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
#	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
#else:   #32bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus

from variation.src import AbstractVariationWorkflow
from variation.src import Stock_250kDB
from variation.src import DefineAssociationLandscapePipeline

class BootstrapAssociationPeakWorkflow(DefineAssociationLandscapePipeline):
	__doc__ = __doc__
	option_default_dict = DefineAssociationLandscapePipeline.option_default_dict.copy()
	for key in [('result_id_ls', 0, )]:
		option_default_dict.pop(key)
	#remove some irrelevant options
	ProcessOptions.removeCertainOptions(option_default_dict, option_dict_to_remove=DefineAssociationLandscapePipeline.my_option_dict)
	
	option_default_dict.update({
						('min_score_ls', 0, ): ['4', 'f', 1, 'minimum score to call a peak'],\
						('min_overlap_ratio_ls', 1, ): ['0.05', '', 1, 'minimum overlap ratio, overlap length/total' ],\
						('peakPadding', 0, int): [0, '', 1, 'the padding around each peak (use only to extend the overlap between peaks)' ],\
						})
	
	def __init__(self,  **keywords):
		"""
		2011-10
		"""
		DefineAssociationLandscapePipeline.__init__(self, **keywords)
		
		
		listArgumentName_data_type_ls = [("min_score_ls", float), ('min_overlap_ratio_ls', float)]
		listArgumentName2hasContent = ProcessOptions.processListArguments(listArgumentName_data_type_ls, emptyContent=[],\
																class_to_have_attr=self)
	
	def addTwoAssociationLocusFileOverlapJob(self, workflow=None, executable=None, \
				inputFileList=None, inputFile=None, outputFile=None, \
				parentJobLs=None, job_max_memory=100, job_max_walltime = 60, \
				extraDependentInputLs=None, \
				transferOutput=False, **keywords):
		"""
		2012.11.28
			job_max_walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
		"""
		extraArgumentList = []
		job = self.addAbstractMatrixFileWalkerJob(workflow=workflow, executable=executable, inputFileList=inputFileList, \
								inputFile=inputFile, outputFile=outputFile, \
					outputFnamePrefix=None, whichColumn=None, whichColumnHeader=None, \
					logY=None, valueForNonPositiveYValue=None, \
					minNoOfTotal=None,\
					samplingRate=None, \
					inputFileFormat=2, outputFileFormat=2,\
					parentJobLs=parentJobLs, \
					extraDependentInputLs=extraDependentInputLs, extraArgumentList=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=job_max_memory, \
					job_max_walltime=job_max_walltime,\
					sshDBTunnel=False, \
					objectWithDBArguments=None, **keywords)
		return job
		
	
	def preReduce(self, workflow=None, outputDirPrefix="", passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		returnData = PassingData()
		returnData.jobDataLs = []
		
		mapDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir="%sMap"%(outputDirPrefix))
		passingData.mapDirJob = mapDirJob
		returnData.mapDirJob = mapDirJob
		
		countAssociationLocusOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, \
											outputDir="%sCountAssociationLocus"%(outputDirPrefix))
		passingData.countAssociationLocusOutputDirJob = countAssociationLocusOutputDirJob
		returnData.countAssociationLocusOutputDirJob = countAssociationLocusOutputDirJob
		
		plotOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, \
											outputDir='%sPlot'%(outputDirPrefix))
		passingData.plotOutputDirJob = plotOutputDirJob
		returnData.plotOutputDirJob = plotOutputDirJob
		
		reduceOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, \
												outputDir="%sReduce"%(outputDirPrefix))
		passingData.reduceOutputDirJob = reduceOutputDirJob
		returnData.reduceOutputDirJob = reduceOutputDirJob
		
		return returnData
	
	def addAllJobs(self, workflow=None, association_result_ls=None, data_dir=None, min_MAF=None, \
				neighbor_distance=None, max_neighbor_distance=None, \
				min_score_ls=None, min_overlap_ratio_ls=None, ground_score=None,\
				peakPadding=None, tax_id=None, \
				outputDirPrefix="", transferOutput=True, job_max_memory=2000, **keywords):
		"""
		2012.11.21
			
		"""
		if workflow is None:
			workflow = self
		sys.stderr.write("Adding jobs for %s association results ..."%(len(association_result_ls)))
		
		returnData = PassingData()
		returnData.jobDataLs = []
		
		passingData = PassingData(fnamePrefix=None, \
					outputDirPrefix=outputDirPrefix, \
					jobData=None,\
					preReduceReturnData=None,\
					)
		
		preReduceReturnData = self.preReduce(workflow=workflow, outputDirPrefix=outputDirPrefix, \
									passingData=passingData, transferOutput=False,\
									**keywords)
		
		mapDirJob = preReduceReturnData.mapDirJob
		plotOutputDirJob = preReduceReturnData.plotOutputDirJob
		countAssociationLocusOutputDirJob = preReduceReturnData.countAssociationLocusOutputDirJob
		reduceOutputDirJob = preReduceReturnData.reduceOutputDirJob
		
		passingData.preReduceReturnData = preReduceReturnData
		
		association_group_key2orderIndex = {}
		association_group_key2resultList = {}
		association_group_key2reduceAssociationPeakJobMatrix = {}
		association_group_key2countAssociationLocusJobList = {}
		
		#the job matrix is a matrix of AssociationPeak2AssociationLocus jobs.
		#	each row is index by association threshold.
		#	each column is index by overlap threshold.
		
		resultID2associationResultFile = {}
		for result in association_result_ls:
			associationResultFile = self.registerOneInputFile(inputFname=result.getFilePath(oldDataDir=self.db.data_dir, newDataDir=self.data_dir), \
										folderName=self.pegasusFolderName)
			resultID2associationResultFile[result.id] = associationResultFile
			
			call_method_id = result.call_method_id
			analysis_method_id = result.analysis_method_id
			association_group_key = (call_method_id, analysis_method_id)
			if association_group_key not in association_group_key2resultList:
				association_group_key2resultList[association_group_key] = []
				association_group_key2reduceAssociationPeakJobMatrix[association_group_key] = []
				for min_score in min_score_ls:
					association_group_key2reduceAssociationPeakJobMatrix[association_group_key].append([])
				association_group_key2countAssociationLocusJobList[association_group_key] = []
			association_group_key2resultList[association_group_key].append(result)
		
		
		for i in xrange(len(min_score_ls)):
			min_score = min_score_ls[i]
			#create a folder here.
			associationMinScoreDir = "association_min_score_%s"%(min_score)
			associationMinScoreDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=associationMinScoreDir)
			
			for j in xrange(len(min_overlap_ratio_ls)):
				min_overlap_ratio = min_overlap_ratio_ls[j]
				#add PlotAssociationLocusFrequencyOnGenome job
				associationLocusFrequencyOnGenomeFnamePrefix = os.path.join(plotOutputDirJob.output, 'frequency_manhattan_min_score_%s_min_overlap_%s'%\
											(min_score, min_overlap_ratio))
				outputFile = File('%s.png'%(associationLocusFrequencyOnGenomeFnamePrefix))
				plotAssociationLocusFrequencyOnGenomeJob = self.addAbstractPlotJob(executable=self.PlotAssociationLocusFrequencyOnGenome, \
					inputFileList=None, inputFile=None, outputFile=outputFile, \
					outputFnamePrefix=None, whichColumn=None, whichColumnHeader="no_of_results", whichColumnPlotLabel="numberOfResults", \
					logX=None, logY=None, valueForNonPositiveYValue=-1, \
					xScaleLog=None, yScaleLog=1,\
					missingDataNotation='NA',\
					xColumnHeader="start", xColumnPlotLabel="genomePosition", \
					minNoOfTotal=1, maxNoOfTotal=None,\
					figureDPI=300, formatString='.', ylim_type=2, samplingRate=1, need_svg=False, \
					inputFileFormat=None, outputFileFormat=None,\
					parentJobLs=[plotOutputDirJob], \
					extraDependentInputLs=None, \
					extraArgumentList=None, extraArguments=None, transferOutput=True,  job_max_memory=2000, \
					sshDBTunnel=self.needSSHDBTunnel, \
					objectWithDBArguments=self)
				
				#for HistogramAssociationLocusAdjacencyDistance,
				#comparing association locus from two different call methods, but same analysis method
				analysis_method_id2AssociationLocusJobList = {}
				
				for association_group_key, result_ls in association_group_key2resultList.iteritems():
					#add a AssociationPeak2AssociationLocus job
					associationLocusFnamePrefix = 'call_%s_analysis_%s_min_score_%s_min_overlap_%s'%\
											(association_group_key[0], association_group_key[1], min_score, min_overlap_ratio)
					associationLocusFile = File(os.path.join(associationMinScoreDirJob.output, '%s_locus.h5'%(associationLocusFnamePrefix)))
					#input to this job is added later
					associationLocusJob = self.addAssociationPeak2LocusJob(executable=self.AssociationPeak2AssociationLocus, \
								inputFile=None, outputFile=associationLocusFile, min_overlap_ratio=min_overlap_ratio, \
								parentJobLs=[associationMinScoreDirJob], job_max_memory=2000, job_max_walltime = 60, \
								extraDependentInputLs=None, \
								transferOutput=False)
					
					call_method_id, analysis_method_id = association_group_key[:2]
					associationLocusJob.call_method_id = call_method_id
					associationLocusJob.analysis_method_id = analysis_method_id
					
					association_group_key2reduceAssociationPeakJobMatrix[association_group_key][i].append(associationLocusJob)
					
					self.addInputToStatMergeJob(statMergeJob=plotAssociationLocusFrequencyOnGenomeJob, parentJobLs=[associationLocusJob])
					
					if analysis_method_id not in analysis_method_id2AssociationLocusJobList:
						analysis_method_id2AssociationLocusJobList[analysis_method_id] = []
					analysis_method_id2AssociationLocusJobList[analysis_method_id].append(associationLocusJob)
					
					for result in result_ls:
						associationResultFile = resultID2associationResultFile.get(result.id)
						#associationResultFile = self.registerOneInputFile(inputFname=result.getFilePath(oldDataDir=self.db.data_dir, newDataDir=self.data_dir), \
						#				folderName=self.pegasusFolderName)
						#add DefineAssociationLandscape job
						outputFnamePrefix = '%s_result_%s'%(associationLocusFnamePrefix, result.id)
						landscapeOutputFile = File(os.path.join(associationMinScoreDirJob.output, '%s_landscape.h5'%(outputFnamePrefix)))
						
						defineLandscapeJob = self.addDefineLandscapeJob(workflow, executable=workflow.DefineAssociationLandscape, \
										result_id=result.id, neighbor_distance=neighbor_distance, \
										max_neighbor_distance=max_neighbor_distance,\
										min_MAF=min_MAF, tax_id=tax_id, \
										data_dir=data_dir, logFile=None,\
										landscapeOutputFile=landscapeOutputFile,\
										extraDependentInputLs=[associationResultFile], \
										parentJobLs=[associationMinScoreDirJob], sshDBTunnel=self.needSSHDBTunnel,\
										job_max_memory=2000, transferOutput=False)
						#add landscape -> peak job
						outputFile = File(os.path.join(associationMinScoreDirJob.output, '%s_peak.h5'%(outputFnamePrefix)))
						associationPeakJob = self.addAssociationLandscape2PeakJob(executable=self.AssociationLandscape2Peak, \
							inputFile=defineLandscapeJob.output, outputFile=outputFile, min_score=min_score, ground_score=ground_score, \
							parentJobLs=[defineLandscapeJob], job_max_memory=2000, job_max_walltime = 60, \
							extraDependentInputLs=None, \
							transferOutput=False)
						self.addInputToStatMergeJob(statMergeJob=associationLocusJob, parentJobLs=[associationPeakJob])
				
				#add HistogramAssociationLocusAdjacencyDistance jobs
				for analysis_method_id, associationLocusJobList in analysis_method_id2AssociationLocusJobList.iteritems():
					#add HistogramAssociationLocusAdjacencyDistance job for each pair of association locus jobs
					for k in xrange(len(associationLocusJobList)-1):
						associationLocusJob1 = associationLocusJobList[k]
						for l in xrange(k+1, len(associationLocusJobList)):
							associationLocusJob2 = associationLocusJobList[l]
							
							outputFnamePrefix = 'distance_from_call_%s_to_call_%s_ana_%s_loci_min_score_%s_min_overlap_%s'%\
										(associationLocusJob1.call_method_id, associationLocusJob2.call_method_id, \
										analysis_method_id, min_score, min_overlap_ratio)
							outputFile = File(os.path.join(associationMinScoreDirJob.output, '%s.h5'%(outputFnamePrefix)))
							
							twoAssociationLocusOverlapJob = self.addTwoAssociationLocusFileOverlapJob(executable=self.TwoAssociationLocusFileOverlap, \
														inputFileList=[associationLocusJob1.output, associationLocusJob2.output], \
														inputFile=None, outputFile=outputFile, \
														parentJobLs=[associationLocusJob1, associationLocusJob2], \
														job_max_memory=1000, job_max_walltime = 60, \
														transferOutput=False)
							#add HistogramAssociationLocusAdjacencyDistance job. the order of input matters
							outputFile = File(os.path.join(plotOutputDirJob.output, 'histogram_of_%s.png'%(outputFnamePrefix)))
							histogramAssociationLocusAdjacencyDistanceJob = self.addDrawHistogramJob(executable=self.DrawHistogram, \
								inputFileList=None, inputFile=twoAssociationLocusOverlapJob.output, outputFile=outputFile, \
								outputFnamePrefix=None, whichColumn=None, whichColumnHeader="fractionCoveredByAssociation2", \
								whichColumnPlotLabel="fractionCall%sCoveredByCall%s"%(associationLocusJob1.call_method_id, associationLocusJob2.call_method_id), \
								logY=1, valueForNonPositiveYValue=-10, \
								missingDataNotation='NA',\
								minNoOfTotal=1, maxNoOfTotal=None,\
								figureDPI=300, formatString='.', ylim_type=2, samplingRate=1, need_svg=False, \
								logCount=True, inputFileFormat=2, \
								parentJobLs=[plotOutputDirJob, twoAssociationLocusOverlapJob], \
								extraDependentInputLs=None, \
								extraArgumentList=None, extraArguments=None, transferOutput=True, job_max_memory=1000)
							
							#now reverse the order
							outputFnamePrefix = 'distance_from_call_%s_to_call_%s_ana_%s_loci_min_score_%s_min_overlap_%s'%\
											(associationLocusJob2.call_method_id, associationLocusJob1.call_method_id, \
											analysis_method_id, min_score, min_overlap_ratio)
							outputFile = File(os.path.join(associationMinScoreDirJob.output, '%s.h5'%(outputFnamePrefix)))
							
							twoAssociationLocusOverlapJob = self.addTwoAssociationLocusFileOverlapJob(executable=self.TwoAssociationLocusFileOverlap, \
														inputFileList=[associationLocusJob2.output, associationLocusJob1.output], \
														inputFile=None, outputFile=outputFile, \
														parentJobLs=[associationLocusJob1, associationLocusJob2], \
														job_max_memory=1000, job_max_walltime = 60, \
														transferOutput=False)
							#add HistogramAssociationLocusAdjacencyDistance job. the order of input matters
							outputFile = File(os.path.join(plotOutputDirJob.output, 'histogram_of_%s.png'%(outputFnamePrefix)))
							histogramAssociationLocusAdjacencyDistanceJob = self.addDrawHistogramJob(executable=self.DrawHistogram, \
								inputFileList=None, inputFile=twoAssociationLocusOverlapJob.output, outputFile=outputFile, \
								outputFnamePrefix=None, whichColumn=None, whichColumnHeader="fractionCoveredByAssociation2", \
								whichColumnPlotLabel="fractionCall%sCoveredByCall%s"%(associationLocusJob2.call_method_id, associationLocusJob1.call_method_id), \
								logY=1, valueForNonPositiveYValue=-10, \
								missingDataNotation='NA',\
								minNoOfTotal=1, maxNoOfTotal=None,\
								figureDPI=300, formatString='.', ylim_type=2, samplingRate=1, need_svg=False, \
								logCount=True, inputFileFormat=2, \
								parentJobLs=[plotOutputDirJob, twoAssociationLocusOverlapJob], \
								extraDependentInputLs=None, \
								extraArgumentList=None, extraArguments=None, transferOutput=True, job_max_memory=1000)
							
		sys.stderr.write("\t %s jobs added before PlotAssociationLocusFrequencyVsAssociationThreshold jobs.\n"%(self.no_of_jobs))
		#add a PlotAssociationLocusFrequencyVsAssociationThreshold job
		for j in xrange(len(min_overlap_ratio_ls)):
			min_overlap_ratio = min_overlap_ratio_ls[j]
			outputFnamePrefix = os.path.join(plotOutputDirJob.output, 'locus_frequency_vs_association_threshold_min_overlap_%s'%\
										(min_overlap_ratio))
			outputFile = File('%s.png'%(outputFnamePrefix))
			plotAssociationLocusFrequencyVsAssociationThresholdJob = self.addAbstractPlotJob(executable=self.PlotAssociationLocusFrequencyVsAssociationThreshold, \
					inputFileList=None, inputFile=None, outputFile=outputFile, \
					outputFnamePrefix=None, whichColumn=None, whichColumnHeader="no_of_association_loci", whichColumnPlotLabel="numberOfPeaks", \
					logX=None, logY=None, valueForNonPositiveYValue=-1, \
					yScaleLog=1,\
					missingDataNotation='NA',\
					xColumnHeader="min_score", xColumnPlotLabel=None, \
					minNoOfTotal=1, maxNoOfTotal=None,\
					figureDPI=300, formatString='.-', ylim_type=2, samplingRate=1, need_svg=False, \
					inputFileFormat=2, outputFileFormat=None,\
					parentJobLs=[plotOutputDirJob], \
					extraDependentInputLs=None, \
					extraArgumentList=None, extraArguments=None, transferOutput=True,  job_max_memory=2000, \
					sshDBTunnel=None, \
					objectWithDBArguments=None)
			
			for association_group_key, countAssociationLocusJobList in association_group_key2countAssociationLocusJobList.iteritems():
				#add a CountAssociationLocus job across different association thresholds (min_score)
				countAssociationLocusFnamePrefix = os.path.join(countAssociationLocusOutputDirJob.output, 'call_%s_analysis_%s_min_overlap_%s'%\
										(association_group_key[0], association_group_key[1], min_overlap_ratio))
				countAssociationLocusFile = File('%s_loci_count.h5'%(countAssociationLocusFnamePrefix))
				#input to this job is added later
				countAssociationLocusJob = self.addCountAssociationLocusJob(executable=self.CountAssociationLocus, \
								inputFileList=None, inputFile=None, outputFile=countAssociationLocusFile, \
								parentJobLs=[countAssociationLocusOutputDirJob], job_max_memory=1000, job_max_walltime = 6, \
								extraDependentInputLs=None, \
								transferOutput=False)
				countAssociationLocusJobList.append(countAssociationLocusJob)
				for i in xrange(len(min_score_ls)):
					min_score = min_score_ls[i]
					associationLocusJob = association_group_key2reduceAssociationPeakJobMatrix[association_group_key][i][j]
					self.addInputToStatMergeJob(statMergeJob=countAssociationLocusJob, parentJobLs=[associationLocusJob])
				
				self.addInputToStatMergeJob(statMergeJob=plotAssociationLocusFrequencyVsAssociationThresholdJob, \
										parentJobLs=[countAssociationLocusJob])
		sys.stderr.write("\t %s jobs added before BoxPlotAssociationLocusAttributeVsOverlapThreshold jobs.\n"%(self.no_of_jobs))
		
		#add the BoxPlotAssociationLocusAttributeVsOverlapThreshold jobs
		for association_group_key, jobMatrix in association_group_key2reduceAssociationPeakJobMatrix.iteritems():
			for i in xrange(len(min_score_ls)):
				min_score = min_score_ls[i]
				#add a BoxPlotAssociationLocusAttributeVsOverlapThreshold job
				outputFnamePrefix = os.path.join(plotOutputDirJob.output, 'call_%s_analysis_%s_boxplot_association_locus_connectivity_vs_overlap_threshold_min_score_%s'%\
							(association_group_key[0], association_group_key[1], min_score))
				outputFile = File('%s.png'%(outputFnamePrefix))
				boxPlotAssociationLocusAttributeVsOverlapThresholdJob = self.addAbstractPlotJob(executable=self.BoxPlotAssociationLocusAttributeVsOverlapThreshold, \
					inputFileList=None, inputFile=None, outputFile=outputFile, \
					outputFnamePrefix=None, whichColumn=None, whichColumnHeader="connectivity", whichColumnPlotLabel="connectivity", \
					logX=None, logY=None, valueForNonPositiveYValue=-1, \
					missingDataNotation='NA',\
					xColumnHeader=None, xColumnPlotLabel=None, title="call-%s-analysis-%s-min_score-%s"%\
							(association_group_key[0], association_group_key[1], min_score), \
					minNoOfTotal=1, maxNoOfTotal=None,\
					figureDPI=300, formatString='.', ylim_type=2, samplingRate=1, need_svg=False, \
					inputFileFormat=2, outputFileFormat=None,\
					parentJobLs=[plotOutputDirJob], \
					extraDependentInputLs=None, \
					extraArgumentList=None, extraArguments=None, transferOutput=True,  job_max_memory=2000, \
					sshDBTunnel=None, \
					objectWithDBArguments=None)
				
				for j in xrange(len(min_overlap_ratio_ls)):
					min_overlap_ratio = min_overlap_ratio_ls[j]
					associationLocusJob = jobMatrix[i][j]
					self.addInputToStatMergeJob(statMergeJob=boxPlotAssociationLocusAttributeVsOverlapThresholdJob, \
										parentJobLs=[associationLocusJob])
		
		
		sys.stderr.write("\t %s jobs.\n"%(self.no_of_jobs))
		return returnData
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2012.11.13
			add more
		2012.2.15
		"""
		DefineAssociationLandscapePipeline.registerCustomExecutables(self, workflow=workflow)
		
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		
		#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		executableClusterSizeMultiplierList = []
		
		BoxPlotAssociationLocusAttributeVsOverlapThreshold = Executable(namespace=namespace, name="BoxPlotAssociationLocusAttributeVsOverlapThreshold", \
						version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		BoxPlotAssociationLocusAttributeVsOverlapThreshold.addPFN(PFN("file://" + \
						os.path.join(self.variationSrcPath, "association_peak/plot/BoxPlotAssociationLocusAttributeVsOverlapThreshold.py"), site_handler))
		executableClusterSizeMultiplierList.append((BoxPlotAssociationLocusAttributeVsOverlapThreshold, 0.1))
		
		PlotAssociationLocusFrequencyOnGenome = Executable(namespace=namespace, name="PlotAssociationLocusFrequencyOnGenome", \
						version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		PlotAssociationLocusFrequencyOnGenome.addPFN(PFN("file://" + \
						os.path.join(self.variationSrcPath, "association_peak/plot/PlotAssociationLocusFrequencyOnGenome.py"), site_handler))
		executableClusterSizeMultiplierList.append((PlotAssociationLocusFrequencyOnGenome, 0.5))
		
		PlotAssociationLocusFrequencyVsAssociationThreshold = Executable(namespace=namespace, name="PlotAssociationLocusFrequencyVsAssociationThreshold", \
						version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		PlotAssociationLocusFrequencyVsAssociationThreshold.addPFN(PFN("file://" + \
						os.path.join(self.variationSrcPath, "association_peak/plot/PlotAssociationLocusFrequencyVsAssociationThreshold.py"), site_handler))
		executableClusterSizeMultiplierList.append((PlotAssociationLocusFrequencyVsAssociationThreshold, 0.2))
		
		TwoAssociationLocusFileOverlap = Executable(namespace=namespace, name="TwoAssociationLocusFileOverlap", \
						version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		TwoAssociationLocusFileOverlap.addPFN(PFN("file://" + \
						os.path.join(self.variationSrcPath, "association_peak/TwoAssociationLocusFileOverlap.py"), site_handler))
		executableClusterSizeMultiplierList.append((TwoAssociationLocusFileOverlap, 0.3))
		
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	
	def run(self):
		"""
		2011-10
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		result_query = self.db.getResultLs(analysis_method_id_ls=self.analysis_method_id_ls, \
						phenotype_method_id_ls=self.phenotype_method_id_ls, call_method_id_ls=self.call_method_id_ls,\
						cnv_method_id=self.cnv_method_id)
		
		workflow = self.initiateWorkflow()
		
		self.registerExecutables()
		self.registerCustomExecutables(workflow)
		association_result_ls = [ row for row in result_query]
		
		self.addAllJobs(workflow=workflow, association_result_ls=association_result_ls, data_dir=self.data_dir, min_MAF=self.min_MAF, \
				neighbor_distance=self.neighbor_distance, max_neighbor_distance=self.max_neighbor_distance, \
				min_score_ls=self.min_score_ls, min_overlap_ratio_ls=self.min_overlap_ratio_ls, ground_score=self.ground_score,\
				peakPadding=self.peakPadding, tax_id=self.tax_id, \
				outputDirPrefix=self.pegasusFolderName, transferOutput=True)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		self.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = BootstrapAssociationPeakWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
