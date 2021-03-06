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
		-o dags/BootstrapAssociationPeak/Bootstrap_call57_75_ana_1_32_min_score_3_4_5_min_overlap_0.05_0.2_0.5_0.8.xml
		-u yh  --db_passwd secret -z banyan -j condorpool -l condorpool

	#2013.1.11, on hoffman2
	%s  --call_method_id_ls 57,75 --analysis_method_id_ls 1,32 --min_score_ls 2,3,4,5 --min_overlap_ratio_ls 0.50,0.80
		-o dags/BootstrapAssociationPeak/Bootstrap_call57_75_ana_1_32_min_score_2_3_4_5_min_overlap_0.5_0.8.xml
		--phenotype_method_id_ls 1-35,39-48,57-82,158-159,161-179,182-186,272-283,314-351,362-380,418-589
		-z localhost -j hcondor -l hcondor
		--data_dir /u/home/eeskin/polyacti/NetworkData/250k/db/
		--db_user yh --db_passwd secret -C 50 --commit --needSSHDBTunnel
	
Description:
	2013.2.6 a workflow that checks association landscape/peak/loci under different threshold settings.
		it'll check with db to see if some landscape/peak/loci have been registered in db. if yes, it'll skip those jobs.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])


#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:	   #64bit
#	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
#	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
#else:   #32bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, utils

from variation.src.association_peak.DefineAssociationLandscapePipeline import DefineAssociationLandscapePipeline

class BootstrapAssociationPeakWorkflow(DefineAssociationLandscapePipeline):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(DefineAssociationLandscapePipeline.option_default_dict)
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
				parentJobLs=None, job_max_memory=100, walltime = 60, \
				extraDependentInputLs=None, \
				transferOutput=False, **keywords):
		"""
		2012.11.28
			walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
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
					walltime=walltime,\
					sshDBTunnel=False, \
					objectWithDBArguments=None, **keywords)
		return job
		
	
	def preReduce(self, workflow=None, outputDirPrefix="", passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		returnData = PassingData()
		returnData.jobDataLs = []
		
		mapDirJob = self.addMkDirJob(outputDir="%sMap"%(outputDirPrefix))
		passingData.mapDirJob = mapDirJob
		returnData.mapDirJob = mapDirJob
		
		countAssociationLocusOutputDirJob = self.addMkDirJob(outputDir="%sCountAssociationLocus"%(outputDirPrefix))
		passingData.countAssociationLocusOutputDirJob = countAssociationLocusOutputDirJob
		returnData.countAssociationLocusOutputDirJob = countAssociationLocusOutputDirJob
		
		plotOutputDirJob = self.addMkDirJob(outputDir='%sPlot'%(outputDirPrefix))
		passingData.plotOutputDirJob = plotOutputDirJob
		returnData.plotOutputDirJob = plotOutputDirJob
		
		reduceOutputDirJob = self.addMkDirJob(outputDir="%sReduce"%(outputDirPrefix))
		passingData.reduceOutputDirJob = reduceOutputDirJob
		returnData.reduceOutputDirJob = reduceOutputDirJob
		
		return returnData
	
	def addDefineAssociationLandscapeJobs(self, db_250k=None, association_result_ls=None, mapDirJob=None, \
										association_landscape_type=None,min_score_ls=None,\
										data_dir=None, tax_id=None, passingData=None):
		"""
		2013.2.7
		"""
		sys.stderr.write("\t Adding define-association-landscape jobs, #jobs=%s..."%(self.no_of_jobs))
		
		for result in association_result_ls:
			associationResultFile = self.registerOneInputFile(inputFname=result.getFileAbsPath(oldDataDir=db_250k.data_dir, \
											newDataDir=self.data_dir), \
										folderName=self.pegasusFolderName)
			
			association_landscape = db_250k.checkAssociationLandscape(result_id=result.id, \
														association_landscape_type_id=association_landscape_type.id)
			if association_landscape:	#2013.2.6 already in db
				landscapeOutputFname = association_landscape.getFileAbsPath(oldDataDir=db_250k.data_dir, \
											newDataDir=data_dir)
				landscapeOutputFile = self.registerOneInputFile(inputFname=landscapeOutputFname, \
										folderName=self.pegasusFolderName)
				defineLandscapeJob = PassingData(output=landscapeOutputFile)
				landscape2DBJob = None
			else:
				#associationResultFile = self.registerOneInputFile(inputFname=result.getFileAbsPath(oldDataDir=self.db_250k.data_dir, \
				#				newDataDir=self.data_dir), \
				#				folderName=self.pegasusFolderName)
				#add DefineAssociationLandscape job
				outputFnamePrefix = 'result_%s_neighbor_%s_max_neighbor_%s_min_MAF_%s_tax_id_%s'%\
						(result.id, association_landscape_type.neighbor_distance, \
						association_landscape_type.max_neighbor_distance, association_landscape_type.min_MAF, tax_id)
				landscapeOutputFile = File(os.path.join(mapDirJob.output, '%s_landscape.h5'%(outputFnamePrefix)))
				defineLandscapeJob = self.addDefineLandscapeJob(executable=self.DefineAssociationLandscape, \
								result_id=result.id, neighbor_distance=association_landscape_type.neighbor_distance, \
								max_neighbor_distance=association_landscape_type.max_neighbor_distance,\
								min_MAF=association_landscape_type.min_MAF, tax_id=tax_id, \
								data_dir=data_dir, logFile=None,\
								landscapeOutputFile=landscapeOutputFile, \
								extraDependentInputLs=[associationResultFile], \
								parentJobLs=[mapDirJob], sshDBTunnel=self.needSSHDBTunnel,\
								job_max_memory=2000, transferOutput=False)
						
				logFile = File(os.path.join(mapDirJob.output, '%s_2db_log.tsv'%\
											(outputFnamePrefix)))
				landscape2DBJob = self.addAssociationLandscape2DBJob(executable=self.AssociationLandscape2DB, inputFile=defineLandscapeJob.output, \
							result_id=result.id, \
							data_dir=self.data_dir, logFile=logFile, commit=self.commit, \
							min_MAF=association_landscape_type.min_MAF, \
							neighbor_distance=association_landscape_type.neighbor_distance, \
							max_neighbor_distance=association_landscape_type.max_neighbor_distance, \
							parentJobLs=[mapDirJob, defineLandscapeJob], \
							extraDependentInputLs=None, transferOutput=True, extraArguments=None, job_max_memory=1000, sshDBTunnel=self.needSSHDBTunnel)
			passingData.resultID2defineLandscapeJobData[result.id] = PassingData(defineLandscapeJob=defineLandscapeJob, landscape2DBJob=landscape2DBJob,\
																association_landscape=association_landscape)

			call_method_id = result.call_method_id
			analysis_method_id = result.analysis_method_id
			association_group_key = (call_method_id, analysis_method_id)
			if association_group_key not in passingData.association_group_key2resultList:
				passingData.association_group_key2resultList[association_group_key] = []
				passingData.association_group_key2reduceAssociationPeakJobMatrix[association_group_key] = []
				for min_score in min_score_ls:
					passingData.association_group_key2reduceAssociationPeakJobMatrix[association_group_key].append([])
				passingData.association_group_key2countAssociationLocusJobList[association_group_key] = []
			passingData.association_group_key2resultList[association_group_key].append(result)
			
		sys.stderr.write(" \t %s total jobs.\n"%(self.no_of_jobs))
		return passingData
	
	def addCheckTwoAssociationLocusFileOverlapJobs(self, analysis_method_id2AssociationLocusJobList=None,\
												min_score=None, min_overlap_ratio=None, associationMinScoreDirJob=None,\
												plotOutputDirJob=None, biologyCategoryID2PhenotypeIDSet=None):
		"""
		2013.2.7
		
		"""
		sys.stderr.write("\t\t Adding jobs handling two-association-locus-files #jobs=%s ..."%(self.no_of_jobs))
		#add TwoAssociationLocusFileOverlap & HistogramAssociationLocusAdjacencyDistance jobs
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
												job_max_memory=1000, walltime = 60, \
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
						figureDPI=200, formatString='.', ylim_type=2, samplingRate=1, need_svg=False, \
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
												job_max_memory=1000, walltime = 60, \
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
						figureDPI=200, formatString='.', ylim_type=2, samplingRate=1, need_svg=False, \
						logCount=True, inputFileFormat=2, \
						parentJobLs=[plotOutputDirJob, twoAssociationLocusOverlapJob], \
						extraDependentInputLs=None, \
						extraArgumentList=None, extraArguments=None, transferOutput=True, job_max_memory=1000)
					
					#2013.2.7 
					outputFnamePrefix = 'compare_association_call_%s_to_call_%s_ana_%s_loci_min_score_%s_min_overlap_%s'%\
									(associationLocusJob1.call_method_id, associationLocusJob2.call_method_id, \
									analysis_method_id, min_score, min_overlap_ratio)
					outputFile = File(os.path.join(associationMinScoreDirJob.output, '%s.h5'%(outputFnamePrefix)))
					
					compareGWAssociationLocusByPhenotypeVectorJob = self.addTwoAssociationLocusFileOverlapJob(\
												executable=self.CompareTwoGWAssociationLocusByPhenotypeVector, \
												inputFileList=[associationLocusJob1.output, associationLocusJob2.output], \
												inputFile=None, outputFile=outputFile, \
												parentJobLs=[associationLocusJob1, associationLocusJob2], \
												job_max_memory=1000, walltime = 60, \
												transferOutput=True)
					
					outputFile = File(os.path.join(plotOutputDirJob.output, '%s_2Dhist.png'%(outputFnamePrefix)))
					self.addDraw2DHistogramOfMatrixJob(inputFileList=None, inputFile=compareGWAssociationLocusByPhenotypeVectorJob.output, \
							outputFile=outputFile, \
							outputFnamePrefix=None, whichColumn=None, whichColumnHeader="fraction_of_input2_only_phenotypes", \
							whichColumnPlotLabel="fraction_of_call%s_ana%s_only_phenotypes"%(associationLocusJob2.call_method_id, \
																						associationLocusJob2.analysis_method_id), \
							logX=False, logY=False, logZ=False, valueForNonPositiveYValue=-1, \
							missingDataNotation='NA',\
							xColumnHeader="fraction_of_input1_only_phenotypes", \
							xColumnPlotLabel="fraction_of_call%s_ana%s_only_phenotypes"%(associationLocusJob1.call_method_id, \
																						associationLocusJob1.analysis_method_id), \
							minNoOfTotal=0,\
							figureDPI=100, formatString='.', samplingRate=1, need_svg=False, \
							inputFileFormat=2, outputFileFormat=None,\
							zColumnHeader=None, \
							parentJobLs=[plotOutputDirJob, compareGWAssociationLocusByPhenotypeVectorJob], \
							extraDependentInputLs=None, \
							extraArgumentList=None, extraArguments=None, transferOutput=True,  job_max_memory=2000)
					#2013.2.23 plot the peak frequency over genome distance for each peak-category (confirmatory, del-gwas-novel, etc.)
					# 0: all; 1: confirmatory; 2: del-gwas-novel; 3:SNP-gwas-novel; 4:specific-in-both; 5: other
					#. also split for each phenotype category 
					for associationLocusCategory in xrange(6):
						#choose  peaks that belong to specific category
						outputFile = File(os.path.join(associationMinScoreDirJob.output, \
											'%s_associationLocusCategory%s.h5'%(outputFnamePrefix, associationLocusCategory)))
						filterTwoGWAssociationLocusComparisonResultJob = self.addGenericJob(executable=self.FilterTwoGWAssociationLocusComparisonResult, \
									inputFile=compareGWAssociationLocusByPhenotypeVectorJob.output, \
									outputFile=outputFile, inputFileList=[], \
									parentJobLs=[compareGWAssociationLocusByPhenotypeVectorJob, associationMinScoreDirJob], \
									extraDependentInputLs=[], extraOutputLs=[], transferOutput=False, \
									extraArguments=None, extraArgumentList=["--associationLocusCategory %s"%(associationLocusCategory)], \
									job_max_memory=1000,  sshDBTunnel=None, \
									key2ObjectForJob=None, no_of_cpus=None, walltime=120)
						#add a plot job
						#add PlotAssociationLocusFrequencyOnGenome job
						plotFilenamePrefix = os.path.join(plotOutputDirJob.output, \
													'frequency_manhattan_%s_associationLocusCategory%s'%(outputFnamePrefix, associationLocusCategory))
						outputFile = File('%s.png'%(plotFilenamePrefix))
						if self.drivername=='mysql':
							genome_dbname = 'genome'
						else:
							genome_dbname = self.dbname
						plotAssociationLocusFrequencyOnGenomeJob = self.addAbstractPlotJob(executable=self.PlotGenomeWideData, \
							inputFileList=None, inputFile=filterTwoGWAssociationLocusComparisonResultJob.output, outputFile=outputFile, \
							outputFnamePrefix=None, whichColumn=None, whichColumnHeader="no_of_total_phenotypes", \
							whichColumnPlotLabel="numberOfPhenotypes", \
							logX=None, logY=None, valueForNonPositiveYValue=-1, \
							xScaleLog=None, yScaleLog=None,\
							missingDataNotation='NA',\
							xColumnHeader="start", xColumnPlotLabel="genomePosition", \
							minNoOfTotal=1, maxNoOfTotal=None,\
							figureDPI=200, formatString='.', ylim_type=2, samplingRate=1, need_svg=False, \
							inputFileFormat=None, outputFileFormat=None,\
							parentJobLs=[plotOutputDirJob, filterTwoGWAssociationLocusComparisonResultJob], \
							extraDependentInputLs=None, \
							extraArgumentList=['--genome_drivername=%s'%self.drivername,\
								'--genome_hostname=%s'%self.hostname,\
								'--genome_dbname=%s'%(genome_dbname),\
								'--genome_schema=genome',\
								'--genome_db_user=%s'%(self.db_user),\
								'--genome_db_passwd=%s'%(self.db_passwd),\
								'--tax_id=%s'%(self.tax_id), '--drawCentromere', '--chromosomeHeader chromosome'], \
							extraArguments=None, transferOutput=True,  job_max_memory=2000, \
							sshDBTunnel=self.needSSHDBTunnel, \
							objectWithDBArguments=self)
						for biologyCategoryID, phenotypeIDSet in biologyCategoryID2PhenotypeIDSet.iteritems():
							phenotypeIDListInStr = utils.getSuccinctStrOutOfList(phenotypeIDSet)
							#. add a filter peak job
							outputFile = File(os.path.join(associationMinScoreDirJob.output, \
									'%s_associationLocusCategory%s_biologyCategory%s.h5'%(outputFnamePrefix, associationLocusCategory,\
																						biologyCategoryID)))
							filterTwoGWAssociationLocusComparisonResultJob = self.addGenericJob(executable=self.FilterTwoGWAssociationLocusComparisonResult, \
									inputFile=compareGWAssociationLocusByPhenotypeVectorJob.output, \
									outputFile=outputFile, inputFileList=[], \
									parentJobLs=[compareGWAssociationLocusByPhenotypeVectorJob, associationMinScoreDirJob], \
									extraDependentInputLs=[], extraOutputLs=[], transferOutput=False, \
									extraArguments=None, extraArgumentList=["--associationLocusCategory %s --phenotypeIDList %s"%\
																		(associationLocusCategory, phenotypeIDListInStr)], \
									job_max_memory=1000,  sshDBTunnel=None, \
									key2ObjectForJob=None, no_of_cpus=None, walltime=120)
							#. add a plot job
							#add PlotAssociationLocusFrequencyOnGenome job
							plotFilenamePrefix = os.path.join(plotOutputDirJob.output, \
											'frequency_manhattan_%s_associationLocusCategory%s_biologyCategory%s'%\
											(outputFnamePrefix, associationLocusCategory, biologyCategoryID))
							outputFile = File('%s.png'%(plotFilenamePrefix))
							if self.drivername=='mysql':
								genome_dbname = 'genome'
							else:
								genome_dbname = self.dbname
							plotAssociationLocusFrequencyOnGenomeJob = self.addAbstractPlotJob(executable=self.PlotGenomeWideData, \
								inputFileList=None, inputFile=filterTwoGWAssociationLocusComparisonResultJob.output, outputFile=outputFile, \
								outputFnamePrefix=None, whichColumn=None, whichColumnHeader="no_of_total_phenotypes", \
								whichColumnPlotLabel="numberOfPhenotypes", \
								logX=None, logY=None, valueForNonPositiveYValue=-1, \
								xScaleLog=None, yScaleLog=None,\
								missingDataNotation='NA',\
								xColumnHeader="start", xColumnPlotLabel="genomePosition", \
								minNoOfTotal=1, maxNoOfTotal=None,\
								figureDPI=200, formatString='.', ylim_type=2, samplingRate=1, need_svg=False, \
								inputFileFormat=None, outputFileFormat=None,\
								parentJobLs=[plotOutputDirJob, filterTwoGWAssociationLocusComparisonResultJob], \
								extraDependentInputLs=None, \
								extraArgumentList=['--genome_drivername=%s'%self.drivername,\
									'--genome_hostname=%s'%self.hostname,\
									'--genome_dbname=%s'%(genome_dbname),\
									'--genome_schema=genome',\
									'--genome_db_user=%s'%(self.db_user),\
									'--genome_db_passwd=%s'%(self.db_passwd),\
									'--tax_id=%s'%(self.tax_id), '--drawCentromere', '--chromosomeHeader chromosome'], \
								extraArguments=None, transferOutput=True,  job_max_memory=2000, \
								sshDBTunnel=self.needSSHDBTunnel, \
								objectWithDBArguments=self)
		sys.stderr.write("%s total jobs.\n"%(self.no_of_jobs))
		
	def addAssociationPeakAndLocusJobs(self, db_250k=None, min_score_ls=None, association_landscape_type=None,\
									data_dir=None, ground_score=None, min_overlap_ratio_ls=None, 
									phenotype_method_id_of_associations_set=None, \
									biologyCategoryID2PhenotypeIDSet=None,\
									plotOutputDirJob=None, passingData=None,\
									**keywords):
		"""
		#2013.2.7 add the association-peak, association-locus jobs
		"""
		sys.stderr.write("\t Adding association-peak & association-locus jobs, #jobs=%s ... \n"%(self.no_of_jobs))
		for i in xrange(len(min_score_ls)):
			min_score = min_score_ls[i]
			#create a folder here.
			associationMinScoreDir = "association_min_score_%s"%(min_score)
			associationMinScoreDirJob = self.addMkDirJob(outputDir=associationMinScoreDir)
			
			resultID2associationPeakJob = {}
			
			association_peak_type = self.db_250k.getAssociationPeakType(association_landscape_type_id=association_landscape_type.id, \
															min_score=min_score)
			for resultID, defineLandscapeJobData in passingData.resultID2defineLandscapeJobData.iteritems():
				defineLandscapeJob = defineLandscapeJobData.defineLandscapeJob
				landscape2DBJob = defineLandscapeJobData.landscape2DBJob
				association_landscape = defineLandscapeJobData.association_landscape
				genome_wide_association_peak = self.db_250k.checkGenomeWideAssociationPeak(result_id=resultID, \
													association_landscape_id=association_landscape.id,\
													association_peak_type_id=association_peak_type.id)
				if genome_wide_association_peak:
					associationPeakFname = genome_wide_association_peak.getFileAbsPath(oldDataDir=db_250k.data_dir, \
												newDataDir=data_dir)
					associationPeakFile = self.registerOneInputFile(inputFname=associationPeakFname, \
											folderName=self.pegasusFolderName)
					associationPeakJob = PassingData(output=associationPeakFile)
					associationPeak2DBJob = None
				else:
					#add landscape -> peak job
					outputFnamePrefix = os.path.splitext(os.path.basename(defineLandscapeJob.output.name))[0]
					outputFnamePrefix = '%s_min_score_%s_ground_score_%s'%(outputFnamePrefix, min_score, ground_score)
					outputFile = File(os.path.join(associationMinScoreDirJob.output, '%s_peak.h5'%(outputFnamePrefix)))
					
					associationPeakJob = self.addAssociationLandscape2PeakJob(executable=self.AssociationLandscape2Peak, \
						inputFile=defineLandscapeJob.output, outputFile=outputFile, \
						min_score=min_score, ground_score=ground_score, \
						parentJobLs=[associationMinScoreDirJob, defineLandscapeJob], job_max_memory=2000, walltime = 60, \
						extraDependentInputLs=None, \
						transferOutput=False)
					
					if landscape2DBJob or association_landscape:	#it needs the landscape already in db_250k
						logFile = File(os.path.join(associationMinScoreDirJob.output, '%s_2db_log.tsv'%\
												(outputFnamePrefix)))
						associationPeak2DBJob = self.addAssociationLandscape2DBJob(executable=self.AssociationPeak2DB, inputFile=associationPeakJob.output, \
								result_id=resultID, \
								data_dir=self.data_dir, logFile=logFile, commit=self.commit, \
								min_MAF=association_landscape_type.min_MAF, \
								neighbor_distance=association_landscape_type.neighbor_distance, \
								max_neighbor_distance=association_landscape_type.max_neighbor_distance, \
								parentJobLs=[associationMinScoreDirJob, associationPeakJob, landscape2DBJob], \
								extraDependentInputLs=None, transferOutput=True, extraArguments="--min_score %s"%(min_score), \
								job_max_memory=1000, sshDBTunnel=self.needSSHDBTunnel)
				resultID2associationPeakJob[resultID] = associationPeakJob
			
			for j in xrange(len(min_overlap_ratio_ls)):
				min_overlap_ratio = min_overlap_ratio_ls[j]
				#for HistogramAssociationLocusAdjacencyDistance,
				#comparing association locus from two different call methods, but same analysis method
				analysis_method_id2AssociationLocusJobList = {}
				
				for association_group_key, result_ls in passingData.association_group_key2resultList.iteritems():
					call_method_id, analysis_method_id = association_group_key[:2]
					associationLocusFnamePrefix = 'call_%s_analysis_%s_min_score_%s_min_overlap_%s'%\
							(call_method_id, analysis_method_id, min_score, min_overlap_ratio)
					
					association_locus_type = db_250k.getAssociationLocusType(association_peak_type_id=association_peak_type.id, \
									min_overlap_ratio=min_overlap_ratio, \
									min_connectivity=None)
					genome_wide_association_locus = db_250k.checkGenomeWideAssociationLocus(association_locus_type_id=association_locus_type.id,\
								call_method_id=call_method_id, analysis_method_id=analysis_method_id,\
								call_method_id_ls='%s'%(call_method_id),\
								analysis_method_id_ls="%s"%(analysis_method_id),\
								phenotype_method_id_ls=utils.getSuccinctStrOutOfList(phenotype_method_id_of_associations_set))
					if genome_wide_association_locus:
						associationLocusFname = genome_wide_association_locus.getFileAbsPath(oldDataDir=db_250k.data_dir, \
													newDataDir=data_dir)
						associationLocusFile = self.registerOneInputFile(inputFname=associationLocusFname, \
												folderName=self.pegasusFolderName)
						associationLocusJob = PassingData(output=associationLocusFile,call_method_id=call_method_id, \
														analysis_method_id=analysis_method_id)
						associationLocus2DBJob = None
					else:
						#add a AssociationPeak2AssociationLocus job
						associationLocusFile = File(os.path.join(associationMinScoreDirJob.output, '%s_locus.h5'%(associationLocusFnamePrefix)))
						#input to this job is added later
						associationLocusJob = self.addAssociationPeak2LocusJob(executable=self.AssociationPeak2AssociationLocus, \
									inputFile=None, outputFile=associationLocusFile, min_overlap_ratio=min_overlap_ratio, \
									parentJobLs=[associationMinScoreDirJob], job_max_memory=2000, walltime = 60, \
									extraDependentInputLs=None, \
									transferOutput=False)
						
						associationLocusJob.call_method_id = call_method_id
						associationLocusJob.analysis_method_id = analysis_method_id
						
						
						#add to the database
						logFile = File(os.path.join(associationMinScoreDirJob.output, '%s_2db_log.tsv'%\
												(associationLocusFnamePrefix)))
						associationLocus2DBJob = self.addAssociationLandscape2DBJob(executable=self.AssociationLocus2DB, \
								inputFile=associationLocusJob.output, \
								result_id=None, \
								data_dir=self.data_dir, logFile=logFile, commit=self.commit, \
								min_MAF=association_landscape_type.min_MAF, \
								neighbor_distance=association_landscape_type.neighbor_distance, \
								max_neighbor_distance=association_landscape_type.max_neighbor_distance, \
								parentJobLs=[associationMinScoreDirJob, associationLocusJob], \
								extraDependentInputLs=None, transferOutput=True, \
								extraArguments="--min_score %s --min_overlap_ratio %s"%(min_score, min_overlap_ratio), \
								job_max_memory=1000, sshDBTunnel=self.needSSHDBTunnel)
					
					passingData.association_group_key2reduceAssociationPeakJobMatrix[association_group_key][i].append(associationLocusJob)
					
					#add PlotAssociationLocusFrequencyOnGenome job
					associationLocusFrequencyOnGenomeFnamePrefix = os.path.join(plotOutputDirJob.output, \
												'frequency_manhattan_%s'%(associationLocusFnamePrefix))
					outputFile = File('%s.png'%(associationLocusFrequencyOnGenomeFnamePrefix))
					if self.drivername=='mysql':
						genome_dbname = 'genome'
					else:
						genome_dbname = self.dbname
					plotAssociationLocusFrequencyOnGenomeJob = self.addAbstractPlotJob(executable=self.PlotAssociationLocusFrequencyOnGenome, \
						inputFileList=None, inputFile=associationLocusFile, outputFile=outputFile, \
						outputFnamePrefix=None, whichColumn=None, whichColumnHeader="no_of_results", \
						whichColumnPlotLabel="numberOfResults", \
						logX=None, logY=None, valueForNonPositiveYValue=-1, \
						xScaleLog=None, yScaleLog=1,\
						missingDataNotation='NA',\
						xColumnHeader="start", xColumnPlotLabel="genomePosition", \
						minNoOfTotal=1, maxNoOfTotal=None,\
						figureDPI=200, formatString='.', ylim_type=2, samplingRate=1, need_svg=False, \
						inputFileFormat=None, outputFileFormat=None,\
						parentJobLs=[plotOutputDirJob, associationLocusJob], \
						extraDependentInputLs=None, \
						extraArgumentList=['--genome_drivername=%s'%self.drivername,\
							'--genome_hostname=%s'%self.hostname,\
							'--genome_dbname=%s'%(genome_dbname),\
							'--genome_schema=genome',\
							'--genome_db_user=%s'%(self.db_user),\
							'--genome_db_passwd=%s'%(self.db_passwd),\
							'--tax_id=%s'%(self.tax_id), '--drawCentromere'], extraArguments=None, transferOutput=True,  job_max_memory=2000, \
						sshDBTunnel=self.needSSHDBTunnel, \
						objectWithDBArguments=self)
					
					#self.addInputToStatMergeJob(statMergeJob=plotAssociationLocusFrequencyOnGenomeJob, parentJobLs=[associationLocusJob])
					
					if analysis_method_id not in analysis_method_id2AssociationLocusJobList:
						analysis_method_id2AssociationLocusJobList[analysis_method_id] = []
					analysis_method_id2AssociationLocusJobList[analysis_method_id].append(associationLocusJob)
					
					if isinstance(associationLocusJob, Job):	#2013.2.6 it could PassingData
						for result in result_ls:
							associationPeakJob = resultID2associationPeakJob.get(result.id)
							self.addInputToStatMergeJob(statMergeJob=associationLocusJob, parentJobLs=[associationPeakJob])
					self.addCheckTwoAssociationLocusFileOverlapJobs(analysis_method_id2AssociationLocusJobList=analysis_method_id2AssociationLocusJobList, \
													min_score=min_score, min_overlap_ratio=min_overlap_ratio, \
													associationMinScoreDirJob=associationMinScoreDirJob, plotOutputDirJob=plotOutputDirJob,\
													biologyCategoryID2PhenotypeIDSet=biologyCategoryID2PhenotypeIDSet)
		sys.stderr.write("\t %s total jobs.\n"%(self.no_of_jobs))
		
		
	
	def addAllJobs(self, workflow=None, db_250k=None, association_result_ls=None, \
				data_dir=None, min_MAF=None, \
				neighbor_distance=None, max_neighbor_distance=None, \
				min_score_ls=None, min_overlap_ratio_ls=None, ground_score=None,\
				peakPadding=None, tax_id=None, \
				outputDirPrefix="", transferOutput=True, job_max_memory=2000, **keywords):
		"""
		2012.11.21
			
		"""
		if workflow is None:
			workflow = self
		#2013.2.24 setup some data strutures here
		phenotype_method_id_of_associations_set = set()
		biologyCategoryID2PhenotypeIDSet = {}
		for association_result in association_result_ls:
			biologyCategoryID = association_result.phenotype_method.biology_category_id
			if biologyCategoryID not in biologyCategoryID2PhenotypeIDSet:
				biologyCategoryID2PhenotypeIDSet[biologyCategoryID] = set()
			biologyCategoryID2PhenotypeIDSet[biologyCategoryID].add(association_result.phenotype_method_id)
			phenotype_method_id_of_associations_set.add(association_result.phenotype_method_id)
		
		sys.stderr.write("Adding jobs for %s association results (%s phenotypes, %s biology categories) #jobs=%s... \n"%\
							(len(association_result_ls),\
							len(phenotype_method_id_of_associations_set), \
							len(biologyCategoryID2PhenotypeIDSet), self.no_of_jobs))
		
		returnData = PassingData()
		returnData.jobDataLs = []
		
		passingData = PassingData(fnamePrefix=None, \
					outputDirPrefix=outputDirPrefix, \
					jobData=None,\
					preReduceReturnData=None,\
					association_group_key2orderIndex = {},\
					association_group_key2resultList = {},\
					association_group_key2reduceAssociationPeakJobMatrix = {},\
					association_group_key2countAssociationLocusJobList = {},\
					resultID2defineLandscapeJobData = {},
					)
		
		preReduceReturnData = self.preReduce(workflow=workflow, outputDirPrefix=outputDirPrefix, \
									passingData=passingData, transferOutput=False,\
									**keywords)
		
		mapDirJob = preReduceReturnData.mapDirJob
		plotOutputDirJob = preReduceReturnData.plotOutputDirJob
		countAssociationLocusOutputDirJob = preReduceReturnData.countAssociationLocusOutputDirJob
		reduceOutputDirJob = preReduceReturnData.reduceOutputDirJob
		
		passingData.preReduceReturnData = preReduceReturnData
		
		
		
		#the job matrix is a matrix of AssociationPeak2AssociationLocus jobs.
		#	each row is index by association threshold.
		#	each column is index by overlap threshold.
		
		association_landscape_type = db_250k.getAssociationLandscapeType(min_MAF=min_MAF, \
											neighbor_distance=neighbor_distance, \
											max_neighbor_distance=max_neighbor_distance)
		#2013.2.7
		self.addDefineAssociationLandscapeJobs(db_250k=db_250k, association_result_ls=association_result_ls, \
											mapDirJob=mapDirJob, association_landscape_type=association_landscape_type, \
											min_score_ls=min_score_ls, data_dir=data_dir, tax_id=tax_id, passingData=passingData)
		
		self.addAssociationPeakAndLocusJobs(db_250k=db_250k, min_score_ls=min_score_ls, association_landscape_type=association_landscape_type,\
										data_dir=data_dir, ground_score=ground_score, min_overlap_ratio_ls=min_overlap_ratio_ls, \
										phenotype_method_id_of_associations_set=phenotype_method_id_of_associations_set, \
										biologyCategoryID2PhenotypeIDSet = biologyCategoryID2PhenotypeIDSet,\
										plotOutputDirJob=plotOutputDirJob, passingData=passingData)
		
		sys.stderr.write("Adding PlotAssociationLocusFrequencyVsAssociationThreshold jobs, #jobs=%s ..."%(self.no_of_jobs))
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
					figureDPI=200, formatString='.-', ylim_type=2, samplingRate=1, need_svg=False, \
					inputFileFormat=2, outputFileFormat=None,\
					parentJobLs=[plotOutputDirJob], \
					extraDependentInputLs=None, \
					extraArgumentList=None, extraArguments=None, transferOutput=True,  job_max_memory=2000, \
					sshDBTunnel=None, \
					objectWithDBArguments=None)
			
			for association_group_key, countAssociationLocusJobList in passingData.association_group_key2countAssociationLocusJobList.iteritems():
				#add a CountAssociationLocus job across different association thresholds (min_score)
				countAssociationLocusFnamePrefix = os.path.join(countAssociationLocusOutputDirJob.output, 'call_%s_analysis_%s_min_overlap_%s'%\
										(association_group_key[0], association_group_key[1], min_overlap_ratio))
				countAssociationLocusFile = File('%s_loci_count.h5'%(countAssociationLocusFnamePrefix))
				#input to this job is added later
				countAssociationLocusJob = self.addCountAssociationLocusJob(executable=self.CountAssociationLocus, \
								inputFileList=None, inputFile=None, outputFile=countAssociationLocusFile, \
								parentJobLs=[countAssociationLocusOutputDirJob], job_max_memory=1000, walltime = 6, \
								extraDependentInputLs=None, \
								transferOutput=False)
				countAssociationLocusJobList.append(countAssociationLocusJob)
				for i in xrange(len(min_score_ls)):
					min_score = min_score_ls[i]
					associationLocusJob = passingData.association_group_key2reduceAssociationPeakJobMatrix[association_group_key][i][j]
					self.addInputToStatMergeJob(statMergeJob=countAssociationLocusJob, parentJobLs=[associationLocusJob])
				
				self.addInputToStatMergeJob(statMergeJob=plotAssociationLocusFrequencyVsAssociationThresholdJob, \
										parentJobLs=[countAssociationLocusJob])
		sys.stderr.write("\t %s total jobs so far.\n"%(self.no_of_jobs))
		
		sys.stderr.write("Adding BoxPlotAssociationLocusAttributeVsOverlapThreshold jobs,#jobs=%s ..."%(self.no_of_jobs))
		#add the BoxPlotAssociationLocusAttributeVsOverlapThreshold jobs
		for association_group_key, jobMatrix in passingData.association_group_key2reduceAssociationPeakJobMatrix.iteritems():
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
					figureDPI=200, formatString='.', ylim_type=2, samplingRate=1, need_svg=False, \
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
		
		
		sys.stderr.write("\t %s total jobs.\n"%(self.no_of_jobs))
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
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.variationSrcPath, \
													"association_peak/TwoAssociationLocusFileOverlap.py"), \
													name="TwoAssociationLocusFileOverlap", clusterSizeMultipler=0.3)
		#2013.2.7	new simpler way of adding an executable
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.variationSrcPath, \
													"association_peak/CompareTwoGWAssociationLocusByPhenotypeVector.py"), \
													name="CompareTwoGWAssociationLocusByPhenotypeVector", clusterSizeMultipler=0)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.variationSrcPath, \
													"association_peak/filters/FilterTwoGWAssociationLocusComparisonResult.py"), \
													name="FilterTwoGWAssociationLocusComparisonResult", clusterSizeMultipler=0)
	
	
	def run(self):
		"""
		2011-10
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		result_query = self.db_250k.getResultLs(analysis_method_id_ls=self.analysis_method_id_ls, \
						phenotype_method_id_ls=self.phenotype_method_id_ls, call_method_id_ls=self.call_method_id_ls,\
						cnv_method_id=self.cnv_method_id)
		
		workflow = self.initiateWorkflow()
		
		self.registerExecutables()
		self.registerCustomExecutables(workflow)
		association_result_ls = [ row for row in result_query]
		
		self.addAllJobs(workflow=workflow, db_250k=self.db_250k, association_result_ls=association_result_ls, \
				data_dir=self.data_dir, min_MAF=self.min_MAF, \
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
