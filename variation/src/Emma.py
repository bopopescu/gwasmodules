#!/usr/bin/env python2.5
"""
Usage: Emma.py [OPTIONS] -o R_FILE SNPS_DATA_FILE PHENOTYPE_DATA_FILE PHENOTYPE_INDEX

Option:

	-o ..., --rFile=...			Emma output.
	-d ..., --delim=...			default is ", "	  
	-m ..., --missingval=...		default is "NA"
	-c ..., --chr=...			default is all chromosomes
	--logTransform				Log transforms the phenotype values, before running Emma.
	--removeOutliers=...			Remove outliers using removeOutliers*IQR as the fence.
	--addConstant=...			Adds the given value to the phenotypic values (before log-transforming) (=0 for sd/10)
	--negate				Negates the values
	--kinshipDatafile=			Datafile which is used to calculate the kinship matrix.		(default is the same as the given genotype file)
	--BoundaryStart=...			Only the region within the boundary is considered in the GWA. (Default is no boundaries)
	--BoundaryEnd=...		   
	--phenotypeRanks			Use ranked phenotype values as the phenotype values.  (Invariant of transformation)
	--minMAF=...				Remove all SNPs which have MAF smaller than the given argument.  (0.0 is set as default).
	--LRT					Use ML and LRT instead of REML and a t-test.
	--phenotypeFileType=...			1 (default) if file has tsv format, 2 if file has csv format and contains accession names (instead of ecotype ID)
	-h, --help				show this help
	--parallel=...				Run Emma on the cluster with standard parameters.  The arguement is used for runid as well as output files.  
	--parallelAll				Run Emma on all phenotypes.
	--onlyMissing				Only missing runs (for parallelAll option only).
	--onlyOriginal96			Run Emma on the original 96 accessions only.
	--onlyOriginal192			Run Emma on the original 192 2010 data accessions only.
	--onlyBelowLatidue=...			Only accessions below the given latitude are left in.
	--complement				Choose the complement set of accessions (works only with onlyOriginal96 or onlyOriginal192).
	--subSample=...				Pick a random subsample of the given size.
	--subsampleTest=num_acc,num_samples	Outputs num_samples files, each with random num_acc accessions.
	
	--testRobustness			Perform a robustness test

	--permutationFilter=...                 Only use a random sample of proportional size equal to the given number.


	--sr							Perform a second run.
	--srSkipFirstRun				Skip first run
	--srInput=pvalFile 		       	Use given results as input. (use with srSkipFirstRun)
	--srOutput=resultFile 		    SOutput new results in given file. 
	--srPar=topQuantile,windowSize	        Default topQuantile is 0.95, and windowSize = 30000 bases. 
	
	
	--useLinearRegress			Use linear regression instead of KW.
	--regressionCofactors=...		Include a pickled dictionary containing factors.
	--FriLerAsCofactor
	--FriColAsCofactor

	--memReq=...				Request memory (on cluster), 5g is default
	--walltimeReq=...			Request time limit (on cluster), 100:00:00 (100 hours) is default

Examples:
	Emma.py -o emma_result.r  250K.csv phenotypes.tsv phenotype_index 
	
Description:

"""

import sys, getopt, traceback
import os, env, phenotypeData, plotResults, tempfile
import gwaResults, SecondStageAnalysis, dataParsers, snpsdata 
import math, time, numpy, util

#Not neccessary if using script based emma.
from numpy import *
from rpy import r


def calcKinship(snps):
	"""
	Requires EMMA to be installed.
	"""
	a = array(snps)
	#r.library("emma")
	print r.source("emma.R")
	print a
	#print r.emma_kinship(a)
	return r.emma_kinship(a)


def runEmma(phed,p_i,k,snps):
	#Assume that the accessions are ordered.
	i = phed.getPhenIndex(p_i)
	#r.library("emma")
	r.source("emma.R")
	phenValues = []
	for vals in phed.phenotypeValues:
		phenValues.append(float(vals[i]))
	phenArray = array([phenValues])
	snpsArray = array(snps)
	res = r.emma_REML_t(phenArray,snpsArray,k)
	#print res
	return res

def run_emma_w_missing_data(snpsd,phed,p_i,chromosome,k,missing_val='NA'):
	"""
	Runs Emma w. missing data.
	"""
	#Assume that the accessions are ordered.
	snps = snpsd.snps
	r.source("emma.R")
	phen_vals = phed.getPhenVals(p_i,noNAs=False)
	if phen_vals.count('NA'):
		print "Coordinating SNPs and phenotypes, removing %d accessions which were missing phenotypic values." % phen_vals.count('NA')
		ai_to_keep = [i for i in range(len(phen_vals)) if phen_vals[i]!='NA']
		#print ai_to_keep
		new_snps = []
		for snp in snps:
			new_snps.append([snp[i] for i in ai_to_keep])
		snps = new_snps
		#print snps[0]
		phen_vals = [phen_vals[i] for i in ai_to_keep]
		#print phen_vals
		new_k = zeros((len(ai_to_keep),len(ai_to_keep)))
		for j1 in range(len(ai_to_keep)):
			for j2 in range(len(ai_to_keep)):
				new_k[j1,j2]=k[ai_to_keep[j1],ai_to_keep[j2]]
		k = new_k
		 
	p_vals = []
	for snp in snps:
		
		if snp.count(missing_val):
			acc_to_keep = []
			new_snp = []
			for j, nt in enumerate(snp):
				if not nt == missing_val:
					acc_to_keep.append(j)
					new_snp.append(nt)
			new_k = zeros((len(acc_to_keep),len(acc_to_keep)))
			new_phen_vals = []
			for j1 in range(len(acc_to_keep)):
				for j2 in range(len(acc_to_keep)):
					new_k[j1,j2]=k[acc_to_keep[j1],acc_to_keep[j2]]
				new_phen_vals.append(phen_vals[acc_to_keep[j1]])
			
			if len(set(snp))==1:
				p_vals.append(1)
				print phenArray
				print snpsArray
				print new_k 
			else:
				phenArray = array([new_phen_vals])
				snpsArray = array([new_snp])
				#print "Running EMMA w. missing data."
				try:
					p_vals.extend(r.emma_REML_t(phenArray,snpsArray,new_k)['ps'])
				except Exception, err_str:
					print err_str
					print phenArray
					print snpsArray
					print new_k 
			  
		else:
			phenArray = array([phen_vals])
			snpsArray = array([snp])
			#print "Running EMMA w.o. missing data."
			try:
				p_vals.extend(r.emma_REML_t(phenArray,snpsArray,k)['ps'])
			except Exception, err_str:
				print err_str
				print phenArray
				print snpsArray
				print k 
	gwas_result = gwaResults.Result(snpsds=[snpsd],name="Emma_"+phed.getPhenotypeName(p_i),phenotypeID=p_i,scores=p_vals,chromosomes=[chromosome]) 
	return gwas_result


def _runEmma_(snps,phenValues,k):
	phenArray = array([phenValues])
	snpsArray = array(snps)
	res = r.emma_REML_t(phenArray,snpsArray,k)
	#print res
	return res


def _sampleSNPs_(snps,n,withReplacement=False):
	l = len(snps)
	i = 0
	snpSample = []
	if withReplacement:
		while i<n:
			r = random.randint(0,l-1)
			snpSample.append(snps[r])
			i += 1
	else:
		while i<n:
			r = random.randint(0,len(snps)-1)
			snpSample.append(snps[r])
			snps.remove(snps[r])
			i += 1		
	return snpSample
	


#For cluster use only!
#emmadir = '/home/cmb-01/bvilhjal/Projects/Python-snps/'
emmadir = '/home/cmbpanfs-01/bvilhjal/src/'
#resultDir = '/home/cmb-01/bvilhjal/results/'
resultDir = '/home/cmbpanfs-01/bvilhjal/results/'




def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["rFile=","chr=", "delim=", "missingval=", "BoundaryStart=", "removeOutliers=", "addConstant=",
			"logTransform", "BoundaryEnd=", "phenotypeFileType=", "help", "parallel=", "parallelAll", "LRT", "minMAF=", 
			"kinshipDatafile=", "phenotypeRanks", "onlyMissing","onlyOriginal96", "onlyOriginal192", "onlyBelowLatidue=", 
			"complement", "negate", "srInput=", "sr","srOutput=", "srPar=","srSkipFirstRun", "testRobustness",
			"permutationFilter=","useLinearRegress", "regressionCofactors=", "FriLerAsCofactor","FriColAsCofactor",
			"memReq=","walltimeReq=",]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:c:d:m:h", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	phenotypeRanks = False
	removeOutliers = None
	addConstant = -1
	phenotypeFileType = 1
	rFile = None
	delim = ","
	missingVal = "NA"
	help = 0
	minMAF=0.0
	boundaries = [-1,-1]
	chr=None
	parallel = None
	logTransform = False
	negate = False
	parallelAll = False
	lrt = False
	kinshipDatafile = None 
	onlyMissing = False
	onlyOriginal96 = False
	onlyOriginal192 = False
	onlyBelowLatidue = None
	complement = False

	sr = False
	srOutput = False
	srInput = False
	srSkipFirstRun = False
	srTopQuantile = 0.95
	srWindowSize = 30000
	
	testRobustness = False
	permutationFilter = 0.002

	useLinearRegress=False
	regressionCofactors=None
	FriLerAsCofactor = False
	FriColAsCofactor = False

	memReq = "5g"
	walltimeReq = "150:00:00"

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
		elif opt in ("-o","--rFile"):
			rFile = arg
		elif opt in ("--phenotypeFileType"):
			phenotypeFileType = int(arg)
		elif opt in ("--BoundaryStart"):
			boundaries[0] = int(arg)
		elif opt in ("--BoundaryEnd"):
			boundaries[1] = int(arg)
		elif opt in ("--addConstant"):
			addConstant = float(arg)
		elif opt in ("--parallel"):
			parallel = arg
		elif opt in ("--minMAF"):
			minMAF = float(arg)
		elif opt in ("--parallelAll"):
			parallelAll = True
		elif opt in ("--onlyMissing"):
			onlyMissing = True
		elif opt in ("--onlyOriginal96"):
			onlyOriginal96 = True
		elif opt in ("--onlyOriginal192"):
			onlyOriginal192 = True
		elif opt in ("--onlyBelowLatidue"):
			onlyBelowLatidue = float(arg)
		elif opt in ("--complement"):
			complement = True
		elif opt in ("--logTransform"):
			logTransform = True
		elif opt in ("--negate"):
			negate = True
		elif opt in ("--removeOutliers"):
			removeOutliers = float(arg)
		elif opt in ("--LRT"):
			lrt = True
		elif opt in ("-c","--chr"):
			chr = int(arg)
		elif opt in ("-d","--delim"):
			delim = arg
		elif opt in ("-m","--missingval"):
			missingVal = arg
		elif opt in ("--kinshipDatafile"):
			kinshipDatafile = arg
		elif opt in ("--phenotypeRanks"):
			phenotypeRanks = True
		elif opt in ("--sr"):
			sr = True
		elif opt in ("--srSkipFirstRun"):
			srSkipFirstRun = True
		elif opt in ("--srInput"):
			srInput = arg
		elif opt in ("--srOutput"):
			srOutput = arg
		elif opt in ("--srPar"):
			vals = arg.split(",")
			srTopQuantile = float(vals[0]) 
			srWindowSize = int(vals[1]) 
		elif opt in ("--testRobustness"):
			testRobustness = True
		elif opt in ("--permutationFilter"):
			permutationFilter = float(arg)
		elif opt in ("--FriLerAsCofactor"):
			FriLerAsCofactor=True
		elif opt in ("--FriColAsCofactor"):
			FriColAsCofactor=True
		elif opt in ("--useLinearRegress"): 
			useLinearRegress=True
		elif opt in ("--regressionCofactors"):
			regressionCofactors=arg
		elif opt in ("--memReq"):
			memReq=arg
		elif opt in ("--walltimeReq"):
			walltimeReq=arg
		else:
			if help==0:
				print "Unkown option!!\n"
				print __doc__
			sys.exit(2)

	if len(args)<3 and not parallel:
		if help==0:
			print "Arguments are missing!!\n"
			print __doc__
		sys.exit(2)

	print "Emma is being set up with the following parameters:"
	print "output:",rFile
	print "phenotypeRanks:",phenotypeRanks
	print "phenotypeFileType:",phenotypeFileType
	print "parallel:",parallel
	print "parallelAll:",parallelAll
	print "minMAF:",minMAF
	print "LRT:",lrt
	print "delim:",delim
	print "missingval:",missingVal
	print "kinshipDatafile:",kinshipDatafile
	print "chr:",chr
	print "boundaries:",boundaries
	print "onlyMissing:",onlyMissing
	print "onlyOriginal96:",onlyOriginal96
	print "onlyOriginal192:",onlyOriginal192
	print "onlyBelowLatidue:",onlyBelowLatidue
	print "complement:",complement
	print "negate:",negate
	print "logTransform:",logTransform
	print "addConstant:",addConstant
	print "removeOutliers:",removeOutliers
	print "sr:",sr
	print "srSkipFirstRun:",srSkipFirstRun
	print "srInput:",srInput
	print "srOutput:",srOutput
	print "srTopQuantile:",srTopQuantile
	print "srWindowSize:",srWindowSize
	print "testRobustness:",testRobustness
	print "permutationFilter:",permutationFilter
	print "useLinearRegress:",useLinearRegress
	print "regressionCofactors:",regressionCofactors
	print "FriLerAsCofactor:",FriLerAsCofactor
	print "FriColAsCofactor:",FriColAsCofactor
	print "walltimeReq:",walltimeReq
	print "memReq:",memReq

	def runParallel(phenotypeIndex,phed):
		#Cluster specific parameters
		print phenotypeIndex
		phenName = phed.getPhenotypeName(phenotypeIndex)
		outFileName = resultDir+"Emma_"+parallel+"_"+phenName

		shstr = "#!/bin/csh\n"
		shstr += "#PBS -l walltime="+walltimeReq+"\n"
		shstr += "#PBS -l mem="+memReq+"\n"
		shstr +="#PBS -q cmb\n"

		shstr += "#PBS -N E"+phenName+"_"+parallel+"\n"
		shstr += "set phenotypeName="+parallel+"\n"
		shstr += "set phenotype="+str(phenotypeIndex)+"\n"
		if useLinearRegress:
			outFileName = resultDir+"LR_"+parallel+"_"+phenName
		shstr += "(python "+emmadir+"Emma.py -o "+outFileName+" "
		if useLinearRegress:
			shstr+=" --useLinearRegress "

		if regressionCofactors:
			shstr+=" --regressionCofactors="+str(regressionCofactors)+" "			
		if FriLerAsCofactor:
			shstr+=" --FriLerAsCofactor "			
		if FriColAsCofactor:
			shstr+=" --FriColAsCofactor "			
		if onlyOriginal96:
			shstr+=" --onlyOriginal96 "			
		elif onlyOriginal192:
			shstr+=" --onlyOriginal192 "
		if onlyBelowLatidue:
			shstr+=" --onlyBelowLatidue="+str(onlyBelowLatidue)+" "
		if logTransform:
			shstr += " --logTransform "
		if negate:
			shstr += " --negate "
		if removeOutliers:
			shstr += " --removeOutliers="+str(removeOutliers)+" "
		if phenotypeRanks:
			shstr += " --phenotypeRanks "
		if testRobustness:
			shstr+=" --testRobustness "

		shstr+=" --permutationFilter="+str(permutationFilter)+" "

		if sr:
			shstr += " --sr "			
			if not srOutput:
				output = resultDir+"Emma_"+parallel+"_"+phenName+".sr.pvals"				
			shstr += " --srOutput="+str(output)+" "
			if srSkipFirstRun:
				if not srInput:
					output = resultDir+"Emma_"+parallel+"_"+phenName+".pvals"
				shstr += " --srInput="+str(output)+" "
				shstr += " --srSkipFirstRun "				
			shstr += " --srPar="+str(srTopQuantile)+","+str(srWindowSize)+" "
			
		if kinshipDatafile:
			shstr += " --kinshipDatafile="+str(kinshipDatafile)+" "			
		shstr += " --addConstant="+str(addConstant)+" "			
		shstr += snpsDataFile+" "+phenotypeDataFile+" "+str(phenotypeIndex)+" "
		shstr += "> "+outFileName+"_job"+".out) >& "+outFileName+"_job"+".err\n"

		f = open(parallel+".sh",'w')
		f.write(shstr)
		f.close()

		#Execute qsub script
		os.system("qsub "+parallel+".sh ")

	snpsDataFile = args[0]
	phenotypeDataFile = args[1]
	if parallel:  #Running on the cluster..
		phed = phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter='\t')  #Get Phenotype data 
		if parallelAll:
			for phenotypeIndex in phed.phenIds:
				if onlyMissing:
					phenName = phed.getPhenotypeName(phenotypeIndex)
					pvalFile = resultDir+"Emma_"+parallel+"_"+phenName+".pvals"
					res = None
					try:
						res = os.stat(pvalFile)

					except Exception:
						print "File",pvalFile,"does not exist."
					if res and res.st_size>0:
						print "File",pvalFile,"already exists, and is non-empty."
						if sr:
							srInput = resultDir+"Emma_"+parallel+"_"+phenName+".sr.pvals"
							srRes = None
							try:
								srRes = os.stat(srInput)
							except Exception:
								print "File",srInput,"does not exist."
							if srRes and srRes.st_size>0:
								print "File",srInput,"already exists, and is non-empty."
							else:
								runParallel(phenotypeIndex,phed)
							
					else:
						print "Setting up the run."
						runParallel(phenotypeIndex,phed)
											
				else:
					runParallel(phenotypeIndex,phed)
		else:
			phenotypeIndex = int(args[2])
			runParallel(phenotypeIndex,phed)
		return
	else:
		phenotypeIndex = int(args[2])


	print "phenotypeIndex:",phenotypeIndex
	print "\nStarting program now!\n"



	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=delim, missingVal=missingVal)

	#Load phenotype file
	phed = phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter='\t')  #Get Phenotype data 
	numAcc = len(snpsds[0].accessions)

	#Removing outliers
	if removeOutliers:
		print "Remoing outliers"
		phed.naOutliers(phenotypeIndex,removeOutliers)
	
	#If onlyOriginal96, then remove all other phenotypes..
	if onlyOriginal96: 
		print "Filtering for the first 96 accessions"
		original_96_ecotypes = phenotypeData._getFirst96Ecotypes_()
		original_96_ecotypes = map(str,original_96_ecotypes)
		keepEcotypes = []
		if complement:
			for acc in phed.accessions:
				if not acc in original_96_ecotypes:
					keepEcotypes.append(acc)
		else:
			keepEcotypes = original_96_ecotypes
		phed.filterAccessions(keepEcotypes)
		print "len(phed.accessions)", len(phed.accessions)
		
	if onlyOriginal192: 
		print "Filtering for the first 192 accessions"
		original_192_ecotypes = phenotypeData._getFirst192Ecotypes_()
		original_192_ecotypes = map(str,original_192_ecotypes)
		keepEcotypes = []
		if complement:
			for acc in phed.accessions:
				if not acc in original_192_ecotypes:
					keepEcotypes.append(acc)
		else:
			keepEcotypes = original_192_ecotypes
		phed.filterAccessions(keepEcotypes)
		print "len(phed.accessions)", len(phed.accessions)
	
	if onlyBelowLatidue:
		print "Filtering for the accessions which orginate below latitude",onlyBelowLatidue
		eiDict = phenotypeData._getEcotypeIdInfoDict_()
		print eiDict
		keepEcotypes = []
		for acc in phed.accessions:
			acc = int(acc)
			if eiDict.has_key(acc) and eiDict[acc][2] and eiDict[acc][2]<onlyBelowLatidue:
				keepEcotypes.append(str(acc))
			elif eiDict.has_key(acc) and eiDict[acc][2]==None:
				keepEcotypes.append(str(acc))
				
		phed.filterAccessions(keepEcotypes)
		print "len(phed.accessions)", len(phed.accessions)
	sys.stdout.write("Finished prefiltering phenotype accessions.\n")
	sys.stdout.flush()

	phenotype = phed.getPhenIndex(phenotypeIndex)

	accIndicesToKeep = []			
	phenAccIndicesToKeep = []
	#Checking which accessions to keep and which to remove .
	for i in range(0,len(snpsds[0].accessions)):
		acc1 = snpsds[0].accessions[i]
		for j in range(0,len(phed.accessions)):
			acc2 = phed.accessions[j]
			if acc1==acc2 and phed.phenotypeValues[j][phenotype]!='NA':
				accIndicesToKeep.append(i)
				phenAccIndicesToKeep.append(j)
				break	

	print "\nFiltering accessions in genotype data:"
	#Filter accessions which do not have the phenotype value (from the genotype data).
	for snpsd in snpsds:
		sys.stdout.write(".")
		sys.stdout.flush()
		snpsd.removeAccessionIndices(accIndicesToKeep)
	print ""
	print numAcc-len(accIndicesToKeep),"accessions removed from genotype data, leaving",len(accIndicesToKeep),"accessions in all."
		

	print "\nNow filtering accessions in phenotype data:"
	phed.removeAccessions(phenAccIndicesToKeep) #Removing accessions that don't have genotypes or phenotype values

	print "Verifying number of accessions: len(phed.accessions)==len(snpsds[0].accessions) is",len(phed.accessions)==len(snpsds[0].accessions)
	if len(phed.accessions)!=len(snpsds[0].accessions):
		raise Exception

	#Filtering monomorphic
	print "Filtering monomorphic SNPs"
	for snpsd in snpsds:
		print "Removed", str(snpsd.filterMonoMorphicSnps()),"Snps"

	#Remove minor allele frequencies
	if minMAF!=0:
		sys.stdout.write("Filterting SNPs with MAF<"+str(minMAF)+".")
		for snpsd in snpsds:
			sys.stdout.write(".")
			sys.stdout.flush()
			snpsd.filterMinMAF(minMAF)

	#Removing SNPs which are outside of boundaries.
	if chr:
		print "\nRemoving SNPs which are outside of boundaries."
		snpsds[chr-1].filterRegion(boundaries[0],boundaries[1])
		snpsds = [snpsds[chr-1]]
	
	#Ordering accessions in genotype data to fit phenotype data.
	print "Ordering genotype data accessions."
	accessionMapping = []
	i = 0
	for acc in phed.accessions:
		if acc in snpsds[0].accessions:
			accessionMapping.append((snpsds[0].accessions.index(acc),i))
			i += 1

	#print zip(accessionMapping,snpsds[0].accessions)
	print "len(snpsds[0].snps)",len(snpsds[0].snps)

	for snpsd in snpsds:
		sys.stdout.write(".")
		sys.stdout.flush()
		snpsd.orderAccessions(accessionMapping)
	print "\nGenotype data has been ordered."
		
	#Converting format to 01
	newSnpsds = []
	sys.stdout.write("Converting data format")
	for snpsd in snpsds:
		sys.stdout.write(".")
		sys.stdout.flush()
		newSnpsds.append(snpsd.getSnpsData(missingVal=missingVal))
	print ""


	
	print "Checking kinshipfile:",kinshipDatafile
	
	if kinshipDatafile:  #Is there a special kinship file?
		kinshipSnpsds = dataParsers.parseCSVData(kinshipDatafile, format=1, deliminator=delim, missingVal=missingVal)

		accIndicesToKeep = []			
		#Checking which accessions to keep and which to remove (genotype data).
		sys.stdout.write("Removing accessions which do not have a phenotype value for "+phed.phenotypeNames[phenotype]+".")
		sys.stdout.flush()
		for i in range(0,len(kinshipSnpsds[0].accessions)):
			acc1 = kinshipSnpsds[0].accessions[i]
			for j in range(0,len(phed.accessions)):
				acc2 = phed.accessions[j]
				if acc1==acc2 and phed.phenotypeValues[j][phenotype]!='NA':
					accIndicesToKeep.append(i)
					break	
		print accIndicesToKeep
	
		for snpsd in kinshipSnpsds:
			sys.stdout.write(".")
			sys.stdout.flush()
			snpsd.removeAccessionIndices(accIndicesToKeep)
		print ""
		print numAcc-len(accIndicesToKeep),"accessions removed from kinship genotype data, leaving",len(accIndicesToKeep),"accessions in all."
	
		print "Ordering kinship data accessions."
		accessionMapping = []
		i = 0
		for acc in snpsds[0].accessions:
			if acc in kinshipSnpsds[0].accessions:
				accessionMapping.append((kinshipSnpsds[0].accessions.index(acc),i))
				i += 1

		print zip(accessionMapping,snpsds[0].accessions)
		print "len(snpsds[0].snps)",len(snpsds[0].snps)
		
		for snpsd in kinshipSnpsds:
			sys.stdout.write(".")
			sys.stdout.flush()
			snpsd.orderAccessions(accessionMapping)
		print "Kinship genotype data has been ordered."

		newKinshipSnpsds = []
		sys.stdout.write("Converting data format")
		for snpsd in kinshipSnpsds:
			sys.stdout.write(".")
			sys.stdout.flush()
			newKinshipSnpsds.append(snpsd.getSnpsData(missingVal=missingVal))  #This data might have NAs
		print ""
		kinshipSnpsds = newKinshipSnpsds

	else:
		kinshipSnpsds = newSnpsds
		

	print "Found kinship data."

	#Ordering accessions according to the order of accessions in the genotype file
#	accessionMapping = []
#	i = 0
#	for acc in snpsds[0].accessions:
#		if acc in phed.accessions:
#			accessionMapping.append((phed.accessions.index(acc),i))
#			i += 1
#	phed.orderAccessions(accessionMapping)

	
	#Negating phenotypic values
	if negate: 
		phed.negateValues(phenotypeIndex)

	if logTransform and not phed.isBinary(phenotypeIndex) and phed.getMinValue(phenotypeIndex)<=0:
		addConstant=0

	#Adding a constant.
	if addConstant!=-1:
		if addConstant==0:
			addConstant = math.sqrt(phed.getVariance(phenotypeIndex))/10
			addConstant = addConstant - phed.getMinValue(phenotypeIndex)
			
		print "Adding a constant to phenotype:",addConstant
		phed.addConstant(phenotypeIndex,addConstant)
	
		
	
	#Log-transforming
	if logTransform:
		print "Log transforming phenotype"
		phed.logTransform(phenotypeIndex)
	#Converting phenotypes to Ranks
	elif phenotypeRanks:
		phed.transformToRanks(phenotypeIndex)
	
	if not chr:
		snpsDataset = snpsdata.SNPsDataSet(newSnpsds,[1,2,3,4,5])
		kinshipSnpsDataset = snpsdata.SNPsDataSet(kinshipSnpsds,[1,2,3,4,5])
	else:
		snpsDataset = snpsdata.SNPsDataSet(newSnpsds,[chr])
		kinshipSnpsDataset = snpsdata.SNPsDataSet(kinshipSnpsds,[chr])
		
	
	phenotypeName = phed.getPhenotypeName(phenotypeIndex)

	sys.stdout.flush()
	
	if testRobustness:
		print "Starting a robustness test"
		allSNPs = []
		for snpsd in snpsDataset.snpsDataList:
			allSNPs += snpsd.snps
		phenVals = phed.getPhenVals(phenotypeIndex)
		_robustness_test_(allSNPs,phenVals,rFile,filter=permutationFilter)
		sys.exit(0)

	if useLinearRegress:
		phenVals = phed.getPhenVals(phenotypeIndex)
		d0 = {}
		d0["phen"]=phenVals
		dh = {}
		dh["phen"]=phenVals
		import rpy, gc
		if regressionCofactors: #Adds ler and col as cofactors
			import pickle
			f = open(regressionCofactors,"r")
			co_factors = pickle.load(f)
			f.close()
			#inserting co factors into model
			for factor in co_factors:
				d[factor] = co_factors[factor]
		import analyzeHaplotype as ah
		(ler_factor, col_factor) = ah.getLerAndColAccessions(newSnpsds,True)
		if FriColAsCofactor:
			d0["col"] = col_factor
			dh["col"] = col_factor
		if FriLerAsCofactor:
			d0["ler"] = ler_factor
			dh["ler"] = ler_factor
		chr_pos_pvals = []
		stats = []
		sys.stdout.write("Applying the linear model")
		sys.stdout.flush()
		for i in range(0,len(newSnpsds)): #[3]:#
			snpsd=newSnpsds[i] 
			sys.stdout.write("|")
			sys.stdout.flush()
			gc.collect()  #Calling garbage collector, in an attempt to clean up memory..
			for j in range(0,len(snpsd.snps)):
				if j%5000==0:
					sys.stdout.write(".")
					sys.stdout.flush()
				#if snpsd.positions[j]>1700000:
				#	break
				snp = snpsd.snps[j]
				d0["snp"]=snp
				try:
					rpy.set_default_mode(rpy.NO_CONVERSION)
					aov0 = rpy.r.aov(r("phen ~ ."), data = d0)
					aovh = rpy.r.aov(r("phen ~ ."), data = dh)
					rpy.set_default_mode(rpy.BASIC_CONVERSION)
					s0 = rpy.r.summary(aov0)
					sh = rpy.r.summary(aovh)
					#print s0,sh
					rss_0 = s0['Sum Sq'][-1]
					if type(sh['Sum Sq'])!=float:
						rss_h = sh['Sum Sq'][-1]
						
					else:
						rss_h = sh['Sum Sq']
					f = (rss_h-rss_0)/(rss_0/(len(phenVals)-len(d0)+1))
					pval = rpy.r.pf(f,1, len(phenVals), lower_tail = False)
				except Exception, err_str:
					print "Calculating p-value failed"#,err_str
					pval =1.0
				#print "dh:",dh
				#print "d0:",d0
				#print "rss_h,rss_0:",rss_h,rss_0
				#print "f,p:",f,pval
				chr_pos_pvals.append([i+1,snpsd.positions[j],pval])
				mafc = min(snp.count(snp[0]),len(snp)-snp.count(snp[0]))
				maf = mafc/float(len(snp))
				stats.append([maf,mafc])
		sys.stdout.write("\n")
		#Write out to a result file
		sys.stdout.write("Writing results to file\n")
		sys.stdout.flush()
		pvalFile = rFile+".pvals"
		f = open(pvalFile,"w")
		f.write("Chromosome,position,p-value,marf,maf\n")
		for i in range(0,len(chr_pos_pvals)):
			chr_pos_pval = chr_pos_pvals[i]
			stat = stats[i]
			f.write(str(chr_pos_pval[0])+","+str(chr_pos_pval[1])+","+str(chr_pos_pval[2])+","+str(stat[0])+","+str(stat[1])+"\n")		
		f.close()
		
		#Plot results
		print "Generating a GW plot."
		phenotypeName=phed.getPhenotypeName(phenotypeIndex)
		res = gwaResults.Result(pvalFile,name="LM_"+phenotypeName, phenotypeID=phenotypeIndex)
		res.negLogTransform()
		pngFile = pvalFile+".png"
		plotResults.plotResult(res,pngFile=pngFile,percentile=90,type="pvals",ylab="$-$log$_{10}(p)$", plotBonferroni=True,usePylab=False)	
		
		
	elif (not sr) or (sr and not srSkipFirstRun):
		sys.stdout.write("Running Primary Emma.\n")
		sys.stdout.flush()
		pvalFile = _runEmmaScript_(snpsDataset, kinshipSnpsDataset, phed, phenotypeIndex, rFile, chr=chr, delim=delim, missingVal=missingVal, boundaries=boundaries, lrt=lrt)
		res = gwaResults.Result(pvalFile,name="EMMA_"+phenotypeName, phenotypeID=phenotypeIndex)
		res.filterMAF()
		res.negLogTransform()
		pngFile = pvalFile+".png"
		plotResults.plotResult(res,pngFile=pngFile,percentile=90,type="pvals",ylab="$-$log$_{10}(p)$", plotBonferroni=True,usePylab=False)	
		srInput = pvalFile

	if sr:
		_secondRun_(srOutput,srInput,srTopQuantile,srWindowSize,newSnpsds,phed,phenotypeIndex,kinshipSnpsDataset)
		print "Generating second run GW plot."
		res = gwaResults.Result(srInput,name="KW_"+phenotypeName, phenotypeID=phenotypeIndex)
		res.filterMAF()
		res.negLogTransform()
		srRes = gwaResults.Result(srOutput,name="EMMA_SR_"+phenotypeName, phenotypeID=phenotypeIndex)
		srRes.filterMAF()
		srRes.negLogTransform()
		srPngFile = pvalFile+".sr.png"
		plotResults.plotResultWithSecondRun(res,srRes,pngFile=srPngFile,ylab="$-$log$_{10}(p)$", plotBonferroni=True)	
		
		
		

def _runEmmaScript_(snpsDataset, kinshipSnpsDataset, phed, p_i, rFile, chr=None, delim=",", missingVal="NA", boundaries=[-1,-1], lrt=False):

	#Writing files
	cluster = "/home/cmb-01/bvilhjal/"==env.home_dir #Am I running on the cluster?
	#tempfile.tempdir = "/home/cmb-01/bvilhjal/tmp/" #(Temporary) debug hack...
	if not cluster:
		tempfile.tempdir='/tmp'
	(fId, phenotypeTempFile) = tempfile.mkstemp()
	os.close(fId)
	(fId, genotypeTempFile) = tempfile.mkstemp()
	os.close(fId)
	(fId, kinshipTempFile) = tempfile.mkstemp()
	os.close(fId)
	
	phed.writeToFile(phenotypeTempFile, [phed.getPhenIndex(p_i)])	
	sys.stdout.write( "Phenotype file written\n")
	sys.stdout.flush()


	#Put into a function....
	snpsDataset.writeToFile(genotypeTempFile, deliminator=delim, missingVal = missingVal)
	sys.stdout.write( "Genotype file written\n")
	sys.stdout.flush()
	kinshipSnpsDataset.writeToFile(kinshipTempFile, deliminator=delim, missingVal = missingVal)
	sys.stdout.write( "Kinship genotype file written\n")
	sys.stdout.flush()

	phenotypeName = phed.getPhenotypeName(p_i)

	rDataFile = rFile+".rData"
	pvalFile = rFile+".pvals"
	rstr = _generateRScript_(genotypeTempFile, phenotypeTempFile, rDataFile, pvalFile, kinshipTempFile, name = phenotypeName, boundaries = boundaries, chr=chr,cluster=cluster, lrt=lrt)
	rFile2 = rFile+".R"
	f = open(rFile2,'w')
	f.write(rstr)
	f.close()
	outRfile = rFile+"_R.out"
	errRfile = rFile+"_R.err"
	print "Running R file:"
	cmdStr = "(R --vanilla < "+rFile2+" > "+outRfile+") >& "+errRfile
	sys.stdout.write(cmdStr+"\n")
	sys.stdout.flush()
	os.system(cmdStr)
	print "Emma output saved in R format in", rDataFile
	return pvalFile
	
	
def _generateRScript_(genotypeFile, phenotypeFile, rDataFile, pvalFile, kinshipDatafile, name=None, boundaries=[-1,-1],chr=None, cluster=False, lrt=False):
	
#	if cluster:
#		rstr = 'library(emma,lib.loc="/home/cmb-01/bvilhjal/Projects/emma");\n'
#	else:
#		rstr = "library(emma);\n"
	rstr = 'source("emma.R");\n'  #Using modified version of Emma (from Yu).
	rstr += 'data250K <- read.csv("'+str(genotypeFile)+'", header=TRUE);\n'
	rstr += "mat250K <- as.matrix(data250K);\n"
	rstr += 'data2010 <- read.csv("'+kinshipDatafile+'", header=TRUE);\n'
	rstr += "mat2010 <- as.matrix(data2010);\n"
	rstr += 'phenotData <- read.csv("'+str(phenotypeFile)+'",header=TRUE);\n'
	rstr += """
phenMat <- t(as.matrix(as.matrix(phenotData)[,2]));
res <- list();
mat250KAll <- mat250K[,3:length(mat250K[1,])]
mat2010All <- mat2010[,3:length(mat2010[1,])]
K2010<- emma.kinship(mat2010All,use="all") #Kinship matrix
"""
	
	rstr += """
pvals <- c();
genotype_var_perc <- c();
remls <- c();
vgs <- c();
ves <- c();
positions <- c();
chrs <- c();
maf <- c();
marf <- c();
"""
	if chr:
		rstr += "for (chr in ("+str(chr)+")){"
	else:
		rstr += "for (chr in (1:5)){"
			
	rstr += """
  mat250K <- as.matrix(data250K);
  mat250K <- mat250K[mat250K[,1]==chr,];

  pos <- mat250K[,2];
  mat250K <- mat250K[,3:length(mat250K[1,])];
"""

	if lrt: 
		rstr += "  res[[chr]] <- emma.ML.LRT(phenMat,mat250K,K2010);"
	else:
		rstr += "  res[[chr]] <- emma.REML.t(phenMat,mat250K,K2010);"
	rstr +="""
  res[[chr]][["pos"]] <- pos;
  res[[chr]][["chr_pos"]] <- pos;

  pvals <- append(pvals,res[[chr]]$ps);
  remls <- append(remls,res[[chr]]$REMLs)
  vgs <- append(vgs,res[[chr]]$vgs)
  ves <- append(ves,res[[chr]]$ves)
  genotype_var_perc <- append(genotype_var_perc,res[[chr]]$genotype_var_perc)
  positions <- append(positions,pos);
  chrs <- append(chrs,rep(chr,length(pos)));
  af <- c();
  arf <- c();
  for(i in (1:length(mat250K[,1]))){
    f = factor(mat250K[i,]);
    freq <- summary(f)[[1]]
    alleleCount <- length(f)
    v <- freq/alleleCount;
    af[i] <- min(freq,alleleCount-freq);
    arf[i] <- min(v,1-v);
  }
  res[[chr]][["maf"]] <- af;
  res[[chr]][["marf"]] <- arf;

  marf <- append(marf,arf);
  maf <- append(maf,af);
}
res[["K"]] <- K2010;

#write to a pvalue-file
l <- list();
l[["Chromasome"]]<-chrs;
l[["Positions"]]<-positions;
l[["Pvalues"]]<-pvals;
l[["MARF"]]<-marf;
l[["MAF"]]<-maf;
l[["genotype_var_perc"]]<-genotype_var_perc;
l[["REMLs"]]<-remls;
l[["vgs"]]<-vgs;
l[["ves"]]<-ves;
dl <- as.data.frame(l);
"""
	rstr +=' write.table(dl,file="'+pvalFile+'", sep=", ", row.names = FALSE);\n'		
	
	#111708 - Following code (R object) was removed. 
#	rstr += """
##Save data as R object.
#res[[2]]$pos <- res[[2]]$pos+res[[1]]$pos[length(res[[1]]$pos)];
#res[[3]]$pos <- res[[3]]$pos+res[[2]]$pos[length(res[[2]]$pos)];
#res[[4]]$pos <- res[[4]]$pos+res[[3]]$pos[length(res[[3]]$pos)];
#res[[5]]$pos <- res[[5]]$pos+res[[4]]$pos[length(res[[4]]$pos)];
#
#res[["ylim"]] <- c(min(min(-log(res[[1]]$ps)),min(-log(res[[2]]$ps)),min(-log(res[[3]]$ps)),min(-log(res[[4]]$ps)),min(-log(res[[5]]$ps))), max(max(-log(res[[1]]$ps)),max(-log(res[[2]]$ps)),max(-log(res[[3]]$ps)),max(-log(res[[4]]$ps)),max(-log(res[[5]]$ps))));
#
#res[["xlim"]] <- c(min(res[[1]]$pos),max(res[[5]]$pos));
#
#res[["FRI"]]<-c(269026+res[[3]]$pos[length(res[[3]]$pos)],271503+res[[3]]$pos[length(res[[3]]$pos)]);
#res[["FLC"]]<-c(3173498+res[[4]]$pos[length(res[[4]]$pos)],3179449+res[[4]]$pos[length(res[[4]]$pos)]);
#
#
#"""
#	if name:
#		rstr += 'res[["lab"]]= "Emma p-values for '+name+'";\n'
#	else:
#		rstr += 'res[["lab"]]="";\n'
#	rstr += 'save(file="'+rDataFile+'",res);\n'
	return rstr	


def _secondRun_(srOutput,srInput,srTopQuantile,srWindowSize,snpsds,phed,p_i,kinshipSnpsDataset):
	"""
	Sets up a second run, using and/or on SNPs.
	"""
	print "Preparing a second stage analyzis."
	sys.stdout.flush()
	#Collect top n% SNPs.
	print "Loading first stage results."
	sys.stdout.flush()
	result = gwaResults.SNPResult(resultFile=srInput,snpsds=snpsds,phenotypeID=p_i)
	# - filter results
	result.negLogTransform()
		
	res = SecondStageAnalysis.retrieveSecondRunSNPs(result,srTopQuantile,srWindowSize)
	snps = res['snps']
	snpTypes = res['snpTypes']
	chromosomes = res['chromosomes']
	positions = res['positions']
	marfs = res['marfs']
	mafs = res['mafs']
	
	
	#Calc kinship matrix:
	print "Calculating kinship matrix."
	sys.stdout.flush()
	k_snps = kinshipSnpsDataset.getSnps()
	#k_snps = _sampleSNPs_(k_snps,10000)
	k = calcKinship(k_snps)
	#print "K:",k
	
	#Run Emma:
	sys.stdout.write( "Running Emma.\n")
	sys.stdout.flush()
	res = runEmma(phed,p_i,k,snps)
	pvals = list(res["ps"])
	genotype_var_perc = list(res["genotype_var_perc"])
	#print zip(positions,pvals,genotype_var_perc)
	
	sys.stdout.write("Writing results to file.\n")
	sys.stdout.flush()
	#Write results to file!
	f = open(srOutput,"w")
	f.write("Chromosome,position,p-value,marf,maf,genotype_var_perc,snpType,second_pos\n")
	for (chr,pos,pval,marf,maf,gvp,snpType) in zip(chromosomes,positions,pvals,marfs,mafs,genotype_var_perc,snpTypes):
		f.write(str(chr)+","+str(pos[0])+","+str(pval[0])+","+str(marf)+","+str(maf)+","+str(gvp[0])+","+str(snpType)+","+str(pos[1])+"\n")
	f.close()

	#Plot results?
	
	
def getKinshipMatrix():
	#snpsDataFile="/Network/Data/250k/dataFreeze_011209/250K_f13_012509.csv"
	snpsDataFile="/home/cmb-01/bvilhjal/Projects/data/250K_f13_012609.csv"
	import dataParsers,snpsdata
	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=",")#,debug=True)
	snps = []
	sys.stdout.write("Converting format")
	for snpsd in snpsds:
		sys.stdout.write(".")
		sys.stdout.flush()
		snps += snpsd.getSnpsData(missingVal="NA").snps
	print ""
	#snps = _sampleSNPs_(snps,100)
	print "Calculating kinship"
	K = calcKinship(snps)
	eDict = phenotypeData._getEcotypeIdToStockParentDict_()
	accessions = map(int,snpsd.accessions)
	#for et in accessions:
	#print eDict[et]
	for i in range(0,len(accessions)):
		et = accessions[i]
		info = eDict[et]
		st = str(et)+", "+str(info[0])+", "+str(info[1])+":"
		st += str(K[i][0])
		for j in range(1,i+1):
			st += ", "+str(K[i][j])
		print st



	
def _robustness_test_(all_snps,phenVals,outputFile,filter=0.1):
	import analyzePhenotype
	import analyzeHaplotype
	import random
	
	new_all_snps = []
	for snp in all_snps:
		if snp.count(0)>1 and snp.count(1)>1:
			new_all_snps.append(snp)
	print "Filtered",len(all_snps)-len(new_all_snps)," with minor allele count <2."
	all_snps = new_all_snps
	

	def getLeaveOneOutK(K,i):
		l = range(0,len(K))
		l.pop(i)
		new_k = numpy.core.take(K,l,0)
		new_k = numpy.core.take(new_k,l,1)
		return new_k
	
	print "Calculating kinship"
	t1 = time.time()
	K = calcKinship(all_snps)
	t2 = time.time()
	print "Took",t2-t1,"seconds."

	
	"""
	Leave one out test..
	"""
	if filter <1.0:
		snps = random.sample(all_snps,int(len(all_snps)*filter))
		print "Number of SNPs:",len(snps)
	else:
		snps = all_snps 
#	K = calcKinship(snps)

	
	print "running EMMA"
	t1 = time.time()
	true_pvals = _runEmma_(snps,phenVals,K)["ps"]
	true_pvals = map(float,true_pvals)
	t2 = time.time()
	print "Took",t2-t1,"seconds."

	
	log_true_pvals = []
	for pval in true_pvals:
		log_true_pvals.append(-math.log(pval,10))
	
	perm_pvalues_list = []
	for i in range(0,len(phenVals)):
		newPhenvals = phenVals[:]
		newPhenvals.pop(i)
				
		newSNPs = []
		for snp in snps:
			newSNP = snp[:]
			newSNP.pop(i)
			newSNPs.append(newSNP)
		
		new_k = getLeaveOneOutK(K,i)
		
		print "running EMMA"
		t1 = time.time()
		pvals = _runEmma_(newSNPs,newPhenvals,new_k)["ps"]
		pvals = map(float,pvals)
		t2 = time.time()
		print "Took",t2-t1,"seconds."

		perm_pvalues_list.append(pvals)
		
		
	
	delta_pvals_list = []
	delta_log_pvals_list = []
	for perm_pvals in perm_pvalues_list:
		log_pvals = []
		delta_pvals = []
		delta_log_pvals = []
		for i in range(0,len(true_pvals)):
			pval = perm_pvals[i]
			true_pval = true_pvals[i]
			delta_pvals.append(true_pval-pval)

			log_true_pval = log_true_pvals[i]
			log_pval = -math.log(pval,10)
			log_pvals.append(log_pval)
			delta_log_pvals.append(log_true_pval-log_pval)
		
		delta_pvals_list.append(delta_pvals)
		delta_log_pvals_list.append(delta_log_pvals)
	
	sd_log_pvals = []
	sd_pvals = []
	t_delta_log_pvals_list = map(list,zip(*delta_log_pvals_list))
	t_delta_pvals_list = map(list,zip(*delta_pvals_list))
	for i in range(0,len(true_pvals)):
		sd_log_pvals.append(util.calcSD(t_delta_log_pvals_list[i]))
		sd_pvals.append(util.calcSD(t_delta_pvals_list[i]))
	
	
	#Write SDs out to file, to be able to replot, or plot together with other methods... etc
	import csv
	sd_log_pval_file = outputFile+".rob.log_pvals_sd"
	f = open(sd_log_pval_file,"w")
	w = csv.writer(f)
	w.writerow(["log_true_pval","sd_log_pvals"])
	l = zip(log_true_pvals,sd_log_pvals)
	w.writerows(l)
	f.close()
	
	
	#Plot things....
	pngFile_log_pvals = outputFile+".rob.log_pval.png"
	pngFile_pval = outputFile+".rob.pval.png"
	pngFile_sd_log_pval = outputFile+".rob.sd_log_pval.png"
	pngFile_sd_pval = outputFile+".rob.sd_pval.png"


	min_val = min(true_pvals)
	max_val = max(true_pvals)
	val_range = max_val-min_val
	
	min_log_val = min(log_true_pvals)
	max_log_val = max(log_true_pvals)
	log_val_range = max_val-min_val


	print "Plotting graphs"
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	plt.figure(figsize=(10,7))
	max_perm_val = 0
	min_perm_val = 0
	for i in range(0,len(perm_pvalues_list)):
		delta_log_pvals = delta_log_pvals_list[i]
		plt.plot(log_true_pvals,delta_log_pvals,"b.")
		max_perm_val = max(max_perm_val,max(delta_log_pvals))
		min_perm_val = min(min_perm_val,min(delta_log_pvals))
	perm_val_range = max_perm_val - min_perm_val
	v = [min_log_val-0.02*log_val_range, max_log_val+0.02*log_val_range, min_perm_val-0.02*perm_val_range, max_perm_val+0.02*perm_val_range]
	plt.axis(v)
	plt.savefig(pngFile_log_pvals, format = "png")

	plt.figure(figsize=(10,7))
	max_perm_val = 0
	min_perm_val = 0
	for i in range(0,len(perm_pvalues_list)):
		delta_pvals = delta_pvals_list[i]
		plt.plot(true_pvals,delta_pvals,"b.")
		max_perm_val = max(max_perm_val,max(delta_pvals))
		min_perm_val = min(min_perm_val,min(delta_pvals))
	perm_val_range = max_perm_val - min_perm_val
	plt.axis([min_val-0.02*val_range, max_val+0.02*val_range, min_perm_val-0.02*perm_val_range, max_perm_val+0.02*perm_val_range])
	plt.savefig(pngFile_pval, format = "png")

	plt.figure(figsize=(10,7))
	max_sd_log_pval = max(sd_log_pvals)
	min_sd_log_pval = min(sd_log_pvals)
	sd_val_range = max_sd_log_pval-min_sd_log_pval
	plt.plot(log_true_pvals,sd_log_pvals,"b.")
	plt.axis([min_log_val-0.02*log_val_range, max_log_val+0.02*log_val_range, min_sd_log_pval-0.02*sd_val_range, max_sd_log_pval+0.02*sd_val_range])
	plt.savefig(pngFile_sd_log_pval, format = "png")

	plt.figure(figsize=(10,7))
	max_sd_pval = max(sd_pvals)
	min_sd_pval = min(sd_pvals)
	sd_val_range = max_sd_pval-min_sd_pval
	plt.plot(true_pvals,sd_pvals,"b.")
	plt.axis([min_val-0.02*val_range, max_val+0.02*val_range, min_sd_pval-0.02*sd_val_range, max_sd_pval+0.02*sd_val_range])
	plt.savefig(pngFile_sd_pval, format = "png")
	print "Done testing robustness"


if __name__ == '__main__':
	_run_()
	#getKinshipMatrix()
	print "Done!"


