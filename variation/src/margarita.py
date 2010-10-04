#!/usr/bin/env python
"""
Usage: margarita.py [OPTIONS] SNPS_DATA_FILE PHENOTYPE_DATA_FILE PHENOTYPE_INDEX

Option:
	-i ..., --id=...		the run ID for this Margarita run.  (temporary filenames will be generated based on this ID).
	-s ..., --scoreFile=...	 	ARG mapping score file
	-d ..., --delim=...		default is ", "	  
	-m ..., --missingval=...	default is "NA"
	-c ..., --chr=...		default is all chromosomes
	--parallel=...			Creates the PBS (portable batch system) to set up the run with the given run ID.  
					(Generates a script which it executes).
	--parallelAll			Runs margarita on all phenotypes (only on cluster).
	--numARG=...			default is 30
	--numMarkers=...		default is 100
	--numPerm=...			default is 0
	--smartCutoff=...		default is 10
	--BoundaryStart=...		Only the region within the boundary is considered in the GWA. (Default is no boundaries)
	--BoundaryEnd=...		   
	--phenotypeFileType=...	 	1 (default) if file has tsv format, 2 if file has csv format and contains accession names (instead of ecotype ID)
	--binary			If the phenotype is binary.
	-a ..., --withArrayId=...   	0 for no array ID info (default), 1 if file has array ID info.
	-b, --debug			enable debug
	-h, --help			show this help

Examples:
	margarita.py   /tmp/250K.csv phenotypes.tsv 2
	
Description:
	An python wrapper for Margarita (blended)

Dependencies:
	Requires MySQLdb to be installed, and uses snpsdata.py and dataParsers.py.
"""

## Change log (started 10/23/08)


#from dataParser import *
import env

class Margarita:
	"""
	A wrapper class for margarita, that aids whole genome ananlysis using margarita.

	Requires both margarita and Java to be installed appropriately.
	"""

	_FLC_ = [5, 3173498, 3179449]
	_FRI_ = [4, 269026, 271503]
	_knownLoci_ = [("FRI", _FRI_), ("FLC", _FLC_)]

	def __init__(self, margFile, outFile, numARG = 100, numMarkers = 200, numPerm = 0, smartCutoff = 12):
		self.margFile = margFile
		self.outfile = outFile	
		self.numARG = numARG
		self.numMarkers = numMarkers
		self.numPerm = numPerm
		self.smartCutoff = smartCutoff
	

	"""
	/System/Library/Frameworks/JavaVM.framework/Versions/1.6.0/Home/bin/javac -Xlint:unchecked 
	-cp /Users/bjarnivilhjalmsson/Projects/gwas_programs/margarita/:/Users/bjarnivilhjalmsson/Projects/gwas_programs/ 
	/Users/bjarnivilhjalmsson/Projects/gwas_programs/margarita/Main.java
	"""
	#Final variables.
	import env
	margaritaCommand = env.java_dir + 'java -cp ' + env.margarita_dir + 'margarita/jsci-core.jar:' + env.margarita_dir + ' margarita/Main '
	def _generateMargaritaFile_():
		"""IMPLEMENT"""

	def _calcVar(self, list):
		mean = sum(list) / float(len(list))
		var = 0
		for i in range(0, len(list)):
			var = var + (list[i] - mean) * (list[i] - mean)
		var = var / float(len(list) - 1)
		return var

		 

	def gwa(self, snpsd, phed, phenotype = 1, chromosome = 1, boundaries = None, binary = False):
		import random, sys

		print "Preparing to run Margarita on phenotype", phed.getPhenotypeName(phenotype), "on chromosome", chromosome
		phenIndex = phed.getPhenIndex(phenotype)
		print "Internal phenotype index is", phenIndex
	  
		runsets = [] #List of run-batches..
		i = 0

		if boundaries:
			start = 0
			while snpsd.positions[start] < boundaries[0]:
				start = start + 1
			end = start
			while end < len(snpsd.positions) - 1 and snpsd.positions[end] < boundaries[1]:
				end = end + 1
			i = start
			while i < end - self.numMarkers:
				l = range(i, i + self.numMarkers)
				runsets.append(l)
				i = i + self.numMarkers
			runsets.append(range(i, end))

		else: #no boundaries
			while i < len(snpsd.positions) - 2 - self.numMarkers:
				l = range(i, i + self.numMarkers)
				runsets.append(l)
				i = i + self.numMarkers
			runsets.append(range(i, len(snpsd.positions)))

		runsetsList = [runsets]

		#Accessions with the phenotype
		accessionsIndices = []
		commonAccessions = []
		for i in range(0, len(snpsd.accessions)):
			acc1 = snpsd.accessions[i]
			for j in range(0, len(phed.accessions)):
				acc2 = phed.accessions[j]
				if acc1 == acc2 and phed.phenotypeValues[j][phenIndex] != 'NA':
					accessionsIndices.append([i, j])
					commonAccessions.append(acc1)

		print len(commonAccessions), "accessions found with phenotype values for", phed.phenotypeNames[phenIndex] + '.'

		# Variables to capture the results
		pos = []
		permPvals = []
		chiPvals = []
		mapScores = []

		chromasomesUnderCut = []
		meanNumCases = []
		caseFrequences = []
		controlFrequences = []

		bestcutScores = [] 
		bestcutVariances = []
		mafs = []
		marfs = []

		decoder = {0:'0', 1:'1', - 1:'M', 2:'M'}

		print "Running margarita:"
		for runsets in runsetsList:
			for runset in runsets:  #Loop over all run batches.
				sys.stdout.write(".")
				sys.stdout.flush()
				#Writing margarita file!
				st = ""
				st += str(len(accessionsIndices)) + " " + str(len(runset)) + "\n"
				for i in runset:
					st += str(snpsd.positions[i]) + "\n"	
				for i in range(0, len(accessionsIndices)):
					st += str(phed.phenotypeValues[accessionsIndices[i][1]][phenIndex]) + "\n"  #The phenotype
				for i in range(0, len(accessionsIndices)):
					for j in runset:
						st += str(decoder[snpsd.snps[j][accessionsIndices[i][0]]])   #Decoding the SNP
					st += "\n"

				f = open(self.margFile, "w")
				f.write(st)
				f.close()

				import os
				done = False
				numTries = 0
				while (numTries < 10):   #Check whether run was ok.  (Some bug in margrita causes it to crash occasionally.)
					numTries = numTries + 1
					if binary:
						command = self.margaritaCommand + " " + self.margFile + ' ' + str(self.numARG) + ' ' + str(self.numPerm) + ' -binary -smart ' + str(self.smartCutoff) + ' > ' + self.outfile
					else:
						command = self.margaritaCommand + " " + self.margFile + ' ' + str(self.numARG) + ' ' + str(self.numPerm) + ' -smart ' + str(self.smartCutoff) + ' > ' + self.outfile
					
					os.system(command)
					try:
						f = open(self.outfile, "r")
						#print self.outfile
						lines = f.readlines()
						f.close()
						line = lines[0]
						i = 0
						while not lines[i].startswith("%INTERPRETATION"):
							i = i + 1
						i = i + 2  
						interpretationResults = []
						while not lines[i].startswith("%MAPPING"):
							interpretationResults.append((lines[i].rstrip()).split())
							i = i + 1
						i = i + 1
						labels2 = (lines[i].rstrip()).split()
						i = i + 1
						mappingResults = []
						while i < len(lines):
							mappingResults.append((lines[i].rstrip()).split())
							i = i + 1
						
						if len(mappingResults) != len(runset):
							print "mappingResults length:", len(mappingResults), ", runset length:", len(runset)
							#raise Exception()
							if numTries == 9:
								break
						else:
							break
					except IndexError:
						print "failed", numTries
				if numTries == 10:
					print "runset:",runset
					print "len(runset):",len(runset)
					print "binary:",binary
					raise Exception

				newChromasomesUnderCut = ["%UNDERCUTS\n"]
				for j in range(0, len(interpretationResults), self.numARG):
					bslist = []
					vals = [0] * len(accessionsIndices)  
					caseSum = 0
					caseFrequency = 0
					controlFrequency = 0
					for k in range(0, self.numARG):
						result = interpretationResults[j + k]
						caseFrequency = caseFrequency + float(result[6])
						controlFrequency = controlFrequency + float(result[7])
						bslist.append(float(result[3]))
						newChromasomesUnderCut.append("MARKER: " + str(j / self.numARG) + " ARG: " + str(k) + " UNDERCUT: " + result[5] + "\n")
						underCut = []
						for m in range(0, len(result[5])):
							val = int(result[5][m])
							underCut.append(val)
							vals[m] = vals[m] + val
						caseSum = caseSum + sum(underCut)
					

					meanNumCases.append(str(caseSum / float(self.numARG)))
					bestcutScores.append(str(max(bslist)))
					bestcutVariances.append(str(self._calcVar(bslist)))
					caseFrequency = caseFrequency / self.numARG
					caseFrequences.append(str(caseFrequency))
					controlFrequency = controlFrequency / self.numARG
					controlFrequences.append(str(controlFrequency))
				
				for result in mappingResults:
					if float(result[3]) == 0.0:
						res = str(float(1.0) / float(self.numPerm * 2))
					else:
						res = result[3]
					permPvals.append(res)
					mapScores.append(result[2])
					chiPvals.append(result[4])
					pos.append(result[1])

				#Calculating MAFs
				for j in runset:
					totalCount = 0 
					alleleCount = 0
					for i in range(0, len(accessionsIndices)):
						if snpsd.snps[j][accessionsIndices[i][0]] == 1:
							alleleCount += 1
							totalCount += 1 
						elif snpsd.snps[j][accessionsIndices[i][0]] == 0:
							totalCount += 1
					arf = 0
					if totalCount > 0:
						arf = alleleCount / float(totalCount)
					marfs.append(min(arf, 1.0 - arf))
					mafs.append(min(alleleCount, totalCount - alleleCount))


				#Loop over all batches ends.
			#Loop over runsets ends.
		
		txtStr = ""  #pos, scores, and MAF text files
		rstr = ""
		for j in range(0, len(pos)):
			txtStr = txtStr + str(chromosome) + ", " + str(int(float(pos[j]))) + ", " + str(mapScores[j]) + ", " + str(marfs[j]) + ", " + str(mafs[j]) + "\n"
				

		return (rstr, txtStr, permPvals)
		


	def gwaWithTrees(self, id, snpsd, phed, phenotype = 0, chromosome = None, boundaries = None):
		import random, sys

		print "Preparing to run Margarita on phenotype", phed.getPhenotypeName(phenotype), "on chromosome", chromosome
		phenIndex = getPhenIndex(phenotype)
		print "Internal phenotype index is", phenIndex
	  
		runsets = [] #List of run-batches..
		i = 0

		if boundaries:
			start = 0
			while snpsd.positions[start] < boundaries[0]:
				start = start + 1
			end = start
			while end < len(snpsd.positions) - 1 and snpsd.positions[end] < boundaries[1]:
				end = end + 1
			i = start
			while i < end - self.numMarkers:
				l = range(i, i + self.numMarkers)
				runsets.append(l)
				i = i + self.numMarkers
			runsets.append(range(i, end))

		else: #no boundaries
			while i < len(snpsd.positions) - self.numMarkers:
				l = range(i, i + self.numMarkers)
				runsets.append(l)
				i = i + self.numMarkers
			runsets.append(range(i, len(snpsd.positions)))

		runsetsList = [runsets]

		#Accessions with the phenotype
		accessionsIndices = []
		commonAccessions = []
		for i in range(0, len(snpsd.accessions)):
			acc1 = snpsd.accessions[i]
			for j in range(0, len(phed.accessions)):
				acc2 = phed.accessions[j]
				if acc1 == acc2 and phed.phenotypeValues[j][phenIndex] != 'NA':
					accessionsIndices.append([i, j])
					commonAccessions.append(acc1)

		print len(commonAccessions), "accessions found with phenotype values for", phed.phenotypeNames[phenIndex] + '.'

		# Variables to capture the results
		pos = []
		permPvals = []
		chiPvals = []
		mapScores = []

		chromasomesUnderCut = []
		meanNumCases = []
		caseFrequences = []
		controlFrequences = []

		bestcutScores = [] 
		bestcutVariances = []
		mafs = []
		marfs = []

		decoder = {0:'0', 1:'1', - 1:'M', 2:'M'}

		print "Running margarita:"
		for runsets in runsetsList:
			for runset in runsets:  #Loop over all run batches.
				sys.stdout.write(".")
				sys.stdout.flush()
				#Writing margarita file!
				st = ""
				st += str(len(accessionsIndices)) + " " + str(len(runset)) + "\n"
				for i in runset:
					st += str(snpsd.positions[i]) + "\n"	
				for i in range(0, len(accessionsIndices)):
					st += str(phed.phenotypeValues[accessionsIndices[i][1]][phenIndex]) + "\n"  #The phenotype
				for i in range(0, len(accessionsIndices)):
					for j in runset:
						st += str(decoder[snpsd.snps[j][accessionsIndices[i][0]]])   #Decoding the SNP
					st += "\n"

				f = open(self.margFile, "w")
				f.write(st)
				f.close()

				print "wrote to", self.margFile, ":", st
				import os
				done = False
				numTries = 0
				while (numTries < 10):				 #Check whether run was ok.  (Some bug in margrita causes it to crash sometimes.)
					numTries = numTries + 1
					if binary:
						command = self.margaritaCommand + " " + self.margFile + ' ' + str(self.numARG) + ' ' + str(self.numPerm) + ' -binary -smart ' + str(self.smartCutoff) + ' > ' + self.outfile
						#print command
						os.system(command)
					else:
						os.system(self.margaritaCommand + " " + self.margFile + ' ' + str(self.numARG) + ' ' + str(self.numPerm) + ' -smart ' + str(self.smartCutoff) + ' > ' + self.outfile) 
					try:
						f = open(outfile, "r")
						lines = f.readlines()
						f.close()
						line = lines[0]
						i = 0
						while not lines[i].startswith("%TREES"):
							i = i + 1
							newTreeList = []
						while not lines[i].startswith("%INTERPRETATION"):
							newTreeList.append(lines[i])
							i = i + 1
						i = i + 2  
						interpretationResults = []
						while not lines[i].startswith("%MAPPING"):
							interpretationResults.append((lines[i].rstrip()).split())
							i = i + 1
						i = i + 1
						labels2 = (lines[i].rstrip()).split()
						i = i + 1
						mappingResults = []
						while i < len(lines):
							mappingResults.append((lines[i].rstrip()).split())
							i = i + 1
						break
					except IndexError:
						print "failed", numTries
				if numTries == 10:
					raise Exception

				newChromasomesUnderCut = ["%UNDERCUTS\n"]
				for j in range(0, len(interpretationResults), numArg):
					bslist = []
					vals = [0] * len(data.individuals)
					caseSum = 0
					caseFrequency = 0
					controlFrequency = 0
					for k in range(0, numArg):
						result = interpretationResults[j + k]
						caseFrequency = caseFrequency + float(result[6])
						controlFrequency = controlFrequency + float(result[7])
						bslist.append(float(result[3]))
						newChromasomesUnderCut.append("MARKER: " + str(j / numArg) + " ARG: " + str(k) + " UNDERCUT: " + result[5] + "\n")
						underCut = []
						for m in range(0, len(result[5])):
							val = int(result[5][m])
							underCut.append(val)
							vals[m] = vals[m] + val
						caseSum = caseSum + sum(underCut)

					meanNumCases.append(str(caseSum / float(numArg)))
					chromasomesUnderCut.append(vals)
					bestcutScores.append(str(max(bslist)))
					bestcutVariances.append(str(self._calcVar(bslist)))
					caseFrequency = caseFrequency / numArg
					caseFrequences.append(str(caseFrequency))
					controlFrequency = controlFrequency / numArg
					controlFrequences.append(str(controlFrequency))
				
				newPos = []  #keeping track of new positions for the tree output
				for result in mappingResults:
					if float(result[3]) == 0.0:
						res = str(float(1.0) / float(numPerm * 2))
					else:
						res = result[3]
					permPvals.append(res)
					mapScores.append(result[2])
					chiPvals.append(result[4])
					pos.append(result[1])
					newPos.append(result[1])
				

				print "Current position: " + pos[ - 1]
				treeList.append('%TREEMARKERS\n')
				for j in range(0, len(newPos)):
					treeList.append("Marker " + str(j) + " at " + newPos[j] + "\n")
				treeList = treeList + newTreeList + newChromasomesUnderCut
		
				#Calculating MAFs
				for j in runset:
					totalCount = 0 
					alleleCount = 0
					for i in range(0, len(accessionsIndices)):
						if snpsd.snps[j][accessionsIndices[i][0]] == 1:
							alleleCount += 1
							totalCount += 1 
						elif snpsd.snps[j][accessionsIndices[i][0]] == 0:
							totalCount += 1
					arf = 0
					if totalCount > 0:
						arf = alleleCount / float(totalCount)
					marfs.append(min(arf, 1.0 - arf))
					mafs.append(min(alleleCount, totalCount - alleleCount))


				#Loop over all batches ends.
			#Loop over runsets ends.

		txtStr = ""  #pos and scores text files
		for j in range(0, len(pos)):
			txtStr = txtStr + str(chromasome + 1) + ", " + str(int(float(pos[j]))) + ", " + mapScores[j] + "\n"
				
		rstr = "bestcutScores <- c(" + (",".join(bestcutScores)) + ")\n" + "bestcutVariances <- c(" + (",".join(bestcutVariances)) + ")\n" + "caseFrequences <- c(" + (",".join(caseFrequences)) + ")\n" + "controlFrequences <- c(" + (",".join(controlFrequences)) + ")\n" + "meanNumCases <- c(" + (",".join(meanNumCases)) + ")\n"
		rstr = rstr + "meanFriFreqG1 <- c(" + (",".join(meanFriFreq[0])) + ")\n" + "meanFriFreqG2 <- c(" + (",".join(meanFriFreq[1])) + ")\n" + "meanFriFreqG3 <- c(" + (",".join(meanFriFreq[2])) + ")\n"

		rstr = rstr + "permPvals <- c(" + (",".join(permPvals)) + ")\n" + "chiPvals <- c(" + (",".join(chiPvals)) + ")\n" + "mapScores <- c(" + (",".join(mapScores)) + ")\n" + "pos <- c(" + (",".join(pos)) + ")\n"
				
		rstr = rstr + "par(mfrow=c(2,1))\n"
		rstr = rstr + 'plot(pos,bestcutScores,pch=20,main="Best-cut Scores")\n'
		rstr = rstr + 'plot(pos,bestcutVariances,pch=20,main="Best-cut Scores Variance Estim.")\n\n'

		rstr = rstr + "par(mfrow=c(3,1))\n"
		rstr = rstr + 'plot(pos,meanNumCases,pch=20,main="Mean size of group 1")\n'
		rstr = rstr + 'plot(pos,caseFrequences,pch=20,main="Mean rank of group 1")\n'
		rstr = rstr + 'plot(pos,controlFrequences,pch=20,main="Mean rank of group 2")\n\n'

		rstr = rstr + "par(mfrow=c(3,1))\n"
		rstr = rstr + 'plot(pos,-log(permPvals)/log(10),pch=20,main="Permutation p-values")\n'
		rstr = rstr + 'plot(pos,mapScores,pch=20,main="ARG mapping scores")\n'
		rstr = rstr + 'plot(pos,-log(chiPvals)/log(10),pch=20,main="Chi-square p-values")\n'
		
		f = open(rfile, "w")
		f.write(rstr)
		f.close()
		os.system("R --vanilla < " + rfile + " > " + outRfile)

		f = open(txtfile, "w")
		f.write(txtStr)
		f.close()

		f = open(treefile, "w")
		f.write(" ".join(treeList))
		f.close()

		print "Done!"
	
 
 

	def getMarginalTrees(self, treeFile, snpsd, position, argBoundaries = None, numArg = 1):
	   
		print "Preparing to run Margarita"

		i = 0

		start = 0
		while snpsd.positions[start] < argBoundaries[0]:
			start = start + 1
		end = start
		while end < len(snpsd.positions) - 1 and snpsd.positions[end] < argBoundaries[1]:
			end = end + 1
		i = start

		runIndices = range(start, end)

		numMarkers = end - start

		# Variables to capture the results
		pos = []
		permPvals = []
		chiPvals = []
		mapScores = []

		treeList = ['NUMARGS: ' + str(numArg) + ' NUMMARKERS: ' + str(numMarkers) + '\n']

		decoder = {0:'0', 1:'1', - 1:'M', 2:'M'}


		print "Running margarita:"
		sys.stdout.write(".")
		sys.stdout.flush()
		
		#Writing margarita file!
		st = ""
		st += str(len(snpsd.accessions)) + " " + str(numMarkers) + "\n"
		for i in runset:
			st += str(snpsd.positions[i]) + "\n"	
		# Construct a dummy phentype.
		st += "0\n" 
		for i in range(1, len(snpsd.accessions)):
			st += "1\n"  #The phenotype
		for i in range(0, len(snpsd.accessions)):
			for j in boundaryIndices:
				st += str(decoder[snpsd.snps[j][i]])   #Decoding the SNP
			st += "\n"

		f = open(self.margFile, "w")
		f.write(st)
		f.close()



		import os
		done = False
		numTries = 0
		while (numTries < 10):				 #Check whether run was ok.  (Some bug in margrita causes it to crash sometimes.)
			numTries = numTries + 1
			os.system(self.margaritaCommand + " " + self.margFile + ' ' + str(numArg) + ' 1 -trees > ' + self.outfile)
			try:
				f = open(outfile, "r")
				lines = f.readlines()
				f.close()
				line = lines[0]
				i = 0
				while not lines[i].startswith("%TREES"):
					i = i + 1
					newTreeList = []
				while not lines[i].startswith("%INTERPRETATION"):
					newTreeList.append(lines[i])
					i = i + 1
				i = i + 2  
				interpretationResults = []
				while not lines[i].startswith("%MAPPING"):
					interpretationResults.append((lines[i].rstrip()).split())
					i = i + 1
				i = i + 1
				labels2 = (lines[i].rstrip()).split()
				i = i + 1
				mappingResults = []
				while i < len(lines):
					mappingResults.append((lines[i].rstrip()).split())
					i = i + 1
					break
			except IndexError:
				print "failed", numTries
			if numTries == 10:
				raise Exception

		bestcutVariances = []
		caseFrequencies = []
		controlFrequencies = []
		newChromasomesUnderCut = ["%UNDERCUTS\n"]
		for j in range(0, len(interpretationResults), numArg):
			bslist = []
			vals = [0] * len(data.individuals)
			caseSum = 0
			caseFrequency = 0
			controlFrequency = 0
			for k in range(0, numArg):
				result = interpretationResults[j + k]
				caseFrequency = caseFrequency + float(result[6])
				controlFrequency = controlFrequency + float(result[7])
				bslist.append(float(result[3]))
				newChromasomesUnderCut.append("MARKER: " + str(j / numArg) + " ARG: " + str(k) + " UNDERCUT: " + result[5] + "\n")
				underCut = []
				for m in range(0, len(result[5])):
					val = int(result[5][m])
					underCut.append(val)
					vals[m] = vals[m] + val
				caseSum = caseSum + sum(underCut)

			meanNumCases.append(str(caseSum / float(numArg)))
			chromasomesUnderCut.append(vals)
			bestcutScores.append(str(max(bslist)))
			bestcutVariances.append(str(self._calcVar(bslist)))
			caseFrequency = caseFrequency / numArg
			caseFrequencies.append(str(caseFrequency))
			controlFrequency = controlFrequency / numArg
			controlFrequencies.append(str(controlFrequency))
				
		newPos = []  #keeping track of new positions for the tree output
		for result in mappingResults:
			if float(result[3]) == 0.0:
				res = str(float(1.0) / float(numPerm * 2))
			else:
				res = result[3]
			permPvals.append(res)
			mapScores.append(result[2])
			chiPvals.append(result[4])
			pos.append(result[1])
			newPos.append(result[1])
				

		print "Current position: " + pos[ - 1]
		treeList.append('%TREEMARKERS\n')
		for j in range(0, len(newPos)):
			treeList.append("Marker " + str(j) + " at " + newPos[j] + "\n")
		treeList = treeList + newTreeList + newChromasomesUnderCut
		
		f = open(treeFile, "w")
		f.write(" ".join(treeList))
		f.close()

		print "Done!"
	
 

	def parseTreeFile(self, treeFileName, rFile, psFile, runNum, argNum, markerNum, accessionMapping = None):
	  
		f = open(treeFileName, "r")
		lines = f.readlines()
		f.close()

		line = (lines[0].rstrip()).split()
		print line
		numRuns = int(line[1])
		numArgs = int(line[3])
		numMarkers = int(line[5])
		runCount = 0
		i = 0
		while i < len(lines) and runCount < runNum:
			i = i + 1
			while not lines[i].startswith(" %TREEMARKERS"):
				i = i + 1
			runCount = runCount + 1
		while not lines[i].startswith(" %TREES"):
			i = i + 1
		argCount = 0
		line = []
		while argCount < argNum:
			i = i + 1
			while not lines[i].startswith(" TREE:"):
				i = i + 1
			line = (lines[i].rstrip()).split()
			argCount = int(line[2]) + 1
		line = (lines[i].rstrip()).split()
		#argCount = int(line[2])+1
		markerCount = int(line[4]) + 1
		while markerCount < markerNum:
			i = i + 1
			while not lines[i].startswith(" TREE:"):
				i = i + 1
			line = (lines[i].rstrip()).split()
			markerCount = int(line[4]) + 1
		i = i + 1
		nodeList = []
		while i < len(lines) and not lines[i].startswith(" TREE:") and not lines[i].startswith(" %TREEMARKERS"):
			line = (lines[i].rstrip()).split()
			nodeList.append([int(line[0]), int(line[1]), int(line[2])])
			i = i + 1
			
		while not lines[i].startswith(" %UNDERCUTS"):  #Find corresponding undercut
			i = i + 1
		i = i + 1
		line = (lines[i].rstrip()).split()
		markerCount = int(line[1]) + 1
		while markerCount < markerNum:
			i = i + numArgs
			line = (lines[i].rstrip()).split()
			markerCount = int(line[1]) + 1
		argCount = int(line[3]) + 1
		while argCount < argNum:
			i = i + 1
			line = (lines[i].rstrip()).split()
			argCount = int(line[3]) + 1
		undercut = []
		line = (lines[i].rstrip()).split()
		for num in line[5]:
			undercut.append(num)
		
		# Now construct distance matrix
		numAccessions = nodeList[0][2]
		distMatrix = []
		for i in range(0, numAccessions):
			l = []
			for j in range(0, numAccessions):
				l.append(0)
			distMatrix.append(l)
		
		# Fill in distance matrix.
		def _kidnap_(node, numMarkers):
			if node < numMarkers:
				return [node]
			else:
				i = node - numMarkers
				node1 = nodeList[i][0]
				node2 = nodeList[i][1]
				return _kidnap_(node1, numMarkers) + _kidnap_(node2, numMarkers)

		distance = 1
		for node in nodeList:
			childs1 = _kidnap_(node[0], numAccessions)
			childs2 = _kidnap_(node[1], numAccessions)
			for child1 in childs1:
				for child2 in childs2:
					if distMatrix[child1][child2] != 0:
						raise Exception("Tree is wrong.. or something")
					distMatrix[child1][child2] = distance
					distMatrix[child2][child1] = distance
			distance = distance + 1
		
		if accessionMapping:
			pass
		#f = open(accessionFile,"r")
		#lines = f.readlines()
		#f.close()

		#rStr = lines[0]
		rStr = rStr + "undercut <- c(" + ", ".join(undercut) + ");\n"
		
		rStr = rStr + "distVec <- c();\n"
		for i in range(0, numAccessions):
			rStr = rStr + "v <- c("
			for j in range(0, numAccessions - 1):
				rStr = rStr + str(distMatrix[i][j]) + ","
			rStr = rStr + str(distMatrix[i][ - 1]) + ");\n" + "distVec <- append(distVec,v);\n"
		#rStr = rStr+"print(distVec);\n"
		rStr = rStr + "distMat<-matrix(distVec, " + str(numAccessions) + ", " + str(numAccessions) + ");\n"
		rStr = rStr + "print(distMat);\n"
		rStr = rStr + 'ag <- hclust(as.dist(distMat),method ="average");\n'
		rStr = rStr + "acc<-acc[ag$order,];\n" + "undercut<-undercut[ag$order];\n" + 'par(mai=c(1,0.75,0.5,0.75));\n'
		rStr = rStr + 'postscript(file="' + psFile + '");\n'						
		rStr = rStr + 'par(las=2);\n' + 'plot(ag,main="",hang =-1 ,sub ="",xlab="",ylab="",axes=F,labels=F);\n'
		rStr = rStr + 'for (j in 1:length(acc[,1])){ \n '
		rStr = rStr + '  mtext ( acc[j,1], side =1, col=undercut[j]*3+1, at=j, line=1, cex=0.8 );\n' + '}\n'
		rStr = rStr + 'for (j in 1:length(acc[,1])){ \n '
		rStr = rStr + '  mtext ( acc[j,2], side =1, col=as.numeric(acc[j,3])+1, at=j, line=-1, cex=0.8 );\n' + '}\n'
		rStr = rStr + "dev.off()\n"
		f = open(rFile, "w")
		f.write(rStr)
		f.close()
		outRfile = rFile + ".out"
		os.system("R --vanilla --file=" + rFile + " > " + outRfile)


		
def _runTest_():
	import dataParsers
	import phenotypeData
	
	#Get phenotype data
	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_transformed_publishable_v2.tsv"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter = '\t')  #Get Phenotype data 

	#Get SNPs data 
	snpsDataFile = "/Network/Data/250k/dataFreeze_080608/250K_f10_080608.csv"
	snpsds = dataParsers.parseCSVData(snpsDataFile) #Get SNPs data 

	psFile = env.homedir + "tmp/tree.ps"
	marg_file = env.homedir + "tmp/test"
	out_file = env.homedir + "tmp/test_out"
	rFile = env.homedir + "tmp/tree_test.r"

	#Run Margarita
	marg = Margarita(marg_file, out_file)
	chr = 4
	snpsd = snpsds[chr - 1].getSnpsData()
	marg.gwaWithTrees(marg_file, snpsd, phed, phenotype = 1, numMarkers = 200, chromosome = chr, boundaries = [200000, 350000], numPerm = 1, cutoff = 16, numArg = 100)
	#(self, id, snpsd, phed, phenotype=0, boundaries = None, numMarkers = 100, numPerm = 500000, cutoff = 16, numArg = 50)

	#which marginal tree
	runNum = 1
	argNum = 1
	markerNum = 1
	
	marg.parseTreeFile(marg_file + ".marg.trees", rFile, psFile, runNum, argNum, markerNum)


import sys, getopt, traceback

def _run_():
	import os
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)	   

	long_options_list = ["id=", "chr=", "numARG=", "numMarkers=", "numPerm=", "smartCutoff=", "BoundaryStart=", "BoundaryEnd=", "binary", "delim=", "missingval=", "withArrayId=", "phenotypeFileType=", "debug", "parallel=", "parallelAll", "help", "scoreFile="]
	
	try:
		opts, args = getopt.getopt(sys.argv[1:], "i:s:c:d:m:a:bh", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
	import tempfile
	tempfile.tempdir = '/tmp'
	(fId, id) = tempfile.mkstemp()
	os.close(fId)		
	scoreFile = None
	chr = None
	numARG = 30
	numMarkers = 100
	numPerm = 0
	smartCutoff = 10
	binary = False
	delim = ","
	missingVal = "NA"
	debug = None
	report = None
	help = 0
	withArrayId = 0
	boundaries = [ - 1, - 1]
	phenotypeFileType = 1
	parallel = None
	parallelAll = False
	snpsDataFile = None

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
		elif opt in ("-i", "--id"):
			id = '/tmp/' + arg
		elif opt in ("-s", "--scoreFile"):
			scoreFile = arg
		elif opt in ("-c", "--chr"):
			chr = int(arg)
		elif opt in ("--numARG"):
			numARG = int(arg)
		elif opt in ("--numMarkers"):
			numMarkers = int(arg)
		elif opt in ("--numPerm"):
			numPerm = int(arg)
		elif opt in ("--BoundaryStart"):
			boundaries[0] = int(arg)
		elif opt in ("--BoundaryEnd"):
			boundaries[1] = int(arg)
		elif opt in ("--smartCutoff"):
			smartCutoff = int(arg)
		elif opt in ("--phenotypeFileType"):
			phenotypeFileType = int(arg)
		elif opt in ("--binary"):
			binary = True
		elif opt in ("--parallel"):
			parallel = arg
		elif opt in ("--parallelAll"):
			parallelAll = True
		elif opt in ("-d", "--delim"):
			delim = arg
		elif opt in ("-m", "--missingval"):
			missingVal = arg	
		elif opt in ("-a", "--withArrayId"):
			withArrayId = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1

	if len(args) < 3 and not parallel:
		if help == 0:
			print "Arguments are missing!!\n"
			print __doc__
		sys.exit(2)

	if boundaries[0] == boundaries[1] and boundaries[0] == - 1:
		boundaries = None

	margFile = id + ".marg"
	outFile = margFile + ".out"		


	def runParallel(phenotypeIndex):
		#Cluster specific parameters
		#margdir = '/home/cmb-01/bvilhjal/Projects/Python-snps/'
		resultDir = env.results_dir #'/home/cmb-01/bvilhjal/results/'
		import phenotypeData
		phed = phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter = '\t')  #Get Phenotype data 
		phenName = phed.getPhenotypeName(phenotypeIndex)
		phenName = phenName.replace("/", "_div_")
		phenName = phenName.replace("*", "_star_")
 
		outFileName = resultDir + "Marg_" + parallel + "_" + phenName
		scoreFile = outFileName + ".score" 

		shstr = """#!/bin/csh
#PBS -l walltime=120:00:00
#PBS -l mem=4g
#PBS -q cmb
"""

		shstr += "#PBS -N M" + phenName + "_" + parallel + "\n"
		#shstr += "(python " + margdir + "margarita.py "
		shstr += "(python " + env.script_dir + "margarita.py "
		if phed.isBinary(phenotypeIndex):
			shstr += " --binary "
		shstr += " -s " + scoreFile
		shstr += " -a " + str(withArrayId) + " "
		shstr += snpsDataFile + " " + phenotypeDataFile + " " + str(phenotypeIndex) + " "
		shstr += "> " + outFileName + ".out) >& " + outFileName + ".err\n"
		
		f = open(parallel + ".sh", 'w')
		f.write(shstr)
		f.close()

		#Execute qsub script
		os.system("qsub " + parallel + ".sh ")
		
	#Nested function ends

	snpsDataFile = args[0]
	phenotypeDataFile = args[1]
	if parallel:  #Running on the cluster..
		if len(args) > 2:
			phenotypeIndex = int(args[2])
			runParallel(phenotypeIndex)
			return
		
		else:
			snpsDataFile = args[0]
			if not parallelAll:
				phenotypeIndex = int(args[1])
				runParallel(phenotypeIndex)
				return

		import phenotypeData
		phed = phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter = '\t')  #Get Phenotype data 
		for phenotypeIndex in phed.phenIds:
			runParallel(phenotypeIndex)
		return

	phenotypeIndex = int(args[2])


	#Print out information about this run...
	print "Preparing a blended margarita...."
	print "Num ARG:", numARG
	print "Num Markers:", numMarkers
	print "Num Permutations:", numPerm
	print "Smart cutoff:", smartCutoff
	print "Binary:", binary
	print "ScoreFile:", scoreFile


	
	import dataParsers, snpsdata, phenotypeData
	#phenotypeFile = "/Users/bjarni/Projects/Python-snps/tinaPhenos_041808.csv"
	if phenotypeFileType == 1:
		phed = phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter = '\t')  #Get Phenotype data 
	elif phenotypeFileType == 2:
		phed = phenotypeData.readPhenotypeFile(phenotypeDataFile, accessionDecoder = dataParsers.accessionName2EcotypeId, type = 2)	

	snpsds = dataParsers.parseCSVData(snpsDataFile, deliminator = delim, missingVal = missingVal, withArrayIds = bool(withArrayId)) #Get SNPs data 



	marg = Margarita(margFile, outFile, numARG, numMarkers, numPerm, smartCutoff)

	if chr:
		snpsd = snpsds[chr - 1].getSnpsData()
		marg.gwa(snpsd, phed, phenotype = phenotypeIndex, boundaries = boundaries, chromosome = chr, binary = binary)
	else:
		scoreStr = ""
		for chr in [0, 1, 2, 3, 4]:
			snpsd = snpsds[chr].getSnpsData()
			(newRStr, newScoreStr, permPvals) = marg.gwa(snpsd, phed, phenotype = phenotypeIndex, boundaries = boundaries, chromosome = chr + 1, binary = binary)
			scoreStr += newScoreStr

		f = open(scoreFile, 'w')
		f.write(scoreStr)
		f.close()


def _generateBinomData_(numIndivid = 10, numDatasets = 100, nseg_sites = 100):
	import random, snpsdata
	snpsds = []
	positions = range(1, nseg_sites + 1)
	for i in range(0, numDatasets):
		snps = []
		for pos in positions:
			snp = []
			for j in range(numIndivid):
				snp.append(round(random.random()))
			snps.append(snp)
		snpsds.append(snpsdata.SnpsData(snps, positions))
	return snpsds

def _testMargarita_(rho = 1, theta = 1, numIndivid = 10, numDatasets = 100, nseg_sites = 100):
	import os, tempfile
	(fId, tempDataFile) = tempfile.mkstemp()
	os.close(fId)

	#"""
	# Generate data using ms.
	msCommand = "~/Projects/programs/ms/ms " + str(numIndivid) + " " + str(numDatasets) + " -s " + str(nseg_sites) + " -t " + str(theta) + " -r " + str(rho) + " " + str(nseg_sites) + " > " + tempDataFile + ""
	print msCommand
	os.system(msCommand)
		

	import dataParsers

	# Parse datasets
	print "parsing ms dataset"
	snpsds = dataParsers.parseMSFile(tempDataFile)
	#"""

	"""
	snpsds = _generateBinomData_()
	"""
	print len(snpsds)
	accessions = range(0, numIndivid)
	for snpsd in snpsds:
		snpsd.accessions = accessions
	print snpsds[0].accessions

	

	import phenotypeData
	#FIXME phenotypeData
	phenValues = []
	for i in range(0, numIndivid / 2):
		phenValues.append([0])
	for i in range(0, numIndivid - numIndivid / 2):
		phenValues.append([1])
		
	print phenValues
	phed = phenotypeData.PhenotypeData(accessions, ["test_phenotype"], phenValues)
	

	(fId, tempMargFile) = tempfile.mkstemp()
	os.close(fId)
	(fId, tempOutFile) = tempfile.mkstemp()
	os.close(fId)

	marg = Margarita(tempMargFile, tempOutFile, 30, nseg_sites, 20000, 20)

	minPvalList = []
	for snpsd in snpsds:	# for all ms datasets 
		#minPval = 1.0
		#while minPval == 1.0: 
		(a, b, permPvals) = marg.gwa(snpsd, phed, phenotype = 0, binary = True)
		pvals = []
		for pval in permPvals:
			pvals.append(float(pval))
		print pvals
		minPval = min(pvals)
		minPvalList.append(minPval)
	print minPvalList
		# run margarita

	# Keep track of num of sign hits for each run..

if __name__ == '__main__':
	_run_()
	#_runTest_()
	#_testMargarita_()
