from snpsdata import *

	
#Statistics........
class Statistic:
	name = None	

	def __init__(self,snpsd):
		raise NotImplementedError(caller + ' must be implemented in subclass')
		
	def saveToFile(self,filename):
		import pickle 
		d = {"stats":self.stats, "lims":self.lims}
		f = open(filename,"w")
		st = pickle.dumps(d)
		f.write(st)
		f.close()
		

	def scaleLimits(self,s):
		for i in xrange(0,len(self.lims)):
				self.lims[i] = self.lims[i]*s

	
	def loadFromFile(self,filename):
		import pickle 
		f = open(filename,"r")
		lines = f.readlines()
		st = ""
		for line in lines:
			st = st+line
		d = pickle.loads(st)
		self.stats = d["stats"]
		self.lims = d["lims"]
		f.close()
	
	def isClose(self,snpsdList): 
		statistic = self.__class__(snpsdList)
		#print statistic
		i=0
		while (i < len(self.lims) and abs(statistic.stats[i]-self.stats[i]) < self.lims[i]):
			i = i+1
		return i == (len(self.lims)-1)



class S(Statistic):
	"""
	Final stat!.
	"""
	name = "S"
		
	#Variance files
	vardir = "/Network/Servers/oak.usc.edu/Volumes/RAID/export/home/bjarni/summer2007/bin/vars/"
	varfiles = ["r2_0-2000.var","r2_2000-10000.var","r2_0-25000.var","D_0-10000.var","D_0-25000.var","P2_5000.var","d_0.2-2.10.var"]
	


	def __init__(self,snpsdList = None):	
		if snpsdList:
			self.stats = [0.0]*7
			self.stats[0] = r2MeanUnbiased(snpsdList,2000)
			self.stats[1] = r2MeanUnbiased(snpsdList,10000,2000)
			self.stats[2] = r2MeanUnbiased(snpsdList,25000)
			
			self.stats[3] = meanDUnbiased(snpsdList,10000)
			self.stats[4] = meanDUnbiased(snpsdList,25000)

			self.stats[5] = patternIIUnbiased(snpsdList,5000)[0]

			self.stats[6] = abs(self.stats[0] - self.stats[1])
			

		#Load variances!!!
		self.means = []
		self.vars = []
		for file in self.varfiles:
			f = open(self.vardir+file,"r")
			import pickle 
			lines = f.readlines()
			st = ""
			for line in lines:
				st = st+line
			d = pickle.loads(st)
			self.means.append(d["means"])
			self.vars.append(d["vars"])
			f.close()

		self.weights = [1,1,1,1,1,0.5,1]
		self.varsOfMeans = [0.00545066049266,0.00663270311036,0.00630405283105,0.00301587115351,0.00571859642099,0.000176389906518,0.00287438181774]
		self.C=0.22 # Constant, accept threshold.



	def saveToFile(self,filename):
		import pickle 
		d = {"stats":self.stats, "C":self.C}
		f = open(filename,"w")
		st = pickle.dumps(d)
		f.write(st)
		f.close()


	def loadFromFile(self,filename):
		import pickle 
		f = open(filename,"r")
		lines = f.readlines()
		st = ""
		for line in lines:
			st = st+line
		d = pickle.loads(st)
		self.stats = d["stats"]
		self.C = d["C"]
		f.close()


	def isClose(self,snpsdList,rhoIndex,tIndex): 
		statistic = self.__class__(snpsdList)
                #print statistic
		d=0
		for i in range(0,len(self.stats)):
			var = self.vars[i][0][rhoIndex][tIndex]+self.varsOfMeans[i]
			d = d + ((statistic.stats[i]-self.stats[i])*(statistic.stats[i]-self.stats[i])*(self.weights[i]))/var
			self.weights[i] = (self.weights[i])**(1.0/2.0)
		return (d**(1.0/2.0)) < self.C*sum(self.weights)

	def __str__(self):
		return "r2 0-2000: "+str(self.stats[0])+"     r2 2000-10000: "+str(self.stats[1])+"     r2 0-25000: "+str(self.stats[2])+"     D 0-10000: "+str(self.stats[3])+"     D 0-25000: "+str(self.stats[4])+"     P_II 0-5000: "+str(self.stats[5])+"     r2 0-2000 - 2000-10000 : "+str(self.stats[6])
	

	
class S2(Statistic):
	"""
	11/16

	Long range stat.
	"""
	name = "S2"
		
	#Variance files
	vardir = "/data/varfiles/"
	varfiles = ["r2_2000-5000.var","r2_5000-25000.var","r2MAF10-25_2000-25000.var","r2MAF25-50_2000-25000.var","D_2000-10000.var","D_10000-25000.var","d_2.5-5.25.var"]
	


	def __init__(self,snpsdList = None):	
		if snpsdList:
			self.stats = [0.0]*6
			self.stats[0] = r2MeanUnbiased(snpsdList,5000,2000)
			self.stats[1] = r2MeanUnbiased(snpsdList,25000,5000)
			
			self.stats[2] = r2MeanUnbiasedMAF(snpsdList,25000,2000,mafs=[0.1,0.25])
			self.stats[3] = r2MeanUnbiasedMAF(snpsdList,25000,2000,mafs=[0.25,0.5])

			self.stats[4] = meanDUnbiased(snpsdList,10000,2000)
			self.stats[5] = meanDUnbiased(snpsdList,25000,2000)

			self.stats[6] = abs(self.stats[0] - self.stats[1])
			

		#Load variances!!!
		self.means = []
		self.vars = []
		for file in self.varfiles:
			f = open(self.vardir+file,"r")
			import pickle 
			lines = f.readlines()
			st = ""
			for line in lines:
				st = st+line
			d = pickle.loads(st)
			self.means.append(d["means"])
			self.vars.append(d["vars"])
			f.close()

		self.weights = [0.5,0.5,0.5,0.5,0.5,0.5,1]
		self.varsOfMeans = [0.00580753875156,0.00468511496677,0.00325359397664,0.0128751955649,0.0128751955649,0.00346848867541,0.00670791114706,0.001198911309]
		self.C=0.2 # Constant, accept threshold.



	def saveToFile(self,filename):
		import pickle 
		d = {"stats":self.stats, "C":self.C}
		f = open(filename,"w")
		st = pickle.dumps(d)
		f.write(st)
		f.close()


	def loadFromFile(self,filename):
		import pickle 
		f = open(filename,"r")
		lines = f.readlines()
		st = ""
		for line in lines:
			st = st+line
		d = pickle.loads(st)
		self.stats = d["stats"]
		self.C = d["C"]
		f.close()


	def isClose(self,snpsdList,rhoIndex,tIndex): 
		statistic = self.__class__(snpsdList)
                #print statistic
		d=0
		for i in range(0,len(self.stats)):
			var = self.vars[i][0][rhoIndex][tIndex]+self.varsOfMeans[i]
			d = d + ((statistic.stats[i]-self.stats[i])*(statistic.stats[i]-self.stats[i])*(self.weights[i]))/var
			self.weights[i] = (self.weights[i])**(1.0/2.0)
		return (d**(1.0/2.0)) < self.C*sum(self.weights)

	def __str__(self):
		return "r2 2000-5000: "+str(self.stats[0])+"     r2 5000-25000: "+str(self.stats[1])+"     r2 10%<MAF<25% 2000-25000: "+str(self.stats[2])+"     r2 MAF>25% 2000-25000: "+str(self.stats[2])+"     D 2000-10000: "+str(self.stats[3])+"     D 2000-25000: "+str(self.stats[4])+"     r2 2000-5000 - 5000-25000 : "+str(self.stats[5])
	



""" -------------------------------------------------- """



def patternInIInew(snpsds,windowSize,threshold=0.1):
	totcount = 0
	pIcount = 0
	pIIcount = 0
	for snpsd in snpsds:
		snpsd.snpsFilterRare(threshold)
		if len(snpsd.snps)>2:
			freqs =	snpsd.calcFreqs(windowSize)	
			for i in xrange(0,len(freqs)-1):
				if len(freqs[i])>1:
					for j in xrange(0,len(freqs[i])-1):
						d12 = D(freqs[i][j])  #snps1,snps2
						for k in xrange(j+1,len(freqs[i])):
							d13 = D(freqs[i][k]) #snps1,snps3
							d23 = D(freqs[i+j+1][k-j-1])  #snps2,snps3
							totcount = totcount +1
							if (d12<d13):
								pIcount = pIcount + 1
								if (d23<d13):
									pIIcount = pIIcount + 1
							else:
								if (d23<d13):
									pIcount = pIcount + 1
			del freqs
			del snpsd.freqs
		"""
		print "pIcount", pIcount
		print "pIIcount", pIIcount
		print "totcount", totcount
		"""
	if totcount > 0:
		return (pIcount/float(totcount),pIIcount/float(totcount),totcount)
	else:
		return (0,0,0)



def patternII(snpsds,windowSize,threshold=0.1):
	totcount = 0
	pIIcount = 0
	for snpsd in snpsds:
		snpsd.snpsFilterRare(threshold)
		if len(snpsd.snps)>2:
			freqs =	snpsd.calcFreqs(windowSize)	
			for i in xrange(0,len(freqs)-1):
				if len(freqs[i])>1:
					for j in xrange(0,len(freqs[i])-1):
						d12 = D(freqs[i][j])  #snps1,snps2
						for k in xrange(j+1,len(freqs[i])):
							d13 = D(freqs[i][k]) #snps1,snps3
							d23 = D(freqs[i+j+1][k-j-1])  #snps2,snps3
							totcount = totcount +1
							if (d12<d13 and d23<d13):
								pIIcount = pIIcount + 1
			del freqs
			del snpsd.freqs
		
	if totcount>0:
		return (pIIcount/float(totcount),totcount)
	else:
		return(0,0)

def patternIIUnbiased(snpsds,windowSize,threshold=0.1):
	totcount = 0
	pIIcount = 0
	for snpsd in snpsds:
		snpsd.snpsFilterRare(threshold)
		if len(snpsd.snps)>2:
			freqs =	snpsd.calcFreqsUnbiased(windowSize)	
			for i in xrange(0,len(freqs)-1):
				if len(freqs[i])>1:
					for j in xrange(0,len(freqs[i])-1):
						if (i+j+1 < len(freqs)):
							d12 = D(freqs[i][j])  #snps1,snps2
							for k in xrange(j+1,len(freqs[i])):
								d13 = D(freqs[i][k]) #snps1,snps3
								d23 = D(freqs[i+j+1][k-j-1])  #snps2,snps3
								totcount = totcount +1
								if (d12<d13 and d23<d13):
									pIIcount = pIIcount + 1
			del freqs
			del snpsd.freqs
		
	if totcount>0:
		return (pIIcount/float(totcount),totcount)
	else:
		return(0,0)


def patternIII(snpsds,windowSize):
	totcount = 0
	pIIIcount = 0
	for snpsd in snpsds:
		if len(snpsd.snps)>2:
			freqs =	snpsd.calcFreqs(windowSize)	
			for i in xrange(0,len(freqs)-1):
				if len(freqs[i])>1:
					for j in xrange(0,len(freqs[i])-1):
						d12 = D(freqs[i][j])  #snps1,snps2
						print "D12 =",d12
						for k in xrange(j+1,len(freqs[i])):
							d23 = D(freqs[i+j+1][k-j-1])  #snps2,snps3
							print "D23 =",d23
							totcount = totcount +1
							if (d12<0.5 and d23<0.5):
								pIIIcount = pIIIcount + 1
			del freqs
			del snpsd.freqs
	if totcount>0:
		return (pIIIcount/float(totcount),totcount)
	else:
		return(0,0)


def patternIVnew(snpsds,windowSize,threshold=0.1):
	totcount = 0
	pIVcount = 0
	for snpsd in snpsds:
		if len(snpsd.snps)!=0 and len(snpsd.snps[0])!=0:
			snpsd.snpsFilterRare(threshold)
			if len(snpsd.snps)>1:
				freqs =	snpsd.calcFreqs(windowSize)	
				for i in xrange(0,len(freqs)-1):
					for j in xrange(0,len(freqs[i])):
						d12 = D(freqs[i][j])  #snps1,snps2
						totcount = totcount +1
						if (d12<1.0):
							pIVcount = pIVcount + 1
				del freqs
				del snpsd.freqs
	if totcount>0:
		return (pIVcount/float(totcount),totcount)
	else:
		return (0,0)


def meanD(snpsdList,windowSize, innerWindowSize=0, threshold=0.1):
	count =	0
	h=0
	for snpsd in snpsdList:
		if len(snpsd.snps)!=0 and len(snpsd.snps[0])!=0:
			snpsd.snpsFilterRare(threshold)
			freqs =	snpsd.calcFreqs(windowSize,innerWindowSize)	
			for i in xrange(0,len(freqs)):
				for j in xrange(0,len(freqs[i])):
					h = D(freqs[i][j]) + h
					count =	count+1
			snpsd.freqs=[]
			del freqs
	if count!= 0:
		h = h/float(count)
	else:
		h=0
	return h


def meanDUnbiased(snpsdList,windowSize, innerWindowSize=0, threshold=0.1):
	count =	0
	h=0
	for snpsd in snpsdList:
		if len(snpsd.snps)!=0 and len(snpsd.snps[0])!=0:
			snpsd.snpsFilterRare(threshold)
			freqs =	snpsd.calcFreqsUnbiased(windowSize,innerWindowSize)	
			for i in xrange(0,len(freqs)):
				for j in xrange(0,len(freqs[i])):
					h = D(freqs[i][j]) + h
					count =	count+1
			snpsd.freqs=[]
			del freqs
	if count!= 0:
		h = h/float(count)
	else:
		h=0
	return h

def calcCAF(snpsdList):
	count = 0
	tot = 0
	for snpsd in snpsdList:
		tot = tot + len(snpsd.positions)
		snpsdd = snpsd.snpsFilter()
		count = count + len(snpsdd[3].positions) + len(snpsdd[4].positions) + len(snpsdd[5].positions) + len(snpsdd[6].positions)
	return count/float(tot)
		
def calcMMAF(snpsdList, threshold=0.15):
	c = 0
	mmaf = 0
	for snpsd in snpsdList:
		snpsd.snpsFilterRare(threshold)
		mmaf = mmaf + snpsd.meanAF()*len(snpsd.snps)
		c = c + len(snpsd.snps)
	if c:
		return mmaf/c    
	else:
		return 0


def r2Frac(snpsdList,threshold,windowSize, innerWindowSize=0):
	count =	0
	h=0
	for snpsd in snpsdList:
		freqs =	snpsd.calcFreqs(windowSize,innerWindowSize)	
		for i in xrange(0,len(freqs)):
			for j in xrange(0,len(freqs[i])):
				r = r2(freqs[i][j])
				count =	count+1
				if r > threshold: 
					h = h +	1
		snpsd.freqs=[]
		del freqs
	if count!= 0:
		h = h/float(count)
	else:
		h=0
	return h
		

def r2MeanNew(snpsdList,windowSize, innerWindowSize=0, threshold=0.1):
	count =	0
	h=0
	for snpsd in snpsdList:
		if len(snpsd.snps)!=0 and len(snpsd.snps[0])!=0:
			snpsd.snpsFilterRare(threshold)
			freqs =	snpsd.calcFreqs(windowSize,innerWindowSize)	
			for i in xrange(0,len(freqs)):
				for j in xrange(0,len(freqs[i])):
					h = r2(freqs[i][j]) + h
					count =	count+1
			snpsd.freqs=[]
			del freqs
	if count!= 0:
		h = h/float(count)
	else:
		h=0
	return h


def r2MeanUnbiased(snpsdList,windowSize, innerWindowSize=0, threshold=0.1):
	count =	0
	h=0
	for snpsd in snpsdList:
		if len(snpsd.snps)!=0 and len(snpsd.snps[0])!=0:
			snpsd.snpsFilterRare(threshold)
			freqs =	snpsd.calcFreqsUnbiased(windowSize,innerWindowSize)	
			for i in xrange(0,len(freqs)):
				for j in xrange(0,len(freqs[i])):
					h = r2(freqs[i][j]) + h
					count =	count+1
			snpsd.freqs=[]
			del freqs
	if count!= 0:
		h = h/float(count)
	else:
		h=0
	return h


def r2MeanUnbiasedMAF(snpsdList,windowSize, innerWindowSize=0, mafs=[0.3,0.5]):
	count =	0
	h=0
	for snpsd in snpsdList:
		if len(snpsd.snps)!=0 and len(snpsd.snps[0])!=0:
			snpsd.snpsFilterMAF(mafs)
			freqs =	snpsd.calcFreqsUnbiased(windowSize,innerWindowSize)	
			for i in xrange(0,len(freqs)):
				for j in xrange(0,len(freqs[i])):
					h = r2(freqs[i][j]) + h
					count =	count+1
			snpsd.freqs=[]
			del freqs
	if count!= 0:
		h = h/float(count)
	else:
		h=0
	return h

def r2MeanSD(snpsdList,windowSize, innerWindowSize=0):
	""" Calculates standard deviation as well.. """
	count =	0
	h=0
	v=0
	for snpsd in snpsdList:
		freqs =	snpsd.calcFreqs(windowSize,innerWindowSize)	
		for i in xrange(0,len(freqs)):
			for j in xrange(0,len(freqs[i])):
				h = r2(freqs[i][j]) + h
				count =	count+1
	if count!= 0:
		h = h/float(count)
		for snpsd in snpsdList:
			for i in xrange(0,len(snpsd.freqs)):
				for j in xrange(0,len(snpsd.freqs[i])):
					l = (r2(freqs[i][j]) - h)
					v = v+l*l
			snpsd.freqs=[]
			del freqs
		if count!=1:
			v=v/(count-1)
		else:
			v=1
	else:
		h=0
		v=1

	return (h,v)

def meanEHH(snpsdList,windowSize, innerWindowSize=0, threshold=0.1):
	ehhcount = 0
	totehh = 0
	i = 0
	for snpsd in snpsdList:
		i = i + 1
		snpsd.snpsFilterRare(threshold)
		res = snpsd.totalEHH(windowSize,innerWindowSize)
		totehh = totehh + res[0]
		ehhcount = ehhcount + res[1]
		if ehhcount == 0 and totehh == 0:			
			ehhcount=1
	return totehh/ehhcount



#-----------------------------------------------------------------------------------#
# Statistics testing

class Stat:
	def __init__(self):
		pass
		
	def calcStat(self,snpsdList = None):
		raise NotImplementedError(caller + ' must be implemented in subclass')




class T1new(Stat):
	name = "Pattern 1 w MAF>0.1"
	def __init__(self,windowSize=5000):
		self.windowSize=windowSize
		self.name = self.name+" ("+str(windowSize)+" bases)"

	def calcStat(self,snpsdList):
		return patternInIInew(snpsdList,self.windowSize)[0]

class T2new(Stat):
	name = "Pattern 2 w MAF>0.1"
	def __init__(self,windowSize=5000):
		self.windowSize=windowSize
		self.name = self.name+" ("+str(windowSize)+" bases)"

	def calcStat(self,snpsdList):
		return patternInIInew(snpsdList,self.windowSize)[1]

class T3(Stat):
	name = "Pattern 2 w MAF>0.1"
	def __init__(self,windowSize=10000):
		self.windowSize=windowSize
		self.name = self.name+" ("+str(windowSize)+" bases)"

	def calcStat(self,snpsdList):
		return patternIIUnbiased(snpsdList,self.windowSize)[0]

class T4new(Stat):
	name = "Pattern 4 w MAF>0.1"
	def __init__(self,windowSize=50000):
		self.windowSize=windowSize
		self.name = self.name+" ("+str(windowSize)+" bases)"

	def calcStat(self,snpsdList):
		return patternIVnew(snpsdList,self.windowSize)[0]


class T5(Stat):
	name = "Mean r2 w MAF>0.1"
	def __init__(self,windowSize=2000,innerWindowSize=0):
		self.windowSize=windowSize
		self.innerWindowSize=innerWindowSize
		self.name = self.name+" (Window size: "+str(windowSize)+" bases.  Inner window size: "+str(innerWindowSize)+" bases)"

	def calcStat(self,snpsdList):
		return r2MeanUnbiased(snpsdList,self.windowSize,self.innerWindowSize)


class T6(Stat):
	name = "Difference between 2000-5000 and 5000-25000 bases"

	def calcStat(self,snpsdList):
		s1 = r2MeanUnbiased(snpsdList,5000,2000)
		s2 = r2MeanUnbiased(snpsdList,25000,5000)
		return s1-s2     #Dif

class T7(Stat):
	name = "Common allele freq."

	def calcStat(self,snpsdList):
		return calcCAF(snpsdList)    


class T8(Stat):
	name = "Mean allele freq."

	def calcStat(self,snpsdList):
		return calcMMAF(snpsdList)    

class T9(Stat):
	name = "Mean D w MAF>0.1"
	def __init__(self,windowSize=2000,innerWindowSize=0):
		self.windowSize=windowSize
		self.innerWindowSize=innerWindowSize
		self.name = self.name+" (Window size: "+str(windowSize)+" bases.  Inner window size: "+str(innerWindowSize)+" bases)"

	def calcStat(self,snpsdList):
		return meanDUnbiased(snpsdList,self.windowSize,self.innerWindowSize)

class T10(Stat):
	name = "Mean r2 w MAFs"
	def __init__(self,windowSize=25000,innerWindowSize=2000,mafs=[0.1,0.25]):
		self.windowSize=windowSize
		self.innerWindowSize=innerWindowSize
		self.mafs = mafs
		self.name = self.name+" (Window sizes: ["+str(innerWindowSize)+","+str(windowSize)+"]; MAFs: "+str(mafs)+".)"

	def calcStat(self,snpsdList):
		return r2MeanUnbiasedMAF(snpsdList,self.windowSize,self.innerWindowSize,self.mafs)

class T11(Stat):
	name = "Mean EHH w MAF>0.1"
	def __init__(self,windowSize=20000,innerWindowSize=2000):
		self.windowSize=windowSize
		self.innerWindowSize=innerWindowSize
		self.name = self.name+" (Window size: "+str(windowSize)+" bases.  Inner window size: "+str(innerWindowSize)+" bases)"

	def calcStat(self,snpsdList):
		return meanEHH(snpsdList,self.windowSize,self.innerWindowSize)


class StatTester:
	#database = "/data/mcmcDataPopExp_bnum50000_self0.95_seg100_num1000/"+"assess_"
	database = "/data/mcmcData_bnum50000_self0.95_seg100_num1000/"+"assess_"

	def __init__(self, stat, thetas=None, rhos=None, ts=None):
		self.stat=stat
		self.rhos = rhos
		self.ts = ts
		self.thetas = thetas
		
	def testStat(self, filterProb=1, sampleNum=500):

		values = []
		for theta in self.thetas:				
			l1 = []
			for rho in self.rhos:
				l2 = []
				for t in self.ts:	
					l2.append([])
				l1.append(l2)
			values.append(l1)

		for j in xrange(0,len(self.thetas)):				
			theta = self.thetas[j]
			v1 = values[j]
			for k in xrange(0,len(self.rhos)):
				v2 = v1[k]
				r = self.rhos[k]
				print 100*k/len(self.rhos),"%"
				for l in xrange(0,len(self.ts)):
					t = self.ts[l]
					filename = self.database+"t"+str(t)+"_r"+str(r)+"_th"+str(theta)+"_dat"
					#print filename
					snpsdList = parseMSDataFilter(filename, baseScale=50000,sampleNum=sampleNum,filterProb=filterProb)
					for snpsd in snpsdList:
						val = self.stat.calcStat([snpsd])
						#if val >1: 
						#	print val
						v2[l].append(val)
						
		self.values = values


	def rstr_rhos_ts(self,lab="",thIndex=0):
		lab = lab + "\n(theta="+str(self.thetas[thIndex])+")"
		st = "v <- c("
		values = self.values[thIndex]
		for i in xrange(0,len(self.rhos)):
			for j in xrange(0,len(self.ts)):
				st = st+str(sum(values[i][j])/float(len(values[i][j])))+","
		st = st[0:-1]+")\nx <- c("
		for t in self.ts:
			if t>0:
				st = st+str(t)+","
			else:
				st = st+"1,"
		st = st[0:-1]+")\ny <- c("
		for rho in self.rhos:
			st = st+str(rho)+","			
		st = st[0:-1]+")\nm <- matrix(v,"+str(len(self.ts))+","+str(len(self.rhos))+')\nimage.plot(x,y,m,xlab="Ts",ylab="Rho",main="'+lab+'")\n'
		return st
		

	def rstr_rhos_thetas(self,lab="",tsIndex=0):
		lab = lab + "\n(ts="+str(self.ts[tsIndex])+")"
		st = "v <- c("
		values = self.values
		for i in xrange(0,len(self.thetas)):
			for j in xrange(0,len(self.rhos)):
				st = st+str(sum(values[i][j][tsIndex])/float(len(values[i][j][tsIndex])))+","
		st = st[0:-1]+")\nx <- c("
		for t in self.thetas:
			st = st+str(t)+","
		st = st[0:-1]+")\ny <- c("
		for rho in self.rhos:
			st = st+str(rho)+","			
		st = st[0:-1]+")\nm <- matrix(v,"+str(len(self.rhos))+","+str(len(self.thetas))+')\nimage.plot(y,x,m,ylab="Theta",xlab="Rho",main="'+lab+'")\n'
		return st
		
	def snr(self):
		""" Signal to noise ratio. """
		values = self.values
		vars = []
		means = []
		for i in xrange(0,len(self.thetas)):
			for j in xrange(0,len(self.rhos)):
				for k in xrange(0,len(self.ts)):
					vals = values[i][j][k]
					mean = sum(vals)/float(len(vals))
					means.append(mean)
					s = 0
					for l in xrange(0,len(vals)):
						s = s+(vals[l]-mean)*(vals[l]-mean)
					vars.append(s/float(len(vals)-1))
		print "Mean variance ", sum(vars)/len(vars)


		mean = sum(means)/len(means)
		s = 0
		for i in xrange(0,len(means)):
			s = s+(means[i]-mean)*(means[i]-mean)

		ovar = 0
		c = 0
		for i in xrange(0,len(self.thetas)):
			for j in xrange(0,len(self.rhos)):
				for k in xrange(0,len(self.ts)):
					vals = values[i][j][k]
					c = c+len(vals)
					for val in vals:
						ovar = ovar + (val-mean)*(val-mean)
		ovar = ovar/(c-1)

		print "Variance of the means ", s/float(len(means)-1)
		print "SNR =",(s/float(len(means)-1))/(sum(vars)/len(vars))
		print "Overall mean =",mean
		print "Overall variance =",ovar

	def saveVarFile(self,filename):
		values = self.values
		vars = []
		means = []
		for i in xrange(0,len(self.thetas)):
			l1m = []
			l1v = []
			for j in xrange(0,len(self.rhos)):
				l2m=[]
				l2v=[]
				for k in xrange(0,len(self.ts)):
					vals = values[i][j][k]
					mean = sum(vals)/float(len(vals))
					s = 0
					for l in xrange(0,len(vals)):
						s = s+(vals[l]-mean)*(vals[l]-mean)
					l2v.append(s/float(len(vals)-1))
					l2m.append(mean)
				l1m.append(l2m)
				l1v.append(l2v)
			means.append(l1m)
			vars.append(l1v)

		#Writing R code
		lab = "\n(theta="+str(self.thetas[0])+")"
		st = "v <- c("
		values = vars[0]
		for i in xrange(0,len(self.rhos)):
			for j in xrange(0,len(self.ts)):
				st = st+str(values[i][j])+","
		st = st[0:-1]+")\nx <- c("
		for t in self.ts:
			if t>0:
				st = st+str(t)+","
			else:
				st = st+"1,"
		st = st[0:-1]+")\ny <- c("
		for rho in self.rhos:
			st = st+str(rho)+","			
		st = st[0:-1]+")\nm <- matrix(v,"+str(len(self.ts))+","+str(len(self.rhos))+')\nimage.plot(x,y,m,xlab="Ts",ylab="Rho",main="'+lab+'")\n'
		print "Variance!:\n"+st


		import pickle 
		d = {"means":means, "vars":vars}
		f = open(filename,"w")
		st = pickle.dumps(d)
		f.write(st)
		f.close()


def runTest():
	import time
	rt = time.time()
	print "\npar(mfrow=c(2,2))"
	sampleNum=500
 	
	rhos=[]
	for i in range(1,25):
		rhos.append(i*2.0)
	ts=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,-1]
	thetas=[250.0]#,150.0,250.0,350.0,450.0]	
	

	t4 = T10()
	tester = StatTester(t4, thetas=thetas, rhos=rhos, ts=ts)
	tester.testStat(filterProb=1.1, sampleNum=sampleNum)
	print ""
	print tester.rstr_rhos_ts(lab = t4.name)
	print tester.rstr_rhos_thetas(lab = t4.name, tsIndex=1)
	tester.snr()
	tester.saveVarFile("r2MAF10-25_2000-25000.var")


	dif = int(time.time() - rt)
	print "It took "+str(dif/60)+" min. and "+str(dif%60)+" sec. to run."

if __name__ == "__main__":
	runTest()
#-----------------------------------------------------------------------------------#



#
#
#
