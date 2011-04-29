"""
Contains statistics object(s) which can be applied in summary statistics calculations.
"""
import util, math
class _Statistic_:
    name = None	
    stats = None # A list which contains the actual statistic values.
    lims = None  # The limits within which it is declared close.
    limit = None  # A single global limit within which a pair of statistics are declared similar.
    scales = None # Scales by which the statistics should be scaled
    
    def __init__(self,lims=None,limit=None,scales=None):
        self.limit = 0.2
        self.lims = [self.limit]
        self.scales = [1] 
        if lims and lims[0]:
            self.lims = lims
        if limit:
            self.limit = limit
        if scales:
            self.scales = scales 
            

    def saveToFile(self,filename):
        import pickle 
        d = {"stats":self.stats, "lims":self.lims}
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
        self.lims = d["lims"]
        f.close()
	
    def setObsData(self,obsData):
        self.stats = self.calcStat(obsData) 

    def calcStat(self,data):  #Calculates the statistics
        raise NotImplementedError(caller + ' must be implemented in subclass')            

    def scaleByVariance(self,datas):
        statisticsList = []
        for data in datas:
            statisticsList.append(self.calcStat(data))

        statLists = util.transposeDoubleLists(statisticsList)
        vars = [1]*len(statLists)
        self.scales = [0]*len(statLists)
        for i in range(0,len(statLists)):  #Now iterating over statistics (instead of datas like before transposition)
            statList = statLists[i]
            #print statList
            if len(statList) >0:
                vars[i] = util.calcVar(statList)
            self.scales[i] = 1.0/math.sqrt(vars[i])
                
            
    def isCloseStandard(self,data):
        """
        Uses the list of limits (self.lims).
        """
        statistic = self.calcStat(data)
	#print statistic
        i=0
        while (i < len(self.lims) and abs(statistic[i]-self.stats[i])*self.scales[i] < self.lims[i]):
            print statistic[i], self.stats[i], self.scales[i], self.lims[i]
            i = i+1
        return i == (len(self.lims))

    def euclidDist(self,data): 
        """
        Uses the single limit value (self.limit).  Returns the actual dist. value
        """
        statistic = self.calcStat(data)
	#print statistic
        i=0
        dist = 0
        for i in xrange(0,len(self.stats)):
            d = (statistic[i]-self.stats[i])*self.scales[i]
            dist += d*d
        dist = math.sqrt(dist)
        return dist


    def isCloseEuclid(self,data): 
        """
        Uses the single limit value (self.limit).
        """
        dist = self.euclidDist(data) 
        return dist < self.limit


class _StatisticsSet_(_Statistic_):
    _statistics_ = []  #The list of statistics.
    
    def __init__(self,statList,limit=None):
        """
        Inheriting classes can call _StatisticSet_.__init__(self) for convenience.
        """
        self.lims = []
        self.scales = []
        self._statistics_ = statList
        self.limit = 0
        for stat in self._statistics_:
            self.limit += stat.limit*stat.limit
            if stat.lims:
                for lim in stat.lims:
                    self.lims.append(lim)
            if stat.scales:
                for scale in stat.scales:
                    self.scales.append(scale)
        if limit:
            self.limit = limit
        else:
            self.limit = math.sqrt(self.limit)
            

    def calcStat(self,data):  #Calculates the statistics
        stats = []
        for statistic in self._statistics_:
            svals = statistic.calcStat(data)
            for v in svals:
                stats.append(v)
        return(stats)
            
class NumSegSites(_Statistic_):

    def __init__(self,lim=None,scale=None):
        self.name="Num_seg_sites"
        _Statistic_.__init__(self,lims=[lim],limit=lim,scales=[scale])  


    def calcStat(self,data): #Calculates the statistics
        return [len(data.positions)]  #Num of seg. sites..


class StupidStat(_Statistic_):

    def __init__(self,lim=None,scale=None):
        self.name="Stupid_stat"
        _Statistic_.__init__(self,lims=[lim],limit=lim,scales=[scale])  

    def calcStat(self,data): #Calculates the statistics
        import random
        return [25.0*random.random()] 


class PairwiseDiff(_Statistic_):

    def __init__(self,lim=None,scale=None):
        self.name="pairwise_diff"
        _Statistic_.__init__(self,lims=[lim],limit=lim,scales=[scale])  

    def calcStat(self,data): #Calculates the statistics
        meanDiff = 0
        if len(data.snps) > 0:
            transposedSnps = util.transposeDoubleLists(data.snps)
            diffCounts = []
            for i in xrange(0,len(data.snps[0])):
                for j in xrange(0,len(data.snps[0])-i-1):
                    hapl1 = transposedSnps[i]
                    hapl2 = transposedSnps[j]
                    diffCount = 0
                    for k in xrange(0,len(hapl1)):
                        if hapl1[k]!=hapl2[k]:
                            diffCount += 1.0
                    diffCounts.append(diffCount)
            meanDiff = sum(diffCounts)/len(diffCounts)
        return [meanDiff] 


class MeanR2(_Statistic_):

    def __init__(self,lim=None,scale=None):
        self.name="Mean_r2"
        _Statistic_.__init__(self,lims=[lim],limit=lim,scales=[scale])  
        self.limit = 0.4
        self.lims = [self.limit]


    def calcStat(self,data): #Calculates the statistics
        meanR2 = 25*r2MeanUnbiased(data,windowSize=0.5,threshold=0)
        return [meanR2] 



class NumHaplotypes(_Statistic_):

    def __init__(self,lim=None,scale=None):
        self.name="Num_haplotypes"
        _Statistic_.__init__(self,lims=[lim],limit=lim,scales=[scale])  

    def calcStat(self,data): #Calculates the statistics
        numHaplotypes = 0
        if len(data.snps)>0:
            transposedSnps = util.transposeDoubleLists(data.snps)
            snpsStringList = []
            for hapl in transposedSnps:
                snpsStringList.append(str(hapl))
            numHaplotypes = len(set(snpsStringList))

        return [numHaplotypes] 


class MaxHaplotypeFreq(_Statistic_):

    def __init__(self,lim=None,scale=None):
        self.name="Max_haplotype_freq"
        _Statistic_.__init__(self,lims=[lim],limit=lim,scales=[scale])  

    def calcStat(self,data): #Calculates the statistics
        maxHaplFreq = 0
        if len(data.snps)>0:
            transposedSnps = util.transposeDoubleLists(data.snps)
            snpsStringList = []
            for hapl in transposedSnps:
                snpsStringList.append(str(hapl))
            haplotypes = set(snpsStringList)
            maxCount = 0
            for hapl in haplotypes:
                c = snpsStringList.count(hapl)
                if c > maxCount:
                    maxCount = c
            maxHaplFreq = float(maxCount)/len(snpsStringList)
        return [maxHaplFreq] 


class NumSingletonHaplotypes(_Statistic_):

    def __init__(self,lim=None,scale=None):
        self.name="Num_singleton_haplotypes"
        _Statistic_.__init__(self,lims=[lim],limit=lim,scales=[scale])  

    def calcStat(self,data): #Calculates the statistics
        singletonCount = 0
        if len(data.snps)>0:
            transposedSnps = util.transposeDoubleLists(data.snps)
            snpsStringList = []
            for hapl in transposedSnps:
                snpsStringList.append(str(hapl))
            haplotypes = set(snpsStringList)
            for hapl in haplotypes:
                if snpsStringList.count(hapl) == 1:
                    singletonCount += 1
        return [singletonCount] 


class StatComb1(_StatisticsSet_):

    def __init__(self,lim=None,limit=None,scale=None):
        statList = [NumSegSites(),PairwiseDiff()]
        _StatisticsSet_.__init__(self,statList=statList,limit=1.0)
        self.name="Num_haplotypes_and_mean_num_pairwise_diff"
        self.lims=[0.1,0.1]

class StatComb2(_StatisticsSet_):

    def __init__(self,lim=None,limit=None,scale=None):
        statList = [NumSegSites(),PairwiseDiff(),NumHaplotypes()]
        _StatisticsSet_.__init__(self,statList=statList,limit=1.0)
        self.name="Num_haplotypes_and_mean_num_pairwise_diff_and_num_haplotypes"
        self.lims=[0.1,0.1,0.1]

class StatComb3(_StatisticsSet_):

    def __init__(self,lim=None,limit=None,scale=None):
        statList = [NumHaplotypes(),MeanR2()]
        _StatisticsSet_.__init__(self,statList=statList,limit=1.0)
        self.name="Num_haplotypes_and_mean_r2"
        self.lims=[0.1,0.1,0.1]



def r2MeanUnbiased(snpsd,windowSize, innerWindowSize=0, threshold=0.1):
    count = 0
    h=0
    if len(snpsd.snps)!=0 and len(snpsd.snps[0])!=0:
        snpsd.snpsFilterRare(threshold)
        freqs =	snpsd.calcFreqsUnbiased(windowSize,innerWindowSize)	
        #print freqs 
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

def r2(freqs):
	f1 = freqs[1]+freqs[3]
	f2 = freqs[2]+freqs[3]
	D = freqs[3]-f1*f2
	divisor = f1*f2*(1-f1)*(1-f2)
	if divisor != 0:
		return D*D/divisor
	else:
		return -1

