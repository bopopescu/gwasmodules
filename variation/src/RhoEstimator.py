#!/usr/bin/env python2.5
"""
Usage: RhoEstimator.py [OPTIONS] OBS_DATA_FILE

Option:

	--simulatedDataFile=...     A file with simulated ms data, to be used for the rejection scheme
	--testStatistics            Run a test on statistics and their combinations.
	--useMCMC                   Use MCMC without likelihoods (Not implemented yet)
        --numSim=...                Number of simulations used to estimate the value(s)
 	-h, --help	            show this help

Examples:
	RhoEstimator.py --simulatedDataFile=sim_file obs_file
	
Description:
        Estimates rho using
"""

import sys, getopt, traceback, os, math
import snpsStatistics
import dataParsers,random,tempfile

def _run_():
    if len(sys.argv) == 1:
        print __doc__
        sys.exit(2)
	
    long_options_list = ["simulatedDataFile=", "useMCMC", "testStatistics", "help", "numSim="]
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", long_options_list)

    except:
        traceback.print_exc()
        print sys.exc_info()
        print __doc__
        sys.exit(2)
	

    obsDataFile = args[0]
    help = 0
    simulatedDataFile = None
    useMCMC=False
    testStatistics=False
    numSim = 100
	
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            help = 1
            print __doc__
        elif opt in ("--simulatedDataFile"):
            simulatedDataFile = arg
        elif opt in ("--numSim"):
            numSim = int(arg)
        elif opt in ("--testStatistics"):
            testStatistics = True
        elif opt in ("--useMCMC"):
            useMCMC = True
        else:
            if help==0:
                print "Unkown option!!\n"
                print __doc__
            sys.exit(2)
    
    

    if testStatistics:
        _testStatistics_()
    else:
        statistic  = snpsStatistics.StatComb1()
        obsData = dataParsers.parseMSFile(obsDataFile)[0]
        _estimateRho_(obsData,statistic,numSimulations=numSim) 


def _estimateRho_(obsData,statistic,acceptRate=0.1,simDatas=None,simRhos=None,numSimulations=100,priorType="Unif",priorParameters=(0,50),theta=5,numHaplotypes=50,type="Median"):
    """
    Estimate rho using a fixed acceptance rate
    """
    rankList = _getEstimateRanks_(obsData,statistic,simDatas=simDatas,simRhos=simRhos,numSimulations=numSimulations,priorType=priorType,priorParameters=priorParameters,theta=theta,numHaplotypes=numHaplotypes)
    goodRhos = []
    lim = int(round(acceptRate*len(rankList)))
    for i in range(0,lim):
        goodRhos.append(rankList[i][1])

    if type=="Mean":
        rhoEst = sum(goodRhos)/len(goodRhos)
    elif type=="Median":
        goodRhos.sort()
        rhoEst = goodRhos[len(goodRhos)/2]
        if len(goodRhos)%2==0:
            rhoEst = 0.5*(rhoEst+goodRhos[len(goodRhos)/2-1])
    
    print "Rho est =", rhoEst, ", num rho accepted =",len(goodRhos)
    return rhoEst

def _getEstimateRanks_(obsData,statistic,simDatas=None,simRhos=None,numSimulations=100,priorType="Unif",priorParameters=(0,50),theta=5,numHaplotypes=50):
    """
    Return a ranked list of rhos.
    """
    import dataParsers,random,tempfile

    # Load observed data
    statistic.setObsData(obsData)  #calculates the obs statistic and saves

    if not simDatas:
        print "Simulating data"
        (simDatas,simRhos) = _getRhoSimData_(numSimulations=numSimulations,priorType=priorType,priorParameters=priorParameters,theta=theta,numHaplotypes=numHaplotypes)

    # Scale statistics accordingly 
    statistic.scaleByVariance(simDatas)

    
    # Perform estimation
    euclidDistList = []
    for i in range(0,len(simDatas)):        
        # Calc statistic
        euclidDistList.append(statistic.euclidDist(simDatas[i]))
    
    rankList = zip(euclidDistList,simRhos)
    rankList.sort()
    return rankList



def _getRhoSimData_(numSimulations=100,priorType="Unif",priorParameters=(0,50),theta=5,numHaplotypes=50,numSamples=1):
    simRhos = []
    simDatas = []

    if priorType=="Unif":
        for i in range(0,numSimulations):  # For each simulation
            # Generate rho
            simRho = random.uniform(*priorParameters)
            simRhos.append(simRho)
    elif priorType=="FixedList":
        simRhos = priorParameters

    for i in range(0,numSimulations):  # For each simulation
        # Generate simulated data
        tmpFile = tempfile.mkstemp()
        os.close(tmpFile[0])        
        tmpFile = tmpFile[1] 
        _runMS_(tmpFile,rho=simRhos[i],theta=theta,numHaplotypes=numHaplotypes,numSamples=numSamples)

        # Load simulated data
        #print tmpFile
        simDatas.append(dataParsers.parseMSFile(tmpFile)[0])
    print "Simulated",len(simDatas),"ms datasets."
    return(simDatas,simRhos)
   

def _runMS_(outFile,rho=5,theta=5,numHaplotypes=50,numSamples=1):
    #./ms 50 1 -t 10 -r 50 10
    import os
    msPath = "/Users/bjarni/Projects/summer2007/ms/"

    runCommand = msPath+"ms "+str(numHaplotypes)+" "+str(numSamples)+" -t "+str(theta)+" "
    if rho:
        runCommand += " -r "+str(rho)+" "+str(100)
    runCommand += " > "+outFile
    
    #print runCommand
    os.system(runCommand)
    

def _testStatistics_():
    resDir = "/Users/bjarni/PM-599/"

    import pylab 
    parameterValues = []
    for i in range(0,100):
        parameterValues.append(0.25+(i/2.0))
    (obsDatas,obsRhos) = _getRhoSimData_(numSimulations=100,priorType="FixedList",priorParameters=parameterValues)    
    (simDatas,simRhos) = _getRhoSimData_(numSimulations=5000,priorType="Unif",priorParameters=(0,50))    

    statisticsList = [snpsStatistics.NumSegSites(),snpsStatistics.StupidStat(),snpsStatistics.NumHaplotypes(),snpsStatistics.MaxHaplotypeFreq(),snpsStatistics.NumSingletonHaplotypes(),snpsStatistics.PairwiseDiff(),snpsStatistics.StatComb1(),snpsStatistics.StatComb3(),snpsStatistics.MeanR2()]

    estRhosList=[] # [stat_index][rho_index]
    rmseList = []
    for statistic in statisticsList:
        estRhos = []
        rmseSum = 0
        print "Trying statistic :",statistic.name
        for (obsData,obsRho) in zip(obsDatas,obsRhos):
            print "Trying rho =",obsRho,":"
            estRho = _estimateRho_(obsData,statistic,simDatas=simDatas,simRhos=simRhos,type="Median")
            error = (estRho-obsRho)
            rmseSum += error*error
            estRhos.append(estRho)
        rmse = math.sqrt(rmseSum/float(len(obsRhos)))
        rmseList.append(rmse)
        estRhosList.append(estRhos)
        print statistic.name,": RMSE =",rmse

        pylab.plot(obsRhos,estRhos,".")
        pylab.plot([0, 50], [0, 50],"r--")
        pylab.axis([0, 50, 0, 50])
        pylab.title(statistic.name)
        pylab.xlabel('Observed rho')
        pylab.ylabel('Estimated rho')
        filename = resDir+'rho_plot1_'+statistic.name+'.pdf'
        print "Saving figure in file"
        pylab.savefig(resDir+'rho_plot1_'+statistic.name+'.pdf',format="pdf")
        pylab.clf()

        rankList = _getEstimateRanks_(obsDatas[20],statistic,simDatas=simDatas,simRhos=simRhos)
        l1 = []
        l2 = []
        for (dist,erho) in rankList:
            l1.append(dist)
            l2.append(erho)
        
        estRhos = []
        for i in range(1,len(l2)):
            estRhos.append(sum(l2[:i])/float(i))
        pylab.plot(l1, l2,".")
        buffer = abs(min(l1)-max(l1))*0.02
        pylab.plot([min(l1)-buffer,max(l1)+buffer], [10.25,10.25],"r--")
        pylab.plot(l1[1:],estRhos,"g-")
        pylab.axis([min(l1)-buffer,max(l1)+buffer, -1, 51])
        pylab.title(statistic.name)
        pylab.xlabel('Distance (euclidian)')
        pylab.ylabel('Simulated rho')
        pylab.savefig(resDir+'rho_plot2_'+statistic.name+'.pdf',format="pdf")
        pylab.clf()

        pylab.plot(range(0,len(l1)), l2,".")
        buffer = len(l1)*0.02
        pylab.plot([0-buffer,len(l1)-1+buffer], [10.25,10.25],"r--")
        pylab.plot(range(1,len(l1)),estRhos,"g-")
        pylab.axis([0-buffer,len(l1)-1+buffer, -1, 51])
        pylab.title(statistic.name)
        pylab.xlabel('Ranks of distance (euclidian)')
        pylab.ylabel('Simulated rho')
        pylab.savefig(resDir+'rho_plot3_'+statistic.name+'.pdf',format="pdf")
        pylab.clf()

if __name__ == '__main__':
    #_test1_()
    _run_()


