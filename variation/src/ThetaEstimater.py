#!/usr/bin/env python2.5
"""
Usage: ThetaEstimator.py [OPTIONS] OBS_DATA_FILE

Option:

	--simulatedDataFile=...     A file with simulated ms data, to be used for the rejection scheme
	--testStatistics            Run a test on statistics and their combinations.
	--useMCMC                   Use MCMC without likelihoods (Not implemented yet)
        --numSim=...                Number of simulations used to estimate the value(s)
 	-h, --help	            show this help

Examples:
	ThetaEstimator.py --simulatedDataFile=sim_file obs_file
	
Description:
        Estimates theta using
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
        _estimateTheta_(obsData,statistic,numSimulations=numSim) 



def _estimateTheta_(obsData,statistic,acceptRate=0.1,simDatas=None,simThetas=None,numSimulations=100,priorType="Unif",priorParameters=(0,50),rho=None,numHaplotypes=50,type="Median"):
    """
    Estimate theta using a fixed acceptance rate
    """
    rankList = _getEstimateRanks_(obsData,statistic,simDatas=simDatas,simThetas=simThetas,numSimulations=numSimulations,priorType=priorType,priorParameters=priorParameters,rho=rho,numHaplotypes=numHaplotypes)
    goodThetas = []
    lim = int(round(acceptRate*len(rankList)))
    for i in range(0,lim):
        goodThetas.append(rankList[i][1])

    if type=="Mean":
        thetaEst = sum(goodThetas)/len(goodThetas)
    if type=="Median":
        goodThetas.sort()
        thetaEst = goodThetas[len(goodThetas)/2]
        if len(goodThetas)%2==0:
            thetaEst = 0.5*(thetaEst+goodThetas[len(goodThetas)/2-1])
    
    print "Theta est =", thetaEst, ", num theta accepted =",len(goodThetas)
    return thetaEst

def _getEstimateRanks_(obsData,statistic,simDatas=None,simThetas=None,numSimulations=100,priorType="Unif",priorParameters=(0,50),rho=None,numHaplotypes=50):
    """
    Return a ranked list of thetas.
    """
    import dataParsers,random,tempfile

    # Load observed data
    statistic.setObsData(obsData)  #calculates the obs statistic and saves

    if not simDatas:
        print "Simulating data"
        (simDatas,simThetas) = _getThetaSimData_(numSimulations=numSimulations,priorType=priorType,priorParameters=priorParameters,rho=rho,numHaplotypes=numHaplotypes)

    # Scale statistics accordingly 
    statistic.scaleByVariance(simDatas)

    
    # Perform estimation
    euclidDistList = []
    for i in range(0,len(simDatas)):        
        # Calc statistic
        euclidDistList.append(statistic.euclidDist(simDatas[i]))
    
    rankList = zip(euclidDistList,simThetas)
    rankList.sort()
    return rankList



def _getThetaSimData_(numSimulations=100,priorType="Unif",priorParameters=(0,50),rho=None,numHaplotypes=50,numSamples=1):
    simThetas = []
    simDatas = []

    if priorType=="Unif":
        for i in range(0,numSimulations):  # For each simulation
            # Generate theta
            simTheta = random.uniform(*priorParameters)
            simThetas.append(simTheta)
    elif priorType=="FixedList":
        simThetas = priorParameters

    for i in range(0,numSimulations):  # For each simulation
        # Generate simulated data
        tmpFile = tempfile.mkstemp()
        os.close(tmpFile[0])        
        tmpFile = tmpFile[1] 
        _runMS_(tmpFile,theta=simThetas[i],rho=rho,numHaplotypes=numHaplotypes,numSamples=numSamples)

        # Load simulated data
        #print tmpFile
        simDatas.append(dataParsers.parseMSFile(tmpFile)[0])
    print "Simulated",len(simDatas),"ms datasets."
    return(simDatas,simThetas)
   

def _runMS_(outFile,theta=5,rho=None,numHaplotypes=50,numSamples=1):
    #./ms 50 1 -t 10 -r 50 10
    import os
    msPath = "/Users/bjarni/Projects/summer2007/ms/"

    runCommand = msPath+"ms "+str(numHaplotypes)+" "+str(numSamples)+" -t "+str(theta)+" "
    if rho:
        runCommand += " -r "+str(rho[0])+" "+str(rho[1])
    runCommand += " > "+outFile
    
    #print runCommand
    os.system(runCommand)
    

def _testStatistics_():
    resDir = "/Users/bjarni/PM-599/"

    import pylab 
    parameterValues = []
    for i in range(0,100):
        parameterValues.append(0.25+(i/2.0))
    (obsDatas,obsThetas) = _getThetaSimData_(numSimulations=100,priorType="FixedList",priorParameters=parameterValues)    
    (simDatas,simThetas) = _getThetaSimData_(numSimulations=5000,priorType="Unif",priorParameters=(0,50))    

    statisticsList = [snpsStatistics.NumSegSites(),snpsStatistics.StupidStat(),snpsStatistics.NumHaplotypes(),snpsStatistics.MaxHaplotypeFreq(),snpsStatistics.NumSingletonHaplotypes(),snpsStatistics.PairwiseDiff(),snpsStatistics.StatComb1(),snpsStatistics.StatComb2(),snpsStatistics.MeanR2()]

    estThetasList=[] # [stat_index][theta_index]
    rmseList = []
    for statistic in statisticsList:
        estThetas = []
        rmseSum = 0
        print "Trying statistic :",statistic.name
        for (obsData,obsTheta) in zip(obsDatas,obsThetas):
            print "Trying theta =",obsTheta,":"
            estTheta = _estimateTheta_(obsData,statistic,simDatas=simDatas,simThetas=simThetas,type="Median")
            error = (estTheta-obsTheta)
            rmseSum += error*error
            estThetas.append(estTheta)
        rmse = math.sqrt(rmseSum/float(len(obsThetas)))
        rmseList.append(rmse)
        estThetasList.append(estThetas)
        print statistic.name,": RMSE =",rmse

        pylab.plot(obsThetas,estThetas,".")
        pylab.plot([0, 50], [0, 50],"r--")
        pylab.axis([0, 50, 0, 50])
        pylab.title(statistic.name)
        pylab.xlabel('Observed theta')
        pylab.ylabel('Estimated theta')
        pylab.savefig(resDir+'theta_plot1_'+statistic.name+'.pdf',format="pdf")
        pylab.clf()

        rankList = _getEstimateRanks_(obsDatas[20],statistic,simDatas=simDatas,simThetas=simThetas)
        l1 = []
        l2 = []
        for (dist,etheta) in rankList:
            l1.append(dist)
            l2.append(etheta)
        
        estThetas = []
        for i in range(1,len(l2)):
            estThetas.append(sum(l2[:i])/float(i))
        pylab.plot(l1, l2,".")
        buffer = abs(min(l1)-max(l1))*0.02
        pylab.plot([min(l1)-buffer,max(l1)+buffer], [10.25,10.25],"r--")
        pylab.plot(l1[1:],estThetas,"g-")
        pylab.axis([min(l1)-buffer,max(l1)+buffer, -1, 51])
        pylab.title(statistic.name)
        pylab.xlabel('Distance (euclidian)')
        pylab.ylabel('Simulated theta')
        pylab.savefig(resDir+'theta_plot2_'+statistic.name+'.pdf',format="pdf")
        pylab.clf()

        pylab.plot(range(0,len(l1)), l2,".")
        buffer = len(l1)*0.02
        pylab.plot([0-buffer,len(l1)-1+buffer], [10.25,10.25],"r--")
        pylab.plot(range(1,len(l1)),estThetas,"g-")
        pylab.axis([0-buffer,len(l1)-1+buffer, -1, 51])
        pylab.title(statistic.name)
        pylab.xlabel('Ranks of distance (euclidian)')
        pylab.ylabel('Simulated theta')
        pylab.savefig(resDir+'theta_plot3_'+statistic.name+'.pdf',format="pdf")
        pylab.clf()

if __name__ == '__main__':
    #_test1_()
    _run_()


