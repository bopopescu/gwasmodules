#!/bin/sh

condor_q -long -attributes Owner,Cmd,JobStatus,RemoteHost,ClusterId,DAGNodeName,DAGParentNodeNames,pegasus_wf_name,JobRunCount,ExitStatus,LastJobStatus,RemoteUserCpu,RemoteSysCpu 2>&1 |less

#condor_q -cputime -dag|less
#condor_q -dag|less
