#!/bin/sh

condor_q -long -attributes Owner,Cmd,JobStatus,RemoteHost,ClusterId,DAGNodeName,DAGParentNodeNames,pegasus_wf_name,JobRunCount,ExitStatus,LastJobStatus |less
