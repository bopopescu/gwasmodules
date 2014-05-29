SUBDIRS 	= 
CXXTARGETS	= 
CXXLDFLAGS 	= $(BoostLib)
CTARGETS	= 

ExtraTargets	= ReduceSameChromosomeAlignmentDepthFiles


SharedLibTargets	= ~/script/pymodule/yhio/AbstractMatrixFileWalkerCC.so

SRCS	= ReduceSameChromosomeAlignmentDepthFiles.cc

ReduceSameChromosomeAlignmentDepthFiles:	%:	%.o $(SharedLibTargets)
	-$(CXXCOMPILER) $< $(SharedLibTargets) $(CXXFLAGS) $(CXXLDFLAGS) -o $@

include ../../Makefile.common
