SUBDIRS 	= 
CXXTARGETS	=
CXXLDFLAGS 	= $(BoostLib)
CTARGETS	= 

ExtraTargets	= CheckTwoVCFOverlapCC
SharedLibTargets	= ~/script/pymodule/yhio/AbstractMatrixFileWalkerCC.so 

CheckTwoVCFOverlapCC:	%:	%.o $(SharedLibTargets)
	-$(CXXCOMPILER) $< ~/script/vcflib/libvcf.a $(SharedLibTargets) $(CXXFLAGS) -o $@ $(CXXLDFLAGS)


SRCS	= CheckTwoVCFOverlapCC.cc

include ../../Makefile.common
#all:	$(CXXTARGETS) $(CTARGETS) CalculateColCorBetweenTwoHDF5
