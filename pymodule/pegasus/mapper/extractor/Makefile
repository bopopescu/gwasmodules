
CXXTARGETS	= 
CXXLDFLAGS	= $(BoostLib)
CTARGETS	= 

ExtraTargets	= SelectRowsFromMatrixCC
SharedLibTargets	= ~/script/pymodule/yhio/AbstractMatrixFileWalkerCC.so

SRCS	= SelectRowsFromMatrixCC.cc

SelectRowsFromMatrixCC:	%:	%.o $(SharedLibTargets)
	-$(CXXCOMPILER) $< $(SharedLibTargets) $(CXXFLAGS) $(CXXLDFLAGS) -o $@


include ../../../Makefile.common
