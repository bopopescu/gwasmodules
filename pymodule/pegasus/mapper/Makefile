SUBDIRS 	= extractor
CXXTARGETS	= CalculateColCorBetweenTwoHDF5 FindMaxLDBetweenPeakAndEachLocus
CXXLDFLAGS 	= $(HDF5Lib)
CTARGETS	= 

ExtraTargets	= CalculateMedianMeanOfInputColumn

CalculateMedianMeanOfInputColumn:	%:	%.o
	-$(CXXCOMPILER) $< $(CXXFLAGS) -lboost_program_options -lboost_iostreams -larmadillo -o $@


SRCS	= CalculateColCorBetweenTwoHDF5.cc FindMaxLDBetweenPeakAndEachLocus.cc CalculateMedianMeanOfInputColumn.cc

include ../../Makefile.common
#all:	$(CXXTARGETS) $(CTARGETS) CalculateColCorBetweenTwoHDF5
