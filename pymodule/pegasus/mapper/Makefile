SUBDIRS 	= extractor
CXXStaticTargets	= CalculateColCorBetweenTwoHDF5 FindMaxLDBetweenPeakAndEachLocus
CXXLDFLAGS 	= $(HDF5Lib)
CTARGETS	= 

ExtraTargets	= CalculateMedianMeanOfInputColumn

CalculateMedianMeanOfInputColumn:	%:	%.o
	-$(CXXCOMPILER) $< $(CXXFLAGS) -o $@ $(StaticStartFlags) -lboost_program_options -lboost_iostreams $(StaticEndFlags) -larmadillo
#2015.03.13 didn't find static armadillo in ubuntu libraries

SRCS	= CalculateColCorBetweenTwoHDF5.cc FindMaxLDBetweenPeakAndEachLocus.cc CalculateMedianMeanOfInputColumn.cc

include ../../Makefile.common
#all:	$(CXXTARGETS) $(CTARGETS) CalculateColCorBetweenTwoHDF5
