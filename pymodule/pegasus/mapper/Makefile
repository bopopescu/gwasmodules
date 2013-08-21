SUBDIRS 	= extractor
CXXTARGETS	= CalculateColCorBetweenTwoHDF5 FindMaxLDBetweenPeakAndEachLocus
CXXLDFLAGS 	= $(HDF5Lib)
CTARGETS	= 

ExtraTargets	= CalculateMedianMeanOfInputColumn

CalculateMedianMeanOfInputColumn:	%:	%.o
	-$(CXXCOMPILER) $< $(CXXFLAGS) -lboost_program_options -lboost_iostreams -larmadillo -o $@



include ../../Makefile.common
#all:	$(CXXTARGETS) $(CTARGETS) CalculateColCorBetweenTwoHDF5
