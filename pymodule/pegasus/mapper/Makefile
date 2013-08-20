
cxx_targets	= CalculateColCorBetweenTwoHDF5 FindMaxLDBetweenPeakAndEachLocus
cxxLinkFlags 	= $(HDF5Lib)
c_targets	= 

extraTargets	= CalculateMedianMeanOfInputColumn

CalculateMedianMeanOfInputColumn:	%:	%.o
	-$(cxx_compiler) $< $(cxx_flags) -lboost_program_options -lboost_iostreams -larmadillo -o $@



include ../../Makefile.common
#all:	$(cxx_targets) $(c_targets) CalculateColCorBetweenTwoHDF5
