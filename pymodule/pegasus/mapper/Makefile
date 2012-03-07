HDF5Lib  = -lhdf5_cpp -l hdf5 

.SUFFIXES:	.i .cc .c .so .o _wrap.o .py
.PHONY:	swig

cxx_targets	= CalculateColCorBetweenTwoHDF5 FindMaxLDBetweenPeakAndEachLocus
cxx_flags	= -O2
cxx_compiler	= g++
c_targets	= 
c_flags		= -O2
c_compiler	= gcc


all:	$(cxx_targets) $(c_targets)


$(cxx_targets):	%:	%.o
	$(cxx_compiler) $(cxx_flags) $(HDF5Lib) $< -o $@

.cc.o:
	$(cxx_compiler) $(cxx_flags) -c $< -o $@

$(c_targets):	%:	%.o
	$(c_compiler) $(c_flags)  $< -o $@

.c.o:
	$(c_compiler) $(c_flags) -c $< -o $@


clean:
	-rm -f *.o $(cxx_targets) $(c_targets)
