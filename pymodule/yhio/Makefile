HDF5Lib  = -lhdf5_cpp -l hdf5 

.SUFFIXES:	.i .cc .c .so .o _wrap.o .py
.PHONY:	swig

BoostLib  = -lboost_program_options -lboost_iostreams

cxx_targets	= AbstractMatrixFileWalker
cxx_flags	= -O2
cxx_compiler	= g++
c_targets	= 
c_flags		= -O2
c_compiler	= gcc


all:	$(cxx_targets) $(c_targets)


$(cxx_targets):	%:	%.o
	-$(cxx_compiler) $< $(cxx_flags) -o $@ $(BoostLib)

.cc.o:
	$(cxx_compiler) $< $(cxx_flags) -c -o $@

$(c_targets):	%:	%.o
	$(c_compiler) $< $(c_flags) -o $@

.c.o:
	$(c_compiler) $< $(c_flags) -c -o $@


clean:
	-rm -f *.o $(cxx_targets) $(c_targets)
