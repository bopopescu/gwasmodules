
BoostPythonTarget	= GADA.so
BoostPythonLib  = -lboost_python
# 2010-6-8 "-L ~/lib64/boost/lib" removed from BoostPythonLib. "-lboost_python-mt" changed to "-lboost_python"

SharedLibFlags	= -shared -fPIC
IncludeDirs	= `python-config --includes`
# 2010-6-8 "-I ~/lib64/boost/include/" removed from IncludeDirs

BoostLib	= -lboost_program_options -lboost_iostreams

CFLAGS	= -O3 
CXXFLAGS	= -fPIC $(CFLAGS)
CXXCOMPILER	= g++
CXXLDFLAGS	= $(BoostLib)

all:	GADA $(BoostPythonTarget)

# for boost.python target(shared libraries)
$(BoostPythonTarget):	%.so:	%.o BaseGADA.o
	-$(CXXCOMPILER) $< $(SharedLibFlags) BaseGADA.o -o $@ $(BoostPythonLib) -lm
	#2013.08.16 -l... should be put in the last. If they are put in front of input files, they will be ignored!

.cc.o:
	-$(CXXCOMPILER) $< $(CXXFLAGS) $(IncludeDirs) $(SharedLibFlags) -c -o $@

.c.o:
	-$(CXXCOMPILER) $< $(CXXFLAGS) $(IncludeDirs) $(SharedLibFlags) -c -o $@

GADA:	GADA.cc BaseGADA.cc BaseGADA.h consoledefines.h
	-$(CXXCOMPILER) GADA.cc BaseGADA.cc $(CXXFLAGS) $(IncludeDirs) -DGADABIN -o GADA -lm $(CXXLDFLAGS)
	# `python-config --libs` and $(BoostPythonLib) would be needed if boost.python is used. However, -DGADABIN renders those code off. 

clean:
	-rm GADA *.o $(BoostPythonTarget)
