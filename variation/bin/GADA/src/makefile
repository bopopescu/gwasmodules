
BoostPythonTarget	= GADA.so
BoostPythonLib  = -lboost_python
# 2010-6-8 "-L ~/lib64/boost/lib" removed from BoostPythonLib. "-lboost_python-mt" changed to "-lboost_python"
CFLAGS	= -O3
CXXFLAGS	= -fPIC $(CFLAGS)
SharedLibFlags	= -shared -fPIC
IncludeDirs = `python-config --includes`
# 2010-6-8 "-I ~/lib64/boost/include/" removed from IncludeDirs
CXXCOMPILER = g++

all:	GADA $(BoostPythonTarget)

# for boost.python target(shared libraries)
$(BoostPythonTarget):	%.so:	%.o BaseGADA.o
	-$(CXXCOMPILER) $(BoostPythonLib) $(SharedLibFlags) -lm $< BaseGADA.o -o $@


.cc.o:
	-$(CXXCOMPILER) $(CXXFLAGS) $(IncludeDirs) $(SharedLibFlags) -c $< -o $@

.c.o:
	-$(CXXCOMPILER) $(CXXFLAGS) $(IncludeDirs) $(SharedLibFlags) -c $< -o $@

GADA:	GADA.cc BaseGADA.c BaseGADA.h consoledefines.h
	-$(CXXCOMPILER) $(CXXFLAGS) $(IncludeDirs) -DGADABIN -o GADA GADA.cc BaseGADA.c -lm
	# `python-config --libs` and $(BoostPythonLib) would be needed if boost.python is used. However, -DGADABIN renders those code off. 

clean:
	-rm GADA *.o $(BoostPythonTarget)
