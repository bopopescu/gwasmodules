
BoostPythonTarget	= GADA.so
BoostPythonLib  = -lboost_python
# 2010-6-8 "-L ~/lib64/boost/lib" removed from BoostPythonLib. "-lboost_python-mt" changed to "-lboost_python"

SharedLibFlags	= -shared -fPIC
IncludeDirs	= `python-config --includes` -I ~/script/
# 2010-6-8 "-I ~/lib64/boost/include/" removed from IncludeDirs

BoostLib	= -lboost_program_options -lboost_iostreams

CFLAGS	= -O3 
CXXFLAGS	= -fPIC $(CFLAGS) -std=c++0x
CXXCOMPILER	= g++
CXXLDFLAGS	= $(BoostLib)

all:	depend GADA $(BoostPythonTarget)

# for boost.python target(shared libraries)
$(BoostPythonTarget):	%.so:	GADA_shared.o BaseGADA_shared.o
	#2013.08.16 -l... should be put in the last. If they are put in front of input files, they will be ignored!
	-$(CXXCOMPILER) $< BaseGADA_shared.o -DGADALib -o $@ $(SharedLibFlags) -lm $(BoostPythonLib)

GADA_shared.o:	GADA.cc
	-$(CXXCOMPILER) $< -c -o $@ $(CXXFLAGS) -DGADALib $(IncludeDirs) $(SharedLibFlags)
BaseGADA_shared.o:	BaseGADA.cc
	-$(CXXCOMPILER) $< -c -o $@ $(CXXFLAGS) -DGADALib $(IncludeDirs) $(SharedLibFlags)

.cc.o:
	-$(CXXCOMPILER) $< -c -o $@ $(CXXFLAGS) $(IncludeDirs) -DGADABIN $(SharedLibFlags)

.c.o:
	-$(CXXCOMPILER) $< -c -o $@ $(CXXFLAGS) $(IncludeDirs) -DGADABIN $(SharedLibFlags)


GADA:	GADA.o BaseGADA.o BaseGADA.h consoledefines.h
	# `python-config --libs` and $(BoostPythonLib) would be needed if boost.python is used. However, -DGADABIN renders those code off. 
	-$(CXXCOMPILER) GADA.o BaseGADA.o -o GADA $(CXXFLAGS) $(IncludeDirs) -DGADABIN -lm $(CXXLDFLAGS)

depend: .depend

.depend: GADA.cc BaseGADA.cc
	rm -f ./.depend
	-$(CXXCOMPILER) $^ $(CXXFLAGS) $(IncludeDirs) -MM >>./.depend;

include .depend

clean:
	-rm GADA *.o $(BoostPythonTarget)
