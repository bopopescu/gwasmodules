
cxx_targets	= 
cxxLinkFlags	= $(BoostLib)
c_targets	= 

extraTargets	= SelectRowsFromMatrix

SelectRowsFromMatrix:	%:	%.o ../../../yhio/AbstractMatrixFileWalker.so
	-$(cxx_compiler) $< ../../../yhio/AbstractMatrixFileWalker.so $(cxx_flags) $(cxxLinkFlags) -o $@


include ../../../Makefile.common
