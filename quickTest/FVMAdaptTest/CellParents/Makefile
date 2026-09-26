CASES = serial parallel
TOOLS = ugrid2vog
SOURCES = test_cell_parents.cc
include $(dir $(lastword $(MAKEFILE_LIST)))../common.mk

ifdef CASE_WORK
tet.vog:
	cp $(INPUT)/tet.ugrid .
	$(SERIAL) $(LOCI_BASE)/bin/ugrid2vog $(LOCI_CMD_FLAGS) -o -m tet

serial: test_cell_parents.exe tet.vog
	$(SERIAL) ./test_cell_parents.exe --no-colors

parallel: test_cell_parents.exe tet.vog
	$(THREE_RANKS) ./test_cell_parents.exe --no-colors

test_cell_parents.exe: test_cell_parents.o $(LOCI_BASE)/lib/libfvmadaptfunc.$(LIB_SUFFIX)
	$(LD) -o $@ test_cell_parents.o -L$(LOCI_BASE)/lib -lfvmadaptfunc $(LIBS) $(LDFLAGS)
endif
