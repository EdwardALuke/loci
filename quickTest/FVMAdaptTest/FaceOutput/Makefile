CASES = empty_ranks
SOURCES = test_face_output.cc
include $(dir $(lastword $(MAKEFILE_LIST)))../common.mk

ifdef CASE_WORK
empty_ranks: test_face_output.exe
	$(THREE_RANKS) ./test_face_output.exe --no-colors

test_face_output.exe: test_face_output.o $(LOCI_BASE)/lib/libfvmadaptfunc.$(LIB_SUFFIX)
	$(LD) -o $@ test_face_output.o -L$(LOCI_BASE)/lib -lfvmadaptfunc $(LIBS) $(LDFLAGS)
endif
