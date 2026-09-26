CASES = plans
SOURCES = test_adaptation.cc test_refinement.cc
include $(dir $(lastword $(MAKEFILE_LIST)))../common.mk

ifdef CASE_WORK
plans: test_adaptation.exe
	$(TIMEOUT) --kill-after=10s $(TIME_LIMIT) ./test_adaptation.exe --no-colors

test_adaptation.exe: $(SOURCES:.cc=.o) $(LOCI_BASE)/lib/libfvmadaptfunc.$(LIB_SUFFIX)
	$(LD) -o $@ $(SOURCES:.cc=.o) -L$(LOCI_BASE)/lib -lfvmadaptfunc $(LIBS) $(LDFLAGS)
endif
