VARS_FILES= $(wildcard *.vars)
TEST_FILES= $(subst .vars,.test, $(VARS_FILES))

include $(LOCI_BASE)/Loci.conf
include $(LOCI_BASE)/version.conf
include $(TEST_BASE)/test.conf

default: TestResults

FRC:

TestResults: $(TEST_FILES) FRC
	cat $(TEST_FILES) > TestResults ; rm $(TEST_FILES)


%.test : %.vars
	@echo Test $*
	@ MPIRUN="$(MPIRUN)" LOCI_FLAGS="$(LOCI_CMD_FLAGS)" ./test.sh $*.vars $*.log
	tail -n 1 $*.log > $*.test
	@grep TEST $*.test


clean:
	rm -f $(TEST_FILES)

distclean: clean
	rm -fr TestResults output debug *~ *.log 


