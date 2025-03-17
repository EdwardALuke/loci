VARS_FILES= $(wildcard *.vars)
TEST_FILES= $(subst .vars,.test, $(VARS_FILES))



#include $(LOCI_BASE)/test.conf
#include ../../test.conf
#include test.conf

default: TestResults

FRC:

TestResults: $(TEST_FILES) FRC
	grep TEST $(TEST_FILES) > TestResults ; rm -f $(TEST_FILES) 


%.test : %.vars
	@echo Test $*
	@./test.sh $*.vars $*.test
	@grep TEST $*.test


clean:
	rm -f $(TEST_FILES)

distclean: clean
	rm -fr TestResults output debug *~


