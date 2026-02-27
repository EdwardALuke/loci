###############################################################################
#
# Copyright 2008-2025, Mississippi State University
#
# This file is part of the Loci Framework.
#
# The Loci Framework is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The Loci Framework is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
#
###############################################################################

include $(LOCI_BASE)/Loci.conf

SOURCE_DIR := $(dir $(realpath $(firstword $(MAKEFILE_LIST))))
vpath %.cc $(SOURCE_DIR)

INCLUDES = -I$(SOURCE_DIR)/../contrib/doctest

TESTS = test_rule
OBJS = $(addsuffix .o,$(TESTS))
TEST_RESULTS = $(addsuffix .test,$(TESTS))

.PHONY: tests $(TESTS)

tests: TestResults

TestResults: $(TEST_RESULTS)
	cat $(TEST_RESULTS) > TestResults; rm $(TEST_RESULTS)

$(TEST_RESULTS): %.test: %.o
	$(LD) -o $*.exe $^ $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)
	@(LD_LIBRARY_PATH=$(LD_LIBRARY_PATH) DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH) ./$*.exe >& $*.log || true)
	tail -n 1 $*.log > $*.log1
	echo -n SystemTests/$*":" > $@
	@if grep -q "SUCCESS!" $*.log1 ; then echo " PASSED!"; else echo " FAILED!"; fi >>$@
	rm $*.log1 $*.exe

$(TESTS): %: %.test

clean:
	rm -fr *.o *.d *.exe *.log *.log1 *.test

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~) $(wildcard crash_dump.*) core debug output

distclean: clean
	rm -fr TestResults $(JUNK) $(DEPEND_FILES) $(TEST_RESULTS)

DEPEND_FILES=$(subst .o,.d,$(OBJS))

# include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean distclean ),)
-include $(DEPEND_FILES)
endif
