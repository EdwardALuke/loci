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

INCLUDES = -I./include -I$(SOURCE_DIR)/include -I$(SOURCE_DIR)/../Containers/include

TESTS = test_options_units test_options_list

.PHONY: test_options_units test_options_list


tests: $(TESTS)

test_options_units.o: $(SOURCE_DIR)/test_options_units.cc
	$(CXX) $(COPT) $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) -o $@ -c $<

test_options_list.o: $(SOURCE_DIR)/test_options_list.cc
	$(CXX) $(COPT) $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) -o $@ -c $<

test_options_units: test_options_units.o
	$(LD) -o $@.exe $^ $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)
	./$@.exe
	rm $@.exe

test_options_list: test_options_list.o
	$(LD) -o $@.exe $^ $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)
	./$@.exe
	rm $@.exe


clean:
	rm -fr *.o *.exe *.d

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~) $(wildcard crash_dump.*) core debug output $(wildcard *.o)

distclean: clean
	rm -fr test_options_units test_options_list $(JUNK) $(DEPEND_FILES)

DEPEND_FILES=$(subst .o,.d,$(OBJS))

# include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean distclean ),)
-include $(DEPEND_FILES)
endif
