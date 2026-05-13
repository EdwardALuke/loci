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
# LOCI_BASE should be set before the Makefile
# Set TARGET to the name of your program
# Set FILES to list '.loci' files that will be compiled into your module, or

TARGET=FVMModUnitTest
FILES = $(AUTOMATIC_FILES)
INCLUDES = -I./include

###########################################################################
# No changes required below this line
###########################################################################

include $(LOCI_BASE)/Loci.conf

AUTOMATIC_FILES = $(call loci_compile_files,)
AUTOMATIC_OBJS = $(call loci_file2objs,$(FILES))

OBJS := $(AUTOMATIC_OBJS)


$(TARGET): $(OBJS)
	$(LD) -o $(TARGET) $(OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

.PHONY: FORCE

tests: $(TARGET) FORCE
	$(MAKE) -C tests TestResults LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)
	mv tests/TestResults .

clean:
	rm -fr $(OBJS) $(TARGET) TestResults
	$(MAKE)  -C tests clean

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~) $(wildcard crash_dump.*)  core debug output $(wildcard *.o) 
# ".cc" files created from .loci files
LOCI_LPP_FILES = $(LOCI_FILES:.loci=.$(LPP_I_SUFFIX))

distclean: clean
	rm -fr  $(JUNK) $(LOCI_LPP_FILES) $(DEPEND_FILES)
	$(MAKE) -C tests distclean

DEPEND_FILES=$(subst .o,.d,$(OBJS))

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean distclean ),)
-include $(DEPEND_FILES)
endif
