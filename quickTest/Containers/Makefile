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

###########################################################################
# No changes required below this line
###########################################################################

include $(LOCI_BASE)/Loci.conf

INCLUDES = -I./include

TESTS = test_vtypes test_VFAD

.PHONY: test_vtypes test_VFAD


tests: $(TESTS)

test_vtypes: test_vtypes.o  
	$(LD) -o $@.exe $^ $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)
	./$@.exe
	rm $@.exe

test_VFAD: test_VFAD.o 
	$(LD) -o $@.exe $^ $(LOCA# L_LIBS) $(LIBS) $(LDFLAGS)
	./$@.exe
	rm $@.exe




#tests: $(TARGET) FORCE
#	$(MAKE) -C tests TestResults LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)

clean:
	rm -fr *.o *.exe *.d

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~) $(wildcard crash_dump.*)  core debug output $(wildcard *.o) 
# ".cc" files created from .loci files
LOCI_LPP_FILES = $(LOCI_FILES:.loci=.$(LPP_I_SUFFIX))

distclean: clean
	rm -fr test_vtypes $(JUNK) $(LOCI_LPP_FILES) $(DEPEND_FILES)

DEPEND_FILES=$(subst .o,.d,$(OBJS))

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean distclean ),)
-include $(DEPEND_FILES)
endif

