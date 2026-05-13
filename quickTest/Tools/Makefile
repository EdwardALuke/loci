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

INCLUDES = -I../contrib/doctest -I$(LOCI_BASE)/include

AUTOMATIC_FILES = $(call loci_compile_files,)
AUTOMATIC_OBJS = $(call loci_file2objs,$(AUTOMATIC_FILES))
AUTOMATIC_TESTS= $(subst .o,.test, $(AUTOMATIC_OBJS))
TESTS = $(AUTOMATIC_TESTS)
#vtypes.test VFAD.test

TestResults: $(TESTS)
	cat $(TESTS) > TestResults; rm $(TESTS)

%.test: %.o
	$(LD) -o $*.exe $^ $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)
	@(./$*.exe >& $*.log || /bin/true)
	tail -n 1 $*.log > $*.log1
	echo -n Containers/$*":" > $@
	@if grep -q "SUCCESS!" $*.log1 ; then echo " PASSED!"; else echo " FAILED!"; fi >>$@
	rm $*.log1 $*.exe

clean:
	rm -fr *.o *.exe *.d *.log *.log1

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~) $(wildcard crash_dump.*)  core debug output $(wildcard *.o) 
# ".cc" files created from .loci files
LOCI_LPP_FILES = $(LOCI_FILES:.loci=.$(LPP_I_SUFFIX))

distclean: clean
	rm -fr TestResults $(JUNK) $(LOCI_LPP_FILES) $(DEPEND_FILES) $(TESTS)

DEPEND_FILES=$(subst .o,.d,$(OBJS))

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean distclean ),)
-include $(DEPEND_FILES)
endif
