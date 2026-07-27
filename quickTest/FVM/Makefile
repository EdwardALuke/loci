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

# List your object files here
OBJS  =  bcFixedMass.o bcImpermeable.o bcOutflow.o gridReader.o gridSetup.o fluxInviscid.o main.o IOgrid.o IOplot.o solverPeriodic.o solverRungeKutta.o solverSetup.o solverTimestep.o solverIC.o solverEoS.o
INCLUDES = -I./include

# List the name of your compiled program here
TARGET = fluidTest


#############################################################################
# No changes needed below this line

include $(LOCI_BASE)/Loci.conf
#COPT= -g -O0

.PHONY: FORCE

tests: fluidTest ndiff FORCE
	$(MAKE) -C tests TestResults LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)
	mv tests/TestResults .


fluidTest: $(OBJS)
	$(LD) -o $(TARGET) $(OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

ndiff: ndiff.o
	$(LD) -o  ndiff ndiff.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

clean:
	rm -fr $(OBJS) $(TARGET) TestResults
	$(MAKE) -C tests clean

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~) $(wildcard crash_dump.*)  core debug output $(wildcard *.o) ndiff
# ".cc" files created from .loci files
LOCI_INTERMEDIATE_FILES = $(subst .loci,.cc, $(wildcard *.loci) )

distclean: clean
	rm -fr $(JUNK) $(LOCI_INTERMEDIATE_FILES) $(DEPEND_FILES)
	$(MAKE) -C tests distclean

DEPEND_FILES=$(subst .o,.d,$(OBJS))

#include automatically generated dependencies
-include $(DEPEND_FILES)
