###############################################################################
#
# Copyright 2008-2026, Mississippi State University
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
include $(LOCI_BASE)/version.conf

INCLUDES = -I../Core -I$(LOCI_BASE)/include
LIB_PATH = $(LOCI_BASE)/lib:$(LD_LIBRARY_PATH)
DY_PATH = $(LOCI_BASE)/lib:$(DYLD_LIBRARY_PATH)
SOURCES = fvmadapt_illustrations.cc ../Core/fixtures.cc
OBJS = fvmadapt_illustrations.o fixtures.o

.PHONY: default examples clean distclean FRC

default: examples

examples: fvmadapt_illustrations.exe
	@rm -fr output
	@mkdir -p output
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) ./fvmadapt_illustrations.exe --write-dir output

fvmadapt_illustrations.exe: $(OBJS) $(LOCI_BASE)/lib/libfvmadaptfunc.$(LIB_SUFFIX)
	$(LD) -o $@ $(OBJS) -L$(LOCI_BASE)/lib -lfvmadaptfunc $(LIBS) $(LDFLAGS)

fvmadapt_illustrations.o: fvmadapt_illustrations.cc
	$(CXX) $(COPT) $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) -o $@ -c $<

fixtures.o: ../Core/fixtures.cc ../Core/fvmadapt_core.h
	$(CXX) $(COPT) $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) -o $@ -c ../Core/fixtures.cc

clean:
	rm -fr *.o *.exe *.d *.log TestResults output

distclean: clean

FRC:
