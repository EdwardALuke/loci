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
#
# Tools Makefile
# 

include $(LOCI_BASE)/Loci.conf

TARGET = libTools.$(LIB_SUFFIX)

INCLUDES = 

OBJS = fpe.o ftrn_reader.o identity.o eventNotify.o debugger.o \
       parse.o options_list.o expr.o  digraph.o unit_type.o \
	fpe_sun.o mpi_stubb.o intervalSet.o variable.o

compile: $(TARGET)

libTools.a: $(OBJS) 
	$(AR) $(ARFLAGS) libTools.a $(OBJS) 
	$(RANLIB) libTools.a

LIB_OBJS=$(OBJS:.o=_lo.o)
libTools.$(LIB_SUFFIX): $(LIB_OBJS)
	$(DYNAMIC_LD) $(DYNAMIC_LD_FLAGS) libTools.$(LIB_SUFFIX) $(LIB_FLAGS) $(LIB_OBJS)  

JUNK =  *~ core ii_files ti_files KCC_files so_locations rii_files

clean:
	rm -fr $(LIB_OBJS) $(OBJS) $(TARGET) $(JUNK) libTools.$(LIB_SUFFIX)

DEPEND_FILES=$(subst .o,.d,$(OBJS))

uninstall:
	rm -f $(LOCI_INSTALL_PATH)/lib/$(TARGET)
	rm -fr $(LOCI_INSTALL_PATH)/include/Tools

distclean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK) $(DEPEND_FILES) libTools.$(LIB_SUFFIX)

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean distclean ),)
-include $(DEPEND_FILES)
endif
