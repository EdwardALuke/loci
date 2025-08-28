##############################################################################
#
# Copyright 2008, Mississippi State University
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
##############################################################################
LOCI_BASE?=..
include $(LOCI_BASE)/Loci.conf

MODULE_NAME  = fvmadapt

LIB_NAME = fvmadaptfunc

MODULE_FILES = $(call loci_compile_files,)
MODULE_OBJS  = $(call loci_file2objs,$(MODULE_FILES))

LIB_FILES = $(call loci_compile_files,library/)
LIB_OBJS  = $(call loci_file2objs,$(LIB_FILES))

LOCAL_LIBS      = $(LIBXML2_LIBS)

JUNK =		*~  core 

LIB_OBJS1=$(MODULE_OBJS:.o=_lo.o)
LIB_OBJS2=$(LIB_OBJS:.o=_lo.o)

INCLUDES=	-I$(LOCI_BASE)/include/FVMAdapt $(LIBXML2_INCLUDE) $(LIBXML2_DEFINE)

all: $(MODULE_NAME)_m.so lib$(LIB_NAME).$(LIB_SUFFIX) 

$(MODULE_NAME)_m.so: $(LIB_OBJS1)
	$(SHARED_LD) $(SHARED_LD_FLAGS) $(MODULE_NAME)_m.so $(LIB_FLAGS) $(LIB_OBJS1) $(LOCAL_LIBS) $(BASE_LIBS)

lib$(LIB_NAME).$(LIB_SUFFIX): $(LIB_OBJS2)
	$(DYNAMIC_LD) $(DYNAMIC_LD_FLAGS) lib$(LIB_NAME).$(LIB_SUFFIX) $(LIB_FLAGS) $(LIB_OBJS2) $(BASE_LIBS)


FRC :

clean:
	rm -fr $(MODULE_OBJS)  $(LIB_OBJS) $(LIB_OBJS1) $(LIB_OBJS2) $(MODULE_NAME)_m.so lib$(LIB_NAME).$(LIB_SUFFIX) $(JUNK)

# intermediate files created from .loci files
LOCI_INTERMEDIATE_FILES = $(subst .loci,.$(LPP_I_SUFFIX), $(wildcard *.loci) )

uninstall:
	rm -f $(LOCI_INSTALL_PATH)/lib/$(MODULE_NAME)_m.so $(LOCI_INSTALL_PATH)/lib/lib$(LIB_NAME).$(LIB_SUFFIX)
	rm -fr $(LOCI_INSTALL_PATH)/include/FVMAdapt

distclean: clean
	rm -fr $(DEPEND_FILES) $(LOCI_INTERMEDIATE_FILES)

DEPEND_FILES=$(subst .o,.d,$(MODULE_OBJS),$(LIB_OBJS))
#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean distclean uninstall),)
-include $(DEPEND_FILES)
endif
