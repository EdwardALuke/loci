###############################################################################
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
###############################################################################
# This is the Loci base directory
LOCI_BASE ?= /usr/tmp/lush/Loci
# Put objects in the module here
OBJS =	geometry.o holeCutting.o interpolation.o motion.o
MODULE_NAME = fvmoverset

#############################################################################
# No changes needed below this line

include $(LOCI_BASE)/Loci.conf

JUNK = *~  core ti_files ii_files rii_files

LIB_OBJS=$(OBJS:.o=_lo.o)
FAD_OBJS=$(OBJS:.o=_FAD.o)

default: $(MODULE_NAME)_m.so $(MODULE_NAME)FAD_m.so

all: $(MODULE_NAME)_m.so $(MODULE_NAME)FAD_m.so

$(MODULE_NAME)_m.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) $(MODULE_NAME)_m.so $(LIB_FLAGS) $(LIB_OBJS)

$(MODULE_NAME)FAD_m.so: $(FAD_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) $(MODULE_NAME)FAD_m.so $(LIB_FLAGS) $(FAD_OBJS)

FRC : 

clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(FAD_OBJS) $(MODULE_NAME)_m.so $(MODULE_NAME)FAD_m.so 

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~)  core 
# intermediate files created from .loci files
LOCI_INTERMEDIATE_FILES = $(subst .loci,.$(LPP_I_SUFFIX), $(wildcard *.loci) )

uninstall:
	rm -f $(LOCI_INSTALL_PATH)/lib/$(MODULE_NAME)_m.so $(LOCI_INSTALL_PATH)/lib/$(MODULE_NAME)FAD_m.so
	rm -fr $(LOCI_INSTALL_PATH)/include/FVMOverset


distclean: 
	rm -fr $(OBJS) $(LIB_OBJS) $(FAD_OBJS) $(MODULE_NAME)_m.so $(MODULE_NAME)FAD_m.so $(JUNK) $(LOCI_INTERMEDIATE_FILES) $(DEPEND_FILES)

DEPEND_FILES=$(subst .o,.d,$(OBJS))

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean distclean uninstall),)
-include $(DEPEND_FILES)
endif
