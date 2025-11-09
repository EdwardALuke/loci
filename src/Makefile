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
LOCI_BASE = $(shell pwd)
LOCI_COMPILE = true
include Loci.conf
include version.conf


DEFINES = $(DEBUG) $(SYSTEM_DEFINES) $(MACHINE_SPECIFIC)

.PHONEY: default library docs install lpp FVMtools FVMMod FVMAdapt FVMOverset
.PHONEY: METIS SPRNG
default: library lpp FVMtools FVMMod FVMAdapt FVMOverset METIS SPRNG

library: Tools/libTools.$(LIB_SUFFIX) System/libLoci.$(LIB_SUFFIX)

docs:
	$(MAKE) -C Tutorial/docs LOCI_BASE="$(LOCI_BASE)"

install: default docs
	INSTALL_DIR=$(INSTALL_DIR) LOCI_INSTALL_DIR=$(LOCI_INSTALL_DIR) \
	INSTALL_METIS=$(INSTALL_METIS) GIT_INFO=$(GIT_INFO) \
	GIT_BRANCH=$(GIT_BRANCH) bash Install.bash

install_minimal: default
	INSTALL_DIR=$(INSTALL_DIR) LOCI_INSTALL_DIR=$(LOCI_INSTALL_DIR) \
	INSTALL_METIS=$(INSTALL_METIS) GIT_INFO=$(GIT_INFO) \
	GIT_BRANCH=$(GIT_BRANCH) bash Install.bash

.PHONY: METIS
METIS :
ifeq ($(INSTALL_METIS),1)
	$(MAKE) -C ParMetis-4.0 LOCI_BASE="$(LOCI_BASE)"
endif


.PHONY : SPRNG
SPRNG : 
	$(MAKE) -C sprng LOCI_BASE="$(LOCI_BASE)"


.PHONY: FORCE

LOCI_MAJOR_VERSION=$(shell echo $(GIT_INFO) | sed "s/v\([0-9][0-9]*\).*/\1/p; t; s/.*/0/")
LOCI_MINOR_VERSION=$(shell echo $(GIT_INFO) | sed "s/v[0-9][0-9]*[.]\([0-9][0-9]*\).*/\1/p; t; s/.*/0/")
include/Config/version_info.h: FORCE
	@echo \#ifndef CONFIG_LOCI_INFO_H > version_info.h
	@echo \#define CONFIG_LOCI_INFO_H >> version_info.h
	@echo \#define LOCI_VERSION_MAJOR $(LOCI_MAJOR_VERSION) >> version_info.h
	@echo \#define LOCI_VERSION_MINOR $(LOCI_MINOR_VERSION) >> version_info.h
	@echo \#endif >> version_info.h
	@# Only update version_info.h if it has changed.
	@if [ -e $@ ]; then if ! cmp version_info.h $@; then mv version_info.h $@ ; else rm version_info.h ; fi ; else mv version_info.h $@ ; fi

include/Config/version_details.h: FORCE
	@echo \#ifndef CONFIG_LOCI_DETAILS_H > version_details.h
	@echo \#define CONFIG_LOCI_DETAILS_H >> version_details.h
	@echo \#define LOCI_VERSION \"$(GIT_INFO)\" >> version_details.h
	@echo \#define LOCI_BRANCH \"$(GIT_BRANCH)\" >> version_details.h
	@echo \#endif >> version_details.h
	@# Only update version_details.h if it has changed.
	@if [ -e $@ ]; then if ! cmp version_details.h $@; then mv version_details.h $@ ; else rm version_details.h ; fi ; else mv version_details.h $@ ; fi

Tools/libTools.$(LIB_SUFFIX): include/Config/version_info.h FORCE
	$(MAKE) -C Tools LOCI_BASE="$(LOCI_BASE)" libTools.$(LIB_SUFFIX)

System/libLoci.$(LIB_SUFFIX): include/Config/version_info.h include/Config/version_details.h FORCE
	$(MAKE) -C System LOCI_BASE="$(LOCI_BASE)" libLoci.$(LIB_SUFFIX)

lpp: Tools/libTools.$(LIB_SUFFIX) include/Config/version_info.h include/Config/version_details.h FORCE
	$(MAKE) -C lpp LOCI_BASE="$(LOCI_BASE)" lpp

FVMtools: library lpp FVMAdapt SPRNG FORCE
	$(MAKE) -C FVMtools LOCI_BASE="$(LOCI_BASE)" default

FVMMod: lpp FORCE
	$(MAKE) -C FVMMod LOCI_BASE="$(LOCI_BASE)" default

FVMAdapt: lpp FORCE
	$(MAKE) -C FVMAdapt LOCI_BASE="$(LOCI_BASE)" all

FVMOverset: lpp FORCE
	$(MAKE) -C FVMOverset LOCI_BASE="$(LOCI_BASE)" all

JUNK =

.PHONY: test
test: default
	$(MAKE) -C quickTest LOCI_BASE="$(LOCI_BASE)"

.PHONY: clean

testclean:
	$(MAKE) -C quickTest distclean LOCI_BASE="$(LOCI_BASE)"

clean:
	rm -f $(JUNK)
	$(MAKE) -C Tools clean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C System clean LOCI_BASE="$(LOCI_BASE)"
ifeq ($(INSTALL_METIS),1)
	$(MAKE) -C ParMetis-4.0 clean LOCI_BASE="$(LOCI_BASE)"
endif
	$(MAKE) -C sprng clean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C FVMtools clean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C FVMMod clean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C FVMAdapt clean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C FVMOverset clean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C quickTest clean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C lpp clean LOCI_BASE="$(LOCI_BASE)"


LOCI_HEADER_FILES= $(wildcard include/*.h include/*.lh include/Config/*.h \
                              include/MPI_stubb/*.h) include/Loci

.PHONY: uninstall
uninstall: FORCE
# Remove top level header files
	@for d in $(LOCI_HEADER_FILES); do rm -f  $(LOCI_INSTALL_PATH)/$$d; done
	@for d in Config MPI_stubb; do rmdir $(LOCI_INSTALL_PATH)/include/$$d; done
# Remove Configuration files
	@rm -f $(LOCI_INSTALL_PATH)/Loci.conf $(LOCI_INSTALL_PATH)/sys.conf \
	      $(LOCI_INSTALL_PATH)/comp.conf $(LOCI_INSTALL_PATH)/version.conf 
ifeq ($(INSTALL_METIS),1)
	rm -f $(LOCI_INSTALL_PATH)/lib/libgk.$(LIB_SUFFIX) 
	rm -f $(LOCI_INSTALL_PATH)/lib/libmetis.$(LIB_SUFFIX) 
	rm -f $(LOCI_INSTALL_PATH)/lib/libparmetis.$(LIB_SUFFIX) 
endif
	$(MAKE) -C sprng uninstall LOCI_BASE="$(LOCI_BASE)" LOCI_INSTALL_PATH="$(LOCI_INSTALL_PATH)"
	$(MAKE) -C FVMMod uninstall LOCI_BASE="$(LOCI_BASE)" LOCI_INSTALL_PATH="$(LOCI_INSTALL_PATH)"
	$(MAKE) -C FVMAdapt uninstall LOCI_BASE="$(LOCI_BASE)" LOCI_INSTALL_PATH="$(LOCI_INSTALL_PATH)"
	$(MAKE) -C FVMOverset uninstall LOCI_BASE="$(LOCI_BASE)" LOCI_INSTALL_PATH="$(LOCI_INSTALL_PATH)"
	$(MAKE) -C FVMtools uninstall LOCI_BASE="$(LOCI_BASE)" LOCI_INSTALL_PATH="$(LOCI_INSTALL_PATH)"
	$(MAKE) -C lpp uninstall LOCI_BASE="$(LOCI_BASE)" LOCI_INSTALL_PATH="$(LOCI_INSTALL_PATH)"
	$(MAKE) -C System uninstall LOCI_BASE="$(LOCI_BASE)" LOCI_INSTALL_PATH="$(LOCI_INSTALL_PATH)"
	$(MAKE) -C Tools uninstall LOCI_BASE="$(LOCI_BASE)" LOCI_INSTALL_PATH="$(LOCI_INSTALL_PATH)"
# remove documentation
	rm -fr $(LOCI_INSTALL_PATH)/docs/1D-Diffusion
	rm -fr $(LOCI_INSTALL_PATH)/docs/Datatypes
	rm -fr $(LOCI_INSTALL_PATH)/docs/heat
	rm -f $(LOCI_INSTALL_PATH)/docs/tutorial.pdf

.PHONY: distclean
distclean:
	rm -f $(JUNK)
ifeq ($(INSTALL_METIS),1)
	$(MAKE) -C ParMetis-4.0 distclean LOCI_BASE="$(LOCI_BASE)"
endif
	$(MAKE) -C sprng distclean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C FVMMod distclean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C FVMAdapt distclean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C FVMOverset distclean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C quickTest distclean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C FVMtools distclean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C lpp distclean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C System distclean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C Tools distclean LOCI_BASE="$(LOCI_BASE)"
