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

install: default #docs
	INSTALL_DIR=$(INSTALL_DIR) LOCI_INSTALL_DIR=$(LOCI_INSTALL_DIR) \
	INSTALL_METIS=$(INSTALL_METIS) GIT_INFO=$(GIT_INFO) \
	GIT_BRANCH=$(GIT_BRANCH) bash Install.bash

ifeq ($(LOCI_INSTALL_DIR),)
LOCI_INSTALLED_DIR=$(INSTALL_DIR)
else
LOCI_INSTALLED_DIR=$(INSTALL_DIR)/$(LOCI_INSTALL_DIR)
endif

.PHONY: METIS
METIS :
ifeq ($(INSTALL_METIS),1)
	$(MAKE) -C ParMetis-4.0 LOCI_BASE="$(LOCI_BASE)"
endif


.PHONY : SPRNG
SPRNG : 
	$(MAKE) -C sprng LOCI_BASE="$(LOCI_BASE)"


.PHONY: FORCE

LOCI_MAJOR_VERSION=$(shell echo $(GIT_INFO) | sed "s/v\([0-9][0-9]*\).*/\1/")
LOCI_MINOR_VERSION=$(shell echo $(GIT_INFO) | sed "s/v[0-9][0-9]*[.]\([0-9][0-9]*\).*/\1/")
include/Config/version_info.h: FORCE
	@echo \#ifndef CONFIG_LOCI_INFO_H > version_info.h
	@echo \#define CONFIG_LOCI_INFO_H >> version_info.h
	@echo \#define LOCI_VERSION \"$(GIT_INFO)\" >> version_info.h
	@echo \#define LOCI_BRANCH \"$(GIT_BRANCH)\" >> version_info.h
	@echo \#define LOCI_VERSION_MAJOR $(LOCI_MAJOR_VERSION) >> version_info.h
	@echo \#define LOCI_VERSION_MINOR $(LOCI_MINOR_VERSION) >> version_info.h
	@echo \#endif >> version_info.h
	@# Only update version_info.h if it has changed.
	@if [ -e $@ ]; then if ! cmp version_info.h $@; then mv version_info.h $@ ; else rm version_info.h ; fi ; else mv version_info.h $@ ; fi

Tools/libTools.$(LIB_SUFFIX): include/Config/version_info.h FORCE
	$(MAKE) -C Tools LOCI_BASE="$(LOCI_BASE)" libTools.$(LIB_SUFFIX)

System/libLoci.$(LIB_SUFFIX): include/Config/version_info.h FORCE
	$(MAKE) -C System LOCI_BASE="$(LOCI_BASE)" libLoci.$(LIB_SUFFIX)

lpp: Tools/libTools.$(LIB_SUFFIX) include/Config/version_info.h FORCE
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


LS_INC_H =$(shell cd $(LOCI_BASE)/include/; ls *.h) Loci
LS_INC_LH =$(shell cd $(LOCI_BASE)/include/; ls *.lh) 

uninstall_include:
	for d in $(LS_INC_H);  do rm -f  $(LOCI_INSTALLED_DIR)/include/$$d; done
	for d in $(LS_INC_LH); do rm -f  $(LOCI_INSTALLED_DIR)/include/$$d; done
	for d in Config MPI_stubb; do rm -rf $(LOCI_INSTALLED_DIR)/include/$$d; done

uninstall_conf:
	rm -f $(LOCI_INSTALLED_DIR)/Loci.conf $(LOCI_INSTALLED_DIR)/sys.conf \
			$(LOCI_INSTALLED_DIR)/comp.conf $(LOCI_INSTALLED_DIR)/version.conf 

.PHONY: uninstall
uninstall: uninstall_include uninstall_conf
ifeq ($(INSTALL_METIS),1)
	$(MAKE) -C ParMetis-4.0 uninstall LOCI_BASE="$(LOCI_BASE)" LOCI_INSTALLED_DIR="$(LOCI_INSTALLED_DIR)"
endif
	$(MAKE) -C sprng uninstall LOCI_BASE="$(LOCI_BASE)" LOCI_INSTALLED_DIR="$(LOCI_INSTALLED_DIR)"
	$(MAKE) -C FVMMod uninstall LOCI_BASE="$(LOCI_BASE)" LOCI_INSTALLED_DIR="$(LOCI_INSTALLED_DIR)"
	$(MAKE) -C FVMAdapt uninstall LOCI_BASE="$(LOCI_BASE)" LOCI_INSTALLED_DIR="$(LOCI_INSTALLED_DIR)"
	$(MAKE) -C FVMOverset uninstall LOCI_BASE="$(LOCI_BASE)" LOCI_INSTALLED_DIR="$(LOCI_INSTALLED_DIR)"
	$(MAKE) -C FVMtools uninstall LOCI_BASE="$(LOCI_BASE)" LOCI_INSTALLED_DIR="$(LOCI_INSTALLED_DIR)"
	$(MAKE) -C lpp uninstall LOCI_BASE="$(LOCI_BASE)" LOCI_INSTALLED_DIR="$(LOCI_INSTALLED_DIR)"
	$(MAKE) -C System uninstall LOCI_BASE="$(LOCI_BASE)" LOCI_INSTALLED_DIR="$(LOCI_INSTALLED_DIR)"
	$(MAKE) -C Tools uninstall LOCI_BASE="$(LOCI_BASE)" LOCI_INSTALLED_DIR="$(LOCI_INSTALLED_DIR)"
# for d in $(LS_INC_H);  do rm -f  ${LOCI_INSTALLED_DIR}/include/$$d; done
# for d in $(LS_INC_LH); do rm -f  ${LOCI_INSTALLED_DIR}/include/$$d; done
# rm -f $(LOCI_INSTALLED_DIR)/Loci.conf $(LOCI_INSTALLED_DIR)/sys.conf $(LOCI_INSTALLED_DIR)/comp.conf 
# for d in Config MPI_stubb; do rm -rf $(LOCI_INSTALLED_DIR)/include/$$d; done

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
