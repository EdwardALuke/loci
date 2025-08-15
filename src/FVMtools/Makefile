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
include $(LOCI_BASE)/Loci.conf

COPT = $(COPTLESS)
#COPT = -g -O0 #-DDEBUG -DBOUNDS_CHECK

# List the name of your compiled program here

APPS = extract plot3d2vog make_periodic ugrid2vog msh2vog cfd++2vog fluent2vog cobalt2vog vogmerge vogcut vogcheck extruder vog2surf ccm2vog refmesh marker refine cgns2vog ugrid2cgns  cgns2ugrid cgns2ensight cgns2surf 

# List your object files here
OBJS =  extract.o extract_ensight.o extract_2dgv.o extract_ascii.o extract_tecplot.o extract_vtk.o extract_fv.o extract_cut.o extract_surf.o extract_mean.o extract_cgns.o

LOCAL_LIBS =

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~)  core 

default: $(APPS)

all: $(APPS)

DEPS = $(LOCI_BASE)/lib/libTools.${LIB_SUFFIX} \
       $(LOCI_BASE)/lib/libLoci.$(LIB_SUFFIX) \
       $(LOCI_BASE)/lib/libsprng.$(LIB_SUFFIX)

$(LOCI_BASE)/lib/libTools.${LIB_SUFFIX}:

$(LOCI_BASE)lib/libLoci.${LIB_SUFFIX}:

$(LOCI_BASE)/lib/libsprng.${LIB_SUFFIX}:

extract: $(OBJS) $(DEPS)
	$(LD) -o extract $(OBJS) $(TEC360LIB) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)
.PHONY: FORCE
libadf/libadf.a: FORCE
	$(MAKE) -C libadf libadf.a

extruder: extruder.o $(DEPS)
	$(LD) -o extruder extruder.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

vogcheck: vogCheck.o $(DEPS)
	$(LD) -o vogcheck vogCheck.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

vog2surf: vog2surf.o $(DEPS)
	$(LD) -o vog2surf vog2surf.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

extract_chem: extract_chem.o  $(DEPS)
	$(LD) -o extract_chem extract_chem.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

ugrid2vog: ugrid2vog.o vogtools.o $(DEPS)
	$(LD) -o ugrid2vog ugrid2vog.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

msh2vog: msh2vog.o vogtools.o $(DEPS)
	$(LD) -o msh2vog msh2vog.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

cfd++2vog: cfd++2vog.o vogtools.o $(DEPS)
	$(LD) -o cfd++2vog cfd++2vog.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

ugridmirror: ugridmirror.o vogtools.o $(DEPS)
	$(LD) -o ugridmirror ugridmirror.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

cobalt2vog: cobalt2vog.o vogtools.o $(DEPS)
	$(LD) -o cobalt2vog cobalt2vog.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

fluent2vog: fluent2vog.o vogtools.o $(DEPS)
	$(LD) -o fluent2vog fluent2vog.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

ccm2vog.o: ccm2vog.cc
	$(CXX) $(COPT) $(EXCEPTIONS) $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) -Ilibadf -c ccm2vog.cc

ccm2vog: ccm2vog.o vogtools.o libadf/libadf.a $(DEPS)
	 $(LD) -o ccm2vog ccm2vog.o vogtools.o libadf/libadf.a $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

vogmerge: vogmerge.o $(DEPS)
	$(LD) -o vogmerge vogmerge.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

vogcut: vogcut.o vogtools.o $(DEPS)
	$(LD) -o vogcut vogcut.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)


plot3d2vog: plot3d2vog.o vogtools.o $(DEPS)
	$(LD) -o plot3d2vog plot3d2vog.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

make_periodic: make_periodic.o predicates.o vogtools.o $(DEPS)
	$(LD) -o make_periodic make_periodic.o predicates.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

cgns2vog: cgns2vog.o hash.o vogtools.o $(DEPS)
	 $(LD) -o cgns2vog cgns2vog.o hash.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

cgns2ugrid: cgns2ugrid.o hash.o binaryio.o getargs.o  $(DEPS)
	 $(LD) -o cgns2ugrid cgns2ugrid.o hash.o binaryio.o getargs.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

cgns2ensight: cgns2ensight.o $(DEPS)
	 $(LD) -o cgns2ensight cgns2ensight.o  $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

cgns2surf: cgns2surf.o hash.o $(DEPS)
	 $(LD) -o cgns2surf cgns2surf.o hash.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)


ugrid2cgns: ugrid2cgns.o hash.o binaryio.o getargs.o $(DEPS)
	 $(LD) -o ugrid2cgns ugrid2cgns.o hash.o binaryio.o getargs.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

predicates.o: predicates.c
	$(CC) $(COPTLESS) -c $(DEFINES) predicates.c

make_periodic.o: make_periodic.cc 
	$(CXX) $(COPTLESS) -c $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) -c make_periodic.cc

refmesh: refmesh.o $(DEPS) $(LOCI_BASE)/lib/libfvmadaptfunc.${LIB_SUFFIX}
	$(LD) -o  refmesh  refmesh.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS) -lfvmadaptfunc
marker: marker.o $(DEPS) $(LOCI_BASE)/lib/libfvmadaptfunc.${LIB_SUFFIX}
	$(LD) -o  marker  marker.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS) -lfvmadaptfunc
refine: refine.o $(DEPS) $(LOCI_BASE)/lib/libfvmadaptfunc.${LIB_SUFFIX}
	$(LD) -o  refine  refine.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS) -lfvmadaptfunc


clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(APPS) $(JUNK) *.o
	$(MAKE) -C libadf clean

# intermediate files created from .loci files
LOCI_INTERMEDIATE_FILES = $(subst .loci,.$(LPP_I_SUFFIX), $(wildcard *.loci) )

uninstall:
	@for d in $(APPS);  do rm -f  $(LOCI_INSTALL_PATH)/bin/$$d; done
	@rm -f $(LOCI_INSTALL_PATH)/bin/extract_movie

distclean: clean 
	rm -fr $(DEPEND_FILES) $(LOCI_INTERMEDIATE_FILES)
	$(MAKE) -C libadf distclean

DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean distclean uninstall),)
-include $(DEPEND_FILES)
endif
