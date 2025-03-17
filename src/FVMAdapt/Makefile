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

MODULE_OBJS =	quadface_rule.o  \
		create_prism_maps.o \
		create_hexcell_maps.o \
		merge_general_face.o \
		merge_general_edge.o \
		merge_hex_face.o \
		merge_prism_face.o \
		merge_mixed_face.o \
		update_general_face.o \
		update_general_edge.o \
		update_quad_face.o \
		update_prism_face.o \
		update_mixed_face.o \
		balance_general_face.o \
		balance_general_face_derefine.o \
		balance_general_edge.o \
		balance_general_edge_derefine.o \
		balance_quad_face.o \
		balance_quad_face_derefine.o \
		balance_prism_face.o \
		balance_prism_face_derefine.o \
		balance_mixed_face.o \
		balance_mixed_face_derefine.o \
		set_hex_nums.o \
		set_general_nums.o \
		set_prism_nums.o \
		set_offset.o \
		reset_hex_nums.o \
		reset_general_nums.o \
		reset_prism_nums.o \
		reset_offset.o \
		make_general_cellplan.o \
		make_hex_cellplan.o \
		make_prism_cellplan.o \
		balance_hex_cell.o \
		balance_prism_cell.o \
		balance_general_cell.o \
		cellplan_io.o \
		check_folded_face.o \
		get_node_remap.o \
		read_tag.o \
		parse_tree.o \
		get_fine_grid.o \
		get_fine_face.o \
		get_hex_fine_grid.o \
		get_hex_fine_face.o \
		get_prism_fine_grid.o \
		get_prism_fine_face.o \
		get_mixed_face_fine_grid.o \
		write_vog_module.o

LIB_OBJS  = 	build_general_cell.o \
		build_hexcell.o \
		build_prismcell.o \
		color_matrix.o \
		diamondcell.o \
		extract_hex_edge.o \
		extract_hex_face.o \
		extract_prism_face.o \
		node_edge.o \
		quadface.o \
		face.o \
		hexcell.o \
		prism.o \
		get_c1_general.o \
		get_c1_hex.o \
		get_c1_prism.o \
		read_par.o \
		transfer_fc.o \
		transfer_plan.o \
		write_vog.o  \
		classify_cell.o globals.o 

LOCAL_LIBS      = $(LIBXML2_LIBS)

JUNK =		*~ \
		core \
		ti_files \
		ii_files \
		rii_files

LIB_OBJS1=$(MODULE_OBJS:.o=_lo.o)
LIB_OBJS2=$(LIB_OBJS:.o=_lo.o)

INCLUDES=	-I$(LOCI_BASE)/include/FVMAdapt $(LIBXML2_INCLUDE) $(LIBXML2_DEFINE)

all: $(MODULE_NAME)_m.so lib$(LIB_NAME).$(LIB_SUFFIX) 

$(MODULE_NAME)_m.so: $(LIB_OBJS1)
	$(SHARED_LD) $(SHARED_LD_FLAGS) $(MODULE_NAME)_m.so $(LIB_FLAGS) $(LIB_OBJS1) $(LOCAL_LIBS)

lib$(LIB_NAME).$(LIB_SUFFIX): $(LIB_OBJS2)
	$(DYNAMIC_LD) $(DYNAMIC_LD_FLAGS) lib$(LIB_NAME).$(LIB_SUFFIX) $(LIB_FLAGS) $(LIB_OBJS2)


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
