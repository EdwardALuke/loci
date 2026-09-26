CASES = hexBase hexMode0 hexMode1 hexMode2 prismBase prismMode0 prismMode1 prismMode2
TOOLS = ugrid2vog marker refmesh vogcheck extract
NEEDS_H5DUMP = 1
include $(dir $(lastword $(MAKEFILE_LIST)))../common.mk

ifdef CASE_WORK
hexBase hexMode0 hexMode1 hexMode2: FIXTURE = thin_hex
hexBase hexMode0 hexMode1 hexMode2: TAG = hex_center_column.tag
prismBase prismMode0 prismMode1 prismMode2: FIXTURE = thin_prism
prismBase prismMode0 prismMode1 prismMode2: TAG = prism_lower.tag
hexBase prismBase: MODE = base
hexMode0 prismMode0: MODE = 0
hexMode1 prismMode1: MODE = 1
hexMode2 prismMode2: MODE = 2

# Both fixtures are one cell thick in z. Modes 0/1 preserve that direction;
# mode 2 splits in all directions. Cell types and coordinate planes distinguish
# the directional outcomes without requiring visualization.
$(CASES):
	cp $(INPUT)/$(FIXTURE).ugrid mesh.ugrid
	$(SERIAL) $(LOCI_BASE)/bin/ugrid2vog $(LOCI_CMD_FLAGS) -o -m mesh
	@if test "$(MODE)" = base; then cp mesh.vog $@.vog; else \
	  $(SERIAL) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g mesh.vog \
	    -ctag $(INPUT)/$(TAG) -mode $(MODE) -o refined.plan && \
	  $(SERIAL) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g mesh.vog -r refined.plan -o $@.vog; \
	fi
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) $@ > mesh.check 2>&1
	$(SERIAL) $(LOCI_BASE)/bin/extract -ascii $@ 0 x y z > mesh.coordinates
	sh $(INPUT)/mode_summary.sh $(FIXTURE) $(MODE) $@.vog mesh.check \
	  output/$@.topo mesh.coordinates > mesh.actual
	diff -u $(INPUT)/dats/$@.dat mesh.actual
endif
