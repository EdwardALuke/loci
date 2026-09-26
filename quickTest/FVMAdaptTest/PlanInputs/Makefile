CASES = xml parameter
TOOLS = marker refmesh vogcheck
NEEDS_H5DUMP = 1
include $(dir $(lastword $(MAKEFILE_LIST)))../common.mk

ifdef CASE_WORK
# Refine the unit-cube mesh inside an XML region using a spacing tolerance.
xml:
	$(PARALLEL) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) \
	  -g $(INPUT)/../RefMesh/testGrid.vog -xml $(INPUT)/region.xml -tol 0.08 -mode 2 -o region.plan
	$(PARALLEL) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) \
	  -g $(INPUT)/../RefMesh/testGrid.vog -r region.plan -o region.vog
	$(PARALLEL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) region > region.check 2>&1
	$(MESH_SUMMARY) region.vog region.check > region.actual
	diff -u $(INPUT)/dats/xml_region.dat region.actual

# Refine a block mesh using the point source in the parameter file.
parameter:
	$(PARALLEL) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) \
	  -g $(TEST_BASE)/FVMModUnitTests/tests/block.vog -par $(INPUT)/sources.par -mode 2 -o source.plan
	$(PARALLEL) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) \
	  -g $(TEST_BASE)/FVMModUnitTests/tests/block.vog -r source.plan -o source.vog
	$(PARALLEL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) source > source.check 2>&1
	$(MESH_SUMMARY) source.vog source.check > source.actual
	diff -u $(INPUT)/dats/spacing_source.dat source.actual
endif
