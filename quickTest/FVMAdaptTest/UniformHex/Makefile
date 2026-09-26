CASES = round_trip final_cell_tag
TOOLS = refine marker refmesh vogcheck
NEEDS_H5DUMP = 1
include $(dir $(lastword $(MAKEFILE_LIST)))../common.mk

ifdef CASE_WORK
block.vog:
	cp $(TEST_BASE)/FVMModUnitTests/tests/block.vog $@

# Refine all 64 cells, then derefine the saved plan back to the base mesh.
round_trip: block.vog
	awk 'BEGIN { for (i=0; i<64; ++i) print 1 }' > refine.tag
	awk 'BEGIN { for (i=0; i<512; ++i) print -1 }' > derefine.tag
	$(SERIAL) $(LOCI_BASE)/bin/refine $(LOCI_CMD_FLAGS) -g block.vog \
	  -ctag refine.tag -mode 2 -o refined.vog -oplan refined.plan
	$(SERIAL) $(LOCI_BASE)/bin/refine $(LOCI_CMD_FLAGS) -g block.vog \
	  -r refined.plan -ctag derefine.tag -mode 2 -o restored.vog -oplan restored.plan
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) refined > refined.check 2>&1
	$(MESH_SUMMARY) refined.vog refined.check > refined.actual
	diff -u $(INPUT)/dats/blockRef.dat refined.actual
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) restored > restored.check 2>&1
	$(MESH_SUMMARY) restored.vog restored.check > restored.actual
	diff -u $(INPUT)/dats/blockRoundTrip.dat restored.actual

# 64 tags divide as 21/21/22 across three ranks. The final tag must not be lost:
# replacing the final hex by eight children increases the cell count to 71.
final_cell_tag: block.vog
	awk 'BEGIN { for (i=1; i<64; ++i) print 0; print 1 }' > final.tag
	$(THREE_RANKS) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g block.vog \
	  -ctag final.tag -mode 2 -o final.plan
	$(THREE_RANKS) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g block.vog -r final.plan -o final.vog
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) final > final.check 2>&1
	$(MESH_SUMMARY) final.vog final.check > final.actual
	diff -u $(INPUT)/dats/finalCellRef.dat final.actual
endif
