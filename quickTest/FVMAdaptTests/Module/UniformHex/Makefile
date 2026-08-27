include $(LOCI_BASE)/Loci.conf
include $(LOCI_BASE)/version.conf
include $(TEST_BASE)/test.conf

SOURCE_GRID = ../../../FVMModUnitTests/tests/block.vog
LIB_PATH = $(LOCI_BASE)/lib:$(LD_LIBRARY_PATH)
DY_PATH = $(LOCI_BASE)/lib:$(DYLD_LIBRARY_PATH)
SERIAL = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) $(SERIALRUN)
THREE_RANKS = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) $(MPI_RUN) -np 3
MESH_BEHAVIOR = sh ../mesh_behavior.sh

.PHONY: default clean distclean FRC

default: TestResults

# Exercise full isotropic refinement and a derefinement restart to the base
# mesh. With 64 tags over three ranks, the final tag must still refine exactly
# one cell.
TestResults: block.vog full_refine.tag derefine.tag final_cell_only.tag FRC
	$(SERIAL) $(LOCI_BASE)/bin/refine $(LOCI_CMD_FLAGS) -g block.vog \
	  -ctag full_refine.tag -mode 2 -o blockRef.vog -oplan block.plan
	$(SERIAL) $(LOCI_BASE)/bin/refine $(LOCI_CMD_FLAGS) -g block.vog \
	  -r block.plan -ctag derefine.tag -mode 2 \
	  -o blockRoundTrip.vog -oplan blockRoundTrip.plan
	$(THREE_RANKS) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g block.vog \
	  -ctag final_cell_only.tag -mode 2 -o finalCell.plan
	$(THREE_RANKS) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g block.vog \
	  -r finalCell.plan -o finalCellRef.vog
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) blockRef > blockRef.check 2>&1
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) blockRoundTrip > blockRoundTrip.check 2>&1
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) finalCellRef > finalCellRef.check 2>&1
	$(MESH_BEHAVIOR) blockRef.vog blockRef.check > blockRef.actual
	diff -u dats/blockRef.dat blockRef.actual
	$(MESH_BEHAVIOR) blockRoundTrip.vog blockRoundTrip.check > blockRoundTrip.actual
	diff -u dats/blockRoundTrip.dat blockRoundTrip.actual
	$(MESH_BEHAVIOR) finalCellRef.vog finalCellRef.check > finalCellRef.actual
	diff -u dats/finalCellRef.dat finalCellRef.actual
	rm -rf output debug
	rm -f block.vog *.plan *.tag *Ref.vog blockRoundTrip.vog
	rm -f *.check *.quality *.actual
	@echo "FVMAdaptTests/Module/UniformHex: PASSED!" > $@
	@cat $@

block.vog: $(SOURCE_GRID)
	cp $< $@

full_refine.tag:
	awk 'BEGIN { for (i = 0; i < 64; ++i) print 1 }' > $@

derefine.tag:
	awk 'BEGIN { for (i = 0; i < 512; ++i) print -1 }' > $@

final_cell_only.tag:
	awk 'BEGIN { for (i = 1; i < 64; ++i) print 0; print 1 }' > $@

clean:
	rm -rf output debug
	rm -f TestResults block.vog *.plan *.tag *Ref.vog blockRoundTrip.vog
	rm -f *.check *.quality *.actual

distclean: clean

FRC:
