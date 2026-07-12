include $(LOCI_BASE)/Loci.conf
include $(LOCI_BASE)/version.conf
include $(TEST_BASE)/test.conf

LIB_PATH = $(LOCI_BASE)/lib:$(LD_LIBRARY_PATH)
DY_PATH = $(LOCI_BASE)/lib:$(DYLD_LIBRARY_PATH)
SERIAL = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) $(SERIALRUN)
MESH_BEHAVIOR = sh ../mesh_behavior.sh

.PHONY: default clean distclean FRC

default: TestResults

# Refine each cell type independently to exercise adaptation across their
# shared prism/hex face.
TestResults: mixed.ugrid refine_prism.tag refine_hex.tag FRC
	$(SERIAL) $(LOCI_BASE)/bin/ugrid2vog $(LOCI_CMD_FLAGS) -o -m mixed > mixed.convert 2>&1
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) mixed > mixed.check 2>&1
	$(SERIAL) $(LOCI_BASE)/bin/refine $(LOCI_CMD_FLAGS) -g mixed.vog \
	  -ctag refine_prism.tag -mode 2 -o mixedRefPrism.vog
	$(SERIAL) $(LOCI_BASE)/bin/refine $(LOCI_CMD_FLAGS) -g mixed.vog \
	  -ctag refine_hex.tag -mode 2 -o mixedRefHex.vog
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) mixedRefPrism > mixedRefPrism.check 2>&1
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) mixedRefHex > mixedRefHex.check 2>&1
	$(MESH_BEHAVIOR) mixed.vog mixed.check > mixed.actual
	diff -u dats/mixed.dat mixed.actual
	$(MESH_BEHAVIOR) mixedRefPrism.vog mixedRefPrism.check > mixedRefPrism.actual
	diff -u dats/mixedRefPrism.dat mixedRefPrism.actual
	$(MESH_BEHAVIOR) mixedRefHex.vog mixedRefHex.check > mixedRefHex.actual
	diff -u dats/mixedRefHex.dat mixedRefHex.actual
	rm -rf output debug
	rm -f mixed.vog mixedRefPrism.vog mixedRefHex.vog *.plan
	rm -f *.convert *.check *.quality *.actual
	@echo "FVMAdaptTests/Module/Mixed: PASSED!" > $@
	@cat $@

clean:
	rm -rf output debug
	rm -f TestResults mixed.vog mixedRefPrism.vog mixedRefHex.vog *.plan
	rm -f *.convert *.check *.quality *.actual

distclean: clean

FRC:
