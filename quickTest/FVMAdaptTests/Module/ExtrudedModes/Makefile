include $(LOCI_BASE)/Loci.conf
include $(LOCI_BASE)/version.conf
include $(TEST_BASE)/test.conf

LIB_PATH = $(LOCI_BASE)/lib:$(LD_LIBRARY_PATH)
DY_PATH = $(LOCI_BASE)/lib:$(DYLD_LIBRARY_PATH)
SERIAL = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) $(SERIALRUN)
MODE_BEHAVIOR = sh mode_behavior.sh

CASES = hexBase hexMode0 hexMode1 hexMode2 \
	prismBase prismMode0 prismMode1 prismMode2

.PHONY: default clean distclean FRC

default: TestResults

# Both fixtures are one cell thick in z. The hex tags select the center column
# and the prism tags select one triangle. The hex cells are also longer in x
# than y. Consequently, mode 0 splits x/y relative to the short z direction,
# mode 1 suppresses z and selects only the longer x direction, and mode 2
# splits isotropically in x/y/z. On the square prism, modes 0 and 1 both split
# in the extruded x/y plane while mode 2 also splits z.
# The coordinate-plane and cell-type counts make those directional outcomes
# part of the checked behavior rather than leaving them as visual inspection.
TestResults: thin_hex.ugrid thin_prism.ugrid \
		hex_center_column.tag prism_lower.tag mode_behavior.sh FRC
	$(SERIAL) $(LOCI_BASE)/bin/ugrid2vog $(LOCI_CMD_FLAGS) \
	  -o -m thin_hex > thin_hex.convert 2>&1
	$(SERIAL) $(LOCI_BASE)/bin/ugrid2vog $(LOCI_CMD_FLAGS) \
	  -o -m thin_prism > thin_prism.convert 2>&1
	cp thin_hex.vog hexBase.vog
	cp thin_prism.vog prismBase.vog
	$(SERIAL) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g thin_hex.vog \
	  -ctag hex_center_column.tag -mode 0 -o hexMode0.plan
	$(SERIAL) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g thin_hex.vog \
	  -ctag hex_center_column.tag -mode 1 -o hexMode1.plan
	$(SERIAL) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g thin_hex.vog \
	  -ctag hex_center_column.tag -mode 2 -o hexMode2.plan
	$(SERIAL) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g thin_prism.vog \
	  -ctag prism_lower.tag -mode 0 -o prismMode0.plan
	$(SERIAL) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g thin_prism.vog \
	  -ctag prism_lower.tag -mode 1 -o prismMode1.plan
	$(SERIAL) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g thin_prism.vog \
	  -ctag prism_lower.tag -mode 2 -o prismMode2.plan
	$(SERIAL) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g thin_hex.vog \
	  -r hexMode0.plan -o hexMode0.vog
	$(SERIAL) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g thin_hex.vog \
	  -r hexMode1.plan -o hexMode1.vog
	$(SERIAL) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g thin_hex.vog \
	  -r hexMode2.plan -o hexMode2.vog
	$(SERIAL) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g thin_prism.vog \
	  -r prismMode0.plan -o prismMode0.vog
	$(SERIAL) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g thin_prism.vog \
	  -r prismMode1.plan -o prismMode1.vog
	$(SERIAL) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g thin_prism.vog \
	  -r prismMode2.plan -o prismMode2.vog
	@for case in $(CASES); do \
	  $(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) $$case \
	    > $$case.check 2>&1 || exit $$?; \
	  $(SERIAL) $(LOCI_BASE)/bin/extract -ascii $$case 0 x y z \
	    > $$case.coordinates || exit $$?; \
	done
	$(MODE_BEHAVIOR) thin_hex base hexBase.vog hexBase.check \
	  output/hexBase.topo hexBase.coordinates > hexBase.actual
	$(MODE_BEHAVIOR) thin_hex 0 hexMode0.vog hexMode0.check \
	  output/hexMode0.topo hexMode0.coordinates > hexMode0.actual
	$(MODE_BEHAVIOR) thin_hex 1 hexMode1.vog hexMode1.check \
	  output/hexMode1.topo hexMode1.coordinates > hexMode1.actual
	$(MODE_BEHAVIOR) thin_hex 2 hexMode2.vog hexMode2.check \
	  output/hexMode2.topo hexMode2.coordinates > hexMode2.actual
	$(MODE_BEHAVIOR) thin_prism base prismBase.vog prismBase.check \
	  output/prismBase.topo prismBase.coordinates > prismBase.actual
	$(MODE_BEHAVIOR) thin_prism 0 prismMode0.vog prismMode0.check \
	  output/prismMode0.topo prismMode0.coordinates > prismMode0.actual
	$(MODE_BEHAVIOR) thin_prism 1 prismMode1.vog prismMode1.check \
	  output/prismMode1.topo prismMode1.coordinates > prismMode1.actual
	$(MODE_BEHAVIOR) thin_prism 2 prismMode2.vog prismMode2.check \
	  output/prismMode2.topo prismMode2.coordinates > prismMode2.actual
	@for case in $(CASES); do \
	  diff -u dats/$$case.dat $$case.actual || exit $$?; \
	done
	rm -rf output debug
	rm -f *.vog *.plan *.convert *.check *.quality *.coordinates *.actual
	@echo "FVMAdaptTests/Module/ExtrudedModes: PASSED!" > $@
	@cat $@

clean:
	rm -rf output debug
	rm -f TestResults *.vog *.plan *.convert *.check *.quality
	rm -f *.coordinates *.actual

distclean: clean

FRC:
