include $(LOCI_BASE)/Loci.conf
include $(LOCI_BASE)/version.conf
include $(TEST_BASE)/test.conf

LIB_PATH = $(LOCI_BASE)/lib:$(LD_LIBRARY_PATH)
DY_PATH = $(LOCI_BASE)/lib:$(DYLD_LIBRARY_PATH)
RUN = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH)
MESH_BEHAVIOR = sh ../mesh_behavior.sh

XML_GRID = ../RefMesh/testGrid.vog
PAR_GRID = $(TEST_BASE)/FVMModUnitTests/tests/block.vog

.PHONY: default clean distclean FRC

default: TestResults

# Exercise an XML refinement region and a parameter-file spacing source. The
# checks use final mesh behavior rather than depending on plan-file encoding.
TestResults: $(XML_GRID) $(PAR_GRID) region.xml sources.par FRC
	$(RUN) $(MPIRUN) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) \
	  -g $(XML_GRID) -xml region.xml -tol 0.08 -mode 2 -o xml_region.plan
	$(RUN) $(MPIRUN) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) \
	  -g $(XML_GRID) -r xml_region.plan -o xml_region.vog > xml_region.log 2>&1
	$(RUN) $(MPIRUN) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) \
	  -g $(PAR_GRID) -par sources.par -mode 2 -o spacing_source.plan
	$(RUN) $(MPIRUN) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) \
	  -g $(PAR_GRID) -r spacing_source.plan -o spacing_source.vog > spacing_source.log 2>&1
	$(RUN) $(MPIRUN) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) xml_region > xml_region.check 2>&1
	$(RUN) $(MPIRUN) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) spacing_source > spacing_source.check 2>&1
	$(MESH_BEHAVIOR) xml_region.vog xml_region.check > xml_region.actual
	diff -u dats/xml_region.dat xml_region.actual
	$(MESH_BEHAVIOR) spacing_source.vog spacing_source.check > spacing_source.actual
	diff -u dats/spacing_source.dat spacing_source.actual
	rm -fr *.plan *.vog *.log *.check *.quality *.actual output debug
	@echo "FVMAdaptTests/Module/PlanInputs: PASSED!" > $@
	@cat $@

clean:
	rm -fr TestResults *.plan *.vog *.log *.check *.quality *.actual output debug

distclean: clean

FRC:
