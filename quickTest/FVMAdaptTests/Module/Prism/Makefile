include $(LOCI_BASE)/Loci.conf
include $(LOCI_BASE)/version.conf
include $(TEST_BASE)/test.conf

LIB_PATH = $(LOCI_BASE)/lib:$(LD_LIBRARY_PATH)
DY_PATH = $(LOCI_BASE)/lib:$(DYLD_LIBRARY_PATH)
SERIAL = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) $(SERIALRUN)
THREE_RANKS = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) $(MPI_RUN) -np 3
MESH_BEHAVIOR = sh ../mesh_behavior.sh

.PHONY: default empty-partition clean distclean FRC

default: TestResults

# Refine one prism isotropically, then apply derefinement tags to the saved cell
# plan and restore the base topology.
TestResults: prism.ugrid refine.tag derefine.tag FRC
	$(SERIAL) $(LOCI_BASE)/bin/ugrid2vog $(LOCI_CMD_FLAGS) -o -m prism > prism.convert 2>&1
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) prism > prism.check 2>&1
	$(SERIAL) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g prism.vog \
	  -ctag refine.tag -mode 2 -o prism.plan
	$(SERIAL) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g prism.vog \
	  -r prism.plan -o prismRef.vog
	$(SERIAL) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g prism.vog \
	  -r prism.plan -ctag derefine.tag -mode 2 -o prismRoundTrip.plan
	$(SERIAL) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g prism.vog \
	  -r prismRoundTrip.plan -o prismRoundTrip.vog
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) prismRef > prismRef.check 2>&1
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) prismRoundTrip > prismRoundTrip.check 2>&1
	$(MESH_BEHAVIOR) prism.vog prism.check > prism.actual
	diff -u dats/prism.dat prism.actual
	$(MESH_BEHAVIOR) prismRef.vog prismRef.check > prismRef.actual
	diff -u dats/prismRef.dat prismRef.actual
	$(MESH_BEHAVIOR) prismRoundTrip.vog prismRoundTrip.check > prismRoundTrip.actual
	diff -u dats/prismRoundTrip.dat prismRoundTrip.actual
	rm -rf output debug
	rm -f prism.vog prismRef.vog prismParallelRef.vog prismRoundTrip.vog derefine.tag
	rm -f *.plan *.convert *.check *.quality *.actual
	@echo "FVMAdaptTests/Module/Prism: PASSED!" > $@
	@cat $@

# A one-cell tag file leaves empty cell ranges on two ranks. This diagnostic
# records the intended behavior while the later plan-output path is being fixed.
empty-partition: prism.ugrid refine.tag FRC
	$(SERIAL) $(LOCI_BASE)/bin/ugrid2vog $(LOCI_CMD_FLAGS) -o -m prism > prism.convert 2>&1
	$(THREE_RANKS) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g prism.vog \
	  -ctag refine.tag -mode 2 -o prismParallel.plan
	$(THREE_RANKS) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g prism.vog \
	  -r prismParallel.plan -o prismParallelRef.vog
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) prismParallelRef > prismParallelRef.check 2>&1
	$(MESH_BEHAVIOR) prismParallelRef.vog prismParallelRef.check > prismParallelRef.actual
	diff -u dats/prismRef.dat prismParallelRef.actual

derefine.tag:
	awk 'BEGIN { for (i = 0; i < 6; ++i) print -1 }' > $@

clean:
	rm -rf output debug
	rm -f TestResults prism.vog prismRef.vog prismParallelRef.vog prismRoundTrip.vog derefine.tag
	rm -f *.plan *.convert *.check *.quality *.actual

distclean: clean

FRC:
