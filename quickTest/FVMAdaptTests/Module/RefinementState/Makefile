include $(LOCI_BASE)/Loci.conf
include $(LOCI_BASE)/version.conf
include $(TEST_BASE)/test.conf

TARGET = test_refinement_state.exe
FILES = test_refinement_state.loci
OBJS = $(call loci_file2objs,$(FILES))
LOCI_INTERMEDIATE_FILES = $(FILES:.loci=.$(LPP_I_SUFFIX))

SOURCE_GRID = ../Prism/prism.ugrid
LIB_PATH = $(LOCI_BASE)/lib:$(LD_LIBRARY_PATH)
DY_PATH = $(LOCI_BASE)/lib:$(DYLD_LIBRARY_PATH)
SERIAL = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) $(SERIALRUN)
THREE_RANKS = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) $(MPI_RUN) -np 3

.PHONY: default clean distclean FRC

default: TestResults

$(TARGET): $(OBJS) $(LOCI_BASE)/lib/libfvmadaptfunc.$(LIB_SUFFIX)
	$(LD) -o $@ $(OBJS) -L$(LOCI_BASE)/lib -lfvmadaptfunc $(LIBS) $(LDFLAGS)

prism.ugrid: $(SOURCE_GRID)
	cp $< $@

prism.vog: prism.ugrid
	$(SERIAL) $(LOCI_BASE)/bin/ugrid2vog $(LOCI_CMD_FLAGS) -o -m prism \
	  > prism.convert 2>&1

# Exercise two online cycles in serial and the first-cycle handoff with empty
# source-cell MPI partitions. A scheduled rule consumes all three public facts
# over geom_cells.
TestResults: $(TARGET) prism.vog FRC
	$(SERIAL) ./$(TARGET) --scheduleoutput prism > serial.log 2>&1
	$(THREE_RANKS) ./$(TARGET) --scheduleoutput prism > parallel.log 2>&1
	rm -rf output debug
	rm -f prism.ugrid prism.vog prism.convert serial.log parallel.log
	rm -f currentPlan nextPlan *.plan *.quality
	@echo "FVMAdaptTests/Module/RefinementState: PASSED!" > $@
	@cat $@

clean:
	rm -rf output debug
	rm -f TestResults $(TARGET) $(OBJS) $(LOCI_INTERMEDIATE_FILES)
	rm -f prism.ugrid prism.vog prism.convert serial.log parallel.log
	rm -f currentPlan nextPlan *.plan *.quality

distclean: clean
	rm -f $(OBJS:.o=.d)

FRC:

DEPEND_FILES = $(OBJS:.o=.d)

ifeq ($(filter $(MAKECMDGOALS),clean distclean),)
-include $(DEPEND_FILES)
endif
