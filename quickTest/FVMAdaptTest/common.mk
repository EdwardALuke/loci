# Group Makefiles name their CASES and provide their recipes inside an
# `ifdef CASE_WORK` guard. Run each recipe in work/<case>, never beside inputs.
GROUP := $(notdir $(abspath $(dir $(firstword $(MAKEFILE_LIST)))))
INPUT := $(TEST_BASE)/FVMAdaptTest/$(GROUP)

.DEFAULT_GOAL := TestResults
.PHONY: TestResults check clean distclean $(CASES)

clean distclean:
	rm -rf work
	rm -f TestResults

ifeq ($(filter clean distclean,$(MAKECMDGOALS)),)
ifneq ($(wildcard $(LOCI_BASE)/Loci.conf),)
include $(LOCI_BASE)/Loci.conf
include $(TEST_BASE)/test.conf
endif

# Use the same build for modules and tools, not a stale shell setting.
# A make-command-line override is still available for comparison runs.
export LOCI_MODULE_PATH = $(LOCI_BASE)/lib
H5DUMP ?= h5dump
TIMEOUT_SECONDS ?= 120
export H5DUMP
LOCI_CMD_FLAGS += --timeout $(TIMEOUT_SECONDS)
SERIAL = LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(LOCI_BASE)/lib DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$(LOCI_BASE)/lib $(SERIALRUN)
PARALLEL = LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(LOCI_BASE)/lib DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$(LOCI_BASE)/lib $(MPIRUN)
THREE_RANKS = LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(LOCI_BASE)/lib DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$(LOCI_BASE)/lib $(MPI_RUN) -np 3
MESH_SUMMARY = sh $(TEST_BASE)/FVMAdaptTest/mesh_summary.sh

ifndef CASE_WORK
check:
	@test -f "$(LOCI_BASE)/Loci.conf" || { echo 'Set LOCI_BASE to a built or installed Loci tree.' >&2; exit 1; }
	@for tool in $(TOOLS); do test -x "$(LOCI_BASE)/bin/$$tool" || { echo "Missing tool: $(LOCI_BASE)/bin/$$tool" >&2; exit 1; }; done
	@for cmd in $(firstword $(MPI_RUN)) $(if $(NEEDS_H5DUMP),$(H5DUMP)) $(if $(SOURCES),$(firstword $(CXX))); do \
	  command -v "$$cmd" >/dev/null || { echo "Missing command: $$cmd" >&2; exit 1; }; done
	@printf '%s\n' '$(GROUP): Loci=$(LOCI_BASE), modules=$(LOCI_MODULE_PATH)'

CASE_RESULTS := $(foreach c,$(CASES),work/$c/TestResults)
.PHONY: FORCE

TestResults: $(CASE_RESULTS)
	@rm -f $@
	@cat $^ > $@

$(CASES): %: work/%/TestResults
	@:

$(CASE_RESULTS): work/%/TestResults: FORCE | check
	@rm -rf work/$*
	@mkdir -p work/$*
	@status=0; $(MAKE) --no-print-directory -C work/$* -f "$(INPUT)/Makefile" CASE_WORK=1 $* \
	  > work/$*/run.log 2>&1 || status=$$?; \
	if test $$status -eq 0; then result=PASSED; else result=FAILED; fi; \
	echo "FVMAdaptTest/$(GROUP)/$*: $$result" > $@; \
	cat $@; \
	if test $$status -ne 0; then \
	  tail -n 35 work/$*/run.log; \
	  echo "Log: $(INPUT)/work/$*/run.log"; \
	fi
else
INCLUDES = -I$(TEST_BASE)/contrib/doctest -I$(LOCI_BASE)/include
%.o: $(INPUT)/%.cc
	$(CXX) $(COPT) $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) -c $< -o $@
endif
endif
