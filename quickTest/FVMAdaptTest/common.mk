# Group Makefiles name their CASES and provide their recipes inside an
# `ifdef CASE_WORK` guard. Run each recipe in work/<case>, never beside inputs.
INPUT := $(abspath $(dir $(firstword $(MAKEFILE_LIST))))
GROUP := $(notdir $(INPUT))
TEST_BASE ?= $(abspath $(INPUT)/../..)
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
TIMEOUT ?= timeout
TIME_LIMIT ?= 120s
export H5DUMP
SERIAL = $(TIMEOUT) --kill-after=10s $(TIME_LIMIT) $(SERIALRUN)
PARALLEL = $(TIMEOUT) --kill-after=10s $(TIME_LIMIT) $(MPIRUN)
THREE_RANKS = $(TIMEOUT) --kill-after=10s $(TIME_LIMIT) $(MPI_RUN) -np 3
MESH_SUMMARY = sh $(INPUT)/../mesh_summary.sh

ifndef CASE_WORK
check:
	@test -f "$(LOCI_BASE)/Loci.conf" || { echo 'Set LOCI_BASE to a built or installed Loci tree.' >&2; exit 1; }
	@for tool in $(TOOLS); do test -x "$(LOCI_BASE)/bin/$$tool" || { echo "Missing tool: $(LOCI_BASE)/bin/$$tool" >&2; exit 1; }; done
	@for cmd in $(firstword $(MPI_RUN)) $(TIMEOUT) $(if $(NEEDS_H5DUMP),$(H5DUMP)) $(if $(SOURCES),$(firstword $(CXX))); do \
	  command -v "$$cmd" >/dev/null || { echo "Missing command: $$cmd" >&2; exit 1; }; done
	@printf '%s\n' '$(GROUP): Loci=$(LOCI_BASE), modules=$(LOCI_MODULE_PATH)'

TestResults:
	@rm -f $@
	@$(MAKE) --no-print-directory check
	@status=0; for case in $(CASES); do \
	  $(MAKE) --no-print-directory $$case || status=1; \
	  cat work/$$case/TestResults >> $@; \
	done; exit $$status

$(CASES): check
	@rm -rf work/$@
	@mkdir -p work/$@
	@status=0; $(MAKE) --no-print-directory -C work/$@ -f "$(INPUT)/Makefile" CASE_WORK=1 $@ \
	  > work/$@/run.log 2>&1 || status=$$?; \
	if test $$status -eq 0; then result=PASSED; else result=FAILED; fi; \
	echo "FVMAdaptTest/$(GROUP)/$@: $$result" > work/$@/TestResults; \
	cat work/$@/TestResults; \
	if test $$status -ne 0; then \
	  tail -n 35 work/$@/run.log; \
	  echo "Log: $(INPUT)/work/$@/run.log"; \
	fi; exit $$status
else
INCLUDES = -I$(TEST_BASE)/contrib/doctest -I$(LOCI_BASE)/include
%.o: $(INPUT)/%.cc
	$(CXX) $(COPT) $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) -c $< -o $@
endif
endif
