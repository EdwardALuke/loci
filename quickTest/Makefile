TEST_BASE ?= $(shell pwd)

SUITES = Containers Tools FVMModUnitTests FVM SystemTests FVMAdaptTest
SUITE_RESULTS := $(addsuffix /TestResults,$(SUITES))
CLEAN_SUITES := $(addprefix clean-,$(SUITES))

.DEFAULT_GOAL := default
DISTCLEAN_SUITES := $(addprefix distclean-,$(SUITES))
.PHONY: default test $(SUITES) $(CLEAN_SUITES) $(DISTCLEAN_SUITES) clean distclean FORCE

default: TestResults
	@echo ''; echo '========== quickTest summary =========='
	@grep PASS TestResults
	@grep FAIL TestResults || true
	@echo '======================================='
	@! grep -q FAIL TestResults

test: default

TestResults: $(SUITE_RESULTS)
	@cat $^ > $@

# Each suite sub-make is allowed to fail; TestResults is still produced.
$(SUITE_RESULTS): %/TestResults: FORCE
	@$(MAKE) --no-print-directory -C $* $(if $(filter FVMModUnitTests,$*),tests) \
	  LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) \
	  || echo "$*: FAILED (make error)" > $@

FORCE:

# FVMModUnitTests uses 'tests' as its entry target (builds the binary first).
FVMModUnitTests: FVMModUnitTests/TestResults
	@:

$(filter-out FVMModUnitTests,$(SUITES)): %: %/TestResults
	@:

clean: $(CLEAN_SUITES)
	rm -f TestResults

$(CLEAN_SUITES): clean-%:
	@$(MAKE) --no-print-directory -C $* LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean

distclean: $(DISTCLEAN_SUITES)
	rm -f TestResults

$(DISTCLEAN_SUITES): distclean-%:
	@$(MAKE) --no-print-directory -C $* LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean
