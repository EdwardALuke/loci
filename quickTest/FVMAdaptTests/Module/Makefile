TESTS = RefMesh UniformHex Prism Mixed PlanInputs
TEST_RESULTS = $(addsuffix /TestResults,$(TESTS))

.PHONY: default $(TESTS) clean distclean FRC

default: TestResults

TestResults: $(TESTS)
	cat $(TEST_RESULTS) > $@
	rm $(TEST_RESULTS)
	@cat $@

$(TESTS): FRC
	$(MAKE) -C $@ LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)

clean:
	rm -f TestResults
	@for test in $(TESTS); do \
	  $(MAKE) -C $$test LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean || exit $$?; \
	done

distclean:
	rm -f TestResults
	@for test in $(TESTS); do \
	  $(MAKE) -C $$test LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean || exit $$?; \
	done

FRC:
