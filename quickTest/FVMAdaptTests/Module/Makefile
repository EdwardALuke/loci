TESTS = RefMesh
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
	$(MAKE) -C RefMesh LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean

distclean:
	rm -f TestResults
	$(MAKE) -C RefMesh LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean

FRC:
