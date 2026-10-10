TEST_FILES=cylinder.test fixedMass.test  wedge.test shock_tube.test 


include $(LOCI_BASE)/Loci.conf
include $(LOCI_BASE)/version.conf
include $(TEST_BASE)/test.conf

default: TestResults

TestResults: $(TEST_FILES) FRC
	tail -n 1 $(TEST_FILES) > TestResults; rm -f $(TEST_FILES)

FRC:

%.test : %.vars
	@FLOWTEST="$(FLOWTEST)" EXTRACT="$(EXTRACT)" NDIFF=../ndiff LOCI_BASE=$(LOCI_BASE) ./tests.sh $*.vars || true
	@if ! test -f $*.test; then echo "Test FVM/$*: FAILED!" > $*.test; fi
	@cat $*.test

clean:
	rm -fr TEST_* TestResults $(TEST_FILES)

distclean:
	rm -fr TEST_* TestResults $(TEST_FILES)


