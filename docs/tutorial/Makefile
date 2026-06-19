.PHONY: default clean distclean FRC

PDF_FILES = docs/tutorial.pdf LociTutorialSlides/LociSlides.pdf

default: ${PDF_FILES}

docs/tutorial.pdf: FRC
	$(MAKE) -C docs LOCI_BASE=$(LOCI_BASE)

LociTutorialSlides/LociSlides.pdf: FRC
	$(MAKE) -C LociTutorialSlides LOCI_BASE=$(LOCI_BASE)

clean:
	$(MAKE) -C docs clean LOCI_BASE=$(LOCI_BASE)
	$(MAKE) -C LociTutorialSlides clean LOCI_BASE=$(LOCI_BASE)

distclean:
	$(MAKE) -C docs distclean LOCI_BASE=$(LOCI_BASE)
	$(MAKE) -C LociTutorialSlides distclean LOCI_BASE=$(LOCI_BASE)
	rm -f *~
