.PHONY: default clean distclean FRC

# Clean only the docs tree where make is running. The generated top-level
# distclean target removes OBJ first, then invokes this Makefile in source docs.
default: FRC
	$(MAKE) -C tutorial   LOCI_BASE=$(LOCI_BASE)
	$(MAKE) -C doxygen    LOCI_BASE=$(LOCI_BASE)
	$(MAKE) -C developer  LOCI_BASE=$(LOCI_BASE)

clean:
	$(MAKE) -C tutorial clean  LOCI_BASE=$(LOCI_BASE)
	$(MAKE) -C doxygen clean   LOCI_BASE=$(LOCI_BASE)
	$(MAKE) -C developer clean LOCI_BASE=$(LOCI_BASE)


distclean: 
	$(MAKE) -C tutorial distclean  LOCI_BASE=$(LOCI_BASE)
	$(MAKE) -C doxygen distclean   LOCI_BASE=$(LOCI_BASE)
	$(MAKE) -C developer distclean LOCI_BASE=$(LOCI_BASE)
	rm -f *~
