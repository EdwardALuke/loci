.PHONY: clean distclean

# Clean only the docs tree where make is running. The generated top-level
# distclean target removes OBJ first, then invokes this Makefile in source docs.

clean:
	@if [ -f tutorial/docs/clean ]; then \
	  (cd tutorial/docs && sh ./clean); \
	fi

distclean: clean
	rm -f tutorial/docs/tutorial.pdf generated-docs-manifest.txt
	rm -rf doxygen/html doxygen/latex
