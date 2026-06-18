.PHONY: clean distclean default


default: html/index.html

# This is not ideal, technically all of the documentation has a dependency
# on all of the source files
html/index.html: Doxyfile
	@if command -v doxygen >/dev/null 2>&1; then \
	  rm -fr html latex; \
	  doxygen Doxyfile;\
	else \
	  echo "Skipping Doxygen API docs: doxygen not found." ; \
	fi

clean: 

distclean: 
	rm -fr html latex
