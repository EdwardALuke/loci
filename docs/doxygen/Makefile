.PHONY: clean distclean default FRC

DOXYGEN_DIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
STAGE_HTML_ASSETS := $(DOXYGEN_DIR)stage_html_assets.sh

default: html/index.html

# This is not ideal, technically all of the documentation has a dependency
# on all of the source files
html/index.html: Doxyfile FRC
	@if command -v doxygen >/dev/null 2>&1; then \
	  rm -fr html latex && \
	  doxygen Doxyfile && \
	  bash "$(STAGE_HTML_ASSETS)" ../.. html; \
	else \
	  echo "Skipping Doxygen API docs: doxygen not found." ; \
	fi

FRC:

clean: 

distclean: 
	rm -fr html latex *~
