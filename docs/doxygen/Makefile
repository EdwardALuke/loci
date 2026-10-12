.PHONY: clean distclean default FRC

DOXYGEN_DIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
BUILD_ROOT := $(abspath ../..)
SOURCE_ROOT := $(abspath $(DOXYGEN_DIR)../..)
STAGE_HTML_ASSETS := $(DOXYGEN_DIR)stage_html_assets.sh
# Figure rendering is memory intensive, so keep serial builds as the default.
DOC_JOBS ?= 1

default: html/index.html

# This is not ideal, technically all of the documentation has a dependency
# on all of the source files
html/index.html: Doxyfile FRC
	@if command -v doxygen >/dev/null 2>&1; then \
	  rm -fr html latex && \
	  doxygen Doxyfile && \
	  bash "$(STAGE_HTML_ASSETS)" "$(BUILD_ROOT)" html "$(SOURCE_ROOT)" "$(DOC_JOBS)"; \
	else \
	  echo "Skipping Doxygen API docs: doxygen not found." ; \
	fi

FRC:

clean: 

distclean: 
	rm -fr html latex *~
