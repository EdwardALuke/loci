.PHONY: clean distclean

# When invoked through the symlinked OBJ/docs tree, also find the source docs
# tree so distclean can remove generated docs packaged by tarball-with-docs.
SOURCE_DOCS_DIR := $(shell target=$$(readlink Makefile 2>/dev/null || true); \
	if [ -n "$$target" ]; then \
	  if [ "$${target#/}" != "$$target" ]; then \
	    dir=$$(dirname "$$target"); \
	  else \
	    dir=$$(cd "$$(dirname "$$target")" && pwd); \
	  fi; \
	else \
	  dir=$$(pwd); \
	fi; \
	cd "$$dir" && pwd)

clean:
	@seen=""; \
	for dir in "$(CURDIR)" "$(SOURCE_DOCS_DIR)"; do \
	  [ -d "$$dir" ] || continue; \
	  case " $$seen " in *" $$dir "*) continue ;; esac; \
	  seen="$$seen $$dir"; \
	  if [ -f "$$dir/tutorial/docs/clean" ]; then \
	    (cd "$$dir/tutorial/docs" && sh ./clean); \
	  fi; \
	done

distclean: clean
	@seen=""; \
	for dir in "$(CURDIR)" "$(SOURCE_DOCS_DIR)"; do \
	  [ -d "$$dir" ] || continue; \
	  case " $$seen " in *" $$dir "*) continue ;; esac; \
	  seen="$$seen $$dir"; \
	  rm -f "$$dir/tutorial/docs/tutorial.pdf" \
	        "$$dir/generated-docs-manifest.txt"; \
	  rm -rf "$$dir/doxygen/html" "$$dir/doxygen/latex"; \
	done
