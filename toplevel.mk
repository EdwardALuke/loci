
.PHONY: default install install_minimal test docs clean uninstall distclean tarball tarball-with-docs

default:
	@$(MAKE) -C ${OBJDIR} $@
install:
	@$(MAKE) -C ${OBJDIR} $@
install_minimal:
	@$(MAKE) -C ${OBJDIR} $@
test:
	@$(MAKE) -C ${OBJDIR} $@
docs:
	@$(MAKE) -C ${OBJDIR} $@
clean:
	@$(MAKE) -C ${OBJDIR} $@
uninstall:
	@$(MAKE) -C ${OBJDIR} $@
distclean:
	rm -fr ${OBJDIR} copy.out config.log
	$(MAKE) -C docs distclean
	rm -f Makefile
tarball:
	version_string=$$(sed -e "s:^GIT_INFO.*:GIT_INFO = $(GIT_INFO):g" \
	                      -e "s:^GIT_BRANCH.*:GIT_BRANCH = $(GIT_BRANCH):g" \
	                      src/conf/version.conf); \
	git archive --format=tgz --prefix=Loci-$(LOCI_REVISION_NAME)/ \
	            --add-virtual-file=Loci-$(LOCI_REVISION_NAME)/src/version.conf:"$$version_string" \
	            -o Loci-$(LOCI_REVISION_NAME).tgz HEAD
tarball-with-docs: dosc
	bash src/scripts/tarball_with_docs.sh "$(LOCI_REVISION_NAME)" "$(GIT_INFO)" "$(GIT_BRANCH)" "${OBJDIR}"
