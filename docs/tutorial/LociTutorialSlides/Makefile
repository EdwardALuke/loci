.PHONY: clean distclean FRC

TEX_FILES = Loci_Intro.tex

LociSlides.pdf: ${TEX_FILES} 
	pdflatex Loci_Intro &&  mv Loci_Intro.pdf LociSlides.pdf

clean: FRC
	rm -f *.aux *.toc *.log *.out *.snm *.nav

distclean: clean
	rm -f Loci_Intro.pdf LociSlides.pdf *~




