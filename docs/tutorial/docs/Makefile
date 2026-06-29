.PHONY: clean distclean FRC

TEX_FILES = 	appendix.tex \
		datatype.tex \
		intro.tex \
		io.tex \
		main.tex 


tutorial.pdf: ${TEX_FILES} 
#	 Convert source files to latex verbatim files
	cat begin.tex > Makefile_ex.tex
	cat ../heat/Makefile >> Makefile_ex.tex
	echo '\\end{verbatim}' >> Makefile_ex.tex
	cat begin.tex > entities_cc.tex
	cat ../Datatypes/entities.cc >> entities_cc.tex
	echo '\\end{verbatim}' >> entities_cc.tex
	cat begin.tex > helpers_cc.tex
	cat ../Datatypes/helpers.cc >> helpers_cc.tex
	echo '\\end{verbatim}' >> helpers_cc.tex
#	 Make document
	pdflatex  main && pdflatex  main && mv main.pdf tutorial.pdf
#	 Remove latex version of source files
	rm Makefile_ex.tex entities_cc.tex helpers_cc.tex

clean: FRC
	rm -f *.aux *.toc *.log main.pdf

distclean: clean
	rm -f *~ tutorial.pdf

