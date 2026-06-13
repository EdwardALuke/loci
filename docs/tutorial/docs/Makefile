TEX_FILES = appendix.tex \
		datatype.tex \
		intro.tex \
		io.tex \
		main.tex 


SOURCE_FILES =  entities_cc.tex \
		helpers_cc.tex \
		Makefile_ex.tex

tutorial.pdf: ${TEX_FILES} ${SOURCE_FILES}
	pdflatex  main && pdflatex  main && mv main.pdf tutorial.pdf
	

clean: FRC
	./clean

FRC:

Makefile_ex.tex: ../heat/Makefile
	cat begin.tex > Makefile_ex.tex
	cat ../heat/Makefile >> Makefile_ex.tex
	echo '\\end{verbatim}' >> Makefile_ex.tex

entities_cc.tex: ../Datatypes/entities.cc
	cat begin.tex > entities_cc.tex
	cat ../Datatypes/entities.cc >> entities_cc.tex
	echo '\\end{verbatim}' >> entities_cc.tex

helpers_cc.tex: ../Datatypes/helpers.cc
	cat begin.tex > helpers_cc.tex
	cat ../Datatypes/helpers.cc >> helpers_cc.tex
	echo '\\end{verbatim}' >> helpers_cc.tex



