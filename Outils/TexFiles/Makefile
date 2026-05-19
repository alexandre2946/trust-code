FILE3=mini_prm
FILE4=mini_jupyter


default:  $(FILE3)

jupyter: $(FILE4)

$(FILE3): $(FILE3).pdf
$(FILE3).pdf: $(FILE3).tex Makefile 
	pdflatex -interaction=nonstopmode $(FILE3).tex > LOG

$(FILE4): $(FILE4).pdf
$(FILE4).pdf: $(FILE4).tex Makefile 
	xelatex -interaction=nonstopmode $(FILE4).tex > LOG
