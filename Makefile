DOCTYPE = DMTN
DOCNUMBER = 069
DOCNAME = $(DOCTYPE)-$(DOCNUMBER)

$(DOCNAME).pdf: $(DOCNAME).tex
	latexmk -bibtex -xelatex $(DOCNAME)
