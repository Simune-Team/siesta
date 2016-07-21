SHELL := /bin/bash
rerun = "(There were undefined references|Rerun to get (cross-references|the bars) right)"
latex = pdflatex -shell-escape
doc = siesta
grep = @-grep -i -s
splitidx = splitindex
docidx = $(doc).idx
bibtex = bibtex
makeidx = makeindex

# Default build both
pdf: siesta tbtrans
	@echo "Done"

# Generic document compilation
doc: $(doc).tex
	$(latex) $(doc).tex

siesta: doc=siesta
siesta: doc

tbtrans: doc=tbtrans
tbtrans: doc

compress:
	gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.5 \
		-dNOPAUSE -dQUIET -dBATCH \
		-sOutputFile=$(doc)-compressed.pdf $(doc).pdf

extra:
	$(bibtex) $(doc)
	$(splitidx) $(docidx)
	$(makeidx) $(doc)

# Generic creation of the final documentation
final-doc: $(doc).tex clean
	$(latex) $(doc).tex
	$(latex) $(doc).tex
	$(bibtex) $(doc)
	$(splitidx) $(docidx)
	$(makeidx) $(doc)
	$(latex) $(doc).tex
	$(latex) $(doc).tex

final-siesta: doc=siesta
final-siesta: final-doc

final-tbtrans: doc=tbtrans
final-tbtrans: final-doc

final: clean final-siesta final-tbtrans

clean:
	rm -f {siesta,tbtrans}.{lot,oct,aux,auxlock,makefile,figlist,bbl,log,ilg,lof,ind,out,blg,gls,glo,glg,idx,toc,bcf,run.xml} {siesta,tbtrans}-*.*

cleanpdf:
	rm -f {siesta,tbtrans}.pdf

check: citations references badness undefined underfull overfull warning

.PHONY: citations
citations:
	@echo "***************"
	@echo Finds citations
	$(grep) "citation" $(doc).log

.PHONY: references
references:
	@echo "***************"
	@echo Finds References
	$(grep) "reference" $(doc).log

.PHONY: warning
warning:
	@echo "***************"
	@echo Finds Warnings
	$(grep) "warning" $(doc).log

.PHONY: badness
badness:
	@echo "***************"
	@echo Finds badness
	$(grep) "badness" $(doc).log

.PHONY: undefined
undefined:
	@echo "***************"
	@echo Finds undefined
	$(grep) "undefined" $(doc).log

.PHONY: underfull
underfull:
	@echo "***************"
	@echo Finds underfull
	$(grep) "underfull" $(doc).log

.PHONY: overfull
overfull:
	@echo "***************"
	@echo Finds overfull
	$(grep) "overfull" $(doc).log

