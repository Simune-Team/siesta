SHELL := /bin/bash
rerun = "(There were undefined references|Rerun to get (cross-references|the bars) right)"
doc   = siesta
latex = pdflatex -shell-escape
grep = @-grep -i -s
docidx = $(doc).idx
bibtex = echo "BIBTEX"
makeidx = makeindex

pdf: $(doc)
	@echo "Done"

$(doc): $(doc).tex
	$(latex) $(doc).tex

compress:
	gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.5 \
		-dNOPAUSE -dQUIET -dBATCH \
		-sOutputFile=$(doc)-compressed.pdf $(doc).pdf

extra:
	$(bibtex) $(doc)
	$(makeidx) $(doc)

final: clean $(doc).tex
	$(MAKE) $(doc)
	$(MAKE) $(doc)
	$(MAKE) extra
	$(MAKE) $(doc)
	$(MAKE) $(doc)
	$(MAKE) extra
	$(MAKE) $(doc)

clean:
	rm -f $(doc).{lot,oct,aux,auxlock,makefile,figlist,bbl,log,ilg,lof,ind,out,blg,gls,glo,glg,idx,toc,bcf,run.xml} $(doc)-blx.bib

cleanpdf:
	rm -f $(doc).pdf

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

