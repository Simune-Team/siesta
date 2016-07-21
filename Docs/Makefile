SHELL := /bin/bash
rerun = "(There were undefined references|Rerun to get (cross-references|the bars) right)"
latex = pdflatex -shell-escape

grep = @-grep -i -s
splitidx = splitindex
docidx = $(doc).idx
bibtex = bibtex
makeidx = makeindex

# Handler for creating the documentation in different
# text-sizes
SIZE ?= 11
input = "\\def\\siestatextsize{$(SIZE)pt}"


# Default build both
pdf: tbtrans siesta

all: final-screen final

doc: $(doc).tex

# Currently there are two manuals
%.pdf: %.tex FORCE
	$(latex) "$(input) \\input{$<}"

%-screen.pdf: %.tex FORCE
	$(latex) "$(input) \\def\\siestascreen{1} \\input{$<}"

# Basic targets
.PHONY: siesta siesta-screen
siesta: siesta.pdf
siesta-screen: siesta-screen.pdf

.PHONY: tbtrans tbtrans-screen
tbtrans: tbtrans.pdf
tbtrans-screen: tbtrans-screen.pdf

extra-siesta:
	$(bibtex) siesta
	$(splitidx) siesta.idx
	$(makeidx) siesta

extra-tbtrans:
	$(bibtex) tbtrans
	$(splitidx) tbtrans.idx
	$(makeidx) tbtrans

# Create final SIESTA manual
final-siesta:
	$(MAKE) siesta
	$(MAKE) siesta
	$(MAKE) extra-siesta
	$(MAKE) siesta
	$(MAKE) siesta

# Create final SIESTA manual for the screen
final-siesta-screen:
	$(MAKE) siesta-screen
	$(MAKE) siesta-screen
	$(MAKE) extra-siesta
	$(MAKE) siesta-screen
	$(MAKE) siesta-screen
	-mv siesta.pdf siesta-screen.pdf

# Create final tbtrans manual
final-tbtrans:
	$(MAKE) tbtrans
	$(MAKE) tbtrans
	$(MAKE) extra-tbtrans
	$(MAKE) tbtrans
	$(MAKE) tbtrans

# Create final tbtrans manual for the screen
final-tbtrans-screen:
	$(MAKE) tbtrans-screen
	$(MAKE) tbtrans-screen
	$(MAKE) extra-tbtrans
	$(MAKE) tbtrans-screen
	$(MAKE) tbtrans-screen
	-mv tbtrans.pdf tbtrans-screen.pdf


final: clean final-siesta final-tbtrans

final-screen: clean final-siesta-screen final-tbtrans-screen


# Clean targets
clean:
	rm -f {siesta,tbtrans}.{lot,oct,aux,auxlock,makefile,figlist,bbl,log,ilg,lof,ind,out,blg,gls,glo,glg,idx,toc,bcf,run.xml} {siesta,tbtrans}-*.[^p]*

clean-all: cleanall
cleanall: clean
	rm -f {siesta,siesta-screen,tbtrans,tbtrans-screen}.pdf


compress:
	gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.5 \
		-dNOPAUSE -dQUIET -dBATCH \
		-sOutputFile=siesta-compressed.pdf siesta.pdf

help:
	@echo "Help for creating the documentation for SIESTA"
	@echo ""
	@echo " Type:"
	@echo "   make"
	@echo " to build the documentation once."
	@echo " Type:"
	@echo "   make final"
	@echo " to build the documentation with references and everything."
	@echo ""
	@echo " The text-size may be controlled by adding SIZE=10|11|12"
	@echo " on the command line. I.e."
	@echo "   make final SIZE=10"
	@echo " will build the SIESTA manual with textsize 10pt."
	@echo ""
	@echo " The manuals may also be build in screen-compatible sizes"
	@echo " To compile for a 4:3 screen type:"
	@echo "   make final-screen"
	@echo " And the resulting manual will be suffixed with '-screen.pdf'"
	@echo ""
	@echo " Lastly, one may build all by using:"
	@echo "   make all"

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

FORCE:
