SHELL     := /bin/bash

TITLEPAGE := title.tex

LATEXC    := lualatex
LXFLAGS   := -interaction=batchmode
SDFLAGS   := --natbib --latex-engine=$(LATEXC)
PRODFLAGS := --include-before-body=$(TITLEPAGE)
BIBTEX    := bibtex

SOURCES   := main.md
TARGET    := paper

BUILDDIR  := build
RESDIR    := latex-common/resources
BIBDIR    := latex-common/bibliography
SETTINGS  := latex-common/scholdoc/kma.yml

all: pdf

html: buildenv
	scholdoc $(SDFLAGS) --citeproc --output=$(BUILDDIR)/$(TARGET).html $(SOURCES)
	@firefox $(BUILDDIR)/$(TARGET).html

latex: buildenv
	scholdoc $(SETTINGS) $(SDFLAGS) $(PRODFLAGS) --output=$(BUILDDIR)/$(TARGET).tex $(SOURCES)

prev: buildenv
	scholdoc $(SETTINGS) $(SDFLAGS) --to=latex --output=$(BUILDDIR)/$(TARGET).pdf $(SOURCES)
	@evince $(BUILDDIR)/$(TARGET).pdf 2> /dev/null &

buildenv: $(BUILDDIR)
$(BUILDDIR):
	@mkdir -p $(BUILDDIR)

pdf: latex links
	@cd $(BUILDDIR) ; \
	$(LATEXC) $(LXFLAGS) -pdf $(TARGET).tex ; \
	$(BIBTEX) $(TARGET).aux ; \
	$(LATEXC) $(LXFLAGS) -pdf $(TARGET).tex ; \
	$(LATEXC) $(LXFLAGS) -pdf $(TARGET).tex
	@evince $(BUILDDIR)/$(TARGET).pdf 2> /dev/null &

links: buildenv
	@ln -sf $(shell pwd)/$(RESDIR) $(BUILDDIR)
	@ln -sf $(shell pwd)/$(BIBDIR) $(BUILDDIR)

clean:
	rm -rf $(BUILDDIR)

.PHONY: html latex prev pdf clean
