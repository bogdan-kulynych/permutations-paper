## Prerequisites

Certain settings assume lualatex is used for PDF generation

```
sudo apt-get install lualatex
sudo apt-get install haskell-platform
```

[Sholdoc](http://scholdoc.scholarlymarkdown.com/) is used to compile Scholarly Markdown to TeX.

```
sudo cabal update
sudo cabal install scholdoc
sudo cabal install scholdoc-citeproc
```

## Build

```
# Produce build/paper.pdf
make pdf

# Produce build/paper.tex
make latex
```
