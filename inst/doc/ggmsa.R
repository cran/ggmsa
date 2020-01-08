## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

CRANpkg <- function(pkg) {
    cran <- "https://cran.r-project.org/package"
    fmt <- "[%s](%s=%s)"
    sprintf(fmt, pkg, cran, pkg)
}

Biocpkg <- function(pkg) {
    sprintf("[%s](http://bioconductor.org/packages/%s)", pkg, pkg)
}

## ----echo=FALSE, results='hide', message=FALSE--------------------------------
library(ggmsa)

## ----eval=FALSE---------------------------------------------------------------
#  ## installing the package
#  install.packages("ggmsa")
#  ## ## loading the package
#  library("ggmsa")

## ----warning=FALSE------------------------------------------------------------
sequences <- system.file("extdata", "sample.fasta", package = "ggmsa")
print(sequences)

## ----fig.height = 2, fig.width = 10, warning=FALSE----------------------------
ggmsa(sequences, 320, 360, color = "Clustal")

## ----fig.height = 2, fig.width = 10, warning=FALSE----------------------------
ggmsa(sequences, 320, 360, color = "Chemistry_AA")

## ----fig.height = 2, fig.width = 10, warning=FALSE----------------------------
ggmsa(sequences, 320, 360, color = "Shapely_AA")

## ----fig.height = 2, fig.width = 10, warning=FALSE----------------------------
ggmsa(sequences, 320, 360, color = "Taylor_AA")

## ----fig.height = 2, fig.width = 10, warning=FALSE----------------------------
ggmsa(sequences, 320, 360, color = "Zappo_AA")

## ----fig.height = 2, fig.width = 10, warning=FALSE----------------------------
ggmsa(sequences, 320, 360, font = NULL, color = "Chemistry_AA")

## ----fig.height=2, fig.width=10, message=FALSE, warning=FALSE-----------------
library(Biostrings)
x <- readAAStringSet(sequences)
d <- as.dist(stringDist(x, method = "hamming")/width(x)[1])
library(ape)
tree <- bionj(d)
library(ggtree)
p <- ggtree(tree ) + geom_tiplab()

data = tidy_msa(x, 164, 213)
p + geom_facet(geom = geom_msa, data = data,  panel = 'msa',
               font = NULL, color = "Chemistry_AA") +
    xlim_tree(1)

## ----fig.height=6, fig.width=13, message=FALSE, warning=FALSE-----------------
f <- system.file("extdata", "LeaderRepeat_All.fa", package = "ggmsa")
s <- readDNAStringSet(f)
strings <- as.character(s)
p1 <- ggmsa(s, font = NULL, color = 'Chemistry_NT')

library(ggseqlogo)
library(cowplot)

p2 <- axis_canvas(p1, axis='x')+ geom_logo(strings, 'probability') 
pp <- insert_xaxis_grob(p1, p2, position="top", grid::unit(.05, "null"))
ggdraw(pp)

