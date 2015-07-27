## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"---------------------------------
BiocStyle::latex(use.unsrturl=FALSE)

## ----setup, include=FALSE, cache=FALSE-------------------------------------------------------
library(knitr)
# set global chunk options for knitr
opts_chunk$set(comment=NA, warning=FALSE, message=FALSE, fig.path='figure/systemPipeR-')
options(formatR.arrow=TRUE, width=95)
unlink("test.db")

## ----eval=FALSE------------------------------------------------------------------------------
#  source("http://bioconductor.org/biocLite.R") # Sources the biocLite.R installation script
#  biocLite("systemPipeR") # Installs systemPipeR from Bioconductor
#  biocLite("tgirke/systemPipeRdata", build_vignettes=TRUE, dependencies=TRUE) # From github

## ----eval=FALSE------------------------------------------------------------------------------
#  genWorkenvir(workflow="varseq")
#  setwd("varseq")

## ----eval=FALSE------------------------------------------------------------------------------
#  vignette("systemPipeR", package = "systemPipeR")

## ----sessionInfo, results='asis'-------------------------------------------------------------
toLatex(sessionInfo())

