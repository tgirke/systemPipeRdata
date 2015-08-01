## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()
options(width=100, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")))

## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE-------------------
suppressPackageStartupMessages({
    library(systemPipeR)
    library(BiocGenerics)
    library(S4Vectors)
})

## ----install, eval=FALSE-------------------------------------------------
## source("http://bioconductor.org/biocLite.R") # Sources the biocLite.R installation script
## biocLite("systemPipeR") # Installs systemPipeR from Bioconductor
## biocLite("tgirke/systemPipeRdata", build_vignettes=TRUE, dependencies=TRUE) # From github

## ----generate_workenvir, eval=FALSE--------------------------------------
## genWorkenvir(workflow="varseq")
## setwd("varseq")

## ----run_workflow, eval=FALSE--------------------------------------------
## vignette("systemPipeR", package = "systemPipeR")

## ----sessionInfo---------------------------------------------------------
sessionInfo()

