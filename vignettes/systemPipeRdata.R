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
})

## ----install, eval=FALSE-------------------------------------------------
## source("http://bioconductor.org/biocLite.R") # Sources the biocLite.R installation script
## biocLite("tgirke/systemPipeRdata", build_vignettes=TRUE, dependencies=TRUE) # Installs from github
## biocLite("systemPipeRdata") # Installs from Bioconductor once available there

## ----load_systemPipeRdata, eval=TRUE-------------------------------------
library("systemPipeRdata") # Loads the package

## ----documentation_systemPipeRdata, eval=FALSE---------------------------
## library(help="systemPipeRdata") # Lists package info
## vignette("systemPipeRdata") # Opens vignette

## ----generate_workenvir, eval=FALSE--------------------------------------
## genWorkenvir(workflow="varseq", mydirname=NULL)
## setwd("varseq")

## ----generate_workenvir_from_shell, eval=FALSE, engine="sh"--------------
## $ echo 'library(systemPipeRdata); genWorkenvir(workflow="varseq", mydirname=NULL)' | R --slave
## $ cd varseq
## $ R

## ----load_systemPipeR, eval=TRUE-----------------------------------------
library("systemPipeR") # Loads systemPipeR which needs to be installed via biocLite() from Bioconductor 

## ----documentation_systemPipeR, eval=FALSE-------------------------------
## vignette("systemPipeR", package = "systemPipeR")

## ----return_samplepaths, eval=TRUE---------------------------------------
pathList()

## ----sessionInfo---------------------------------------------------------
sessionInfo()

