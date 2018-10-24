## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()
options(width=100, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")))

## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE-------------------
suppressPackageStartupMessages({
    library(systemPipeR)
    library(systemPipeRdata)
    library(BiocGenerics)
})

## ----install, eval=FALSE-------------------------------------------------
## if (!requireNamespace("BiocManager", quietly=TRUE))
##     install.packages("BiocManager")
## BiocManager::install("systemPipeRdata") # Installs from Bioconductor once
##                                           # available there
## BiocManager::install("tgirke/systemPipeR", build_vignettes=TRUE,
##                      dependencies=TRUE)  # Installs from github

## ----load_systemPipeRdata, eval=TRUE-------------------------------------
library("systemPipeRdata") # Loads the package

## ----documentation_systemPipeRdata, eval=FALSE---------------------------
## library(help="systemPipeRdata") # Lists package info
## vignette("systemPipeRdata") # Opens vignette

## ----generate_workenvir, eval=FALSE--------------------------------------
## genWorkenvir(workflow="varseq", mydirname=NULL)
## setwd("varseq")

## ----workflow_template_structure, eval=FALSE-----------------------------
## workflow_name/ # *.Rnw/*.Rmd scripts and targets file
## param/ # parameter files for command-line software
## data/ # inputs e.g. FASTQ, reference, annotations
## results/ # analysis result files

## ----load_systemPipeR, eval=TRUE-----------------------------------------
library("systemPipeR") 
# Loads systemPipeR which needs to be installed via BiocManager from Bioconductor

## ----documentation_systemPipeR, eval=FALSE-------------------------------
## vignette("systemPipeR", package = "systemPipeR")

## ----return_samplepaths, eval=TRUE---------------------------------------
pathList()

## ----sessionInfo---------------------------------------------------------
sessionInfo()

