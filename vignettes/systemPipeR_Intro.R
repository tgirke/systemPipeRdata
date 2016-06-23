## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()
options(width=100, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")),
    warning=FALSE, message=FALSE)

## ----download_latest, eval=FALSE-----------------------------------------
## download.file("https://raw.githubusercontent.com/tgirke/systemPipeRdata/master/vignettes/systemPipeR_Intro.Rmd", "systemPipeR_Intro.Rmd")

## ----install, eval=FALSE-------------------------------------------------
## source("http://bioconductor.org/biocLite.R") # Sources the biocLite.R installation script
## biocLite("systemPipeR") # Installs systemPipeR
## biocLite("systemPipeRdata") # Installs systemPipeRdata

## ----documentation, eval=FALSE-------------------------------------------
## library("systemPipeR") # Loads the package
## library(help="systemPipeR") # Lists package info
## vignette("systemPipeR") # Opens vignette

## ----genRna_workflow, eval=FALSE-----------------------------------------
## library(systemPipeRdata)
## genWorkenvir(workflow="riboseq", bam=TRUE)
## setwd("riboseq")

## ----targetsSE, eval=TRUE------------------------------------------------
library(systemPipeR)
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR") 
read.delim(targetspath, comment.char = "#")

## ----targetsPE, eval=TRUE------------------------------------------------
targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
read.delim(targetspath, comment.char = "#")[1:2,1:6]

## ----comment_lines, eval=TRUE--------------------------------------------
readLines(targetspath)[1:4]

## ----targetscomp, eval=TRUE----------------------------------------------
readComp(file=targetspath, format="vector", delim="-")

## ----param_structure, eval=TRUE------------------------------------------
parampath <- system.file("extdata", "tophat.param", package="systemPipeR")

## ----param_import, eval=TRUE---------------------------------------------
args <- suppressWarnings(systemArgs(sysma=parampath, mytargets=targetspath))
args

## ----sysarg_access, eval=TRUE--------------------------------------------
names(args)

## ----sysarg_access2, eval=TRUE-------------------------------------------
sysargs(args)[1]
modules(args)
cores(args)
outpaths(args)[1]

## ----sysarg_json, eval=TRUE----------------------------------------------
systemArgs(sysma=parampath, mytargets=targetspath, type="json")

## ----sessionInfo---------------------------------------------------------
sessionInfo()

