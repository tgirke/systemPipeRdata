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
    library(BiocParallel)
    library(Biostrings)
    library(Rsamtools)
    library(GenomicRanges)
    library(ggplot2)
    library(GenomicAlignments)
    library(ShortRead)
})

## ----install, eval=FALSE-------------------------------------------------
## source("http://bioconductor.org/biocLite.R") # Sources the biocLite.R installation script
## biocLite("systemPipeR") # Installs systemPipeR from Bioconductor
## biocLite("tgirke/systemPipeRdata", build_vignettes=TRUE, dependencies=TRUE) # From github

## ----documentation, eval=FALSE-------------------------------------------
## library("systemPipeR") # Loads the package
## library(help="systemPipeR") # Lists package info
## vignette("systemPipeR") # Opens vignette

## ----targetsSE, eval=TRUE------------------------------------------------
library(systemPipeR)
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
read.delim(targetspath, comment.char = "#")

## ----targetsPE, eval=TRUE------------------------------------------------
targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
read.delim(targetspath, comment.char = "#")[1:2,1:6]

## ----targetscomp, eval=TRUE----------------------------------------------
readComp(file=targetspath, format="vector", delim="-")

## ----param_structure, eval=TRUE------------------------------------------
parampath <- system.file("extdata", "tophat.param", package="systemPipeR")
read.delim(parampath, comment.char = "#")

## ----param_import, eval=TRUE---------------------------------------------
args <- suppressWarnings(systemArgs(sysma=parampath, mytargets=targetspath))
args

## ----sysarg_access, eval=TRUE--------------------------------------------
names(args)
modules(args)
cores(args)
outpaths(args)[1]
sysargs(args)[1]

## ----load_package, eval=FALSE--------------------------------------------
## library(systemPipeR)

## ----construct_sysargs, eval=FALSE---------------------------------------
## args <- systemArgs(sysma="trim.param", mytargets="targets.txt")

## ----preprocessing, eval=FALSE-------------------------------------------
## preprocessReads(args=args, Fct="trimLRPatterns(Rpattern='GCCCGGGTAA', subject=fq)",
##                 batchsize=100000, overwrite=TRUE, compress=TRUE)
## writeTargetsout(x=args, file="targets_trim.txt")

## ----custom_preprocessing, eval=FALSE------------------------------------
## args <- systemArgs(sysma="trimPE.param", mytargets="targetsPE.txt")
## filterFct <- function(fq, cutoff=20, Nexceptions=0) {
##     qcount <- rowSums(as(quality(fq), "matrix") <= cutoff)
##     fq[qcount <= Nexceptions] # Retains reads where Phred scores are >= cutoff with N exceptions
## }
## preprocessReads(args=args, Fct="filterFct(fq, cutoff=20, Nexceptions=0)", batchsize=100000)
## writeTargetsout(x=args, file="targets_PEtrim.txt")

## ----fastq_quality, eval=FALSE-------------------------------------------
## fqlist <- seeFastq(fastq=infile1(args), batchsize=10000, klength=8)
## pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
## seeFastqPlot(fqlist)
## dev.off()

## ----sessionInfo---------------------------------------------------------
sessionInfo()

