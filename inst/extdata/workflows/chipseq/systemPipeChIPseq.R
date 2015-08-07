## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"---------------------------------
BiocStyle::latex(use.unsrturl=FALSE)

## ----setup, include=FALSE, cache=FALSE-------------------------------------------------------
library(knitr)
# set global chunk options for knitr
opts_chunk$set(comment=NA, warning=FALSE, message=FALSE, fig.path='figure/systemPipeR-')
options(formatR.arrow=TRUE, width=95)
unlink("test.db")

## ----eval=TRUE-------------------------------------------------------------------------------
library(systemPipeR)

## ----eval=FALSE------------------------------------------------------------------------------
#  source("systemPipeChIPseq_Fct.R")

## ----eval=TRUE-------------------------------------------------------------------------------
targets <- read.delim("targets_chip.txt", comment.char = "#")[,1:4]
targets

## ----eval=FALSE------------------------------------------------------------------------------
#  args <- systemArgs(sysma="param/bowtieSE.param", mytargets="targets_chip.txt")
#  fqlist <- seeFastq(fastq=infile1(args), batchsize=100000, klength=8)
#  pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
#  seeFastqPlot(fqlist)
#  dev.off()

## ----eval=FALSE------------------------------------------------------------------------------
#  args <- systemArgs(sysma="param/bowtieSE.param", mytargets="targets_chip.txt")
#  sysargs(args)[1] # Command-line parameters for first FASTQ file

## ----eval=FALSE------------------------------------------------------------------------------
#  moduleload(modules(args))
#  system("bowtie2-build ./data/tair10.fasta ./data/tair10.fasta")
#  resources <- list(walltime="20:00:00", nodes=paste0("1:ppn=", cores(args)), memory="10gb")
#  reg <- clusterRun(args, conffile=".BatchJobs.R", template="torque.tmpl", Njobs=18, runid="01",
#                    resourceList=resources)
#  runCommandline(args)

## ----eval=FALSE------------------------------------------------------------------------------
#  file.exists(outpaths(args))

## ----eval=FALSE------------------------------------------------------------------------------
#  read_statsDF <- alignStats(args=args)
#  write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")

## ----eval=FALSE------------------------------------------------------------------------------
#  read.delim("results/alignStats.xls")

## ----eval=FALSE------------------------------------------------------------------------------
#  symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"),
#              urlbase="http://biocluster.ucr.edu/~tgirke/",
#  	    urlfile="./results/IGVurl.txt")

## ----eval=FALSE------------------------------------------------------------------------------
#  writeTargetsout(x=args, file="targets_bam.txt", overwrite=TRUE)
#  args <- systemArgs(sysma="param/macs2_noinput.param", mytargets="targets_bam.txt")
#  sysargs(args)[1] # Command-line parameters for first FASTQ file
#  runCommandline(args)
#  file.exists(outpaths(args))

## ----eval=FALSE------------------------------------------------------------------------------
#  writeTargetsout(x=args, file="targets_bam.txt", overwrite=TRUE)
#  writeTargetsRef(infile="targets_bam.txt", outfile="targets_bam_ref.txt", silent=FALSE, overwrite=FALSE)
#  args <- systemArgs(sysma="param/macs2.param", mytargets="targets_bam_ref.txt")
#  sysargs(args)[1] # Command-line parameters for first FASTQ file
#  runCommandline(args)
#  file.exists(outpaths(args))

## ----sessionInfo, results='asis'-------------------------------------------------------------
toLatex(sessionInfo())

