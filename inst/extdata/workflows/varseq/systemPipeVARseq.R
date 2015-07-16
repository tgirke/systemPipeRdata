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

## ----eval=TRUE-------------------------------------------------------------------------------
source("systemPipeVARseq_Fct.R")

## ----eval=TRUE-------------------------------------------------------------------------------
targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")[,1:5]
targets

## ----eval=TRUE, messages=FALSE, warning=FALSE------------------------------------------------
library(BiocParallel); library(BatchJobs)
args <- systemArgs(sysma="param/trimPE.param", mytargets="targetsPE.txt")
f <- function(x) {
    library(systemPipeR)
    source("systemPipeVARseq_Fct.R")
    args <- systemArgs(sysma="param/trimPE.param", mytargets="targetsPE.txt")
    preprocessReads(args=args[x], Fct="filterFct(fq, cutoff=20, Nexceptions=0)", batchsize=100000)
}
funs <- makeClusterFunctionsTorque("torque.tmpl")
param <- BatchJobsParam(length(args), resources=list(walltime="20:00:00", nodes="1:ppn=1", memory="16gb"), cluster.functions=funs)
register(param)
d <- bplapply(seq(along=args), f)
writeTargetsout(x=args, file="targets_PEtrim.txt", overwrite=TRUE)

## ----eval=TRUE-------------------------------------------------------------------------------
library(BiocParallel); library(BatchJobs)
args <- systemArgs(sysma="param/bwa.param", mytargets="targets_PEtrim.txt")
f <- function(x) {
    library(systemPipeR)
    args <- systemArgs(sysma="param/bwa.param", mytargets="targets_PEtrim.txt")
    seeFastq(fastq=infile1(args)[x], batchsize=100000, klength=8)
}
funs <- makeClusterFunctionsTorque("torque.tmpl")
param <- BatchJobsParam(length(args), resources=list(walltime="20:00:00", nodes="1:ppn=1", memory="16gb"), cluster.functions=funs)
register(param)
fqlist <- bplapply(seq(along=args), f)
pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
seeFastqPlot(unlist(fqlist, recursive=FALSE))
dev.off()

## ----eval=TRUE-------------------------------------------------------------------------------
args <- systemArgs(sysma="param/bwa.param", mytargets="targets_PEtrim.txt")
sysargs(args)[1] # Command-line parameters for first FASTQ file

## ----eval=FALSE------------------------------------------------------------------------------
#  bampaths <- runCommandline(args=args)

## ----eval=TRUE-------------------------------------------------------------------------------
moduleload(modules(args))
system("bwa index -a bwtsw ./data/tair10.fasta")
resources <- list(walltime="20:00:00", nodes=paste0("1:ppn=", cores(args)), memory="4gb")
reg <- clusterRun(args, conffile=".BatchJobs.R", template="torque.tmpl", Njobs=18, runid="01", 
                  resourceList=resources)
waitForJobs(reg)
writeTargetsout(x=args, file="targets_bam.txt", overwrite=TRUE)

## ----eval=FALSE------------------------------------------------------------------------------
#  file.exists(outpaths(args))

## ----eval=TRUE-------------------------------------------------------------------------------
library(gmapR); library(BiocParallel); library(BatchJobs)
gmapGenome <- GmapGenome(reference(args), directory="data", name="gmap_tair10chr", create=TRUE)
args <- systemArgs(sysma="param/gsnap.param", mytargets="targets_PEtrim.txt")
f <- function(x) {
    library(gmapR); library(systemPipeR)
    args <- systemArgs(sysma="param/gsnap.param", mytargets="targets_PEtrim.txt")
    gmapGenome <- GmapGenome(reference(args), directory="data", name="gmap_tair10chr", create=FALSE)
    p <- GsnapParam(genome=gmapGenome, unique_only=TRUE, molecule="DNA", max_mismatches=3)
    o <- gsnap(input_a=infile1(args)[x], input_b=infile2(args)[x], params=p, output=outfile1(args)[x])
}
funs <- makeClusterFunctionsTorque("torque.tmpl")
param <- BatchJobsParam(length(args), resources=list(walltime="20:00:00", nodes="1:ppn=1", memory="6gb"), cluster.functions=funs)
register(param)
d <- bplapply(seq(along=args), f)
writeTargetsout(x=args, file="targets_gsnap_bam.txt", overwrite=TRUE)
file.exists(outpaths(args))

## ----eval=TRUE-------------------------------------------------------------------------------
args <- systemArgs(sysma="param/bwa.param", mytargets="targets_PEtrim.txt")
read_statsDF <- alignStats(args=args) 
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")

## ----eval=FALSE------------------------------------------------------------------------------
#  symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"),
#              urlbase="http://biocluster.ucr.edu/~tgirke/",
#              urlfile="./results/IGVurl.txt")

## ----eval=TRUE-------------------------------------------------------------------------------
#system("java -jar CreateSequenceDictionary.jar R=./data/tair10.fasta O=./data/tair10.dict")
system("java -jar /opt/picard/1.81/CreateSequenceDictionary.jar R=./data/tair10.fasta O=./data/tair10.dict")
args <- systemArgs(sysma="param/gatk.param", mytargets="targets_bam.txt")
resources <- list(walltime="20:00:00", nodes=paste0("1:ppn=", 1), memory="10gb")
reg <- clusterRun(args, conffile=".BatchJobs.R", template="torque.tmpl", Njobs=18, runid="01",
                  resourceList=resources)
waitForJobs(reg)
writeTargetsout(x=args, file="targets_gatk.txt", overwrite=TRUE)

## ----eval=TRUE-------------------------------------------------------------------------------
args <- systemArgs(sysma="param/sambcf.param", mytargets="targets_bam.txt")
resources <- list(walltime="20:00:00", nodes=paste0("1:ppn=", 1), memory="10gb")
reg <- clusterRun(args, conffile=".BatchJobs.R", template="torque.tmpl", Njobs=18, runid="01",
                  resourceList=resources)
waitForJobs(reg)
writeTargetsout(x=args, file="targets_sambcf.txt", overwrite=TRUE)

