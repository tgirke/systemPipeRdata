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
#  library(systemPipeRdata)
#  genWorkenvir(workflow="chipseq")
#  setwd("chipseq")

## ----eval=FALSE------------------------------------------------------------------------------
#  source("systemPipeChIPseq_Fct.R")

## ----eval=TRUE-------------------------------------------------------------------------------
targets <- read.delim("targets_chip.txt", comment.char = "#")
targets[1:4,-c(5,6)]

## ----eval=FALSE------------------------------------------------------------------------------
#  args <- systemArgs(sysma="param/bowtieSE.param", mytargets="targets_chip_trim.txt")
#  fqlist <- seeFastq(fastq=infile1(args), batchsize=100000, klength=8)
#  pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
#  seeFastqPlot(fqlist)
#  dev.off()

## ----eval=FALSE------------------------------------------------------------------------------
#  args <- systemArgs(sysma="param/bowtieSE.param", mytargets="targets_chip_trim.txt")
#  sysargs(args)[1] # Command-line parameters for first FASTQ file
#  moduleload(modules(args)) # Skip if a module system is not used
#  system("bowtie2-build ./data/tair10.fasta ./data/tair10.fasta") # Indexes reference genome
#  runCommandline(args)
#  writeTargetsout(x=args, file="targets_bam.txt", overwrite=TRUE)

## ----eval=FALSE------------------------------------------------------------------------------
#  file.exists(outpaths(args))

## ----eval=FALSE------------------------------------------------------------------------------
#  read_statsDF <- alignStats(args=args)
#  write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
#  read.delim("results/alignStats.xls")

## ----eval=FALSE------------------------------------------------------------------------------
#  symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"),
#              urlbase="http://biocluster.ucr.edu/~tgirke/",
#              urlfile="./results/IGVurl.txt")

## ----eval=FALSE------------------------------------------------------------------------------
#  args <- systemArgs(sysma=NULL, mytargets="targets_bam.txt")
#  args_merge <- mergeBamByFactor(args, overwrite=TRUE)
#  writeTargetsout(x=args_merge, file="targets_mergeBamByFactor.txt", overwrite=TRUE)

## ----eval=FALSE------------------------------------------------------------------------------
#  args <- systemArgs(sysma="param/macs2_noinput.param", mytargets="targets_mergeBamByFactor.txt")
#  sysargs(args)[1] # Command-line parameters for first FASTQ file
#  runCommandline(args)
#  file.exists(outpaths(args))
#  writeTargetsout(x=args, file="targets_macs.txt", overwrite=TRUE)

## ----eval=FALSE------------------------------------------------------------------------------
#  writeTargetsRef(infile="targets_mergeBamByFactor.txt", outfile="targets_bam_ref.txt", silent=FALSE, overwrite=TRUE)
#  args <- systemArgs(sysma="param/macs2.param", mytargets="targets_bam_ref.txt")
#  sysargs(args)[1] # Command-line parameters for first FASTQ file
#  runCommandline(args)
#  file.exists(outpaths(args))
#  writeTargetsout(x=args, file="targets_macs.txt", overwrite=TRUE)

## ----eval=FALSE------------------------------------------------------------------------------
#  library(ChIPpeakAnno); library(GenomicFeatures)
#  args <- systemArgs(sysma="param/annotate_peaks.param", mytargets="targets_macs.txt")
#  txdb <- loadDb("./data/tair10.sqlite")
#  ge <- genes(txdb, columns=c("tx_name", "gene_id", "tx_type"))
#  for(i in seq(along=args)) {
#      peaksGR <- as(read.delim(infile1(args)[i], comment="#"), "GRanges")
#      annotatedPeak <- annotatePeakInBatch(peaksGR, AnnotationData=genes(txdb))
#      df <- data.frame(as.data.frame(annotatedPeak), as.data.frame(values(ge[values(annotatedPeak)$feature,])))
#      write.table(df, outpaths(args[i]), quote=FALSE, row.names=FALSE, sep="\t")
#  }
#  writeTargetsout(x=args, file="targets_peakanno.txt", overwrite=TRUE)

## ----eval=FALSE, include=FALSE---------------------------------------------------------------
#  ## Perform previous step with full genome annotation from Biomart
#  # txdb <- makeTxDbFromBiomart(biomart="ENSEMBL_MART_PLANT", dataset="athaliana_eg_gene")
#  # tx <- transcripts(txdb, columns=c("tx_name", "gene_id", "tx_type"))
#  # ge <- genes(txdb, columns=c("tx_name", "gene_id", "tx_type")) # works as well
#  # seqlevels(ge) <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrC", "ChrM")
#  # table(mcols(tx)$tx_type)
#  # tx <- tx[!duplicated(unstrsplit(values(tx)$gene_id, sep=",")) # Keeps only first transcript model for each gene]
#  # annotatedPeak <- annotatePeakInBatch(macsOutput, AnnotationData = tx)

## ----eval=FALSE------------------------------------------------------------------------------
#  library(ChIPseeker)
#  txdb <- loadDb("./data/tair10.sqlite")
#  for(i in seq(along=args)) {
#      peakAnno <- annotatePeak(infile1(args)[i], TxDb=txdb, verbose=FALSE)
#      df <- as.data.frame(peakAnno)
#      write.table(df, outpaths(args[i]), quote=FALSE, row.names=FALSE, sep="\t")
#  }
#  writeTargetsout(x=args, file="targets_peakanno.txt", overwrite=TRUE)

## ----eval=FALSE------------------------------------------------------------------------------
#  peak <- readPeakFile(infile1(args)[1])
#  covplot(peak, weightCol="X.log10.pvalue.")
#  peakHeatmap(outpaths(args)[1], TxDb=txdb, upstream=1000, downstream=1000, color="red")
#  plotAvgProf2(outpaths(args)[1], TxDb=txdb, upstream=1000, downstream=1000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

## ----eval=FALSE------------------------------------------------------------------------------
#  library(GenomicRanges)
#  args <- systemArgs(sysma="param/count_rangesets.param", mytargets="targets_macs.txt")
#  args_bam <- systemArgs(sysma=NULL, mytargets="targets_bam.txt")
#  bfl <- BamFileList(outpaths(args_bam), yieldSize=50000, index=character())
#  countDFnames <- countRangeset(bfl, args, mode="Union", ignore.strand=TRUE)
#  writeTargetsout(x=args, file="targets_countDF.txt", overwrite=TRUE)

## ----eval=FALSE------------------------------------------------------------------------------
#  args_diff <- systemArgs(sysma="param/rundiff.param", mytargets="targets_countDF.txt")
#  cmp <- readComp(file=args_bam, format="matrix")
#  dbrlist <- runDiff(args=args_diff, diffFct=run_edgeR, targets=targetsin(args_bam),
#                      cmp=cmp[[1]], independent=TRUE, dbrfilter=c(Fold=2, FDR=1))
#  writeTargetsout(x=args_diff, file="targets_rundiff.txt", overwrite=TRUE)

## ----eval=FALSE------------------------------------------------------------------------------
#  args <- systemArgs(sysma="param/macs2.param", mytargets="targets_bam_ref.txt")
#  args_anno <- systemArgs(sysma="param/annotate_peaks.param", mytargets="targets_macs.txt")
#  annofiles <- outpaths(args_anno)
#  gene_ids <- sapply(names(annofiles), function(x) unique(as.character(read.delim(annofiles[x])[,"gene_id"])))
#  load("data/GO/catdb.RData")
#  BatchResult <- GOCluster_Report(catdb=catdb, setlist=gene_ids, method="all", id_type="gene", CLSZ=2, cutoff=0.9, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)

## ----eval=FALSE------------------------------------------------------------------------------
#  library(Biostrings); library(seqLogo); library(BCRANK)
#  args <- systemArgs(sysma="param/annotate_peaks.param", mytargets="targets_macs.txt")
#  rangefiles <- infile1(args)
#  for(i in seq(along=rangefiles)) {
#      df <- read.delim(rangefiles[i], comment="#")
#      peaks <- as(df, "GRanges")
#      names(peaks) <- paste0(as.character(seqnames(peaks)), "_", start(peaks), "-", end(peaks))
#      peaks <- peaks[order(values(peaks)$X.log10.pvalue, decreasing=TRUE)]
#      pseq <- getSeq(FaFile("./data/tair10.fasta"), peaks)
#      names(pseq) <- names(peaks)
#      writeXStringSet(pseq, paste0(rangefiles[i], ".fasta"))
#  }

## ----eval=FALSE------------------------------------------------------------------------------
#  set.seed(0)
#  BCRANKout <- bcrank(paste0(rangefiles[1], ".fasta"), restarts=25, use.P1=TRUE, use.P2=TRUE)
#  toptable(BCRANKout)
#  topMotif <- toptable(BCRANKout, 1)
#  weightMatrix <- pwm(topMotif, normalize = FALSE)
#  weightMatrixNormalized <- pwm(topMotif, normalize = TRUE)
#  pdf("results/seqlogo.pdf")
#  seqLogo(weightMatrixNormalized)
#  dev.off()

## ----sessionInfo, results='asis'-------------------------------------------------------------
toLatex(sessionInfo())

