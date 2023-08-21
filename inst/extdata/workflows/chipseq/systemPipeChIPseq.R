## pre code {

## white-space: pre !important;

## overflow-x: scroll !important;

## word-break: keep-all !important;

## word-wrap: initial !important;

## }


## ----style, echo = FALSE, results = 'asis'----------------
BiocStyle::markdown()
options(width = 60, max.print = 1000)
knitr::opts_chunk$set(
    eval = as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache = as.logical(Sys.getenv("KNITR_CACHE", "TRUE")), 
    tidy.opts = list(width.cutoff = 60), tidy = TRUE)


## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE----
suppressPackageStartupMessages({
    library(systemPipeR)
})


## ----genNew_wf, eval=FALSE--------------------------------
## systemPipeRdata::genWorkenvir(workflow = "chipseq", mydirname = "chipseq")
## setwd("chipseq")


## ----load_targets_file, eval=TRUE-------------------------
targetspath <- system.file("extdata", "targetsPE_chip.txt", package = "systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")
targets[1:4,-c(5,6)]


## ----create_workflow, message=FALSE, eval=FALSE-----------
## library(systemPipeR)
## sal <- SPRproject()
## sal


## ----load_SPR, message=FALSE, eval=FALSE, spr=TRUE--------
## cat(crayon::blue$bold("To use this workflow, following R packages are expected:\n"))
## cat(c("'ggbio", "ChIPseeker", "GenomicFeatures", "GenomicRanges", "Biostrings",
##        "seqLogo", "BCRANK", "readr'\n"), sep = "', '")
## targetspath <- system.file("extdata", "targetsPE_chip.txt", package = "systemPipeR")
## ###pre-end
## appendStep(sal) <- LineWise(code = {
##                             library(systemPipeR)
##                             }, step_name = "load_SPR")


## ----fastq_report, eval=FALSE, message=FALSE, spr=TRUE----
## appendStep(sal) <- LineWise(
##     code = {
##         targets <- read.delim(targetspath, comment.char = "#")
##         updateColumn(sal, step = "load_SPR", position = "targetsWF") <- targets
##         fq_files <- getColumn(sal, "load_SPR", "targetsWF", column = 1)
##         fqlist <- seeFastq(fastq = fq_files, batchsize = 10000, klength = 8)
##         png("./results/fastqReport.png", height = 162, width = 288 * length(fqlist))
##         seeFastqPlot(fqlist)
##         dev.off()
##     },
##     step_name = "fastq_report",
##     dependency = "load_SPR"
## )


## ----preprocessing, message=FALSE, eval=FALSE, spr=TRUE----
## appendStep(sal) <- SYSargsList(
##     step_name = "preprocessing",
##     targets = targetspath, dir = TRUE,
##     wf_file = "preprocessReads/preprocessReads-pe.cwl",
##     input_file = "preprocessReads/preprocessReads-pe.yml",
##     dir_path = system.file("extdata/cwl", package = "systemPipeR"),
##     inputvars = c(
##         FileName1 = "_FASTQ_PATH1_",
##         FileName2 = "_FASTQ_PATH2_",
##         SampleName = "_SampleName_"
##     ),
##     dependency = c("fastq_report")
## )


## ----custom_preprocessing_function, eval=FALSE------------
## appendStep(sal) <- LineWise(
##     code = {
##         filterFct <- function(fq, cutoff = 20, Nexceptions = 0) {
##             qcount <- rowSums(as(quality(fq), "matrix") <= cutoff, na.rm = TRUE)
##             # Retains reads where Phred scores are >= cutoff with N exceptions
##             fq[qcount <= Nexceptions]
##         }
##         save(list = ls(), file = "param/customFCT.RData")
##     },
##     step_name = "custom_preprocessing_function",
##     dependency = "preprocessing"
## )


## ----editing_preprocessing, message=FALSE, eval=FALSE-----
## yamlinput(sal, "preprocessing")$Fct
## yamlinput(sal, "preprocessing", "Fct") <- "'filterFct(fq, cutoff=20, Nexceptions=0)'"
## yamlinput(sal, "preprocessing")$Fct ## check the new function
## cmdlist(sal, "preprocessing", targets = 1) ## check if the command line was updated with success


## ----bowtie2_index, eval=FALSE, spr=TRUE------------------
## appendStep(sal) <- SYSargsList(
##     step_name = "bowtie2_index",
##     dir = FALSE, targets = NULL,
##     wf_file = "bowtie2/bowtie2-index.cwl",
##     input_file = "bowtie2/bowtie2-index.yml",
##     dir_path = system.file("extdata/cwl", package = "systemPipeR"),
##     inputvars = NULL,
##     dependency = c("preprocessing")
## )


## ----bowtie2_alignment, eval=FALSE, spr=TRUE--------------
## appendStep(sal) <- SYSargsList(
##     step_name = "bowtie2_alignment",
##     dir = TRUE,
##     targets = targetspath,
##     wf_file = "workflow-bowtie2/workflow_bowtie2-pe.cwl",
##     input_file = "workflow-bowtie2/workflow_bowtie2-pe.yml",
##     dir_path = system.file("extdata/cwl", package = "systemPipeR"),
##     inputvars = c(
##         FileName1 = "_FASTQ_PATH1_",
##         FileName2 = "_FASTQ_PATH2_",
##         SampleName = "_SampleName_"
##     ),
##     dependency = c("bowtie2_index")
## )


## ----bowtie2_alignment_check, eval=FALSE------------------
## cmdlist(sal, step="bowtie2_alignment", targets=1)


## ----align_stats, eval=FALSE, spr=TRUE--------------------
## appendStep(sal) <- LineWise(
##     code = {
##         fqpaths <- getColumn(sal, step = "bowtie2_alignment", "targetsWF", column = "FileName1")
##         bampaths <- getColumn(sal, step = "bowtie2_alignment", "outfiles", column = "samtools_sort_bam")
##         read_statsDF <- alignStats(args = bampaths, fqpaths = fqpaths, pairEnd = TRUE)
##         write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
##         },
##     step_name = "align_stats",
##     dependency = "bowtie2_alignment")


## ----bam_IGV, eval=FALSE, spr=TRUE------------------------
## appendStep(sal) <- LineWise(
##     code = {
##         bampaths <- getColumn(sal, step = "bowtie2_alignment", "outfiles",
##                               column = "samtools_sort_bam")
##         symLink2bam(
##             sysargs = bampaths, htmldir = c("~/.html/", "somedir/"),
##             urlbase = "http://cluster.hpcc.ucr.edu/~tgirke/",
##             urlfile = "./results/IGVurl.txt")
##     },
##     step_name = "bam_IGV",
##     dependency = "bowtie2_alignment",
##     run_step = "optional"
## )


## ----rle_object, eval=FALSE-------------------------------
## bampaths <- getColumn(sal, step = "bowtie2_alignment", "outfiles", column = "samtools_sort_bam")
## aligns <- readGAlignments(bampaths[1])
## cov <- coverage(aligns)
## cov


## ----resize_align, eval=FALSE-----------------------------
## trim(resize(as(aligns, "GRanges"), width = 200))


## ----rle_slice, eval=FALSE--------------------------------
## islands <- slice(cov, lower = 15)
## islands[[1]]


## ----plot_coverage, eval=FALSE----------------------------
## library(ggbio)
## myloc <- c("Chr1", 1, 100000)
## ga <- readGAlignments(bampaths[1], use.names=TRUE,
##                       param=ScanBamParam(which=GRanges(myloc[1],
##                         IRanges(as.numeric(myloc[2]), as.numeric(myloc[3])))))
## autoplot(ga, aes(color = strand, fill = strand), facets = strand ~ seqnames, stat = "coverage")


## ----merge_bams, eval=FALSE, spr=TRUE---------------------
## appendStep(sal) <- LineWise(
##     code = {
##         bampaths <- getColumn(sal, step = "bowtie2_alignment", "outfiles", column = "samtools_sort_bam")
##         merge_bams <- mergeBamByFactor(args=bampaths, targetsDF = targetsWF(sal)[["bowtie2_alignment"]], out_dir = file.path("results", "merge_bam") ,overwrite=TRUE)
##         updateColumn(sal, step = "merge_bams", position = "targetsWF") <- merge_bams
##     },
##     step_name = "merge_bams",
##     dependency = "bowtie2_alignment"
## )


## ----call_peaks_macs_noref, eval=FALSE, spr=TRUE----------
## cat("Running preprocessing for call_peaks_macs_noref\n")
## # Previous Linewise step is not run at workflow building time,
## # but we need the output as input for this sysArgs step. So
## # we use some preprocess code to predict the output paths
## # to update the output targets of merge_bams, and then
## # them into this next step during workflow building phase.
## mergebam_out_dir = file.path("results", "merge_bam") # make sure this is the same output directory used in merge_bams
## targets_merge_bam <- targetsWF(sal)$bowtie2_alignment
## targets_merge_bam <- targets_merge_bam[, -which(colnames(targets_merge_bam) %in% c("FileName1", "FileName2", "FileName"))]
## targets_merge_bam <- targets_merge_bam[!duplicated(targets_merge_bam$Factor), ]
## targets_merge_bam <- cbind(FileName = file.path(mergebam_out_dir, paste0(targets_merge_bam$Factor, "_merged.bam")), targets_merge_bam)
## updateColumn(sal, step = "merge_bams", position = "targetsWF") <- targets_merge_bam
## # write it out as backup, so you do not need to use preprocess code above again
## writeTargets(sal, step = "merge_bams", file = "targets_merge_bams.txt", overwrite = TRUE)
## 
## ###pre-end
## appendStep(sal) <- SYSargsList(
##     step_name = "call_peaks_macs_noref",
##     targets = "targets_merge_bams.txt",
##     wf_file = "MACS2/macs2-noinput.cwl",
##     input_file = "MACS2/macs2-noinput.yml",
##     dir_path = system.file("extdata/cwl", package = "systemPipeR"),
##     inputvars = c(
##         FileName = "_FASTQ_PATH1_",
##         SampleName = "_SampleName_"
##     ),
##     dependency = c("merge_bams")
## )


## ----call_peaks_macs_withref, eval=FALSE, spr=TRUE--------
## cat("Running preprocessing for call_peaks_macs_withref\n")
## # To generate the reference targets file for the next step, use `writeTargetsRef`,
## # this file needs to be present at workflow building time
## # Use following preprocess code to do so:
## writeTargetsRef(infile = "targets_merge_bams.txt", outfile = "targets_bam_ref.txt", silent = FALSE, overwrite = TRUE)
## 
## ###pre-end
## appendStep(sal) <- SYSargsList(
##     step_name = "call_peaks_macs_withref",
##     targets = "targets_bam_ref.txt",
##     wf_file = "MACS2/macs2-input.cwl",
##     input_file = "MACS2/macs2-input.yml",
##     dir_path = system.file("extdata/cwl", package = "systemPipeR"),
##     inputvars = c(
##         FileName1 = "_FASTQ_PATH1_",
##         FileName2 = "_FASTQ_PATH2_",
##         SampleName = "_SampleName_"
##     ),
##     dependency = c("merge_bams")
## )


## ----consensus_peaks, eval=FALSE, spr=TRUE----------------
## appendStep(sal) <- LineWise(
##     code = {
##         peaks_files <- getColumn(sal, step = "call_peaks_macs_noref", "outfiles", column = "peaks_xls")
##         peak_M1A <- peaks_files["M1A"]
##         peak_M1A <- as(read.delim(peak_M1A, comment="#")[,1:3], "GRanges")
##         peak_A1A <- peaks_files["A1A"]
##         peak_A1A <- as(read.delim(peak_A1A, comment="#")[,1:3], "GRanges")
##         (myol1 <- subsetByOverlaps(peak_M1A, peak_A1A, minoverlap=1))
##             # Returns any overlap
##         myol2 <- olRanges(query=peak_M1A, subject=peak_A1A, output="gr")
##             # Returns any overlap with OL length information
##         myol2[values(myol2)["OLpercQ"][,1]>=50]
##             # Returns only query peaks with a minimum overlap of 50%
##     },
##     step_name = "consensus_peaks",
##     dependency = "call_peaks_macs_noref"
## )


## ----annotation_ChIPseeker, eval=FALSE, spr=TRUE----------
## appendStep(sal) <- LineWise(
##     code = {
##         library(ChIPseeker); library(GenomicFeatures)
##         peaks_files <- getColumn(sal, step = "call_peaks_macs_noref", "outfiles", column = "peaks_xls")
##         txdb <- suppressWarnings(makeTxDbFromGFF(file="data/tair10.gff", format="gff", dataSource="TAIR",
##                                 organism="Arabidopsis thaliana"))
##         for(i in seq(along=peaks_files)) {
##             peakAnno <- annotatePeak(peaks_files[i], TxDb=txdb, verbose=FALSE)
##             df <- as.data.frame(peakAnno)
##             outpaths <- paste0("./results/", names(peaks_files), "_ChIPseeker_annotated.xls")
##             names(outpaths) <- names(peaks_files)
##             write.table(df, outpaths[i], quote=FALSE, row.names=FALSE, sep="\t")
##         }
##         updateColumn(sal, step = "annotation_ChIPseeker", position = "outfiles") <- data.frame(outpaths)
##     },
##     step_name = "annotation_ChIPseeker",
##     dependency = "call_peaks_macs_noref"
## )


## ----ChIPseeker_plots, eval=FALSE, spr=TRUE---------------
## appendStep(sal) <- LineWise(
##     code = {
##         peaks_files <- getColumn(sal, step = "call_peaks_macs_noref", "outfiles", column = "peaks_xls")
##         peak <- readPeakFile(peaks_files[1])
##         png("results/peakscoverage.png")
##         covplot(peak, weightCol="X.log10.pvalue.")
##         dev.off()
##         png("results/peaksHeatmap.png")
##         peakHeatmap(peaks_files[1], TxDb=txdb, upstream=1000, downstream=1000,
##                     color="red")
##         dev.off()
##         png("results/peaksProfile.png")
##         plotAvgProf2(peaks_files[1], TxDb=txdb, upstream=1000, downstream=1000,
##                      xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency",
##                      conf=0.05)
##         dev.off()
##     },
##     step_name = "ChIPseeker_plots",
##     dependency = "annotation_ChIPseeker"
## )


## ----annotation_ChIPpeakAnno, eval=FALSE, spr=TRUE--------
## appendStep(sal) <- LineWise(
##     code = {
##         library(ChIPpeakAnno); library(GenomicFeatures)
##         peaks_files <- getColumn(sal, step = "call_peaks_macs_noref", "outfiles", column = "peaks_xls")
##         txdb <- suppressWarnings(makeTxDbFromGFF(file="data/tair10.gff", format="gff", dataSource="TAIR",
##                                 organism="Arabidopsis thaliana"))
##         ge <- genes(txdb, columns=c("tx_name", "gene_id", "tx_type"))
##         for(i in seq(along=peaks_files)) {
##             peaksGR <- as(read.delim(peaks_files[i], comment="#"), "GRanges")
##             annotatedPeak <- annotatePeakInBatch(peaksGR, AnnotationData=genes(txdb))
##             df <- data.frame(as.data.frame(annotatedPeak),
##                              as.data.frame(values(ge[values(annotatedPeak)$feature,])))
##             df$tx_name <- as.character(lapply(df$tx_name, function(x) paste(unlist(x), sep='', collapse=', ')))
##             df$tx_type <- as.character(lapply(df$tx_type, function(x) paste(unlist(x), sep='', collapse=', ')))
##             outpaths <- paste0("./results/", names(peaks_files), "_ChIPpeakAnno_annotated.xls")
##             names(outpaths) <- names(peaks_files)
##             write.table(df, outpaths[i], quote=FALSE, row.names=FALSE, sep="\t")
##             }
##     },
##     step_name = "annotation_ChIPpeakAnno",
##     dependency = "call_peaks_macs_noref",
##     run_step = "optional"
## )


## ----count_peak_ranges, eval=FALSE, spr=TRUE--------------
## appendStep(sal) <- LineWise(
##     code = {
##         library(GenomicRanges)
##         bam_files <- getColumn(sal, step = "bowtie2_alignment", "outfiles", column = "samtools_sort_bam")
##         args <- getColumn(sal, step = "call_peaks_macs_noref", "outfiles", column = "peaks_xls")
##         outfiles <- paste0("./results/", names(args), "_countDF.xls")
##         bfl <- BamFileList(bam_files, yieldSize=50000, index=character())
##         countDFnames <- countRangeset(bfl, args, outfiles, mode="Union", ignore.strand=TRUE)
##         updateColumn(sal, step = "count_peak_ranges", position = "outfiles") <- data.frame(countDFnames)
##     },
##     step_name = "count_peak_ranges",
##     dependency = "call_peaks_macs_noref",
## )


## ----diff_bind_analysis, eval=FALSE, spr=TRUE-------------
## appendStep(sal) <- LineWise(
##     code = {
##         countDF_files <- getColumn(sal, step = "count_peak_ranges", "outfiles")
##         outfiles <- paste0("./results/", names(countDF_files), "_peaks_edgeR.xls")
##         names(outfiles) <- names(countDF_files)
##         cmp <- readComp(file =stepsWF(sal)[["bowtie2_alignment"]],
##                         format="matrix")
##         dbrlist <- runDiff(args=countDF_files, outfiles = outfiles, diffFct=run_edgeR,
##                            targets=targetsWF(sal)[["bowtie2_alignment"]], cmp=cmp[[1]],
##                            independent=TRUE, dbrfilter=c(Fold=2, FDR=1))
##     },
##     step_name = "diff_bind_analysis",
##     dependency = "count_peak_ranges",
## )


## ----go_enrich, eval=FALSE, spr=TRUE----------------------
## appendStep(sal) <- LineWise(
##     code = {
##         annofiles <- getColumn(sal, step = "annotation_ChIPseeker", "outfiles")
##         gene_ids <- sapply(annofiles,
##                            function(x) unique(as.character
##                                               (read.delim(x)[,"geneId"])), simplify=FALSE)
##         load("data/GO/catdb.RData")
##         BatchResult <- GOCluster_Report(catdb=catdb, setlist=gene_ids, method="all",
##                                         id_type="gene", CLSZ=2, cutoff=0.9,
##                                         gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
##         write.table(BatchResult, "results/GOBatchAll.xls", quote=FALSE, row.names=FALSE, sep="\t")
##     },
##     step_name = "go_enrich",
##     dependency = "annotation_ChIPseeker",
## )


## ----parse_peak_sequences, eval=FALSE, spr=TRUE-----------
## appendStep(sal) <- LineWise(
##     code = {
##         library(Biostrings); library(seqLogo); library(BCRANK)
##         rangefiles <- getColumn(sal, step = "call_peaks_macs_noref", "outfiles")
##         for(i in seq(along=rangefiles)) {
##             df <- read.delim(rangefiles[i], comment="#")
##             peaks <- as(df, "GRanges")
##             names(peaks) <- paste0(as.character(seqnames(peaks)), "_", start(peaks),
##                                    "-", end(peaks))
##             peaks <- peaks[order(values(peaks)$X.log10.pvalue., decreasing=TRUE)]
##             pseq <- getSeq(FaFile("./data/tair10.fasta"), peaks)
##             names(pseq) <- names(peaks)
##             writeXStringSet(pseq, paste0(rangefiles[i], ".fasta"))
##             }
##         },
##     step_name = "parse_peak_sequences",
##     dependency = "call_peaks_macs_noref",
## )


## ----bcrank_enrich, eval=FALSE, spr=TRUE------------------
## appendStep(sal) <- LineWise(
##     code = {
##         library(Biostrings); library(seqLogo); library(BCRANK)
##         rangefiles <- getColumn(sal, step = "call_peaks_macs_noref", "outfiles")
##         set.seed(0)
##         BCRANKout <- bcrank(paste0(rangefiles[1], ".fasta"), restarts=25,
##                             use.P1=TRUE, use.P2=TRUE)
##         toptable(BCRANKout)
##         topMotif <- toptable(BCRANKout, 1)
##         weightMatrix <- pwm(topMotif, normalize = FALSE)
##         weightMatrixNormalized <- pwm(topMotif, normalize = TRUE)
##         png("results/seqlogo.png")
##         seqLogo(weightMatrixNormalized)
##         dev.off()
##         },
##     step_name = "bcrank_enrich",
##     dependency = "call_peaks_macs_noref",
## )


## ----sessionInfo, eval=FALSE, spr=TRUE--------------------
## appendStep(sal) <- LineWise(
##     code = {
##         sessionInfo()
##         },
##     step_name = "sessionInfo",
##     dependency = "bcrank_enrich"
## )


## ----runWF, eval=FALSE------------------------------------
## sal <- runWF(sal)


## ----runWF_cluster, eval=FALSE----------------------------
## # wall time in mins, memory in MB
## resources <- list(conffile=".batchtools.conf.R",
##                   template="batchtools.slurm.tmpl",
##                   Njobs=18,
##                   walltime=120,
##                   ntasks=1,
##                   ncpus=4,
##                   memory=1024,
##                   partition = "short"
##                   )
## sal <- addResources(sal, c("bowtie2_alignment"), resources = resources)
## sal <- runWF(sal)


## ----plotWF, eval=FALSE-----------------------------------
## plotWF(sal, rstudio = TRUE)


## ----statusWF, eval=FALSE---------------------------------
## sal
## statusWF(sal)


## ----logsWF, eval=FALSE-----------------------------------
## sal <- renderLogs(sal)


## ----list_tools-------------------------------------------
if(file.exists(file.path(".SPRproject", "SYSargsList.yml"))) {
    local({
        sal <- systemPipeR::SPRproject(resume = TRUE)
        systemPipeR::listCmdTools(sal)
        systemPipeR::listCmdModules(sal)
    })
} else {
    cat(crayon::blue$bold("Tools and modules required by this workflow are:\n"))
    cat(c("BLAST 2.14.0+"), sep = "\n")
}


## ----report_session_info, eval=TRUE-----------------------
sessionInfo()

