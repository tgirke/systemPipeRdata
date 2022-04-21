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
## systemPipeRdata::genWorkenvir(workflow = "rnaseq", mydirname = "rnaseq")
## setwd("rnaseq")


## ----load_targets_file, eval=TRUE-------------------------
targetspath <- system.file("extdata", "targetsPE.txt", package = "systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")
targets[1:4,-c(5,6)]


## ----create_workflow, message=FALSE, eval=FALSE-----------
## library(systemPipeR)
## sal <- SPRproject()
## sal


## ----load_SPR, message=FALSE, eval=FALSE, spr=TRUE--------
## appendStep(sal) <- LineWise(code = {
##                 library(systemPipeR)
##                 }, step_name = "load_SPR")


## ----preprocessing, message=FALSE, eval=FALSE, spr=TRUE----
## appendStep(sal) <- SYSargsList(
##     step_name = "preprocessing",
##     targets = "targetsPE.txt", dir = TRUE,
##     wf_file = "preprocessReads/preprocessReads-pe.cwl",
##     input_file = "preprocessReads/preprocessReads-pe.yml",
##     dir_path = system.file("extdata/cwl", package = "systemPipeR"),
##     inputvars = c(
##         FileName1 = "_FASTQ_PATH1_",
##         FileName2 = "_FASTQ_PATH2_",
##         SampleName = "_SampleName_"
##     ),
##     dependency = c("load_SPR"))


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


## ----trimming, eval=FALSE, spr=TRUE-----------------------
## appendStep(sal) <- SYSargsList(
##     step_name = "trimming",
##     targets = "targetsPE.txt",
##     wf_file = "trimmomatic/trimmomatic-pe.cwl", input_file = "trimmomatic/trimmomatic-pe.yml",
##     dir_path = system.file("extdata/cwl", package = "systemPipeR"),
##     inputvars=c(FileName1="_FASTQ_PATH1_", FileName2="_FASTQ_PATH2_", SampleName="_SampleName_"),
##     dependency = "load_SPR",
##     run_step = "optional")


## ----fastq_report, eval=FALSE, message=FALSE, spr=TRUE----
## appendStep(sal) <- LineWise(code = {
##                 fastq <- getColumn(sal, step = "preprocessing", "targetsWF", column = 1)
##                 fqlist <- seeFastq(fastq = fastq, batchsize = 10000, klength = 8)
##                 pdf("./results/fastqReport.pdf", height = 18, width = 4*length(fqlist))
##                 seeFastqPlot(fqlist)
##                 dev.off()
##                 }, step_name = "fastq_report",
##                 dependency = "preprocessing")


## ----hisat2_index, eval=FALSE, spr=TRUE-------------------
## appendStep(sal) <- SYSargsList(
##     step_name = "hisat2_index",
##     dir = FALSE,
##     targets=NULL,
##     wf_file = "hisat2/hisat2-index.cwl",
##     input_file="hisat2/hisat2-index.yml",
##     dir_path="param/cwl",
##     dependency = "load_SPR"
## )


## ----hisat2_mapping, eval=FALSE, spr=TRUE-----------------
## appendStep(sal) <- SYSargsList(
##     step_name = "hisat2_mapping",
##     dir = TRUE,
##     targets ="preprocessing",
##     wf_file = "workflow-hisat2/workflow_hisat2-pe.cwl",
##     input_file = "workflow-hisat2/workflow_hisat2-pe.yml",
##     dir_path = "param/cwl",
##     inputvars = c(preprocessReads_1 = "_FASTQ_PATH1_", preprocessReads_2 = "_FASTQ_PATH2_",
##                   SampleName = "_SampleName_"),
##     rm_targets_col = c("FileName1", "FileName2"),
##     dependency = c("preprocessing", "hisat2_index")
## )


## ----bowtie2_alignment, eval=FALSE------------------------
## cmdlist(sal, step="hisat2_mapping", targets=1)


## ----align_stats, eval=FALSE, spr=TRUE--------------------
## appendStep(sal) <- LineWise(
##     code = {
##         fqpaths <- getColumn(sal, step = "preprocessing", "targetsWF", column = "FileName1")
##         bampaths <- getColumn(sal, step = "hisat2_mapping", "outfiles", column = "samtools_sort_bam")
##         read_statsDF <- alignStats(args = bampaths, fqpaths = fqpaths, pairEnd = TRUE)
##         write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
##         },
##     step_name = "align_stats",
##     dependency = "hisat2_mapping")


## ----bam_IGV, eval=FALSE, spr=TRUE------------------------
## appendStep(sal) <- LineWise(
##     code = {
##         bampaths <- getColumn(sal, step = "hisat2_mapping", "outfiles",
##                   column = "samtools_sort_bam")
##         symLink2bam(
##             sysargs = bampaths, htmldir = c("~/.html/", "somedir/"),
##             urlbase = "http://cluster.hpcc.ucr.edu/~tgirke/",
##             urlfile = "./results/IGVurl.txt")
##     },
##     step_name = "bam_IGV",
##     dependency = "hisat2_mapping",
##     run_step = "optional"
## )


## ----create_db, eval=FALSE, spr=TRUE----------------------
## appendStep(sal) <- LineWise(
##     code = {
##         library(GenomicFeatures)
##         txdb <- suppressWarnings(makeTxDbFromGFF(file="data/tair10.gff", format="gff", dataSource="TAIR", organism="Arabidopsis thaliana"))
##         saveDb(txdb, file="./data/tair10.sqlite")
##         },
##     step_name = "create_db",
##     dependency = "hisat2_mapping")


## ----read_counting, eval=FALSE, spr=TRUE------------------
## appendStep(sal) <- LineWise(
##     code = {
##         library(GenomicFeatures); library(BiocParallel)
##         txdb <- loadDb("./data/tair10.sqlite")
##         outpaths <- getColumn(sal, step = "hisat2_mapping", "outfiles", column = "samtools_sort_bam")
##         eByg <- exonsBy(txdb, by = c("gene"))
##         bfl <- BamFileList(outpaths, yieldSize = 50000, index = character())
##         multicoreParam <- MulticoreParam(workers = 4); register(multicoreParam); registered()
##         counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode = "Union",
##                                                                  ignore.strand = TRUE,
##                                                                  inter.feature = FALSE,
##                                                                  singleEnd = FALSE,
##                                                                  BPPARAM = multicoreParam))
##         countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
##         rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- names(bfl)
##         rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
##         write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
##         write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
##         ## Creating a SummarizedExperiment object
##         colData <- data.frame(row.names=SampleName(sal, "hisat2_mapping"),
##                               condition=getColumn(sal, "hisat2_mapping", position = "targetsWF", column = "Factor"))
##         colData$condition <- factor(colData$condition)
##         countDF_se <- SummarizedExperiment::SummarizedExperiment(assays = countDFeByg,
##                                                                  colData = colData)
##         ## Add results as SummarizedExperiment to the workflow object
##         SE(sal, "read_counting") <- countDF_se
##         },
##     step_name = "read_counting",
##     dependency = "create_db")


## ----sample_tree, eval=FALSE, spr=TRUE--------------------
## appendStep(sal) <- LineWise(
##     code = {
##         library(DESeq2, quietly=TRUE); library(ape, warn.conflicts=FALSE)
##         ## Extracting SummarizedExperiment object
##         se <- SE(sal, "read_counting")
##         dds <- DESeqDataSet(se, design = ~ condition)
##         d <- cor(assay(rlog(dds)), method="spearman")
##         hc <- hclust(dist(1-d))
##         pdf("results/sample_tree.pdf")
##         plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
##         dev.off()
##         },
##     step_name = "sample_tree",
##     dependency = "read_counting")


## ----run_edger, eval=FALSE, spr=TRUE----------------------
## appendStep(sal) <- LineWise(
##     code = {
##         library(edgeR)
##         countDF <- read.delim("results/countDFeByg.xls", row.names=1, check.names=FALSE)
##         cmp <- readComp(stepsWF(sal)[['hisat2_mapping']], format="matrix", delim="-")
##         edgeDF <- run_edgeR(countDF=countDF, targets=targetsWF(sal)[['hisat2_mapping']], cmp=cmp[[1]], independent=FALSE, mdsplot="")
##         },
##     step_name = "run_edger",
##     dependency = "read_counting")


## ----custom_annot, eval=FALSE, spr=TRUE-------------------
## appendStep(sal) <- LineWise(
##     code = {
##         library("biomaRt")
##         m <- useMart("plants_mart", dataset="athaliana_eg_gene", host="https://plants.ensembl.org")
##         desc <- getBM(attributes=c("tair_locus", "description"), mart=m)
##         desc <- desc[!duplicated(desc[,1]),]
##         descv <- as.character(desc[,2]); names(descv) <- as.character(desc[,1])
##         edgeDF <- data.frame(edgeDF, Desc=descv[rownames(edgeDF)], check.names=FALSE)
##         write.table(edgeDF, "./results/edgeRglm_allcomp.xls", quote=FALSE, sep="\t", col.names = NA)
##         },
##     step_name = "custom_annot",
##     dependency = "run_edger")


## ----filter_degs, eval=FALSE, spr=TRUE--------------------
## appendStep(sal) <- LineWise(
##     code = {
##         edgeDF <- read.delim("results/edgeRglm_allcomp.xls", row.names=1, check.names=FALSE)
##         pdf("results/DEGcounts.pdf")
##         DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=20))
##         dev.off()
##         write.table(DEG_list$Summary, "./results/DEGcounts.xls", quote=FALSE, sep="\t", row.names=FALSE)
##         },
##     step_name = "filter_degs",
##     dependency = "custom_annot")


## ----venn_diagram, eval=FALSE, spr=TRUE-------------------
## appendStep(sal) <- LineWise(
##     code = {
##         vennsetup <- overLapper(DEG_list$Up[6:9], type="vennsets")
##         vennsetdown <- overLapper(DEG_list$Down[6:9], type="vennsets")
##         pdf("results/vennplot.pdf")
##         vennPlot(list(vennsetup, vennsetdown), mymain="", mysub="", colmode=2, ccol=c("blue", "red"))
##         dev.off()
##         },
##     step_name = "venn_diagram",
##     dependency = "filter_degs")


## ----get_go_annot, eval=FALSE, spr=TRUE-------------------
## appendStep(sal) <- LineWise(
##     code = {
##         library("biomaRt")
##         # listMarts() # To choose BioMart database
##         # listMarts(host="plants.ensembl.org")
##         m <- useMart("plants_mart", host="https://plants.ensembl.org")
##         #listDatasets(m)
##         m <- useMart("plants_mart", dataset="athaliana_eg_gene", host="https://plants.ensembl.org")
##         # listAttributes(m) # Choose data types you want to download
##         go <- getBM(attributes=c("go_id", "tair_locus", "namespace_1003"), mart=m)
##         go <- go[go[,3]!="",]; go[,3] <- as.character(go[,3])
##         go[go[,3]=="molecular_function", 3] <- "F"; go[go[,3]=="biological_process", 3] <- "P"; go[go[,3]=="cellular_component", 3] <- "C"
##         go[1:4,]
##         if(!dir.exists("./data/GO")) dir.create("./data/GO")
##         write.table(go, "data/GO/GOannotationsBiomart_mod.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
##         catdb <- makeCATdb(myfile="data/GO/GOannotationsBiomart_mod.txt", lib=NULL, org="", colno=c(1,2,3), idconv=NULL)
##         save(catdb, file="data/GO/catdb.RData")
##         },
##     step_name = "get_go_annot",
##     dependency = "filter_degs")


## ----go_enrich, eval=FALSE, spr=TRUE----------------------
## appendStep(sal) <- LineWise(
##     code = {
##         library("biomaRt")
##         load("data/GO/catdb.RData")
##         DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=50), plot=FALSE)
##         up_down <- DEG_list$UporDown; names(up_down) <- paste(names(up_down), "_up_down", sep="")
##         up <- DEG_list$Up; names(up) <- paste(names(up), "_up", sep="")
##         down <- DEG_list$Down; names(down) <- paste(names(down), "_down", sep="")
##         DEGlist <- c(up_down, up, down)
##         DEGlist <- DEGlist[sapply(DEGlist, length) > 0]
##         BatchResult <- GOCluster_Report(catdb=catdb, setlist=DEGlist, method="all", id_type="gene", CLSZ=2, cutoff=0.9, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
##         m <- useMart("plants_mart", dataset="athaliana_eg_gene", host="https://plants.ensembl.org")
##         goslimvec <- as.character(getBM(attributes=c("goslim_goa_accession"), mart=m)[,1])
##         BatchResultslim <- GOCluster_Report(catdb=catdb, setlist=DEGlist, method="slim", id_type="gene", myslimv=goslimvec, CLSZ=10, cutoff=0.01, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
##         write.table(BatchResultslim, "results/GOBatchSlim.xls", row.names=FALSE, quote=FALSE, sep="\t")
##         },
##     step_name = "go_enrich",
##     dependency = "get_go_annot")


## ----go_plot, eval=FALSE, spr=TRUE------------------------
## appendStep(sal) <- LineWise(
##     code = {
##         gos <- BatchResultslim[grep("M6-V6_up_down", BatchResultslim$CLID), ]
##         gos <- BatchResultslim
##         pdf("results/GOslimbarplotMF.pdf", height=8, width=10)
##         goBarplot(gos, gocat="MF")
##         goBarplot(gos, gocat="BP")
##         goBarplot(gos, gocat="CC")
##         dev.off()
##         },
##     step_name = "go_plot",
##     dependency = "go_enrich")


## ----heatmap, eval=FALSE, spr=TRUE------------------------
## appendStep(sal) <- LineWise(
##     code = {
##         library(pheatmap)
##         geneids <- unique(as.character(unlist(DEG_list[[1]])))
##         y <- assay(rlog(dds))[geneids, ]
##         pdf("results/heatmap1.pdf")
##         pheatmap(y, scale="row", clustering_distance_rows="correlation", clustering_distance_cols="correlation")
##         dev.off()
##         },
##     step_name = "heatmap",
##     dependency = "go_enrich")


## ----sessionInfo, eval=FALSE, spr=TRUE--------------------
## appendStep(sal) <- LineWise(
##     code = {
##         sessionInfo()
##         },
##     step_name = "sessionInfo",
##     dependency = "heatmap")


## ----runWF, eval=FALSE------------------------------------
## sal <- runWF(sal)


## ----runWF_cluster, eval=FALSE----------------------------
## resources <- list(conffile=".batchtools.conf.R",
##                   template="batchtools.slurm.tmpl",
##                   Njobs=18,
##                   walltime=120, ## minutes
##                   ntasks=1,
##                   ncpus=4,
##                   memory=1024, ## Mb
##                   partition = "short"
##                   )
## sal <- addResources(sal, c("hisat2_mapping"), resources = resources)
## sal <- runWF(sal)


## ----plotWF, eval=FALSE-----------------------------------
## plotWF(sal, rstudio = TRUE)


## ----statusWF, eval=FALSE---------------------------------
## sal
## statusWF(sal)


## ----logsWF, eval=FALSE-----------------------------------
## sal <- renderLogs(sal)

