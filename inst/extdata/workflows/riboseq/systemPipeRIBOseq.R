## pre code {

## white-space: pre !important;

## overflow-x: scroll !important;

## word-break: keep-all !important;

## word-wrap: initial !important;

## }


## ----style, echo = FALSE, results = 'asis'----------------
BiocStyle::markdown()
options(width=60, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")), 
    tidy.opts=list(width.cutoff=60), tidy=TRUE)


## ----setup, echo=FALSE, message=FALSE, warning=FALSE, eval=FALSE----
## suppressPackageStartupMessages({
##     library(systemPipeR)
## })


## ----load_systempiper, eval=TRUE, message=FALSE, warning=FALSE----
library(systemPipeR)


## ----genNew_wf, eval=FALSE--------------------------------
## systemPipeRdata::genWorkenvir(workflow = "riboseq", mydirname = "riboseq")
## setwd("riboseq")


## ----load_targets, eval=TRUE------------------------------
targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")[,1:4]
targets


## ----create_workflow, message=FALSE, eval=FALSE-----------
## library(systemPipeR)
## sal <- SPRproject()
## sal


## ----load_SPR, message=FALSE, eval=FALSE, spr=TRUE--------
## appendStep(sal) <- LineWise(
##     code = {
##         library(systemPipeR)
##         library(rtracklayer)
##         library(GenomicFeatures)
##         library(ggplot2)
##         library(grid)
##         library(BiocParallel)
##         library(DESeq2, quietly=TRUE)
##         library(ape, warn.conflicts=FALSE)
##         library(edgeR)
##         library(biomaRt)
##         library(BBmisc) # Defines suppressAll()
##         library(pheatmap)
##         library(BiocParallel)
##     }, step_name = "load_SPR")


## ----preprocessing, message=FALSE, eval=FALSE, spr=TRUE----
## appendStep(sal) <- SYSargsList(
##     step_name = "preprocessing",
##     targets = "targetsPE.txt", dir = TRUE,
##     wf_file = "preprocessReads/preprocessReads-pe.cwl",
##     input_file = "preprocessReads/preprocessReads-pe_riboseq.yml",
##     dir_path = system.file("extdata/cwl", package = "systemPipeR"),
##     inputvars = c(
##         FileName1 = "_FASTQ_PATH1_",
##         FileName2 = "_FASTQ_PATH2_",
##         SampleName = "_SampleName_"
##     ),
##     dependency = c("load_SPR"))


## ----preprocessing_check, message=FALSE, eval=FALSE-------
## yamlinput(sal, step="preprocessing")$Fct
## # [1] "'trimbatch(fq, pattern=\"ACACGTCT\", internalmatch=FALSE, minpatternlength=6, Nnumber=1, polyhomo=50, minreadlength=16, maxreadlength=101)'"
## cmdlist(sal, step = "preprocessing", targets = 1 )


## ----fastq_report, eval=FALSE, message=FALSE, spr=TRUE----
## appendStep(sal) <- LineWise(
##     code = {
##         fq_files <- getColumn(sal, "preprocessing", "targetsWF", column = 1)
##         fqlist <- seeFastq(fastq = fq_files, batchsize = 10000, klength = 8)
##         png("./results/fastqReport.png", height = 162, width = 288 * length(fqlist))
##         seeFastqPlot(fqlist)
##         dev.off()
##     },
##     step_name = "fastq_report",
##     dependency = "preprocessing"
## )


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
##     targets ="targetsPE.txt",
##     wf_file = "workflow-hisat2/workflow_hisat2-pe.cwl",
##     input_file = "workflow-hisat2/workflow_hisat2-pe.yml",
##     dir_path = "param/cwl",
##     inputvars = c(FileName1 = "_FASTQ_PATH1_", FileName2 = "_FASTQ_PATH2_",
##                   SampleName = "_SampleName_"),
##     dependency = c("hisat2_index")
## )


## ----bowtie2_alignment, eval=FALSE------------------------
## cmdlist(sal, step="hisat2_mapping", targets=1)


## ----align_stats, eval=FALSE, spr=TRUE--------------------
## appendStep(sal) <- LineWise(
##     code = {
##         fqpaths <- getColumn(sal, step = "hisat2_mapping", "targetsWF", column = "FileName1")
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


## ----genFeatures, eval=FALSE, spr=TRUE--------------------
## appendStep(sal) <- LineWise(
##     code = {
##         txdb <- suppressWarnings(makeTxDbFromGFF(file="data/tair10.gff", format="gff3", dataSource="TAIR", organism="Arabidopsis thaliana"))
##         feat <- genFeatures(txdb, featuretype="all", reduce_ranges=TRUE, upstream=1000,
##                             downstream=0, verbose=TRUE)
##     },
##     step_name = "genFeatures",
##     dependency = "hisat2_mapping",
##     run_step = "mandatory"
## )


## ----featuretypeCounts, eval=FALSE, spr=TRUE--------------
## appendStep(sal) <- LineWise(
##     code = {
##         outpaths <- getColumn(sal, step = "hisat2_mapping", "outfiles", column = "samtools_sort_bam")
##         fc <- featuretypeCounts(bfl=BamFileList(outpaths, yieldSize=50000), grl=feat,
##                                 singleEnd=FALSE, readlength=NULL, type="data.frame")
##         p <- plotfeaturetypeCounts(x=fc, graphicsfile="results/featureCounts.png",
##                                    graphicsformat="png", scales="fixed", anyreadlength=TRUE,
##                                    scale_length_val=NULL)
##     },
##     step_name = "featuretypeCounts",
##     dependency = "genFeatures",
##     run_step = "mandatory"
## )


## ----featuretypeCounts_length, eval=FALSE, spr=TRUE-------
## appendStep(sal) <- LineWise(
##     code = {
##         fc2 <- featuretypeCounts(bfl=BamFileList(outpaths, yieldSize=50000), grl=feat,
##                                  singleEnd=TRUE, readlength=c(74:76,99:102), type="data.frame")
##         p2 <- plotfeaturetypeCounts(x=fc2, graphicsfile="results/featureCounts2.png",
##                                     graphicsformat="png", scales="fixed", anyreadlength=FALSE,
##                             scale_length_val=NULL)
##     },
##     step_name = "featuretypeCounts_length",
##     dependency = "featuretypeCounts",
##     run_step = "mandatory"
## )


## ----pred_ORF, eval=FALSE, spr=TRUE-----------------------
## appendStep(sal) <- LineWise(
##     code = {
##         txdb <- suppressWarnings(makeTxDbFromGFF(file="data/tair10.gff", format="gff3", organism="Arabidopsis"))
##         futr <- fiveUTRsByTranscript(txdb, use.names=TRUE)
##         dna <- extractTranscriptSeqs(FaFile("data/tair10.fasta"), futr)
##         uorf <- predORF(dna, n="all", mode="orf", longest_disjoint=TRUE, strand="sense")
##     },
##     step_name = "pred_ORF",
##     dependency = "featuretypeCounts_length"
## )


## ----scale_ranges, eval=FALSE, spr=TRUE-------------------
## appendStep(sal) <- LineWise(
##     code = {
##         grl_scaled <- scaleRanges(subject=futr, query=uorf, type="uORF", verbose=TRUE)
##         export.gff3(unlist(grl_scaled), "results/uorf.gff")
##     },
##     step_name = "scale_ranges",
##     dependency = "pred_ORF"
## )


## ----translate, eval=FALSE, spr=TRUE----------------------
## appendStep(sal) <- LineWise(
##     code = {
##         translate(unlist(getSeq(FaFile("data/tair10.fasta"), grl_scaled[[7]])))
##     },
##     step_name = "translate",
##     dependency = "scale_ranges"
## )


## ----add_features, eval=FALSE, spr=TRUE-------------------
## appendStep(sal) <- LineWise(
##     code = {
##         feat <- genFeatures(txdb, featuretype="all", reduce_ranges=FALSE)
##         feat <- c(feat, GRangesList("uORF"=unlist(grl_scaled)))
##     },
##     step_name = "add_features",
##     dependency = c("genFeatures", "scale_ranges")
## )


## ----pred_sORFs, eval=FALSE, spr=TRUE---------------------
## appendStep(sal) <- LineWise(
##     code = {
##         feat <- genFeatures(txdb, featuretype="intergenic", reduce_ranges=TRUE)
##         intergenic <- feat$intergenic
##         strand(intergenic) <- "+"
##         dna <- getSeq(FaFile("data/tair10.fasta"), intergenic)
##         names(dna) <- mcols(intergenic)$feature_by
##         sorf <- suppressWarnings(predORF(dna, n="all", mode="orf", longest_disjoint=TRUE, strand="both"))
##         sorf <- sorf[width(sorf) > 60] # Remove sORFs below length cutoff, here 60bp
##         intergenic <- split(intergenic, mcols(intergenic)$feature_by)
##         grl_scaled_intergenic <- scaleRanges(subject=intergenic, query=sorf, type="sORF", verbose=TRUE)
##         export.gff3(unlist(grl_scaled_intergenic), "sorf.gff")
##         translate(getSeq(FaFile("data/tair10.fasta"), unlist(grl_scaled_intergenic)))
##     },
##     step_name = "pred_sORFs",
##     dependency = c("add_features")
## )


## ----binned_CDS_coverage, eval=FALSE, spr=TRUE------------
## appendStep(sal) <- LineWise(
##     code = {
##         grl <- cdsBy(txdb, "tx", use.names=TRUE)
##         fcov <- featureCoverage(bfl=BamFileList(outpaths[1:2]), grl=grl[1:4],
##                                 resizereads=NULL, readlengthrange=NULL, Nbins=20, method=mean,
##                                 fixedmatrix=FALSE, resizefeatures=TRUE, upstream=20,
##                                 downstream=20, outfile="results/featureCoverage.xls",
##                                 overwrite=TRUE)
##     },
##     step_name = "binned_CDS_coverage",
##     dependency = c("add_features")
## )


## ----coverage_upstream_downstream, eval=FALSE, spr=TRUE----
## appendStep(sal) <- LineWise(
##     code = {
##         fcov <- featureCoverage(bfl=BamFileList(outpaths[1:4]), grl=grl[1:12], resizereads=NULL,
##                                 readlengthrange=NULL, Nbins=NULL, method=mean, fixedmatrix=TRUE,
##                                 resizefeatures=TRUE, upstream=20, downstream=20,
##                                 outfile="results/featureCoverage.xls", overwrite=TRUE)
##         png("./results/coverage_upstream_downstream.png", height=12, width=24, units="in", res=72)
##         plotfeatureCoverage(covMA=fcov, method=mean, scales="fixed", extendylim=2,
##                             scale_count_val=10^6)
##         dev.off()
##     },
##     step_name = "coverage_upstream_downstream",
##     dependency = c("binned_CDS_coverage")
## )


## ----coverage_combined, eval=FALSE, spr=TRUE--------------
## appendStep(sal) <- LineWise(
##     code = {
##         fcov <- featureCoverage(bfl=BamFileList(outpaths[1:4]), grl=grl[1:4],
##                                 resizereads=NULL, readlengthrange=NULL, Nbins=20, method=mean,
##                                 fixedmatrix=TRUE, resizefeatures=TRUE, upstream=20,
##                                 downstream=20,outfile="results/featureCoverage.xls",
##                                 overwrite=TRUE)
##         png("./results/featurePlot.png", height=12, width=24, units="in", res=72)
##         plotfeatureCoverage(covMA=fcov, method=mean, scales="fixed", extendylim=2,
##                             scale_count_val=10^6)
##         dev.off()
##     },
##     step_name = "coverage_combined",
##     dependency = c("binned_CDS_coverage", "coverage_upstream_downstream")
## )


## ----coverage_nuc_level, eval=FALSE, spr=TRUE-------------
## appendStep(sal) <- LineWise(
##     code = {
##         fcov <- featureCoverage(bfl=BamFileList(outpaths[1:2]), grl=grl[1],
##                                 resizereads=NULL, readlengthrange=NULL, Nbins=NULL, method=mean,
##                                 fixedmatrix=FALSE, resizefeatures=TRUE, upstream=20,
##                                 downstream=20, outfile=NULL)
##     },
##     step_name = "coverage_nuc_level",
##     dependency = c("coverage_combined")
## )


## ----read_counting, eval=FALSE, spr=TRUE------------------
## appendStep(sal) <- LineWise(
##     code = {
##         txdb <- loadDb("./data/tair10.sqlite")
##         eByg <- exonsBy(txdb, by=c("gene"))
##         bfl <- BamFileList(outpaths, yieldSize=50000, index=character())
##         multicoreParam <- MulticoreParam(workers = 8); register(multicoreParam); registered()
##         counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union",
##                                                                  ignore.strand=TRUE,
##                                                                  inter.feature=FALSE,
##                                                                  singleEnd=FALSE,
##                                                                  BPPARAM = multicoreParam))
##         countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
##         rownames(countDFeByg) <- names(rowRanges(counteByg[[1]]))
##         colnames(countDFeByg) <- names(bfl)
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
##     },
##     step_name = "read_counting",
##     dependency = c("featuretypeCounts")
## )


## ----read_counting_view, eval=TRUE------------------------
read.delim(system.file("extdata/countDFeByg.xls", package = "systemPipeR"),
           row.names=1, check.names=FALSE)[1:4,1:5]


## ----read_rpkm_view, eval=FALSE---------------------------
## read.delim(system.file("extdata/rpkmDFeByg.xls", package = "systemPipeR"),
##            row.names=1, check.names=FALSE)[1:4,1:5]


## ----sample_tree, eval=FALSE, eval=FALSE, spr=TRUE--------
## appendStep(sal) <- LineWise(
##     code = {
##         ## Extracting SummarizedExperiment object
##         se <- SE(sal, "read_counting")
##         dds <- DESeqDataSet(se, design = ~ condition)
##         d <- cor(assay(rlog(dds)), method="spearman")
##         hc <- hclust(dist(1-d))
##         png("results/sample_tree.png")
##         plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
##         dev.off()
##         },
##     step_name = "sample_tree",
##     dependency = "read_counting")


## ----run_edgeR, eval=FALSE, spr=TRUE----------------------
## appendStep(sal) <- LineWise(
##     code = {
##         countDF <- read.delim("results/countDFeByg.xls", row.names=1, check.names=FALSE)
##         cmp <- readComp(stepsWF(sal)[['hisat2_mapping']], format="matrix", delim="-")
##         edgeDF <- run_edgeR(countDF=countDF, targets=targetsWF(sal)[['hisat2_mapping']], cmp=cmp[[1]], independent=FALSE, mdsplot="")
##         },
##     step_name = "run_edgeR",
##     dependency = "read_counting")


## ----custom_annot, eval=FALSE, spr=TRUE-------------------
## appendStep(sal) <- LineWise(
##     code = {
##         m <- useMart("plants_mart", dataset="athaliana_eg_gene", host="https://plants.ensembl.org")
##         desc <- getBM(attributes=c("tair_locus", "description"), mart=m)
##         desc <- desc[!duplicated(desc[,1]),]
##         descv <- as.character(desc[,2]); names(descv) <- as.character(desc[,1])
##         edgeDF <- data.frame(edgeDF, Desc=descv[rownames(edgeDF)], check.names=FALSE)
##         write.table(edgeDF, "./results/edgeRglm_allcomp.xls", quote=FALSE, sep="\t", col.names = NA)
##         },
##     step_name = "custom_annot",
##     dependency = "run_edgeR")


## ----filter_degs, eval=FALSE, spr=TRUE--------------------
## appendStep(sal) <- LineWise(
##     code = {
##         edgeDF <- read.delim("results/edgeRglm_allcomp.xls", row.names=1, check.names=FALSE)
##         png("results/DEGcounts.png")
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
##         png("results/vennplot.png")
##         vennPlot(list(vennsetup, vennsetdown), mymain="", mysub="", colmode=2, ccol=c("blue", "red"))
##         dev.off()
##         },
##     step_name = "venn_diagram",
##     dependency = "filter_degs")


## ----get_go_annot, eval=FALSE, spr=TRUE-------------------
## appendStep(sal) <- LineWise(
##     code = {
##         # listMarts() # To choose BioMart database
##         # listMarts(host="plants.ensembl.org")
##         # m <- useMart("plants_mart", host="https://plants.ensembl.org")
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
##         },
##     step_name = "go_enrich",
##     dependency = "get_go_annot")


## ----go_plot, eval=FALSE, spr=TRUE------------------------
## appendStep(sal) <- LineWise(
##     code = {
##         gos <- BatchResultslim[grep("M6-V6_up_down", BatchResultslim$CLID), ]
##         gos <- BatchResultslim
##         png("results/GOslimbarplotMF.png", height=8, width=10)
##         goBarplot(gos, gocat="MF")
##         goBarplot(gos, gocat="BP")
##         goBarplot(gos, gocat="CC")
##         dev.off()
##         },
##     step_name = "go_plot",
##     dependency = "go_enrich")


## ----diff_loading, eval=FALSE, spr=TRUE-------------------
## appendStep(sal) <- LineWise(
##     code = {
##         countDFeByg <- read.delim("results/countDFeByg.xls", row.names=1, check.names=FALSE)
##         coldata <- S4Vectors::DataFrame(assay=factor(rep(c("Ribo","mRNA"), each=4)),
##                                         condition= factor(rep(as.character(targetsWF(sal)[['hisat2_mapping']]$Factor[1:4]), 2)),
##                                         row.names=as.character(targetsWF(sal)[['hisat2_mapping']]$SampleName)[1:8])
##         coldata
##         },
##     step_name = "diff_loading",
##     dependency = "go_plot")


## ----diff_translational_eff, eval=FALSE, spr=TRUE---------
## appendStep(sal) <- LineWise(
##     code = {
##         dds <- DESeq2::DESeqDataSetFromMatrix(countData=as.matrix(countDFeByg[,rownames(coldata)]),
##                                               colData = coldata,
##                                               design = ~ assay + condition + assay:condition)
##         # model.matrix(~ assay + condition + assay:condition, coldata) # Corresponding design matrix
##         dds <- DESeq2::DESeq(dds, test="LRT", reduced = ~ assay + condition)
##         res <- DESeq2::results(dds)
##         head(res[order(res$padj),],4)
##         write.table(res, file="transleff.xls", quote=FALSE, col.names = NA, sep="\t")
##         },
##     step_name = "diff_translational_eff",
##     dependency = "diff_loading")


## ----heatmap, eval=FALSE, spr=TRUE------------------------
## appendStep(sal) <- LineWise(
##     code = {
##         geneids <- unique(as.character(unlist(DEG_list[[1]])))
##         y <- assay(rlog(dds))[geneids, ]
##         y <- y[rowSums(y[])>0,]
##         png("results/heatmap1.png")
##         pheatmap(y, scale="row", clustering_distance_rows="correlation", clustering_distance_cols="correlation")
##         dev.off()
##         },
##     step_name = "heatmap",
##     dependency = "diff_translational_eff")


## ----sessionInfo, eval=FALSE, spr=TRUE--------------------
## appendStep(sal) <- LineWise(
##     code = {
##         sessionInfo()
##         },
##     step_name = "sessionInfo",
##     dependency = "heatmap")


## ----runWF, eval=FALSE------------------------------------
## sal <- runWF(sal, run_step = "mandatory")


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
## sal <- runWF(sal, run_step = "mandatory")


## ----plotWF, eval=FALSE-----------------------------------
## plotWF(sal, rstudio = TRUE)


## ----statusWF, eval=FALSE---------------------------------
## sal
## statusWF(sal)


## ----logsWF, eval=FALSE-----------------------------------
## sal <- renderLogs(sal)


## ----sessionInfo_final, eval=TRUE-------------------------
sessionInfo()

