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
    tidy.opts = list(width.cutoff = 60), tidy = TRUE
)


## ----setup, echo=FALSE, message=FALSE, warning=FALSE, eval=TRUE----
suppressPackageStartupMessages({
    library(systemPipeR)
})


## ----genNew_wf, eval=FALSE--------------------------------
## systemPipeRdata::genWorkenvir(workflow = "varseq", mydirname = "varseq")
## setwd("varseq")


## ----create_sal, message=FALSE, eval=FALSE----------------
## sal <- SPRproject()


## ----load_SPR, message=FALSE, eval=FALSE, spr=TRUE--------
## appendStep(sal) <- LineWise(
##     code = {
##         library(systemPipeR)
##         },
##     step_name = "load_SPR")


## ----fastq_report_pre, eval=FALSE, message=FALSE, spr=TRUE----
## appendStep(sal) <- LineWise(
##     code = {
##         targets <- read.delim("targetsPE_varseq.txt", comment.char = "#")
##         updateColumn(sal, step = "load_SPR", position = "targetsWF") <- targets
##         fq_files <- getColumn(sal, "load_SPR", "targetsWF", column = 1)
##         fqlist <- seeFastq(fastq = fq_files, batchsize = 10000, klength = 8)
##         pdf("./results/fastqReport_pre.pdf", height = 18, width = 4 * length(fqlist))
##         seeFastqPlot(fqlist)
##         dev.off()
##     },
##     step_name = "fastq_report_pre",
##     dependency = "load_SPR"
## )


## ----trimmomatic, eval=FALSE, spr=TRUE--------------------
## appendStep(sal) <- SYSargsList(
##     step_name = "trimmomatic",
##     targets = "targetsPE_varseq.txt",
##     wf_file = "trimmomatic/trimmomatic-pe.cwl",
##     input_file = "trimmomatic/trimmomatic-pe.yml",
##     dir_path = "param/cwl",
##     inputvars = c(
##         FileName1 = "_FASTQ_PATH1_",
##         FileName2 = "_FASTQ_PATH2_",
##         SampleName = "_SampleName_"
##     ),
##     dependency = c("fastq_report_pre"),
##     run_step = "optional"
## )


## ----preprocessing, message=FALSE, eval=FALSE, spr=TRUE----
## appendStep(sal) <- SYSargsList(
##     step_name = "preprocessing",
##     targets = "targetsPE_varseq.txt", dir = TRUE,
##     wf_file = "preprocessReads/preprocessReads-pe.cwl",
##     input_file = "preprocessReads/preprocessReads-pe.yml",
##     dir_path = "param/cwl",
##     inputvars = c(
##         FileName1 = "_FASTQ_PATH1_",
##         FileName2 = "_FASTQ_PATH2_",
##         SampleName = "_SampleName_"
##     ),
##     dependency = c("fastq_report_pre"),
##     run_step = "optional"
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


## ----fastq_report_pos, eval=FALSE, message=FALSE, spr=TRUE----
## appendStep(sal) <- LineWise(
##     code = {
##         fq_files <- getColumn(sal, "preprocessing", "outfiles", column = 1) ## get outfiles path
##         fqlist <- seeFastq(fastq = fq_files, batchsize = 10000, klength = 8)
##         pdf("./results/fastqReport_pos.pdf", height = 18, width = 4 * length(fqlist))
##         seeFastqPlot(fqlist)
##         dev.off()
##     },
##     step_name = "fastq_report_pos",
##     dependency = "trimmomatic",
##     run_step = "optional"
## )


## ----bwa_index, eval=FALSE, spr=TRUE----------------------
## appendStep(sal) <- SYSargsList(
##     step_name = "bwa_index",
##     dir = FALSE, targets = NULL,
##     wf_file = "gatk/workflow_bwa-index.cwl",
##     input_file = "gatk/gatk.yaml",
##     dir_path = "param/cwl",
##     dependency = "load_SPR"
## )


## ----fasta_index, eval=FALSE, spr=TRUE--------------------
## appendStep(sal) <- SYSargsList(
##     step_name = "fasta_index",
##     dir = FALSE, targets = NULL,
##     wf_file = "gatk/workflow_fasta_dict.cwl",
##     input_file = "gatk/gatk.yaml",
##     dir_path = "param/cwl",
##     dependency = "bwa_index"
## )


## ----faidx_index, eval=FALSE, spr=TRUE--------------------
## appendStep(sal) <- SYSargsList(
##     step_name = "faidx_index",
##     dir = FALSE, targets = NULL,
##     wf_file = "gatk/workflow_fasta_faidx.cwl",
##     input_file = "gatk/gatk.yaml",
##     dir_path = "param/cwl",
##     dependency = "fasta_index"
## )


## ----bwa_alignment, eval=FALSE, spr=TRUE------------------
## appendStep(sal) <- SYSargsList(
##     step_name = "bwa_alignment",
##     targets = "targetsPE_varseq.txt",
##     wf_file = "gatk/workflow_bwa-pe.cwl",
##     input_file = "gatk/gatk.yaml",
##     dir_path = "param/cwl",
##     inputvars = c(
##         FileName1 = "_FASTQ_PATH1_",
##         FileName2 = "_FASTQ_PATH2_",
##         SampleName = "_SampleName_"
##     ),
##     dependency = c("faidx_index")
## )


## ----align_stats, eval=FALSE, spr=TRUE--------------------
## appendStep(sal) <- LineWise(
##     code = {
##         bampaths <- getColumn(sal, step = "bwa_alignment", "outfiles", column = "samtools_sort_bam")
##         fqpaths <- getColumn(sal, step = "bwa_alignment", "targetsWF", column = "FileName1")
##         read_statsDF <- alignStats(args = bampaths, fqpaths = fqpaths, pairEnd = TRUE)
##         write.table(read_statsDF, "results/alignStats.xls", row.names = FALSE, quote = FALSE, sep = "\t")
##     },
##     step_name = "align_stats",
##     dependency = "bwa_alignment",
##     run_step = "optional"
## )


## ----bam_urls, eval=FALSE, spr=TRUE-----------------------
## appendStep(sal) <- LineWise(
##     code = {
##         bampaths <- getColumn(sal, step = "bwa_alignment", "outfiles", column = "samtools_sort_bam")
##         symLink2bam(
##             sysargs = bampaths, htmldir = c("~/.html/", "somedir/"),
##             urlbase = "http://cluster.hpcc.ucr.edu/~tgirke/",
##             urlfile = "./results/IGVurl.txt"
##         )
##     },
##     step_name = "bam_urls",
##     dependency = "bwa_alignment",
##     run_step = "optional"
## )


## ----fastq2ubam, eval=FALSE, spr=TRUE---------------------
## appendStep(sal) <- SYSargsList(
##     step_name = "fastq2ubam",
##     targets = "targetsPE_varseq.txt",
##     wf_file = "gatk/workflow_gatk_fastq2ubam.cwl",
##     input_file = "gatk/gatk.yaml",
##     dir_path = "param/cwl",
##     inputvars = c(
##         FileName1 = "_FASTQ_PATH1_",
##         FileName2 = "_FASTQ_PATH2_",
##         SampleName = "_SampleName_"
##     ),
##     dependency = c("faidx_index")
## )


## ----merge_bam, eval=FALSE, spr=TRUE----------------------
## appendStep(sal) <- SYSargsList(
##     step_name = "merge_bam",
##     targets = c("bwa_alignment", "fastq2ubam"),
##     wf_file = "gatk/workflow_gatk_mergebams.cwl",
##     input_file = "gatk/gatk.yaml",
##     dir_path = "param/cwl",
##     inputvars = c(
##         bwa_men_sam = "_bwasam_",
##         ubam = "_ubam_",
##         SampleName = "_SampleName_"
##     ),
##     rm_targets_col = c("preprocessReads_1", "preprocessReads_2"),
##     dependency = c("bwa_alignment", "fastq2ubam")
## )


## ----sort, eval=FALSE, spr=TRUE---------------------------
## appendStep(sal) <- SYSargsList(
##     step_name = "sort",
##     targets = "merge_bam",
##     wf_file = "gatk/workflow_gatk_sort.cwl",
##     input_file = "gatk/gatk.yaml",
##     dir_path = "param/cwl",
##     inputvars = c(merge_bam = "_mergebam_", SampleName = "_SampleName_"),
##     rm_targets_col = c(
##         "bwa_men_sam", "ubam", "SampleName_fastq2ubam",
##         "Factor_fastq2ubam", "SampleLong_fastq2ubam",
##         "Experiment_fastq2ubam", "Date_fastq2ubam"
##     ),
##     dependency = c("merge_bam")
## )


## ----mark_dup, eval=FALSE, spr=TRUE-----------------------
## appendStep(sal) <- SYSargsList(
##     step_name = "mark_dup",
##     targets = "sort",
##     wf_file = "gatk/workflow_gatk_markduplicates.cwl",
##     input_file = "gatk/gatk.yaml",
##     dir_path = "param/cwl",
##     inputvars = c(sort_bam = "_sort_", SampleName = "_SampleName_"),
##     rm_targets_col = c("merge_bam"),
##     dependency = c("sort")
## )


## ----fix_tag, eval=FALSE, spr=TRUE------------------------
## appendStep(sal) <- SYSargsList(
##     step_name = "fix_tag",
##     targets = "mark_dup",
##     wf_file = "gatk/workflow_gatk_fixtag.cwl",
##     input_file = "gatk/gatk.yaml",
##     dir_path = "param/cwl",
##     inputvars = c(mark_bam = "_mark_", SampleName = "_SampleName_"),
##     rm_targets_col = c("sort_bam"),
##     dependency = c("mark_dup")
## )


## ----hap_caller, eval=FALSE, spr=TRUE---------------------
## appendStep(sal) <- SYSargsList(
##     step_name = "hap_caller",
##     targets = "fix_tag",
##     wf_file = "gatk/workflow_gatk_haplotypecaller.cwl",
##     input_file = "gatk/gatk.yaml",
##     dir_path = "param/cwl",
##     inputvars = c(fixtag_bam = "_fixed_", SampleName = "_SampleName_"),
##     rm_targets_col = c("mark_bam"),
##     dependency = c("fix_tag")
## )


## ----import, eval=FALSE, spr=TRUE-------------------------
## appendStep(sal) <- SYSargsList(
##     step_name = "import",
##     targets = NULL, dir = FALSE,
##     wf_file = "gatk/workflow_gatk_genomicsDBImport.cwl",
##     input_file = "gatk/gatk.yaml",
##     dir_path = "param/cwl",
##     dependency = c("hap_caller")
## )


## ----call_variants, eval=FALSE, spr=TRUE------------------
## appendStep(sal) <- SYSargsList(
##     step_name = "call_variants",
##     targets = NULL, dir = FALSE,
##     wf_file = "gatk/workflow_gatk_genotypeGVCFs.cwl",
##     input_file = "gatk/gatk.yaml",
##     dir_path = "param/cwl",
##     dependency = c("import")
## )


## ----filter, eval=FALSE, spr=TRUE-------------------------
## appendStep(sal) <- SYSargsList(
##     step_name = "filter",
##     targets = NULL, dir = FALSE,
##     wf_file = "gatk/workflow_gatk_variantFiltration.cwl",
##     input_file = "gatk/gatk.yaml",
##     dir_path = "param/cwl",
##     dependency = c("call_variants")
## )


## ----create_vcf, eval=FALSE, spr=TRUE---------------------
## appendStep(sal) <- SYSargsList(
##     step_name = "create_vcf",
##     targets = "hap_caller",
##     wf_file = "gatk/workflow_gatk_select_variant.cwl",
##     input_file = "gatk/gatk.yaml",
##     dir_path = "param/cwl",
##     inputvars = c(SampleName = "_SampleName_"),
##     dependency = c("hap_caller", "filter")
## )


## ----create_vcf_BCFtool, eval=FALSE, spr=TRUE-------------
## appendStep(sal) <- SYSargsList(
##     step_name = "create_vcf_BCFtool",
##     targets = "bwa_alignment", dir = TRUE,
##     wf_file = "workflow-bcftools/workflow_bcftools.cwl",
##     input_file = "workflow-bcftools/bcftools.yml",
##     dir_path = "param/cwl",
##     inputvars = c(bwa_men_sam = "_bwasam_", SampleName = "_SampleName_"),
##     rm_targets_col = c("preprocessReads_1", "preprocessReads_2"),
##     dependency = "bwa_alignment",
##     run_step = "optional"
## )


## ----inspect_vcf, eval=FALSE------------------------------
## library(VariantAnnotation)
## vcf_raw <- getColumn(sal, "create_vcf")
## vcf <- readVcf(vcf_raw[1], "A. thaliana")
## vcf
## vr <- as(vcf, "VRanges")
## vr


## ----filter_vcf, eval=FALSE, spr=TRUE---------------------
## appendStep(sal) <- LineWise(
##     code = {
##         vcf_raw <- getColumn(sal, "create_vcf")
##         library(VariantAnnotation)
##         filter <- "totalDepth(vr) >= 2 & (altDepth(vr) / totalDepth(vr) >= 0.8)"
##         vcf_filter <- suppressWarnings(filterVars(vcf_raw, filter, organism = "A. thaliana", out_dir = "results/vcf_filter"))
##         # updateColumn(sal, 'create_vcf', data.frame(vcf_filter=vcf_filter))
##     },
##     step_name = "filter_vcf",
##     dependency = "create_vcf",
##     run_step = "optional"
## )


## ----filter_vcf_BCFtools, eval=FALSE, spr=TRUE------------
## appendStep(sal) <- LineWise(
##     code = {
##         vcf_raw <- getColumn(sal, "create_vcf_BCFtool")
##         library(VariantAnnotation)
##         filter <- "totalDepth(vr) >= 2 & (altDepth(vr) / totalDepth(vr) >= 0.8)"
##         vcf_filter <- suppressWarnings(filterVars(vcf_raw, filter, organism = "A. thaliana", out_dir = "results/vcf_filter_BCFtools"))
##         # updateColumn(sal, 'create_vcf', data.frame(vcf_filter=vcf_filter))
##     },
##     step_name = "filter_vcf_BCFtools",
##     dependency = "create_vcf_BCFtool",
##     run_step = "optional"
## )


## ----check_filter, eval=FALSE-----------------------------
## copyEnvir(sal, "vcf_raw", globalenv())
## copyEnvir(sal, "vcf_filter", globalenv())
## length(as(readVcf(vcf_raw[1], genome = "Ath"), "VRanges")[, 1])
## length(as(readVcf(vcf_filter[1], genome = "Ath"), "VRanges")[, 1])


## ----annotate_basics, eval=FALSE--------------------------
## library("GenomicFeatures")
## # comment the next line if optional step "filter_vcf" is included
## vcf_filter <- getColumn(sal, "create_vcf")
## # uncomment the next line if optional step "filter_vcf" is included
## # copyEnvir(sal, "vcf_filter", globalenv())
## txdb <- loadDb("./data/tair10.sqlite")
## vcf <- readVcf(vcf_filter[1], "A. thaliana")
## locateVariants(vcf, txdb, CodingVariants())


## ----annotate_basics_non-synon, eval=FALSE----------------
## fa <- FaFile("data/tair10.fasta")
## predictCoding(vcf, txdb, seqSource = fa)


## ----annotate_vcf, eval=FALSE, spr=TRUE-------------------
## appendStep(sal) <- LineWise(
##     code = {
##         # comment the next line if optional step "filter_vcf" is included
##         vcf_filter <- getColumn(sal, "create_vcf")
##         # uncomment the next line if optional step "filter_vcf" is included
##         # copyEnvir(sal, "vcf_filter", globalenv())
##         library("GenomicFeatures")
##         txdb <- loadDb("./data/tair10.sqlite")
##         fa <- FaFile("data/tair10.fasta")
##         vcf_anno <- suppressMessages(suppressWarnings(variantReport(vcf_filter, txdb = txdb, fa = fa, organism = "A. thaliana", out_dir = "results/vcf_anno")))
##     },
##     step_name = "annotate_vcf",
##     dependency = "create_vcf"
## )
## 


## ----view_annotation, eval=FALSE--------------------------
## copyEnvir(sal, "vcf_anno", globalenv())
## read.delim(vcf_anno[1])[38:40, ]


## ----combine_var, eval=FALSE, spr=TRUE--------------------
## appendStep(sal) <- LineWise(
##     code = {
##         combineDF <- combineVarReports(vcf_anno, filtercol = c(Consequence = "nonsynonymous"))
##         write.table(combineDF, "./results/combineDF_nonsyn.tsv", quote = FALSE, row.names = FALSE, sep = "\t")
##     },
##     step_name = "combine_var",
##     dependency = "annotate_vcf"
## )


## ----summary_var, eval=FALSE, spr=TRUE--------------------
## appendStep(sal) <- LineWise(
##     code = {
##         write.table(varSummary(vcf_anno), "./results/variantStats.tsv", quote = FALSE, col.names = NA, sep = "\t")
##     },
##     step_name = "summary_var",
##     dependency = "combine_var"
## )


## ----venn_diagram, eval=FALSE, spr=TRUE-------------------
## appendStep(sal) <- LineWise(
##     code = {
##         ## make a list of first three samples
##         varlist <- sapply(names(vcf_anno[1:3]), function(x) as.character(read.delim(vcf_anno[x])$VARID))
##         vennset <- overLapper(varlist, type = "vennsets")
##         pdf("./results/vennplot_var.pdf")
##         vennPlot(list(vennset), mymain = "Venn Plot of First 3 Samples", mysub = "", colmode = 2, ccol = c("red", "blue"))
##         dev.off()
##     },
##     step_name = "venn_diagram",
##     dependency = "annotate_vcf"
## )


## ----plot_variant, eval=FALSE, spr=TRUE-------------------
## appendStep(sal) <- LineWise(
##     code = {
##         # comment the next line if optional step "filter_vcf" is included
##         vcf_filter <- getColumn(sal, "create_vcf")
##         # uncomment the next line if optional step "filter_vcf" is included
##         # copyEnvir(sal, "vcf_filter", globalenv())
##         library(ggbio)
##         library(VariantAnnotation)
##         mychr <- "ChrM"
##         mystart <- 19000
##         myend <- 21000
##         bams <- getColumn(sal, "fix_tag")
##         vcf <- suppressWarnings(readVcf(vcf_filter["M6B"], "A. thaliana"))
##         ga <- readGAlignments(bams["M6B"], use.names = TRUE, param = ScanBamParam(which = GRanges(mychr, IRanges(mystart, myend))))
##         p1 <- autoplot(ga, geom = "rect")
##         p2 <- autoplot(ga, geom = "line", stat = "coverage")
##         p3 <- autoplot(vcf[seqnames(vcf) == mychr], type = "fixed") +
##             xlim(mystart, myend) +
##             theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())
##         p4 <- autoplot(loadDb("./data/tair10.sqlite"), which = GRanges(mychr, IRanges(mystart, myend)), names.expr = "gene_id")
##         p1_4 <- tracks(Reads = p1, Coverage = p2, Variant = p3, Transcripts = p4, heights = c(0.3, 0.2, 0.1, 0.35)) + ylab("")
##         ggbio::ggsave(p1_4, filename = "./results/plot_variant.png", units = "in")
##     },
##     step_name = "plot_variant",
##     dependency = "create_vcf"
## )


## ----sessionInfo------------------------------------------
sessionInfo()


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

