## pre code {
## white-space: pre !important;
## overflow-x: scroll !important;
## word-break: keep-all !important;
## word-wrap: initial !important;
## }

## ----style, echo = FALSE, results = 'asis'--------------------------------------------------------------------------------------------------------------
BiocStyle::markdown()
options(width = 60, max.print = 1000)
knitr::opts_chunk$set(
    eval = as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache = as.logical(Sys.getenv("KNITR_CACHE", "TRUE")),
    tidy.opts = list(width.cutoff = 60), tidy = TRUE
)


## ----setup, echo=FALSE, message=FALSE, warning=FALSE, eval=TRUE-----------------------------------------------------------------------------------------
suppressPackageStartupMessages({
    library(systemPipeR)
})


## ----download_commands, eval=TRUE-----------------------------------------------------------------------------------------------------------------------
targets <- read.delim("targetsPE_varseq.txt", comment.char = "#")

build_wget <- function(url_base, file_path) {
    url_clean <- sub("/+$$", "", url_base)
    src <- paste0(url_clean, "/", basename(file_path))
    dest <- file_path
    sprintf("wget %s -O %s", src, dest)
}

commands <- unlist(Map(function(url, f1, f2) {
    c(build_wget(url, f1), build_wget(url, f2))
}, targets$url, targets$FileName1, targets$FileName2))

cat(paste0(paste0(commands, " &\n", collapse = ""), "wait\n"))


## # hg38 reference gnome
## wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
## gunzip hg38.fa.gz
## 
## # Download SnpEff
## wget https://snpeff.odsp.astrazeneca.com/versions/snpEff_latest_core.zip
## unzip snpEff_latest_core.zip
## rm snpEff_latest_core.zip
## # the tool path can be used with `java -jar snpEff/snpEff.jar`

## ----generate_workenvir, eval=FALSE---------------------------------------------------------------------------------------------------------------------
# library(systemPipeRdata)
# genWorkenvir(workflow = "varseq", mydirname = "varseq")
# setwd("varseq")


## ----load_targets_file, eval=TRUE-----------------------------------------------------------------------------------------------------------------------
targetspath <- system.file("extdata", "workflows", "varseq", "targetsPE_varseq.txt", package = "systemPipeRdata")
targets <- read.delim(targetspath, comment.char = "#")
targets[1:4, -(5:8)]


## ----project_varseq, eval=FALSE-------------------------------------------------------------------------------------------------------------------------
# library(systemPipeR)
# sal <- SPRproject()
# sal <- importWF(sal, file_path = "systemPipeVARseq.Rmd", verbose = FALSE)
# sal


## ----run_varseq, eval=FALSE-----------------------------------------------------------------------------------------------------------------------------
# sal <- runWF(sal)


## ----plot_varseq, eval=FALSE----------------------------------------------------------------------------------------------------------------------------
# plotWF(sal)


## ----varseq-toplogy, eval=TRUE, warning= FALSE, echo=FALSE, out.width="100%", fig.align = "center", fig.cap= "Toplogy graph of VAR-Seq workflow.", warning=FALSE----
knitr::include_graphics("results/plotwf_varseq.png")


## ----report_varseq, eval=FALSE--------------------------------------------------------------------------------------------------------------------------
# # Scientific report
# sal <- renderReport(sal)
# rmarkdown::render("systemPipeVARseq.Rmd", clean = TRUE, output_format = "BiocStyle::html_document")
# 
# # Technical (log) report
# sal <- renderLogs(sal)


## ----status_varseq, eval=FALSE--------------------------------------------------------------------------------------------------------------------------
# statusWF(sal)


## ----load_SPR, message=FALSE, eval=FALSE, spr=TRUE------------------------------------------------------------------------------------------------------
# cat(crayon::blue$bold("To use this workflow, following R packages are expected:\n"))
# cat(c("'ggplot2', 'dplyr'\n"), sep = "', '")
# ###pre-end
# appendStep(sal) <- LineWise(
#     code = {
#         library(systemPipeR)
#         },
#     step_name = "load_SPR"
# )


## ----fastqc, eval=FALSE, spr=TRUE-----------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "fastqc",
#     targets = "targetsPE_varseq.txt",
#     wf_file = "fastqc/workflow_fastqc.cwl",
#     input_file = "fastqc/fastqc.yml",
#     dir_path = "param/cwl",
#     inputvars = c(
#         FileName1  = "_FASTQ_PATH1_",
#         FileName2  = "_FASTQ_PATH2_"
#     ),
#     dependency = "load_SPR"
# )


## ----trimmomatic, eval=FALSE, spr=TRUE------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "trimmomatic",
#     targets = "targetsPE_varseq.txt",
#     wf_file = "trimmomatic/trimmomatic-pe.cwl",
#     input_file = "trimmomatic/trimmomatic-pe.yml",
#     dir_path = "param/cwl",
#     inputvars = c(
#         FileName1 = "_FASTQ_PATH1_",
#         FileName2 = "_FASTQ_PATH2_",
#         SampleName = "_SampleName_"
#     ),
#     dependency = c("load_SPR"),
#     run_step = "optional"
# )


## ----preprocessing, message=FALSE, eval=FALSE, spr=TRUE-------------------------------------------------------------------------------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "preprocessing",
#     targets = "targetsPE_varseq.txt", dir = TRUE,
#     wf_file = "preprocessReads/preprocessReads-pe.cwl",
#     input_file = "preprocessReads/preprocessReads-pe.yml",
#     dir_path = "param/cwl",
#     inputvars = c(
#         FileName1 = "_FASTQ_PATH1_",
#         FileName2 = "_FASTQ_PATH2_",
#         SampleName = "_SampleName_"
#     ),
#     dependency = c("load_SPR"),
#     run_step = "optional"
# )


## ----custom_preprocessing_function, eval=FALSE----------------------------------------------------------------------------------------------------------
# appendStep(sal) <- LineWise(
#     code = {
#         filterFct <- function(fq, cutoff = 20, Nexceptions = 0) {
#             qcount <- rowSums(as(quality(fq), "matrix") <= cutoff, na.rm = TRUE)
#             # Retains reads where Phred scores are >= cutoff with N exceptions
#             fq[qcount <= Nexceptions]
#         }
#         save(list = ls(), file = "param/customFCT.RData")
#     },
#     step_name = "custom_preprocessing_function",
#     dependency = "preprocessing"
# )


## ----editing_preprocessing, message=FALSE, eval=FALSE---------------------------------------------------------------------------------------------------
# yamlinput(sal, "preprocessing")$Fct
# yamlinput(sal, "preprocessing", "Fct") <- "'filterFct(fq, cutoff=20, Nexceptions=0)'"
# yamlinput(sal, "preprocessing")$Fct ## check the new function
# cmdlist(sal, "preprocessing", targets = 1) ## check if the command line was updated with success


## ----bwa_index, eval=FALSE, spr=TRUE--------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "bwa_index",
#     dir = FALSE, targets = NULL,
#     wf_file = "gatk/workflow_bwa-index.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     dependency = "load_SPR"
# )


## ----fasta_index, eval=FALSE, spr=TRUE------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "fasta_index",
#     dir = FALSE, targets = NULL,
#     wf_file = "gatk/workflow_fasta_dict.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     dependency = "bwa_index"
# )


## ----faidx_index, eval=FALSE, spr=TRUE------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "faidx_index",
#     dir = FALSE, targets = NULL,
#     wf_file = "gatk/workflow_fasta_faidx.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     dependency = "fasta_index"
# )


## ----bwa_alignment, eval=FALSE, spr=TRUE----------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "bwa_alignment",
#     targets = "targetsPE_varseq.txt",
#     wf_file = "gatk/workflow_bwa-pe.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     inputvars = c(
#         FileName1 = "_FASTQ_PATH1_",
#         FileName2 = "_FASTQ_PATH2_",
#         SampleName = "_SampleName_"
#     ),
#     dependency = c("faidx_index")
# )


## ----align_stats, eval=FALSE, spr=TRUE------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- LineWise(
#     code = {
#         bampaths <- getColumn(sal, step = "bwa_alignment", "outfiles", column = "samtools_sort_bam")
#         fqpaths <- getColumn(sal, step = "bwa_alignment", "targetsWF", column = "FileName1")
#         read_statsDF <- alignStats(args = bampaths, fqpaths = fqpaths, pairEnd = TRUE)
#         write.table(read_statsDF, "results/alignStats.xls", row.names = FALSE, quote = FALSE, sep = "\t")
#     },
#     step_name = "align_stats",
#     dependency = "bwa_alignment",
#     run_step = "optional"
# )


## ----bam_urls, eval=FALSE, spr=TRUE---------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- LineWise(
#     code = {
#         bampaths <- getColumn(sal, step = "bwa_alignment", "outfiles", column = "samtools_sort_bam")
#         symLink2bam(
#             sysargs = bampaths, htmldir = c("~/.html/", "somedir/"),
#             urlbase = "http://cluster.hpcc.ucr.edu/~tgirke/",
#             urlfile = "./results/IGVurl.txt"
#         )
#     },
#     step_name = "bam_urls",
#     dependency = "bwa_alignment",
#     run_step = "optional"
# )


## ----fastq2ubam, eval=FALSE, spr=TRUE-------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "fastq2ubam",
#     targets = "targetsPE_varseq.txt",
#     wf_file = "gatk/workflow_gatk_fastq2ubam.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     inputvars = c(
#         FileName1 = "_FASTQ_PATH1_",
#         FileName2 = "_FASTQ_PATH2_",
#         SampleName = "_SampleName_"
#     ),
#     dependency = c("faidx_index")
# )


## ----merge_bam, eval=FALSE, spr=TRUE--------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "merge_bam",
#     targets = c("bwa_alignment", "fastq2ubam"),
#     wf_file = "gatk/workflow_gatk_mergebams.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     inputvars = c(
#         bwa_men_sam = "_bwasam_",
#         ubam = "_ubam_",
#         SampleName = "_SampleName_"
#     ),
#     rm_targets_col = c("preprocessReads_1", "preprocessReads_2"),
#     dependency = c("bwa_alignment", "fastq2ubam")
# )


## ----sort, eval=FALSE, spr=TRUE-------------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "sort",
#     targets = "merge_bam",
#     wf_file = "gatk/workflow_gatk_sort.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     inputvars = c(merge_bam = "_mergebam_", SampleName = "_SampleName_"),
#     rm_targets_col = c(
#         "bwa_men_sam", "ubam", "SampleName_fastq2ubam",
#         "Factor_fastq2ubam", "SampleLong_fastq2ubam",
#         "Experiment_fastq2ubam", "Date_fastq2ubam"
#     ),
#     dependency = c("merge_bam")
# )


## ----mark_dup, eval=FALSE, spr=TRUE---------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "mark_dup",
#     targets = "sort",
#     wf_file = "gatk/workflow_gatk_markduplicates.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     inputvars = c(sort_bam = "_sort_", SampleName = "_SampleName_"),
#     rm_targets_col = c("merge_bam"),
#     dependency = c("sort")
# )


## ----fix_tag, eval=FALSE, spr=TRUE----------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "fix_tag",
#     targets = "mark_dup",
#     wf_file = "gatk/workflow_gatk_fixtag.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     inputvars = c(mark_bam = "_mark_", SampleName = "_SampleName_"),
#     rm_targets_col = c("sort_bam"),
#     dependency = c("mark_dup")
# )


## ----hap_caller, eval=FALSE, spr=TRUE-------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "hap_caller",
#     targets = "fix_tag",
#     wf_file = "gatk/workflow_gatk_haplotypecaller.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     inputvars = c(fixtag_bam = "_fixed_", SampleName = "_SampleName_"),
#     rm_targets_col = c("mark_bam"),
#     dependency = c("fix_tag")
# )


## ----import, eval=FALSE, spr=TRUE-----------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "import",
#     targets = NULL, dir = FALSE,
#     wf_file = "gatk/workflow_gatk_genomicsDBImport.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     dependency = c("hap_caller")
# )


## ----call_variants, eval=FALSE, spr=TRUE----------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "call_variants",
#     targets = NULL, dir = FALSE,
#     wf_file = "gatk/workflow_gatk_genotypeGVCFs.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     dependency = c("import")
# )


## ----filter, eval=FALSE, spr=TRUE-----------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "filter",
#     targets = NULL, dir = FALSE,
#     wf_file = "gatk/workflow_gatk_variantFiltration.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     dependency = c("call_variants")
# )


## ----create_vcf, eval=FALSE, spr=TRUE-------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "create_vcf",
#     targets = "hap_caller",
#     wf_file = "gatk/workflow_gatk_select_variant.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     inputvars = c(SampleName = "_SampleName_"),
#     dependency = c("hap_caller", "filter")
# )


## ----create_vcf_BCFtool, eval=FALSE, spr=TRUE-----------------------------------------------------------------------------------------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "create_vcf_BCFtool",
#     targets = "bwa_alignment", dir = TRUE,
#     wf_file = "workflow-bcftools/workflow_bcftools.cwl",
#     input_file = "workflow-bcftools/bcftools.yml",
#     dir_path = "param/cwl",
#     inputvars = c(bwa_men_sam = "_bwasam_", SampleName = "_SampleName_"),
#     rm_targets_col = c("preprocessReads_1", "preprocessReads_2"),
#     dependency = "bwa_alignment",
#     run_step = "optional"
# )


## ----inspect_vcf, eval=FALSE----------------------------------------------------------------------------------------------------------------------------
# library(VariantAnnotation)
# vcf_raw <- getColumn(sal, "create_vcf")
# vcf <- readVcf(vcf_raw[1], "Homo sapiens")
# vcf
# vr <- as(vcf, "VRanges")
# vr


## ----filter_vcf, eval=FALSE, spr=TRUE-------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- LineWise(
#     code = {
#         vcf_raw <- getColumn(sal, "create_vcf")
#         library(VariantAnnotation)
#         filter <- "totalDepth(vr) >= 20 & (altDepth(vr) / totalDepth(vr) >= 0.8)"
#         vcf_filter <- suppressWarnings(filterVars(vcf_raw, filter, organism = "Homo sapiens", out_dir = "results/vcf_filter"))
#     },
#     step_name = "filter_vcf",
#     dependency = "create_vcf",
#     run_step = "optional"
# )


## ----filter_vcf_BCFtools, eval=FALSE, spr=TRUE----------------------------------------------------------------------------------------------------------
# appendStep(sal) <- LineWise(
#     code = {
#         vcf_raw <- getColumn(sal, step = "create_vcf_BCFtool",
#                              position = "outfiles", column = "bcftools_call")
#         library(VariantAnnotation)
#         filter <- "rowSums(vr) >= 2 & (rowSums(vr[,3:4])/rowSums(vr[,1:4]) >= 0.8)"
#         vcf_filter_bcf <- suppressWarnings(filterVars(vcf_raw, filter, organism = "Homo sapiens", out_dir = "results/vcf_filter_BCFtools", varcaller = "bcftools"))
# 
#         updateColumn(sal, 'create_vcf', "outfiles") <- data.frame(vcf_filter_bcf=vcf_filter_bcf)
#     },
#     step_name = "filter_vcf_BCFtools",
#     dependency = "create_vcf_BCFtool",
#     run_step = "optional"
# )


## ----check_filter, eval=FALSE---------------------------------------------------------------------------------------------------------------------------
# copyEnvir(sal, "vcf_raw", globalenv())
# copyEnvir(sal, "vcf_filter", globalenv())
# length(as(readVcf(vcf_raw[1], genome = "Homo sapiens"), "VRanges")[, 1])
# length(as(readVcf(vcf_filter[1], genome = "Homo sapiens"), "VRanges")[, 1])


## ----summary_filter, eval=FALSE, spr=TRUE---------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- LineWise(
#     code = {
#           # read in the cohort VCF file
#           vcf_all <- suppressWarnings(VariantAnnotation::readVcf("./results/samples_filter.vcf.gz", "Homo sapiens"))
# 
#           filter_values <- VariantAnnotation::filt(vcf_all)
#           filter_values[is.na(filter_values)] <- "" # ensure character comparisons work
#           overall_counts <- table(filter_values) |>
#              dplyr::as_tibble() |>
#              dplyr::arrange(dplyr::desc(n))
# 
#           vcf_all_ft <- VariantAnnotation::geno(vcf_all)$FT
#           passes_per_sample <- apply(vcf_all_ft, 2, function(x) sum(x == "PASS", na.rm = TRUE))
#           fails_per_sample <- apply(vcf_all_ft, 2, function(x) sum(x != "PASS", na.rm = TRUE))
#           sample_filter_summary <- dplyr::tibble(
#               sample = names(passes_per_sample),
#               passed = passes_per_sample,
#               filtered = fails_per_sample
#           )
# 
#           write.table(overall_counts, file = "results/summary_filter_overall.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
#           write.table(sample_filter_summary, file = "results/summary_filter_per_sample.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
# 
#           p_filter <- ggplot2::ggplot(sample_filter_summary, ggplot2::aes(x = sample)) +
#             ggplot2::geom_bar(ggplot2::aes(y = passed, fill = "Passed"), stat = "identity") +
#             ggplot2::geom_bar(ggplot2::aes(y = -filtered, fill = "Filtered"), stat = "identity") +
#             ggplot2::coord_flip() +
#             ggplot2::labs(y = "Number of Variants", fill = "Variant Status", title = "Variant Filtering Summary per Sample") +
#             ggplot2::theme_minimal() +
#             ggplot2::scale_fill_manual(values = c("Passed" = "steelblue", "Filtered" = "salmon")) +
#             ggplot2::theme(
#                 plot.title = ggplot2::element_text(hjust = 0.5),
#                 axis.text.y = ggplot2::element_text(size = 8)
#             )
#           png("results/summary_filter_plot.png", width = 800, height = 600)
#           print(p_filter)
#           dev.off()
# 
#           # clean up RAM
#           try(rm(vcf_all, vcf_all_ft, filter_values, passes_per_sample, fails_per_sample), silent = TRUE)
#           invisible(gc())
# 
#     },
#     step_name = "summary_filter",
#     dependency = "filter"
# )


## ----annotate_vcf, eval=FALSE, spr=TRUE-----------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- SYSargsList(
#   step_name = "annotate_vcf",
#   targets = "create_vcf", dir = TRUE,
#   wf_file = "gatk/snpeff.cwl",
#   input_file = "gatk/gatk.yaml",
#   dir_path = "param/cwl",
#   inputvars = c(SampleName = "_SampleName_", vcf_raw = "_vcf_raw_"),
#   dependency = c("create_vcf")
# )


## ----combine_var, eval=FALSE, spr=TRUE------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- LineWise(
#   code = {
#     vcf_anno <- getColumn(sal, "annotate_vcf", position = "outfiles", column = "ann_vcf")
# 
#     clean_vcf_file <- function(path) {
#       lines <- readLines(path, warn = FALSE)
#       if (!length(lines)) return(path)
#       header_start <- which(grepl("^##fileformat=", lines, perl = TRUE))[1]
#       if (is.na(header_start) || header_start <= 1) return(path)
#       writeLines(lines[header_start:length(lines)], path)
#       path
#     }
# 
#     vcf_anno <- vapply(vcf_anno, clean_vcf_file, FUN.VALUE = character(1))
#     vr_from_vcf <- function(path) {
#       message("Importing annotated VCF: ", path)
#       vcf <- suppressWarnings(VariantAnnotation::readVcf(path, "Homo sapiens"))
# 
#       ft <- VariantAnnotation::geno(vcf)$FT
#       if (!is.null(ft) && !isTRUE(is.character(ft))) {
#         ft_vec <- as.character(ft)
#         ft_clean <- matrix(ft_vec, ncol = 1)
#         if (!is.null(dim(ft))) {
#           ft_clean <- matrix(ft_vec, nrow = nrow(ft), dimnames = dimnames(ft))
#         }
#         VariantAnnotation::geno(vcf)$FT <- ft_clean
#       }
#       suppressWarnings(as(vcf, "VRanges"))
#     }
# 
#     vcf_vranges <- lapply(vcf_anno, vr_from_vcf)
#   },
#   step_name = "combine_var",
#   dependency = "annotate_vcf"
# )


## ----summary_var, eval=FALSE, spr=TRUE------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- LineWise(
#     code = {
#         ann_cols <- c(
#             "allele", "consequence", "effect", "gene", "gene_id", "feature_type",
#             "feature_id", "transcript_biotype", "rank", "hgvs_c", "hgvs_p",
#             "cdna", "cds", "aa", "distance", "warnings"
#         )
# 
#         extract_ann_strings <- function(ann_obj, n) {
#             if (is.null(ann_obj)) {
#                 return(rep(NA_character_, n))
#             }
#             if (inherits(ann_obj, "CompressedList") || is.list(ann_obj)) {
#                 ann_list <- as.list(ann_obj)
#                 return(vapply(ann_list, function(x) {
#                     if (length(x)) x[1] else NA_character_
#                 }, character(1)))
#             }
#             if (is.character(ann_obj)) {
#                 return(vapply(ann_obj, function(x) {
#                     if (is.na(x) || !nzchar(x)) return(NA_character_)
#                     strsplit(x, ",", fixed = TRUE)[[1]][1]
#                 }, character(1)))
#             }
#             rep(NA_character_, n)
#         }
# 
#         expand_ann_fields <- function(strings) {
#             if (!length(strings)) {
#                 return(matrix(NA_character_, nrow = 0, ncol = length(ann_cols), dimnames = list(NULL, ann_cols)))
#             }
#             mats <- lapply(strings, function(entry) {
#                 row <- rep(NA_character_, length(ann_cols))
#                 if (!is.na(entry) && nzchar(entry)) {
#                     tokens <- strsplit(entry, "\\|", fixed = FALSE)[[1]]
#                     row[seq_len(min(length(tokens), length(ann_cols)))] <- tokens[seq_len(min(length(tokens), length(ann_cols)))]
#                 }
#                 row
#             })
#             mat <- do.call(rbind, mats)
#             colnames(mat) <- ann_cols
#             mat
#         }
# 
#         summarize_sample <- function(sample_id, vr) {
#             if (!length(vr)) return(NULL)
#             ann_vec <- extract_ann_strings(S4Vectors::mcols(vr)$ANN, length(vr))
#             ann_df <- as.data.frame(expand_ann_fields(ann_vec), stringsAsFactors = FALSE)
#             if (!nrow(ann_df)) return(NULL)
#             data.frame(
#                 sample = sample_id,
#                 seqnames = as.character(GenomicRanges::seqnames(vr)),
#                 start = BiocGenerics::start(vr),
#                 end = BiocGenerics::end(vr),
#                 ref = as.character(VariantAnnotation::ref(vr)),
#                 alt = as.character(VariantAnnotation::alt(vr)),
#                 gene = ann_df$gene,
#                 consequence = ann_df$consequence,
#                 effect = ann_df$effect,
#                 stringsAsFactors = FALSE
#             )
#         }
# 
#         variant_tables <- Map(summarize_sample, names(vcf_vranges), vcf_vranges)
#         variant_tables <- Filter(function(x) !is.null(x) && nrow(x), variant_tables)
# 
#         summary_var <- dplyr::bind_rows(variant_tables) |>
#           dplyr::as_tibble() |>
#           dplyr::mutate(location = paste0(seqnames, ":", start)) |>
#           dplyr::relocate(location, .after = "sample")
#         utils::write.table(summary_var, file = "results/variant_summary_long.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
#     },
#     step_name = "summary_var",
#     dependency = "combine_var"
# )


## ----plot_var_stats, eval=FALSE, spr=TRUE---------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- LineWise(
#     code = {
#         library(ggplot2)
#         plot_summary_data <- summary_var |>
#           dplyr::filter(
#             consequence %in% c(
#               "nonsynonymous_variant", "stop_gained", "frameshift_variant",
#               "splice_acceptor_variant", "splice_donor_variant",
#               "start_lost", "stop_lost"
#             ),
#             effect == "HIGH"
#           ) |>
#           dplyr::mutate(
#             sample = factor(sample, levels = getColumn(sal, step = "annotate_vcf", position = "targetsWF", column = "SampleName")),
#             sex = getColumn(sal, step = "annotate_vcf", position = "targetsWF", column = "Factor")[sample]
#           )
#                 assign("plot_summary_data", plot_summary_data, envir = .GlobalEnv)
# 
#         png("./results/var_summary.png")
#         p_var_summary <- ggplot(plot_summary_data) +
#           geom_bar(aes(x = sample, fill = sex), alpha = 0.75) +
#           scale_fill_brewer(palette = "Set2") +
#           theme_minimal()
#         print(p_var_summary)
#         dev.off()
#     },
#     step_name = "plot_var_stats",
#     dependency = "summary_var"
# )


## ----plot_var_boxplot, eval=FALSE, spr=TRUE-------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- LineWise(
#     code = {
#         if (!exists("plot_summary_data", inherits = FALSE)) {
#             stop("'plot_summary_data' not found. Please run 'plot_var_stats' first.")
#         }
#         library(ggplot2)
#         boxplot_data <- plot_summary_data |>
#             dplyr::count(sample, sex, name = "n_variants")
# 
#         p_label <- tryCatch({
#             if (dplyr::n_distinct(boxplot_data$sex) < 2) return("Wilcoxon test not applicable")
#             p_val <- stats::wilcox.test(n_variants ~ sex, data = boxplot_data)$p.value
#             paste0("Wilcoxon p = ", signif(p_val, 3))
#         }, error = function(...) "Wilcoxon test failed")
# 
#         label_y <- max(boxplot_data$n_variants, na.rm = TRUE) * 1.1
# 
#         png("./results/var_summary_boxplot.png")
#         p_summary_boxplot <- ggplot(boxplot_data, aes(x = sex, y = n_variants, fill = sex)) +
#             geom_boxplot(alpha = 0.6, outlier.shape = NA) +
#             geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
#             labs(
#                 title = "High-impact variants by sex",
#                 x = "Sex",
#                 y = "Variant count"
#             ) +
#             scale_fill_brewer(palette = "Set2", guide = "none") +
#             annotate("text", x = 1.5, y = label_y, label = p_label, fontface = "bold") +
#             theme_minimal() +
#             expand_limits(y = label_y * 1.05)
#         print(p_summary_boxplot)
#         dev.off()
#     },
#     step_name = "plot_var_boxplot",
#     dependency = "plot_var_stats"
# )


## ----venn_diagram, eval=FALSE, spr=TRUE-----------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- LineWise(
#     code = {
#         top_n <- min(3, length(variant_tables))
#         selected_tables <- variant_tables[seq_len(top_n)]
#         if (is.null(names(selected_tables)) || any(!nzchar(names(selected_tables)))) {
#             names(selected_tables) <- paste0("Sample_", seq_along(selected_tables))
#         }
# 
#         variant_sets <- lapply(selected_tables, function(df) {
#             if (!nrow(df)) return(character(0))
#             unique(paste0(df$seqnames, ":", df$start, "_", df$ref, "/", df$alt))
#         })
# 
#         vennset <- overLapper(variant_sets, type = "vennsets")
#         png("./results/vennplot_var.png")
#         vennPlot(vennset, mymain = "Venn Plot of First 3 Samples", mysub = "", colmode = 2, ccol = c("red", "blue"))
#         dev.off()
#     },
#     step_name = "venn_diagram",
#     dependency = "summary_var"
# )


## ----plot_variant, eval=FALSE, spr=TRUE-----------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- LineWise(
#     code = {
# 
#         first_high <- summary_var |>
#           dplyr::filter(effect == "HIGH") |>
#           head(n = 1)
# 
#         library(ggbio)
#         library(VariantAnnotation)
#         mychr <- as.character(first_high$seqnames)
#         mystart <- as.numeric(first_high$start) - 500
#         myend <- as.numeric(first_high$end) + 500
#         bam_path <- getColumn(sal, "fix_tag")[first_high$sample]
#         vcf_path <- getColumn(sal, step = "create_vcf")[first_high$sample]
# 
#         vcf <- suppressWarnings(readVcf(vcf_path, "Homo sapiens"))
#         ga <- readGAlignments(bam_path, use.names = TRUE, param = ScanBamParam(which = GRanges(mychr, IRanges(mystart, myend))))
#         simplify_info <- function(vcf_obj) {
#             # Drop list-like INFO fields so VRanges coercion receives plain vectors
#             info_df <- VariantAnnotation::info(vcf_obj)
#             if (!ncol(info_df)) return(vcf_obj)
#             keep_idx <- vapply(as.list(info_df), function(col) is.atomic(col) && !is.list(col), logical(1))
#             if (any(keep_idx)) {
#                 info(vcf_obj) <- info_df[, keep_idx, drop = FALSE]
#             } else {
#                 info(vcf_obj) <- S4Vectors::DataFrame()
#             }
#             vcf_obj
#         }
#         normalize_ft <- function(vcf_obj) {
#             # Ensure genotype FT matrix stores plain character strings per sample
#             ft <- VariantAnnotation::geno(vcf_obj)$FT
#             if (is.null(ft) || !length(ft)) return(vcf_obj)
#             ft_vec <- as.character(ft)
#             dims <- dim(ft)
#             if (is.null(dims)) {
#                 ft_mat <- matrix(ft_vec, ncol = 1)
#                 colnames(ft_mat) <- colnames(vcf_obj)
#             } else {
#                 ft_mat <- matrix(ft_vec, nrow = dims[1], dimnames = dimnames(ft))
#             }
#             VariantAnnotation::geno(vcf_obj)$FT <- ft_mat
#             vcf_obj
#         }
#         vcf_chr <- normalize_ft(simplify_info(vcf[seqnames(vcf) == mychr]))
#         vr <- suppressWarnings(as(vcf_chr, "VRanges"))
#         vr_region <- vr[start(vr) >= mystart & end(vr) <= myend]
#         if (!length(vr_region)) {
#             vr_region <- vr
#         }
#         p1 <- autoplot(ga, geom = "rect")
#         p2 <- autoplot(ga, geom = "line", stat = "coverage")
#         p3 <- autoplot(vr_region, type = "fixed") +
#             xlim(mystart, myend) +
#             theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())
#         p1_3 <- tracks(
#             place_holder = ggplot2::ggplot(),
#             Reads = p1,
#             Coverage = p2,
#             Variant = p3,
#             heights = c(0, 0.3, 0.2, 0.1)
#         ) + ylab("")
#         ggbio::ggsave(p1_3, filename = "./results/plot_variant.png", units = "in")
#     },
#     step_name = "plot_variant",
#     dependency = "summary_var"
# )


## ----sessionInfo, eval=FALSE, spr=TRUE------------------------------------------------------------------------------------------------------------------
# appendStep(sal) <- LineWise(
#     code = {
#         sessionInfo()
#         },
#     step_name = "sessionInfo",
#     dependency = "plot_variant")


## ----runWF, eval=FALSE----------------------------------------------------------------------------------------------------------------------------------
# sal <- runWF(sal)


## ----runWF_cluster, eval=FALSE--------------------------------------------------------------------------------------------------------------------------
# # wall time in mins, memory in MB
# resources <- list(conffile=".batchtools.conf.R",
#                   template="batchtools.slurm.tmpl",
#                   Njobs=8,
#                   walltime=120,
#                   ntasks=1,
#                   ncpus=4,
#                   memory=1024,
#                   partition = "short"
#                   )
# sal <- addResources(sal, c("hisat2_mapping"), resources = resources)
# sal <- runWF(sal)


## ----plotWF, eval=FALSE---------------------------------------------------------------------------------------------------------------------------------
# plotWF(sal, rstudio = TRUE)


## ----statusWF, eval=FALSE-------------------------------------------------------------------------------------------------------------------------------
# sal
# statusWF(sal)


## ----logsWF, eval=FALSE---------------------------------------------------------------------------------------------------------------------------------
# sal <- renderLogs(sal)


## ----list_tools, eval=FALSE-----------------------------------------------------------------------------------------------------------------------------
# if(file.exists(file.path(".SPRproject", "SYSargsList.yml"))) {
#     local({
#         sal <- systemPipeR::SPRproject(resume = TRUE)
#         systemPipeR::listCmdTools(sal)
#         systemPipeR::listCmdModules(sal)
#     })
# } else {
#     cat(crayon::blue$bold("Tools and modules required by this workflow are:\n"))
#     cat(c("trimmomatic/0.39", "samtools/1.14", "gatk/4.2.0.0", "bcftools/1.15",
#           "bwa/0.7.17", "snpEff/5.3"), sep = "\n")
# }


## ----report_session_info, eval=TRUE---------------------------------------------------------------------------------------------------------------------
sessionInfo()

