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


## ----generate_workenvir, eval=FALSE-----------------------
# library(systemPipeRdata)
# genWorkenvir(workflow = "varseq", mydirname = "varseq")
# setwd("varseq")


## ----download_commands, eval=FALSE------------------------
# targets <- read.delim(system.file("extdata", "workflows", "varseq", "targetsPE_varseq.txt", package = "systemPipeRdata"), comment.char = "#")
# source("VARseq_helper.R") # defines helper functions
# options(timeout = 3600) # increase time limit for downloads
# commands <- varseq_example_fastq(targets)
# for (cmd in commands) eval(parse(text = cmd))
# download_ref(ref="hg38.fa.gz")
# download_tool(tool="snpEff_latest")


## ----load_targets_file, eval=TRUE-------------------------
targetspath <- system.file("extdata", "workflows", "varseq", "targetsPE_varseq.txt", package = "systemPipeRdata")
targets <- read.delim(targetspath, comment.char = "#")
targets[1:4, -(5:8)]


## ----project_varseq, eval=FALSE---------------------------
# library(systemPipeR)
# sal <- SPRproject()
# # sal <- SPRproject(resume=TRUE, load.envir=TRUE) # to resume workflow if needed
# sal <- importWF(sal, file_path = "systemPipeVARseq.Rmd", verbose = FALSE)
# sal


## ----run_varseq, eval=FALSE-------------------------------
# sal <- runWF(sal)


## ----plot_varseq, eval=FALSE------------------------------
# plotWF(sal)


## ----varseq-toplogy, eval=TRUE, warning= FALSE, echo=FALSE, out.width="100%", fig.align = "center", fig.cap= "Topology graph of VAR-Seq workflow.", warning=FALSE----
knitr::include_graphics("results/plotwf_varseq.png")


## ----report_varseq, eval=FALSE----------------------------
# # Scientific report
# sal <- renderReport(sal)
# rmarkdown::render("systemPipeVARseq.Rmd", clean = TRUE, output_format = "BiocStyle::html_document")
# 
# # Technical (log) report
# sal <- renderLogs(sal)


## ----status_varseq, eval=FALSE----------------------------
# statusWF(sal)


## ----save_sal, eval=FALSE---------------------------------
# # sal <- write_SYSargsList(sal)


## ----load_SPR, message=FALSE, eval=FALSE, spr=TRUE--------
# cat(crayon::blue$bold("To use this workflow, following R packages are expected:\n"))
# cat(c("'ggplot2', 'dplyr'\n"), sep = "', '")
# ###pre-end
# appendStep(sal) <- LineWise(
#     code = {
#         library(systemPipeR)
#         },
#     step_name = "load_SPR"
# )


## ----fastqc, eval=FALSE, spr=TRUE-------------------------
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


## ----fastq_report, eval=FALSE, message=FALSE, spr=TRUE----
# appendStep(sal) <- LineWise(code = {
#   fastq1 <- getColumn(sal, step = "fastqc", "targetsWF", column = 1)
#   fastq2 <- getColumn(sal, step = "fastqc", "targetsWF", column = 2)
#   fastq <- setNames(c(rbind(fastq1, fastq2)), c(rbind(names(fastq1), names(fastq2))))
#   fqlist <- seeFastq(fastq = fastq, batchsize = 1000, klength = 8)
#   png("./results/fastqReport_varseq.png", height = 650, width = 288 * length(fqlist))
#   seeFastqPlot(fqlist)
#   dev.off()
# }, step_name = "fastq_report",
# dependency = "fastqc")


## ----trimmomatic, eval=FALSE, spr=TRUE--------------------
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


## ----preprocessing, message=FALSE, eval=FALSE, spr=TRUE----
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


## ----bwa_index, eval=FALSE, spr=TRUE----------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "bwa_index",
#     dir = FALSE, targets = NULL,
#     wf_file = "gatk/workflow_bwa-index.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     dependency = "load_SPR"
# )


## ----fasta_index, eval=FALSE, spr=TRUE--------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "fasta_index",
#     dir = FALSE, targets = NULL,
#     wf_file = "gatk/workflow_fasta_dict.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     dependency = "bwa_index"
# )


## ----faidx_index, eval=FALSE, spr=TRUE--------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "faidx_index",
#     dir = FALSE, targets = NULL,
#     wf_file = "gatk/workflow_fasta_faidx.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     dependency = "fasta_index"
# )


## ----bwa_alignment, eval=FALSE, spr=TRUE------------------
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


## ----align_stats, eval=FALSE, spr=TRUE--------------------
# appendStep(sal) <- LineWise(
#     code = {
#         bampaths <- getColumn(sal, step = "bwa_alignment", "outfiles", column = "samtools_sort_bam")
#         fqpaths <- getColumn(sal, step = "bwa_alignment", "targetsWF", column = "FileName1")
#         read_statsDF <- alignStats(args = bampaths, fqpaths = fqpaths, pairEnd = TRUE)
#         write.table(read_statsDF, "results/alignStats_varseq.xls", row.names = FALSE, quote = FALSE, sep = "\t")
#     },
#     step_name = "align_stats",
#     dependency = "bwa_alignment",
#     run_step = "optional"
# )


## ----align_stats_view, eval=TRUE--------------------------
read.table("results/alignStats_varseq.xls", header = TRUE)[1:4,]


## ----bam_urls, eval=FALSE, spr=TRUE-----------------------
# appendStep(sal) <- LineWise(
#     code = {
#         bampaths <- getColumn(sal, step = "bwa_alignment", "outfiles", column = "samtools_sort_bam")
#         symLink2bam(
#             sysargs = bampaths, htmldir = c("~/.html/", "<somedir>/"),
#             urlbase = "<base_url>/~<username>/",
#             urlfile = "./results/IGVurl.txt"
#         )
#     },
#     step_name = "bam_urls",
#     dependency = "bwa_alignment",
#     run_step = "optional"
# )


## ----fastq2ubam, eval=FALSE, spr=TRUE---------------------
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


## ----merge_bam, eval=FALSE, spr=TRUE----------------------
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


## ----sort, eval=FALSE, spr=TRUE---------------------------
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


## ----mark_dup, eval=FALSE, spr=TRUE-----------------------
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


## ----fix_tag, eval=FALSE, spr=TRUE------------------------
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


## ----hap_caller, eval=FALSE, spr=TRUE---------------------
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


## ----import, eval=FALSE, spr=TRUE-------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "import",
#     targets = NULL, dir = FALSE,
#     wf_file = "gatk/workflow_gatk_genomicsDBImport.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     dependency = c("hap_caller")
# )


## ----call_variants, eval=FALSE, spr=TRUE------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "call_variants",
#     targets = NULL, dir = FALSE,
#     wf_file = "gatk/workflow_gatk_genotypeGVCFs.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     dependency = c("import")
# )


## ----filter, eval=FALSE, spr=TRUE-------------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "filter",
#     targets = NULL, dir = FALSE,
#     wf_file = "gatk/workflow_gatk_variantFiltration.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     dependency = c("call_variants")
# )


## ----create_vcf, eval=FALSE, spr=TRUE---------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "create_vcf",
#     targets = "hap_caller",
#     wf_file = "gatk/workflow_gatk_select_variant.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     inputvars = c(SampleName = "_SampleName_"),
#     dependency = c("hap_caller", "filter")
# )


## ----inspect_vcf, eval=FALSE------------------------------
# library(VariantAnnotation)
# vcf_raw <- getColumn(sal, "create_vcf")
# vcf <- readVcf(vcf_raw[1], "Homo sapiens")
# vcf
# vr <- as(vcf, "VRanges")
# vr


## ----filter_vcf, eval=FALSE, spr=TRUE---------------------
# appendStep(sal) <- LineWise(
#     code = {
#         source("VARseq_helper.R")
#         vcf_raw <- getColumn(sal, "create_vcf")
#         library(VariantAnnotation)
#         filter <- "totalDepth(vr) >= 20 & (altDepth(vr) / totalDepth(vr) >= 0.8)"
#         vcf_filter <- suppressWarnings(filter_vars(vcf_raw, filter, organism = "Homo sapiens", out_dir = "results/vcf_filter"))
#     },
#     step_name = "filter_vcf",
#     dependency = "create_vcf",
#     run_step = "optional"
# )


## ----check_filter, eval=FALSE-----------------------------
# copyEnvir(sal, "vcf_raw", globalenv())
# copyEnvir(sal, "vcf_filter", globalenv())
# length(as(readVcf(vcf_raw[1], genome = "Homo sapiens"), "VRanges")[, 1])
# length(as(readVcf(vcf_filter[1], genome = "Homo sapiens"), "VRanges")[, 1])


## ----summary_filter, eval=FALSE, spr=TRUE-----------------
# appendStep(sal) <- LineWise(
#     code = {
# 
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
#           library(ggplot2)
#           source("VARseq_helper.R")
#           png("results/summary_filter_plot.png", width = 800, height = 600)
#           p_filter <- plot_summary_filter_plot(sample_filter_summary)
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


## ----annotate_vcf, eval=FALSE, spr=TRUE-------------------
# appendStep(sal) <- SYSargsList(
#     step_name = "annotate_vcf",
#     targets = "create_vcf", dir = TRUE,
#     wf_file = "gatk/snpeff.cwl",
#     input_file = "gatk/gatk.yaml",
#     dir_path = "param/cwl",
#     inputvars = c(SampleName = "_SampleName_", vcf_raw = "_vcf_raw_"),
#     dependency = c("create_vcf")
# )


## ----combine_var, eval=FALSE, spr=TRUE--------------------
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


## ----summary_var, eval=FALSE, spr=TRUE--------------------
# appendStep(sal) <- LineWise(
#     code = {
#         source("VARseq_helper.R")
#         summary_var <- extract_var_table(vcf_vranges)
#         utils::write.table(summary_var, file = "results/variant_summary_long.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
#     },
#     step_name = "summary_var",
#     dependency = "combine_var"
# )


## ----plot_var_consequence, eval=FALSE, spr=TRUE-----------
# appendStep(sal) <- LineWise(
#     code = {
#         library(ggplot2)
#         effect_counts <- summary_var |>
#           dplyr::count(sample, effect, name = "n")
# 
#         png("./results/var_consequence_log.png", width = 1000, height = 600)
#         p_var_consequence_log <- ggplot(effect_counts, aes(sample, n, fill = effect)) +
#         geom_col(position = "dodge", alpha = 0.8) +
#         geom_text(aes(label = scales::comma(n)), position = position_dodge(width = 0.9), vjust = -0.2, size = 3) +
#         scale_y_continuous(trans = "log10", labels = scales::comma_format()) +
#         labs(y = "Variant count (log10 scale)", title = "Variant consequences (log scale)") +
#         theme_minimal() +
#         theme(axis.text.x = element_text(angle = 45, hjust = 1))
#         print(p_var_consequence_log)
#         dev.off()
#     },
#     step_name = "plot_var_consequence",
#     dependency = "summary_var"
# )


## ----plot_var_stats, eval=FALSE, spr=TRUE-----------------
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


## ----plot_var_boxplot, eval=FALSE, spr=TRUE---------------
# appendStep(sal) <- LineWise(
#     code = {
#         source("VARseq_helper.R")
#         # change the sample and group columns as needed
#         p_summary_boxplot <- plot_summary_boxplot(plot_summary_data, sample_col = "sample", group_col = "sex")
#         library(ggplot2)
#         png("./results/var_summary_boxplot.png")
#         print(p_summary_boxplot)
#         dev.off()
#     },
#     step_name = "plot_var_boxplot",
#     dependency = "plot_var_stats"
# )


## ----venn_diagram, eval=FALSE, spr=TRUE-------------------
# appendStep(sal) <- LineWise(
#     code = {
#         unique_samples <- summary_var |> dplyr::distinct(sample) |> dplyr::pull(sample)
#         if (!length(unique_samples)) {
#             stop("No samples available in `summary_var`; cannot draw Venn diagram.")
#         }
# 
#         top_n <- min(3, length(unique_samples))
#         selected_samples <- unique_samples[seq_len(top_n)]
# 
#         variant_df <- summary_var |>
#             dplyr::filter(sample %in% selected_samples) |>
#             dplyr::distinct(sample, seqnames, start, ref, alt) |>
#             dplyr::mutate(variant_id = paste0(seqnames, ":", start, "_", ref, "/", alt))
# 
#         variant_sets <- split(variant_df$variant_id, variant_df$sample)
#         vennset <- overLapper(variant_sets, type = "vennsets")
#         png("./results/vennplot_var.png")
#         vennPlot(vennset, mymain = "Venn Plot of First 3 Samples", mysub = "", colmode = 2, ccol = c("red", "blue"))
#         dev.off()
#     },
#     step_name = "venn_diagram",
#     dependency = "summary_var"
# )


## ----plot_variant, eval=FALSE, spr=TRUE-------------------
# appendStep(sal) <- LineWise(
#     code = {
#         source("VARseq_helper.R")
#         first_high <- summary_var |>
#           dplyr::filter(effect == "HIGH") |>
#           head(n = 1)
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


## ----non_syn_vars, eval=FALSE, spr=TRUE-------------------
# appendStep(sal) <- LineWise(
#     code = {
#         vardf <- read.delim("results/variant_summary_long.tsv")
#         source("VARseq_helper.R")
#         common_nonsyn_entrez <- filterNonSyn(df=vardf)
#         writeLines(common_nonsyn_entrez, "results/common_nonsyn_entrez")
#     },
#     step_name = "non_syn_vars",
#     dependency = "plot_var_consequence"
# )


## ----pathenrich, eval=FALSE, spr=TRUE---------------------
# appendStep(sal) <- LineWise(
#     code = {
#         common_nonsyn_entrez <- readLines("results/common_nonsyn_entrez")
#         source("VARseq_helper.R")
#         reacdb <- load_reacList(org="R-HSA")
# 
#         library(fgsea); library(data.table); library(ggplot2)
#         foraRes <- fora(genes=common_nonsyn_entrez, universe=unique(unlist(reacdb)), pathways=reacdb)
#         if (!dir.exists("results/fea")) dir.create("results/fea", recursive = TRUE)
#         foraRes$overlapGenes <- vapply(foraRes$overlapGenes, toString, FUN.VALUE = character(1))
#         write.table(foraRes, file = "results/fea/foraRes.xls", row.names = FALSE, sep = "\t", quote = FALSE)
#         foraRes$pathway <- gsub("\\(.*\\) ", "", foraRes$pathway)
#         foraRes$pathway <- factor(foraRes$pathway, levels = rev(foraRes$pathway))
#         png("./results/fea/pathenrich.png", width = 680)
#         ggplot(head(foraRes, 15), aes(pathway, overlap, fill = padj)) +
#             geom_bar(position="dodge", stat="identity") + coord_flip() +
#             scale_fill_distiller(palette = "RdBu", direction=-1, limits = range(head(foraRes$padj, 15))) +
#             theme(axis.text=element_text(angle=0, hjust=1, size=12), axis.title = element_text(size = 14))
#         dev.off()
#     },
#     step_name = "pathenrich",
#     dependency = "non_syn_vars"
# )


## ----drug_target_analysis, eval=FALSE, spr=TRUE-----------
# appendStep(sal) <- LineWise(
#     code = {
#         ## Configure paths for drugTargetInteractions. Under chembldb provide path to chembl_xx.db on your system
#         # chembldb <- system.file("extdata", "chembl_sample.db", package="drugTargetInteractions")
#         resultsPath <- "results/drug_target/"
#         config <- drugTargetInteractions::genConfig(chemblDbPath=chembldb, resultsPath=resultsPath)
#         downloadUniChem(config=config)
#         cmpIdMapping(config=config)
#         foraRes <- read.delim("results/fea/foraRes.xls")
#         entrez_ids <- unlist(strsplit(foraRes[13,7], ", ")) # select here pathway of interest
#         source("VARseq_helper.R")
#         drugMap <- runGeneTargetDrug(entrez=entrez_ids)[[1]]
#         drugMap <- drugMap[!grepl("Query_", drugMap$GeneName), c("GeneName", "UniProt_ID", "Target_Desc", "Drug_Name", "CHEMBL_CMP_ID", "MOA", "Mesh_Indication")]
#         drugMap <- drugMap[!is.na(drugMap$CHEMBL_CMP_ID),]
#         if (!dir.exists("results/drug_target")) dir.create("results/drug_target", recursive = TRUE)
#         write.table(drugMap, file="results/drug_target/drug_target.xls", row.names=FALSE, sep="\t", quote=FALSE)
#     },
#     step_name = "drug_target",
#     dependency = "pathenrich"
# )


## ----read_drug_table, eval=TRUE---------------------------
drugMap <- read.delim("results/drug_target.xls")
DT::datatable(drugMap)


## ----sessionInfo, eval=FALSE, spr=TRUE--------------------
# appendStep(sal) <- LineWise(
#     code = {
#         sessionInfo()
#         },
#     step_name = "sessionInfo",
#     dependency = "drug_target")


## ----runWF, eval=FALSE------------------------------------
# sal <- runWF(sal)


## ----runWF_cluster, eval=FALSE----------------------------
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
# sal <- addResources(sal, c("bwa_alignment"), resources = resources)
# sal <- runWF(sal)


## ----plotWF, eval=FALSE-----------------------------------
# plotWF(sal, rstudio = TRUE)


## ----statusWF, eval=FALSE---------------------------------
# sal
# statusWF(sal)


## ----logsWF, eval=FALSE-----------------------------------
# sal <- renderLogs(sal)


## ----list_tools, eval=TRUE--------------------------------
if(file.exists(file.path(".SPRproject", "SYSargsList.yml"))) {
    local({
        sal <- systemPipeR::SPRproject(resume = TRUE)
        systemPipeR::listCmdTools(sal)
        systemPipeR::listCmdModules(sal)
    })
} else {
    cat(crayon::blue$bold("Tools and modules required by this workflow are:\n"))
    cat(c("trimmomatic/0.39", "samtools/1.14", "gatk/4.2.0.0", "bcftools/1.15", 
          "bwa/0.7.17", "snpEff/5.3"), sep = "\n")
}


## ----report_session_info, eval=TRUE-----------------------
sessionInfo()

