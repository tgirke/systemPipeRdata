###########################################
## Helper Functions for VAR-Seq Workflow ##
###########################################
## Author: Le Zhang (w. some mods by Thomas Girke)
## Last update: Dec 2025

## Download example FASTQ files specified in targets file (e.g. targetsPE_varseq.txt)
varseq_example_fastq <- function(targets) {
	build_wget <- function(url_base, file_path) {
		url_clean <- sub("/+$$", "", url_base)
		src <- paste0(url_clean, "/", basename(file_path))
		sprintf("download.file('%s', '%s')", src, file_path)
	}
	download_pairs <- Map(function(url, f1, f2) {
		c(build_wget(url, f1), build_wget(url, f2))
	}, targets$url, targets$FileName1, targets$FileName2)
	commands <- unlist(download_pairs, use.names = FALSE)
    return(commands)
}

## Download reference genome
download_ref <- function(ref="hg38.fa.gz") {
	# Download hg38 reference genome
	if(ref=="hg38.fa.gz") {
	    download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz", "./data/hg38.fa.gz")
	    R.utils::gunzip("./data/hg38.fa.gz")
	} else {
        stop("Unsupported ref value.")
    }
}

## Download software (here snpEff)
download_tool <- function(tool="snpEff_latest") {
    # Download SnpEff
	# the tool path can be used with `java -jar snpEff/snpEff.jar`"
	if(tool=="snpEff_latest") {
		download.file("https://snpeff.odsp.astrazeneca.com/versions/snpEff_latest_core.zip", "snpEff_latest_core.zip")
		unzip("snpEff_latest_core.zip")
		unlink("snpEff_latest_core.zip")
    } else {
        stop("Unsupported tool value.")
    }
}

extract_var_table <- function(vcf_vranges) {
    ann_cols <- c(
        "allele", "consequence", "effect", "gene", "gene_id", "feature_type",
        "feature_id", "transcript_biotype", "rank", "hgvs_c", "hgvs_p",
        "cdna", "cds", "aa", "distance", "warnings"
    )

    extract_ann_strings <- function(ann_obj, n) {
        if (is.null(ann_obj)) {
            return(rep(NA_character_, n))
        }
        if (inherits(ann_obj, "CompressedList") || is.list(ann_obj)) {
            ann_list <- as.list(ann_obj)
            return(vapply(ann_list, function(x) {
                if (length(x)) x[1] else NA_character_
            }, character(1)))
        }
        if (is.character(ann_obj)) {
            return(vapply(ann_obj, function(x) {
                if (is.na(x) || !nzchar(x)) return(NA_character_)
                strsplit(x, ",", fixed = TRUE)[[1]][1]
            }, character(1)))
        }
        rep(NA_character_, n)
    }

    expand_ann_fields <- function(strings) {
        if (!length(strings)) {
            return(matrix(NA_character_, nrow = 0, ncol = length(ann_cols), dimnames = list(NULL, ann_cols)))
        }
        mats <- lapply(strings, function(entry) {
            row <- rep(NA_character_, length(ann_cols))
            if (!is.na(entry) && nzchar(entry)) {
                tokens <- strsplit(entry, "\\|", fixed = FALSE)[[1]]
                row[seq_len(min(length(tokens), length(ann_cols)))] <- tokens[seq_len(min(length(tokens), length(ann_cols)))]
            }
            row
        })
        mat <- do.call(rbind, mats)
        colnames(mat) <- ann_cols
        mat
    }

    summarize_sample <- function(sample_id, vr) {
        if (!length(vr)) return(NULL)
        ann_vec <- extract_ann_strings(S4Vectors::mcols(vr)$ANN, length(vr))
        ann_df <- as.data.frame(expand_ann_fields(ann_vec), stringsAsFactors = FALSE)
        if (!nrow(ann_df)) return(NULL)
        data.frame(
            sample = sample_id,
            seqnames = as.character(GenomicRanges::seqnames(vr)),
            start = BiocGenerics::start(vr),
            end = BiocGenerics::end(vr),
            ref = as.character(VariantAnnotation::ref(vr)),
            alt = as.character(VariantAnnotation::alt(vr)),
            gene = ann_df$gene,
            consequence = ann_df$consequence,
            effect = ann_df$effect,
            stringsAsFactors = FALSE
        )
    }

    variant_tables <- Map(summarize_sample, names(vcf_vranges), vcf_vranges)
    variant_tables <- Filter(function(x) !is.null(x) && nrow(x), variant_tables)

    dplyr::bind_rows(variant_tables) |>
        dplyr::as_tibble() |>
        dplyr::mutate(location = paste0(seqnames, ":", start)) |>
        dplyr::relocate(location, .after = "sample")
}

plot_summary_filter_plot <- function(sample_filter_summary) {
	ggplot(sample_filter_summary, aes(x = sample)) +
		geom_bar(aes(y = passed, fill = "Passed"), stat = "identity") +
		geom_bar(aes(y = -filtered, fill = "Filtered"), stat = "identity") +
		coord_flip() +
		labs(y = "Number of Variants", fill = "Variant Status", title = "Variant Filtering Summary per Sample") +
		theme_minimal() +
		scale_fill_manual(values = c("Passed" = "steelblue", "Filtered" = "salmon")) +
		theme(
			plot.title = element_text(hjust = 0.5),
			axis.text.y = element_text(size = 8)
		)
}

plot_summary_boxplot <- function(plot_summary_data, sample_col = "sample", group_col = "sex") {
    boxplot_data <- plot_summary_data |>
        dplyr::count(dplyr::across(dplyr::all_of(c(sample_col, group_col))), name = "n_variants")

    p_label <- tryCatch({
        if (dplyr::n_distinct(boxplot_data[[group_col]]) < 2) return("Wilcoxon test not applicable")
        p_val <- stats::wilcox.test(n_variants ~ get(group_col), data = boxplot_data)$p.value
        paste0("Wilcoxon p = ", signif(p_val, 3))
    }, error = function(...) "Wilcoxon test failed")

    label_y <- max(boxplot_data$n_variants, na.rm = TRUE) * 1.1

    ggplot(boxplot_data, aes(x = get(group_col), y = n_variants, fill = get(group_col))) +
        geom_boxplot(alpha = 0.6, outlier.shape = NA) +
        geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
        labs(
            title = paste0("High-impact variants by ", group_col),
            x = group_col,
            y = "Variant count"
        ) +
        scale_fill_brewer(palette = "Set2", guide = "none") +
        annotate("text", x = 1.5, y = label_y, label = p_label, fontface = "bold") +
        theme_minimal() +
        expand_limits(y = label_y * 1.05)
}

simplify_info <- function(vcf_obj) {
	# Drop list-like INFO fields so VRanges coercion receives plain vectors
	info_df <- VariantAnnotation::info(vcf_obj)
	if (!ncol(info_df)) return(vcf_obj)
	keep_idx <- vapply(as.list(info_df), function(col) is.atomic(col) && !is.list(col), logical(1))
	if (any(keep_idx)) {
		info(vcf_obj) <- info_df[, keep_idx, drop = FALSE]
	} else {
		info(vcf_obj) <- S4Vectors::DataFrame()
	}
	vcf_obj
}
normalize_ft <- function(vcf_obj) {
	# Ensure genotype FT matrix stores plain character strings per sample
	ft <- VariantAnnotation::geno(vcf_obj)$FT
	if (is.null(ft) || !length(ft)) return(vcf_obj)
	ft_vec <- as.character(ft)
	dims <- dim(ft)
	if (is.null(dims)) {
		ft_mat <- matrix(ft_vec, ncol = 1)
		colnames(ft_mat) <- colnames(vcf_obj)
	} else {
		ft_mat <- matrix(ft_vec, nrow = dims[1], dimnames = dimnames(ft))
	}
	VariantAnnotation::geno(vcf_obj)$FT <- ft_mat
	vcf_obj
}

