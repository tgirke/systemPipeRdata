## External systemPipeRdata resources
##
## The large example resources formerly distributed under inst/extdata
## are stored externally and retrieved through BiocFileCache.

## Zenodo record:
## DOI: 10.5281/zenodo.22983364
## Version: 1.0.0

.SPRDATA_VERSION <- "1.0.0"

## Zenodo URLs 
.SPRDATA_URL <- c(
    data = paste0(
        "https://zenodo.org/records/22983364/files/",
        "systemPipeRdata-data-1.0.0.tar.gz?download=1"
    ),
    bam = paste0(
        "https://zenodo.org/records/22983364/files/",
        "systemPipeRdata-bam-1.0.0.tar.gz?download=1"
    )
)

## Published MD5 checksums for the Zenodo archives
.SPRDATA_MD5 <- c(
    data = "2c1bee287632937e968d53160e276d06",
    bam = "6720efe1b99f674a19303f651d670b7d"
)


#' Retrieve external systemPipeRdata resources
#'
#' Retrieves and caches the external example resources used by
#' \pkg{systemPipeRdata}. FASTQ and annotation resources are always made
#' available. Precomputed BAM files are retrieved only when
#' \code{bam = TRUE}.
#'
#' Downloaded archives are managed by \pkg{BiocFileCache}. Extracted
#' resources are stored in a versioned systemPipeRdata cache directory.
#'
#' @param bam Logical. If \code{TRUE}, also retrieve the precomputed BAM
#'   files. Default is \code{FALSE}.
#'
#' @return A list with paths to the FASTQ and annotation directories.
#'   If \code{bam = TRUE}, the list also contains the BAM directory.
#'
#' @examples
#' ## The first call downloads ~175MB (~260MB with bam = TRUE) from Zenodo;
#' ## later calls use the BiocFileCache copy.
#' if (interactive()) {
#'     paths <- getSPRdata()
#'     paths$fastqdir
#'     paths$annotationdir
#'
#'     paths <- getSPRdata(bam = TRUE)
#'     paths$bamdir
#' }
#'
#' @export
getSPRdata <- function(bam = FALSE) {

    if (!is.logical(bam) || length(bam) != 1L || is.na(bam)) {
        stop("'bam' must be TRUE or FALSE.", call. = FALSE)
    }

    if (!requireNamespace("BiocFileCache", quietly = TRUE)) {
        stop(
            "Package 'BiocFileCache' is required. ",
            "Install it with BiocManager::install('BiocFileCache').",
            call. = FALSE
        )
    }

    data_root <- .spr_get_resource("data")

    out <- list(
        fastqdir = file.path(data_root, "fastq"),
        annotationdir = file.path(data_root, "annotation"),
        bamdir = NULL
    )

    if (bam) {
        bam_root <- .spr_get_resource("bam")
        out$bamdir <- file.path(bam_root, "bam")
    }

    required <- unlist(out, use.names = FALSE)
    required <- required[!is.na(required)]

    if (!all(dir.exists(required))) {
        stop(
            "One or more systemPipeRdata resource directories ",
            "are missing after extraction.",
            call. = FALSE
        )
    }

    out
}


.spr_get_resource <- function(resource) {

    if (!resource %in% c("data", "bam")) {
        stop("Unknown systemPipeRdata resource: ", resource, call. = FALSE)
    }

    archive <- paste0(
        "systemPipeRdata-",
        resource,
        "-",
        .SPRDATA_VERSION,
        ".tar.gz"
    )

    ## During development, local archives can be supplied through:
    ##
    ## options(
    ##     systemPipeRdata.resource.dir =
    ##         "/path/to/systemPipeRdata/zenodo"
    ## )
    ##
    resource_dir <- getOption("systemPipeRdata.resource.dir")

    if (!is.null(resource_dir)) {

        source <- file.path(resource_dir, archive)

        if (!file.exists(source)) {
            stop(
                "Local systemPipeRdata resource archive not found: ",
                source,
                call. = FALSE
            )
        }

        rtype <- "local"

    } else {

        source <- .SPRDATA_URL[[resource]]

        if (is.na(source)) {
            stop(
                "The Zenodo URL for the '",
                resource,
                "' resource has not yet been configured.",
                call. = FALSE
            )
        }

        rtype <- "web"
    }

    bfc <- BiocFileCache::BiocFileCache(ask = FALSE)

    rname <- paste(
        "systemPipeRdata",
        .SPRDATA_VERSION,
        resource,
        sep = "::"
    )

    hit <- BiocFileCache::bfcquery(
        bfc,
        rname,
        field = "rname",
        exact = TRUE
    )

    if (nrow(hit) == 0L) {

        archive_path <- BiocFileCache::bfcadd(
            bfc,
            rname = rname,
            fpath = source,
            rtype = rtype,
            action = "copy",
            download = TRUE
        )

    } else {

        archive_path <- BiocFileCache::bfcrpath(
            bfc,
            rids = hit$rid[[1L]]
        )
    }

    extract_root <- file.path(
        tools::R_user_dir(
            "systemPipeRdata",
            which = "cache"
        ),
        .SPRDATA_VERSION,
        resource
    )

    complete_file <- file.path(
        extract_root,
        ".complete"
    )

    if (!file.exists(complete_file)) {

        ## Verify archive integrity before extraction. This is only done
        ## when extraction is required, avoiding repeated checksum
        ## calculation for already extracted resources.
        md5 <- unname(tools::md5sum(archive_path))

        if (!identical(md5, .SPRDATA_MD5[[resource]])) {
            stop(
                "Checksum verification failed for systemPipeRdata resource '",
                resource,
                "'.",
                call. = FALSE
            )
        }

        ## Remove an incomplete extraction from an interrupted
        ## previous operation.
        if (dir.exists(extract_root)) {
            unlink(
                extract_root,
                recursive = TRUE,
                force = TRUE
            )
        }

        dir.create(
            extract_root,
            recursive = TRUE,
            showWarnings = FALSE
        )

        utils::untar(
            archive_path,
            exdir = extract_root
        )

        if (!file.create(complete_file)) {
            stop(
                "Unable to mark systemPipeRdata resource extraction ",
                "as complete.",
                call. = FALSE
            )
        }
    }

    extract_root
}
