##################################################
## Initialize Workflow Environments from GitHub ##
##################################################

#########################################
## Clone workflow with gert::git_clone ##
########################################
## Optional: user can use just 'gert::git_clone' directly or 'git clone' from command-line
## The advantage of 'gert' over 'git clone' is that it eliminates the git dependency  
#' @keywords internal
.ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed. Install it via install.packages('%s').", pkg, pkg),
         call. = FALSE)
  }
}

#' Clone a workflow repository using gert and enforce version requirements
#'
#' Convenience wrapper around \code{gert::git_clone()} that clones a workflow
#' repository and enforces minimum package version requirements specified in
#' the workflow's \code{manifest.yml} file.
#'
#' If \code{manifest.yml} contains a block like:
#'
#' \preformatted{
#' systemPipeR:
#'   bioc_min: "2.17.1"
#'
#' systemPipeRdata:
#'   bioc_min: "2.15.4"
#' }
#'
#' then the installed versions of the corresponding packages are checked.
#'
#' By default, version mismatches cause an error and the cloned directory
#' is removed to avoid leaving a partially initialized workflow.
#'
#' @param url GitHub HTTPS URL of the workflow repository to clone.
#' @param mydirname Target directory to clone into.
#' @param force_min_version Logical; if \code{FALSE} (default), enforce
#'   \code{bioc_min} requirements and stop on version mismatch. If \code{TRUE},
#'   continue with a warning (at your own risk). Missing required packages
#'   always cause an error.
#' @param cleanup_on_fail Logical; if \code{TRUE} (default) and a strict
#'   version check fails, remove the cloned directory before stopping.
#' @param ... Additional arguments passed to \code{gert::git_clone()}.
#'
#' @return Invisibly returns the normalized path of the cloned repository.
#'
#' @examples
#' \dontrun{
#' ## Strict (default)
#' genWorkenvir_gh(
#'   "https://github.com/systemPipeR/sprwf-new.git",
#'   "sprwf-new"
#' )
#'
#' ## Allow version mismatch (not recommended)
#' genWorkenvir_gh(
#'   "https://github.com/systemPipeR/sprwf-new.git",
#'   "sprwf-new",
#'   force_min_version = TRUE
#' )
#' }
genWorkenvir_gh <- function(url,
                            mydirname,
                            force_min_version = FALSE,
                            cleanup_on_fail = TRUE,
                            ...) {

  .ensure_pkg("gert")

  if (dir.exists(mydirname) &&
      length(list.files(mydirname, all.files = TRUE, no.. = TRUE)) > 0) {
    stop(sprintf("Target mydirname '%s' exists and is not empty.",
                 mydirname),
         call. = FALSE)
  }

  # Clone first
  gert::git_clone(url = url, path = mydirname, ...)
  repo_path <- normalizePath(mydirname, winslash = "/", mustWork = TRUE)

  # Now enforce min versions transactionally
  manifest <- file.path(repo_path, "manifest.yml")

  result <- tryCatch({

    mins <- .read_bioc_min_from_manifest(
      manifest = manifest,
      pkgs = c("systemPipeR", "systemPipeRdata")
    )

    .check_min_pkg_versions(
      mins,
      context = sprintf("'%s'", manifest),
      force = force_min_version,
      quiet = FALSE
    )

    TRUE  # success flag

  }, error = function(e) {

    if (!isTRUE(force_min_version) && isTRUE(cleanup_on_fail)) {
      message("Version requirement failed - removing cloned directory: ", repo_path)
      try(unlink(repo_path, recursive = TRUE, force = TRUE),
          silent = TRUE)
    }

    stop(e$message, call. = FALSE)
  })

  invisible(repo_path)
}

#################################################
## Download param dir/files with 'getParam_gh' ##
#################################################
## Helper function for safe destination handling to avoid unintentional overwrite
#' @keywords internal
.handle_existing_dest <- function(dest, tag, mode = c("backup","side_by_side","ask","abort","overwrite"),
                                  force = TRUE, quiet = FALSE) {
  mode <- match.arg(mode)
  say <- function(...) if (!quiet) message(...)

  is_nonempty_dir <- dir.exists(dest) &&
    length(list.files(dest, all.files = TRUE, no.. = TRUE)) > 0

  ts <- format(Sys.time(), "%Y%m%d-%H%M%S")

  if (!dir.exists(dest)) {
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    return(list(dest = dest, backup = NA_character_, action = "created"))
  }

  if (!is_nonempty_dir) {
    return(list(dest = dest, backup = NA_character_, action = "reuse_empty"))
  }

  if (mode == "abort") {
    stop(sprintf("Destination '%s' exists and is not empty. Aborting (mode='abort').", dest), call. = FALSE)
  }

  if (mode == "side_by_side") {
    new_dest <- sprintf("%s.NEW-%s-%s", dest, tag, ts)
    dir.create(new_dest, recursive = TRUE, showWarnings = FALSE)
    say("Destination exists; installing side-by-side into: ", new_dest)
    return(list(dest = new_dest, backup = NA_character_, action = "side_by_side"))
  }

  if (mode == "ask") {
    if (!interactive()) {
      stop(sprintf(
        "Destination '%s' exists and is not empty, and session is non-interactive.\nRefusing to modify it. Use mode='backup' or mode='side_by_side' explicitly.",
        dest
      ), call. = FALSE)
    }
    cat(sprintf("\nDestination '%s' already exists and is not empty.\n", dest))
    cat("Choose action:\n")
    cat("  [1] Backup existing and install fresh (recommended)\n")
    cat("  [2] Install side-by-side (no changes to existing)\n")
    cat("  [3] Abort\n")
    cat("Selection [1]: ")
    ans <- trimws(readline())
    if (ans == "" || ans == "1") mode <- "backup"
    else if (ans == "2") mode <- "side_by_side"
    else stop("Aborted by user.", call. = FALSE)
  }

  if (mode == "backup") {
    backup <- sprintf("%s.BAK-%s", dest, ts)
    say("Destination exists; moving current directory to: ", backup)
    ok <- file.rename(dest, backup)
    if (!ok) {
      stop(sprintf(
        "Failed to rename '%s' to '%s'. Check permissions or whether files are in use.",
        dest, backup
      ), call. = FALSE)
    }
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    return(list(dest = dest, backup = backup, action = "backup_replace"))
  }

  if (mode == "overwrite") {
    if (!isTRUE(force)) {
      stop(sprintf("Refusing to overwrite '%s' because force=FALSE (mode='overwrite').", dest), call. = FALSE)
    }
    say("Destination exists; overwriting contents (mode='overwrite').")
    existing <- list.files(dest, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    if (length(existing)) unlink(existing, recursive = TRUE, force = TRUE)
    return(list(dest = dest, backup = NA_character_, action = "overwrite"))
  }

  stop("Unhandled mode.", call. = FALSE)
}

#' @keywords internal
read_block_from_manifest <- function(block_name, manifest = "manifest.yml", allow_missing = FALSE) {
  if (!file.exists(manifest)) {
    if (allow_missing) return(list(name=NA_character_, repo=NA_character_, tag=NA_character_))
    stop(sprintf("Cannot find '%s' in current working directory.", manifest), call. = FALSE)
  }
  x <- readLines(manifest, warn = FALSE)

  i0 <- grep(sprintf("^\\s*%s\\s*:\\s*$", block_name), x)
  if (!length(i0)) {
    if (allow_missing) return(list(name=NA_character_, repo=NA_character_, tag=NA_character_))
    stop(sprintf("manifest.yml has no '%s:' block.", block_name), call. = FALSE)
  }
  i0 <- i0[1]

  i1 <- i0
  if (i0 < length(x)) {
    for (i in (i0 + 1):length(x)) {
      if (grepl("^\\S", x[i])) break
      i1 <- i
    }
  }
  block <- if (i1 >= i0 + 1) x[(i0 + 1):i1] else character()

  get_key <- function(key) {
    m <- grep(sprintf("^\\s+%s\\s*:\\s*", key), block)
    if (!length(m)) return(NA_character_)
    val <- sub(sprintf("^\\s+%s\\s*:\\s*", key), "", block[m[1]])
    val <- sub("\\s+#.*$", "", val)
    val <- trimws(val)
    val <- sub('^"(.*)"$', "\\1", val)
    val <- sub("^'(.*)'$", "\\1", val)
    if (!nzchar(val)) return(NA_character_)
    val
  }

  list(name = get_key("name"), repo = get_key("repo"), tag = get_key("tag"))
}

#' @keywords internal
.install_gh_archive_dir <- function(repo,
                                    tag,
                                    source_subdir,
                                    dest,
                                    force = TRUE,
                                    mode = c("backup","side_by_side","ask","abort","overwrite"),
                                    quiet = FALSE) {

  mode <- match.arg(mode)
  say <- function(...) if (!quiet) message(...)

  if (is.null(repo) || !nzchar(repo)) stop("repo is empty.", call. = FALSE)
  if (is.null(tag)  || !nzchar(tag))  stop("tag is empty.", call. = FALSE)

  # Resolve owner/repo and a clean base GitHub URL
  base <- repo
  base <- sub("\\.git$", "", base)
  base <- sub("/+$", "", base)

  if (grepl("^https?://", base)) {
    m <- regexec("^https?://github\\.com/([^/]+/[^/#?]+)", base)
    reg <- regmatches(base, m)[[1]]
    if (length(reg) < 2) stop("Could not parse GitHub repo from URL: ", repo, call. = FALSE)
    orgrepo <- reg[2]
    repo_url <- sprintf("https://github.com/%s", orgrepo)
  } else {
    orgrepo <- base
    repo_url <- sprintf("https://github.com/%s", orgrepo)
  }

  # GitHub tag archive zip URL
  zip_url <- sprintf("%s/archive/refs/tags/%s.zip", repo_url, utils::URLencode(tag, reserved = TRUE))

  say("Installing GitHub archive directory")
  say("  repo: ", orgrepo)
  say("  tag:  ", tag)
  say("  url:  ", zip_url)
  say("  dest: ", dest)
  say("  src:  ", source_subdir)
  say("  mode: ", mode)

  # Destination handling (SAFE by default)
  dest_info <- .handle_existing_dest(dest = dest, tag = tag, mode = mode, force = force, quiet = quiet)
  dest <- dest_info$dest

  # Temp workspace
  tmp_root <- tempfile("gh_pack_")
  dir.create(tmp_root)
  on.exit(unlink(tmp_root, recursive = TRUE, force = TRUE), add = TRUE)

  zip_path <- file.path(tmp_root, "pack.zip")

  # Download (base R). On Windows, mode="wb" is important.
  say("Downloading tag archive...")
  tryCatch(
    utils::download.file(zip_url, destfile = zip_path, mode = "wb", quiet = quiet),
    error = function(e) {
      stop(sprintf("Failed to download tag archive. URL: %s\nError: %s", zip_url, conditionMessage(e)),
           call. = FALSE)
    }
  )

  # Extract
  say("Extracting...")
  utils::unzip(zip_path, exdir = tmp_root)

  # GitHub expands to: <repo>-<tag>/
  top_dirs <- list.dirs(tmp_root, recursive = FALSE, full.names = TRUE)
  top_dirs <- top_dirs[basename(top_dirs) != "__MACOSX"]
  top_dir <- top_dirs[dir.exists(top_dirs)][1]
  if (is.na(top_dir) || !dir.exists(top_dir)) {
    stop("Could not locate extracted top-level directory in archive.", call. = FALSE)
  }

  # Source directory inside archive
  if (is.null(source_subdir) || !nzchar(source_subdir)) {
    stop("source_subdir is empty.", call. = FALSE)
  }
  src_dir <- file.path(top_dir, source_subdir)
  if (!dir.exists(src_dir)) {
    stop(sprintf("source_subdir '%s' does not exist in archive for tag '%s'.", source_subdir, tag),
         call. = FALSE)
  }

  # Copy files into dest
  files <- list.files(src_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  if (!length(files)) stop("Source directory appears empty; nothing to copy.", call. = FALSE)

  say("Installing files into ", dest, " ...")
  ok <- file.copy(from = files, to = dest, recursive = TRUE, overwrite = TRUE, copy.mode = TRUE)
  if (!all(ok)) stop("Some files failed to copy into destination. Check permissions/locked files.", call. = FALSE)

  # Provenance
  writeLines(
    c(
      paste0("repo=", orgrepo),
      paste0("tag=", tag),
      paste0("archive_url=", zip_url),
      paste0("source_subdir=", source_subdir),
      paste0("date=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
      paste0("mode=", mode),
      if (!is.na(dest_info$backup)) paste0("backup_dir=", dest_info$backup) else NULL
    ),
    file.path(dest, "INSTALLED_FROM.txt")
  )

  if (!quiet && !is.na(dest_info$backup)) {
    message("Backup saved as: ", dest_info$backup)
  }

  say("Done. Installed into ", dest, "/")
  invisible(list(dest = normalizePath(dest, winslash = "/", mustWork = TRUE),
                 tag = tag,
                 repo = orgrepo,
                 backup = dest_info$backup,
                 mode = mode))
}

#' Install a shared CWL parameter directory from a GitHub tag archive
#'
#' Downloads a GitHub tag archive (\code{.zip}) for a parameter repository and installs
#' the contents of a subdirectory (default: \code{"param"}) into the workflow repository
#' (default destination: \code{"param"}).
#'
#' The parameter repository URL and tag are typically specified in the workflow repository's
#' \code{manifest.yml} under the \code{param:} block (keys \code{repo} and \code{tag}).
#'
#' @param param_repo Parameter repository identifier. Either \code{"ORG/repo"} or a GitHub URL
#'   like \code{"https://github.com/ORG/repo"}. If left at the default placeholder
#'   (\code{"ORG/cwl-params"}) and \code{manifest.yml} exists, the repo is read from the manifest.
#' @param dest Destination directory in the workflow repository (default \code{"param"}).
#' @param tag Git tag to install. If \code{NULL}, the tag is resolved in priority:
#'   \enumerate{
#'     \item environment variable \code{PARAM_TAG}
#'     \item \code{manifest.yml} \code{param: tag:} (and \code{param: repo:} if \code{param_repo} is left at default)
#'     \item \code{default_tag}
#'   }
#' @param default_tag Fallback tag used if no other tag source is available.
#' @param source_subdir Subdirectory within the parameter repo archive to copy (default \code{"param"}).
#' @param force Logical; only relevant when \code{mode="overwrite"}. If \code{FALSE}, refuse destructive overwrite.
#' @param mode Behavior when \code{dest} exists and is non-empty:
#'   \describe{
#'     \item{\code{"backup"}}{Rename existing \code{dest} to \code{dest.BAK-<timestamp>} and install fresh (default).}
#'     \item{\code{"side_by_side"}}{Install into \code{dest.NEW-<tag>-<timestamp>} without modifying existing \code{dest}.}
#'     \item{\code{"ask"}}{Prompt interactively; refuses in non-interactive sessions.}
#'     \item{\code{"abort"}}{Stop with an error if \code{dest} exists and is non-empty.}
#'     \item{\code{"overwrite"}}{Delete contents of \code{dest} and install fresh (requires \code{force=TRUE}).}
#'   }
#' @param quiet Logical; suppress progress messages.
#'
#' @return Invisibly returns a list with elements:
#' \describe{
#'   \item{dest}{Normalized destination path where parameters were installed.}
#'   \item{tag}{Resolved tag used for installation.}
#'   \item{repo}{Resolved \code{ORG/repo} used for archive download.}
#'   \item{backup}{Backup directory path (if a backup occurred), otherwise \code{NA}.}
#'   \item{mode}{Effective mode used.}
#' }
#' @details
#' \strong{Versioned snapshots via Git tags:}
#' This function installs files from a GitHub \emph{tag archive} URL of the form
#' \code{https://github.com/<ORG>/<REPO>/archive/refs/tags/<TAG>.zip}.
#' To publish an updated parameter (or data) pack, developers must create and push
#' a new Git tag in the corresponding repository (and update the workflow \code{manifest.yml}
#' \code{tag:} field accordingly). Creating a GitHub Release for the tag is optional but
#' recommended for discoverability.
#'
#' Users can override the tag without editing the manifest via the environment variable
#' \code{PARAM_TAG} (or \code{DATA_TAG} for sample data).
#'
#' @examples
#' \dontrun{
#' ## Run from inside a workflow repo, usually generated by genWorkenvir_gh, 
#' ## that has manifest.yml:
#' getParam_gh()
#'
#' ## Explicit repo and tag:
#' getParam_gh(param_repo = "systemPipeR/sprwfcmp-param", tag = "v1.1.1")
#'
#' ## Safer install without touching existing param/:
#' getParam_gh(mode = "side_by_side")
#'
#' ## Developer workflow (run in param/data repo):
#' ## git tag -a v1.2.0 -m "Update CWL params"
#' ## git push --tags
#' ## Then update manifest.yml in each workflow repo to tag: "v1.2.0"
#' }

getParam_gh <- function(param_repo = "ORG/cwl-params",
                        dest = "param",
                        tag = NULL,
                        default_tag = "v1.0.0",
                        source_subdir = "param",
                        force = TRUE,
                        mode = c("backup","side_by_side","ask","abort","overwrite"),
                        quiet = FALSE) {

  mode <- match.arg(mode)

  # Resolve tag/repo: explicit args > env var > manifest.yml > default_tag
  if (is.null(tag) || !nzchar(tag)) {
    env_tag <- Sys.getenv("PARAM_TAG", unset = "")
    if (nzchar(env_tag)) {
      tag <- trimws(env_tag)
    } else if (file.exists("manifest.yml")) {
      mf <- read_block_from_manifest("param", "manifest.yml")
      if (!is.na(mf$tag) && nzchar(mf$tag)) tag <- mf$tag
      if (identical(param_repo, "ORG/cwl-params") && !is.na(mf$repo) && nzchar(mf$repo)) {
        param_repo <- mf$repo
      }
    } else {
      tag <- default_tag
    }
  }
  if (!nzchar(tag)) stop("Resolved param tag is empty.", call. = FALSE)

  .install_gh_archive_dir(
    repo = param_repo,
    tag = tag,
    source_subdir = source_subdir,
    dest = dest,
    force = force,
    mode = mode,
    quiet = quiet
  )
}

#' Install a shared workflow data directory from a GitHub tag archive
#'
#' Downloads a GitHub tag archive (\code{.zip}) for a data repository and installs
#' the contents of the \code{data/} directory into the workflow repository (default
#' destination: \code{"data"}).
#'
#' The data repository URL and tag are typically specified in the workflow repository's
#' \code{manifest.yml} under the \code{data:} block (keys \code{repo} and \code{tag}).
#'
#' @param data_repo Data repository identifier. Either \code{"ORG/repo"} or a GitHub URL
#'   like \code{"https://github.com/ORG/repo"}. If left at the default placeholder
#'   (\code{"ORG/cwl-data"}) and \code{manifest.yml} exists, the repo is read from the manifest.
#' @param dest Destination directory in the workflow repository (default \code{"data"}).
#' @param tag Git tag to install. If \code{NULL}, the tag is resolved in priority:
#'   \enumerate{
#'     \item environment variable \code{DATA_TAG}
#'     \item \code{manifest.yml} \code{data: tag:} (and \code{data: repo:} if \code{data_repo} is left at default)
#'     \item \code{default_tag}
#'   }
#' @param default_tag Fallback tag used if no other tag source is available.
#' @param source_subdir Subdirectory within the data repo archive to copy (default \code{"data"}).
#' @param force Logical; only relevant when \code{mode="overwrite"}. If \code{FALSE}, refuse destructive overwrite.
#' @param mode Behavior when \code{dest} exists and is non-empty: \code{"backup"} (default),
#'   \code{"side_by_side"}, \code{"ask"}, \code{"abort"}, or \code{"overwrite"}.
#' @param quiet Logical; suppress progress messages.
#'
#' @return Invisibly returns a list with elements \code{dest}, \code{tag}, \code{repo}, \code{backup}, and \code{mode}.
#'
#' @details
#' \strong{Versioned snapshots via Git tags:}
#' This function installs files from a GitHub \emph{tag archive} URL of the form
#' \code{https://github.com/<ORG>/<REPO>/archive/refs/tags/<TAG>.zip}.
#' To publish an updated parameter (or data) pack, developers must create and push
#' a new Git tag in the corresponding repository (and update the workflow \code{manifest.yml}
#' \code{tag:} field accordingly). Creating a GitHub Release for the tag is optional but
#' recommended for discoverability.
#'
#' Users can override the tag without editing the manifest via the environment variable
#' \code{PARAM_TAG} (or \code{DATA_TAG} for sample data).
#'
#' @examples
#' \dontrun{
#' ## Run from inside a workflow repo that has manifest.yml:
#' getData_gh()
#'
#' ## Developer workflow (run in param/data repo):
#' ## git tag -a v1.2.0 -m "Update CWL params"
#' ## git push --tags
#' ## Then update manifest.yml in each workflow repo to tag: "v1.2.0"
#' }
getData_gh <- function(data_repo = "ORG/cwl-data",
                       dest = "data",
                       tag = NULL,
                       default_tag = "v1.0.0",
                       source_subdir = "data",
                       force = TRUE,
                       mode = c("backup","side_by_side","ask","abort","overwrite"),
                       quiet = FALSE) {

  mode <- match.arg(mode)
  say <- function(...) if (!quiet) message(...)

  # If user explicitly provided repo/tag, just proceed
  explicit_repo <- !identical(data_repo, "ORG/cwl-data") && nzchar(data_repo)
  explicit_tag  <- !is.null(tag) && nzchar(tag)

  mf <- read_block_from_manifest("data", "manifest.yml", allow_missing = TRUE)

  # Resolve repo from manifest only if placeholder and manifest provides it
  if (!explicit_repo && !is.na(mf$repo) && nzchar(mf$repo)) data_repo <- mf$repo

  # Resolve tag: explicit > env > manifest > default
  if (!explicit_tag) {
    env_tag <- Sys.getenv("DATA_TAG", unset = "")
    if (nzchar(env_tag)) tag <- trimws(env_tag)
    else if (!is.na(mf$tag) && nzchar(mf$tag)) tag <- mf$tag
    else tag <- default_tag
  }

  # If still placeholder repo and manifest had no repo, treat as "no sample data configured"
  if (identical(data_repo, "ORG/cwl-data") && (is.na(mf$repo) || !nzchar(mf$repo)) && !explicit_repo) {
    say("No 'data:' entry in manifest.yml; skipping sample data install.")
    return(invisible(NULL))
  }

  if (!nzchar(tag)) stop("Resolved data tag is empty.", call. = FALSE)

  .install_gh_archive_dir(
    repo = data_repo,
    tag = tag,
    source_subdir = source_subdir,
    dest = dest,
    force = force,
    mode = mode,
    quiet = quiet
  )
}

###############################################
## Internal helper functions for GitHub tags ##
###############################################

#' @keywords internal
available_gh_tags <- function(repo,
                              per_page = 100,
                              max_pages = 3,
                              quiet = TRUE) {

  base <- repo
  base <- sub("\\.git$", "", base)
  base <- sub("/+$", "", base)

  if (grepl("^https?://", base)) {
    m <- regexec("^https?://github\\.com/([^/]+/[^/#?]+)", base)
    reg <- regmatches(base, m)[[1]]
    if (length(reg) < 2) stop("Could not parse GitHub repo from URL: ", repo, call. = FALSE)
    orgrepo <- reg[2]
  } else {
    orgrepo <- base
  }

  fetch_page <- function(page) {
    api <- sprintf("https://api.github.com/repos/%s/tags?per_page=%d&page=%d",
                   orgrepo, per_page, page)
    tmp <- tempfile("tags_")
    on.exit(unlink(tmp, force = TRUE), add = TRUE)

    ok <- tryCatch({
      utils::download.file(api, tmp, mode = "wb", quiet = quiet)
      TRUE
    }, error = function(e) FALSE)

    if (!ok) stop("Failed to query GitHub tags via API: ", api, call. = FALSE)

    txt <- paste(readLines(tmp, warn = FALSE), collapse = "\n")
    m <- gregexpr('"name"\\s*:\\s*"[^"]+"', txt, perl = TRUE)
    hits <- regmatches(txt, m)[[1]]
    if (!length(hits)) return(character())
    tags <- sub('^"name"\\s*:\\s*"([^"]+)"$', "\\1", hits, perl = TRUE)
    unique(tags[nzchar(tags)])
  }

  tags_all <- character()
  for (p in seq_len(max_pages)) {
    tags <- fetch_page(p)
    if (!length(tags)) break
    tags_all <- c(tags_all, tags)
    if (length(tags) < per_page) break
  }

  unique(tags_all)
}

#' @keywords internal
choose_gh_tag <- function(repo, default = 1) {
  tags <- available_gh_tags(repo, quiet = TRUE)
  if (!length(tags)) stop("No tags found.", call. = FALSE)

  cat("\nAvailable tags for ", repo, ":\n", sep = "")
  for (i in seq_along(tags)) cat(sprintf("  [%d] %s\n", i, tags[i]))
  cat("\nSelect a tag number [default ", default, "]: ", sep = "")

  ans <- readline()
  if (!nzchar(ans)) return(tags[default])

  i <- suppressWarnings(as.integer(ans))
  if (is.na(i) || i < 1 || i > length(tags)) stop("Invalid selection.", call. = FALSE)
  tags[i]
}

#' @keywords internal
.read_yaml_block_lines <- function(key, x) {
  # returns the indented lines belonging to a top-level key block:  key:
  i0 <- grep(sprintf("^\\s*%s\\s*:\\s*$", key), x)
  if (!length(i0)) return(character())
  i0 <- i0[1]

  i1 <- i0
  if (i0 < length(x)) {
    for (i in (i0 + 1):length(x)) {
      if (grepl("^\\S", x[i])) break   # next top-level key
      i1 <- i
    }
  }
  if (i1 >= i0 + 1) x[(i0 + 1):i1] else character()
}

#' @keywords internal
.read_bioc_min_from_manifest <- function(manifest = "manifest.yml",
                                        pkgs = c("systemPipeR", "systemPipeRdata")) {
  if (!file.exists(manifest)) return(setNames(rep(NA_character_, length(pkgs)), pkgs))

  x <- readLines(manifest, warn = FALSE)
  out <- setNames(rep(NA_character_, length(pkgs)), pkgs)

  for (pkg in pkgs) {
    block <- .read_yaml_block_lines(pkg, x)
    if (!length(block)) next
    m <- grep("^\\s+bioc_min\\s*:\\s*", block)
    if (!length(m)) next

    val <- sub("^\\s+bioc_min\\s*:\\s*", "", block[m[1]])
    val <- sub("\\s+#.*$", "", val)
    val <- trimws(val)
    val <- sub('^"(.*)"$', "\\1", val)
    val <- sub("^'(.*)'$", "\\1", val)
    if (nzchar(val)) out[pkg] <- val
  }

  out
}

#' @keywords internal
.check_min_pkg_versions <- function(min_versions,
                                    context = "manifest.yml",
                                    force = FALSE,
                                    quiet = FALSE) {
  say <- function(...) if (!quiet) message(...)

  pkgs <- names(min_versions)
  for (pkg in pkgs) {
    req <- min_versions[[pkg]]
    if (is.na(req) || !nzchar(req)) next

    # Missing package is always a hard stop
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf(
        "Package '%s' is required (>= %s) per %s, but is not installed.\nInstall/update via BiocManager::install('%s').",
        pkg, req, context, pkg
      ), call. = FALSE)
    }

    cur <- utils::packageVersion(pkg)
    reqv <- as.package_version(req)

    if (cur < reqv) {
      msg <- sprintf(
        "Package '%s' version %s is installed, but version >= %s is required per %s.\nPlease update via BiocManager::install('%s'). Alternatively, proceed at your own risk by setting force_min_version=TRUE",
        pkg, as.character(cur), req, context, pkg
      )
      if (!isTRUE(force)) stop(msg, call. = FALSE)

      say("WARNING: ", msg)
      say("Proceeding anyway because force_min_version=TRUE (at your own risk).")
    }
  }

  invisible(TRUE)
}



