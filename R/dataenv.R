#################################
## Obtain paths to sample data ##
#################################
pathList <- function() {
  list(
    targets=system.file("extdata/param", "targets.txt", package="systemPipeRdata", mustWork=TRUE),
    targetsPE=system.file("extdata/param", "targetsPE.txt", package="systemPipeRdata", mustWork=TRUE),
    annotationdir=system.file("extdata/annotation", "", package="systemPipeRdata", mustWork=TRUE),
    fastqdir=system.file("extdata/fastq", "", package="systemPipeRdata", mustWork=TRUE),
    bamdir=system.file("extdata/bam", "", package="systemPipeRdata", mustWork=TRUE),
    paramdir=system.file("extdata/param", "", package="systemPipeRdata", mustWork=TRUE),
    workflows=system.file("extdata/workflows", "", package="systemPipeRdata", mustWork=TRUE),
    chipseq=system.file("extdata/workflows/chipseq", "", package="systemPipeRdata", mustWork=TRUE),
    rnaseq=system.file("extdata/workflows/rnaseq", "", package="systemPipeRdata", mustWork=TRUE),
    riboseq=system.file("extdata/workflows/riboseq", "", package="systemPipeRdata", mustWork=TRUE),
    varseq=system.file("extdata/workflows/varseq", "", package="systemPipeRdata", mustWork=TRUE),
    new=system.file("extdata/workflows/new", "", package="systemPipeRdata", mustWork=TRUE)
  )
}

#########################################
## Generate environments for workflows ##
#########################################
genWorkenvir <- function(workflow, mydirname=NULL, bam=FALSE, ref="master", subdir=NULL, url=NULL, urlname=NULL) {
  ## Input validity check
  if(grepl("/", workflow)==TRUE) {
      ## if the package is defined it
      mydirname2 <- .genWorkenvirPKG(package_repo=workflow, ref=ref, subdir=subdir, mydirname=mydirname, build_env=TRUE, silent=TRUE)
      genWorkdata(path = mydirname2)
  } else {
    check_workflow <- dir(system.file("extdata/workflows", "", package="systemPipeRdata", mustWork=TRUE))
    if(!workflow %in% check_workflow) stop(paste("workflow can only be assigned one of:", paste(check_workflow, collapse=", ")))
    if(workflow=="new" && is.null(mydirname)) warning("It is recommended to specify the workflow directory name, using 'mydirname' argument")
    if(all(!c(is.null(mydirname), is.character(mydirname)))) stop("mydirname can only be assigned 'NULL' or a character vector of length 1")
    ## If 'mydirname' is NULL (default) use value assigned to 'workflow' as directory name
    if(is.null(mydirname)) {
      mydirname2 <- workflow 
    } else {
      mydirname2 <- mydirname
    }
    ## Generate temp workflow directory
    mydirname2temp <- paste0(mydirname2, "_temp", paste(sample(0:9, 4), collapse=""))
    if(dir.exists(mydirname2) | dir.exists(mydirname2temp)) stop("Directory name assigned to 'mydirname' or its temp variant exists already. Please assign to 'mydirname' a different name or rename/delete existing one.")
    dir.create(mydirname2temp)
    ## Move workflow templates into workflow directory
    if(workflow=="varseq") {
      file.copy(pathList()$varseq, mydirname2temp, recursive=TRUE)
      file.rename(paste0(normalizePath(mydirname2temp), "/", workflow), mydirname2) # generates final dir
      unlink(mydirname2temp, recursive=TRUE) # removes temp dir
      # if we desire to rename *Rmd files... not sure if it is important...
      # if(!is.null(mydirname)){
      #   rename.files <- list.files(mydirname2, pattern = "systemPipeR*", full.names = T)
      #   sapply(rename.files, function(x) { file.rename(from=x, to=sub('.*\\.', paste0(mydirname2, "/systemPipeR_", mydirname2, "."), x))})
      # }
      file.copy(c(paste0(pathList()$paramdir, "gatk_run.sh"), paste0(pathList()$paramdir, "sambcf_run.sh")), paste0(mydirname2, "/"))
      file.copy(c(paste0(pathList()$paramdir, "Makefile_varseq")), paste0(mydirname2, "/Makefile"))
    } else if(workflow=="rnaseq") {
      file.copy(pathList()$rnaseq, mydirname2temp, recursive=TRUE)
      file.rename(paste0(normalizePath(mydirname2temp), "/", workflow), mydirname2) # generates final dir
      unlink(mydirname2temp, recursive=TRUE) # removes temp dir
      file.copy(c(paste0(pathList()$paramdir, "Makefile_rnaseq")), paste0(mydirname2, "/Makefile"))
    } else if(workflow=="riboseq") {
      file.copy(pathList()$riboseq, mydirname2temp, recursive=TRUE)
      file.rename(paste0(normalizePath(mydirname2temp), "/", workflow), mydirname2) # generates final dir
      unlink(mydirname2temp, recursive=TRUE) # removes temp dir
      file.copy(c(paste0(pathList()$paramdir, "Makefile_riboseq")), paste0(mydirname2, "/Makefile"))
      # file.copy(c(paste0(pathList()$paramdir, "bibtex.bib")), paste0(mydirname2, "/bibtex.bib"), overwrite=TRUE)
    } else if(workflow=="chipseq") {
      file.copy(pathList()$chipseq, mydirname2temp, recursive=TRUE)
      file.rename(paste0(normalizePath(mydirname2temp), "/", workflow), mydirname2) # generates final dir
      unlink(mydirname2temp, recursive=TRUE) # removes temp dir
      file.copy(c(paste0(pathList()$paramdir, "targetsPE_chip.txt"), paste0(pathList()$paramdir, "targets_chip.txt")), paste0(mydirname2, "/"))
      file.copy(c(paste0(pathList()$paramdir, "Makefile_chipseq")), paste0(mydirname2, "/Makefile"))
      # file.copy(c(paste0(pathList()$paramdir, "bibtex.bib")), paste0(mydirname2, "/bibtex.bib"), overwrite=TRUE)
    } else if(workflow=="new"){
      file.copy(pathList()$new, mydirname2temp, recursive=TRUE)
      file.rename(paste0(normalizePath(mydirname2temp), "/", workflow), mydirname2) # generates final dir
      unlink(mydirname2temp, recursive=TRUE) # removes temp dir
    } else {
      stop(paste("workflow can only be assigned one of:", paste(check_workflow, collapse=", ")))
    }
    ## Moving data and param common files
    file.copy(normalizePath(list.files(pathList()$fastqdir, "*", full.names=TRUE)), paste0(mydirname2, "/data"), overwrite=TRUE, recursive=TRUE)
    file.copy(normalizePath(list.files(pathList()$annotationdir, "*", full.names=TRUE)), paste0(mydirname2, "/data"), overwrite=TRUE, recursive=TRUE)
    file.copy(c(normalizePath(paste0(pathList()$paramdir, "/targetsPE.txt")), normalizePath(paste0(pathList()$paramdir, "/targets.txt"))), paste0(mydirname2, "/"))
    file.copy(c(paste0(pathList()$paramdir, "bibtex.bib")), paste0(mydirname2, "/bibtex.bib"), overwrite=TRUE)
    if(bam==TRUE) file.copy(normalizePath(list.files(pathList()$bamdir, "*", full.names=TRUE)), paste0(mydirname2, "/results"), overwrite=TRUE, recursive=TRUE)
    file.copy(pathList()$paramdir, paste0(mydirname2, "/"), recursive=TRUE)
    file.copy(c(paste0(pathList()$paramdir, "batchtools.slurm.tmpl"), paste0(pathList()$paramdir, ".batchtools.conf.R")), paste0(mydirname2, "/"))
  }
  if(all(c(!is.null(url), is.null(urlname)))) stop("argument 'urlname' missing. The argument can only be assigned 'NULL' when no url is provided. The argument should be assigned as a character vector of length 1")
  if(!is.null(url)){
    download.file(url=url, destfile = paste0("./", mydirname2, "/", urlname), quiet = TRUE)
  }
  print(paste("Generated", mydirname2, "directory. Next run in", workflow, "directory, the R code from *.Rmd template interactively. Alternatively, workflows can be exectued with a single command as instructed in the vignette."))
}
## Usage:
# genWorkenvir(workflow="varseq", mydirname=NULL)
# genWorkenvir(workflow="systemPipeR/systemPipeVARseq", mydirname=NULL)
# genWorkenvir(workflow="systemPipeR/systemPipeRNAseq", mydirname="rnaseq",
# url = "https://raw.githubusercontent.com/systemPipeR/systemPipeRNAseq/cluster/vignettes/systemPipeRNAseq.Rmd",
# urlname = "rnaseq_V-cluster.Rmd")

############################################
## Generate data/param for SPR_packages ##
############################################
genWorkdata <- function(path=getwd(), data=TRUE, param=TRUE){
  ## TODO: we can add a check if it is a sysproject
  suppressWarnings({
    data.path <- normalizePath(paste0(path, "/data/"))
    param.path <- normalizePath(paste0(path, "/param/"))
  })
  if(data==TRUE){
    if (dir.exists(data.path) == FALSE) (dir.create(data.path))
    file.copy(normalizePath(list.files(pathList()$fastqdir, "*", full.names=TRUE)), data.path, overwrite=TRUE, recursive=TRUE)
    file.copy(normalizePath(list.files(pathList()$annotationdir, "*", full.names=TRUE)), data.path, overwrite=TRUE, recursive=TRUE)
    print("The 'demo data'was successfully copied to your project.")
  }
  if (param==TRUE){
    if (dir.exists(param.path) == FALSE) (dir.create(param.path))
    file.copy(pathList()$paramdir, normalizePath(path), recursive=TRUE)
    print("The 'param' directory was successfully copied to your project.")
  }
}

## Usage:
# genWorkdata()

############################################
## Generate environments for SPR_packages ##
############################################
# This function installs packages from GitHub. Internally, it is used remotes::install_github(), which also requires remotes to be installed. 
## Arguments
# package_repo Repository Github address in the format `username/repo`. 
# ref Desired GitHub reference for the `"branch"` name. Default to `"master"`.
# subdir subdirectory within GitHub repo that contains the R package.
# mydirname Specifies the name of the workflow directory and the *Rmd file. Default is NULL, using the workflow/package name.
# build_env Internally the function uses rmarkdown::draft function to populate the directory with all the files,  the same way if it is used the Rstudio create a new file. Default is TRUE.
# silent If TRUE, suppress output.
# Details
# For an interactive() session, the readline() function is used to choose between `yes` or `no`.
# For non-interactive use, if there is no package install, the option YES will be chosen.

.genWorkenvirPKG <- function(package_repo, ref="master", subdir=NULL, mydirname=NULL, build_env=TRUE, silent=FALSE){
  ## Check if pkg is installed
  pkgs <- c("remotes", sub('.*/', '', package_repo))
  pkg_name <- sub('.*/', '', package_repo)
  pkg_req <- sapply(pkgs, function(x) !requireNamespace(x, quietly=TRUE))
  ## If 'mydirname' is NULL (default) use value assigned to 'pkg_name' as directory/*Rmd name
  if(is.null(mydirname)) {
    mydirname2 <- pkg_name 
  } else {
    mydirname2 <- mydirname
  }
  ## Question
  if(any(pkg_req)) {
    ## For an interactive() session
    if(interactive()){
      pkg_install <- readline(cat("There is no package called", paste0(pkgs, collapse = " OR "), "\n",
                                  "Would you like to install this package now? Type a number: \n 1. Yes \n 2. No \n"))
    } else {
      ## For an non-interactive session
      pkg_install <-"1"
    }
    print(paste("Installing package", pkg_name))
    for(i in which(pkg_req)){
      if (pkg_install == "1"){
        if (!requireNamespace("remotes", quietly=TRUE))
          install.packages("remotes", quiet = TRUE)
        remotes::install_github(repo=package_repo, ref=ref, subdir=subdir, quiet = TRUE, upgrade = "always", build_vignettes=TRUE, dependencies=TRUE)
      } else if (pkg_install == 2)
        stop (paste("The", pkg_name, "package is required, you can install the package and come back later!"))
    }
  }
  print(paste("Package", pkg_name, "is successfully installed!"))
  if(build_env==TRUE){
    rmarkdown::draft(paste0(mydirname2, ".Rmd"), template = pkg_name, package = pkg_name, edit=FALSE)
    if(silent != TRUE) print(paste("Generated", mydirname2, "directory. Next run in", pkg_name, "directory, the R code from *.Rmd template interactively. Alternatively, workflows can be exectued with a single command as instructed in the vignette."))
  } else {
    dump <- "do nothing"
  }
  return(mydirname2)
}

## Usage:
# .genWorkenvirPKG(package_repo="systemPipeR/systemPipeRIBOseq")

############################################
## Check Workflow Templates Availability ##
############################################
# This function checks the workflow templates availability from systemPipeRdata package and 
# also from GitHub. 
## Arguments
# github if TRUE, it will return the available Workflow Templates packages from https://github.com/systemPipeR
## Note:
# We are assuming that workflow templates repositories under [systemPipeR Organization](https://github.com/systemPipeR/)
# content the keywords in the description "Workflow Template" and in "Topics" we expected "systempiper" and "release" or "development" words.

availableWF <- function(github = FALSE){
  check_workflow <- dir(system.file("extdata/workflows", "", package="systemPipeRdata", mustWork=TRUE))
  if(github==TRUE){
    ## get all the repos
    tryCatch(get_repos <- jsonlite::fromJSON("https://api.github.com/orgs/systemPipeR/repos"),
             error = function(e){
               stop('You have reached the limit for GitHub API requests. The limit is 60 requests per hour.')
             } 
    )  
    WFrepos <- data.frame(workflow=get_repos$name, html=get_repos$html_url, description=get_repos$description)
    ## Keep repos with the right description
    WFrepos <- WFrepos[WFrepos$description %in% "Workflow Template", ]
    ## branches
    for(i in seq_along(WFrepos$workflow)){
      branches <- jsonlite::fromJSON(paste0("https://api.github.com/repos/systemPipeR/", WFrepos$workflow[i], "/branches"))$name
      WFrepos$branches[i] <- list(branches)
    }
    ## Searching for topics
    release <- data.frame(workflow=jsonlite::fromJSON("https://api.github.com/search/repositories?q=topic:systempiper+topic:release")$items$name)
    devel <- data.frame(workflow=jsonlite::fromJSON("https://api.github.com/search/repositories?q=topic:systempiper+topic:development")$items$name)
    ## Merge topics and version
    WFrepos_devel <- cbind(merge(WFrepos, devel, by="workflow"), version=rep("devel"))
    WFrepos_release <- cbind(merge(WFrepos, release, by="workflow"), version=rep("release"))
    ## final
    WFrepos <- rbind(WFrepos_release, WFrepos_devel)
    WFrepos$workflow <- paste0("systemPipeR/", WFrepos$workflow)
    WFrepos <- WFrepos[c("workflow", "branches", "version", "html", "description")]
    list <- list(systemPipeRdata=check_workflow, github=WFrepos)
    return(list)
  }
  return(list(systemPipeRdata=check_workflow))
}

## Usage:
# availableWF()
# availableWF(github = TRUE)


