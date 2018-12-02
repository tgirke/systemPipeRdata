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
        varseq=system.file("extdata/workflows/varseq", "", package="systemPipeRdata", mustWork=TRUE)
    )
}

#########################################
## Generate environments for workflows ##
#########################################
genWorkenvir <- function(workflow, mydirname=NULL, bam=FALSE) {
    ## Input validity check
    check_workflow <- c("varseq", "rnaseq", "riboseq", "chipseq")
    if(!workflow %in% check_workflow) stop(paste("workflow can only be assigned one of:", paste(check_workflow, collapse=", ")))
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
        file.copy(normalizePath(list.files(pathList()$fastqdir, "*", full.names=TRUE)), "varseq/data", overwrite=TRUE, recursive=TRUE)
        if(bam==TRUE) file.copy(normalizePath(list.files(pathList()$bamdir, "*", full.names=TRUE)), "varseq/results", overwrite=TRUE, recursive=TRUE)
        file.copy(normalizePath(list.files(pathList()$annotationdir, "*", full.names=TRUE)), "varseq/data", overwrite=TRUE, recursive=TRUE)
        file.copy(pathList()$paramdir, "varseq/", recursive=TRUE)
        file.copy(c("varseq/param/batchtools.slurm.tmpl", "varseq/param/.batchtools.conf.R"), "./varseq")
        file.copy(c("varseq/param/gatk_run.sh", "varseq/param/sambcf_run.sh"), "./varseq")
        file.copy(c("varseq/param/targetsPE.txt", "varseq/param/targets.txt"), "./varseq")
        file.copy(c("varseq/param/Makefile_varseq"), "./varseq/Makefile")
        file.copy(c("varseq/param/bibtex.bib"), "./varseq/bibtex.bib")
    } else if(workflow=="rnaseq") {
        file.copy(pathList()$rnaseq, mydirname2temp, recursive=TRUE)
        file.rename(paste0(normalizePath(mydirname2temp), "/", workflow), mydirname2) # generates final dir
        unlink(mydirname2temp, recursive=TRUE) # removes temp dir
        file.copy(normalizePath(list.files(pathList()$fastqdir, "*", full.names=TRUE)), "rnaseq/data", overwrite=TRUE, recursive=TRUE)
        if(bam==TRUE) file.copy(normalizePath(list.files(pathList()$bamdir, "*", full.names=TRUE)), "rnaseq/results", overwrite=TRUE, recursive=TRUE)
        file.copy(normalizePath(list.files(pathList()$annotationdir, "*", full.names=TRUE)), "rnaseq/data", overwrite=TRUE, recursive=TRUE)
        file.copy(pathList()$paramdir, "rnaseq/", recursive=TRUE)
        file.copy(c("rnaseq/param/batchtools.slurm.tmpl", "rnaseq/param/.batchtools.conf.R"), "./rnaseq")
        file.copy(c("rnaseq/param/targetsPE.txt", "rnaseq/param/targets.txt"), "./rnaseq")
        file.copy(c("rnaseq/param/Makefile_rnaseq"), "./rnaseq/Makefile")
        file.copy(c("rnaseq/param/bibtex.bib"), "./rnaseq/bibtex.bib")
    } else if(workflow=="riboseq") {
        file.copy(pathList()$riboseq, mydirname2temp, recursive=TRUE)
        file.rename(paste0(normalizePath(mydirname2temp), "/", workflow), mydirname2) # generates final dir
        unlink(mydirname2temp, recursive=TRUE) # removes temp dir
        file.copy(normalizePath(list.files(pathList()$fastqdir, "*", full.names=TRUE)), "riboseq/data", overwrite=TRUE, recursive=TRUE)
        if(bam==TRUE) file.copy(normalizePath(list.files(pathList()$bamdir, "*", full.names=TRUE)), "riboseq/results", overwrite=TRUE, recursive=TRUE)
        file.copy(normalizePath(list.files(pathList()$annotationdir, "*", full.names=TRUE)), "riboseq/data", overwrite=TRUE, recursive=TRUE)
        file.copy(pathList()$paramdir, "riboseq/", recursive=TRUE)
        file.copy(c("riboseq/param/batchtools.slurm.tmpl", "riboseq/param/.batchtools.conf.R"), "./riboseq")
        file.copy(c("riboseq/param/targetsPE.txt", "riboseq/param/targets.txt"), "./riboseq")
        file.copy(c("riboseq/param/Makefile_riboseq"), "./riboseq/Makefile")
        file.copy(c("riboseq/param/bibtex.bib"), "./riboseq/bibtex.bib")
    } else if(workflow=="chipseq") {
        file.copy(pathList()$chipseq, mydirname2temp, recursive=TRUE)
        file.rename(paste0(normalizePath(mydirname2temp), "/", workflow), mydirname2) # generates final dir
        unlink(mydirname2temp, recursive=TRUE) # removes temp dir
        file.copy(normalizePath(list.files(pathList()$fastqdir, "*", full.names=TRUE)), "chipseq/data", overwrite=TRUE, recursive=TRUE)
        if(bam==TRUE) file.copy(normalizePath(list.files(pathList()$bamdir, "*", full.names=TRUE)), "chipseq/results", overwrite=TRUE, recursive=TRUE)
        file.copy(normalizePath(list.files(pathList()$annotationdir, "*", full.names=TRUE)), "chipseq/data", overwrite=TRUE, recursive=TRUE)
        file.copy(pathList()$paramdir, "chipseq/", recursive=TRUE)
        file.copy(c("chipseq/param/batchtools.slurm.tmpl", "chipseq/param/.batchtools.conf.R"), "./chipseq")
        file.copy(c("chipseq/param/targetsPE_chip.txt", "chipseq/param/targets_chip.txt"), "./chipseq")
        file.copy(c("chipseq/param/Makefile_chipseq"), "./chipseq/Makefile")
        file.copy(c("chipseq/param/bibtex.bib"), "./chipseq/bibtex.bib")
    } else {
        stop(paste("workflow can only be assigned one of:", paste(check_workflow, collapse=", ")))
    }
    print(paste("Generated", mydirname2, "directory. Next run in", workflow, "directory, the R code from *.Rmd (*.Rnw) template interactively. Alternatively, workflows can be exectued with a single command as instructed in the vignette."))
}
## Usage:
# genWorkenvir(workflow="varseq", mydirname=NULL)
