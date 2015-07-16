pathlist <- function() {
    list(
        targets=system.file("extdata/param", "targets.txt", package="systemPipeRdata", mustWork=TRUE),
        targetsPE=system.file("extdata/param", "targetsPE.txt", package="systemPipeRdata", mustWork=TRUE),
        annotationdir=system.file("extdata/annotation", "", package="systemPipeRdata", mustWork=TRUE),
        fastqdir=system.file("extdata/fastq", "", package="systemPipeRdata", mustWork=TRUE),
        paramdir=system.file("extdata/param", "", package="systemPipeRdata", mustWork=TRUE),
        workflows=system.file("extdata/workflows", "", package="systemPipeRdata", mustWork=TRUE),
        chipseq=system.file("extdata/workflows/chipseq", "", package="systemPipeRdata", mustWork=TRUE),
        rnaseq=system.file("extdata/workflows/rnaseq", "", package="systemPipeRdata", mustWork=TRUE),
        varseq=system.file("extdata/workflows/varseq", "", package="systemPipeRdata", mustWork=TRUE)
    )
}

