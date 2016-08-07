\name{genWorkenvir}
\alias{genWorkenvir}
\title{
Generate workflow templates
}
\description{
Generates workflow templates for \code{systemPipeR} package. The template environemts contain a predifined directory structure along with run parameter files and sample data. The strucutre of the workflow templates and the sample data are described in all details in the Overview Vignette of the systemPipeR package here: 
http://bioconductor.org/packages/systemPipeR.
}
\usage{
genWorkenvir(workflow, mydirname=NULL, bam=FALSE)
}
\arguments{
  \item{workflow}{
Can be assigned workflow type as character string. Supported values are:
\code{rnaseq}, \code{riboseq}, \code{varseq} or \code{chipseq} 
}
  \item{mydirname}{
Specifies the name of the workflow directory. The default \code{NULL} uses the name of the chosen workflow. An error is issued if a directory of the same name and path exists already.
}
  \item{bam}{
If \code{bam=TRUE} pregenerated short read alignment (BAM) files will be included in the \code{results} directory of the workflow environment. Note, these BAM files have been generated with the HISAT2 aligner using the FASTQ files provided in the \code{data} directory. 
The default \code{bam=FALSE} omits this step meaning no BAM files will be copied into the \code{results} directory.
}
}
\value{
Workflow directory containing sample data and parameter files along with the following subdirectories:
    \item{param/}{stores parameter files}
    \item{data/}{stores input data}
    \item{results/}{stores output results}

For more details, please consult the Overview Vignette (HTML) of the systemPipeR package (http://bioconductor.org/packages/systemPipeR).
}
\author{
Thomas Girke
}
\examples{
## Return location of sample data
samplepaths <- pathList()
\dontrun{
## Generate varseq workflow environment
genWorkenvir(workflow="varseq", mydirname=NULL)
setwd("varseq")
}
}
\keyword{ utilities }
