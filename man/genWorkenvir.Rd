\name{genWorkenvir}
\alias{genWorkenvir}
\title{
Generate workflow templates
}
\description{
Generates workflow templates for \code{systemPipeR} package. The template environemts
contain a predifined directory structure along with run parameter files and sample data.
The strucutre of the workflow templates and the sample data are described in all
details in the Overview Vignette of the \href{http://bioconductor.org/packages/systemPipeR}{systemPipeR package}.
}
\usage{
genWorkenvir(workflow, mydirname=NULL, bam=FALSE)
}
\arguments{
  \item{workflow}{
  \code{character string} of workflow templates to be generated. Supported values
  can be checked with the \verb{\link{availableWF}()} function.
}
  \item{mydirname}{
  Specifies the name of the workflow directory. The default \code{NULL} uses the
  name of the chosen workflow. An error is issued if a directory of the same name
  and path exists already.
}
  \item{bam}{
  If \code{bam=TRUE} pregenerated short read alignment (BAM) files will be included
  in the \code{results} directory of the workflow environment. Note, these BAM files
  have been generated with the HISAT2 aligner using the FASTQ files provided in the
  \code{data} directory. The default \code{bam=FALSE} omits this step meaning no
  BAM files will be copied into the \code{results} directory.
}
}
\value{
Workflow directory containing sample data and parameter files along with the following subdirectories:
    \item{param/}{stores parameter files}
    \item{data/}{stores input data}
    \item{results/}{stores output results}

For more details, please consult the Overview Vignette (HTML) of the systemPipeR package (http://bioconductor.org/packages/systemPipeR).
}
\details{
Check the output of \verb{\link{availableWF}()} to the current workflow templates
available on \code{systemPipeR} Organization.

For an \code{interactive()} session, the \code{readline()} function provides the option
choose between proceeding or not, through options: \code{yes} or \code{no}.
For non-interactive use, if there is no package install, the option \code{yes} will be selected.
}
\author{
Thomas Girke and Daniela Cassol
}
\examples{
## Return location of sample data
samplepaths <- pathList()
\dontrun{
## Generate varseq workflow environment
genWorkenvir(workflow="varseq", mydirname=NULL, url=NULL, urlname=NULL)
setwd("varseq")
}
}
\keyword{ utilities }
