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
genWorkenvir(workflow, mydirname=NULL, bam=FALSE, ref="master", subdir=NULL, url=NULL, urlname=NULL)
}
\arguments{
  \item{workflow}{
  \code{character string} of workflow templates to be generated. Supported values 
  can be checked with the \verb{\link{availableWF}()} function. Workflow name 
  containing a '/' are treated as GitHub repositories and installed using 
  \verb{BiocManager::\link{install_github}()}. For Github workflow templates, the format
   \code{systemPipeR/repo} is required. See Details.
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
  \item{ref}{ 
  Desired GitHub reference for the \code{branch} name. Default to \code{master} branch.
}
  \item{subdir}{
  subdirectory within GitHub repo that contains the R package, if it is required.
}
  \item{url}{
  \code{character string} of a specifies the \code{URL} of a different version of 
  the \code{R Markdown} workflow template or other file to download. 
  The default \code{NULL} copies the current version available 
  in the `systemPipeRdata` or the workflow template on GitHub being selected. 
}
  \item{urlname}{
  \code{character string} with the name where the downloaded file is saved.
  This is argument is required when the \code{url} is provided. 
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
When installing GitHub Workflow Packages under \href{https://github.com/systemPipeR}{systemPipeR Organization}, 
it is requested to the remotes package function \verb{BiocManager::\link{install_github}()}
to build vignettes and also install all the dependencies, via \code{dependencies=TRUE}, \code{build_vignettes=TRUE}.

Check the output of \verb{\link{availableWF}()} to the current workflow templates
available on \code{systemPipeR} Organization. The argument \code{}

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

## List Workflow Templates from \code{systemPipeRdata} package and \code{systemPipeR} Organization 
availableWF(github = TRUE)
## Generate 'systemPipeR/systemPipeChIPseq' workflow environment 
genWorkenvir(workflow="systemPipeR/systemPipeChIPseq", mydirname=NULL, ref="master", subdir=NULL)
setwd("systemPipeChIPseq")

## Download a specific R Markdown file
genWorkenvir(workflow="systemPipeR/systemPipeRNAseq", mydirname="rnaseq", url = "https://raw.githubusercontent.com/systemPipeR/systemPipeRNAseq/cluster/vignettes/systemPipeRNAseq.Rmd", urlname = "rnaseq_V-cluster.Rmd")
}
}
\keyword{ utilities }
