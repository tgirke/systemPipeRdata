\name{genWorkenvir}
\alias{genWorkenvir}
\title{
Generate workflow environments
}
\description{
Generates recommended directory structure for sample workflows populated
with sample data and parameter files.
}
\usage{
genWorkenvir(workflow)
}
\arguments{
  \item{workflow}{
Can be assigned workflow type as character string. Currently, supported are:
\code{rnaseq}, \code{varseq} or \code{chipseq} 
}
}
\value{
Directory structure with sample data and parameter files.
}
\author{
Thomas Girke
}
\examples{
\dontrun{
## Generate varseq workflow environment
genWorkenvir(workflow="varseq")
}
}
\keyword{ utilities }
