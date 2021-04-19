\name{getSubsetReads}
\alias{getSubsetReads}
\title{Subsetting fastq data}
\description{
Returns subsets of fastq files data based on specific mapping regions or
list of genes or GRanges object.
}
\usage{
getSubsetReads(args,
  geneList = NULL,
  gr = NULL,
  MappingRegion = 1:1e+05,
  sample_range = 90000:1e+05,
  truncate_refs = TRUE,
  id_read_number = TRUE,
  annotation = "data/tair10.gff",
  reference = "data/tair10.fasta",
  annot_outname = "tair10_sub.gff",
  ref_outname = "tair10_sub.fasta",
  outdir = "data/subset/",
  silent = FALSE
)
}
\arguments{
\item{args}{object of class \code{SYSargs2}.}

\item{geneList}{selected genes list to retrieve the reads from the fastq file.}

\item{gr}{an object containing genomic ranges to retrieve the reads from the fastq file.}

\item{MappingRegion}{integers ranges of start and end of chromosome position to retrieve the reads from the fastq file.}

\item{sample_range}{random range to subsetted the fastq file.}

\item{truncate_refs}{logical. If TRUE it will generate reference genome and annotation subset file.}

\item{id_read_number}{if fastq file contains sequence name with read number (\verb{$ri} - \verb{--defline-seq '@$sn[_$rn]/$ri'}).}

\item{annotation}{path to annotation file.}

\item{reference}{path to reference genome.}

\item{annot_outname}{character name of the annotation output file.}

\item{ref_outname}{character name of the reference genome output file.}

\item{outdir}{path to output directory.}

\item{silent}{if set to TRUE, all messages returned by the function will be suppressed.}
}
\value{
Workflow directory containing sample data and parameter files along with the following subdirectories:
    \item{param/}{stores parameter files}
    \item{data/}{stores input data}
    \item{results/}{stores output results}

For more details, please consult the Overview Vignette (HTML) of the systemPipeR package (http://bioconductor.org/packages/systemPipeR).
}
\author{
Thomas Girke, Shiyuan Guo and Daniela Cassol
}
\examples{
\dontrun{
getSubsetReads(args, MappingRegion = 1:900, sample_range = 800:900, outdir = "data/subset/", silent = FALSE)
getSubsetReads(args, MappingRegion = 1:900, sample_range = NULL, outdir = "data/subset/", silent = FALSE)
}
}
\keyword{ utilities }
