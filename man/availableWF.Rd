\name{availableWF}
\alias{availableWF}
\title{
List Available Workflows Templates at \code{systemPipeRdata} package
}
\description{
This function checks the workflow templates availability from \code{systemPipeRdata} package and 
also from \href{https://github.com/systemPipeR/}{systemPipeR Organization} on GitHub.
}
\usage{
availableWF(github = FALSE)
}
\arguments{
  \item{github}{
  logical. If \code{TRUE}, it will return current workflow templates available on \code{systemPipeR} Organization.
}
}
\value{
Return a list with the names of the workflows templates availabe at \code{systemPipeRdata} 
package. If \code{github = TRUE}, it will return an addicional \code{data.frame}
with current workflow templates available on \code{systemPipeR} Organization.
}
\author{
Daniela Cassol
}
\details{
Internally, this function uses the GitHub API, and there is an access limit per hour. 
For more details, please check: \code{system("curl -i https://api.github.com/users/<username>")}.

}
\note{
We are assuming that workflow templates repositories under
\href{https://github.com/systemPipeR/}{systemPipeR Organization} content the keyword 
"Workflow Template" on the Description section and "Topics" section, we expected 
"systempiper" and "release" or "development" words.
}
\seealso{
\code{\link{genWorkenvir}}.
}
\examples{
availableWF()
\dontrun{
## List Workflow Templates from \code{systemPipeR} Organization
availableWF(github = TRUE)
}
}
\keyword{ utilities }