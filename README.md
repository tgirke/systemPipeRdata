
# systemPipeRdata: Workflow templates and sample data <img src="https://github.com/tgirke/systemPipeR/raw/gh-pages/images/systemPipeR.png" align="right" height="139" />

[![platforms](http://www.bioconductor.org/images/shields/availability/all.svg)](http://www.bioconductor.org/packages/devel/data/experiment/html/systemPipeRdata.html#archives)
[![rank](http://www.bioconductor.org/shields/downloads/devel/systemPipeRdata.svg)](http://bioconductor.org/packages/stats/data-experiment/systemPipeRdata/)
[![posts](http://www.bioconductor.org/shields/posts/systemPipeRdata.svg)](https://support.bioconductor.org/t/systempiperdata/)
[![build](http://www.bioconductor.org/shields/build/devel/data-experiment/systemPipeRdata.svg)](http://bioconductor.org/checkResults/devel/data-experiment-LATEST/systemPipeRdata/)
[![updated](http://www.bioconductor.org/shields/lastcommit/devel/data-experiment/systemPipeRdata.svg)](http://bioconductor.org/checkResults/devel/data-experiment-LATEST/systemPipeRdata/)
[![dependencies](http://www.bioconductor.org/shields/dependencies/devel/systemPipeRdata.svg)](http://www.bioconductor.org/packages/devel/data/experiment/html/systemPipeRdata.html#since)

![R-CMD-check](https://github.com/tgirke/systemPipeRdata/workflows/R-CMD-check/badge.svg)

[_systemPipeRdata_](http://bioconductor.org/packages/devel/systemPipeRdata) is a helper 
package to generate with a single command NGS workflow templates that are intended to
be used by its parent package [_systemPipeR_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html). 
The latter is an environment for building *end-to-end* analysis pipelines with
automated report generation for next generation sequence (NGS) applications
such as RNA-Seq, ChIP-Seq, VAR-Seq, Ribo-Seq and many others. 

#### Installation 

To install the package, please use the _`BiocManager::install`_ command:
```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("systemPipeRdata")
```

To obtain the most recent updates immediately, one can install it directly from
github as follow:
```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("tgirke/systemPipeRdata", build_vignettes=TRUE, dependencies=TRUE)
```

Due to the large size of the sample data (~320 MB) provided by _systemPipeRdata_, its download/install may take some time.

To install the parent package _systemPipeR_ itself, please use the _BiocManager::install_ method as instructed
[_here_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html).

#### Usage

Detailed user manuals are available here: 

+ [_systemPipeRdata_ Vignette](http://www.bioconductor.org/packages/devel/data/experiment/vignettes/systemPipeRdata/inst/doc/systemPipeRdata.html)
+ [_systemPipeR_ Overview Vignette](http://www.bioconductor.org/packages/devel/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html)

Additional information can be found on the corresponding Bioconductor package overview 
[_page_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html).
