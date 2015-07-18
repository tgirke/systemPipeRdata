### systemPipeRdata: sample data for the systemPipeR workflow environment

[_systemPipeR_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html)
is an R/Bioconductor package for building *end-to-end* analysis pipelines with
automated report generation for next generation sequence (NGS) applications
such as RNA-Seq, ChIP-Seq, VAR-Seq and many others. An important feature is
support for running command-line software, such as NGS aligners, on both single
machines or compute clusters. This includes both interactive job submissions or
batch submissions to queuing systems of clusters.  Efficient handling of
complex sample sets and experimental designs is facilitated by a consistently
designed sample annotation infrastructure which improves reproducibility and
user-friendliness of many typical analysis workflows in the NGS area. The
input data for _systemPipeR's_ sample workflows are provided by _systemPipeRdata_.
Its vignette is available [_here_](https://github.com/tgirke/systemPipeRdata/blob/master/vignettes/systemPipeRdata.pdf?raw=true).

The presentation for the Bioc2015 event is available [_here_](https://htmlpreview.github.io/?https://github.com/tgirke/systemPipeRdata/blob/master/vignettes/systemPipeR_Presentation.html).

#### Installation 
The easiest way to install _systemPipeRdata_ is via the [_devtools_](http://cran.r-project.org/web/packages/devtools/index.html) 
package from CRAN.
```s
devtools::install_github("tgirke/systemPipeRdata")
```

Alternatively, one can install the package via _biocLite_.
```s
source("http://bioconductor.org/biocLite.R")
biocLite("tgirke/systemPipeRdata", build_vignettes=TRUE, dependencies=TRUE)
```

Due to the large size of the sample data (~280MB), the download/install may take some time.

To install _systemPipeR_ itself, please use the _biocLite_ method as instructed 
[_here_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html).

#### Usage
The data provided by _systemPipeRdata_ are explained in _systemPipeR's_ overview
[_vignette_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html) (manual).
