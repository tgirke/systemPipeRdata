### systemPipeRdata: NGS workflow templates and sample data

[_systemPipeRdata_](https://github.com/tgirke/systemPipeRdata) is an
R/Bioconductor data package that facilitates the generation of NGS workflow
templates for its parent package 
[_systemPipeR_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html). 
The latter is an environment for building *end-to-end* analysis pipelines with
automated report generation for next generation sequence (NGS) applications
such as RNA-Seq, ChIP-Seq, VAR-Seq and many others. A short manual (vignette) for _systemPipeRdata_ 
is available [_here_](https://github.com/tgirke/systemPipeRdata/blob/master/vignettes/systemPipeRdata.pdf?raw=true),
and the main manual for _systemPipeR_ is available
[_here_](https://htmlpreview.github.io/?https://github.com/tgirke/systemPipeR/blob/master/vignettes/systemPipeR.html).

#### Installation 
_systemPipeRdata_ can be installed directly from GitHub using the [_devtools_](http://cran.r-project.org/web/packages/devtools/index.html) 
package from CRAN.
```s
devtools::install_github("tgirke/systemPipeRdata")
```

Alternatively, one can install it using the _biocLite_ command from Bioconductor.
```s
source("http://bioconductor.org/biocLite.R")
biocLite("tgirke/systemPipeRdata", build_vignettes=TRUE, dependencies=TRUE)
```

Due to the large size of the sample data (~250MB) provided by _systemPipeRdata_, its download/install may take some time.

To install the parent package _systemPipeR_ itself, please use the _biocLite_ method as instructed 
[_here_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html).

#### Usage
Examples for using _systemPipeRdata_ are given in _systemPipeR's_ overview
[_vignette_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html).
