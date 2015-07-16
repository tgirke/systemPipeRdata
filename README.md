### systemPipeRdata: sample data for the systemPipeR NGS workflow environment

[_systemPipeR_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html)
is an R/Bioconductor package for building *end-to-end* analysis pipelines with
automated report generation for next generation sequence (NGS) applications
such as RNA-Seq, ChIP-Seq, VAR-Seq and many others. An important feature is
support for running command-line software, such as NGS aligners, on both single
machines or compute clusters. This includes both interactive job submissions or
batch submissions to queuing systems of clusters.  Efficient handling of
complex sample sets and experimental designs is facilitated by a well-defined
sample annotation infrastructure which improves reproducibility and
user-friendliness of many typical analysis workflows in the NGS area.

#### Installation 
The easiest way to install _systemPipeRdata_ is via the _devtools_ package from
CRAN.
```s
devtools::install_github("tgirke/systemPipeRdata")
```
Due to the large size of the sample data (~250MB), the install may take some time.

To install _systemPipeR_ itself, please use the _biocLite_ method as instructed 
[_here_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html).

#### Usage
The data provided by _systemPipeRdata_ are explained in the overview
[_vignette_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html) (manual).
