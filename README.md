
### systemPipeRdata: NGS workflow templates and sample data

[_systemPipeRdata_](http://bioconductor.org/packages/devel/systemPipeRdata) is a helper 
package to generate with a single command NGS workflow templates that are intended to
be used by its parent package [_systemPipeR_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html). 
The latter is an environment for building *end-to-end* analysis pipelines with
automated report generation for next generation sequence (NGS) applications
such as RNA-Seq, ChIP-Seq, VAR-Seq, Ribo-Seq and many others. 

#### Installation 
_systemPipeRdata_ can be installed directly from GitHub using the [_devtools_](http://cran.r-project.org/web/packages/devtools/index.html) 
package from CRAN.
```s
devtools::install_github("tgirke/systemPipeRdata")
```
Alternatively, one can install it with the Bioconductor _`BiocManager::install`_ command.
```s
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("systemPipeRdata")
```
Due to the large size of the sample data (~320 MB) provided by _systemPipeRdata_, its download/install may take some time.

To install the parent package _systemPipeR_ itself, please use the _`BiocManager`_ package as instructed 
[_here_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html).

#### Usage
Detailed user manuals are available here: 
+ [_systemPipeRdata_ Vignette](https://htmlpreview.github.io/?https://github.com/tgirke/systemPipeRdata/blob/master/vignettes/systemPipeRdata.html)
+ [_systemPipeR_ Overview Vignette](https://htmlpreview.github.io/?https://github.com/tgirke/systemPipeR/blob/master/vignettes/systemPipeR.html)

Additional information can be found on the corresponding Bioconductor package overview 
[_page_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html).



