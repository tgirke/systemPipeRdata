
### Material for Bioc2016 Tutorial 
+ [Overview Slide Show](https://docs.google.com/presentation/d/175aup31LvnbIJUAvEEoSkpGsKgtBJ2RpQYd0Gs23dLo/embed?start=false&loop=false&delayms=60000)
+ Introduction 
    + [HTML](https://htmlpreview.github.io/?https://raw.githubusercontent.com/tgirke/systemPipeRdata/master/vignettes/systemPipeR_Intro.html)
    + [Rmd](https://raw.githubusercontent.com/tgirke/systemPipeRdata/master/vignettes/systemPipeR_Intro.Rmd)
+ Demo Workflow - RIBO-Seq 
    + [HTML](https://htmlpreview.github.io/?https://raw.githubusercontent.com/tgirke/systemPipeRdata/master/inst/extdata/workflows/riboseq/systemPipeRIBOseq.html)
    + [Rmd](https://raw.githubusercontent.com/tgirke/systemPipeRdata/master/inst/extdata/workflows/riboseq/systemPipeRIBOseq.Rmd)

#### Installation 
Please install _systemPipeR_ along with its dependencies from Bioconductor (devel branch) and _systemPipeRdata_ from this GitHub repository as shown below.
_systemPipeRdata_ has not many dependencies but provides the data and `Rmd` tutorials. However, the parent package _systemPipeR_ has many dependencies. 

```s
source("http://bioconductor.org/biocLite.R")
biocLite("tgirke/systemPipeRdata", build_vignettes=TRUE, dependencies=TRUE)
biocLite("systemPipeR") # Normal BioC install is sufficient
```

Due to the large size of the sample data (~320 MB) provided by _systemPipeRdata_, its download/install may take some time.

