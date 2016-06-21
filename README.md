
### Material for Bioc2016 Tutorial 
+ [Intro Slide Show](https://docs.google.com/presentation/d/175aup31LvnbIJUAvEEoSkpGsKgtBJ2RpQYd0Gs23dLo/embed?start=false&loop=false&delayms=60000)
+ Tutorial Material 
    + Introduction to `systemPipeR` 
        + [HTML](https://htmlpreview.github.io/?https://raw.githubusercontent.com/tgirke/systemPipeRdata/master/vignettes/systemPipeR_Intro.html)
        + [Rmd](https://raw.githubusercontent.com/tgirke/systemPipeRdata/master/vignettes/systemPipeR_Intro.Rmd)
    + Demo Workflow: RIBO-Seq 
        + [HTML](https://htmlpreview.github.io/?https://raw.githubusercontent.com/tgirke/systemPipeRdata/master/inst/extdata/workflows/riboseq/systemPipeRIBOseq.html)
        + [Rmd](https://raw.githubusercontent.com/tgirke/systemPipeRdata/master/inst/extdata/workflows/riboseq/systemPipeRIBOseq.Rmd)

#### Installation 
Please install _systemPipeRdata_ from this GitHub repository as shown below. This package provides the data and `Rmd` files for the tutorial. 
Its parent package _systemPipeR_ is a dependency and should install along with all its own dependencies automatically. If it doesn't then please also
run the last `biocLite` install command given below.

```s
source("http://bioconductor.org/biocLite.R")
biocLite("tgirke/systemPipeRdata", build_vignettes=TRUE, dependencies=TRUE)
biocLite("systemPipeR")
```

Note, due to the relative large size of the sample data (~320 MB) provided by _systemPipeRdata_, its download/install may take some time.

