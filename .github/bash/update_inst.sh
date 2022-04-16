#!/bin/bash
set -e
# git_link is the http git link of a repo
git_link=${1}
# update_dir is what is the folder path in website you want to copy to
update_dir=${2}
down_folder="/tmp/${update_dir}"

# print info
echo ${update_dir}
echo ${git_link}

# clean temp
echo "Clean dwonload folder ${down_folder}"
rm -rf ${down_folder}

# clone repo
echo "Downloading ${update_dir} to \"${down_folder}\""
git clone ${git_link} ${down_folder}

# copy
echo "Update inst files"
cp -r ${down_folder}/vignettes/systemPipeRIBOseq.Rmd ${update_dir}/extdata/workflows/riboseq/
cp -r ${down_folder}/vignettes/bibtex.bib ${update_dir}/extdata/workflows/riboseq/
cd ${update_dir}/extdata/workflows/riboseq/
Rscript -e "rmarkdown::render('systemPipeRIBOseq.Rmd', c('BiocStyle::html_document'), clean=T); knitr::knit('systemPipeRIBOseq.Rmd', tangle=TRUE)"
rm -rf *cache/
cd ../../../../
echo "done ribo"

cp -r ${down_folder}/vignettes/systemPipeChIPseq.Rmd ${update_dir}/extdata/workflows/chipseq/
cp -r ${down_folder}/vignettes/bibtex.bib ${update_dir}/extdata/workflows/chipseq/
cd ${update_dir}/extdata/workflows/chipseq/
Rscript -e "rmarkdown::render('systemPipeChIPseq.Rmd', c('BiocStyle::html_document'), clean=T); knitr::knit('systemPipeChIPseq.Rmd', tangle=TRUE)"
rm -rf *cache/
cd ../../../../
echo "done chip"

cp -r ${down_folder}/vignettes/systemPipeRNAseq.Rmd ${update_dir}/extdata/workflows/rnaseq/
cp -r ${down_folder}/vignettes/bibtex.bib ${update_dir}/extdata/workflows/rnaseq/
cd ${update_dir}/extdata/workflows/rnaseq/
Rscript -e "rmarkdown::render('systemPipeRNAseq.Rmd', c('BiocStyle::html_document'), clean=T); knitr::knit('systemPipeRNAseq.Rmd', tangle=TRUE)"
rm -rf *cache/
cd ../../../../
echo "done rna"

cp -r ${down_folder}/vignettes/systemPipeVARseq.Rmd ${update_dir}/extdata/workflows/varseq/
cp -r ${down_folder}/vignettes/bibtex.bib ${update_dir}/extdata/workflows/varseq/
cd ${update_dir}/extdata/workflows/varseq/
Rscript -e "rmarkdown::render('systemPipeVARseq.Rmd', c('BiocStyle::html_document'), clean=T); knitr::knit('systemPipeVARseq.Rmd', tangle=TRUE)"
rm -rf *cache/
cd ../../../../
echo "done var"

# cleaning
echo "Clean update folder \"${down_folder}\""
rm -rf ${down_folder}