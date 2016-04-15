# Usage

The package builds on StefansExpressionLib - please install that first.

## Install

source("https://bioconductor.org/biocLite.R")

biocLite(c("DESeq", 'Rsubread') )

library(devtools)

install_github('stela2502/RFclust.SGE')

install_github('stela2502/StefansExpressionSet')

install_github('stela2502/NGSexpressionSet')

or

install_github('stela2502/NGSexpressionSet', build_vignettes=FALSE)