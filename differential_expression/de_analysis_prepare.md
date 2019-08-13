---
title: "Prepare Differential Expression"
author: "Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---

### Create a new RStudio project

Open RStudio and create a new project, for more info see <https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects>

* File > New Project > New Directory > New Project (name the new directory, Ex. Differential Expression) and check "use packrat with this project" if present.

Learn more about packrat see <https://rstudio.github.io/packrat/>

Set some options and make sure the packages edgeR, gplots, RColorBrewer, topGO, KEGGREST, Rgraphviz and org.Hs.eg.db are installed (if not install it), and then load

In the R console run the following commands

```r
if (!any(rownames(installed.packages()) == "edgeR")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("edgeR", version = "3.8")
}
library(edgeR)
```

```
## Loading required package: limma
```

```r
if (!any(rownames(installed.packages()) == "topGO")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("topGO", version = "3.8")
}
```

```
## Bioconductor version 3.8 (BiocManager 1.30.4), R 3.5.2 (2018-12-20)
```

```
## Installing package(s) 'topGO'
```

```
## also installing the dependencies 'graph', 'SparseM'
```

```
## Update old packages: 'BiocParallel', 'caTools', 'class', 'codetools',
##   'dplyr', 'evaluate', 'flexmix', 'forcats', 'formatR', 'GenomeInfoDb',
##   'gplots', 'haven', 'Hmisc', 'igraph', 'irlba', 'knitr', 'later',
##   'Matrix', 'metap', 'mgcv', 'modelr', 'mvtnorm', 'openssl', 'pbapply',
##   'processx', 'proxy', 'purrr', 'R.utils', 'R6', 'RCurl', 'readxl',
##   'reticulate', 'Rsamtools', 'stringi', 'stringr', 'sys', 'tidyr',
##   'vegan', 'xfun'
```

```r
library(topGO)
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following object is masked from 'package:limma':
## 
##     plotMA
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind,
##     colMeans, colnames, colSums, dirname, do.call, duplicated,
##     eval, evalq, Filter, Find, get, grep, grepl, intersect,
##     is.unsorted, lapply, lengths, Map, mapply, match, mget, order,
##     paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind,
##     Reduce, rowMeans, rownames, rowSums, sapply, setdiff, sort,
##     table, tapply, union, unique, unsplit, which, which.max,
##     which.min
```

```
## Loading required package: graph
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Loading required package: GO.db
```

```
## Loading required package: AnnotationDbi
```

```
## Loading required package: stats4
```

```
## Loading required package: IRanges
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```
## 
```

```
## Loading required package: SparseM
```

```
## 
## Attaching package: 'SparseM'
```

```
## The following object is masked from 'package:base':
## 
##     backsolve
```

```
## 
## groupGOTerms: 	GOBPTerm, GOMFTerm, GOCCTerm environments built.
```

```
## 
## Attaching package: 'topGO'
```

```
## The following object is masked from 'package:IRanges':
## 
##     members
```

```r
if (!any(rownames(installed.packages()) == "KEGGREST")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("KEGGREST", version = "3.8")
}
```

```
## Bioconductor version 3.8 (BiocManager 1.30.4), R 3.5.2 (2018-12-20)
```

```
## Installing package(s) 'KEGGREST'
```

```
## Update old packages: 'BiocParallel', 'caTools', 'class', 'codetools',
##   'dplyr', 'evaluate', 'flexmix', 'forcats', 'formatR', 'GenomeInfoDb',
##   'gplots', 'haven', 'Hmisc', 'igraph', 'irlba', 'knitr', 'later',
##   'Matrix', 'metap', 'mgcv', 'modelr', 'mvtnorm', 'openssl', 'pbapply',
##   'processx', 'proxy', 'purrr', 'R.utils', 'R6', 'RCurl', 'readxl',
##   'reticulate', 'Rsamtools', 'stringi', 'stringr', 'sys', 'tidyr',
##   'vegan', 'xfun'
```

```r
library(KEGGREST)

if (!any(rownames(installed.packages()) == "Rgraphviz")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("Rgraphviz", version = "3.8")
}
```

```
## Bioconductor version 3.8 (BiocManager 1.30.4), R 3.5.2 (2018-12-20)
```

```
## Installing package(s) 'Rgraphviz'
```

```
## Update old packages: 'BiocParallel', 'caTools', 'class', 'codetools',
##   'dplyr', 'evaluate', 'flexmix', 'forcats', 'formatR', 'GenomeInfoDb',
##   'gplots', 'haven', 'Hmisc', 'igraph', 'irlba', 'knitr', 'later',
##   'Matrix', 'metap', 'mgcv', 'modelr', 'mvtnorm', 'openssl', 'pbapply',
##   'processx', 'proxy', 'purrr', 'R.utils', 'R6', 'RCurl', 'readxl',
##   'reticulate', 'Rsamtools', 'stringi', 'stringr', 'sys', 'tidyr',
##   'vegan', 'xfun'
```

```r
library(Rgraphviz)
```

```
## Loading required package: grid
```

```
## 
## Attaching package: 'grid'
```

```
## The following object is masked from 'package:topGO':
## 
##     depth
```

```
## 
## Attaching package: 'Rgraphviz'
```

```
## The following objects are masked from 'package:IRanges':
## 
##     from, to
```

```
## The following objects are masked from 'package:S4Vectors':
## 
##     from, to
```

```r
if (!any(rownames(installed.packages()) == "org.Hs.eg.db")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db", version = "3.8")
}
```

```
## Bioconductor version 3.8 (BiocManager 1.30.4), R 3.5.2 (2018-12-20)
```

```
## Installing package(s) 'org.Hs.eg.db'
```

```
## installing the source package 'org.Hs.eg.db'
```

```
## Update old packages: 'BiocParallel', 'caTools', 'class', 'codetools',
##   'dplyr', 'evaluate', 'flexmix', 'forcats', 'formatR', 'GenomeInfoDb',
##   'gplots', 'haven', 'Hmisc', 'igraph', 'irlba', 'knitr', 'later',
##   'Matrix', 'metap', 'mgcv', 'modelr', 'mvtnorm', 'openssl', 'pbapply',
##   'processx', 'proxy', 'purrr', 'R.utils', 'R6', 'RCurl', 'readxl',
##   'reticulate', 'Rsamtools', 'stringi', 'stringr', 'sys', 'tidyr',
##   'vegan', 'xfun'
```

```r
library(org.Hs.eg.db)
```

```
## 
```

```r
library(edgeR)
if (!any(rownames(installed.packages()) == "gplots")){
install.packages("gplots")
}
library(gplots)
```

```
## 
## Attaching package: 'gplots'
```

```
## The following object is masked from 'package:IRanges':
## 
##     space
```

```
## The following object is masked from 'package:S4Vectors':
## 
##     space
```

```
## The following object is masked from 'package:stats':
## 
##     lowess
```

```r
if (!any(rownames(installed.packages()) == "RColorBrewer")){
install.packages("RColorBrewer")
}
library(RColorBrewer)
```


### Download the template Markdown workshop document and open it

In the R console run the following command

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019_August_UCD_mRNAseq_Workshop/master/differential_expression/DE_Analysis.Rmd", "DE_Analysis.Rmd")
```

### Download the data file for the workshop document and preview/open it

This is the the counts file generated after running [Generating Summarized Counts](https://ucdavis-bioinformatics-training.github.io/2019_August_UCD_mRNAseq_Workshop/data_reduction/counts.html). I've also uploaded to the github pages 

In the R console run the following command.

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019_August_UCD_mRNAseq_Workshop/master/differential_expression/de_data/rnaseq_workshop_counts.txt", "rnaseq_workshop_counts.txt")
```

### Edit the file YAML portion

The top YAML (YAML ain't markup language) portion of the doc tells RStudio how to parse the document.

<pre><code>---
title: "Data_in_R"
author: your_name
date: current_date
output:
    html_notebook: default
    html_document: default
---</code></pre>

