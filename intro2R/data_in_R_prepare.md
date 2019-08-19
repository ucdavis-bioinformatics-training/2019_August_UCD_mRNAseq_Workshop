---
title: "Prepare Data_in_R"
author: "Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---

### Create a new RStudio project

Open RStudio and create a new project, for more info see <https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects>

* File > New Project > New Directory > New Project (name the new directory, Ex. Data_in_R) and check "use packrat with this project" if present.

Learn more about packrat see <https://rstudio.github.io/packrat/>

Set some options and make sure the packages 'knitr', 'tidyverse', 'reshape2', and 'gridExtra' are installed (if not install it), and then load

In the R console run the following commands

```r
if (!any(rownames(installed.packages()) == "knitr")){
  install.packages("knitr")
}
library(knitr)

if (!any(rownames(installed.packages()) == "tidyverse")){
  install.packages("tidyverse")
}
library(tidyverse)
```

<div class="output">── Attaching packages ─────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──
</div>

<div class="output">✔ ggplot2 3.1.0     ✔ purrr   0.2.5
✔ tibble  2.0.1     ✔ dplyr   0.7.8
✔ tidyr   0.8.2     ✔ stringr 1.3.1
✔ readr   1.3.1     ✔ forcats 0.3.0
</div>

<div class="output">── Conflicts ────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
✖ dplyr::filter() masks stats::filter()
✖ dplyr::lag()    masks stats::lag()
</div>

```r
if (!any(rownames(installed.packages()) == "reshape2")){
  install.packages("reshape2")
}
library(reshape2)
```

<div class="output">
Attaching package: 'reshape2'
</div>

<div class="output">The following object is masked from 'package:tidyr':

    smiths
</div>

```r
if (!any(rownames(installed.packages()) == "gridExtra")){
  install.packages("gridExtra")
}
library(gridExtra)
```

<div class="output">
Attaching package: 'gridExtra'
</div>

<div class="output">The following object is masked from 'package:dplyr':

    combine
</div>

Learn more about the tidyverse see <https://www.tidyverse.org>.

### Open a new R Notebook

An R notebook is an R Markdown document with chunks that can be executed independently and interactively, with output visible immediately beneath the input. More info see <https://rmarkdown.rstudio.com/r_notebooks.html>

* File -> New File -> R Notebook
* Save the Notebook (Ex. test)

### R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **preview** or **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed R code and plots in chunks like this:

<pre><code>```{r chunk_name}
print('hello world!')
```</code></pre>

Review the [R Markdown page]<https://rmarkdown.rstudio.com> and [R Markdown cheat sheets] <https://rmarkdown.rstudio.com/lesson-15.html>.

Try 'knitting' to html, pdf, and doc as well as previewing the notebook. Open the resulting documents.

Try executing the code chunks in the R Notebook.


### Download the data file for the workshop document and preview/open it

This is the stats file generated after running samtools stats on a bam file generated from running BWA MEM.

In the R console run the following command.

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019_August_UCD_mRNAseq_Workshop/master/intro2R/Data_in_R_files/bwa_mem_Stats.log", "bwa_mem_Stats.log")
```

### Download the template Markdown workshop document and open it

In the R console run the following command

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019_August_UCD_mRNAseq_Workshop/master/intro2R/data_in_R.Rmd", "data_in_R.Rmd")
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


### What are we going to do?

We will recreate some of the plots generated with plot-bamstats on the same file

You can view the output of plot-bamstats -> [bwa_mem_stats.html](Data_in_R_files/bwa_mem_Stats/bwa_mem_Stats.html)
