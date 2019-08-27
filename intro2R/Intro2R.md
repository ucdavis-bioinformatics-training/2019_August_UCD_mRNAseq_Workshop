---
title: "R for Biologist - An Introduction to R"
author: "Bioinformatics Core"
date: "2019-08-27"
output:
  html_document:
    keep_md: TRUE
---

<style type="text/css">
.colsel {
background-color: lightyellow;
}
</style>



R and RStudio
========================================================

### What is R?
[R](http://r-project.org/) is a language and environment for statistical computing and graphics developed in 1993. It provides a wide variety of statistical and graphical techniques (linear and nonlinear modeling, statistical tests, time series analysis, classification, clustering, ...), and is highly extensible, meaning that the user community can write new R tools. It is a GNU project (Free and Open Source).

The R language has its roots in the S language and environment which was developed at Bell Laboratories (formerly AT&T, now Lucent Technologies) by John Chambers and colleagues. R was created by Ross Ihaka and Robert Gentleman at the University of Auckland, New Zealand, and now, R is developed by the R Development Core Team, of which Chambers is a member. R is named partly after the first names of the first two R authors (Robert Gentleman and Ross Ihaka), and partly as a play on the name of S. R can be considered as a different implementation of S. There are some important differences, but much code written for S runs unaltered under R.

Some of R's strengths:

* The ease with which well-designed publication-quality plots can be produced, including mathematical symbols and formulae where needed. Great care has been taken over the defaults for the minor design choices in graphics, but the user retains full control.
* It compiles and runs on a wide variety of UNIX platforms and similar systems (including FreeBSD and Linux), Windows and MacOS.
* R can be extended (easily) via packages.
* R has its own LaTeX-like documentation format, which is used to supply comprehensive documentation, both on-line in a number of formats and in hardcopy.
* It has a vast community both in academia and in business.
* It's FREE!   

### The R environment
R is an integrated suite of software facilities for data manipulation, calculation and graphical display. It includes

* an effective data handling and storage facility,
* a suite of operators for calculations on arrays, in particular matrices,
* a large, coherent, integrated collection of intermediate tools for data analysis,
* graphical facilities for data analysis and display either on-screen or on hardcopy, and
* a well-developed, and effective programming language which includes conditionals, loops, user-defined recursive functions and input and output facilities.

The term "environment" is intended to characterize it as a fully planned and coherent system, rather than an incremental accretion of very specific and inflexible tools, as is frequently the case with other data analysis software.

R, like S, is designed around a true computer language, and it allows users to add additional functionality by defining new functions. Much of the system is itself written in the R dialect of S, which makes it easy for users to follow the algorithmic choices made. For computationally-intensive tasks, C, C++ and Fortran code can be linked and called at run time. Advanced users can write C code to manipulate R objects directly.

Many users think of R as a statistics system. The R group prefers to think of it of an environment within which statistical techniques are implemented.   

### The R Homepage
The R homepage has a wealth of information on it,

[R-project.org](http://r-project.org/)

On the homepage you can:

* Learn more about R
* Download R
* Get Documentation (official and user supplied)
* Get access to CRAN 'Comprehensive R archival network'

### Interface for R

There are many ways one can interface with R language. Here are a few popular ones:

* [RStudio](https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019-March-Bioinformatics-Prerequisites/master/wednesday/Intro2R/rstudio.png)
* [RGui](https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019-March-Bioinformatics-Prerequisites/master/wednesday/Intro2R/rgui.png)
* Jupyter and R notebooks
* text editors, such as vi(m), Emacs...


### RStudio

[RStudio](http://rstudio.com/) started in 2010, to offer R a more full featured integrated development environment (IDE) and modeled after matlab's IDE.

RStudio has many features:

* syntax highlighting
* code completion
* smart indentation
* "Projects"
* workspace browser and data viewer
* embedded plots
* Markdown notebooks, Sweave authoring and knitr with one click pdf or html
* runs on all platforms and over the web
* etc. etc. etc.

RStudio and its team have contributed to many R packages.[13] These include:

* Tidyverse – R packages for data science, including ggplot2, dplyr, tidyr, and purrr
* Shiny – An interactive web technology
* RMarkdown – Insert R code into markdown documents
* knitr – Dynamic reports combining R, TeX, Markdown & HTML
* packrat – Package dependency tool
* devtools – Package development tool

RStudio Cheat Sheets: [rstudio-ide.pdf](https://github.com/rstudio/cheatsheets/raw/master/rstudio-ide.pdf)

---


Topics covered in this introduction to R
====================================================
1. Basic concepts
2. Basic data types in R
3. Import and export data in R
4. Functions in R
5. Basic statistics in R
6. Simple data visulization in R
7. Install packages in R
8. Save data in R session


---

Topic 1. Basic concepts
====================================================

There are three concepts that we should be familiar with before working in R:

* Operators

<table class="table table-striped" style="width: auto !important; ">
<caption>Assignment Operators in R</caption>
 <thead>
  <tr>
   <th style="text-align:center;"> Operator </th>
   <th style="text-align:center;"> Description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> &lt;-, = </td>
   <td style="text-align:center;"> Assignment </td>
  </tr>
</tbody>
</table>

<table class="table table-striped" style="width: auto !important; ">
<caption>Arithmetic Operators in R</caption>
 <thead>
  <tr>
   <th style="text-align:center;"> Operator </th>
   <th style="text-align:center;"> Description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> \+ </td>
   <td style="text-align:center;"> Addition </td>
  </tr>
  <tr>
   <td style="text-align:center;"> \- </td>
   <td style="text-align:center;"> Subtraction </td>
  </tr>
  <tr>
   <td style="text-align:center;"> \* </td>
   <td style="text-align:center;"> Multiplication </td>
  </tr>
  <tr>
   <td style="text-align:center;"> / </td>
   <td style="text-align:center;"> Division </td>
  </tr>
  <tr>
   <td style="text-align:center;"> ^ </td>
   <td style="text-align:center;"> Exponent </td>
  </tr>
  <tr>
   <td style="text-align:center;"> %% </td>
   <td style="text-align:center;"> Modulus </td>
  </tr>
  <tr>
   <td style="text-align:center;"> %/% </td>
   <td style="text-align:center;"> Integer Division </td>
  </tr>
</tbody>
</table>

<table class="table table-striped" style="width: auto !important; ">
<caption>Relational Operators in R</caption>
 <thead>
  <tr>
   <th style="text-align:center;"> Operator </th>
   <th style="text-align:center;"> Description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> &lt; </td>
   <td style="text-align:center;"> Less than </td>
  </tr>
  <tr>
   <td style="text-align:center;"> &gt; </td>
   <td style="text-align:center;"> Greater than </td>
  </tr>
  <tr>
   <td style="text-align:center;"> &lt;= </td>
   <td style="text-align:center;"> Less than or equal to </td>
  </tr>
  <tr>
   <td style="text-align:center;"> &gt;= </td>
   <td style="text-align:center;"> Greater than or equal to </td>
  </tr>
  <tr>
   <td style="text-align:center;"> == </td>
   <td style="text-align:center;"> Equal to </td>
  </tr>
  <tr>
   <td style="text-align:center;"> != </td>
   <td style="text-align:center;"> Not equal to </td>
  </tr>
</tbody>
</table>

<table class="table table-striped" style="width: auto !important; ">
<caption>Logical Operators in R</caption>
 <thead>
  <tr>
   <th style="text-align:center;"> Operator </th>
   <th style="text-align:center;"> Description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> ! </td>
   <td style="text-align:center;"> Logical NOT </td>
  </tr>
  <tr>
   <td style="text-align:center;"> &amp; </td>
   <td style="text-align:center;"> Element-wise logical AND </td>
  </tr>
  <tr>
   <td style="text-align:center;"> &amp;&amp; </td>
   <td style="text-align:center;"> Logical AND </td>
  </tr>
  <tr>
   <td style="text-align:center;"> | </td>
   <td style="text-align:center;"> Element-wise logical OR </td>
  </tr>
  <tr>
   <td style="text-align:center;"> || </td>
   <td style="text-align:center;"> Logical OR </td>
  </tr>
</tbody>
</table>

* Functions

Functions are essential in all programming languages. A function takes zero or more parameters and return a result. The way to use a function in R is:

function.name(parameter1=value1, ...)


* Variables

A variable is a named storage. The name of a variable can have letters, numbers, dot and underscore. However, a valid variable name cannot start with a underscore or a number, or start with a dot that is followed by a number.


---

Topic 2. Basic data types in R
====================================================

### Simple variables: variables that have a numeric value, a character value (such as a string), or a logical value (True or False)

Examples of numeric values.

```{.r .colsel}
# assign number 150 to variable a.
a <- 150
a
```

```
## [1] 150
```

```{.r .colsel}
# assign a number in scientific format to variable b.
b <- 3e-2
b
```

```
## [1] 0.03
```

Examples of character values.

```{.r .colsel}
# assign a string "Professor" to variable title
title <- "Professor"
title
```

```
## [1] "Professor"
```

```{.r .colsel}
# assign a string "Hello World" to variable hello
hello <- "Hello World"
hello
```

```
## [1] "Hello World"
```

Examples of logical values.

```{.r .colsel}
# assign logical value "TRUE" to variable is_female
is_female <- TRUE
is_female
```

```
## [1] TRUE
```

```{.r .colsel}
# assign logical value "FALSE" to variable is_male
is_male <- FALSE
is_male
```

```
## [1] FALSE
```

```{.r .colsel}
# assign logical value to a variable by logical operation
age <- 20
is_adult <- age > 18
is_adult
```

```
## [1] TRUE
```

To find out the type of variable.

```{.r .colsel}
class(is_female)
```

```
## [1] "logical"
```

```{.r .colsel}
# To check whether the variable is a specific type
is.numeric(hello)
```

```
## [1] FALSE
```

```{.r .colsel}
is.numeric(a)
```

```
## [1] TRUE
```

```{.r .colsel}
is.character(hello)
```

```
## [1] TRUE
```

The rule to convert a logical variable to numeric: TRUE > 1, FALSE > 0

```{.r .colsel}
as.numeric(is_female)
```

```
## [1] 1
```

```{.r .colsel}
as.numeric(is_male)
```

```
## [1] 0
```

R does not know how to convert a numeric variable to a character variable.

```{.r .colsel}
b
```

```
## [1] 0.03
```

```{.r .colsel}
as.character(b)
```

```
## [1] "0.03"
```


### Vectors: a vector is a combination of multiple values(numeric, character or logical) in the same object. A vector is created using the function c() (for concatenate).


```{.r .colsel}
friend_ages <- c(21, 27, 26, 32)
friend_ages
```

```
## [1] 21 27 26 32
```

```{.r .colsel}
friend_names <- c("Mina", "Ella", "Anna", "Cora")
friend_names
```

```
## [1] "Mina" "Ella" "Anna" "Cora"
```

One can give names to the elements of a vector.

```{.r .colsel}
# assign names to a vector by specifying them
names(friend_ages) <- c("Mina", "Ella", "Anna", "Carla")
friend_ages
```

```
##  Mina  Ella  Anna Carla 
##    21    27    26    32
```

```{.r .colsel}
# assign names to a vector using another vector
names(friend_ages) <- friend_names
friend_ages
```

```
## Mina Ella Anna Cora 
##   21   27   26   32
```

Or One may create a vector with named elements from scratch.

```{.r .colsel}
friend_ages <- c(Mina=21, Ella=27, Anna=26, Cora=32)
friend_ages
```

```
## Mina Ella Anna Cora 
##   21   27   26   32
```

To find out the length of a vector:

```{.r .colsel}
length(friend_ages)
```

```
## [1] 4
```

To access elements of a vector: by index, or by name if it is a named vector.

```{.r .colsel}
friend_ages[2]
```

```
## Ella 
##   27
```

```{.r .colsel}
friend_ages["Ella"]
```

```
## Ella 
##   27
```

```{.r .colsel}
friend_ages[c(1,3)]
```

```
## Mina Anna 
##   21   26
```

```{.r .colsel}
friend_ages[c("Mina", "Anna")]
```

```
## Mina Anna 
##   21   26
```

```{.r .colsel}
# selecting elements of a vector by excluding some of them.
friend_ages[-3]
```

```
## Mina Ella Cora 
##   21   27   32
```

To select a subset of a vector can be done by logical vector.

```{.r .colsel}
my_friends <- c("Mina", "Ella", "Anna", "Cora")
my_friends
```

```
## [1] "Mina" "Ella" "Anna" "Cora"
```

```{.r .colsel}
has_child <- c("TRUE", "TRUE", "FALSE", "TRUE")
has_child
```

```
## [1] "TRUE"  "TRUE"  "FALSE" "TRUE"
```

```{.r .colsel}
my_friends[has_child == "TRUE"]
```

```
## [1] "Mina" "Ella" "Cora"
```


##### <font color='red'>NOTE: a vector can only hold elements of the same type.</font>  



### Matrices: A matrix is like an Excel sheet containing multiple rows and columns. It is used to combine vectors of the same type.


```{.r .colsel}
col1 <- c(1,3,8,9)
col2 <- c(2,18,27,10)
col3 <- c(8,37,267,19)

my_matrix <- cbind(col1, col2, col3)
my_matrix
```

```
##      col1 col2 col3
## [1,]    1    2    8
## [2,]    3   18   37
## [3,]    8   27  267
## [4,]    9   10   19
```

```{.r .colsel}
rownames(my_matrix) <- c("row1", "row2", "row3", "row4")
my_matrix
```

```
##      col1 col2 col3
## row1    1    2    8
## row2    3   18   37
## row3    8   27  267
## row4    9   10   19
```

```{.r .colsel}
t(my_matrix)
```

```
##      row1 row2 row3 row4
## col1    1    3    8    9
## col2    2   18   27   10
## col3    8   37  267   19
```

To find out the dimension of a matrix:

```{.r .colsel}
ncol(my_matrix)
```

```
## [1] 3
```

```{.r .colsel}
nrow(my_matrix)
```

```
## [1] 4
```

```{.r .colsel}
dim(my_matrix)
```

```
## [1] 4 3
```


Accessing elements of a matrix is done in similar ways to accessing elements of a vector.


```{.r .colsel}
my_matrix[1,3]
```

```
## [1] 8
```

```{.r .colsel}
my_matrix["row1", "col3"]
```

```
## [1] 8
```


```{.r .colsel}
my_matrix[1,]
```

```
## col1 col2 col3 
##    1    2    8
```

```{.r .colsel}
my_matrix[,3]
```

```
## row1 row2 row3 row4 
##    8   37  267   19
```


```{.r .colsel}
my_matrix[col3 > 20,]
```

```
##      col1 col2 col3
## row2    3   18   37
## row3    8   27  267
```


Calculations with matrices.


```{.r .colsel}
my_matrix * 3
```

```
##      col1 col2 col3
## row1    3    6   24
## row2    9   54  111
## row3   24   81  801
## row4   27   30   57
```

```{.r .colsel}
log10(my_matrix)
```

```
##           col1     col2     col3
## row1 0.0000000 0.301030 0.903090
## row2 0.4771213 1.255273 1.568202
## row3 0.9030900 1.431364 2.426511
## row4 0.9542425 1.000000 1.278754
```

Total of each row.

```{.r .colsel}
rowSums(my_matrix)
```

```
## row1 row2 row3 row4 
##   11   58  302   38
```

Total of each column.

```{.r .colsel}
colSums(my_matrix)
```

```
## col1 col2 col3 
##   21   57  331
```

It is also possible to use the function apply() to apply any statistical functions to rows/columns of matrices. The advantage of using apply() is that it can take a function created by user.

The simplified format of apply() is as following:

apply(X, MARGIN, FUN)

X: data matrix
MARGIN: possible values are 1 (for rows) and 2 (for columns)
FUN: the function to apply on rows/columns

To calculate the mean of each row.


```{.r .colsel}
apply(my_matrix, 1, mean)
```

```
##       row1       row2       row3       row4 
##   3.666667  19.333333 100.666667  12.666667
```

To calculate the median of each row


```{.r .colsel}
apply(my_matrix, 1, median)
```

```
## row1 row2 row3 row4 
##    2   18   27   10
```


### Factors: a factor represents categorical or groups in data. The function factor() can be used to create a factor variable.


```{.r .colsel}
friend_groups <- factor(c(1,2,1,2))
friend_groups
```

```
## [1] 1 2 1 2
## Levels: 1 2
```

In R, categories are called factor levels. The function levels() can be used to access the factor levels.


```{.r .colsel}
levels(friend_groups)
```

```
## [1] "1" "2"
```

Change the factor levels.


```{.r .colsel}
levels(friend_groups) <- c("best_friend", "not_best_friend")
friend_groups
```

```
## [1] best_friend     not_best_friend best_friend     not_best_friend
## Levels: best_friend not_best_friend
```

Change the order of levels.


```{.r .colsel}
levels(friend_groups) <- c("not_best_friend", "best_friend")
friend_groups
```

```
## [1] not_best_friend best_friend     not_best_friend best_friend    
## Levels: not_best_friend best_friend
```

By default, the order of factor levels is taken in the order of numeric or alphabetic.


```{.r .colsel}
friend_groups <- factor(c("not_best_friend", "best_friend", "not_best_friend", "best_friend"))
friend_groups
```

```
## [1] not_best_friend best_friend     not_best_friend best_friend    
## Levels: best_friend not_best_friend
```

The factor levels can be specified when creating the factor, if the order does not follow the default rule.


```{.r .colsel}
friend_groups <- factor(c("not_best_friend", "best_friend", "not_best_friend", "best_friend"), levels=c("not_best_friend", "best_friend"))
friend_groups
```

```
## [1] not_best_friend best_friend     not_best_friend best_friend    
## Levels: not_best_friend best_friend
```

If you want to know the number of individuals at each levels, there are two functions.


```{.r .colsel}
summary(friend_groups)
```

```
## not_best_friend     best_friend 
##               2               2
```


```{.r .colsel}
table(friend_groups)
```

```
## friend_groups
## not_best_friend     best_friend 
##               2               2
```


### Data frames: a data frame is like a matrix but can have columns with different types (numeric, character, logical).

A data frame can be created using the function data.frame().


```{.r .colsel}
# creating a data frame using previously defined vectors
friends <- data.frame(name=friend_names, age=friend_ages, child=has_child)
friends
```

```
##      name age child
## Mina Mina  21  TRUE
## Ella Ella  27  TRUE
## Anna Anna  26 FALSE
## Cora Cora  32  TRUE
```

To check whether a data is a data frame, use the function is.data.frame().


```{.r .colsel}
is.data.frame(friends)
```

```
## [1] TRUE
```


```{.r .colsel}
is.data.frame(my_matrix)
```

```
## [1] FALSE
```

One can convert a object to a data frame using the function as.data.frame().


```{.r .colsel}
class(my_matrix)
```

```
## [1] "matrix"
```


```{.r .colsel}
my_data <- as.data.frame(my_matrix)
class(my_data)
```

```
## [1] "data.frame"
```

A data frame can be transposed in the similar way as a matrix.


```{.r .colsel}
my_data
```

```
##      col1 col2 col3
## row1    1    2    8
## row2    3   18   37
## row3    8   27  267
## row4    9   10   19
```


```{.r .colsel}
t(my_data)
```

```
##      row1 row2 row3 row4
## col1    1    3    8    9
## col2    2   18   27   10
## col3    8   37  267   19
```

To obtain a subset of a data frame can be done in similar ways as we have discussed: by index, by row/column names, or by logical values.


```{.r .colsel}
friends["Mina",]
```

```
##      name age child
## Mina Mina  21  TRUE
```


```{.r .colsel}
# The columns of a data frame can be referred to by the names of the columns
friends
```

```
##      name age child
## Mina Mina  21  TRUE
## Ella Ella  27  TRUE
## Anna Anna  26 FALSE
## Cora Cora  32  TRUE
```

```{.r .colsel}
friends$age
```

```
## [1] 21 27 26 32
```

```{.r .colsel}
friends[friends$age > 26,]
```

```
##      name age child
## Ella Ella  27  TRUE
## Cora Cora  32  TRUE
```


```{.r .colsel}
friends[friends$child == "TRUE",]
```

```
##      name age child
## Mina Mina  21  TRUE
## Ella Ella  27  TRUE
## Cora Cora  32  TRUE
```

Function subset() can also be used to get a subset of a data frame.


```{.r .colsel}
# select friends that are older than 26
subset(friends, age > 26)
```

```
##      name age child
## Ella Ella  27  TRUE
## Cora Cora  32  TRUE
```

```{.r .colsel}
# select the information of the ages of friends
subset(friends, select=age)
```

```
##      age
## Mina  21
## Ella  27
## Anna  26
## Cora  32
```


A data frame can be extended.


```{.r .colsel}
# add a column that has the information on the marrital status of friends
friends$married <- c("YES", "YES", "NO", "YES")
friends
```

```
##      name age child married
## Mina Mina  21  TRUE     YES
## Ella Ella  27  TRUE     YES
## Anna Anna  26 FALSE      NO
## Cora Cora  32  TRUE     YES
```

A data frame can also be extended using the functions cbind() and rbind().


```{.r .colsel}
# add a column that has the information on the salaries of friends
cbind(friends, salary=c(4000, 8000, 2000, 6000))
```

```
##      name age child married salary
## Mina Mina  21  TRUE     YES   4000
## Ella Ella  27  TRUE     YES   8000
## Anna Anna  26 FALSE      NO   2000
## Cora Cora  32  TRUE     YES   6000
```

### Lists: a list is an ordered collection of objects, which can be any type of R objects (vectors, matrices, data frames).

A list can be created using the function list().


```{.r .colsel}
my_list <- list(mother="Sophia", father="John", sisters=c("Anna", "Emma"), sister_age=c(5, 10))
my_list
```

```
## $mother
## [1] "Sophia"
## 
## $father
## [1] "John"
## 
## $sisters
## [1] "Anna" "Emma"
## 
## $sister_age
## [1]  5 10
```


```{.r .colsel}
# names of elements in the list
names(my_list)
```

```
## [1] "mother"     "father"     "sisters"    "sister_age"
```


```{.r .colsel}
# number of elements in the list
length(my_list)
```

```
## [1] 4
```

To access elements of a list can be done using its name or index.


```{.r .colsel}
my_list$mother
```

```
## [1] "Sophia"
```


```{.r .colsel}
my_list[["mother"]]
```

```
## [1] "Sophia"
```


```{.r .colsel}
my_list[[1]]
```

```
## [1] "Sophia"
```


```{.r .colsel}
my_list[[3]]
```

```
## [1] "Anna" "Emma"
```


```{.r .colsel}
my_list[[3]][2]
```

```
## [1] "Emma"
```

---

Topic 3. Import and export data in R
====================================================

R base function read.table() is a general funciton that can be used to read a file in table format. The data will be imported as a data frame.


```{.r .colsel}
# To read a local file. If you have downloaded the raw_counts.txt file to your local machine, you may use the following command to read it in, by providing the full path for the file location. The way to specify the full path is the same as taught in the command line session. Here we assume raw_counts.txt is in our current working directory
data <- read.table(file="./Intro2R_files/raw_counts.txt", sep="\t", header=T, stringsAsFactors=F)

# There is a very convenient way to read files from the internet.
data <- read.table(file="https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019_August_UCD_mRNAseq_Workshop/master/intro2R/Intro2R_files/raw_counts.txt", sep="\t", header=T, stringsAsFactors=F)
```

Take a look at the beginning part of the data frame.

```{.r .colsel}
head(data)
```

```
##            C61  C62  C63  C64  C91  C92  C93 C94 I561 I562 I563 I564 I591
## AT1G01010  322  346  256  396  372  506  361 342  638  488  440  479  770
## AT1G01020  149   87  162  144  189  169  147 108  163  141  119  147  182
## AT1G01030   15   32   35   22   24   33   21  35   18    8   54   35   23
## AT1G01040  687  469  568  651  885  978  794 862  799  769  725  715  811
## AT1G01046    1    1    5    4    5    3    0   2    4    3    1    0    2
## AT1G01050 1447 1032 1083 1204 1413 1484 1138 938 1247 1516  984 1044 1374
##           I592 I593 I594 I861 I862 I863 I864 I891 I892 I893 I894
## AT1G01010  430  656  467  143  453  429  206  567  458  520  474
## AT1G01020  156  153  177   43  144  114   50  161  195  157  144
## AT1G01030    8   16   24   42   17   22   39   26   28   39   30
## AT1G01040  567  831  694  345  575  605  404  735  651  725  591
## AT1G01046    8    8    1    0    4    0    3    5    7    0    5
## AT1G01050 1355 1437 1577  412 1338 1051  621 1434 1552 1248 1186
```


Depending on the format of the file, several variants of read.table() are available to make reading a file easier.

read.csv(): for reading "comma separated value" files (.csv).

read.csv2(): variant used in countries that use a comma "," as decimal point and a semicolon ";" as field separators.

read.delim(): for reading "tab separated value" files (".txt"). By default, point(".") is used as decimal point.

read.delim2(): for reading "tab separated value" files (".txt"). By default, comma (",") is used as decimal point.


```{.r .colsel}
# We are going to read a file over the internet by providing the url of the file.
data2 <- read.csv(file="https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019_August_UCD_mRNAseq_Workshop/master/intro2R/Intro2R_files/raw_counts.txt", stringsAsFactors=F)

# To look at the file:
head(data2)
```

```
##                                         C61.C62.C63.C64.C91.C92.C93.C94.I561.I562.I563.I564.I591.I592.I593.I594.I861.I862.I863.I864.I891.I892.I893.I894
## 1                     AT1G01010\t322\t346\t256\t396\t372\t506\t361\t342\t638\t488\t440\t479\t770\t430\t656\t467\t143\t453\t429\t206\t567\t458\t520\t474
## 2                        AT1G01020\t149\t87\t162\t144\t189\t169\t147\t108\t163\t141\t119\t147\t182\t156\t153\t177\t43\t144\t114\t50\t161\t195\t157\t144
## 3                                               AT1G01030\t15\t32\t35\t22\t24\t33\t21\t35\t18\t8\t54\t35\t23\t8\t16\t24\t42\t17\t22\t39\t26\t28\t39\t30
## 4                     AT1G01040\t687\t469\t568\t651\t885\t978\t794\t862\t799\t769\t725\t715\t811\t567\t831\t694\t345\t575\t605\t404\t735\t651\t725\t591
## 5                                                                     AT1G01046\t1\t1\t5\t4\t5\t3\t0\t2\t4\t3\t1\t0\t2\t8\t8\t1\t0\t4\t0\t3\t5\t7\t0\t5
## 6 AT1G01050\t1447\t1032\t1083\t1204\t1413\t1484\t1138\t938\t1247\t1516\t984\t1044\t1374\t1355\t1437\t1577\t412\t1338\t1051\t621\t1434\t1552\t1248\t1186
```



R base function write.table() can be used to export data to a file.


```{.r .colsel}
# To write to a file called "output.txt" in your current working directory.
write.table(data2[1:20,], file="output.txt", sep="\t", quote=F, row.names=T, col.names=T)
```

It is also possible to export data to a csv file.

write.csv()

write.csv2()


\newpage

---

Topic 4. Functions in R
====================================================
### Invoking a function by its name, followed by the parenthesis and zero or more arguments.


```{.r .colsel}
# to find out the current working directory
getwd()
```

```
## [1] "/Users/mattsettles/projects/src/github.com-ucdavis-bioinformatics-training/2019_August_UCD_mRNAseq_Workshop/intro2R"
```

```{.r .colsel}
# to set a different working directory, use setwd
#setwd("/Users/jli/Desktop")

# to list all variables in the environment
ls()
```

```
##  [1] "a"             "age"           "b"             "col1"         
##  [5] "col2"          "col3"          "colFmt"        "data"         
##  [9] "data2"         "friend_ages"   "friend_groups" "friend_names" 
## [13] "friends"       "has_child"     "hello"         "is_adult"     
## [17] "is_female"     "is_male"       "my_data"       "my_friends"   
## [21] "my_list"       "my_matrix"     "title"
```

```{.r .colsel}
# to create a vector from 2 to 3, using increment of 0.1
seq(2, 3, by=0.1)
```

```
##  [1] 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0
```

```{.r .colsel}
# to create a vector with repeated elements
rep(1:3, times=3)
```

```
## [1] 1 2 3 1 2 3 1 2 3
```

```{.r .colsel}
rep(1:3, each=3)
```

```
## [1] 1 1 1 2 2 2 3 3 3
```

```{.r .colsel}
# to get help information on a function in R: ?function.name
?seq
?sort
?rep
```

#### <font color='red'>One useful function to find out information on a variable: str().</font>  



```{.r .colsel}
str(data2)
```

```
## 'data.frame':	33602 obs. of  1 variable:
##  $ C61.C62.C63.C64.C91.C92.C93.C94.I561.I562.I563.I564.I591.I592.I593.I594.I861.I862.I863.I864.I891.I892.I893.I894: chr  "AT1G01010\t322\t346\t256\t396\t372\t506\t361\t342\t638\t488\t440\t479\t770\t430\t656\t467\t143\t453\t429\t206\t"| __truncated__ "AT1G01020\t149\t87\t162\t144\t189\t169\t147\t108\t163\t141\t119\t147\t182\t156\t153\t177\t43\t144\t114\t50\t161\t195\t157\t144" "AT1G01030\t15\t32\t35\t22\t24\t33\t21\t35\t18\t8\t54\t35\t23\t8\t16\t24\t42\t17\t22\t39\t26\t28\t39\t30" "AT1G01040\t687\t469\t568\t651\t885\t978\t794\t862\t799\t769\t725\t715\t811\t567\t831\t694\t345\t575\t605\t404\t"| __truncated__ ...
```


### Two special functions: lapply() and sapply()
#### lapply() is to apply a given function to every element of a vector and obtain a list as results.
#### The difference between lapply() and apply() is that lapply() can be applied on objects like dataframes, lists or vectors. Function apply() only works on an array of dimension 2 or a matrix.

##### To check the syntax of using lapply():

```{.r .colsel}
#?lapply()

data <- as.data.frame(matrix(rnorm(49), ncol=7), stringsAsFactors=F)
dim(data)
```

```
## [1] 7 7
```

```{.r .colsel}
lapply(1:dim(data)[1], function(x){sum(data[x,])})
```

```
## [[1]]
## [1] 2.696972
## 
## [[2]]
## [1] 0.1402984
## 
## [[3]]
## [1] -3.111006
## 
## [[4]]
## [1] 0.1516168
## 
## [[5]]
## [1] 0.01941972
## 
## [[6]]
## [1] -2.913826
## 
## [[7]]
## [1] -1.392908
```

```{.r .colsel}
apply(data, MARGIN=1, sum)
```

```
## [1]  2.69697170  0.14029844 -3.11100645  0.15161683  0.01941972 -2.91382587
## [7] -1.39290826
```

```{.r .colsel}
lapply(1:dim(data)[1], function(x){log10(sum(data[x,]))})
```

```
## Warning in FUN(X[[i]], ...): NaNs produced

## Warning in FUN(X[[i]], ...): NaNs produced

## Warning in FUN(X[[i]], ...): NaNs produced
```

```
## [[1]]
## [1] 0.4308764
## 
## [[2]]
## [1] -0.8529472
## 
## [[3]]
## [1] NaN
## 
## [[4]]
## [1] -0.8192526
## 
## [[5]]
## [1] -1.711757
## 
## [[6]]
## [1] NaN
## 
## [[7]]
## [1] NaN
```

#### The function sapply() works like function lapply(), but tries to simplify the output to the most elementary data structure that is possible. As a matter of fact, sapply() is a "wrapper" function for lapply(). By default, it returns a vector.


```{.r .colsel}
# To check the syntax of using sapply():
#?sapply()

sapply(1:dim(data)[1], function(x){log10(sum(data[x,]))})
```

```
## Warning in FUN(X[[i]], ...): NaNs produced

## Warning in FUN(X[[i]], ...): NaNs produced

## Warning in FUN(X[[i]], ...): NaNs produced
```

```
## [1]  0.4308764 -0.8529472        NaN -0.8192526 -1.7117570        NaN
## [7]        NaN
```

#### If the "simplify" parameter is turned off, sapply() will produced exactly the same results as lapply(), in the form of a list. By default, "simplify" is turned on.

```{.r .colsel}
sapply(1:dim(data)[1], function(x){log10(sum(data[x,]))}, simplify=FALSE)
```

```
## Warning in FUN(X[[i]], ...): NaNs produced

## Warning in FUN(X[[i]], ...): NaNs produced

## Warning in FUN(X[[i]], ...): NaNs produced
```

```
## [[1]]
## [1] 0.4308764
## 
## [[2]]
## [1] -0.8529472
## 
## [[3]]
## [1] NaN
## 
## [[4]]
## [1] -0.8192526
## 
## [[5]]
## [1] -1.711757
## 
## [[6]]
## [1] NaN
## 
## [[7]]
## [1] NaN
```

---

Topic 5. Basic statistics in R
====================================================

<table class="table table-striped table-hover table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:center;"> Description </th>
   <th style="text-align:center;"> R_function </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> Mean </td>
   <td style="text-align:center;"> mean() </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Standard deviation </td>
   <td style="text-align:center;"> sd() </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Variance </td>
   <td style="text-align:center;"> var() </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Minimum </td>
   <td style="text-align:center;"> min() </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Maximum </td>
   <td style="text-align:center;"> max() </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Median </td>
   <td style="text-align:center;"> median() </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Range of values: minimum and maximum </td>
   <td style="text-align:center;"> range() </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Sample quantiles </td>
   <td style="text-align:center;"> quantile() </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Generic function </td>
   <td style="text-align:center;"> summary() </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Interquartile range </td>
   <td style="text-align:center;"> IQR() </td>
  </tr>
</tbody>
</table>

Calculate the mean expression for each sample.


```{.r .colsel}
apply(data, 2, mean)
```

```
##          V1          V2          V3          V4          V5          V6 
## -0.64883337  0.41640650 -0.26557523 -0.06964068 -0.27446079  0.27862684 
##          V7 
## -0.06644240
```

Calculate the range of expression for each sample.


```{.r .colsel}
apply(data, 2, range)
```

```
##              V1         V2         V3        V4         V5        V6
## [1,] -1.8921134 -0.5167689 -0.9437552 -1.037911 -1.1690234 -0.680934
## [2,]  0.9787187  1.3551789  1.0341524  1.365358  0.7155543  1.286370
##             V7
## [1,] -1.027381
## [2,]  1.845798
```

Calculate the quantiles of each samples.


```{.r .colsel}
apply(data, 2, quantile)
```

```
##              V1         V2          V3          V4          V5         V6
## 0%   -1.8921134 -0.5167689 -0.94375525 -1.03791140 -1.16902338 -0.6809340
## 25%  -1.4150087 -0.0871726 -0.84136018 -0.57622027 -0.73619527 -0.2789512
## 50%  -1.1195325  0.2051671 -0.34303930 -0.08251126 -0.07999668  0.2095695
## 75%   0.1605554  1.0228068  0.03816795  0.21001020  0.04231540  0.8466422
## 100%  0.9787187  1.3551789  1.03415239  1.36535801  0.71555429  1.2863702
##              V7
## 0%   -1.0273810
## 25%  -0.5771931
## 50%  -0.4228686
## 75%   0.1468706
## 100%  1.8457978
```

---

Topic 6. Simple data visulization in R
====================================================

Scatter plot and line plot can be produced using the function plot().


```{.r .colsel}
x <- c(1:50)
y <- 1 + sqrt(x)/2
plot(x,y)
```

![](Intro2R_files/figure-html/unnamed-chunk-68-1.png)<!-- -->

```{.r .colsel}
plot(x,y, type="l")
```

![](Intro2R_files/figure-html/unnamed-chunk-68-2.png)<!-- -->

```{.r .colsel}
# plot both the points and lines
## first plot points
plot(x,y)
lines(x,y, type="l")
```

![](Intro2R_files/figure-html/unnamed-chunk-68-3.png)<!-- -->

```{.r .colsel}
## lines() can only be used to add information to a graph, while it cannot produce a graph on its own.
```


boxplot() can be used to summarize data.


```{.r .colsel}
boxplot(data, xlab="Sample ID", ylab="Raw Counts")
```

![](Intro2R_files/figure-html/unnamed-chunk-69-1.png)<!-- -->

add more details to the plot.


```{.r .colsel}
boxplot(data, xlab="Sample ID", ylab="Raw Counts", main="Expression levels", col="blue", border="black")
```

![](Intro2R_files/figure-html/unnamed-chunk-70-1.png)<!-- -->



```{.r .colsel}
x <- rnorm(1000)
boxplot(x)
```

![](Intro2R_files/figure-html/unnamed-chunk-71-1.png)<!-- -->

hist() can be used to create histograms of data.

```{.r .colsel}
hist(x)
```

![](Intro2R_files/figure-html/unnamed-chunk-72-1.png)<!-- -->

```{.r .colsel}
# use user defined break points
hist(x, breaks=seq(range(x)[1]-1, range(x)[2]+1, by=0.5))
```

![](Intro2R_files/figure-html/unnamed-chunk-72-2.png)<!-- -->


```{.r .colsel}
# clear plotting device/area
dev.off()
```

```
## null device 
##           1
```

---

Topic 7. Install packages in R
====================================================

Starting from Bioconductor version 3.8, the installation of packages is recommended to use BiocManager.


```{.r .colsel}
if (!requireNamespace("BiocManager"))
	install.packages("BiocManager")
## install core packages
BiocManager::install()
## install specific packages
BiocManager::install(c("ggplot2", "ShortRead"))
```

* Bioconductor has a repository and release schedule that differ from R (Bioconductor has a ‘devel’ branch to which new packages and updates are introduced, and a stable ‘release’ branch emitted once every 6 months to which bug fixes but not new features are introduced). This mismatch causes that the version detected by install.packages() is sometimes not the most recent 'release'.

* A consequence of the distince 'devel' branch is that install.packages() sometimes points only to the 'release' repository, while users might want to have access to the leading-edge features in the develop version.

* An indirect consequence of Bioconductor's structured release is that packages generally have more extensive dependences with one another.


##### <font color='red'>It is always recommended to update to the most current version of R and Bioconductor. If it is not possible and R < 3.5.0, please use the legacy approach to install Bioconductor packages</font>   


```{.r .colsel}
source("http://bioconductor.org/biocLite.R")
## install core packages
biocLite()
## install specific packages
biocLite("RCircos")
biocLite(c("IdeoViz", "devtools"))
```

The R function install.packages() can be used to install packages that are not part of Bioconductor.


```{.r .colsel}
install.packages("ggplot2", repos="http://cran.us.r-project.org")
```

Install from source of github.

```{.r .colsel}
library(devtools)
#install_github("stephenturner/qqman")
```


---

Topic 8. Save data in R session
====================================================

#### To save history in R session


```{.r .colsel}
#savehistory(file="March27.history")

#loadhistory(file="March27.history")
```

#### To save objects in R session


```{.r .colsel}
save(list=c("x", "data"), file="March27.RData")

#load("March27.RData")
```

Challenge
====================================================

Read in the reference fasta file for PhiX that you have downloaded this morning and find out the length of PhiX genome and GC content.

<font color='blue'>Hint: please look into bioconductor package ShortRead.</font>
<font color='blue'>The PhiX genome can be downloaded by following the instruction in Intro to Command-Line session. The link for the data is at ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz</font>
