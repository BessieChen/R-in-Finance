# display documentation on function "getwd"
help(getwd)
?getwd  # equivalent to "help(getwd)"
help.start()  # open the hypertext documentation
# "<-" and "=" are valid assignment operators
my_var <- 3

# typing a symbol or expression evaluates it
my_var

# text in quotes is interpreted as a string
my_var <- "Hello World!"

# typing a symbol or expression evaluates it
my_var

my_var  # text after hash is treated as comment
getwd()  # get cwd
setwd("C:/Develop/R")  # set cwd
getwd()  # get cwd
Sys.time()  # get date and time

Sys.Date()  # get date only
rm(list=ls())
setwd("C:/Develop/R/lecture_slides/data")
var1 <- 3  # define new object
ls()  # list all objects in workspace
# list objects starting with "v"
ls(pattern=glob2rx("v*"))
# remove all objects starting with "v"
rm(list=ls(pattern=glob2rx("v*")))
save.image()  # save workspace to file .RData in cwd
rm(var1)  # remove object
ls()  # list objects
load(".RData")
ls()  # list objects
var2 <- 5  # define another object
save(var1, var2,  # save selected objects
     file="C:/Develop/R/lecture_slides/data/my_data.RData")
rm(list=ls())  # remove all objects
ls()  # list objects
load_ed <- load(file="C:/Develop/R/lecture_slides/data/my_data.RData")
load_ed
ls()  # list objects
  q()  # quit R session
history(5)  # display last 5 commands
savehistory(file="myfile")  # default is ".Rhistory"
loadhistory(file="myfile")  # default is ".Rhistory"
sessionInfo()  # get R version and other session info
Sys.getenv()[5:7]  # list some environment variables

Sys.getenv("Home")  # get R user HOME directory

Sys.setenv(Home="C:/Develop/data")  # set HOME directory

Sys.getenv("Home")  # get user HOME directory

Sys.getenv("R_home")  # get R_HOME directory

R.home()  # get R_HOME directory

R.home("etc")  # get "etc" sub-directory of R_HOME
# ?options  # long list of global options
# interpret strings as characters, not factors
getOption("stringsAsFactors")  # display option
options("stringsAsFactors")  # display option
options(stringsAsFactors=FALSE)  # set option
# number of digits printed for numeric values
options(digits=3)
# control exponential scientific notation of print method
# positive "scipen" values bias towards fixed notation
# negative "scipen" values bias towards scientific notation
options(scipen=100)
# maximum number of items printed to console
options(max.print=30)
# warning levels options
# negative - warnings are ignored
options(warn=-1)
# zero - warnings are stored and printed after top-level function has completed
options(warn=0)
# one - warnings are printed as they occur
options(warn=1)
# two or larger - warnings are turned into errors
options(warn=2)
# save all options in variable
op_tions <- options()
# restore all options from variable
options(op_tions)
# single numbers are vectors of length 1
1
# character strings are vectors of length 1
"a"
# strings without quotes are variable names
a  # variable "a" doesn't exist
# list elements can have different mode
list(aa=c('a', 'b'), bb=1:5)
data.frame(aa=c('a', 'b'), bb=1:66)
is.atomic(data.frame(aa=c('a', 'b'), bb=1:2))
is.recursive(data.frame(aa=c('a', 'b'), bb=1:2))
my_var <- "hello"
c(typeof(my_var), mode(my_var), class(my_var))

my_var <- 1:5
c(typeof(my_var), mode(my_var), class(my_var))

my_var <- runif(5)
c(typeof(my_var), mode(my_var), class(my_var))

my_var <- matrix(1:10, 2, 5)
c(typeof(my_var), mode(my_var), class(my_var))

my_var <- matrix(runif(10), 2, 5)
c(typeof(my_var), mode(my_var), class(my_var))

my_var <- list(aa=c('a', 'b'), bb=1:5)
c(typeof(my_var), mode(my_var), class(my_var))

my_var <- data.frame(aa=c('a', 'b'), bb=1:2)
c(typeof(my_var), mode(my_var), class(my_var))
# a simple vector has no attributes
attributes(5:10)
my_var <- c(pi=pi, euler=exp(1), gamma=-digamma(1))
# named vector has "names" attribute
attributes(my_var)
my_var <- 1:10
is.vector(my_var)  # is the object a vector?
attributes(my_var) <- list(my_attr="foo")
my_var
is.vector(my_var)  # is the object a vector?
my_var <- 0
attributes(my_var) <- list(class="Date")
my_var  # "Date" object
structure(0, class="numeric")  # "Date" object
structure(0, class="Date")  # "Date" object

my_var <- matrix(runif(10), 2, 5)
class(my_var)  # has implicit class
# but no explicit "class" attribute
attributes(my_var)
c(typeof(my_var), mode(my_var), class(my_var))
# assign explicit "class" attribut
 
class(my_var) <- c("my_class","my_kris")
class(my_var)  # has explicit "class"
# has explicit "class" attribute
attributes(my_var)
is.matrix(my_var)  # is the object a matrix?
is.vector(my_var)  # is the object a vector?
attributes(unclass(my_var))
# integer implicit class derived from type
my_var <- vector(mode="integer", length=10)
c(typeof(my_var), mode(my_var), class(my_var))
# numeric implicit class derived from mode
my_var <- vector(mode="numeric", length=10)
c(typeof(my_var), mode(my_var), class(my_var))
# adding dim attribute changes implicit class to matrix
dim(my_var) <- c(5, 2)
c(typeof(my_var), mode(my_var), class(my_var))
# data frames have implicit dim attribute
my_var <- data.frame(aa=c('a', 'b'), bb=1:28)
c(typeof(my_var), mode(my_var), class(my_var))
attributes(my_var)
dim(my_var)
my_var <- 1:5
c(typeof(my_var), mode(my_var), class(my_var))
mode(my_var) <- "character"  # coerce to "character"
my_var
c(typeof(my_var), mode(my_var), class(my_var))
# explicitly coerce to "character"
my_var <- as.character(1:5)
c(typeof(my_var), mode(my_var), class(my_var))
mat_rix <- matrix(1:10, 2, 5)  # create matrix
# explicitly coerce to "character"
mat_rix <- as.character(mat_rix)
c(typeof(mat_rix), mode(mat_rix), class(mat_rix))
# coercion converted matrix to vector
c(is.matrix(mat_rix), is.vector(mat_rix))
as.logical(0:3)  # explicit coercion to "logical"
as.numeric(c(FALSE, TRUE, TRUE, TRUE))
c(1:3, 'a')  # implicit coercion to "character"
# explicit coercion to "numeric"
as.numeric(c(1:3, 'a'))
"Hello World!"  # type some text
# hello is a variable name, because it's not in quotes
hello  # R interprets "hello" as a variable name
is.vector(1)  # single number is a vector
is.vector("a")  # string is a vector
4:8  # create a vector
# create vector using c() combine function
c(1, 2, 3, 4, 5)
# create vector using c() combine function
c('a', 'b', 'c')
# create vector using c() combine function
c(1, 'b', 'c')
str_var <- "Some string"
str_var
str_var[1]
str_var[2]

length(str_var)  # length of vector
nchar(str_var)  # length of string

# concatenate and echo to console
cat("Hello", "World!")
cat("Enter\ttab")
cat("Enter\nnewline")
cat("Enter\\backslash")
str_var1 <- "Hello"  # define a character string
str_var2 <- "World!"  # define a character string
paste(str_var1, str_var2, sep=' ')  # concatenate and return value
cat(str_var1, str_var2)  # concatenate and echo to console
paste('a', 1:4, sep='-')  # convert, recycle and concatenate
paste('a', 1:4, collapse='-')
paste(c("a1", "a2", "a3"), collapse="+") 
paste(c("a1", "a2", "a3"), sep="+") # collapse vector to string
paste(list("a1", "a2", "a3"), collapse="+")
paste(list("a1", "a2", "a3"), sep="+")
paste("Today is", Sys.time())  # coerce and concatenate strings
paste("Today is", format(Sys.time(), "%B-%d-%Y"))
strsplit("Hello World", split='r')  # split string
strsplit("Hello.World", split='[.]')  # split string
strsplit("Hello.World", split='.', fixed=TRUE)  # split string
substring("Hello World", 3, 6)  # extract characters from 3 to 6
gsub("is", "XX", "is this gratis?")  # replace "is" with "XX"

grep("b", c("abc", "xyz", "cba d", "bbb"))  # get indexes

grep("b", c("abc", "xyz", "cba d", "bbb"), value=TRUE)  # get values

glob2rx("abc.*")  # convert globs into regex
glob2rx("*.doc")
is.vector(1)  # single number is a vector
is.vector("a")  # string is a vector
vec_tor <- c(8, 6, 5, 7)  # create vector
vec_tor
vec_tor[2]  # extract second element
# extract all elements, except the second element
vec_tor[-2]
# create Boolean vector
c(FALSE, TRUE, TRUE)
# extract second and third elements
vec_tor[c(FALSE, TRUE, TRUE)]
letters[5:10]  # vector of letters
c('a', letters[5:10])  # combine two vectors of letters
0:10  # vector of integers from 0 to 10
vector()  # create empty vector
vector(mode="numeric", length=10)  # numeric vector of zeros
seq(10)  # sequence from 1 to 10
seq(along=(-5:5))  # instead of 1:length(obj)
seq_along(c("a", "b", "c"))  # instead of 1:length(obj)
seq(from=0, to=1, len=11)  # decimals from 0 to 1.0
seq(from=0, to=1, by=0.1)  # decimals from 0 to 1.0
seq(-2,2, len=11)  # 10 numbers from -2 to 2
rep(100, times=5)  # replicate a number
character(5)  # create empty character vector
numeric(5)  # create empty numeric vector
numeric(0)  # create zero-length vector
2*4:8  # multiply a vector
2*(4:8)  # multiply a vector
4:8/2  # divide a vector
(0:10)/10  # divide vector - decimals from 0 to 1.0
vec_tor <- c(8, 6, 5, 7)  # create vector
vec_tor
# Boolean vector TRUE if element is equal to second one
vec_tor == vec_tor[2]
# Boolean vector TRUE for elements greater than six
vec_tor > 6
2*vec_tor  # multiply all elements by 2
vec_tor^2  # square all elements
c(11, 5:10)  # combine two vectors
c(vec_tor, 2.0)  # append number to vector
vec_tor <- # create named vector
  c(pi_const=pi, euler=exp(1), gamma=-digamma(1))
vec_tor
names(vec_tor)  # get names of elements
vec_tor["euler"]  # get element named "euler"
names(vec_tor) <- c("pie","eulery","gammy")  # rename elements
vec_tor
unname(vec_tor)  # remove names attribute
letters[5:10]  # vector of letters
c('a', letters[5:10])  # combine two vectors of letters
# create named vector
structure(sample(1:5), names=paste0("el", 1:5))
vec_tor  # named vector
# extract second element
vec_tor[2]
# extract all elements, except the second element
vec_tor[-2]
# extract zero elements - returns zero-length vector
vec_tor[0]
# extract second and third elements
vec_tor[c(FALSE, TRUE, TRUE)]
# extract elements using their names
vec_tor["eulery"]
# extract elements using their names
vec_tor[c("pie", "gammy")]
# subset whole vector
vec_tor[] <- 0
vec_tor <- runif(5)
vec_tor
vec_tor > 0.5  # Boolean vector
# Boolean vector of elements equal to the second one
vec_tor == vec_tor[2]
# extract all elements equal to the second one
vec_tor[vec_tor == vec_tor[2]]
vec_tor < 1  # Boolean vector of elements less than one
# extract all elements greater than one
vec_tor[vec_tor > 1]
vec_tor[vec_tor > 0.5]  # filter elements > 0.5
which(vec_tor > 0.5)  # index of elements > 0.5
mat_rix <- matrix(5:10, nrow=2, ncol=3)  # create a matrix
mat_rix  # by default matrices are constructed column-wise
# create a matrix row-wise
matrix(5:10, nrow=2, byrow=TRUE)
mat_rix[2, 3]  # extract third element from second row
mat_rix[2, ]  # extract second row
mat_rix[, 3]  # extract third column
mat_rix[, c(1,3)]  # extract first and third column
mat_rix[, -2]  # remove second column
# subset whole matrix
mat_rix[] <- 0
# get the number of rows or columns
nrow(vec_tor); ncol(vec_tor)
NROW(vec_tor); NCOL(vec_tor)
nrow(mat_rix); ncol(mat_rix)
NROW(mat_rix); NCOL(mat_rix)
attributes(mat_rix)  # get matrix attributes
dim(mat_rix)  # get dimension attribute
class(mat_rix)  # get class attribute
rownames(mat_rix) <- c("row1", "row2")  # rownames attribute
colnames(mat_rix) <- c("col1", "col2", "col3")  # colnames attribute
mat_rix
mat_rix["row2", "col3"]  # third element from second row
names(mat_rix)  # get the names attribute
dimnames(mat_rix)  # get dimnames attribute
attributes(mat_rix)  # get matrix attributes
mat_rix  # matrix with column names
mat_rix[1, ]  # subset rows by index
mat_rix[, "col1"]  # subset columns by name
mat_rix[, c(TRUE, FALSE, TRUE)]  # subset columns Boolean vector
mat_rix[1, ]  # subsetting can produce a vector!
class(mat_rix); class(mat_rix[1, ])
is.matrix(mat_rix[1, ]); is.vector(mat_rix[1, ])
mat_rix <- matrix(5:10, nrow=2, byrow=TRUE)
mat_rix[1, , drop=FALSE] 
mat_rix[1, , drop=TRUE] # drop=FALSE preserves matrix
class(mat_rix[1, , drop=TRUE])
is.matrix(mat_rix[1, , drop=TRUE]); is.vector(mat_rix[1, , drop=TRUE])
rm(list=ls())
TRUE | FALSE
TRUE | NA
vec_tor1 <- c(2, 4, 6)
vec_tor1 < 5  # element-wise comparison
(vec_tor1 < 5) & (vec_tor1 > 3)
vec_tor1[(vec_tor1 < 5) & (vec_tor1 > 3)]
vec_tor2 <- c(-10, 0, 10)
vec_tor1 < vec_tor2
vec_tor2 <- c(-10, 0, 10,11)
vec_tor1 < vec_tor2
c(FALSE, TRUE, FALSE) & c(TRUE, TRUE, FALSE)
c(FALSE, TRUE, FALSE) | c(TRUE, TRUE, FALSE)
rm(list=ls())
c(FALSE, TRUE, FALSE) && c(TRUE, TRUE, FALSE)
c(FALSE, TRUE, FALSE) || c(TRUE, TRUE, FALSE)
echo_true <- function() {cat("echo_true\t"); TRUE}
echo_false <- function() {cat("echo_false\t"); FALSE}
echo_true() | echo_false()
echo_true() || echo_false()
echo_false() || echo_true()
# echo_false() isn't evaluated at all!
vec_tor <- c(2, 4, 6)
# works (does nothing) using '&&'
if (is.matrix(vec_tor) && (vec_tor[2, 3] > 0)) {
  vec_tor[2, 3] <- 1
}
# no short-circuit so fails (produces an error)
if (is.matrix(vec_tor) & (vec_tor[2, 3] > 0)) {
  vec_tor[2, 3] <- 1
}
?Arithmetic
4.7 * 0.5  # multiplication
4.7 / 0.5  # division
# exponentiation
2**3
2^3
num_var <- 2
num_var==2
identical(num_var, 2)

identical(num_var, NULL)
# this doesn't work:
num_var==NULL
is.null(num_var)

vec_tor <- c(2, 4, 6)
vec_tor==2
identical(vec_tor, 2)

# num_ber is equal to "1.0" within machine precision
num_ber <- 1.0 + 2*sqrt(.Machine$double.eps)
all.equal(num_ber, 1.0)

# info machine precision of computer R is running on
# ?.Machine
# machine precision
.Machine$double.eps
vec_tor <- sample(1:6, 21, replace=TRUE)
mat_rix <- matrix(vec_tor, ncol=3)
vec_tor
which(vec_tor == 5)
# equivalent but slower than above
(1:length(vec_tor))[vec_tor == 5]
which(vec_tor > 5)
# find indices of TRUE elements of Boolean matrix
which((mat_rix == 5)|(mat_rix == 6),
      arr.ind=TRUE)
# equivalent but slower than above
arrayInd(which((mat_rix == 5)|(mat_rix == 6)),dim(mat_rix), dimnames(mat_rix)
 )
which.max(vec_tor)
?which.max
x <- c(1:4, 0:5, 11,5)
which.min(x)
which.max(x)
# equivalent but slower than above
which(vec_tor == max(vec_tor))
which.min(vec_tor)
match(5, vec_tor)
match(53, vec_tor)
# more general but slower than above
which(vec_tor == 5)
match(-5, vec_tor)
5 %in% vec_tor
vec_tor %in% 5
# equivalent to above
match(5, vec_tor, nomatch=0) > 0
-5 %in% vec_tor
c(5, -5) %in% vec_tor
# equivalent to "5 %in% vec_tor"
any(vec_tor == 5)
# equivalent to "-5 %in% vec_tor"
any(vec_tor == (-5))
if (any(vec_tor < 0))
  cat("vector contains negative values\n")
# partial matching of strings
pmatch("med", c("mean", "median", "mode"))
str(findInterval)
# get index of the element of "vec" that matches 5
findInterval(x=5, vec=c(3, 5, 7))
match(5, c(3, 5, 7))
# no exact match
findInterval(x=6, vec=c(3, 5, 7))
match(6, c(3, 5, 7))
# indices of "vec" that match elements of "x"
findInterval(x=1:8, vec=c(3, 5, 7))
# return only indices of inside intervals
findInterval(x=1:8, vec=c(3, 5, 7),
       all.inside=TRUE)
# make rightmost interval inclusive
findInterval(x=1:8, vec=c(3, 5, 7),
       rightmost.closed=TRUE)
num_var1 <- 3  # "<-" and "=" are valid assignment operators
num_var1
num_var1 = 3
num_var1
2<-3  # "<" operator confused with "<-"
2 < -3  # add space or brackets to avoid confusion
# "=" assignment within argument list
median(a=1:10)
a  # x doesn't exist outside the function
# "<-" assignment within argument list
median(a <- 1:10)
a  # x exists outside the function
rm(list=ls())
# expressions enclosed in parenthesis are less ambiguous
-2:5
(-2):5
-(2:5)
# expressions enclosed in parenthesis are less ambiguous
-2*3+5
-2*(3+5)

# expressions can be separated by semicolons or by lines
{1+2;  1:5;2*3}
# or
{1+2
2*3
1:5}

mat_rix <- matrix(nr=3, nc=4)
mat_rix <- 0
# subset whole matrix
mat_rix[] <- 0

# parenthesis and braces require a little additional processing time
library(microbenchmark)
summary(microbenchmark(
  ba_se=sqrt(rnorm(10000)^2),
  pa_ren=sqrt(((((rnorm(10000)^2))))),
  bra_ce=sqrt({{{{rnorm(10000)^2}}}}),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary
rm(list=ls())
num_var1 <- 1

if (num_var1) {  # numeric zero is FALSE, all other numbers are TRUE
  num_var2 <- 4
} else if (num_var1 == 0) {  # 'else if' together on same line
  num_var2 <- 0
} else {  # 'else' together with curly braces
  num_var2 <- -4
}  # end if

num_var2
switch("a", a="aaahh", b="bee", c="see", d=2,
       "else this")
switch("c", a="aaahh", b="bee", c="see", d=2,
       "else this")
switch(3, a="aaahh", b="bee", c="see", d=2,
       "else this")
switch("cc", a="aaahh", b="bee", c="see", d=2,
       "else this")
ccc <- c("b","QQ","a","A","bb")
# note: cat() produces no output for NULL
for(ch in ccc)
  cat(ch,":", switch(EXPR = ch, a = 1, b = 2:3), "\n")
for(ch in ccc)
  cat(ch,":", switch(EXPR = ch, a =, A = 1, b = 2:3, "Otherwise: last"),"\n")
centre <- function(x, type) {
  switch(type,
         mean = mean(x),
         median = median(x),
         trimmed = mean(x, trim = .1))
}
x <- rcauchy(10)
centre(x, "mean")
centre(x, "median")
centre(x, "trimmed")

?switch
# measure of central tendency
centra_lity <- function(in_put,
    meth_od=c("mean", "mean_narm", "median")) {
# validate "meth_od" argument
  #meth_od <- match.arg(meth_od)
  switch(meth_od,
 mean=mean(in_put),
 mean_narm=mean(in_put, na.rm=TRUE),
 median=median(in_put))
}  # end centra_lity
#?match.arg
my_var <- rnorm(100, mean=2)
centra_lity(my_var, "mean")
centra_lity(my_var, "mean_narm")
centra_lity(my_var, "median")
for (in_dex in vec_tor) {ex_pressions}
rm(list=ls())
color_list <- list("red", "white", "blue")
# loop over list
for (some_color in color_list) {
  print(some_color)
}  # end for
# loop over vector
for (in_dex in 1:3) {
  print(color_list[[in_dex]])
}  # end for

# while loops require initialization
in_dex <- 1
# while loop
while (in_dex < 4) {
  print(color_list[[in_dex]])
  in_dex <- in_dex + 1
}  # end while
rm(list=ls())
# loop over a vector and overwrite it
vec_tor <- integer(7)
for (i in 1:7) {
  cat("Changing element:", i, "\n")
  vec_tor[i] <- i^2
}  # end for
# equivalent way (without cat side effect)
for (i in seq_along(vec_tor)) 
  vec_tor[i] <- i^2

# sapply() loop returns vector of values
vec_tor <- sapply(seq_along(vec_tor), 
          function(x) (x^2))
rm(list=ls())
# fib_seq <- numeric()  # zero length numeric vector
# pre-allocate vector instead of "growing" it
fib_seq <- numeric(10)
fib_seq[1] <- 0  # initialize
fib_seq[2] <- 1  # initialize
for (i in 3:10) {  # perform recurrence loop
  fib_seq[i] <- fib_seq[i-1] + fib_seq[i-2]
}  # end for
fib_seq
# allocate character vector
character()
character(5)
is.character(character(5))
# allocate integer vector
integer()
integer(5)
is.integer(integer(5))
is.numeric(integer(5))
# allocate numeric vector
numeric()
numeric(5)
is.integer(numeric(5))
is.numeric(numeric(5))
# allocate Boolean vector
vector()
vector(length=5)
# allocate numeric vector
vector(length=5, mode="numeric")
is.null(vector())
# allocate Boolean matrix
matrix()
is.null(matrix())
# allocate integer matrix
matrix(NA_integer_, nrow=3, ncol=2)
is.integer(matrix(NA_integer_, nrow=3, ncol=2))
# allocate numeric matrix
matrix(NA_real_, nrow=3, ncol=2)
is.numeric(matrix(NA_real_, nrow=3, ncol=2))
vec_tor <- sample(1:9)
vec_tor
vec_tor < 5  # element-wise comparison
vec_tor == 5  # element-wise comparison
mat_rix <- matrix(vec_tor, ncol=3)
mat_rix
mat_rix < 5  # element-wise comparison
mat_rix == 5  # element-wise comparison
mat_rix <- 1:6  # create a vector
class(mat_rix)  # get its class
# is it vector or matrix?
c(is.vector(mat_rix), is.matrix(mat_rix))
structure(mat_rix, dim=c(2, 3))  # matrix object
# adding dimension attribute coerces into matrix
dim(mat_rix) <- c(2, 3)
class(mat_rix)  # get its class
# is it vector or matrix?
c(is.vector(mat_rix), is.matrix(mat_rix))
# assign dimnames attribute
dimnames(mat_rix) <- list(rows=c("row1", "row2"),
                  columns=c("col1", "col2", "col3"))
dimnames(mat_rix) <- list(columns=c("row1", "row2"),
                          rows=c("col1", "col2", "col3"))
mat_rix
mat_rix <- matrix(1:10, 2, 5)  # create matrix
mat_rix
# as.numeric strips dim attribute from matrix
as.numeric(mat_rix)
# explicitly coerce to "character"
mat_rix <- as.character(mat_rix)
c(typeof(mat_rix), mode(mat_rix), class(mat_rix))
# coercion converted matrix to vector
c(is.matrix(mat_rix), is.vector(mat_rix))
vec_tor1 <- 1:3  # define vector
vec_tor2 <- 6:4  # define vector
# bind vectors into columns
cbind(vec_tor1, vec_tor2)
# bind vectors into rows
rbind(vec_tor1, vec_tor2)
# extend to four elements
vec_tor2 <- c(vec_tor2, 7)
# recycling rule applied
cbind(vec_tor1, vec_tor2)
# another example of recycling rule
1:6 + c(10, 20)
# replicate a single element
rep("a", 5)
# replicate the whole vector several times
rep(c("a", "b"), 5)
rep(c("a", "b"), times=5)
# replicate the first element, then the second, etc.
rep(c("a", "b"), each=5)
# replicate to specified length
rep(c("a", "b"), length.out=5)
# define vector and matrix
vec_tor1 <- c(2, 4, 3)
mat_rix <- matrix(sample(1:12), ncol=3)
# multiply matrix by vector column-wise
vec_tor1 * mat_rix
mat_rix * vec_tor1
# multiply matrix by vector row-wise
t(vec_tor1 * t(mat_rix))
vec_tor1
vec_tor2 <- 6:4  # define vector
# multiply two vectors element-by-element
vec_tor1 * vec_tor2
# calculate inner product
vec_tor1 %*% vec_tor2
# calculate inner product and drop dimensions
drop(vec_tor1 %*% vec_tor2)
# multiply columns of matrix by vector
 mat_rix %*% vec_tor1  # single column matrix
 vec_tor1 %*%mat_rix
drop(mat_rix %*% vec_tor1)  # vector
rowSums(t(vec_tor1 * t(mat_rix)))
# using rowSums() and t() is 10 times slower than %*%
library(microbenchmark)
summary(microbenchmark(
  in_ner=drop(mat_rix %*% vec_tor1),
  row_sums=rowSums(t(vec_tor1 * t(mat_rix))),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary
library(microbenchmark)
# multiply matrix by vector fails because dimensions aren't conformable
vec_tor1 %*% mat_rix
# works after transpose
drop(vec_tor1 %*% t(mat_rix))
# calculate inner product
crossprod(vec_tor1, vec_tor2)
# create matrix and vector
mat_rix <- matrix(1:3000, ncol=3)
tmat_rix <- t(mat_rix)
vec_tor <- 1:3

my_v1 <- 1:4  #4x1
my_v2 <- 1:3  #3x1
my_m1 <- matrix(1:12,ncol=3)  #4x3
my_m2 <- matrix(1:9,ncol=3)  #3x3

my_m1*my_m2
my_m2*my_m1

my_m1%*%my_m2
my_m2%*%my_m1

my_v1*my_m1  #4x1 * 4*3 = 4x3
my_m1*my_v1  #4x3 * 4*1 = 4x3 //same as above

my_v1%*%my_m1  #4x1 %*% 4x3 = 1x3
my_m1%*%my_v1  #4x3 %*% 4x1  //error

my_v2*my_m1  #3x1 * 4x3 = 4x3
my_m1*my_v2  #4x3 * 3x1 = 4x3 //same as above

my_v2%*%my_m1  #3x1 %*% 4x3 error
my_m1%*%my_v2  #4x3 %*% 3x1 = 4x1

my_v2*my_m2  #3x1 * 3x3 = 3x3
my_m2*my_v2  #3x3 * 3x3 = 3x3 //same as above

my_v2%*%my_m2  #3x1 %*% 3*3 = 1x3
my_m2%*%my_v2  #3x3 %*% 3*1 = 3x1

# crossprod is slightly faster than "%*%" operator
summary(microbenchmark(
  cross_prod=crossprod(tmat_rix, vec_tor),
  inner_prod=mat_rix %*% vec_tor,
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary
structure(microbenchmark)
# define named vectors
vec_tor1 <- sample(1:4)
names(vec_tor1) <-
  paste0("row", 1:4, "=", vec_tor1)
vec_tor1
vec_tor2 <- sample(1:4)
names(vec_tor2) <-
  paste0("col", 1:3, "=", vec_tor2)
vec_tor2
# calculate outer product of two vectors
mat_rix <- outer(vec_tor1, vec_tor2)
mat_rix
mat_rix2 <- outer(vec_tor2, vec_tor1)
mat_rix2
# calculate vectorized function spanned over two vectors
mat_rix <- outer(vec_tor1, vec_tor2,
           FUN=function(x1, x2) x2*sin(x1))
mat_rix
mat_rix <- outer(vec_tor1, vec_tor2,
                 FUN=function(x1, x2) x2)
mat_rix
# define a function with two arguments
test_func <- function(first_arg, second_arg) {  # body
  first_arg + second_arg  # returns last evaluated statement
}  # end test_func

test_func(1, 2)  # apply the function
args(test_func)  # display argument

# define function that uses variable from enclosure environment
test_func <- function(first_arg, second_arg) {
  first_arg + second_arg + glob_var
}  # end test_func

test_func(3, 2)  # error - glob_var doesn't exist yet!
glob_var <- 10  # create glob_var
test_func(3, 2)  # now works
test_func <- function(first_arg, second_arg) {
# last statement of function is return value
  first_arg + 2*second_arg
}  # end test_func
test_func(first_arg=3, second_arg=2)  # bind by name
test_func(first=3, second=2)  # partial name binding
test_func(3, 2)  # bind by position
test_func(second_arg=2, 3)  # mixed binding
test_func(3, 2, 1)  # too many arguments
test_func(2)  # not enough arguments
# function "paste" has two arguments with default values
str(paste)
# default values of arguments can be specified in argument list
test_func <- function(first_arg, fac_tor=1) {
  fac_tor*first_arg
}  # end test_func
test_func(3)  # default value used for second argument
test_func(3, 10)  # default value over-ridden
# default values can be a vector of strings
test_func <- function(in_put=c("first_val", "second_val")) {
  in_put <- match.arg(in_put)  # match to arg list
  in_put
}  # end test_func
test_func("second_val")
test_func("se")  # partial name binding
test_func("some_val")  # invalid string
# define function that returns NULL for non-numeric argument
test_func <- function(in_put) {
  if (!is.numeric(in_put)) {
    warning(paste("argument", in_put, "isn't numeric"))
    return(NULL)
  }
  2*in_put
}  # end test_func

test_func(2)
test_func("hello")
library(quantmod)
mtcars
sample(32,10)
some_cars <- mtcars[sample(NROW(mtcars), 10), ]
# plot scatterplot horsepower vs miles per gallon
plot(some_cars[, "hp"], some_cars[, "mpg"],
     xlab="horsepower", ylab="miles per gallon",
     main="miles per gallon vs horsepower")
# add a solid red point (pch=16) for the last car
points(x=some_cars[NROW(some_cars), "hp"],
 y=some_cars[NROW(some_cars), "mpg"],
 col="red", pch=16)
# add labels with the car names
text(x=some_cars[, "hp"], y=some_cars[, "mpg"],labels=rownames(some_cars),
     pos=1, cex=0.8)
rownames(some_cars)
rownames(some_cars[, ])
# or add labels using wordcloud, to prevent overlaps
library(wordcloud)
textplot(x=some_cars[, "hp"], y=some_cars[, "mpg"],
   words=rownames(some_cars))
# plot the tree Height
plot(trees[, "Height"],
     type="l",
     lwd=2,
     col="blue",
     main="Tree heights and volumes",
     xlab="tree number", ylab="",
     ylim=c(min(trees[, c("Height", "Volume")]),
      max(trees[, c("Height", "Volume")])))

aa <- trees[, c("Height", "Volume")]
aa
# plot the tree Volume
lines(trees[, "Volume"], lwd=2, col="green")
# add legend
legend(x="left", legend=c("Height", "Volume"),
 inset=0.1, cex=1.0, bg="white",
 lwd=2, lty=c(1, 1), col=c("blue", "green"))
x_var <- seq(-2*pi, 2*pi, len=100)  # x values

# open Windows graphics device
x11(width=11, height=7, title="simple plot")


install.packages(XQuartz)
# plot a sine function using basic line plot
plot(x=x_var, y=sin(x_var), xlab="x-val",
     ylab="y-val", type="l", lwd=2, col="red")
# add a cosine function
lines(x=x_var, y=cos(x_var), lwd=2, col="blue")
# add title
title(main="sine and cosine functions", line=0.1)
# add legend
legend(x="topright", legend=c("sine", "cosine"),
 title="legend", inset=0.1, cex=1.0, bg="white",
 lwd=2, lty=c(1, 1), col=c("red", "blue"))
graphics.off()  # close all graphics devices
par(mar=c(7, 2, 1, 2), mgp=c(2, 1, 0), cex.lab=0.8, cex.axis=0.8, cex.main=0.8, cex.sub=0.5)
# plot a Normal probability distribution
curve(expr=dnorm, type="l", xlim=c(-3, 3),
xlab="", ylab="", lwd=2, col="blue")
# add shifted Normal probability distribution
curve(expr=dnorm(x, mean=1), add=TRUE,
type="l", lwd=2, col="red")

# add title
title(main="Normal probability distribution functions",
line=0.1)
# add legend
legend(x="topright", legend=c("Normal", "shifted"),
 title="legend", inset=0.05, cex=0.8, bg="white",
 lwd=2, lty=c(1, 1), col=c("blue", "red"))
graph_params <- par()  # get existing parameters
par("mar")  # get plot margins
par(mar=c(2, 1, 2, 1))  # set plot margins
par(oma=c(1, 1, 1, 1))  # set outer margins
par(mgp=c(2, 1, 0))  # set title and label margins
par(cex.lab=0.8,  # set font scales
    cex.axis=0.8, cex.main=0.8, cex.sub=0.5)
par(las=1)  # set axis labels to horizontal
par(ask=TRUE)  # pause, ask before plotting
par(mfrow=c(2, 2))  # plot on 2x2 grid by rows
for (i in 1:4) {  # plot 4 panels
  barplot(sample(1:6), main=paste("panel", i),
    col=rainbow(6), border=NA, axes=FALSE)
  box()
}  # end for
par(ask=FALSE)  # restore automatic plotting
par(new=TRUE)  # allow new plot on same chart
par(graph_params)  # restore original parameters
x_var <- seq(-5, 7, length=1000)
y_var <- dnorm(x_var, mean=1.0, sd=2.0)
plot(x_var, y_var, type="l", lty="solid",
     xlab="", ylab="")
title(main="Normal Density Function", line=0.5)
star_t <- 3; fin_ish <- 5  # set lower and upper bounds
# set polygon base
are_a <- ((x_var >= star_t) & (x_var <= fin_ish))
polygon(c(star_t, x_var[are_a], fin_ish),  # draw polygon
  c(-1, y_var[are_a], -1), col="red")
par(mar=c(7, 2, 1, 2), mgp=c(2, 1, 0), cex.lab=0.8, cex.axis=0.8, cex.main=0.8, cex.sub=0.5)
sig_mas <- c(0.5, 1, 1.5, 2)  # sigma values
# create plot colors
col_ors <- c("red", "black", "blue", "green")
# create legend labels
lab_els <- paste("sigma", sig_mas, sep="=")
for (in_dex in 1:4) {  # plot four curves
curve(expr=dnorm(x, sd=sig_mas[in_dex]),
type="l", xlim=c(-4, 4),
xlab="", ylab="", lwd=2,
col=col_ors[in_dex],
add=as.logical(in_dex-1)
)
}  # end for
# add title
title(main="Normal Distributions", line=0.5)
# add legend
legend("topright", inset=0.05, title="Sigmas",
 lab_els, cex=0.8, lwd=2, lty=c(1, 1, 1, 1),
 col=col_ors)
rm(list=ls())
par(mar=c(7, 2, 1, 2), mgp=c(2, 1, 0), cex.lab=0.8, cex.axis=0.8, cex.main=0.8, cex.sub=0.5)
x_var <- seq(-4, 4, length=100)
sig_mas <- c(0.5, 1, 1.5, 2)  # sigma values
# create plot colors
col_ors <- c("red", "black", "blue", "green")
# create legend labels
lab_els <- paste("sigma", sig_mas, sep="=")
# plot the first chart
plot(x_var, dnorm(x_var, sd=sig_mas[1]),
     type="n", xlab="", ylab="",
     main="Normal Distributions")
# add lines to plot
for (in_dex in 1:4) {
  lines(x_var, dnorm(x_var, sd=sig_mas[in_dex]),
  lwd=2, col=col_ors[in_dex])
}  # end for
# add legend
legend("topright", inset=0.05, title="Sigmas",
 lab_els, cex=0.8, lwd=2, lty=c(1, 1, 1, 1),
 col=col_ors)
x11(width=6, height=5)
par(mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
deg_free <- c(2, 5, 8, 11)  # df values
# create plot colors
col_ors <- c("red", "black", "blue", "green")
# create legend labels
lab_els <- paste("df", deg_free, sep="=")
for (in_dex in 1:4) {  # plot four curves
curve(expr=dchisq(x, df=deg_free[in_dex]),
      type="l", xlim=c(0, 20), ylim=c(0, 0.3),
      xlab="", ylab="", lwd=2,
      col=col_ors[in_dex],
      add=as.logical(in_dex-1))
}  # end for
# add title
title(main="Chi-squared Distributions", line=0.5)
# add legend
legend("topright", inset=0.05,
       title="Degrees of freedom", lab_els,
       cex=0.8, lwd=6, lty=c(1, 1, 1, 1),
       col=col_ors)
x11(width=6, height=5)
par(mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
deg_free <- c(2, 5, 8, 11)  # df values
# create plot colors
col_ors <- c("red", "black", "blue", "green")
# create legend labels
lab_els <- paste("df", deg_free, sep="=")
# plot an empty chart
x_var <- seq(0, 20, length=100)
plot(x_var, dchisq(x_var, df=deg_free[1]),
     type="n", xlab="", ylab="", ylim=c(0, 0.3))
# add lines to plot
for (in_dex in 1:4) {
  lines(x_var, dchisq(x_var, df=deg_free[in_dex]),
lwd=2, col=col_ors[in_dex])
}  # end for
# add title
title(main="Chi-squared Distributions", line=0.5)
# add legend
legend("topright", inset=0.05,
       title="Degrees of freedom", lab_els,
       cex=0.8, lwd=6, lty=c(1, 1, 1, 1),
       col=col_ors)
x11(width=6, height=5)
par(mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
deg_free <- c(3, 6, 9)  # df values
col_ors <- c("black", "red", "blue", "green")
lab_els <- c("normal", paste("df", deg_free, sep="="))
# plot a Normal probability distribution
curve(expr=dnorm, type="l", xlim=c(-4, 4),
      xlab="", ylab="", lwd=2)
for (in_dex in 1:3) {  # plot three t-distributions
curve(expr=dt(x, df=deg_free[in_dex]),
      type="l", xlab="", ylab="", lwd=2,
      col=col_ors[in_dex+1], add=TRUE)
}  # end for
# add title
title(main="t-distributions", line=0.5)
# add legend
legend("topright", inset=0.05,
       title="Degrees\n of freedom", lab_els,
       cex=0.8, lwd=6, lty=c(1, 1, 1, 1),
       col=col_ors)
x11(width=6, height=5)
par(mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
x_var <- seq(-4, 4, length=100)
deg_free <- c(3, 6, 9)  # df values
col_ors <- c("black", "red", "blue", "green")
lab_els <- c("normal", paste("df", deg_free, sep="="))
# plot chart of normal distribution
plot(x_var, dnorm(x_var), type="l",
     lwd=2, xlab="", ylab="")
for (in_dex in 1:3) {  # add lines for t-distributions
  lines(x_var, dt(x_var, df=deg_free[in_dex]),
lwd=2, col=col_ors[in_dex+1])
}  # end for
# add title
title(main="t-distributions", line=0.5)
# add legend
legend("topright", inset=0.05,
       title="Degrees\n of freedom", lab_els,
       cex=0.8, lwd=6, lty=c(1, 1, 1, 1),
       col=col_ors)
x11(width=6, height=5)
par(mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
deg_free <- c(3, 5, 9)  # df values
col_ors <- c("black", "red", "blue", "green")
lab_els <- paste0("df1=", deg_free, ", df2=3")
for (in_dex in 1:NROW(deg_free)) {  # plot four curves
curve(expr=df(x, df1=deg_free[in_dex], df2=3),
      type="l", xlim=c(0, 4),
      xlab="", ylab="", lwd=2,
      col=col_ors[in_dex],
      add=as.logical(in_dex-1))
}  # end for
# add title
title(main="F-Distributions", line=0.5)
# add legend
legend("topright", inset=0.05, title="degrees of freedom",
       lab_els, cex=0.8, lwd=2, lty=1,
       col=col_ors)
rm(list=ls())
par(mar=c(7, 2, 1, 2), mgp=c(2, 1, 0), cex.lab=0.8, cex.axis=0.8, cex.main=0.8, cex.sub=0.5)
poisson_events <- 0:11  # Poisson events
poisson_freq <- dpois(poisson_events, lambda=4)
names(poisson_freq) <- as.character(poisson_events)
# Poisson function
poisson_func <- function(x, lambda)
              {exp(-lambda)*lambda^x/factorial(x)}
curve(expr=poisson_func(x, lambda=4), xlim=c(0, 11), main="Poisson distribution",
xlab="No. of events", ylab="Frequency of events", lwd=2, col="red")
legend(x="topright", legend="Poisson density", title="",
 inset=0.05, cex=0.8, bg="white", lwd=4, lty=1, col="red")
set.seed(1121)  # reset random number generator
runif(3)  # three numbers from uniform distribution
runif(3)  # produce another three numbers
set.seed(1121)  # reset random number generator
runif(3)  # produce another three numbers

# produce random number from standard normal distribution
rnorm(1)
# produce five random numbers from standard normal distribution
rnorm(5)
# produce five random numbers from the normal distribution
rnorm(n=5, mean=1, sd=2)  # match arguments by name
# calculate cumulative standard normal distribution
c(pnorm(-2), pnorm(2))
# calculate inverse cumulative standard normal distribution
c(qnorm(0.75), qnorm(0.25))
# define logistic map function
log_map <- function(x, r=4) r*x*(1-x)
log_map(0.25, 4)
# plot logistic map
x11(width=6, height=5)
curve(expr=log_map, type="l", xlim=c(0, 1),
xlab="x[n-1]", ylab="x[n]", lwd=2, col="blue",
main="logistic map")
lines(x=c(0, 0.25), y=c(0.75, 0.75), lwd=2, col="blue")
lines(x=c(0.25, 0.25), y=c(0, 0.75), lwd=2, col="blue")
# calculate uniformly distributed pseudo-random
# sequence using logistic map function
uni_form <- function(see_d, len_gth=10) {
  # pre-allocate vector instead of "growing" it
  out_put <- numeric(len_gth)
  # initialize
  out_put[1] <- see_d
  # perform loop
  for (i in 2:len_gth) {
    out_put[i] <- 4*out_put[i-1]*(1-out_put[i-1])
  }  # end for
  acos(1-2*out_put)/pi
}  # end uni_form
uni_form(see_d=0.1, len_gth=15)
plot(
  density(uni_form(see_d=runif(1), len_gth=1e5)),
  xlab="", ylab="", lwd=2, col="blue",
  main="uniform pseudo-random number density")
set.seed(1121)  # reset random number generator
# flip unbiased coin once, 20 times
rbinom(n=20, size=1, 0.5)
# number of heads after flipping twice, 20 times
rbinom(n=20, size=2, 0.5)
# number of heads after flipping thrice, 20 times
rbinom(n=2, size=3, 0.5)
# number of heads after flipping biased coin thrice, 20 times
rbinom(n=20, size=3, 0.8)
# number of heads after flipping biased coin thrice, 20 times
rbinom(n=20, size=3, 0.2)
# flip unbiased coin once, 20 times
sample(x=0:1, size=20, replace=TRUE)  # fast
as.numeric(runif(20) < 0.5)  # slower
# permutation of five numbers
sample(x=5)
# permutation of four strings
sample(x=c("apple", "grape", "orange", "peach"))
# sample of size three
sample(x=5, size=3)
# sample with replacement
sample(x=5, replace=TRUE)
sample(  # sample of strings
  x=c("apple", "grape", "orange", "peach"),
  size=12,
  replace=TRUE)
# binomial sample: flip coin once, 20 times
sample(x=0:1, size=20, replace=TRUE)
# flip unbiased coin once, 20 times
as.numeric(runif(20) > 0.5)  # slower
rm(list=ls())
set.seed(1121)  # reset random number generator
# sample from Standard Normal Distribution
sam_ple <- rnorm(1000)

mean(sam_ple)  # sample mean

median(sam_ple)  # sample median

sd(sam_ple)  # sample standard deviation
rm(list=ls())
# DAX returns
re_turns <- diff(log(EuStockMarkets[, 1]))
# number of observations
len_gth <- NROW(re_turns)
# mean of DAX returns
mean_rets <- mean(re_turns)
# standard deviation of DAX returns
sd_rets <- sd(re_turns)
# skew of DAX returns
len_gth/((len_gth-1)*(len_gth-2))*
  sum(((re_turns - mean_rets)/sd_rets)^3)
# kurtosis of DAX returns
len_gth*(len_gth+1)/((len_gth-1)^3)*
  sum(((re_turns - mean_rets)/sd_rets)^4)
# random normal returns
re_turns <- rnorm(len_gth, sd=2)
# mean and standard deviation of random normal returns
mean_rets <- mean(re_turns)
sd_rets <- sd(re_turns)
# skew of random normal returns
len_gth/((len_gth-1)*(len_gth-2))*
  sum(((re_turns - mean_rets)/sd_rets)^3)
# kurtosis of random normal returns
len_gth*(len_gth+1)/((len_gth-1)^3)*
  sum(((re_turns - mean_rets)/sd_rets)^4)
set.seed(1121)  # reset random number generator
# sample from Standard Normal Distribution
len_gth <- 1000
sam_ple <- rnorm(len_gth)
# sample mean
mean(sam_ple)
# sample standard deviation
sd(sam_ple)
#Perform two-tailed test that sample is
#from Standard Normal Distribution (mean=0, SD=1)
# generate vector of samples and store in data frame
test_frame <- data.frame(samples=rnorm(1e4))
# get p-values for all the samples
test_frame$p_values <- sapply(test_frame$samples, 
        function(x) 2*pnorm(-abs(x)))
# significance level, two-tailed test, critical value=2*SD
signif_level <- 2*(1-pnorm(2))
# compare p_values to significance level
test_frame$result <-
  test_frame$p_values > signif_level
# number of null rejections
sum(!test_frame$result) / NROW(test_frame)
# show null rejections
head(test_frame[!test_frame$result, ])
x11(width=6, height=5)
par(mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
# plot the Normal probability distribution

curve(expr=dnorm(x, sd=1), type="l", xlim=c(-4, 4),
xlab="", ylab="", lwd=3, col="blue")
title(main="Two-tailed Test", line=0.5)
# plot tails of the distribution using polygons
star_t <- 2; e_nd <- 4
# plot right tail using polygon
x_var <- seq(star_t, e_nd, length=100)
y_var <- dnorm(x_var, sd=1)
y_var[1] <- (-1)
y_var[NROW(y_var)] <- (-1)
polygon(x=x_var, y=y_var, col="green")
# plot left tail using polygon
y_var <- dnorm(-x_var, sd=1)
y_var[1] <- (-1)
y_var[NROW(y_var)] <- (-1)
polygon(x=(-x_var), y=y_var, col="red")
rm(list=ls())
par(oma=c(1, 1, 1, 1), mgp=c(2, 0.5, 0), mar=c(5, 1, 1, 1), cex.lab=0.8, cex.axis=0.8, cex.main=0.8, cex.sub=0.5)
library(ggplot2)  # load ggplot2

qplot(  # simple ggplot2
    main="Standard Normal Distribution",
    c(-4, 4),
    stat="function",
    fun=dnorm,
    geom="line",
    xlab=NULL, ylab=NULL
    ) +  # end qplot

theme(  # modify plot theme
    plot.title=element_text(vjust=-1.0),
    plot.background=element_blank()
    ) +  # end theme

geom_vline(  # add vertical line
  aes(xintercept=c(-2.0, 2.0)),
  colour="red",
  linetype="dashed"
  )  # end geom_vline
rm(list=ls())
par(oma=c(1, 1, 1, 1), mgp=c(2, 0.5, 0), mar=c(5, 1, 1, 1), cex.lab=0.8, cex.axis=0.8, cex.main=0.8, cex.sub=0.5)
#create ggplot2 with shaded area
x_var <- -400:400/100
norm_frame <- data.frame(x_var=x_var,
                 d.norm=dnorm(x_var))
norm_frame$shade <- ifelse(
            abs(norm_frame$x_var) >= 2,
            norm_frame$d.norm, NA)
ggplot(  # main function
  data=norm_frame,
  mapping=aes(x=x_var, y=d.norm)
  ) +  # end ggplot
# plot line
  geom_line() +
# plot shaded area
  geom_ribbon(aes(ymin=0, ymax=shade), fill="red") +
# no axis labels
  xlab("") + ylab("") +
# add title
  ggtitle("Standard Normal Distribution") +
# modify plot theme
  theme(
  plot.title=element_text(vjust=-1.0),
  plot.background=element_blank()
  )  # end theme
# t-test for single sample
t.test(rnorm(100))
# t-test for two samples
t.test(rnorm(100),
       rnorm(100, mean=1))
# Wilcoxon test for normal distribution
wilcox.test(rnorm(100))
# Wilcoxon test for two normal distributions
wilcox.test(rnorm(100), rnorm(100, mean=0.1))
# Wilcoxon test for two normal distributions
wilcox.test(rnorm(100), rnorm(100, mean=1.0))
# Wilcoxon test for a uniform versus normal distribution
wilcox.test(runif(100), rnorm(100))
# KS test for normal distribution
ks.test(rnorm(100), pnorm)
# KS test for uniform distribution
ks.test(runif(100), pnorm)
# KS test for two similar normal distributions
ks.test(rnorm(100), rnorm(100, mean=0.1))
# KS test for two different normal distributions
ks.test(rnorm(100), rnorm(100, mean=1.0))
# calculate DAX percentage returns
dax_rets <- diff(log(EuStockMarkets[, 1]))

# Shapiro-Wilk test for normal distribution
shapiro.test(rnorm(NROW(dax_rets)))

# Shapiro-Wilk test for DAX returns
shapiro.test(dax_rets)

# Shapiro-Wilk test for uniform distribution
shapiro.test(runif(NROW(dax_rets)))
dax_rets <- diff(log(EuStockMarkets[, 1]))
library(tseries)  # load package tseries

# Jarque-Bera test for normal distribution
jarque.bera.test(rnorm(NROW(dax_rets)))

# Jarque-Bera test for DAX returns
jarque.bera.test(dax_rets)

# Jarque-Bera test for uniform distribution
jarque.bera.test(runif(NROW(dax_rets)))
# verify that rtools are working properly:
devtools::find_rtools()
devtools::has_devel()

# load package Rcpp
library(Rcpp)
# get documentation for package Rcpp
# get short description
packageDescription("Rcpp")
# load help page
help(package="Rcpp")
# list all datasets in "Rcpp"
data(package="Rcpp")
# list all objects in "Rcpp"
ls("package:Rcpp")
# remove Rcpp from search path
detach("package:Rcpp")

# define Rcpp function
Rcpp::cppFunction("
  int times_two(int x)
    { return 2 * x;}
  ")  # end cppFunction
# run Rcpp function
times_two(3)
# source Rcpp functions from file
Rcpp::sourceCpp(file="C:/Develop/R/lecture_slides/scripts/rcpp_mult.cpp")
# multiply two numbers
rcpp_mult(2, 3)
rcpp_mult(1:3, 6:4)
# multiply two vectors
rcpp_mult_vec(2, 3)
rcpp_mult_vec(1:3, 6:4)
# define Rcpp function with loop
Rcpp::cppFunction("
double inner_mult(NumericVector x, NumericVector y) {
int x_size = x.size();
int y_size = y.size();
if (x_size != y_size) {
    return 0;
  } else {
    double total = 0;
    for(int i = 0; i < x_size; ++i) {
total += x[i] * y[i];
  }
  return total;
  }
}")  # end cppFunction
# run Rcpp function
inner_mult(1:3, 6:4)
inner_mult(1:3, 6:3)
# define Rcpp Sugar function with loop
Rcpp::cppFunction("
double inner_mult_sugar(NumericVector x, NumericVector y) {
  return sum(x * y);
}")  # end cppFunction
# run Rcpp Sugar function
inner_mult_sugar(1:3, 6:4)
inner_mult_sugar(1:3, 6:3)
# define R function with loop
inner_mult_r <- function(x, y) {
    to_tal <- 0
    for(i in 1:NROW(x)) {
to_tal <- to_tal + x[i] * y[i]
    }
    to_tal
}  # end inner_mult_r
# run R function
inner_mult_r(1:3, 6:4)
inner_mult_r(1:3, 6:3)
# compare speed of Rcpp and R
library(microbenchmark)
summary(microbenchmark(
  pure_r=inner_mult_r(1:10000, 1:10000),
  inner_r=1:10000 %*% 1:10000,
  r_cpp=inner_mult(1:10000, 1:10000),
  r_cpp_sugar=inner_mult_sugar(1:10000, 1:10000),
  times=10))[, c(1, 4, 5)]
# calculate uniformly distributed pseudo-random sequence
uni_form <- function(see_d, len_gth=10) {
  out_put <- numeric(len_gth)
  out_put[1] <- see_d
  for (i in 2:len_gth) {
    out_put[i] <- 4*out_put[i-1]*(1-out_put[i-1])
  }  # end for
  acos(1-2*out_put)/pi
}  # end uni_form

# source Rcpp functions from file
Rcpp::sourceCpp(file="C:/Develop/R/lecture_slides/scripts/uni_form.cpp")
# microbenchmark Rcpp code
library(microbenchmark)
summary(microbenchmark(
  pure_r=runif(1e5),
  r_loop=uni_form(0.3, 1e5),
  r_cpp=uniform_rcpp(0.3, 1e5),
  times=10))[, c(1, 4, 5)]

