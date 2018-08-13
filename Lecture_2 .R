# create factor vector
fac_tor <- factor(c('b', 'c', 'd', 'a', 'c', 'b', 'z','zz'))
fac_tor

fac_tor <- factor(c('b', 'c', 'd', 'a', 'c', 'b'))
fac_tor

fac_tor[5]
attributes(fac_tor)  # get factor attributes
levels(fac_tor) 
is.vector(levels(fac_tor) )# get allowed values
as.numeric(fac_tor)  # get encoding vector number underlyng these letter. a=1, b=2, z=5.
# like 25 
is.vector(fac_tor)
as.factor(1:5)  # coerce vector to factor
# coerce factor to character vector
as.vector(as.factor(1:5))
fac_tor
levels(fac_tor)  # get allowed values
unique(fac_tor)  # get unique elements
# get contingency (frequency) table
table(fac_tor)
# get contingency table using sapply
ee <- sapply(levels(fac_tor),
 function(le_vel) {
   sum(fac_tor==le_vel)
 })  # end sapply
ee
sum(fac_tor=="a")
le_vel <- "b"
ss <- fac_tor == le_vel
ss
sum(fac_tor == le_vel)


library(microbenchmark)
str(findInterval)
# get index of the element of "vec" that matches 5
findInterval(x=5, vec=c(3, 5, 7)) #????????????????????????5??????[5,7)???????????????
match(3, c(3, 5, 7)) #2.
match(4, c(3, 5, 7)) #NA
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
# named numeric vector of breakpoints
brea_ks <- c(freezing=0, very_cold=30,
       cold=50, pleasant=60,
       warm=80, hot=90)
brea_ks
tempe_ratures <- runif(10, min=10, max=100)
feels_like <- names(
  brea_ks[findInterval(x=tempe_ratures,
                 vec=brea_ks)])
names(tempe_ratures) <- feels_like
tempe_ratures
library(microbenchmark)
foo <- sample(0:6) + 0.1
foo

cut(x=foo, breaks=c(2, 4, 6, 8))
rbind(foo, cut(x=foo, breaks=c(2, 4, 6, 8)))
# cut() replicates findInterval()
cut(x=1:8, breaks=c(3, 5, 7), labels=1:2,
    right=FALSE)
findInterval(x=1:8, vec=c(3, 5, 7))
# findInterval() is a compiled function, so it's faster than cut()
vec_tor <- rnorm(1000)
summary(microbenchmark(
  find_interval=
    findInterval(x=vec_tor, vec=c(3, 5, 7)),
  cuut=
    cut(x=vec_tor, breaks=c(3, 5, 7)),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary
setwd("C:/Develop/R/lecture_slides/data")
# define a function that returns invisibly
return_invisible <- function(in_put) {
  invisible(in_put)
}  # end return_invisible

return_invisible(2) #????????????return??????????????????

glob_var <- return_invisible(2)
glob_var

return_invisible <- function(in_put) {
  in_put <- 33*in_put
  invisible(in_put)#??????return????????????invisible???????????????????????????????????????????????????
}  # end return_invisible

glob_var <- return_invisible(2)
glob_var

rm(list=ls())  # remove all objects
# load objects from file
loaded <- load(file="/Users/bessie.chan/Desktop/nyu/courses/6871 R/my_data.RData")
loaded  # vector of loaded objects #?????????????????????Rdata??????????????????????????????????????????load???????????????????????????
ls()  # list objects
calc_skew  # show the function code

getAnywhere(calc_skew)  # display function
str(plot)  # dots for additional plot parameters
bind_dots <- function(in_put, ...) {
  paste0("in_put=", in_put,
 ", dots=", paste(..., sep=", "))
}  # end bind_dots
bind_dots(1, 2, 3)  # "in_put" bound by position
bind_dots(2, in_put=1, 3)  # "in_put" bound by name
bind_dots(2, in_put="c", 3)
bind_dots(2, in_put="1,3", 3)
bind_dots(2, in_put <- "1,3", 3) #[1] "in_put=2, dots=1,3, 3" ?????? <- R?????????in_put ???1???3???
bind_dots(1, 2, 3, foo=10)  # named argument bound to dots???...????????????????????????????????????foo=10
list(firs_t = 'a', se_cond = 'b')
list(firs_t <-  'a', se_cond <- 'b')
c(firs_t = 'a', se_cond = 'b')
c(firs_t <-  'a', se_cond <- 'b')
bind_dots <- function(arg1, arg2, ...) {
  arg1 + 2*arg2 + sum(...)
}  # end bind_dots
bind_dots(3, 2)  # bind arguments by position
bind_dots(3, 2, 5, 8)  # extra arguments bound to dots
str(sum)  # dots before other arguments
sum(1, 2, 3)  # dots bind before other arguments
sum(1, 2, NA, 3, na.rm=TRUE)
bind_dots <- function(..., in_put) {
  paste0("in_put=", in_put,
 ", dots=", paste(..., sep=", "))
}  # end bind_dots
# arguments after dots must be bound by full name
bind_dots(1, 2, 3, in_put=10)
bind_dots(1, 2, 3, in_put=10, foo=4)  # dots bound
bind_dots(1, 2, 3)  # "in_put" not bound
bind_dots <- function(..., in_put=10) {
  paste0("in_put=", in_put,
 ", dots=", paste(..., sep=", "))
}  # end bind_dots
bind_dots(1, 2, 3)  # "in_put" not bound, but has default
bind_dots(1, 2, 3,in_put = 11) 
# wrapper for mean() with default na.rm=TRUE
my_mean <- function(x, na.rm=TRUE, ...) {
  mean(x=x, na.rm=na.rm, ...)
}  # end my_mean
foo <- sample(c(1:10, NA, rep(0.1, t=5)))
mean(c(foo, NA))
mean(c(foo, NA), na.rm=TRUE)
my_mean(c(foo, NA))
my_mean(c(foo, NA), trim=0.4)  # pass extra argument
# wrapper for saving data into default directory
save_data <- function(...,
              file=stop("error: no file name"),
              my_dir="/Users/bessie.chan/Desktop/nyu/courses/6871 R") {
# create file path
  file <- file.path(my_dir, file)
  save(..., file=file)
}  # end save_data
foo <- 1:10
save_data(foo, file="scratch.RData")
save_data(foo, file="scratch.RData", my_dir="/Users/bessie.chan/Desktop/nyu/courses/6871 R")
# wrapper for testing negative arguments
stop_if_neg <- function(in_put) {
  if (!is.numeric(in_put) || in_put<0)
    stop("argument not numeric or negative")
}  # end stop_if_neg
# wrapper for sqrt()
my_sqrt <- function(in_put) {
  stop_if_neg(in_put)
  sqrt(in_put)
}  # end my_sqrt
sqrt(2)
sqrt(-2)
my_sqrt(2)
my_sqrt(-2)
my_sqrt(NA)
setwd("/Users/bessie.chan/Desktop/nyu/courses/6871 R/my_data.RData")
rm(list=ls())  # remove all objects
ls()  # list objects
# load objects from file (side effect)
load(file="/Users/bessie.chan/Desktop/nyu/courses/6871 R/my_data.RData")
ls()  # list objects
glob_var <- 1  # define a global variable
# explore function scope and side effects
side_effect <- function() {
  cat("global glob_var:\t", glob_var, "\n")
# define local "glob_var" variable
  glob_var <- 10
# re-define the global "glob_var"
  glob_var <<- 2  #?????????global scope????????????
  cat("local glob_var:\t", glob_var, "\n")
}  # end side_effect
glob_var <- 1000
side_effect()#?????????????????????????????????????????????global scope???????????????
glob_var

side_effect2 <- function() {
  glob_var <- 10
  cat("global glob_var:\t", glob_var, "\n")
  # define local "glob_var" variable
  
  # re-define the global "glob_var"
  glob_var <<- 2  #?????????global scope????????????
  cat("local glob_var:\t", glob_var, "\n")
}
glob_var <- 1000
side_effect2()
glob_var
# global variable was modified as side effect
glob_var
# create functional that accepts a function as input argument
func_tional <- function(func_name) {
# calculates statistic on random numbers
  set.seed(1)
  func_name(runif(1e4))  # apply the function name
}  # end func_tional
func_tional(mean)
func_tional(sd)
# func_tional accepts function name and additional argument
func_tional <- function(func_name, in_put) {
# produce function name from argument
  func_name <- match.fun(func_name)
# execute function call
  func_name(in_put)
}  # end func_tional
func_tional(sqrt, 4)
# string also works because match.fun() converts it to a function
func_tional("sqrt", 4)
str(sum)  # sum() accepts multiple arguments
# func_tional can't accept indefinite number of arguments
func_tional(sum, 1, 2, 3)
# func_tional accepts function name and dots '...' argument
func_tional <- function(func_name, ...) {
  func_name <- match.fun(func_name)
  func_name(...)  # execute function call
}  # end func_tional

func_tional(sum, 1, 2, 3)
func_tional(sum, 1, 2, NA, 4, 5)
func_tional(sum, 1, 2, NA, 4, 5, na.rm=TRUE) #na.rm=TRUE?????????sum??????????????????????????????func_tional??????bind????????????
# function with three arguments and dots '...' arguments
my_func <- function(in_put, param1, param2, ...) {
  c(input=in_put, param1=param1, param2=param2,
dots=c(...))
}  # end my_func
my_func(1, 2, 3, param2=4, param1=5)
func_tional(my_func, 1, 2, 3, param2=4, param1=5)
func_tional(my_func, 1, 2, 3, 4, 5)
# simple anonymous function
(function(x) (x + 3)) (10)
# anonymous function passed to func_tional
func_tional(func_name=(function(x) (x + 3)), 5)
# anonymous function is default value
func_tional <-
  function(..., func_name=function(x, y, z) {x+y+z}) {
    func_name <- match.fun(func_name)
    func_name(...)  # execute function call
}  # end func_tional
func_tional(2, 3, 4)  # use default func_name
func_tional(2, 3, 4, 5)
# func_name bound by name
func_tional(func_name=sum, 2, 3, 4, 5)
# pass anonymous function to func_name
func_tional(func_name=function(x, y, z) {x*y*z},
    2, 3, 4)
str(sum)  # sum() accepts multiple arguments
# sum() can't accept list of arguments
sum(list(1, 2, 3))
str(do.call)  # "what" argument is a function
# do.call passes list elements into "sum" individually
do.call(sum, list(1, 2, 3))
do.call(sum, list(1, 2, NA, 3))
do.call(sum, list(1, 2, NA, 3, na.rm=TRUE))
# func_tional() accepts list with function name and arguments
func_tional <- function(list_arg) {
# produce function name from argument
  func_name <- match.fun(list_arg[[1]])
# execute function call uing do.call()
  do.call(func_name, list_arg[-1])
}  # end func_tional
arg_list <- list("sum", 1, 2, 3)
qq <- func_tional(arg_list)
rm(list=ls())
str(apply)  # get list of arguments
# create a matrix
mat_rix <- matrix(6:1, nrow=2, ncol=3)
mat_rix
# sum the rows and columns
row_sums <- apply(mat_rix, 1, sum)
row_sums
col_sums <- apply(mat_rix, 2, sum)
col_sums

NCOL(mat_rix)
foo <- numeric(NCOL(mat_rix))
foo
for(i in 1:NCOL(mat_rix)) {foo[i] <- sum(mat_rix[,i])}
foo
col_sums <- apply(mat_rix, 2, sum)
col_sums

c(sum(row_sums), row_sums)
rbind(col_sums, mat_rix)

mat_rix <- cbind(c(sum(row_sums), row_sums),
          rbind(col_sums, mat_rix))
mat_rix
dimnames(mat_rix) <- list(c("col_sums", "row1", "row2"),
                 c("row_sums", "col1", "col2", "col3"))
mat_rix
str(apply)  # get list of arguments
mat_rix <- matrix(sample(12), nrow=3, ncol=4)  # create a matrix
mat_rix
apply(mat_rix, 2, sort)  # sort matrix columns
apply(mat_rix, 2, sort, decreasing=TRUE)  # sort decreasing order #decreasing=TRUE???sort???????????????????????????
mat_rix[2, 2] <- NA  # introduce NA value
mat_rix
# calculate median of columns
apply(mat_rix, 2, median)
# calculate median of columns with na.rm=TRUE
apply(mat_rix, 2, median, na.rm=TRUE)
rm(list=ls())
# DAX percent returns
dax_rets <- 100*diff(log(EuStockMarkets[, 1]))
library(moments)  # load package moments
str(moment)  # get list of arguments
# apply moment function
moment(x=dax_rets, order=3)
# 4x1 matrix of moment orders
moment_orders <- as.matrix(1:4)
moment_orders
# anonymous function allows looping over function parameters
a <- apply(X=moment_orders, MARGIN=1,
      FUN=function(moment_order) {
  moment(x=dax_rets, order=moment_order)
}  # end anonymous function
      )  # end apply
a
a <- sapply(X=moment_orders,
           FUN=function(moment_order) {
             moment(x=dax_rets, order=moment_order)
           }  # end anonymous function
)  # end apply
a
# another way of passing parameters into moment() function
apply(X=moment_orders, MARGIN=1, FUN=moment,
      x=dax_rets)
str(moment)#function (x, order = 1, central = FALSE, absolute = FALSE, na.rm = FALSE)  
#????????????by name???x=dax_rets???????????????????????????????????????order??????bind by position?????????X=moment_orders????????????

# function with three arguments
my_func <- function(arg1, arg2, arg3) {
  c(arg1=arg1, arg2=arg2, arg3=arg3)
}  # end my_func
my_func(1, 2, 3)
da_ta <- as.matrix(1:4)
# pass da_ta to arg1
apply(X=da_ta, MAR=1, FUN=my_func, arg2=2, arg3=3)
?apply()
# pass da_ta to arg2
apply(X=da_ta, MAR=1, FUN=my_func, arg1=1, arg3=3)
# pass da_ta to arg3
apply(X=da_ta, MAR=1, FUN=my_func, arg1=1, arg2=2)
# vector of means of numeric columns
dim(iris)
head(iris)
class(iris)
sapply(iris, class)
apply(iris,2, class) #??????????????????
nrow(iris)
iris[, -5] #??????????????????
iris[, -4]
sapply(iris[, -5], mean)
# list of means of numeric columns
lapply(iris[, -5], mean)
# lapply using anonymous function
do.call(cbind,lapply(iris[, -5], mean)) #cbind. ???????????????list bind?????????????????????????????????vector
class(lapply(iris[, -5], mean))
unlist(lapply(iris[, -5], mean))
apply(iris[, -5], 2, mean)

unlist(lapply(iris,
      function(col_umn) {
        if (is.numeric(col_umn)) mean(col_umn)
      }  # end anonymous function
      )  # end lapply
       )  # end unlist
unlist(sapply(iris, function(col_umn) {
  if (is.numeric(col_umn)) mean(col_umn)}))
sapply(6:10, sqrt)  # sapply on vector
sapply(list(6, 7, 8, 9, 10), sqrt)  # sapply on list
lapply(6:10, sqrt)
lapply(list(6, 7, 8, 9, 10), sqrt)
# calculate means of iris data frame columns
sapply(iris, mean)  # returns NA for Species

# create a matrix
 
mat_rix <- matrix(sample(100), ncol=4)
# calculate column means using apply
apply(mat_rix, 2, mean)

# calculate column means using sapply, with anonymous function
sapply(1:NCOL(mat_rix),
       function(col_index) {  # anonymous function
 mean(mat_rix[, col_index])
  }  # end anonymous function
)  # end sapply
# vectors form columns of matrix returned by sapply
sapply(2:4, function(num) c(el1=num, el2=2*num))
# vectors of different lengths returned as list
sapply(2:4, function(num) 1:num)
# vapply is similar to sapply
vapply(2:4, function(num) c(el1=num, el2=2*num),
       FUN.VALUE=c(row1=0, row2=0))
# vapply produces an error if it can't simplify
vapply(2:4, function(num) 1:num,
       FUN.VALUE=c(row1=0, row2=0))
vapply(runif(3), function(num) c(el1=num, el2=2*num),
       FUN.VALUE=c(row1=0, row2=0))
vapply(c("a","d"), function(num) c(el1=num, el2=2*num),
       FUN.VALUE=c(row1=0, row2=0))
vec_tor1 <- rnorm(1000000)
vec_tor2 <- rnorm(1000000)
big_vector <- numeric(1000000)
system.time(  # sum vectors using "for" loop
  for (i in 1:NROW(vec_tor1)) {
    big_vector[i] <- vec_tor1[i] + vec_tor2[i]
  }  # end for
)  # end system.time
# sum vectors using vectorized "+"
system.time(big_vector <- vec_tor1 + vec_tor2) #???????????????
# allocate memory for cumulative sum
cum_sum <- numeric(NROW(big_vector))
# cumulative sum using "for" loop
cum_sum[1] <- big_vector[1]
system.time(
  for (i in 2:NROW(big_vector)) {
    cum_sum[i] <- cum_sum[i-1] + big_vector[i]
  }  # end for
)  # end system.time
# cumulative sum using "cumsum"
system.time(cum_sum <- cumsum(big_vector)) #???????????????????????????0??????
# calculate row sums two different ways
summary(microbenchmark(
  row_sums=rowSums(big_matrix),
  ap_ply=apply(big_matrix, 1, sum),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary
library(microbenchmark)
str(pmax)
# calculate row maximums two different ways
summary(microbenchmark(
  p_max=
    do.call(pmax.int,
lapply(seq_along(big_matrix[1, ]),
  function(in_dex) big_matrix[, in_dex])),
  l_apply=unlist(
    lapply(seq_along(big_matrix[, 1]),
  function(in_dex) max(big_matrix[in_dex, ]))),
  times=10))[, c(1, 4, 5)]
# create two numeric vectors
vec_tor1 <- sin(0.25*pi*1:10)
vec_tor2 <- cos(0.25*pi*1:10)
# create third vector using 'ifelse'
vec_tor3 <- ifelse(vec_tor1 > vec_tor2,
          vec_tor1, vec_tor2)
# cbind all three together
vec_tor4 <- cbind(vec_tor1, vec_tor2, vec_tor3)

# set plotting parameters
par(mar=c(7, 2, 1, 2), mgp=c(2, 1, 0),
    cex.lab=0.8, cex.axis=0.8, cex.main=0.8,
    cex.sub=0.5)
# plot matrix
matplot(vec_tor4, type="l", lty="solid",
col=c("green", "blue", "red"),
lwd=c(2, 2, 2), xlab="", ylab="")
# add legend
legend(x="bottomright", legend=colnames(vec_tor4),
       title="", inset=0.05, cex=0.8, lwd=2,
       lty=c(1, 1, 1), col=c("green", "blue", "red"))
library(parallel)  # load package parallel
# get short description
packageDescription("parallel")
# load help page
help(package="parallel")
# list all objects in "parallel"
ls("package:parallel")
library(parallel)  # load package parallel
# calculate number of available cores
num_cores <- detectCores() - 1
# define function that pauses execution
paws <- function(x, sleep_time) {
  Sys.sleep(sleep_time)
  x
}  # end paws
# perform parallel loop under Mac-OSX or Linux
paw_s <- mclapply(1:10, paws, mc.cores=num_cores,
          sleep_time=10)
# initialize compute cluster under Windows
clus_ter <- makeCluster(num_cores)
# perform parallel loop under Windows
paw_s <- parLapply(clus_ter, 1:10, paws,
           sleep_time=0.01)
library(microbenchmark)  # load package microbenchmark
# compare speed of lapply versus parallel computing
summary(microbenchmark(
  l_apply=lapply(1:10, paws, sleep_time=0.01),
  parl_apply=
    parLapply(clus_ter, 1:10, paws, sleep_time=0.01),
  times=10)
)[, c(1, 4, 5)]
# stop R processes over cluster under Windows
stopCluster(clus_ter)
library(parallel)  # load package parallel
# calculate number of available cores
num_cores <- detectCores() - 1
# initialize compute cluster under Windows
clus_ter <- makeCluster(num_cores)
# define function that pauses execution
paws <- function(x, sleep_time) {
  Sys.sleep(sleep_time)
  x
}  # end paws
# compare speed of lapply with parallel computing
iter_ations <- 3:10
compute_times <- sapply(iter_ations,
  function(max_iterations, sleep_time) {
    out_put <- summary(microbenchmark(
lapply=lapply(1:max_iterations, paws,
              sleep_time=sleep_time),
parallel=parLapply(clus_ter, 1:max_iterations,
        paws, sleep_time=sleep_time),
times=10))[, c(1, 4)]
    structure(out_put[, 2],
        names=as.vector(out_put[, 1]))
    }, sleep_time=0.01)
compute_times <- t(compute_times)
rownames(compute_times) <- iter_ations
library(parallel)  # load package parallel
plot(x=rownames(compute_times),
     y=compute_times[, "lapply"],
     type="l", lwd=2, col="blue",
     main="Compute times",
     xlab="number of iterations in loop", ylab="",
     ylim=c(0, max(compute_times[, "lapply"])))
lines(x=rownames(compute_times),
y=compute_times[, "parallel"], lwd=2, col="green")
legend(x="topleft", legend=colnames(compute_times),
 inset=0.1, cex=1.0, bg="white",
 lwd=2, lty=c(1, 1), col=c("blue", "green"))
library(parallel)  # load package parallel
# calculate number of available cores
num_cores <- detectCores() - 1
# initialize compute cluster under Windows
clus_ter <- makeCluster(num_cores)
# define large matrix
mat_rix <- matrix(rnorm(7*10^5), ncol=7)
# define aggregation function over column of matrix
agg_regate <- function(col_umn) {
  out_put <- 0
  for (in_dex in 1:NROW(col_umn))
    out_put <- out_put + col_umn[in_dex]
  out_put
}  # end agg_regate
# perform parallel aggregations over columns of matrix
agg_regations <-
  parCapply(clus_ter, mat_rix, agg_regate)
# compare speed of apply with parallel computing
summary(microbenchmark(
  ap_ply=apply(mat_rix, MARGIN=2, agg_regate),
  parl_apply=
    parCapply(clus_ter, mat_rix, agg_regate),
  times=10)
)[, c(1, 4, 5)]
# stop R processes over cluster under Windows
stopCluster(clus_ter)
library(parallel)  # load package parallel
# calculate number of available cores
num_cores <- detectCores() - 1
# initialize compute cluster under Windows
clus_ter <- makeCluster(num_cores)
ba_se <- 2
# fails because child processes don't know ba_se:
parLapply(clus_ter, 2:4,
    function(exponent) ba_se^exponent)
# ba_se passed to child via dots ... argument:
parLapply(clus_ter, 2:4,
    function(exponent, ba_se) ba_se^exponent,
    ba_se=ba_se)
# ba_se passed to child via clusterExport:
clusterExport(clus_ter, "ba_se")
parLapply(clus_ter, 2:4,
    function(exponent) ba_se^exponent)
# fails because child processes don't know zoo::index():
parSapply(clus_ter, c("VTI", "IEF", "DBC"),
    function(sym_bol)
      NROW(index(get(sym_bol, envir=rutils::env_etf))))
# zoo function referenced using "::" in child process:
parSapply(clus_ter, c("VTI", "IEF", "DBC"),
    function(sym_bol)
      NROW(zoo::index(get(sym_bol, envir=rutils::env_etf))))
# package zoo loaded in child process:
parSapply(clus_ter, c("VTI", "IEF", "DBC"),
    function(sym_bol) {
      stopifnot("package:zoo" %in% search() || require("zoo", quietly=TRUE))
      NROW(index(get(sym_bol, envir=rutils::env_etf)))
    })  # end parSapply
# stop R processes over cluster under Windows
stopCluster(clus_ter)
library(parallel)  # load package parallel
# calculate number of available cores
num_cores <- detectCores() - 1
# initialize compute cluster under Windows
clus_ter <- makeCluster(num_cores)
# set seed for cluster under Windows
# doesn't work: set.seed(1121)
clusterSetRNGStream(clus_ter, 1121)
# perform parallel loop under Windows
out_put <- parLapply(clus_ter, 1:70, rnorm, n=100)
sum(unlist(out_put))
# stop R processes over cluster under Windows
stopCluster(clus_ter)
# perform parallel loop under Mac-OSX or Linux
out_put <- mclapply(1:10, rnorm, mc.cores=num_cores, n=100)
# calculate random default probabilities
num_assets <- 100
default_probs <- runif(num_assets, max=0.2)
mean(default_probs)
# calculate number of defaults
uni_form <- runif(num_assets)
sum(uni_form < default_probs)
# simulate average number of defaults
num_simu <- 1000
de_faults <- numeric(num_simu)
# simulate using for() loop (inefficient way)
for (i in 1:num_simu) {  # perform loop
  uni_form <- runif(num_assets)
  de_faults[i] <- sum(uni_form < default_probs)
}  # end for
# calculate average number of defaults
mean(de_faults)
# simulate using vectorized functions  (efficient way)
uni_form <- matrix(runif(num_simu*num_assets),
             ncol=num_simu)
sum(uni_form < default_probs)/num_simu
# plot Standard Normal distribution
curve(expr=dnorm(x),
type="l", xlim=c(-4, 4),
xlab="asset value", ylab="", lwd=2,
col="blue", main="Distribution of Asset Values")
abline(v=qnorm(0.025), col="red", lwd=2)
text(x=qnorm(0.025)-0.1, y=0.15,
 labels="default threshold",
 lwd=2, srt=90, pos=3)
# define correlation parameters
rh_o <- 0.2
rho_sqrt <- sqrt(rh_o) ; rho_sqrtm <- sqrt(1-rh_o)
num_assets <- 5 ; num_simu <- 10000
# calculate vector of systematic factors
system_atic <- rnorm(num_simu)
# simulate asset values using vectorized functions (efficient way)
asset_values <- rho_sqrt*system_atic +
  rho_sqrtm*rnorm(num_simu*num_assets)
dim(asset_values) <- c(num_simu, num_assets)
# calculate correlations between asset values
cor(asset_values)
# simulate asset values using for() loop (inefficient way)
# allocate matrix of assets
asset_values <- matrix(nr=num_simu, nc=num_assets)
# simulate asset values using for() loop
for (i in 1:num_simu) {  # perform loop
  asset_values[i, ] <-
    rho_sqrt*system_atic[i] +
    rho_sqrtm*rnorm(num_assets)
}  # end for
cor(asset_values)
# benchmark the speed of the two methods
library(microbenchmark)
summary(microbenchmark(
  for_loop={for (i in 1:num_simu) {
    rho_sqrt*system_atic[i] +
    rho_sqrtm*rnorm(num_assets)}},
  vector_ized={rho_sqrt*system_atic +
        rho_sqrtm*rnorm(num_simu*num_assets)},
  times=10))[, c(1, 4, 5)]
# calculate random default probabilities
num_assets <- 5
default_probs <- runif(num_assets, max=0.2)
mean(default_probs)
# calculate default thresholds
default_thresh <- qnorm(default_probs)
# calculate number of defaults using vectorized functions (efficient way)
# calculate vector of number of defaults
dim(default_thresh)
dim(asset_values)
de_faults <-
  colSums(t(asset_values < default_thresh))
de_faults / num_simu
aa <- matrix(1:20,nr=4,nc=5)
aa
bb <- c(1,6,3,4)
bb
aa>bb
t(aa)>bb ##so that the bb is placed vertically and then compared with the aa.
default_probs
# calculate number of defaults using for() loop (inefficient way)
# allocate matrix of de_faults
de_faults <- matrix(nr=num_simu, nc=num_assets)
# simulate asset values using for() loop
for (i in 1:num_simu) {  # perform loop
  de_faults[i, ] <-
    (asset_values[i, ] < default_thresh)
}  # end for
colSums(de_faults) / num_simu
default_probs
# calculate correlations between defaults
cor(de_faults)
# define default probabilities
num_assets <- 2
default_prob <- 0.2
default_thresh <- qnorm(default_prob)
# define correlation parameters
rh_o <- 0.2
rho_sqrt <- sqrt(rh_o) ; rho_sqrtm <- sqrt(1-rh_o)
# calculate vector of systematic factors
num_simu <- 1000
system_atic <- rnorm(num_simu)
# simulate asset values using vectorized functions
asset_values <- rho_sqrt*system_atic +
  rho_sqrtm*rnorm(num_simu*num_assets)
dim(asset_values) <- c(num_simu, num_assets)
# calculate number of defaults using vectorized functions
de_faults <- t(t(asset_values) < default_thresh)
# calculate correlations between defaults
cor(de_faults)
# calculate averaage number of defaults and compare to default_prob
colSums(de_faults) / num_simu
default_prob
# define cumulative default probability function
def_prob <- function(x, def_thresh=qnorm(0.1), rh_o=0.1)
  pnorm((sqrt(1-rh_o)*qnorm(x) - def_thresh)/sqrt(rh_o))
def_prob(x=0.2, def_thresh=qnorm(0.2), rh_o=0.2)
# plot cumulative default probability function
curve(expr=def_prob(x, def_thresh=qnorm(0.4), rh_o=0.05),
xlim=c(0, 0.999), lwd=3,
xlab="percent default", ylab="probability",
col="green", main="Cumulative Default Probabilities")
# plot default distribution with higher correlation
curve(expr=def_prob(x, def_thresh=qnorm(0.4), rh_o=0.2),
xlim=c(0, 0.999), add=TRUE, lwd=3,
col="blue", main="")
# add legend
legend(x="topleft",
 legend=c("high correlation", "low correlation"),
 title=NULL, inset=0.05, cex=0.8, bg="white",
 bty="n", lwd=6, lty=c(1, 1), col=c("blue", "green"))
# add unconditional default probability
abline(v=0.4, col="red", lwd=3)
text(x=0.4, y=0.0,
 labels="default probability",
 lwd=2, srt=90, pos=4)
# define default probability density function
vasi_cek <- function(x, def_thresh=-2, rh_o=0.1)
  sqrt((1-rh_o)/rh_o)*exp(-(sqrt(1-rh_o)*qnorm(x) -
  def_thresh)^2/(2*rh_o) + qnorm(x)^2/2)
vasi_cek(0.03, def_thresh=qnorm(0.025), rh_o=0.1)
# plot probability distribution of defaults
curve(expr=vasi_cek(x, def_thresh=qnorm(0.025), rh_o=0.02),
xlim=c(0, 0.1), lwd=3,
xlab="percentage of defaults", ylab="density",
col="green", main="Distribution of Defaults")
# plot default distribution with higher correlation
curve(expr=vasi_cek(x, def_thresh=qnorm(0.025), rh_o=0.1),
xlab="default percentage", ylab="",
add=TRUE, lwd=3, col="blue", main="")
# add legend
legend(x="topright",
 legend=c("high correlation", "low correlation"),
 title=NULL, inset=0.05, cex=0.8, bg="white",
 bty="n", lwd=6, lty=c(1, 1), col=c("blue", "green"))
# add unconditional default probability
abline(v=0.025, col="red", lwd=3)
text(x=0.023, y=8,
 labels="default probability",
 lwd=2, srt=90, pos=3)
# plot default distribution with low correlation
curve(expr=vasi_cek(x, def_thresh=qnorm(0.1), rh_o=0.01),
xlab="default percentage", ylab="", lwd=2,
col="green", main="Distribution of Defaults")
# plot default distribution with high correlation
curve(expr=vasi_cek(x, def_thresh=qnorm(0.1), rh_o=0.99),
xlab="percentage of defaults", ylab="density",
add=TRUE, lwd=2, n=10001, col="blue", main="")
# add legend
legend(x="top",
 legend=c("high correlation", "low correlation"),
 title=NULL, inset=0.1, cex=0.8, bg="white",
 bty="n", lwd=6, lty=c(1, 1), col=c("blue", "green"))
# add unconditional default probability
abline(v=0.1, col="red", lwd=2)
text(x=0.1, y=10, lwd=2, pos=4,
 labels="default probability")
# define Vasicek loss distribution density function
portf_loss <- function(x, def_thresh=-2, rh_o=0.1, l_gd=0.4)
  sqrt((1-rh_o)/rh_o)*exp(-(sqrt(1-rh_o)*qnorm(x/l_gd) - def_thresh)^2/(2*rh_o) + qnorm(x/l_gd)^2/2)/l_gd
integrate(portf_loss, low=0, up=0.3,
  def_thresh=-2, rh_o=0.1, l_gd=0.4)
# plot probability distribution of losses
curve(expr=portf_loss(x, def_thresh=qnorm(0.06), rh_o=0.1),
type="l", xlim=c(0, 0.06),
xlab="loss percentage", ylab="density", lwd=3,
col="orange", main="Distribution of Losses")
# add line for expected loss
abline(v=0.02, col="red", lwd=3)
text(x=0.02-0.001, y=10, labels="expected loss",
 lwd=2, srt=90, pos=3)
# add lines for unexpected loss
abline(v=0.04, col="blue", lwd=3)
arrows(x0=0.02, y0=35, x1=0.04, y1=35,
 code=3, lwd=3, cex=0.5)
text(x=0.03, y=36, labels="unexpected loss",
     lwd=2, pos=3)
# add lines for VaR
abline(v=0.055, col="red", lwd=3)
arrows(x0=0.0, y0=25, x1=0.055, y1=25,
 code=3, lwd=3, cex=0.5)
text(x=0.03, y=26, labels="VaR", lwd=2, pos=3)
text(x=0.055-0.001, y=10, labels="VaR",
 lwd=2, srt=90, pos=3)
# plot probability distribution of losses
curve(expr=portf_loss(x, def_thresh=qnorm(0.1), rh_o=0.1),
type="l", xlim=c(0, 0.06),
xlab="loss percentage", ylab="density", lwd=3,
col="orange", main="Conditional Value at Risk")
# add line for expected loss
abline(v=0.02, col="red", lwd=3)
text(x=0.02-0.001, y=10, labels="expected loss",
 lwd=2, srt=90, pos=3)
# add lines for VaR
abline(v=0.04, col="red", lwd=3)
text(x=0.04-0.001, y=10, labels="VaR",
 lwd=2, srt=90, pos=3)
# add shading for CVaR
va_r <- 0.04; var_max <- 0.07
var_s <- seq(va_r, var_max, length=100)
dens_ity <- sapply(var_s, portf_loss,
  def_thresh=qnorm(0.1), rh_o=0.1)
# draw shaded polygon
polygon(c(va_r, var_s, var_max),
  c(-1, dens_ity, -1), col="red", border=NA)
text(x=0.045, y=0, labels="CVaR", lwd=2, pos=3)
# VaR (quantile of the loss distribution)
var_func <- function(x, def_thresh=qnorm(0.1), rh_o=0.1, l_gd=0.4)
  l_gd*pnorm((sqrt(rh_o)*qnorm(x) + def_thresh)/sqrt(1-rh_o))
var_func(x=0.99, def_thresh=qnorm(0.1), rh_o=0.2, l_gd=0.4)
# plot VaR
curve(expr=var_func(x, def_thresh=qnorm(0.1), rh_o=0.1, l_gd=0.4),
type="l", xlim=c(0, 0.999),
xlab="confidence level", ylab="VaR", lwd=3,
col="orange", main="VaR versus Confidence Level")
# add line for expected loss
abline(h=0.04, col="red", lwd=3)
text(x=0.2, y=0.04, labels="expected loss",
     lwd=2, pos=3)
# integrate portf_loss() over full range
integrate(portf_loss, low=0.0, up=0.3,
    def_thresh=qnorm(0.1), rh_o=0.1, l_gd=0.4)
# calculate expected losses using portf_loss()
integrate(function(x, ...) x*portf_loss(x, ...),
    low=0.0, up=0.3,
    def_thresh=qnorm(0.1), rh_o=0.1, l_gd=0.4)
# calculate confidence levels corresponding to VaR values
var_s <- seq(0.07, 0.12, 0.001)
conf_levels <- sapply(var_s, function(va_r, ...) {
  integrate(portf_loss, low=va_r, up=0.3, ...)
}, def_thresh=qnorm(0.1), rh_o=0.1, l_gd=0.4)  # end sapply
conf_levels <- cbind(as.numeric(t(conf_levels)[, 1]), var_s)
colnames(conf_levels) <- c("conf_levels", "VaRs")
# calculate 95% confidence level VaR value
conf_levels[
  match(TRUE, conf_levels[, "conf_levels"] < 0.05), "VaRs"]
plot(x=1-conf_levels[, "conf_levels"],
     y=conf_levels[, "VaRs"], lwd=2,
     xlab="conf_levels", ylab="VaRs",
     t="l", main="VaR values and confidence levels")
# calculate CVaR values
cvar_s <- sapply(var_s, function(va_r, ...) {
  integrate(function(x, ...) x*portf_loss(x, ...),
      low=va_r, up=0.3, ...)
}, def_thresh=qnorm(0.1), rh_o=0.1, l_gd=0.4)  # end sapply
conf_levels <- cbind(conf_levels, as.numeric(t(cvar_s)[, 1]))
colnames(conf_levels)[3] <- "CVaRs"
# divide CVaR by confidence level
conf_levels[, "CVaRs"] <-
  conf_levels[, "CVaRs"]/conf_levels[, "conf_levels"]
# calculate 95% confidence level CVaR value
conf_levels[match(TRUE,
  conf_levels[, "conf_levels"] < 0.05), "CVaRs"]
# plot CVaRs
plot(x=1-conf_levels[, "conf_levels"],
     y=conf_levels[, "CVaRs"],
     t="l", col="red", lwd=2,
     ylim=range(conf_levels[, c("VaRs", "CVaRs")]),
     xlab="conf_levels", ylab="CVaRs",
     main="CVaR values and confidence levels")
# add VaRs
lines(x=1-conf_levels[, "conf_levels"],
y=conf_levels[, "VaRs"], lwd=2)
# add legend
legend(x="topleft", legend=c("CVaRs", "VaRs"),
 title="default probability = 10%
correlation = 10%
loss given default = 40%",
 inset=0.1, cex=0.8, bg="white", bty="n",
 lwd=6, lty=c(1, 1), col=c("red", "black"))
# Define model parameters
num_assets <- 300
num_simu <- 1000
l_gd <- 0.4
# define correlation parameters
rh_o <- 0.2
rho_sqrt <- sqrt(rh_o)
rho_sqrtm <- sqrt(1-rh_o)
# calculate default probabilities and thresholds
set.seed(1121)
default_probs <- runif(num_assets, max=0.2)
default_thresh <- qnorm(default_probs)
# calculate vector of systematic factors
system_atic <- rnorm(num_simu)
# simulate losses under Vasicek model
asset_values <- matrix(rnorm(num_simu*num_assets), ncol=num_simu)
asset_values <- t(rho_sqrt*system_atic + t(rho_sqrtm*asset_values))
loss_es <-
  l_gd*colSums(asset_values < default_thresh)/num_assets
# calculate VaRs
conf_levels <- seq(0.93, 0.99, 0.01)
var_s <- quantile(loss_es, probs=conf_levels)
plot(x=conf_levels, y=var_s, t="l", lwd=2,
     main="Simulated VaR and confidence levels")
# calculate CVaRs
cvar_s <- sapply(var_s, function(va_r) {
  mean(loss_es[loss_es>va_r])
})  # end sapply
cvar_s <- cbind(cvar_s, var_s)
# alternative CVaR calculation using frequency table
# first calculate frequency table of loss_es
table_losses <- table(loss_es)/num_simu
# calculate CVaRs from frequency table
cvar_s <- sapply(var_s, function(va_r) {
  tai_l <- table_losses[names(table_losses) > va_r]
  tai_l %*% as.numeric(names(tai_l)) / sum(tai_l)
})  # end sapply
# plot CVaRs
plot(x=rownames(cvar_s), y=cvar_s[, "cvar_s"],
     t="l", col="red", lwd=2,
     ylim=range(cvar_s),
     xlab="conf_levels", ylab="CVaRs",
     main="Simulated CVaR and confidence levels")
# add VaRs
lines(x=rownames(cvar_s), y=cvar_s[, "var_s"], lwd=2)
# add legend
legend(x="topleft", legend=c("CVaRs", "VaRs"), bty="n",
 title=NULL, inset=0.05, cex=0.8, bg="white",
 lwd=6, lty=c(1, 1), col=c("red", "black"))
calc_var <- function(default_thresh,
               l_gd=0.6,
               rho_sqrt,
               rho_sqrtm,
               num_simu=1000,
               conf_levels=seq(0.93, 0.99, 0.01)) {
  # Define model parameters
  num_assets <- NROW(default_thresh)
  # Simulate losses under Vasicek model
  system_atic <- rnorm(num_simu)
  asset_values <- matrix(rnorm(num_simu*num_assets), ncol=num_simu)
  asset_values <- t(rho_sqrt*system_atic + t(rho_sqrtm*asset_values))
  loss_es <- l_gd*colSums(asset_values < default_thresh)/num_assets
  # Calculate VaRs and CVaRs
  var_s <- quantile(loss_es, probs=conf_levels)
  cvar_s <- sapply(var_s, function(va_r) {
    mean(loss_es[loss_es>va_r])
  })  # end sapply
  names(cvar_s) <- names(var_s)
  c(var_s, cvar_s)
}  # end calc_var
# define number of bootstrap simulations
num_boot <- 500
num_assets <- NROW(default_probs)
# perform bootstrap of calc_var
set.seed(1121)
boot_strap <- sapply(rep(l_gd, num_boot),
  calc_var,
  default_thresh=qnorm(default_probs),
  rho_sqrt=rho_sqrt,
  rho_sqrtm=rho_sqrtm,
  num_simu=num_simu,
  conf_levels=conf_levels)  # end sapply
boot_strap <- t(boot_strap)
# calculate vectors of standard errors of VaR and CVaR from boot_strap data
std_error_var <- apply(boot_strap[, 1:7], MARGIN=2,
    function(x) c(mean=mean(x), sd=sd(x)))
std_error_cvar <- apply(boot_strap[, 8:14], MARGIN=2,
    function(x) c(mean=mean(x), sd=sd(x)))
# scale the standard errors of VaRs and CVaRs
std_error_var[2, ] <- std_error_var[2, ]/std_error_var[1, ]
std_error_cvar[2, ] <- std_error_cvar[2, ]/std_error_cvar[1, ]
# plot the standard errors of VaRs and CVaRs
plot(x=colnames(std_error_cvar),
  y=std_error_cvar[2, ], t="l", col="red", lwd=2,
  ylim=range(c(std_error_var[2, ], std_error_cvar[2, ])),
  xlab="conf_levels", ylab="CVaRs",
  main="Scaled standard errors of CVaR and VaR")
lines(x=colnames(std_error_var), y=std_error_var[2, ], lwd=2)
legend(x="topleft", legend=c("CVaRs", "VaRs"), bty="n",
 title=NULL, inset=0.05, cex=0.8, bg="white",
 lwd=6, lty=c(1, 1), col=c("red", "black"))
library(parallel)  # load package parallel
num_cores <- detectCores() - 1  # number of cores
clus_ter <- makeCluster(num_cores)  # initialize compute cluster
# perform bootstrap of calc_var for Windows
set.seed(1121)
boot_strap <- parLapply(clus_ter, rep(l_gd, num_boot),
  fun=calc_var, default_probs=default_probs,
  rh_o=rh_o, num_simu=num_simu,
  conf_levels=conf_levels)  # end parLapply
# bootstrap under Mac-OSX or Linux
boot_strap <- mclapply(rep(l_gd, num_boot),
  FUN=calc_var, default_probs=default_probs,
  rh_o=rh_o, num_simu=num_simu,
  conf_levels=conf_levels)  # end mclapply
boot_strap <- rutils::do_call(rbind, boot_strap)
stopCluster(clus_ter)  # stop R processes over cluster
# calculate vectors of standard errors of VaR and CVaR from boot_strap data
std_error_var <- apply(boot_strap[, 1:7], MARGIN=2,
    function(x) c(mean=mean(x), sd=sd(x)))
std_error_cvar <- apply(boot_strap[, 8:14], MARGIN=2,
    function(x) c(mean=mean(x), sd=sd(x)))
# scale the standard errors of VaRs and CVaRs
std_error_var[2, ] <- std_error_var[2, ]/std_error_var[1, ]
std_error_cvar[2, ] <- std_error_cvar[2, ]/std_error_cvar[1, ]
# plot the standard errors of VaRs and CVaRs
plot(x=colnames(std_error_cvar),
  y=std_error_cvar[2, ], t="l", col="red", lwd=2,
  ylim=range(c(std_error_var[2, ], std_error_cvar[2, ])),
  xlab="conf_levels", ylab="CVaRs",
  main="Scaled standard errors of CVaR and VaR")
lines(x=colnames(std_error_var), y=std_error_var[2, ], lwd=2)
legend(x="topleft", legend=c("CVaRs", "VaRs"), bty="n",
 title=NULL, inset=0.05, cex=0.8, bg="white",
 lwd=6, lty=c(1, 1), col=c("red", "black"))
