#################################
### FRE6871 Test #2 February 5, 2018
#################################
# Max score 150pts

# Please write in this file the R code needed to perform the tasks below, 
# rename the file to your_name_test2.R
# and upload it to NYU Classes,


############## Part I
# Summary: Create functions which calculate the indices
# of the TRUE elements of a Boolean vector, similar to
# the function which().
# Implement the functions using several different methods,
# first using a for() loop over vectors, second using an
# lapply() loop, third using vectorized operations over 
# the vectors.

# 1. (20pts) First method:
# Create a function called which_for(), that produces
# the same result as function which(), when applied to
# Boolean vectors.
# You must perform a for() loop.
# hint: When performing the for() loop, you can first
# allocate an integer vector, and then perform a
# for() loop to populate it with the index values.
# You can use functions integer(), seq_along(), and c().

### write your code here
which_for <- function(boolean_vector)
{
  order <- integer()
  for(i in seq_along(boolean_vector))
  {
    if(boolean_vector[i]==TRUE)
    {order <- c(order,i)} #like pop in the vector,leijia.
  }
  return(order)
}
# Verify that which_for() applied to vec_tor produces 
# the same result as function which().
# Use function all.equal().

vec_tor <- sample(c(TRUE, FALSE), size=10, replace=TRUE)

### write your code here
vec_tor <- sample(c(TRUE, FALSE), size=10, replace=TRUE)
all.equal(which(vec_tor),
          which_for(vec_tor))

# 2. (20pts) Second method:
# Create a function called which_apply(), that produces
# the same result as function which(), when applied to
# Boolean vectors.
# You must perform an lapply() loop.
# You can use functions lapply(), seq_along(),
# and unlist().

### write your code here
which_apply <- function(vec_tor) {
  order <- seq_along(vec_tor)
  unlist(lapply(order, function(x) if(vec_tor[x]) x))}
# Verify that which_apply() applied to vec_tor produces 
# the same result as function which().
# Use function all.equal().

vec_tor <- sample(c(TRUE, FALSE), size=10, replace=TRUE)

### write your code here
vec_tor <- sample(c(TRUE, FALSE), size=10, replace=TRUE)
all.equal(which(vec_tor),
          which_apply(vec_tor))

# 3. (20pts) Third method:
# Create a function called which_vec(), that produces
# the same result as function which(), when applied to
# Boolean vectors.
# You cannot perform any type of loop, only vectorized functions.
# hint: You can use functions NROW() or seq_along(),
# and then perform vector subsetting.

### write your code here
which_vac <- function(x)
{ seq_along(x)[x]}

# Verify that which_vec() applied to vec_tor produces 
# the same result as function which().
# Use function all.equal().

vec_tor <- sample(c(TRUE, FALSE), size=10, replace=TRUE)

### write your code here
vec_tor <- sample(c(TRUE, FALSE), size=10, replace=TRUE)
all.equal(which_vac(vec_tor),which(vec_tor))

# 4. (10pts) Benchmark the speed of which_for(), which_apply(), 
# and which_vec(), versus which(), using the function microbenchmark().

library(microbenchmark)

### write your code here
library(microbenchmark)
which_for <- which_for(vec_tor)
which_vec <- which_vac(vec_tor)
which_apply <- which_apply(vec_tor)
which_r <- which(vec_tor)
summary(microbenchmark(
  which_for,
  which_apply,
  which_vec,
  which_r,
  times=10))[, c(1, 4, 5)] 

# You should get a result similar to this (only 
# the ratios of the times matter):
#          expr     mean median
# 1   which_for  3744.35   3422
# 4 which_apply 30497.25  29813
# 5   which_vec  1696.68   1467
# 6     which_r  2283.00   1467



############## Part II
# Summary: Create a function called find_interval(),
# that reproduces the function findInterval(), using 
# two different methods.

# 1. (20pts)
# find_interval() should accept a vector argument called vec_tor,
# containing numeric values, which should each be classified into 
# intervals, according to the break_points.
# It should also accept a vector of breakpoints called break_points,
# which determines the intervals.
# find_interval() should return an integer vector of length equal
# to vec_tor, specifying the intervals to which the numeric values
# contained in vec_tor belong.
# hint: You can perform a for() loop over break_points.
# You can use functions integer(), NROW(), seq_along(),
# logical operators, and a for() loop.

### write your code here
find_interval <- function(x,v)
{
  vec <- sort(v);
  test <- integer(length(x));
  for(i in 1:length(vec))
  {
    boolean <- x>=v[i]
    test <- test+boolean #one by one add in.
  }
  test 
}
# Call find_interval() and findInterval() on a vector of numbers 
# and a vector of breakpoints, to verify that they produce exactly 
# the same outputs.
# Use function all.equal().

vec_tor <- 1:10
break_points <- c(3, 5, 7, 9)

### write your code here

find_interval(vec_tor,break_points)
findInterval(vec_tor,break_points)
all.equal(find_interval(vec_tor,break_points),
          findInterval(vec_tor,break_points))

# 2. (20pts)
# Create a function called find_interval2(), that reproduces 
# the function findInterval().
# You cannot perform a for() loop over break_points.
# hint: You can perform a sapply() loop over vec_tor.
# You can use functions match(), is.na(), NROW(),
# logical operators, and an sapply() loop.

### write your code here
find_interval2 <- function(x,v)
{
  vec <- sort(v);
  test <- integer(length(x));
  in_terval <- seq_along(vec);
  rowSums(sapply(in_terval,function(a){boolean <- x>=vec[a];
  test <- test+boolean}))
}

# Call find_interval2() and findInterval() on a vector of numbers 
# and a vector of breakpoints, to verify that they produce exactly 
# the same outputs.
# Use function all.equal().

vec_tor <- 1:10
break_points <- c(3, 5, 7, 9)

### write your code here
find_interval2(vec_tor,break_points)
findInterval(vec_tor,break_points)
all.equal(find_interval2(vec_tor,break_points),
          findInterval(vec_tor,break_points))

# Benchmark the speed of findInterval() versus find_interval(), 
# and find_interval2(), using the function microbenchmark(),

vec_tor <- rnorm(1000)
break_points <- c(-4, -2, 0, 2, 4)
library(microbenchmark)

### write your code here

library(microbenchmark)
findInterval <- findInterval(vec_tor,break_points)
find_interval <- find_interval(vec_tor,break_points)
find_interval2 <- find_interval2(vec_tor,break_points)
summary(microbenchmark(findInterval,find_interval,find_interval2,times=10))[, c(1, 4, 5)] 
# You should get a result similar to this (only 
# the ratios of the times matter):
#             expr      mean    median
# 1   findInterval   16.7153   15.1515
# 2  find_interval   48.2382   47.1625
# 3 find_interval2 1494.6746 1493.0620



############## Part III
# Summary: Extract (find) the matrix column containing its
# largest element.

# 1. (10pts) Create a matrix called mat_rix, with 20 random
# integer elements, and with 4 columns and 5 rows.
# Assign the names: "col1, ..., col4" to the columns of mat_rix.
# You can use functions matrix(), sample(), colnames() and paste0().
# You can't use function c().

set.seed(1121)

### write your code here
mat_rix <- matrix(sample(20),nrow=5,ncol=4)
colnames(mat_rix) <- paste0("col",1:4)
mat_rix
# You should get the following result:
# > mat_rix
#      col1 col2 col3 col4
# [1,]   12   15    3    2
# [2,]    6   18    5   10
# [3,]   13   17   16   11
# [4,]   19   20    7    1
# [5,]    9    4    8   14

# 2. (20pts) Extract (find) the matrix column containing its
# largest element.
# Be sure to produce a column matrix, not a vector.
# You can use functions max() and which() with the "arr.ind"
# argument, and matrix subsetting with the "drop" argument.

### write your code here

a <- integer()
for (i in seq_along(mat_rix[,1]))
{
  b <- max(mat_rix[i,])
  a <- c(a,b)
}
a
m <- matrix(a, ncol=1)
m
colnames(m) <- list("col2")
m

# You should get the following result:
#      col2
# [1,]   15
# [2,]   18
# [3,]   17
# [4,]   20
# [5,]    4


# 3. (10pts) Find the row and column numbers of all the elements 
# of mat_rix that are greater than 15,
# You can use the function which() with the "arr.ind" argument.

### write your code here
which(mat_rix>15,arr.ind=TRUE)
# You should get the following result:
#      row col
# [1,]   4   1
# [2,]   2   2
# [3,]   3   2
# [4,]   4   2
# [5,]   3   3

# Extract the elements of mat_rix that are greater than 15,

### write your code here
Matr_ix <- mat_rix[(mat_rix>15)];
Matr_ix;

# You should get the following result:
# [1] 19 18 17 20 16


# 4. (10pts) 
# Replace the elements of mat_rix that are greater than 15 
# with NA,

### write your code here
mat_rix[(mat_rix>15)] <- NA;
mat_rix;

# Calculate the column means of mat_rix using the function 
# apply(), and omit the NA values.
# hint: pass the parameter "na.rm=TRUE" to function mean().

### write your code here
column_mean <- apply(mat_rix,2,mean,na.rm=TRUE);
column_mean

# You should get the following result:
#  col1  col2  col3  col4 
# 10.00  9.50  5.75  7.60

