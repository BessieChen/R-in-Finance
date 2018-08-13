#################################
### FRE6871 Test #1 January 29, 2018
#################################
# Max score 120pts

# Please write in this file the R code needed to perform the tasks below, 
# rename the file to your_name_test1.R
# and upload it to NYU Classes,


############## Part I
# 1. (10pts) 
# Coerce the vector called mat_rix into a matrix with two columns.
# The vector should be transformed row-wise.
# Perform it in two different ways.
# You can use the functions matrix(), dim(), attributes(), 
# t(), and structure(),

# First method:
mat_rix <- 1:10

### write your code here
mat_rix <- 1:10
mat_rix <- matrix(mat_rix, ncol = 2, byrow = TRUE)
mat_rix

# Second method:
mat_rix <- 1:10

### write your code here
mat_rix <- 1:10
mat_rix <- matrix(mat_rix, nrow = 2, ncol = 5)
mat_rix <- t(mat_rix)
mat_rix

# Assign to mat_rix the column names "col1" and "col2".
# You can use the function colnames().

### write your code here
colnames(mat_rix) <- c("col1","col2")
mat_rix


# You should get the following result:
# > mat_rix
#      col1 col2
# [1,]    1    2
# [2,]    3    4
# [3,]    5    6
# [4,]    7    8
# [5,]    9   10


############## Part II
# Summary: Extract and multiply vector elements.
# Create a vector:

set.seed(1121)
vec_tor <- sample(1:10)

# 1. (10pts) Extract every third element of vec_tor, 
# starting with the first one, and coerce it into 
# a single row matrix called mat_rix.
# You can use functions seq() and matrix(),

### write your code here
set.seed(1121)
vec_tor <- sample(1:10)
ma_trix <- t(matrix(vec_tor[seq(1,length(vec_tor),3)]))
ma_trix

# You should get the following result:
# > mat_rix
#       [,1] [,2] [,3] [,4]
# [1,]    6    9    8    1


# 2. (10pts) 
# Extract a vector of odd index elements of vec_tor
# (first, third, etc), and a vector of even index elements
# (second, fourth, etc.).
# Calculate the scalar ("inner") product of the two vectors.
# The value should be a vector with a single element, not
# a matrix.
# hint: you can use the "%*%" inner multiplication operator,
# and the functions seq(), drop(), or function matrix(),
# with the proper "byrow" argument,

### write your code here
set.seed(1121)
vec_tor <- sample(1:10)
vec_tor1 <- vec_tor[seq(1, length(vec_tor),2)]
vec_tor2 <- vec_tor[seq(2, length(vec_tor),2)]
drop(vec_tor1 %*% vec_tor2)

# You should get the following result:
# [1] 151


############## Part III
# Summary: Perform for() and sapply() loops over a vector,
# and then perform the equivalent vectorized operations
# over the vector, and measure the increase in speed.

# First create a vector of random numbers as follows:

set.seed(1121)
vec_tor <- sample(1:10)


# 1. (20pts) Perform a for() loop to replace those elements
# of vec_tor that are greater than "5" with the number "5".
# You can use functions for() and seq_along(),

### write your code here
set.seed(1121)
vec_tor <- sample(1:10)
for (i in seq_along(vec_tor))
{if(vec_tor[i]>5) vec_tor[i] <- 5 }
vec_tor

# You should get the following result:
# > vec_tor
# [1] 5 3 5 5 4 5 5 2 5 1


# 2. (20pts) Perform exactly the same calculations as 
# in p.1. using an sapply().
# You must use either functions apply() lapply(), or sapply().
# You can also use functions NROW() and seq_along(),
# and an anonymous function.
# This is an example of sapply() with an anonymous function:

sapply(vec_tor, function(x) (x^2))

set.seed(1121)
vec_tor <- sample(1:10)

### write your code here
set.seed(1121)
vec_tor <- sample(1:10)
vec_tor <- sapply(vec_tor, function(x) (if(x>5){ x <- 5} else x))
vec_tor


# 3. (20pts) Perform the same calculations as in p.1,
# but only using vectorized operations (logical
# operators and subsetting).
# You cannot use any for() or apply() loops.

set.seed(1121)
vec_tor <- sample(1:10)

### write your code here
set.seed(1121)
vec_tor <- sample(1:10)
vec_tor[vec_tor>5] <- 5
vec_tor

# 4. (10pts) Benchmark the CPU time used by the code
# from p.2 with the code from p.3, using the function
# microbenchmark().
# Assign the names "s_apply" and "vector_ized" to each method.

library(microbenchmark)

### write your code here
library(microbenchmark)
summary(microbenchmark(
  s_apply=sapply(vec_tor, function(x) (if(x>5){ x <- 5}else x)),
  vector_ized=vec_tor[vec_tor>5] <- 5,
  times=10))[, c(1, 4, 5)] 

# You should get a result similar to this:
#          expr     mean median
# 1     s_apply 500114.9  39099
# 2 vector_ized   3910.4   2933



############## Part IV
# 1. (20pts) Create the character string "y = x1 + x2 - x3 - x4" from 
# the characters "x", "y", "=", "+", "-", and the vectors 1:2 and 3:4, 
# using the functions paste0(), and paste() with the "collapse" and 
# "sep" string arguments. 
# You can also use these characters with spaces around them, say " x ",
# hint: you must call paste0() and paste() several times and nest them. 
# You cannot use strings like "x3 - x4", etc.

### write your code here
paste("y ", paste(paste0( " x", 1:2, collapse = " +"),paste0( " x", 3:4, collapse = " -"), sep = " -"), sep = "=")

