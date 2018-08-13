#################################
### FRE6871 Test #5 March 5, 2018
#################################
# Max score 200pts

# Please write in this file the R code needed to perform 
# the tasks below, rename the file to your_name_test5.R
# and upload it to NYU Classes,

############## Part I
# Summary: Modify the incorrect for() and sapply() loops.

# Create a vector variable:

vec_tor <- 1:5

# 1. (10pts) Modify the following for() loop so 
# that the values of vec_tor are changed to: 
# [1] 1 2 3 3 3
# You can write any for() loop, as long as it
# performs the task.

for (i in vec_tor) {
  if (vec_tor[i]>3)
   3
  else
    NULL
}  # end for


### write your code here
vec_tor <- 1:5
for (i in vec_tor) {
  if (vec_tor[i]>3)
    vec_tor[i]=3
  else
    NULL
}  # end for
vec_tor


# 2. (10pts) Modify the following sapply() loop so 
# that the values of vec_tor are changed to: 
# [1] 1 2 3 3 3
# YOu can write any sapply() loop, as long as it
# performs the task.

vec_tor <- 1:5

sapply(vec_tor, function(x) {
  if (vec_tor[i] > 3)
    vec_tor[i] <- 3
  else
    vec_tor[i] <- NULL
})  # end sapply

### write your code here
vec_tor <- 1:5
vec_tor <- sapply(vec_tor, function(x) {
  if (x > 3)
    3
  else
    x
}) 
vec_tor


############## Part II
# Summary: Calculate the p-values of the t-test, Wilcoxon-test, 
# and KS-test as a function of the mean separation between two 
# samples.

# 1. (20pts) Create two vectors of random data, and 
# a vector of means as follows:

mean_s <- seq(0.04, 0.1, 0.001)
set.seed(1121)
sample_1 <- rnorm(1000)
sample_2 <- rnorm(1000)

# Perform an sapply loop over mean_s, and in the loop
# add the values of mean_s to sample_2.
# Calculate the p-values of the t-test, Wilcoxon-test, 
# and the KS-test for sample_1 and sample_2.

### write your code here
p_values <- t(sapply(mean_s,function(x){
  a1 <- t.test(sample_1, sample_2+x)$p.value
  a2 <- wilcox.test(sample_1, sample_2+x)$p.value
  a3 <- ks.test(sample_1, sample_2+x)$p.value
  test <- c(t_test = a1, wilcox = a2, ks = a3)
  return(test)
}))
tail(p_values)

# You should get the following results:
# > tail(p_values)
#           t_test     wilcox        ks
# [56,] 0.01257580 0.00851545 0.0869054
# [57,] 0.01181424 0.00800208 0.0776215
# [58,] 0.01109390 0.00751139 0.0776215
# [59,] 0.01041287 0.00704131 0.0776215
# [60,] 0.00976932 0.00659628 0.0776215
# [61,] 0.00916148 0.00617379 0.0691903


# 2. (10pts) Plot the p-values.
# You can use functions x11(), plot(), lines(), and legend().

### write your code here
x11()
p_values <- t(p_values)
plot(mean_s,p_values[1,], col="orange", type="l",lwd=1,
     xlab="mean values", ylab="test value")
lines(mean_s,p_values[2,],
      type="l", lwd=2, col="blue")
lines(mean_s,p_values[3,],
      type="l", lwd=2, col="green")
title(main="P_value versus differences in means", line=0.5)
# add legend
legend("topright", title="Tests",
       paste("test", rownames(p_values), sep="="),
       inset=0.05, cex=0.6, lwd=2,
       lty=c(1, 1, 1), col=c("orange","blue","green"))



# Your plot should be similar to test_wilcox_ks_pvals.png



############## Part III
# Summary: Perform logistic regression to forecast
# if a subject is a student or not based on default,
# balance, and income data.
# Use the Default dataset from package ISLR::Default.
# hint: You can adapt the lecture code, and replace
# default status with student status.

# 1. (20pts) Plot a scatterplot of students and non-students
# with balance on x-axis and income on y-axis.
# hint: You can create a Boolean vector called is_student, 
# which is TRUE if student=="Yes".
# Your plot should be similar to islr_student_scatter.png.

library(ISLR)  # load package ISLR
attach(Default)  # attach credit default data

### write your code here

x_lim <- range(balance)
y_lim <- range(income)
# plot data points for non-defaulters
is_student <- (student=="Yes")
plot(income ~ balance,
     main="Income and balance for students and non-students",
     xlim=x_lim, ylim=y_lim,
     data=Default[!is_student, ],
     pch=1, col="blue")
# plot data points for defaulters
points(income ~ balance,
       data=Default[is_student, ],
       pch=1, col="red")
legend(x="topright", bty="n",
                               legend=c("non-students", "students"),
                               col=c("blue", "red"), lty=1, pch=4)

# 2. (10pts) Perform Wilcoxon test for balance, to
# test if the mean of balance for students is 
# different from that of non-students.
# hint: You can subset balance using is_student.

# Wilcoxon test for balance predictor

### write your code here
wilcox.test(balance[is_student], balance[!is_student])

# You should get the following results:
#
# W = 12994000, p-value < 2.2e-16

# Perform the same Wilcoxon test for income predictor.

### write your code here
wilcox.test(income[is_student], income[!is_student])

# You should get the following results:
#
# W = 474440, p-value < 2.2e-16


# 3. (20pts) Create a formula object from the column 
# names of Default.
# You cannot just type the formula!
# You can use functions colnames(), as.formula(), 
# and paste().

### write your code here
colnames(Default)
col_names <- colnames(Default)
for_mula <- as.formula(paste(col_names[2],
                             paste(col_names[-2], collapse="+"), sep=" ~ "))
for_mula
# You should get the following result:
# > for_mula
# student ~ default + balance + income

# Perform logistic regression to forecast if 
# a subject is a student or not based on default,
# balance, and income data.
# You should use function glm().

# fit multifactor logistic regression model

### write your code here
log_it <- glm(for_mula, data=Default,
              family=binomial(logit))
summary(log_it)

# You should get the following results:
# > summary(log_it)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)
# (Intercept)  8.526e+00  2.366e-01  36.029   <2e-16 ***
# defaultYes  -5.699e-01  2.615e-01  -2.179   0.0293 *
# balance      1.024e-03  1.019e-04  10.045   <2e-16 ***
# income      -3.933e-04  9.638e-06 -40.805   <2e-16 ***


# 4. (20pts)
# Assume that the null hypothesis is that the
# subject is a student.
# A type I error is the incorrect rejection of a
# TRUE null hypothesis (i.e. a "false positive").
# That is, a type I error is when a subject is
# a student but it's classified as a non-student.
# A type II error is the incorrect acceptance of
# a FALSE null hypothesis (i.e. a "false negative").
# That is, a type II error is when a subject is
# not a student but it's classified as a student.

# Calculate a confusion matrix based on the 
# fitted.values versus is_student, and call it 
# confu_sion.
# You can use function table().

# define probability threshold
thresh_old <- 0.5

### write your code here
fore_casts <- predict(log_it, type="response")
fore_casts[1:6]
identical(log_it$fitted.values, fore_casts)
# calculate confusion matrix
confu_sion <- table(is_student, (fore_casts>thresh_old))
#another way to calculate confusion matrix
confu_sion <- table(is_student, (log_it$fitted.values>thresh_old))
confu_sion
apply(confu_sion, 1, sum)
sum(is_student)
sum(!is_student)
all.equal(unname(apply(confu_sion, 1, sum)), 
          c(sum(!is_student), sum(is_student)))
# You should get the following results:
# > confu_sion
# is_student  FALSE TRUE
#      FALSE   6690  366
#      TRUE     268 2676
#
# > apply(confu_sion, 1, sum)
# FALSE  TRUE
#  7056  2944
#
# > sum(is_student)
#   [1] 2944
# > sum(!is_student)
#   [1] 7056
# 
# all.equal(unname(apply(confu_sion, 1, sum)), 
#   c(sum(!is_student), sum(is_student)))
#   [1] TRUE

# Scale confu_sion so that its rows add up to 1.
# You can use function rowSums().

### write your code here
row_Sums <- rowSums(confu_sion)
confu_sion <- confu_sion/row_Sums
confu_sion

# You should get the following results:
# > confu_sion
# is_student       FALSE       TRUE
#       FALSE 0.94812925 0.05187075
#       TRUE  0.09103261 0.90896739

# The off-diagonal elements of the confu_sion
# contain the type I error (row 2, column 1), and
# the type II error (row 1, column 2).
type_I_error_rates <- confu_sion[2,1]
type_II_error_rates <- confu_sion[1,2]
type_I_error_rates
type_II_error_rates


# 5. (20pts)
# Create a function called con_fuse(), which
# calculates the confusion matrix, scales its
# rows from counts to frequencies (so rows
# add up to 1), and returns a named vector of
# type I and type II error rates.
# 
# The function con_fuse() should accept multiple
# arguments as follows:
#
# con_fuse <- function(
#   res_ponse, # Boolean vector of TRUE response values
#   fore_casts, # forecast probabilities from logistic regression model
#   thresh_old) # probability threshold

### write your code here
con_fuse <- function(
  res_ponse, # Boolean vector of TRUE response values
  fore_casts, # forecast probabilities from logistic regression model
  thresh_old)
{
  confu_sion <- table(res_ponse, (fore_casts>thresh_old))
  row_Sums <- rowSums(confu_sion)
  confu_sion <- confu_sion/row_Sums
  c(typeI = confu_sion[2,1],
    typeII = confu_sion[1,2])
}

con_fuse(is_student, log_it$fitted.values, thresh_old)

# You should get the following results:
# > con_fuse(is_student, log_it$fitted.values, thresh_old)
#         typeI    typeII
#     0.0910326 0.0518707


# 6. (20pts)
# Define a vector of probability thresholds.

threshold_s <- seq(0.01, 0.99, len=11)^2

# Calculate a two-column matrix called
# error_rates containing the type I and
# type II error rates for the vector threshold_s
# of probability threshold values.
# Use functions sapply() and con_fuse().

### write your code here
col_names <- colnames(Default)
for_mula <- as.formula(paste(col_names[2],
                             paste(col_names[-2], collapse="+"), sep=" ~ "))
log_it <- glm(for_mula, data=Default,
              family=binomial(logit))
fore_casts <- predict(log_it, type="response")
is_student <- (student=="Yes")
error_rates <- t(sapply(threshold_s, function(x)
{con_fuse(res_ponse=is_student, fore_casts=fore_casts,thresh_old=x)}
))
head(error_rates)
#another way:
error_rates <- t(sapply(threshold_s, con_fuse,
                        res_ponse=is_student, fore_casts=fore_casts
))
#rownames(error_rates) <- paste("threshold",threshold_s,sep="=")
head(error_rates)


# You should get the following result:
# > head(error_rates)
#            typeI      typeII
# [1,] 0.000000000 0.757227891
# [2,] 0.000000000 0.308956916
# [3,] 0.001019022 0.198979592
# [4,] 0.004755435 0.146258503
# [5,] 0.012907609 0.111678005
# [6,] 0.029551630 0.090702948


# 7. (10pts) Plot the sum of error_rates columns.
# Your plot should be similar to islr_student_error_rates.png.

### write your code here

plot(threshold_s,error_rates[,1]+error_rates[,2],col="black", type="l",lwd=2,
     main="Sum of typeI plus typeII error rates",
     xlab="thresholds", ylab="sum of error rates")


# 8. (10pts) Find the threshold value which has the 
# lowest sum of type I plus type II error rates.
# You can use functions which.min() and rowSums().

### write your code here
optim_threshold <- threshold_s[which.min(error_rates[,1]+error_rates[,2])]

# You should get the following threshold value:
# > optim_threshold
# [1] 0.357604


# 9. (20pts) Plot an ROC curve: the true positive
# rate (y-axis) as a function of the false positive
# rate (type I error, x-axis).
# The true positive rate is equal to 1 minus the
# type II error rate.
# Add a green point for the threshold value which
# minimizes the sum of type I plus type II error rates.
# Your plot should be similar to islr_student_roc.png.
# You can use functions plot(), points(), and text().

### write your code here
plot(x=error_rates[, 1],
     y=1-error_rates[, 2],
     xlab="false positive rate",
     ylab="true positive rate",
     main="ROC curve for defaults",
     type="l", lwd=2, col="blue")
min_index <- which.min(error_rates[,1]+error_rates[,2])
y <- 1-error_rates[, 2][min_index]
x <- error_rates[, 1][min_index]
points(x,y,col="green")
text(x+0.13,y+0.01,labels=paste("minimum of typeI plus \n typeII error rates",threshold_s[min_index],sep="= \n"))

# detach the Default credit default data
detach(Default)

