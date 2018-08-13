#################################
### FRE6871 Test #4 February 26, 2018
#################################
# Max score 180pts

# Please write in this file the R code needed to perform 
# the tasks below, rename the file to your_name_test4.R
# and upload it to NYU Classes,

############## Part I
# Summary: Estimate the probability of a Brownian motion
# path crossing a barrier by performing multiple simulations.
# Estimate the crossing probabilities for several different
# barrier levels.
# Calculate the standard error of the crossing probability
# estimate using bootstrap simulation.

# Define the simulation parameters:
# number of steps in each simulation
len_gth <- 1000
# barrier level
bar_rier <- 20
# number of simulations
n_simu <- 500


# 1. (20pts) Simulate multiple Brownian paths, by performing
# an sapply() loop starting from 1 to n_simu, and inside the
# loop perform a vectorized simulation of Brownian motion
# crossing the barrier bar_rier.
# The sapply() loop should return a Boolean vector of length
# n_simu (called cross_ed), that should be equal to TRUE if
# the given path did cross the bar_rier, and otherwise it
# should be FALSE.

# hint: you can adapt code from the slide titled
# "Simulating Brownian Motion Using Vectorized Functions".
# hint: you can use an anonymous function that accepts an
# integer argument (the loop count) and returns a Boolean
# value.
# You can compare the simulated path values to the bar_rier 
# level, to determine if at any point the path reached above 
# the bar_rier. The comparison of the path with bar_rier 
# produces a Boolean vector, whose sum is zero only if path 
# never crossed bar_rier, and is greater than zero if they did.
# You can use functions sapply(), sum(), cumsum(), and
# rnorm().

# reset random number generator
set.seed(1121)

### write your code here
simulate_path <- function(len_gth,bar_rier)
{
  pa_th <- cumsum(rnorm(len_gth))
  cro_ss <- which(pa_th > bar_rier)
  if (NROW(cro_ss)>0) {
    return(TRUE)
  } 
  return(FALSE)
}
cross_ed <- sapply(seq(n_simu),function(x,...)
{
  simulate_path(len_gth,bar_rier)
},len_gth=len_gth,bar_rier=bar_rier)
cross_ed

#another way to do this:
set.seed(1121)
cross_ed <- lapply(seq(n_simu),function(x,...)
{
  simulate_path(len_gth,bar_rier)
},len_gth=len_gth,bar_rier=bar_rier) %>% do.call(rbind,.)
mean(cross_ed)

# Calculate the probability of crossing the bar_rier as the
# number of TRUE elements of cross_ed, divided by n_simu.

### write your code here
sum(cross_ed)/n_simu

# You should get the following result:
# [1] 0.516


# 2. (20pts) Perform the same simulation as in p.1 but
# without using an sapply() loop, only using vectorized
# functions.
# Create a matrix of normal random numbers, with dimensions
# equal to n_simu columns by len_gth rows, and call it path_s.
# Apply function colCumsums() from package matrixStats to
# path_s, to calculate the cumulative sums of its columns.
# The columns of path_s represent vectors of simulated paths.
# You can use functions matrix(), colCumsums(), and rnorm().

# load package matrixStats
library(matrixStats)
# reset random number generator
set.seed(1121)

### write your code here
path_s <- matrix(rnorm(len_gth*n_simu),ncol = n_simu)
path_s <- colCumsums(path_s)

# Following the methodology of p.1, compare the simulated
# paths (columns of path_s) to the bar_rier level, and produce a
# Boolean matrix.
# Sum up the columns of the Boolean matrix to determine the
# the paths that crossed the bar_rier, and call the resulting
# Boolean vector cross_ed.
# cross_ed should be a vector of length n_simu.
# You can use function colSums().

### write your code here

cross_ed <- colSums(path_s>bar_rier_many)>0
cross_ed
# Calculate the probability of crossing the bar_rier as the
# number of TRUE elements of cross_ed, divided by n_simu.

### write your code here
mean(cross_ed)

# You should get the following result:
# [1] 0.516


# 3. (20pts) Estimate the crossing probabilities for
# a vector of different barrier levels.
# Create a named numeric vector called bar_riers with
# values from=5, to=60, by=5.
# You can use functions seq(), structure(), paste0(),
# and names().

### write your code here
bar_riers <- seq(from=5, to=60, by=5)
names(bar_riers) <- paste0("level = ", seq(from=5, to=60, by=5))
bar_riers


# You should get the following result:
# bar_riers
# level=5 level=10 level=15 level=20 level=25 level=30 level=35 level=40 
#       5       10       15       20       25       30       35       40 
# level=45 level=50 level=55 level=60 
#       45       50       55       60 

# Perform an sapply() loop over bar_riers.
# Inside the loop calculate the probabilities of crossing
# the bar_rier, and call the resulting vector prob_s.
# To receive full credit you shouldn't recalculate the
# path_s for each loop (bar_rier), but instead use the
# path_s already calculated in p.2.
# You can use functions sapply(), sum(), colSums(), and
# an anonymous function.

### write your code here
prob_s <- sapply(bar_riers, function(x)
{bar_rier_many <- rep(x,n_simu)
cross_ed <- colSums(path_s>bar_rier_many)>0
mean(cross_ed)
})
prob_s


# Create a scatterplot of prob_s versus bar_riers.
# You can use functions plot() and title().

### write your code here
plot(bar_riers,prob_s,xlab="barrier levels", ylab="probabilities",type = "l",col="blue")
title(main="Barrier Crossing Probabilities", line=-1)


# Your plot should be similar to simu_barrier_probs.png


# 4. (20pts) Calculate the standard error of the crossing
# probability estimate using parallel bootstrap simulation.
# Perform a loop over 1:n_boots, and inside the loop 
# calculate the crossing probabilities by repeating the 
# calculation from p.2.
# You must resample from the vector path_s defined below.
# You cannot use function rnorm() inside the loop.
# Instead you should use function sample.int() to resample 
# from the vector path_s.
# 
# Call the resulting vector of probabilities prob_s (it
# should have length n_boots).
# Calculate the standard error of the crossing
# probability estimate from the vector prob_s.
# Calculate the standard error of the mean(prob_s) of 
# vector prob_s.
# You can use functions sd(), mean(), detectCores(), 
# makeCluster(), clusterSetRNGStream(), clusterExport(), 
# sample.int(), and either parSapply() or mclapply().

# Define the simulation parameters:
len_gth <- 1000
# barrier level
bar_rier <- 20
# number of simulations
n_simu <- 500
# number of bootstrap simulations
n_boots <- 1e2
# define the vector path_s
set.seed(1121)
path_s <- rnorm(n_simu*len_gth)

library(parallel)  # load package parallel

num_cores <- detectCores() - 1
clus_ter <- makeCluster(num_cores)
clusterExport(clus_ter, varlist=c("len_gth","bar_rier","n_simu","path_s"))
clusterSetRNGStream(clus_ter, 1121)

prob_s <- parSapply(clus_ter,seq(n_boots),
                    function(x) {
                      path_s <- matrixStats::colCumsums(matrix(path_s, nc=n_simu))
                      cross_ed <- colSums(path_s[,sample.int(n_simu, replace=TRUE)]>bar_rier)>0
                      sum(cross_ed)/n_simu
                    }
)
prob_s <- rutils::do_call(rbind, prob_s)
stopCluster(clus_ter)

# Calculate the standard error

sd(prob_s)

# You should get a result similar to this:
# [1] 0.02320616

# Calculate the standard error of mean(prob_s):

### write your code here
sd(prob_s)/sqrt(length(prob_s))

# 5. (20pts) Plot a histogram of prob_s.
# Add a curve with the Normal density, with the same mean and sd.
# Add a legend.
# You can use functions hist(), lines(), density(), curve(), and
# legend().

x11()

### write your code here

# Your plot should be similar to simu_barrier_probs_hist.png
hist(prob_s, col="lightgrey", xlab="Option prices", freq=FALSE,
     main="Bootstrap of Barrier Option Simulation")
lines(density(prob_s, adjust=3), type="l", lwd=2, col="blue")
curve(expr=dnorm(x, mean=mean(prob_s), sd=sd(prob_s)),
      add=TRUE, type="l", lwd=2, col="red")
# add legend
legend(x="topright", legend=c("density", "Normal density"),
       title=paste0("mean price=", format(mean(prob_s), digits=3),
                    "\n", "sd=", format(sd(prob_s), digits=3)),
       inset=0.05, cex=1.0, bg="white", bty="n",
       lwd=6, lty=1, col=c("blue", "red"))



############## Part II
# Summary: Calculate the standard errors of regression
# coefficients using parallel bootstrap simulation.

# 1. (20pts) 
# Specify de_sign as a matrix with two columns, containing 
# the returns of the XLF and VTI ETFs, as follows (run all 
# the below code):

library(HighFreq)
de_sign <- rutils::env_etf$re_turns[, c("XLF", "VTI")]
de_sign <- na.omit(de_sign)
de_sign <- coredata(de_sign)
len_gth <- NROW(de_sign)

# Create a formula object called reg_formula, from the 
# column names of de_sign.
# You cannot just type the formula!
# You can use functions colnames(), as.formula(), 
# and paste() with the "collapse" arguments.

### write your code here
reg_formula <- as.formula(paste(colnames(de_sign),collapse = "~"))
reg_formula

# You should get the following result:
# > reg_formula
# XLF ~ VTI

# Perform the regression specified by reg_formula and de_sign.
# Use function lm() as follows:

reg_model <- lm(reg_formula, data=as.data.frame(de_sign))

# Extract the regression coefficients from the list reg_model.

### write your code here
reg_model$coefficients

# You should get the following output:
#   (Intercept)           VTI
# -0.0001815245  1.3650531384


# Extract the standard errors of the regression coefficients 
# from the regression summary list produced by function summary().

### write your code here
reg_model_sum <- summary(reg_model)
attributes(reg_model_sum)$names
reg_model_sum$coefficients[,"Std. Error"]


# You should get the following output:
#   (Intercept)          VTI 
#  0.0001518561 0.0128118890 


# 2. (20pts) Calculate the regression coefficients 
# and their standard errors from the matrix de_sign, 
# using matrix algebra.
# You should get the same values as in p.1 above.
# You can use functions t(), sd(), sum(), and colMeans().
# You cannot use functions lm() or summary(). 

# calculate the de-meaned design matrix

### write your code here
de_XLF <- t(t(de_sign) - colMeans(de_sign))[,1]
de_VTI <- t(t(de_sign) - colMeans(de_sign))[,2]

# calculate the regression beta

### write your code here
be_ta <- sum(de_XLF*de_VTI) / sum(de_VTI^2)

# calculate the regression alpha

### write your code here
al_pha <- mean(de_sign[,1]) - be_ta*mean(de_sign[,2])

# calculate the standard error of regression beta

### write your code here
residue <- de_sign[,1] - (al_pha + be_ta*de_sign[,2])
sum_of_residue_square <- sum(residue^2)   #SSE
sum_of_de_meaned_response_square <- sum((de_sign[,2]-mean(de_sign[,2]))^2)
se_of_beta <- sqrt(sum_of_residue_square/((len_gth-2)*sum_of_de_meaned_response_square))

# Calculate the standard error of regression alpha.
# hint: the standard error of alpha is equal to the 
# standard deviation of the mean of the residuals.

### write your code here
se_of_alpha <- se_of_beta*sqrt(sum(de_sign[,2]^2)/len_gth)

# 3. (20pts) Calculate the standard errors of regression
# coefficients using parallel bootstrap simulation 
# and matrix algebra.
# Inside the loop you must resample from the rows of the 
# matrix de_sign, and recalculate the regression coefficients 
# from the resampled de_sign matrix using matrix algebra.
# 
# The bootstrap calculation should produce a matrix called 
# boot_strap, with dimension of 2 rows and num_boot columns. 
# 
# You should use functions parSapply() or mclapply(), 
# and sample.int().
# You can use functions t(), sd(), sum(), colMeans(), 
# rutils::do_call(), rbind(), and clusterExport().
# You cannot use functions lm() or summary(). 

# Define the simulation parameters:
num_boot <- 1e4

# Perform parallel bootstrap under Windows

### write your code here

# OR
# Perform parallel bootstrap under Mac-OSX or Linux

### write your code here
library(parallel)
num_cores <- detectCores() - 1
clus_ter <- makeCluster(num_cores)
clusterExport(clus_ter,
              varlist=c("de_sign", "num_boot","len_gth","n_simu"))
clusterSetRNGStream(clus_ter, 1121)

beta_boot <- mclapply(seq(num_boot),
                      function(x,...){
                        new_design <- de_sign[c(sample.int(n_simu, replace=TRUE)),]
                        de_XLF <- t(t(new_design) - colMeans(new_design))[,1]
                        de_VTI <- t(t(new_design) - colMeans(new_design))[,2]
                        
                        be_ta <- sum(de_XLF*de_VTI) / sum(de_VTI^2)
                        al_pha <- mean(new_design[,1]) - be_ta*mean(new_design[,2])
                        
                        return(c(alpha = al_pha,beta = be_ta))
                      },
                      mc.cores=num_cores)  # end mclapply
beta_boot <- rutils::do_call(cbind, beta_boot)

# Calculate the standard errors of regression alpha 
# and beta from the boot_strap matrix.
# You can use functions apply(), mean(), and sd().

### write your code here
aa <- apply(beta_boot,MARGIN = 1,function(x){
  c(mean = mean(x),std_error = sd(x))
})
aa


# You should get output similar to this:
#                   alpha       beta
# mean      -0.0001818862 1.36518920
# std_error  0.0001480540 0.03515898

# Write in one sentence the reason why the standard 
# error from bootstrap is about 2.7 times bigger 
# than the standard error from lm()?
# 
# Answer: the bootstrape may not have the whole sample, so it might loose some sample or even repeatedly get a same sample, therefore the 
#the sample maybe more scatter and cause larger sd.


# 4. (20pts) Plot a histogram of the boot-strapped betas.
# Add a curve with the Normal density, with its mean equal 
# to the regression beta and its sd equal to the standard 
# error of beta (from p.1 above).
# Add a legend.
# You can use functions hist(), lines(), density(), 
# curve(), and legend().

x11()

### write your code here
hist(beta_boot[2,], col="lightgrey", xlab="Betas", freq=FALSE,
     main="Boot-straped Betas")
lines(density(beta_boot[2,], adjust=3), type="l", lwd=2, col="blue")
curve(expr=dnorm(x, mean=aa[1,2], sd=aa[2,2]),
      add=TRUE, type="l", lwd=2, col="red")
# add legend
legend(x="topright", legend=c("density", "Normal density"),
       title=paste0("mean beta=", format(aa[1,2], digits=3),
                    "\n", "sd=", format(aa[2,2], digits=3)),
       inset=0.05, cex=1.0, bg="white", bty="n",
       lwd=6, lty=1, col=c("blue", "red"))

# Your plot should be similar to boot_betas.png


