#################################
### FRE6871 Test #3 February 19, 2018
#################################
# Max score 160pts

# Please write in this file the R code needed to perform 
# the tasks below, rename the file to your_name_test3.R
# and upload it to NYU Classes,


############## Part I
# Summary: Calculate the default correlation between 
# Boolean vectors representing default events (instead of 
# correlations between assets), under the Vasicek model.
# Plot the default correlation as a function of asset 
# correlation and default probability. 

# Calculate vectors of systematic and idiosyncratic factors:
num_assets <- 2
num_simu <- 1e5
set.seed(1121)
system_atic <- rnorm(num_simu)
idio_syncratic <- matrix(rnorm(num_simu*num_assets), nc=num_assets)


# 1. (20pts)
# Create a named vector of default thresholds and 
# a named vector of asset correlations.
# You can use functions seq(), qnorm(), and names().

default_probs <- seq(0.01, 0.3, 0.05)

### write your code here
default_threshs <- qnorm(default_probs)
rho_s <- seq(from=0.01, to=0.3, by=0.05)
names(default_threshs) <- paste0("dp = ",seq(from=0.01, to=0.3, by=0.05))
names(rho_s) <- paste0("rho = ",seq(from=0.01, to=0.3, by=0.05))
round(default_threshs, 2)
round(rho_s, 2)


# You should get the following results:
# 
# round(default_threshs, 2)
# dp=0.01 dp=0.06 dp=0.11 dp=0.16 dp=0.21 dp=0.26 
#   -2.33   -1.55   -1.23   -0.99   -0.81   -0.64
# 
# round(rho_s, 2)
# rho=0.01 rho=0.06 rho=0.11 rho=0.16 rho=0.21 rho=0.26
#     0.01     0.06     0.11     0.16     0.21     0.26


# There are two assets and between them one default 
# correlation.
# Calculate the default correlation for every pair of 
# values of the asset correlation and the default 
# probability.
# Calculate a matrix of default correlations for the
# vectors of asset correlations and default 
# probabilities. 
# You must perform two nested sapply() loops, first 
# over default_threshs and then over rho_s.
# You should calculate asset_values from system_atic
# and idio_syncratic, without generating random numbers
# inside the sapply() loops.
# 
# You can use functions sapply(), sqrt(), t(), and 
# cor().

### write your code here
correlation_s <- sapply(default_threshs, function(x) {sapply(rho_s, function(y) {cor(t(t(sqrt(y)*system_atic+sqrt(1-y)*idio_syncratic)< x))[1,2]})})
round(correlation_s, 4)


# You should get the following result:
# > round(correlation_s, 4)
#          dp=0.01 dp=0.06 dp=0.11 dp=0.16 dp=0.21 dp=0.26
# rho=0.01 -0.0043 -0.0055 -0.0011  0.0006  0.0028  0.0071
# rho=0.06  0.0009  0.0083  0.0179  0.0252  0.0316  0.0374
# rho=0.11  0.0104  0.0237  0.0370  0.0502  0.0607  0.0680
# rho=0.16  0.0177  0.0390  0.0593  0.0778  0.0920  0.0948
# rho=0.21  0.0247  0.0545  0.0872  0.1063  0.1205  0.1227
# rho=0.26  0.0314  0.0762  0.1122  0.1326  0.1497  0.1552



############## Part II
# Summary: Calculate the standard errors of VaR and CVaR 
# under the Vasicek model, for different numbers of 
# simulations, using parallel bootstrap simulation.

# 1. (20pts) Create a function called calc_var()
# which simulates losses under the Vasicek model,
# and calculates a vector of VaR and CVaR values.
# The function calc_var() should accept the following
# arguments:
#  num_simu - a vector with the numbers of simulations,
#  default_thresh - vector of default thresholds,
#  rho_sqrt, rho_sqrtm - correlation parameters,
#  l_gd - loss given default,
#  conf_level - a single confidence level.
# 
# The function calc_var() should return a vector that 
# is double the length of num_simu, with its first 
# elements equal to the VaR values, and the last elements 
# equal to the CVaR values.
# 
# The function calc_var() should call rnorm() only twice.
# The function calc_var() should not reset the random 
# number generator.
# 
# Hint: Modify the function calc_var() from the lecture 
# notes, so that it accepts a vector of num_simu values, 
# instead of a vector of conf_levels values.

### write your code here
calc_var <- function(num_simu,
                     default_thresh,
                     rho_sqrt,
                     rho_sqrtm,
                     l_gd, 
                     conf_level) 
{
  num_assets <- NROW(default_thresh)
  vec_tor <- rnorm(max(num_simu))
  mat_rix <- matrix(rnorm(max(num_simu)*num_assets), ncol=max(num_simu))
  asset_values <- sapply(num_simu,function(x){t(rho_sqrt*vec_tor[1:x] + t(rho_sqrtm*mat_rix[,c(1:x)]))}) 
  loss_es <- sapply(asset_values,function(x)l_gd*colSums(x < default_thresh)/num_assets)
  var_s <- sapply(loss_es, quantile, probs=conf_level)
  cvar_s <- 
  sapply(seq_along(num_simu),function(i){mean(loss_es[[i]][loss_es[[i]]>=var_s[i]])})
  names(var_s) <- paste0("num_sim = ",num_simu)
  names(cvar_s) <- names(var_s)
  c(var_s, cvar_s)
}

# Define the model parameters
num_assets <- 30
num_simu <- seq(5e1, 2e2, 2e1)
l_gd <- 0.4
# Define correlation parameters
rh_o <- 0.08
rho_sqrt <- sqrt(rh_o)
rho_sqrtm <- sqrt(1-rh_o)
# Calculate default probabilities and thresholds
set.seed(1121)
default_probs <- runif(num_assets, max=0.1)
default_thresh <- qnorm(default_probs)
conf_level <- 0.99

# Run calc_var() as follows:

set.seed(1121)
aa <- calc_var(num_simu=num_simu,
         default_thresh=default_thresh,
         rho_sqrt=rho_sqrt,
         rho_sqrtm=rho_sqrtm,
         l_gd=l_gd, 
         conf_level=conf_level)
aa <- as.vector(aa)
aa

# You should get the following result:
#  [1] 0.09333333 0.09333333 0.09333333 0.09093333 0.08560000 0.08026667
#  [7] 0.07493333 0.06960000 0.09333333 0.09333333 0.09333333 0.09333333
# [13] 0.09333333 0.09333333 0.09333333 0.09333333


# 2. (20pts) Calculate the standard errors of the VaR and 
# CVaR using bootstrap simulation, by performing an sapply() 
# loop over the number of bootstrap simulations, called num_boot.
# 
# The sapply() loop should produce a matrix called boot_strap, 
# of dimensions 2*NROW(num_simu) rows by num_boot columns. 
# Transpose boot_strap to obtain a matrix with num_boot rows by 
# 2*NROW(num_simu) columns. 
# Add column names to boot_strap: "nsim=50", "nsim=70", ...
# You can use functions sapply(), rep(), t(), and colnames().

# Define number of bootstrap simulations.
num_boot <- 500

# Perform bootstrap of calc_var.
set.seed(1121)

### write your code here
boot_strap <- sapply(rep(l_gd, num_boot),
                     calc_var,
                     num_simu=num_simu,
                     default_thresh=default_thresh,
                     rho_sqrt=rho_sqrt,
                     rho_sqrtm=rho_sqrtm,
                     conf_level=0.95)
boot_strap <- t(boot_strap)
colnames(boot_strap) <- paste0("nsim=",rep(num_simu,2))
boot_strap[1:4, 1:4]


# You should get the following result:
# > boot_strap[1:4, 1:4]
#         nsim=50    nsim=70    nsim=90   nsim=110
# [1,] 0.09333333 0.09333333 0.09333333 0.09093333
# [2,] 0.06013333 0.05746667 0.05480000 0.05333333
# [3,] 0.06666667 0.07080000 0.08000000 0.08000000
# [4,] 0.06013333 0.07493333 0.06960000 0.06546667

# Calculate vectors of the standard errors of VaR and 
# CVaR from the boot_strap data, scaled by their means.
# Call them std_error_var and std_error_cvar.
# Then perform an apply() loop over the columns of boot_strap 
# to get the means and standard errors of VaR and CVaR.
# You can use functions NROW(), apply(), mean(), and sd().

### write your code here
std_error_var_my <- apply(boot_strap[, 1:(NCOL(boot_strap)/2)], MARGIN=2,
                       function(x) c(mean(x), sd(x)))
std_error_cvar_my <- apply(boot_strap[, ((NCOL(boot_strap)/2)+1):NCOL(boot_strap)], MARGIN=2,
                        function(x) c(mean(x), sd(x)))
# scale the standard errors of VaRs and CVaRs
std_error_var <- std_error_var_my[2, ]/std_error_var_my[1, ]
std_error_cvar <- std_error_cvar_my[2, ]/std_error_cvar_my[1, ]

round(std_error_var, 5)
round(std_error_cvar, 5)

# You should get the following result:
# > round(std_error_var, 5)
# nsim=50  nsim=70  nsim=90 nsim=110 nsim=130 nsim=150 nsim=170 nsim=190 
# 0.18501  0.16849  0.15625  0.14999  0.13371  0.12643  0.12393  0.12372 
#
# > round(std_error_cvar, 5)
# nsim=50  nsim=70  nsim=90 nsim=110 nsim=130 nsim=150 nsim=170 nsim=190 
# 0.22580  0.21686  0.20539  0.17075  0.16212  0.15509  0.15288  0.14814


# 3. (10pts) Plot std_error_var and std_error_cvar.
# You should use functions plot() and lines().

### write your code here
plot(x=seq(5e1, 2e2, 2e1),
     y=std_error_cvar, t="l", col="red", lwd=2,
     ylim=range(c(std_error_var, std_error_cvar)),
     xlab="number of simulation", ylab="Scaled standard errors",
     main="Scaled standard errors of CVaR and VaR")
lines(x=seq(5e1, 2e2, 2e1), y=std_error_var, lwd=2, col="blue")
legend(x="topright", legend=c("CVaRs", "VaRs"), bty="n",
       title=NULL, inset=0.05, cex=0.8, bg="white",
       lwd=6, lty=c(1, 1), col=c("red", "black"))

# Your plot should be similar to cvar_std_error_boot.png


# 4. (20pts) Perform the bootstrap simulation from p.2 
# above, but using parallel computing.
# You should use functions parLapply() or mclapply(), 
# and rutils::do_call(), and rbind().

library(parallel)  # load package parallel
num_cores <- detectCores() - 1  # number of cores
clus_ter <- makeCluster(num_cores)  # initialize compute cluster

### write your code here


boot_strap <- mclapply(rep(l_gd, num_boot),
                        FUN=calc_var, num_simu=num_simu,
                        default_thresh=default_thresh,
                        rho_sqrt=rho_sqrt,
                        rho_sqrtm=rho_sqrtm,
                        conf_level = conf_level)  # end mclapply
boot_strap <- rutils::do_call(rbind, boot_strap)


############## Part III
# Summary: Perform Monte Carlo simulation of Brownian motion, 
# to estimate the price of a barrier option.
# (Actual asset prices follow Geometric Brownian motion, but
# in this exercise we will simplify and use Brownian motion, 
# and accept that prices can become negative.)
# Perform bootstrap simulation to estimate the standard error 
# of the estimated price.
# Apply antithetic sampling to the barrier option simulation 
# to reduce the standard error.

# About barrier options:
# One type of barrier option is a knock-out (down and out) 
# call option.
# The knock-out call option pays out the positive difference 
# between the final price minus the strike, but only if the 
# intermediate price never falls bellow the knock-out price.

# Define the simulation parameters:
# number of steps in each simulation
len_gth <- 1000
# strike price
strik_e <- 10
# knock-out barrier level
bar_rier <- (-5)
# number of Brownian motion simulations
n_simu <- 1e4


# 1. (20pts) Simulate multiple Brownian paths, by performing
# an sapply() loop starting from 1 to n_simu.
# Inside the loop perform a vectorized simulation of Brownian 
# motion, and calculate the option pay-out, equal to the 
# positive difference between the final price minus the strike, 
# but only if the intermediate price never falls bellow the 
# knock-out price.
# The sapply() loop should return a numeric vector of length
# n_simu, called pay_outs.
# Hint: you can use an anonymous function that accepts an
# integer argument (the loop count) and returns a numeric
# value.
# You can compare the simulated path values to the bar_rier 
# level, to determine if at any point the path reached below 
# the bar_rier.
# The comparison of the path with bar_rier produces a Boolean
# vector, whose sum is zero only if the path never crossed 
# the bar_rier, and is greater than zero if it did.
# You can use functions sapply(), sum(), cumsum(), and
# rnorm().

# reset random number generator
set.seed(1121)

### write your code here
pay_outs <- sapply(seq(n_simu <- 1e4),function(x){
  #pa_th <- numeric(len_gth)
  pa_th <- cumsum(rnorm(len_gth))
  cro_ss <- which(pa_th < bar_rier)
  if (NROW(cro_ss)>0) {
    pa_th[(cro_ss[1]+1):len_gth] <-
      pa_th[cro_ss[1]]
  }  
  pay_outs <- ifelse(pa_th[len_gth]>strik_e,pa_th[len_gth]-strik_e,0)
  #pay_outs <- round(pay_outs,2)
})

mean(pay_outs)

# Calculate the price of the barrier option as the 
# average value of pay_outs.
# You should get the following result:
# > mean(pay_outs)
# [1] 3.423402


# 2. (20pts) Perform a bootstrap of the simulation 
# from p.1 above, by repeating it in a parallel loop. 
# Inside the bootstrap loop calculate the price of the 
# barrier option as the average value of the pay_outs.
# The bootstrap loop should return a numeric vector of 
# length n_boot, called boot_strap.
# You can use functions detectCores(), makeCluster().
# clusterSetRNGStream(), clusterExport(),  
# and either parSapply() or mcmapply().

# Load package parallel
library(parallel)
# Get the number of cores to use in parallel under Windows
num_cores <- detectCores() - 1
# Initialize compute cluster for parallel bootstrap under Windows
clus_ter <- makeCluster(num_cores)
# number of bootstrap simulations
n_boot <- 10

# Export variables to CPU cores

### write your code here
clusterExport(clus_ter,
              varlist=c("len_gth", "strik_e","bar_rier","n_simu"))

# Reset random number generator for all cores
clusterSetRNGStream(clus_ter, 1121)

### write your code here
boot_strap <- mclapply(seq(n_boot),
                       function(x,...){mean(...)},
                       pay_outs[sample.int(n_simu, replace=TRUE)]
                       , mc.cores=num_cores)  # end mclapply
boot_strap <- rutils::do_call(rbind, boot_strap)
mean(boot_strap)
sd(boot_strap)

stopCluster(clus_ter)

# Calculate the price of the barrier option as the 
# average value of boot_strap.
# Calculate the standard error as the standard deviation 
# of boot_strap.
# You should get results similar to the following:
# > mean(boot_strap)
# [1] 3.529987
# > sd(boot_strap)
# [1] 0.1504767

# Calculate the standard error of the mean of the 
# bootstrap average (that is of mean(boot_strap)):

### write your code here
boot_strap2 <- sapply(seq(1,100,1),function(y)
{mclapply(seq(n_boot),
          function(x,...){mean(...)},
          sample <- 
            pay_outs[sample.int(n_simu, replace=TRUE)]
          , mc.cores=num_cores)})  # end mclapply
boot_strap2 <- rutils::do_call(cbind, boot_strap2)

me_an <- apply(boot_strap2, MARGIN = 2, function(x){mean(x)})

sd(me_an)
# bootstrap average

# 3. (20pts) Apply antithetic sampling to the simulation 
# from p.1 above, and perform a bootstrap simulation for 
# it, as in p.2 above. 

# Reset random number generator for all cores
clusterSetRNGStream(clus_ter, 1121)

### write your code here

### write your code here
boot_strap <- mclapply(seq(n_boot),
                       function(x,...){mean(...)},
                       pay_outs[n_simu-(sample.int(n_simu, replace=TRUE))]
                       , mc.cores=num_cores)  # end mclapply
boot_strap <- rutils::do_call(rbind, boot_strap)
mean(boot_strap)
sd(boot_strap)

stopCluster(clus_ter)


# Calculate the price of the barrier option as the 
# average value of boot_strap.
# Calculate the standard error as the standard deviation 
# of boot_strap.
# You should get results similar to the following:
# > mean(boot_strap)
# [1] 3.540663
# > sd(boot_strap)
# [1] 0.08667923


# 4. (10pts) Plot a histogram of boot_strap.
# Add a curve with the Normal density, with the same mean and sd.
# Add a legend.
# You can use functions hist(), lines(), density(), curve(), and
# legend().

x11(width=6, height=5)

### write your code here
hist(boot_strap, ylim = c(0,10),xlim = c(2.9,3.7))
lines(density(boot_strap),add=TRUE, col="red", lwd=2,add=TRUE)
curve(dnorm(x, mean(boot_strap), sd(boot_strap)),add=TRUE, col="blue", lwd=2)


legend(x="topright",
       legend=c(paste0("mean price=",mean(boot_strap)),paste0("sd=",sd(boot_strap)),"boot_strap histogram", "density","Normal density"),
       title=NULL, inset=0.05, cex=0.8, bg="white",
       bty="n", lwd=5, lty=c(1, 1), col=c("black", "black","black","blue","red"))
# Your plot should be similar to boot_barrier.png


