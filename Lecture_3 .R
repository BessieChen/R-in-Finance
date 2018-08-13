set.seed(1121)  # reset random number generator
# sample from Standard Normal Distribution
len_gth <- 1000
sam_ple <- rnorm(len_gth)
# sample mean - MC estimate
mean(sam_ple)
# sample standard deviation - MC estimate
sd(sam_ple)
# MC estimate of cumulative probability
sam_ple <- sort(sam_ple)
pnorm(1)
sum(sam_ple<1)/len_gth
# MC estimate of quantile
qnorm(0.75)
install(XQuartz)
sam_ple[0.75*len_gth]
x11(width=6, height=5)
par(oma=c(1, 1, 1, 1), mar=c(2, 2, 2, 1), mgp=c(2, 1, 0), cex.lab=0.8, cex.axis=0.8, cex.main=0.8, cex.sub=0.5)
set.seed(1121)  # reset random number generator
lev_el <- 20  # barrier level
len_gth <- 1000  # number of simulation steps
pa_th <- numeric(len_gth)  # allocate path vector
pa_th[1] <- 0  # initialize path
in_dex <- 2  # initialize simulation index
while ((in_dex <= len_gth) &&
 (pa_th[in_dex - 1] < lev_el)) {
# simulate next step
  pa_th[in_dex] <-
    pa_th[in_dex - 1] + rnorm(1)
  in_dex <- in_dex + 1  # advance in_dex
}  # end while
# fill remaining pa_th after it crosses lev_el
if (in_dex <= len_gth)
  pa_th[in_dex:len_gth] <- pa_th[in_dex - 1]
# create daily time series starting 2011
ts_path <- ts(data=pa_th, frequency=20, start=c(2011, 1))
plot(ts_path, type="l", col="black",
     lty="solid", lwd=2, xlab="", ylab="")
abline(h=lev_el, lwd=2, col="red")
title(main="Brownian motion crossing a barrier level",
      line=0.5)
x11(width=6, height=5)
par(oma=c(1, 1, 1, 1), mar=c(2, 2, 2, 1), mgp=c(2, 1, 0), cex.lab=0.8, cex.axis=0.8, cex.main=0.8, cex.sub=0.5)
set.seed(1121)  # reset random number generator
lev_el <- 20  # barrier level
len_gth <- 1000  # number of simulation steps
# simulate path of Brownian motion
pa_th <- cumsum(rnorm(len_gth))
# find index when pa_th crosses lev_el
cro_ss <- which(pa_th > lev_el)
# fill remaining pa_th after it crosses lev_el
if (NROW(cro_ss)>0) {
  pa_th[(cro_ss[1]+1):len_gth] <-
    pa_th[cro_ss[1]]
}  # end if
# create daily time series starting 2011
ts_path <- ts(data=pa_th, frequency=365,
     start=c(2011, 1))
# create plot with horizontal line
plot(ts_path, type="l", col="black",
     lty="solid", lwd=2, xlab="", ylab="")
abline(h=lev_el, lwd=2, col="red")
title(main="Brownian motion crossing a barrier level",
      line=0.5)
set.seed(1121)  # reset random number generator
# sample from Standard Normal Distribution
len_gth <- 1000
sam_ple <- rnorm(len_gth)
# sample mean
mean(sam_ple)
# sample standard deviation
sd(sam_ple)
# bootstrap of sample mean and median
boot_strap <- sapply(1:10000, function(x) {
  boot_sample <- sam_ple[sample.int(len_gth,
                              replace=TRUE)]
  c(mean=mean(boot_sample),
    median=median(boot_sample))
})  # end sapply
boot_strap[, 1:3]
# standard error from formula
sd(sam_ple)/sqrt(len_gth)
# standard error of mean from bootstrap
sd(boot_strap["mean", ])
# standard error of median from bootstrap
sd(boot_strap["median", ])
library(parallel)  # load package parallel
num_cores <- detectCores() - 1  # number of cores
clus_ter <- makeCluster(num_cores)  # initialize compute cluster under Windows
set.seed(1121)  # reset random number generator
# sample from Standard Normal Distribution
len_gth <- 1000
sam_ple <- rnorm(len_gth)
# bootstrap mean and median under Windows
boot_strap <- parLapply(clus_ter, 1:10000,
  function(x, sam_ple, len_gth) {
  boot_sample <- sam_ple[sample.int(len_gth, replace=TRUE)]
  c(mean=mean(boot_sample), median=median(boot_sample))
  }, sam_ple=sam_ple, len_gth=len_gth)  # end parLapply
# bootstrap mean and median under Mac-OSX or Linux
boot_strap <- mclapply(1:10000,
  function(x) {
  boot_sample <- sam_ple[sample.int(len_gth, replace=TRUE)]
  c(mean=mean(boot_sample), median=median(boot_sample))
  }, mc.cores=num_cores)  # end mclapply
boot_strap <- rutils::do_call(rbind, boot_strap)
# means and standard errors from bootstrap
apply(boot_strap, MARGIN=2,
function(x) c(mean=mean(x), sd=sd(x)))
# standard error from formula
sd(sam_ple)/sqrt(len_gth)
stopCluster(clus_ter)  # stop R processes over cluster under Windows
set.seed(1121)  # reset random number generator
# sample from Standard Normal Distribution
len_gth <- 1000
sam_ple <- rnorm(len_gth)
# estimate the 95% quantile
boot_strap <- sapply(1:10000, function(x) {
  boot_sample <- sam_ple[sample.int(len_gth,
                              replace=TRUE)]
  quantile(boot_sample, 0.95)
})  # end sapply
sd(boot_strap)
# estimate the 95% quantile using antithetic sampling
boot_strap <- sapply(1:10000, function(x) {
  boot_sample <- sam_ple[sample.int(len_gth,
                              replace=TRUE)]
  quantile(c(boot_sample, -boot_sample), 0.95)
})  # end sapply
# standard error of mean from bootstrap
sd(boot_strap)
sqrt(2)*sd(boot_strap)
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
de_faults <-
  colSums(t(t(asset_values) < default_thresh))
de_faults / num_simu
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
num_assets <- 3
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
clusterExport(clus_ter,
              varlist=c("default_probs", "rh_o","num_simu","conf_levels"))
boot_strap2 <- parLapply(clus_ter, rep(l_gd, num_boot),
  fun=calc_var, default_thresh=qnorm(default_probs),
  rho_sqrt=rho_sqrt,
  rho_sqrtm=rho_sqrtm,
  num_simu=num_simu,
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
