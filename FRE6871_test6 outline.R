#################################
### FRE6871 Test #6 Outline March 24, 2018
#################################
# Max score 430pts

# This file is an outline for test #6, designed to help prepare for it. 
# The actual test #6 file will be more detailed.

############## Part I
# Summary: Re-order the rows of a data frame according
# to the order of one of its columns.

# 1. (20pts) Create a named vector of student scores 
# between 1 and 10, and call it stu_dents.
# You can use functions sample((), round(), runif(), 
# names(), and c().

set.seed(1121)

### write your code here
n <- 5 #number of score
stu_dents <- round(sample(runif(10)*10,n,replace = FALSE))
stu_dents


# 2. (20pts) Create a data frame of student scores and 
# their score ranks called ros_ter.
# Assign rank strings according to the student score, 
# with the string "first" assigned to the highest student 
# score, and the string "fifth" assigned to the lowest 
# student score.

ra_nks <- c("first", "second", "third", "fourth", "fifth")

# You should permute the ra_nks vector so it corresponds 
# to the student names and scores.
# You should not permute the student names and scores 
# (keep them in their original order).
# You can use functions data.frame() and order().

### write your code here
# L3 <- LETTERS[1:3]
# fac <- sample(L3, 10, replace = TRUE)
# (d <- data.frame(x = 1, y = 1:10, fac = fac))
# ## The "same" with automatic column names:
# data.frame(1, 1:10, sample(L3, 10, replace = TRUE))

rank <- ra_nks[rev(order(stu_dents))]
data.frame(score = stu_dents, rank = rank )


############## Part II
# Summary: Download data from FRED for rates and 
# currencies, and perform PCA analysis.

# Load rutils package.
library(rutils)

# Define symbols for Treasury rates and exchange rates:
sym_bols <- c("DGS1", "DGS5", "DGS20", "DSG100", "DEXJPUS", "DEXUSEU", "DEXUSUK")

# 1. (20pts)
# Create a new environment for data called env_rates.
# Use new.env()

### write your code here
env_rates <- new.env()

# Download data from FRED for sym_bols into env_rates.
# Use quantmod::getSymbols()

### write your code here
library(quantmod) 
quantmod::getSymbols(sym_bols, env=env_rates, src="FRED")


# You will get an error message, because one of the 
# symbols is wrong.
# List the names of all the xts series in env_rates.
# Use ls()

### write your code here
ls(env_rates)

# You should get the following result:
# [1] "DGS1"  "DGS20" "DGS5" 

# One of the symbols in sym_bols is wrong, which caused 
# an error condition in getSymbols(), which aborted the 
# loading of data for the remaining symbols in sym_bols.
# 
# Calculate a Boolean vector indicating the symbols 
# that have been downloaded already.
# Use function ls() and the %in% operator.

### write your code here
sym_bols[sym_bols %in% ls(env_rates)]

# You have to write a loop to download the data one 
# symbol at a time using getSymbols(), for the 
# remaining symbols in sym_bols that were not downloaded.
# Your loop should download data only for those symbols
# that were not downloaded yet.
# You can either use an sapply() or a for() loop.
# You have to call getSymbols() wrapped inside 
# tryCatch(), with error and warning handlers, and 
# a finally function.
# You can also use functions print() and paste().

### write your code here
#tryCatch(1, finally = print("Hello"))
sapply(seq_along(sym_bols), function(index) {  # anonymous function
  tryCatch(  # with error handler
    {  # body
      stopifnot(!sym_bols[index] %in% ls(env_rates))  # check for error
      # broadcast message to console
      #cat("(cat) num_var =", num_var, "\t")
      # # return a value
      quantmod::getSymbols(sym_bols[index], env=env_rates, src="FRED")
      #cat("(cat) num_var =", index, "\t")
      #paste("(return) num_var =", index)
      
    },
    error=function(error_cond) {
      print(paste("error handler: ", error_cond))
    },
    warning=function(warning_cond) {
      print(paste("warning handler: ", warning_cond))
    },
    finally=#quantmod::getSymbols(sym_bols, env=env_rates, src="FRED")
      print(paste("already downloaded", sym_bols[index]))
  )  # end tryCatch
}  # end anonymous function
)

# List the names of all the xts series in env_rates.
# Use ls()

### write your code here
ls(env_rates)

# Combine (cbind) all the xts time series in env_rates
# into a single xts time series called rate_s.
# Calculate the class of re_turns.
# You can use functions rutils::do_call(), cbind(),
# class(), and as.list(env_rates).

### write your code here
rate_s <- rutils::do_call(cbind,as.list(env_rates))
class(rate_s)

# Calculate the number of NA values in rate_s.

### write your code here
sum(sapply(seq_along(ls(env_rates)),function(x)
{sum(is.na(as.list(rate_s)[[x]]))}))

# Carry forward non-NA values in rate_s, and then 
# remove any rows with remaining NA values.
# You can use functions zoo::na.locf() and na.omit().

### write your code here
rate_s <- na.omit(zoo::na.locf(rate_s))


# Calculate the daily changes in rates (differences), 
# and call it re_turns.
# You can use function rutils::diff_it().

### write your code here
re_turns <- rutils::diff_it(rate_s)

# Center and scale the columns of re_turns so that their 
# means are zero and their standard deviations are equal to 1.
# You can use functions apply() and scale(), or functions 
# t(), colMeans, colSums(), and sqrt().

### write your code here
#method 1: using scale()
date_s <- index(re_turns)
re_turns2 <- apply(re_turns,2,scale)
re_turns2 <- xts(re_turns2, date_s)

#method 2:using t(),...
re_turns <- t(t(re_turns) - colMeans(re_turns))
re_turns <- t(t(re_turns) / sqrt(colSums(re_turns^2)/(NROW(re_turns)-1)))


# Calculate the means and standard deviations of the 
# columns of re_turns to verify that they are all 
# equal to 0 and 1 (within machine precision).
# You can use functions colMeans, apply(), and sd().

### write your code here
.Machine$double.eps
apply(re_turns,2,function(x){c(means = mean(x), sd = sd(x))})


# 3. (20pts) Perform principal component analysis PCA
# on re_turns, and call it pc_a.
# Use function prcomp().

### write your code here
pc_a <- prcomp(re_turns, scale=TRUE)

# Calculate the time series of principal components from 
# re_turns and pc_a$rotation, and call it pca_rets.
# You can use the "%*%" operator.

### write your code here

pca_rets <- re_turns %*% pc_a$rotation


# Demonstrate that pca_rets is equal to pc_a$x.
# Use function all.equal()

### write your code here
all.equal(pc_a$x,pca_rets)
# #I just want to change it to time series.
# date_s <- index(re_turns)
# pca_rets <- xts(re_turns %*% pc_a$rotation,
#                 order.by=date_s)


# Calculate the correlations between the pca_rets.
# Use function cor().

### write your code here
cov_ret <- cor(pca_rets)
cov_ret

# 4. (20pts) 
# Invert the matrix pc_a$rotation, and call it inv_rotation.
# You can use function solve().

### write your code here
inv_rotation <- solve(pc_a$rotation)
inv_rotation

# Calculate the time series of returns by multiplying the 
# principal components time series pca_rets by inv_rotation,
# and call it sol_ved.
# Use the %*% operator.

### write your code here
sol_ved <- pca_rets %*% inv_rotation
all.equal(re_turns, sol_ved)
# [1] TRUE

# You should get the following result:
# > all.equal(re_turns, sol_ved)
# [1] TRUE

# Calculate the number of principal components with the 
# largest eigenvalues, which sum up to at least 80% of 
# the total variance.
# Remember that pc_a$sdev are the standard deviations.
# You can use functions which(), cumsum(), and sum().

### write your code here
n_comps_largest <- which((cumsum(pc_a$sdev^2)/sum(pc_a$sdev^2)>0.8))[1]
n_comps_largest

# Calculate the time series of returns from the time 
# series of principal components with the largest 
# eigenvalues which sum up to at least 80% of the total 
# variance, and call it sol_ved.
# hint: First select columns from pca_rets and 
# inv_rotation, and then multiply them together.

### write your code here
date_s <- index(rate_s)
weight_s <- pc_a$rotation[, 1:n_comps_largest]
sol_ved <- pca_rets %*% weight_s
sol_ved <- xts::xts(sol_ved, date_s)
sol_ved

# 5. (20pts) Calculate the cumulative returns of 
# re_turns and sol_ved.
# Use function xts:::cumsum.xts().

### write your code here

# #choice 1:
# weight_s <- pc_a$rotation[,1:n_comps_largest]  #pc_a$rotation[,3]  means under the third component, the weight of each 6  stock. 
# sol_ved <- re_turns %*% weight_s
# dim(sol_ved)
# # > dim(sol_ved)
# # [1] 5014    3
# 
# #choice 2:
# sol_ved <- pca_rets[,1:n_comps_largest] %*% inv_rotation[1:n_comps_largest, ]
# dim(sol_ved)
# # > dim(sol_ved)
# # [1] 5014    6
# all.equal(sol_ved,re_turns)

#choice 3:
sol_ved <- pca_rets[,1:n_comps_largest] %*% inv_rotation[1:n_comps_largest, ]
#dim(sol_ved)
#rownames(sol_ved) <- (re_turns)
# head(sol_ved)
# sol_ved

cum_sol_ved <- xts:::cumsum.xts(sol_ved)
cum_returns <- xts:::cumsum.xts(re_turns)
# head(cum_sol_ved)
# head(cum_returns)
# cum_returns[, 1]
# cum_sol_ved[, 1]
# Plot the cumulative returns of re_turns and sol_ved
# in 6 panels, with 3 rows and 2 columns.
# Each panel should have a plot with two lines.
# Use functions x11(), par(), for(), zoo::plot.zoo(), 
# cbind(), legend(), and paste0().

x11(width=6, height=7)

### write your code here
list <- ls(env_rates)
par(mfrow=c(3, 2))
par(mar=c(2, 2, 0, 1), oma=c(0, 0, 0, 0))
for (sym_bol in list) {
  plot.zoo(
    cbind(cum_returns[, sym_bol], cum_sol_ved[, sym_bol]),
    plot.type="single", col=c("black", "blue"), xlab="", ylab="")
  legend(x="topleft", bty="n",
         legend=paste0(sym_bol, c("", " solved")),
         title=NULL, inset=0.05, cex=1.0, lwd=6,
         lty=1, col=c("black", "blue"))
}  # end for




############## Part III
# Summary: Perform aggregations on a data frame of panel 
# data using the split-apply-combine procedure. 

# 1. (30pts) Download the file CRSPpanel.txt from NYU 
# Classes. 
# The file CRSPpanel.txt contains a data frame with 
# a single day of panel data. The panel data contains 
# fundamental financial data for 265 S&P500 stocks. 
# Read the file into a data frame called panel_data using 
# read.table(), with the header and sep arguments. 

panel_data <- read.table(file="/Users/bessie.chan/Desktop/nyu/courses/6871 R/resourses/CRSPpanel.txt", 
                         header=TRUE, sep="\t")


# The Industry column has 22 unique elements, while
# the Sector column has 10 unique elements. 
# 
# Each unique Industry belongs to a single Sector. 
# But each Sector may have several Industries that 
# belong to it.

# Calculate a data frame with the unique Industries 
# in the first column, and the Sector to which they 
# belong in the second column, and call it ind_sec. 
# You can use the functions match(), levels(), with() 
# and data.frame().
# You cannot use tapply(), aggregate(), or any loops. 

### write your code here
ls(panel_data)
ind <- levels(panel_data$Industry)
sec <- levels(panel_data$Sector)
aa <- seq_along(ind)
bb <- seq_along(sec)
ind_index <- match(ind[aa],panel_data$Industry) #get the first index of Ind level from panel_data
belonged_sector <- panel_data$Sector[ind_index]
ind_sec <- data.frame(Industries=ind,Sectors=belonged_sector)
ind_sec

#another method using with()
get_levels <- function(x){levels(factor(x))}
belonged_sector <- with(panel_data,
                        Sector[
                          match(get_levels(Industry)[seq_along(get_levels(Industry))],Industry)
                          ])
ind_sec <- data.frame(Industries=ind,Sectors=belonged_sector)
ind_sec

#another method from kw:
table <- with(panel_data, data.frame(Industry, Sector))
indus_tries <- levels(panel_data$Industry)
ind_sec <- data.frame(table[match(indus_tries, table$Industry),1],table[match(indus_tries, table$Industry),2])
rownames(ind_sec) <- ind_sec[,1]
colnames(ind_sec) <- c("Industries", "Sectors")

# 2. (20pts) Calculate the same data frame as ind_sec,
# but using function sapply().  Call it ind_sec_2.
# 
# The sapply() loop will return a named vector of strings.
# Call it ind_sec_2.
# The string elements of ind_sec_2 should be the Sectors 
# and the vector names should be the Industries that 
# belong to the Sectors.
# So ind_sec_2 will have many elements with the same 
# Sector, but the vector names should be all different.
# 
# hint: You can either perform an sapply() loop over the 
# levels of panel_data$Industry, or over unique elements 
# of indus_tries. 
# 
# You can also use the functions levels(), match(), 
# unique(), and an anonymous function.

### write your code here
ind <- levels(panel_data$Industry)
sec <- levels(panel_data$Sector)
ind_sec_2 <- sapply(seq_along(ind),function(x){
  aa <- with(panel_data,Sector[match(ind[x],Industry)])
  names(aa) <- ind[x]
  return(aa)
  }      )
ind_sec_2

# Verify that the names(ind_sec_2) match exactly ind_sec[, 1],
# and if they don't then appply a permutation (re-order) to 
# ind_sec_2 so that they do. 
# You can use the functions as.vector(), match(), 
# names(), and all.equal().

### write your code here
a <- as.vector(names(ind_sec_2))
b <- as.vector(ind_sec[, 1])
all.equal(a,b)


# Convert the named vector ind_sec_2 into a data frame equal 
# to as ind_sec.
# Verify that ind_sec_2 is equal to ind_sec using all.equal().
# You can use the functions data.frame(), rownames(), 
# names(), and all.equal().

### write your code here
ind_sec_2 <- data.frame(Industries = names(ind_sec_2),Sectors = as.vector(ind_sec_2))
all.equal(ind_sec_2,ind_sec)

# 3. (20pts) Calculate the same data frame as ind_sec,
# but using function aggregate().  Call it ind_sec_3.
# Rename the columns of ind_sec_3 to: "Industries" and "Sectors".
# Verify that ind_sec_3 is equal to ind_sec using all.equal().
# You can also use the functions c(), list(), colnames(), 
# with(), and all.equal().


### write your code here
ind_sec_3 <- aggregate(x=panel_data$Sector,by=list(panel_data$Industry),FUN=unique)
colnames(ind_sec_3) <- c("Industries","Sectors")
ind_sec_3

all.equal(ind_sec_3,ind_sec)

# 4. (30pts) Calculate the same data frame as ind_sec,
# but using function tapply().  Call it ind_sec_4.
# Perform tapply() on the Sector and Industry columns.
# hint: You can convert the Sector column into a vector.

### write your code here
get_sec <- function(x)
{
  y <- unique(x)
  levels(factor(panel_data$Sector))[y]
}
#as.vector(panel_data$Sector)
#ind_sec_4 <- do.call(cbind,tapply(X=panel_data$Sector, INDEX=panel_data$Industry, FUN=get_sec))
ind_sec_4 <- tapply(X=as.vector(panel_data$Sector), INDEX=panel_data$Industry, FUN=unique)


# tapply() returns an array which you must convert into 
# a data frame. 
# You can use the functions data.frame(), rownames(), 
# names(), and as.vector().

### write your code here
ind_sec_4 <- data.frame(Industries = names(ind_sec_4),Sectors = as.vector(ind_sec_4))



# Verify that ind_sec_2 is equal to ind_sec using all.equal().

### write your code here
all.equal(ind_sec_4,ind_sec)


# 5. (20pts) Each Sector has one or more Industries 
# that belong to it. 
# Calculate a named list of vectors (not factors) of 
# strings with the Industries belonging to the Sectors, 
# and call it sec_ind.
# There are at least three ways of doing this.

# First you should use the function split(). 
# You can also use functions as.vector(), sapply(), 
# unique(), and with(). 

### write your code here
sec_ind <- with(panel_data,sapply(split(Industry, Sector),function(x){as.vector(unique(x))}))
sec_ind

# In the second method you should perform an sapply() 
# loop over levels(Sector), and call the output sec_ind_2.
# You can also use functions as.vector(), unique(), 
# levels(), with(), and an anonymous function. 

### write your code here
sec_ind_2 <- with(panel_data, sapply(levels(Sector),function(x)
  {as.vector(unique(Industry[Sector==x]))}))
sec_ind_2

# Verify that sec_ind_2 is equal to sec_ind using all.equal().

### write your code here
all.equal(sec_ind, sec_ind_2)


# In the third method you should use function tapply(). 
# and call the output sec_ind_3.
# You can also use functions as.vector(), unique(), 
# drop(), as.matrix(), and with(). 

### write your code here

indus_tries <- levels(panel_data$Industry)
sec_ind_3 <- drop(as.matrix(with(panel_data, tapply(Industry, Sector, function(x) indus_tries[unique(x)]))))

# Verify that sec_ind_3 is equal to sec_ind using all.equal().

### write your code here
all.equal(sec_ind, sec_ind_3)


# 6. (20pts) Calculate a named list with the stock tickers 
# of companies in each industry, and call it industry_tickers. 
# 
# You must perform the calculation in two different ways. 
# In the first method you must use function split(), and 
# you can also use functions as.vector() and with().
# hint: Use the columns TICKER and Industry.

### write your code here
#names(panel_data)
industry_tickers <- with(panel_data,as.vector(split(TICKER, Industry)))
#class(industry_tickers)

# In the second method you must use function sapply(). 
# You cannot use tapply(). 
# You can also use the functions levels(), as.vector(), 
# with(), and an anonymous function. 

### write your code here
industry_tickers2 <- with(panel_data,sapply(levels(Industry),function(x)
  {
  TICKER[Industry==x]
})
     )


# Verify that both methods produce the same result.
# Use function all.equal().

### write your code here
all.equal(industry_tickers,industry_tickers2)

# 7. (20pts) Calculate a named vector with the number of 
# companies in each Industry. 
# You can use functions sapply() and NROW(). 
# hint: you can use the vector industry_tickers.

### write your code here
#method 1:
number_of_ind <- sapply(industry_tickers,NROW)
# #method 2:
# number_of_ind2 <- with(panel_data,sapply(levels(Industry),function(x)
# {
#   sum(Industry==x)
# })
# )



# Calculate a named list with the indices of companies in 
# each industry, and call it industry_indices. 
# The index of a company is its row number in panel_data.
# You can use functions sapply() and match(). 
# hint: you can use the vector industry_tickers.

### write your code here
industry_indices <- with(panel_data,sapply(industry_tickers, function(x){match(x,TICKER)}))
industry_indices
class(industry_indices) #"list"



# 8. (30pts) Calculate a named vector (not an array!) with 
# the average "NET.INCOME" of all the companies in each Sector. 
# You can use functions with(), split(), sapply(), and mean(). 

### write your code here
#with(panel_data,split(NET.INCOME,Sector))
average_net_inc <- with(panel_data,sapply(split(NET.INCOME,Sector),mean))
average_net_inc
class(average_net_inc)
#as.vector(average_net_inc)



# 9. (20pts) Calculate a vector of tickers of the 
# companies which have the highest ROE in each 
# Industry, and call it ticker_s.
# You must use functions sapply(), split(), 
# as.vector(), which.max(), and with().

### write your code here
#with(panel_data,split(ROE,Industry))
ticker_s <- with(panel_data,sapply(split(ROE,Industry),function(x){TICKER[match(x[which.max(x)],ROE)]}))
ticker_s
class(ticker_s)

# Subset the panel_data data frame and extract 
# the rows corresponding to ticker_s, and call it 
# max_roes.
# You can use function match(). 

### write your code her
max_roes <- with(panel_data,panel_data[sort(match(ticker_s,TICKER)),])
max_roes


############## Part IV
# Summary: List the class and dimension attributes 
# of objects in the environment rutils::env_etf.
# The environment rutils::env_etf contains ETF time 
# series and other ETF data.

# 1. (10pts) List the names of all the objects in 
# the environment rutils::env_etf, and save the names 
# in a vector of strings called name_s.
# You can use function ls().

### write your code here
name_s <- as.vector(ls(rutils::env_etf))
name_s
class(name_s)

# Create a list with the class attributes of all 
# the objects in the environment rutils::env_etf, 
# and call it class_es.
# You can use functions eapply() and class().

### write your code here
class_es <- eapply(rutils::env_etf,class)
class_es


# 2. (20pts) List the names of all the objects of
# class xts in the environment rutils::env_etf, and 
# save the names in a vector of strings called name_s.
# hint: first calculate a named Boolean vector that
# is TRUE for objects of class xts, then extract the
# names of those objects.
# You can use functions eapply(), is.xts(), unlist(),
# and names().
# Or you can use sapply() instead of eapply(), and
# an anonymous function, and the "%in%" operator.

### write your code here
# x[is.xts(x)]
# aa <- eapply(rutils::env_etf,function(x){strsplit([1], split='.')[1]})
# bb <- unlist(aa)
# ls(rutils::env_etf)[aa]
sapply(class_es,function(x){if(is.xts(class_es[[3]])){names(x)}})
bool <- unname(unlist(eapply(rutils::env_etf,is.xts)))
names(class_es)[bool]


# 3. (20pts) Create a matrix called dimen_sions,
# containing the dimensions of all the xts objects
# in rutils::env_etf.
# You can use the functions eapply(), dim(), rbind(),
# and do.call().
# Or you can use lapply() and an anonymous function,
# instead of eapply().
# hint: eapply() returns a list, and you must flatten
# the list into a matrix using rbind() and do.call().

### write your code here
dimen_sions <- do.call(rbind,eapply(rutils::env_etf,dim)[bool])
dimen_sions
