#################################
### FRE6871 Test #6 March 24, 2018
#################################
# Max score 430pts

# Please write in this file the R code needed to perform
# the tasks below, rename the file to your_name_test6.R
# and upload it to NYU Classes,

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note the extra assignments marked with *extra*
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

############## Part I
# Summary: Re-order the rows of a data frame according
# to the order of one of its columns.

# 1. (20pts) Create a named vector of student scores 
# between 1 and 10, and call it stu_dents.
# You can use functions sample(), round(), runif(), 
# names(), and c().

set.seed(1121)

### write your code here
n <- 5 #number of score
stu_dents <- round(sample(runif(10)*10,n,replace = FALSE),2)
names(stu_dents) <- c("Angie", "Chris", "Suzie",  "Matt",  "Liz")
stu_dents

# You should get the following output:
# > stu_dents
# Angie Chris Suzie  Matt   Liz 
#  6.05  3.87  3.49  7.02  6.02


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
rank <- ra_nks[rev(order(stu_dents))]
ros_ter <- data.frame(score = stu_dents, rank = rank )
ros_ter

# You should get the following output:
# Notice that the student names and scores aren't permuted 
# (they are in the original order).
# > ros_ter
#          score   rank
# Angie     6.05 second
# Chris     3.87 fourth
# Suzie     3.49  fifth
# Matt      7.02  first
# Liz       6.02  third



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

# You should get the following result:
# [1] "DEXJPUS" "DEXUSEU" "DEXUSUK" "DGS1"    "DGS20"   "DGS5"


# extra

# 2. (20pts)
# Calculate a vector with the number of NA values 
# in all the xts series in env_rates.
# You can use functions unlist(), eapply(), sum(), 
# and is.na().

### write your code here
number_na <- function(x)
{sum(is.na(x))}
unlist(eapply(env_rates,number_na))

# You should get output similar to the following 
# (depending on when you downloaded the data):
# DEXUSUK DEXJPUS    DGS1   DGS20 DEXUSEU    DGS5 
#     469     475     627     262     185     627

# Calculate a vector of strings representing the 
# start dates of all the xts series in env_rates.
# You can use functions unlist(), eapply(), format(), 
# and start().

### write your code here
get_start <- function(x)
{format(index(x[1]))}

get_end <- function(x)
{format(index(x[NROW(x)]))}
unlist(eapply(env_rates,get_end))


# You should get output similar to the following 
# (depending on when you downloaded the data):
#      DEXUSUK      DEXJPUS         DGS1        DGS20      DEXUSEU         DGS5 
# "1971-01-04" "1971-01-04" "1962-01-02" "1993-10-01" "1999-01-04" "1962-01-02"

# extra end

# Combine (cbind) all the xts time series in env_rates
# into a single xts time series called rate_s.
# Calculate the class of re_turns.
# You can use functions rutils::do_call(), cbind(),
# class(), and as.list().

### write your code here
rate_s <- rutils::do_call(cbind,as.list(env_rates)[sym_bols])
class(rate_s)
ls(rate_s)

# Calculate the number of NA values in rate_s.

### write your code here
sum(sapply(seq_along(ls(rate_s)),function(x)
{sum(is.na(as.list(rate_s)[[x]]))}))

# You should get output similar to the following 
# (depending on when you downloaded the data):
# [1] 25292

# Carry forward non-NA values in rate_s, and then 
# remove any rows with remaining NA values.
# You can use functions zoo::na.locf() and na.omit().

### write your code here
rate_s <- na.omit(zoo::na.locf(rate_s))
dim(rate_s)
round(head(rate_s), 2)

# You should get output similar to the following 
# (depending on when you downloaded the data):
# > dim(rate_s)
# [1] 5004    6
# 
# > round(head(rate_s), 2)
#            DGS1 DGS5 DGS20 DEXJPUS DEXUSEU DEXUSUK
# 1999-01-04 4.58 4.57  5.42  112.15    1.18    1.66
# 1999-01-05 4.56 4.62  5.48  111.15    1.18    1.66
# 1999-01-06 4.53 4.61  5.42  112.78    1.16    1.65
# 1999-01-07 4.51 4.62  5.48  111.69    1.17    1.65
# 1999-01-08 4.57 4.72  5.57  111.52    1.16    1.64
# 1999-01-11 4.62 4.76  5.61  108.83    1.15    1.64

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
re_turns <- t(t(re_turns) - colMeans(re_turns))
re_turns <- t(t(re_turns) / sqrt(colSums(re_turns^2)/(NROW(re_turns)-1)))
round(head(re_turns), 2)

# You should get the following result:
# > round(head(re_turns), 2)
#       DGS1  DGS5 DGS20 DEXJPUS DEXUSEU DEXUSUK
# [1,]  0.01  0.01  0.01    0.00    0.00    0.01
# [2,] -0.50  0.84  1.12   -1.45   -0.70   -0.16
# [3,] -0.75 -0.16 -1.10    2.37   -1.66   -0.20
# [4,] -0.50  0.17  1.12   -1.58    0.48   -0.55
# [5,]  1.55  1.68  1.67   -0.25   -1.58   -0.96
# [6,]  1.29  0.67  0.75   -3.91   -0.27   -0.32


# Calculate the means and standard deviations of the 
# columns of re_turns to verify that they are all 
# equal to 0 and 1 (within machine precision).
# You can use functions colMeans, apply(), and sd().

### write your code here
.Machine$double.eps
apply(re_turns,2,function(x){c(means = mean(x), sd = sd(x))})


# extra

# Copy the xts series rate_s and re_turns into env_rates.
# Do it in two different ways.  First use the function assign().

### write your code here
class(env_rates)
assign("rate_s", rate_s, envir=env_rates)
assign("re_turns", re_turns, envir=env_rates)

# Copy the xts series using the "$" operator.

### write your code here
env_rates$rate_s <- rate_s
env_rates$re_turns <- re_turns
ls(env_rates)


# Save env_rates to the file rates_data2.RData
# You can use function save().

### write your code here
save(env_rates, file='rates_data2.RData')

# extra end


# 3. (20pts) Perform principal component analysis PCA
# on re_turns, and call it pc_a.
# Use function prcomp().

### write your code here
pc_a <- prcomp(re_turns, scale=TRUE)

# Calculate the time series or matrix of principal components 
# from re_turns and pc_a$rotation, and call it pca_rets.
# You can use the "%*%" operator.

### write your code here
pc_a$sdev
pc_a$x
pc_a$rotation
pca_rets <- re_turns %*% pc_a$rotation


# Demonstrate that pca_rets is equal to pc_a$x.
# Use function all.equal(), and if needed zoo::coredata().

### write your code here
all.equal(pc_a$x,pca_rets)

# Calculate the correlations between the pca_rets.
# Use function cor().

### write your code here
cov_ret <- cor(pca_rets)
cov_ret
round(cov_ret, 6)

# You should get the following result:
# > round(cor_mat, 6)
#     PC1 PC2 PC3 PC4 PC5 PC6
# PC1   1   0   0   0   0   0
# PC2   0   1   0   0   0   0
# PC3   0   0   1   0   0   0
# PC4   0   0   0   1   0   0
# PC5   0   0   0   0   1   0
# PC6   0   0   0   0   0   1


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
sum(pc_a$sdev^2)

# You should get the following result:
# [1] 3
# 
# > sum(pc_a$sdev^2)
# [1] 6

# Calculate the time series of returns from the time 
# series of principal components with the largest 
# eigenvalues which sum up to at least 80% of the total 
# variance, and call it sol_ved.
# hint: First select columns from pca_rets and 
# inv_rotation, and then multiply them together.

### write your code here
sol_ved <- pca_rets[,1:n_comps_largest] %*% inv_rotation[1:n_comps_largest, ]
sol_ved
round(head(sol_ved), 3)

# You should get output similar to this:
# > round(head(sol_ved), 3)
#        DGS1   DGS5  DGS20 DEXJPUS DEXUSEU DEXUSUK
# [1,]  0.008  0.010  0.009   0.002   0.002   0.003
# [2,]  0.599  0.530  0.447  -1.254  -0.180  -0.472
# [3,] -0.675 -0.647 -0.637   2.485  -1.162  -0.691
# [4,]  0.361  0.262  0.223  -1.572   0.203  -0.137
# [5,]  1.637  1.747  1.513  -0.156  -1.197  -1.298
# [6,]  1.108  0.860  0.719  -3.841   0.141  -0.715


# 5. (20pts) Calculate the cumulative returns of 
# re_turns and sol_ved.
# Use function xts:::cumsum.xts().

### write your code here
cum_sol_ved <- xts:::cumsum.xts(sol_ved)
cum_returns <- xts:::cumsum.xts(re_turns)


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
    plot.type="single", col=c("black", "blue"), xlab="", ylab="",
    xy.labels = TRUE)
  legend(x="topleft", bty="n",
         legend=paste0(sym_bol, c("", " solved")),
         title=NULL, inset=0.05, cex=1.0, lwd=6,
         lty=1, col=c("black", "blue"))
}  #

# Your plot should be similar to pca_rates_series_solved.png



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

#panel_data <- read.table(file="C:/Develop/R/lecture_slides/data/CRSPpanel.txt", header=TRUE, sep="\t")
panel_data <- read.table(file="/Users/bessie.chan/Desktop/nyu/courses/6871 R/resourses/CRSPpanel.txt", 
                         header=TRUE, sep="\t")

dim(panel_data)
panel_data[1:3, 1:4]
# You should get the following output for panel_data:
# > dim(panel_data)
# [1] 265  48
# 
# > panel_data[1:3, 1:4]
#       DATE PERMNO    CUSIP     COMPANY.NAME
# 1 20031231  26403 25468710   DISNEY WALT CO
# 2 20031231  89525 20030N10 COMCAST CORP NEW
# 3 20031231  66181 43707610   HOME DEPOT INC

# extra

# Calculate the names of the columns of panel_data 
# that are factors.
# You can use the functions colnames(), sapply(), 
# and is.factor().

### write your code here
colnames(panel_data)[sapply(panel_data,is.factor)]


# You should get the following output:
# [1] "CUSIP"   "COMPANY.NAME"  "TICKER"  "Sector"  "Industry"    

# Calculate the number of numeric columns of panel_data.
# You can use the functions NROW(), colnames(), sapply(), 
# and is.numeric(). 

### write your code here
NROW(colnames(panel_data)[sapply(panel_data,is.numeric)])


# You should get the following output:
# [1] 43

# Note the Industry and Sector columns of panel_data. 
# Calculate the class of the Industry and Sector columns 
# of panel_data.

### write your code here
class(panel_data$Industry)
class(panel_data$Sector)

# The Industry and Sector columns are factors.
# Calculate the number of unique elements of the Industry 
# and Sector columns. 
# Use the functions NROW() and levels(). 

### write your code here
NROW(levels(panel_data$Industry))
NROW(levels(panel_data$Sector))



# extra end

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
head(panel_data)
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

# You should get the following output:
# > ind_sec
#                                  Industries                    Sectors
# 1                  Automobiles & Components     Consumer Discretionary
# 2                             Capital Goods                Industrials
# 3        Commercial & Professional Services                Industrials
# 4               Consumer Durables & Apparel     Consumer Discretionary
# 5                    Diversified Financials                 Financials
# 6                                    Energy                     Energy
# 7                     Food & Drug Retailing           Consumer Staples
# 8                  Food, Beverage & Tobacco           Consumer Staples
# 9          Health Care Equipment & Services                Health Care
# 10             Hotels Restaurants & Leisure     Consumer Discretionary
# 11            Household & Personal Products           Consumer Staples
# 12                                Insurance                 Financials
# 13                                Materials                  Materials
# 14                                    Media     Consumer Discretionary
# 15          Pharmaceuticals & Biotechnology                Health Care
# 16                                Retailing     Consumer Discretionary
# 17 Semiconductors & Semiconductor Equipment     Information Technology
# 18                      Software & Services     Information Technology
# 19          Technology Hardware & Equipment     Information Technology
# 20               Telecommunication Services Telecommunication Services
# 21                           Transportation                Industrials
# 22                                Utilities                  Utilities


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
# hint: You can perform an sapply() loop over the 
# levels of panel_data$Industry. 
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


# You should get the following output, with names of Sectors 
# as list element names, and names of Industries with quotes: 
# 
# > sec_ind
# $`Consumer Discretionary`
# [1] "Media"                        "Retailing"                   
# [3] "Hotels Restaurants & Leisure" "Automobiles & Components"    
# [5] "Consumer Durables & Apparel" 
# 
# $`Consumer Staples`
# [1] "Household & Personal Products" "Food, Beverage & Tobacco"     
# [3] "Food & Drug Retailing"        
# 
# $Energy
# [1] "Energy"
# 
# $Financials
# [1] "Insurance"              "Diversified Financials"
# 
# $`Health Care`
# [1] "Pharmaceuticals & Biotechnology"  "Health Care Equipment & Services"
# 
# $Industrials
# [1] "Commercial & Professional Services" "Capital Goods"                     
# [3] "Transportation"                    
# 
# $`Information Technology`
# [1] "Technology Hardware & Equipment"         
# [2] "Software & Services"                     
# [3] "Semiconductors & Semiconductor Equipment"
# 
# $Materials
# [1] "Materials"
# 
# $`Telecommunication Services`
# [1] "Telecommunication Services"
# 
# $Utilities
# [1] "Utilities"


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
industry_tickers <- with(panel_data,as.vector(split(TICKER, Industry)))
industry_tickers

# You should get the following output, with names 
# of Industries as the list element names:
# 
# > industry_tickers
# $`Automobiles & Components`
# [1] "JCI" "BWA"
# 
# $`Capital Goods`
# [1] "BA"   "MMM"  "HON"  "EMR"  "DHR"  "CMI"  "ETN"  "LMT" 
# [9] "PCP"  "ITW"  "GD"   "RTN"  "NOC"  "GWW"  "PH"   "ROK" 
# [17] "DOV"  "IR"   "FLR"  "TYC"  "FAST" "PNR"  "ROP"  "COL" 
# [25] "PLL"  "LLL"  "MAS"  "JEC" 
# 
# $`Commercial & Professional Services`
# [1] "MNST" "EFX"  "SRCL" "RSG"  "RHI"  "CTAS" "PBI"  "AVY" 
# 
# $`Consumer Durables & Apparel`
# [1] "VFC"  "COH"  "WHR"  "PVH"  "NWL"  "FOSL" "HAS"  "GRMN"
# [9] "LEG"  "DFS"  "SNA"  "MSI" 
# 
# $`Diversified Financials`
# [1] "MCO" "LUK"
# 
# $Energy
# [1] "XOM" "CVX" "SLB" "OXY" "COP" "PXD" "EOG" "APC" "HAL" "APA"
# [11] "NBL" "MRO" "BHI" "VLO" "DVN" "FTI" "RRC" "TSO" "CHK" "MUR"
# [21] "DNR" "ESV" "RDC" "NE"  "NBR" "DO"  "HP"  "NFX"
# 
# $`Food & Drug Retailing`
# [1] "CVS"  "COST" "WAG"  "KR"   "SYY"  "SWY"  "SE"   "WIN" 
# 
# $`Food, Beverage & Tobacco`
# [1] "KO"  "PEP" "GIS" "ADM" "K"   "HSY" "SJM" "STZ" "CCE" "TSN"
# [11] "CPB" "HRL"
# 
# $`Health Care Equipment & Services`
# [1] "ESRX" "MDT"  "BAX"  "WLP"  "AET"  "SYK"  "BDX"  "HUM" 
# [9] "BSX"  "STJ"  "ISRG" "CERN" "ZMH"  "LIFE" "DGX"  "DVA" 
# [17] "LH"   "VAR"  "EW"   "XRAY" "THC" 
# 
# $`Hotels Restaurants & Leisure`
# [1] "MCD"  "YUM"  "WYNN" "DRI"  "IGT" 
# 
# $`Household & Personal Products`
# [1] "PG"  "KMB" "EL"  "CLX" "AVP"
# 
# $Insurance
# [1] "MMC"
# 
# $Materials
# [1] "BTU" "DD"  "DOW" "PX"  "PPG" "APD" "IP"  "NUE" "NEM" "EMN"
# [11] "FMC" "AA"  "ARG" "VMC" "MWV" "IFF" "SEE" "OI"  "BMS" "CLF"
# [21] "ATI" "X"  
# 
# $Media
# [1] "DIS"   "CMCSA" "OMC"   "IPG"   "GCI"   "CVC"   "WPO"  
# 
# $`Pharmaceuticals & Biotechnology`
# [1] "JNJ"  "PFE"  "MRK"  "GILD" "AMGN" "BMY"  "BIIB" "AGN" 
# [9] "REGN" "ALXN" "PRGO" "FRX" 
# 
# $Retailing
# [1] "HD"   "LOW"  "PCLN" "TJX"  "DG"   "BBBY" "AZO"  "LTD" 
# [9] "ORLY" "GPC"  "GPS"  "KMX"  "KSS"  "TIF"  "PETM" "FDO" 
# [17] "GME"  "JCP" 
# 
# $`Semiconductors & Semiconductor Equipment`
# [1] "INTC" "TXN"  "AMAT" "MU"   "ADI"  "BRCM" "ALTR" "KLAC"
# [9] "LRCX" "LSI"  "TER"  "AMD" 
# 
# $`Software & Services`
# [1] "MSFT" "ORCL" "ADP"  "YHOO" "CTSH" "ADBE" "SYMC" "CTXS"
# [9] "PAYX" "CA"   "ADSK"
# 
# $`Technology Hardware & Equipment`
# [1] "TMO"  "WAT"  "PKI"  "IBM"  "CSCO" "QCOM" "EMC"  "GLW" 
# [9] "DELL" "SNDK" "NTAP" "STX"  "JNPR" "FFIV" "HRS"  "JBL" 
# [17] "FLIR"
# 
# $`Telecommunication Services`
# [1] "T"   "VZ"  "CTL"
# 
# $Transportation
# [1] "UNP"  "UPS"  "FDX"  "CSX"  "NSC"  "DAL"  "KSU"  "LUV" 
# [9] "CHRW" "EXPD" "R"   
# 
# $Utilities
# [1] "WMB" "KMI" "EQT" "DUK" "SO"  "EXC" "AEP" "PCG" "ED"  "EIX"
# [11] "XEL" "NU"  "DTE" "CNP" "WEC" "CMS" "SCG" "GAS" "POM" "TE"


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
number_of_ind <- sapply(industry_tickers,NROW)
number_of_ind

# You should get the following output, with names 
# of Industries as the element names:
# 
# Automobiles & Components                            Capital Goods 
#                       2                                       28 
# Commercial & Professional Services              Consumer Durables & Apparel 
#                       8                                       12 
# Diversified Financials                                   Energy 
#                       2                                       28 
# Food & Drug Retailing                 Food, Beverage & Tobacco 
#                       8                                       12 
# Health Care Equipment & Services             Hotels Restaurants & Leisure 
#                       21                                        5 
# Household & Personal Products                                Insurance 
#                       5                                        1 
# Materials                                    Media 
#                       22                                        7 
# Pharmaceuticals & Biotechnology                                Retailing 
#                       12                                       18 
# Semiconductors & Semiconductor Equipment                      Software & Services 
#                       12                                       11 
# Technology Hardware & Equipment               Telecommunication Services 
#                       17                                        3 
# Transportation                                Utilities 
#                       11                                       20


# Calculate a named list with the indices of companies in 
# each industry, and call it industry_indices. 
# The index of a company is its row number in panel_data.
# You can use functions sapply() and match(). 
# hint: you can use the vector industry_tickers.

### write your code here
industry_indices <- with(panel_data,sapply(industry_tickers, function(x){match(x,TICKER)}))
industry_indices
class(industry_indices) #"list"

# You should get the following output, with names 
# of Industries as the element names: 
# 
# > industry_indices
# 
# $`Automobiles & Components`
# [1]  9 21
# 
# $`Capital Goods`
# [1] 139 141 142 144 145 146 147 148 150 151 153 155 156 157 158 160 161 162 163 164 165 166 167 168 170 177 178
# [28] 179
# 
# $`Commercial & Professional Services`
# [1]  64 172 173 174 180 181 184 202
# 
# $`Consumer Durables & Apparel`
# [1]  10  14  23  25  27  33  36  37  38  99 182 228
# 
# $`Diversified Financials`
# [1] 101 102
# 
# $Energy
# [1] 66 67 68 69 70 71 72 73 74 75 77 79 80 81 82 83 85 86 87 88 90 91 92 93 94 95 96 98
# 
# $`Food & Drug Retailing`
# [1]  44  46  47  51  52  63  78 243
# 
# $`Food, Beverage & Tobacco`
# [1] 43 45 49 50 54 55 56 57 59 61 62 65
# 
# $`Health Care Equipment & Services`
# [1] 110 111 112 114 116 119 120 121 122 123 124 125 126 127 130 131 132 134 135 136 138
# 
# $`Hotels Restaurants & Leisure`
# [1]  4  8 20 31 35
# 
# $`Household & Personal Products`
# [1] 42 48 53 58 60
# 
# $Insurance
# [1] 100
# 
# $Materials
# [1]  97 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 203 204 205 206
# 
# $Media
# [1]  1  2 12 29 34 39 40
# 
# $`Pharmaceuticals & Biotechnology`
# [1] 103 104 105 106 107 108 109 115 117 118 128 129
# 
# $Retailing
# [1]  3  5  6  7 11 13 15 16 17 18 19 22 24 26 28 30 32 41
# 
# $`Semiconductors & Semiconductor Equipment`
# [1] 214 216 224 226 230 231 234 237 240 246 247 248
# 
# $`Software & Services`
# [1] 207 212 217 218 219 220 225 232 233 238 239
# 
# $`Technology Hardware & Equipment`
# [1] 113 133 137 208 211 213 215 221 223 227 229 235 236 241 242 244 245
# 
# $`Telecommunication Services`
# [1] 209 210 222
# 
# $Transportation
# [1] 140 143 149 152 154 159 169 171 175 176 183
# 
# $Utilities
# [1]  76  84  89 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265


# 8. (30pts) Calculate a named vector (not an array!) with 
# the average "NET.INCOME" of all the companies in each Sector. 
# You can use functions with(), split(), sapply(), and mean(). 

### write your code here
average_net_inc <- with(panel_data,sapply(split(NET.INCOME,Sector),mean))
average_net_inc
is.vector(average_net_inc) #[1] TRUE

# You should get the following output: 
# Consumer Discretionary    Consumer Staples    Energy    Financials
#   477.3523                  1033.5947       1542.6779     666.9847
# Health Care   Industrials     Information Technology    Materials
# 898.7432        490.6571            861.5079            216.4407
# Telecommunication Services      Utilities
#         1762.2357               229.8667

# extra

# Perform the same calculation as above, but using function tapply().
# You can also use functions with(), mean(), drop(), and as.matrix(). 
# You must produce a named vector, not an array!

### write your code here
average_net_inc2 <- tapply(X=as.vector(panel_data$NET.INCOME), INDEX=panel_data$Sector, FUN=mean)
average_net_inc2



# Perform a similar calculation for the average "NET.INCOME" 
# as above, but calculate a data frame using function aggregate().
# You can use functions with(), tapply(), mean(), drop(), 
# and as.matrix(). 

### write your code here
average_net_inc3 <- aggregate(x=panel_data$NET.INCOME,by=list(panel_data$Sector),FUN=mean)
average_net_inc3

# You should get the following output: 
#                       Group.1         x
# 1      Consumer Discretionary  477.3523
# 2            Consumer Staples 1033.5947
# 3                      Energy 1542.6779
# 4                  Financials  666.9847
# 5                 Health Care  898.7432
# 6                 Industrials  490.6571
# 7      Information Technology  861.5079
# 8                   Materials  216.4407
# 9  Telecommunication Services 1762.2357
# 10                  Utilities  229.8667

# extra end

# 9. (20pts) Calculate a vector of tickers of the 
# companies which have the highest ROE in each 
# Industry, and call it ticker_s.
# You must use functions sapply(), split(), 
# as.vector(), which.max(), and with().

### write your code here
ticker_s <- with(panel_data,sapply(split(ROE,Industry),function(x){TICKER[match(x[which.max(x)],ROE)]}))
unname(ticker_s)



# You should get the following output:
# > unname(ticker_s)
#  [1] "JCI"  "COL"  "PBI"  "COH"  "LUK"  "XOM"  "SYY"  "CPB"  "DVA" 
# [10] "YUM"  "AVP"  "MMC"  "IFF"  "OMC"  "MRK"  "AZO"  "INTC" "ORCL"
# [19] "STX"  "T"    "DAL"  "CNP"


# Subset the panel_data data frame and extract 
# the rows corresponding to ticker_s, and call it 
# max_roes.
# You can use function match(). 

### write your code here
max_roes <- with(panel_data,panel_data[match(ticker_s,TICKER),])
max_roes
max_roes[, c("COMPANY.NAME", "TICKER", "Industry")]

#if I use sort() to sort the sequence.
max_roes2 <- with(panel_data,panel_data[sort(match(ticker_s,TICKER)),])
max_roes2
max_roes2[, c("COMPANY.NAME", "TICKER", "Industry")]


# You should get the following output:
# 
# > max_roes[, c("COMPANY.NAME", "TICKER", "Industry")]
# 
#                         COMPANY.NAME TICKER                                 Industry
# 9               JOHNSON CONTROLS INC    JCI                 Automobiles & Components
# 168             ROCKWELL COLLINS INC    COL                            Capital Goods
# 184                 PITNEY BOWES INC    PBI       Commercial & Professional Services
# 14                         COACH INC    COH              Consumer Durables & Apparel
# 102           LEUCADIA NATIONAL CORP    LUK                   Diversified Financials
# 66                  EXXON MOBIL CORP    XOM                                   Energy
# 52                        SYSCO CORP    SYY                    Food & Drug Retailing
# 62                  CAMPBELL SOUP CO    CPB                 Food, Beverage & Tobacco
# 131                       DAVITA INC    DVA         Health Care Equipment & Services
# 8                     YUM BRANDS INC    YUM             Hotels Restaurants & Leisure
# 60                 AVON PRODUCTS INC    AVP            Household & Personal Products
# 100         MARSH & MCLENNAN COS INC    MMC                                Insurance
# 199 INTERNATIONAL FLAVORS & FRAG INC    IFF                                Materials
# 12                 OMNICOM GROUP INC    OMC                                    Media
# 105                   MERCK & CO INC    MRK          Pharmaceuticals & Biotechnology
# 15                      AUTOZONE INC    AZO                                Retailing
# 214                       INTEL CORP   INTC Semiconductors & Semiconductor Equipment
# 212                      ORACLE CORP   ORCL                      Software & Services
# 235               SEAGATE TECHNOLOGY    STX          Technology Hardware & Equipment
# 209                     A T & T CORP      T               Telecommunication Services
# 159              DELTA AIR LINES INC    DAL                           Transportation
# 259           CENTERPOINT ENERGY INC    CNP                                Utilities



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
dimen_sions <- do.call(rbind,eapply(rutils::env_etf,rm_dim)[bool])
dimen_sions

# You should get the following output:
#   dimen_sions
#          [,1] [,2]
# price_s  2437   20
# DBC      2437    6
# IEF      2437    6
# VTI      2437    6
# XLB      2437    6
# etc.


# extra

# 4. (20pts) Create a copy of rutils::env_etf in 
# your R Workspace, and call it env_etf.

### write your code here
 
ls(rutils::env_etf)
env_etf <- rutils::env_etf

# List the names of all the objects in env_etf 
# whose names start with "X*".
# You can use the functions glob2rx() and ls() with
# the "pattern" argument.

### write your code here
ls(env_etf)
ls(env_etf,pattern=glob2rx("X*"))


# You should get the following output:
# [1] "XLB" "XLE" "XLF" "XLI" "XLK" "XLP" "XLU" "XLV" "XLY"


# 5. (10pts) Remove all the objects in env_etf
# whose names start with "X*".
# You can use the functions rm(), glob2rx() and
# ls() with the "pattern" argument.

### write your code here
rm(envir = env_etf,list=ls(env_etf,pattern=glob2rx("X*")))


# List the names of all the objects in env_etf 

### write your code here
ls(env_etf)

# You should get the following output:
#  [1] "capm_stats"  "DBC"         "etf_list"    "IEF"         "IVW"        
#  [6] "IWB"         "IWD"         "IWF"         "price_s"     "re_turns"   
# [11] "risk_return" "sym_bols"    "VEU"         "VNQ"         "VTI"        
# [16] "VXX"         "VYM"

# extra end


