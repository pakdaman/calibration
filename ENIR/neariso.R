#! /usr/bin/Rscript

library(neariso)

args <- commandArgs(trailingOnly = TRUE)

#Extract directory name
dname = args[1]
setwd(dname)

#Extract base name
bname = args[2]

# Loading z vector by reading the associated data file
zname = paste(bname,"_z.csv", sep = "")
tmp <- read.csv(zname, header=FALSE)
z <- tmp$V1

# loading y vector by reading the associated data file
yname = paste(bname,"_y.csv", sep = "")
tmp <- read.csv(yname, header=FALSE)
y <- tmp$V1

# Calling neariso to find a solution path of near isotonic regression mappings
#res <- neariso(z, lambda=NULL)
res <- neariso2(z, y, lambda=NULL)
lambda =  res$lambda
fname = paste(bname,"_lambda.csv", sep = "")
write.table(lambda, file=fname, row.names=FALSE, col.names=FALSE, sep=",")

beta = res$beta
fname = paste(bname,"_beta.csv", sep = "")
write.table(beta, file=fname, row.names=FALSE, col.names=FALSE, sep=",")

df = res$df
fname = paste(bname,"_df.csv", sep = "")
write.table(df, file=fname, row.names=FALSE, col.names=FALSE, sep=",")


# Garbage
#zname = 'data_z.csv'
#setwd('/Users/mahdi/Research/Calibration/Thesis/Code/R_Matlab/CSV/')
#zname = 'data_z.csv'
#yname = 'data_y.csv'
#tmp <- read.csv(yname, header=FALSE)
#y <- tmp$V1
#tmp <- read.csv(zname, header=FALSE)
#z <- tmp$V1



