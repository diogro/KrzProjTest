source("~/projects/lem-usp-R/matrix.func.r")
require(plyr)
require(ggplot2)
require(reshape2)
library(doMC)
library(foreach)
registerDoMC(10)


p = 39 #Traits
n = 1 #Matrices
sizematrices = qr.Q(qr(matrix(c(rep(1, p), rnorm(p*(p-1))), p, p)))
for(i in 1:(n-1)){
    current.mat = qr.Q(qr(matrix(c(rep(1, p), rnorm(p*(p-1))), p, p)))
    sizematrices = rbind(sizematrices, current.mat)
}
nosizematrices = qr.Q(qr(matrix(rnorm(p*p), p, p)))
for(i in 1:(n-1)){
    current.mat = qr.Q(qr(matrix(rnorm(p*p), p, p)))
    nosizematrices = rbind(nosizematrices, current.mat)
}
#write.csv(sizematrices, "size-39t-eigen-vectors.csv")
#write.csv(nosizematrices, "nosize-39t-eigen-vectors.csv")

file.list = dir("./random_matrices/", "matrix-")
random.matrices = llply(file.list, function(x) as.matrix(read.table(paste("./random_matrices", x, sep = "/"), sep=',')))
CompareToN <- function(n) laply(random.matrices[-n], function(x) {KrzProjection(x, random.matrices[[n]])[1]}, .parallel = T)
comparisons  <- alply(1:100, 1, CompareToN)
comparisons  <-  melt(comparisons)
qplot(value, data=comparisons, geom="histogram")

