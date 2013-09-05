source("~/projects/lem-usp-R/matrix.func.r")
require(plyr)
require(gdata)
require(ggplot2)
require(reshape2)
library(doMC)
library(foreach)
registerDoMC(10)

#file.list = dir("./random_matrices/", "matrix-.+.csv")
#random.matrices = llply(file.list,
                        #function(x) as.matrix(read.table(paste("./random_matrices",
                                                               #x,
                                                               #sep = "/"),
                                                         #sep=',')))

file.list = dir("./random_matrices/", "matrix-.+.xls")
random.matrices = llply(file.list,
                        function(x) as.matrix(read.xls(paste("./random_matrices",
                                                               x,
                                                               sep = "/"), header=F)))

CompareToNProj <- function(n) laply(random.matrices[-n], function(x) {KrzProjection(x, random.matrices[[n]])[1]}, .parallel = T)
comparisons.proj  <- alply(1:100, 1, CompareToNProj, .progress="text")
comparisons.proj  <-  melt(comparisons.proj)
qplot(value, data=comparisons.proj, geom="histogram")

CompareToNCor <- function(n) laply(random.matrices[-n], function(x) {KrzCor(x, random.matrices[[n]])}, .parallel = T)
comparisons.cor  <- alply(1:100, 1, CompareToNCor, .progress="text")
comparisons.cor  <-  melt(comparisons.cor)
qplot(value, data=comparisons.cor, geom="histogram")
