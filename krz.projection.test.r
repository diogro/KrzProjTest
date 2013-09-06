source("~/projects/lem-usp-R/matrix.func.r")
require(plyr)
require(ggplot2)
require(gdata)

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

comparisons.proj <- MultiKrzProjection(random.matrices, num.cores = 12)
qplot(value, data=comparisons.proj, geom="histogram")

comparisons.cor <- MultiKrzProjection(random.matrices, num.cores = 12)
qplot(value, data=comparisons.cor, geom="histogram")

out.comparisons = data.frame(Kzr_Correlation = comparisons.cor[[1]], Kzr_Projection = comparisons.proj[[1]])
write.csv(out.comparisons, "KzrCorr_Projection.csv")
