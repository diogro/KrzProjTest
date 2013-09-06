source("~/projects/lem-usp-R/matrix.func.r")
require(plyr)
require(gdata)
require(ggplot2)
require(reshape2)
library(doMC)
library(foreach)
registerDoMC(12)

load("~/Dropbox/labbio/cov_bayes_data/nwm.Rdata")
cov.matrices = lapply(nwm.matrices, function(x) x[[1]])
raw.hipot = read.table("./hipoteses.funcionais.csv", header = F, sep = '\t', as.is = T)
hipoteses = list(Neuroface = as.matrix(raw.hipot[1:39,1:39]), All = as.matrix(raw.hipot[40:78,1:39]))
mat.list = c(cov.matrices, hipoteses)

MultiKrzProjection  <- function(mat.list, ret.dim.1 = NULL, ret.dim.2 = NULL, num.cores = 1){
    require(plyr)
    require(reshape2)
    if (num.cores > 1) {
        library(doMC)
        library(foreach)
        registerDoMC(num.cores)
        parallel = TRUE
    }
    else
        parallel = FALSE
    CompareToNProj <- function(n) laply(mat.list[-n],
                                        function(x) {KrzProjection(x,
                                                                   mat.list[[n]],
                                                                   ret.dim.1, ret.dim.2)[1]},
                                        .parallel = parallel)
    comparisons.proj  <- alply(1:length(mat.list), 1,
                               CompareToNProj,
                               .progress="text", .parallel = parallel)
    comparisons.proj  <-  melt(comparisons.proj)
    return(comparisons.proj)
}

comparisons.proj = MultiKrzProjection(mat.list, num.cores = 12 )
qplot(value, data = comparisons.proj, geom = "histogram")
