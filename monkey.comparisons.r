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

MultiKrzProjection <- function(mat.list,
                               ret.dim.1 = NULL, ret.dim.2 = NULL,
                               num.cores = 1,
                               full.results = F){
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
    if(full.results){
    CompareToNProj <- function(n) llply(mat.list,
                                        function(x) {KrzProjection(x,
                                                                   mat.list[[n]],
                                                                   ret.dim.1, ret.dim.2)},
                                        .parallel = parallel)
    }
    else{
    CompareToNProj <- function(n) llply(mat.list[names(mat.list) != n],
                                        function(x) {KrzProjection(x,
                                                                   mat.list[[n]],
                                                                   ret.dim.1, ret.dim.2)[1]},
                                        .parallel = parallel)
    }
    comparisons.proj  <- llply(names(mat.list),
                               CompareToNProj,
                               .progress="text", .parallel = parallel)
    if(full.results){
        names(comparisons.proj) = names(mat.list)
        return(comparisons.proj)
    }
    else{
        comparisons.proj <- melt(comparisons.proj)
        comparisons.proj[,4] = names(mat.list)[(comparisons.proj[,4])]
    }
    return(comparisons.proj[,-2])
}


comparisons.cor = MultiKrzCor(cov.matrices)
comparisons.proj = MultiKrzProjection(cov.matrices, num.cores = 12, full.results = F)
qplot(value, data = comparisons.proj, geom = "histogram")
