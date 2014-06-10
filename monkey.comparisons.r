require(gdata)
require(ggplot2)
require(reshape2)
library(doMC)
registerDoMC(12)
library(Morphometrics)

load("~/Dropbox/labbio/cov_bayes_data/Rdatas/nwm.Rdata")
cov.matrices = lapply(nwm.matrices, function(x) x[[1]])
raw.hipot = read.table("./hipoteses.funcionais.csv", header = F, sep = '\t', as.is = T)
hipoteses = list(Neuroface = as.matrix(raw.hipot[1:39,1:39]), All = as.matrix(raw.hipot[40:78,1:39]))
mat.list = c(cov.matrices, hipoteses)

comparisons.cor = KrzCor(cov.matrices)
comparisons.proj = KrzProjection(cov.matrices, num.cores = 12, full.results = F)
comp = melt(data.frame(comparisons.proj))
qplot(value, data = comp, geom = "histogram")

###########################
# Significance Testing
###########################

## Null hipotesis is dissimilarity

## 1) Random Matrices

x <- RandomMatrix(10)
y <- RandomMatrix(10)
x <- cov.matrices[[1]]
y <- cov.matrices[[2]]
KPS_Random <- function(x, y, iterations = 100){
    obs_sim <- KrzProjection(x, y)[[1]]
    dx <- diag(x)
    dy <- diag(y)
    rand_mats_x <- RandomMatrix(dim(x)[[1]], iterations, min(dx), max(dx))
    rand_mats_y <- RandomMatrix(dim(y)[[1]], iterations, min(dy), max(dy))
    random_comps <- unlist(Map(function(x, y) KrzProjection(x, y)[[1]], rand_mats_x, rand_mats_y))
    significance <- sum(random_comps > obs_sim)/iterations
    return(c(SharedVariance = obs_sim, Prob = significance))
}

