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
KPS_RandomMatrix <- function(x, y, iterations = 100){
    obs_sim <- KrzProjection(x, y)[[1]]
    rand_mats_x <- RandomMatrix(dim(x)[[1]], iterations, variance = diag(x))
    rand_mats_y <- RandomMatrix(dim(y)[[1]], iterations, variance = diag(y))
    random_comps <- unlist(Map(function(x, y) KrzProjection(x, y)[[1]], rand_mats_x, rand_mats_y))
    significance <- sum(random_comps > obs_sim)/iterations
    return(c(SharedVariance = obs_sim, Prob = significance))
}

## 2) Random EigenVectors


## 3) Shuffled eigenvectors

KPS_ShuffleEigenVector <- function(x, y, iterations = 100){
    n = dim(x)[1]
    obs_sim <- KrzProjection(x, y)[[1]]
    eigen_x <- eigen(x)
    eigen_y <- eigen(y)
    randomCompFunc = function(index){
        shuffle_x = sample(1:n, n)
        shuffle_y = sample(1:n, n)
        rand_x = t(eigen_x$vectors[,shuffle_x]) %*% diag(eigen_x$values) %*% eigen_x$vectors[,shuffle_x]
        rand_y = t(eigen_y$vectors[,shuffle_y]) %*% diag(eigen_y$values) %*% eigen_y$vectors[,shuffle_y]
        KrzProjection(rand_x, rand_y)[[1]]
    }
    random_comps = aaply(1:iterations, 1, randomCompFunc)
    significance <- sum(random_comps > obs_sim)/iterations
    return(c(SharedVariance = obs_sim, Prob = significance))
}


## Null hipotesis is similarity

## 1) Random populations are shuffled

KPS_ShufflePop = function(x, y, sample.x = 100, sample.y = 100, iterations = 100) {
    n = dim(x)[1]
    obs_sim <- KrzProjection(x, y)[[1]]
    randomCompFunc = function(index){
        pop_x = data.frame(mvtnorm::rmvnorm(n = sample.x, sigma = x))
        pop_x$'pop' = 'x'
        pop_y = data.frame(mvtnorm::rmvnorm(n = sample.y, sigma = y))
        pop_y$'pop' = 'y'
        pop = rbind(pop_x, pop_y)
        shuffle = sample(dim(pop)[1])
        pop$'pop' = pop$'pop'[shuffle]
        rand.P = dlply(pop, .(pop), function(x) cov(x[-dim(x)[2]]))
        KrzProjection(rand.P)[1,2]
    }
    random_comps = aaply(1:iterations, 1, randomCompFunc)
    significance <- sum(random_comps < obs_sim)/iterations
    return(c(SharedVariance = obs_sim, Prob = significance))
}
