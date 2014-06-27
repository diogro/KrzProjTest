require(gdata)
require(ggplot2)
require(reshape2)
library(doMC)
registerDoMC(12)
library(Morphometrics)

load("~/Dropbox/labbio/cov_bayes_data/Rdatas/nwm.Rdata")
cov.matrices = lapply(nwm.matrices, function(x) x[[1]])

###########################
# Significance Testing
###########################

funcSignificance <- function(x, y, randomCompFunc, ..., compFunc=KrzProjection, iterations = 100) {
    obs_sim <- compFunc(x, y)[[1]]
    random_comps = aaply(1:iterations, 1, function (i) randomCompFunc(x, y, compFunc, ...), .parallel = TRUE)
    significance <- sum(random_comps > obs_sim)/iterations
    return (c(SharedVariance = obs_sim, Prob = significance))
}
## Null hipotesis is dissimilarity


## 1) Random Matrices

RandomMat <- function(x, y, compFunc=KrzProjection, ...){
    n = dim(x)[1]
    rand_x <- RandomMatrix(n, variance = diag(x))
    rand_y <- RandomMatrix(n, variance = diag(y))
    return (compFunc(rand_x, rand_y, ...)[[1]])
}

## 2) Random EigenVectors

RandomEV <- function(x, y, compFunc=KrzProjection, evals=NULL, size = TRUE){
    n = dim(y)[1]

    if (is.null(evals))
        evals = diag(eigen(x)$values)

    if (size)
        random_eVecs = qr.Q(qr(matrix(c(rep(1, n), rnorm(n*(n-1))), n, n)))
    else
        random_eVecs = qr.Q(qr(matrix(rnorm(n*n), n, n)))

    rand_x =  t(random_eVecs) %*% evals %*% random_eVecs
    return (compFunc(rand_x, y)[[1]])
}

## 3) Shuffled eigenvectors

ShuffleEV <- function(x, y, compFunc=KrzProjection, eigen_x=NULL, eigen_y=NULL){
    n = dim(x)[1]

    if (is.null(eigen_x))
        eigen_x <- eigen(x)

    if (is.null(eigen_y))
        eigen_y <- eigen(y)

    shuffle_x = sample(1:n, n)
    shuffle_y = sample(1:n, n)
    rand_x = t(eigen_x$vectors[,shuffle_x]) %*% diag(eigen_x$values) %*% eigen_x$vectors[,shuffle_x]
    rand_y = t(eigen_y$vectors[,shuffle_y]) %*% diag(eigen_y$values) %*% eigen_y$vectors[,shuffle_y]
    return(compFunc(rand_x, rand_y)[[1]])
}

## Null hypotesis is similarity

## 1) Random populations are shuffled

ShufflePop <- function(x, y, compFunc=KrzProjection, sample.x = 100, sample.y = 100){
    pop_x = data.frame(mvtnorm::rmvnorm(n = sample.x, sigma = x))
    pop_x$'pop' = 'x'
    pop_y = data.frame(mvtnorm::rmvnorm(n = sample.y, sigma = y))
    pop_y$'pop' = 'y'
    pop = rbind(pop_x, pop_y)
    shuffle = sample(dim(pop)[1])
    pop$'pop' = pop$'pop'[shuffle]
    rand.P = dlply(pop, .(pop), function(x) cov(x[-dim(x)[2]]))
    return (compFunc(rand.P)[1,2])
}

rand_x = RandomMatrix(15, 2, 1, 10)
rand_y = RandomMatrix(15, 2, 1, 10)

sig_dist_RandMat = Map(function(x, y) funcSignificance(x, y, RandomMat, compFunc=PCAsimilarity), rand_x, rand_y)
kr_sig_dist_RandMat = Map(function(x, y) funcSignificance(x, y, RandomMat), rand_x, rand_y)

sig_dist_RandEV = Map(function(x, y) funcSignificance(x, y, RandomEV, diag(eigen(x)$values), compFunc=PCAsimilarity), rand_x, rand_y)
kr_sig_dist_RandEV = Map(function(x, y) funcSignificance(x, y, RandomEV, diag(eigen(x)$values)), rand_x, rand_y)

sig_dist_ShuffleEV = Map(function(x, y) funcSignificance(x, y, ShuffleEV, eigen(x), eigen(y), compFunc=PCAsimilarity), rand_x, rand_y)
kr_sig_dist_ShuffleEV = Map(function(x, y) funcSignificance(x, y, ShuffleEV, eigen(x), eigen(y)), rand_x, rand_y)

sig_dist_ShufflePop = Map(function(x, y) funcSignificance(x, y, ShufflePop, compFunc=PCAsimilarity), rand_x, rand_x)
kr_sig_dist_ShufflePop = Map(function(x, y) funcSignificance(x, y, ShufflePop), rand_x, rand_x)

save.image("./power.Rdata")
