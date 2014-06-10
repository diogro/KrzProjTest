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

## Null hipotesis is dissimilarity

## 1) Random Matrices

KPS_RandomMatrix <- function(x, y, iterations = 100){
    obs_sim <- KrzProjection(x, y)[[1]]
    randomCompFunc = function(index){
        rand_x <- RandomMatrix(dim(x)[[1]], variance = diag(x))
        rand_y <- RandomMatrix(dim(y)[[1]], variance = diag(y))
        KrzProjection(rand_x, rand_y)[[1]]
    }
    random_comps = aaply(1:iterations, 1, randomCompFunc)
    significance <- sum(random_comps > obs_sim)/iterations
    return(c(SharedVariance = obs_sim, Prob = significance))
}

## 2) Random EigenVectors

KPS_RandomEV <- function(x, y, iterations = 100, size = TRUE){
    obs_sim <- KrzProjection(x, y)[[1]]
    n = dim(x)[1]
    evals <- diag(eigen(x)$values)
    randomCompFunc = function(index){
        if(size) random_eVecs = qr.Q(qr(matrix(c(rep(1, n), rnorm(n*(n-1))), n, n)))
        else     random_eVecs = qr.Q(qr(matrix(rnorm(n*n), n, n)))
        rand_x =  t(random_eVecs) %*% evals %*% random_eVecs
        KrzProjection(rand_x, y)[[1]]
    }
    random_comps = aaply(1:iterations, 1, randomCompFunc)
    significance <- sum(random_comps > obs_sim)/iterations
    return(c(SharedVariance = obs_sim, Prob = significance))
}

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

KPS_ShufflePop = function(x, y, sample.x = 100, sample.y = 100, iterations = 100){
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
