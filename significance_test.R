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

KrzProjSignificance <- function(x, y, randomCompFunc = c("RandomMat", "RandomEV", "ShuffleEV", "ShufflePop"), ..., iterations = 100){
    randomCompFunc = match.arg(randomCompFunc)
    n = dim(x)[1]
    if(randomCompFunc == "ShuffleEV"){
        eigen_x <- eigen(x)
        eigen_y <- eigen(y)
    }
    if(randomCompFunc == "RandomEV"){
        evals <- diag(eigen(x)$values)
    }
    obs_sim <- KrzProjection(x, y)[[1]]
    randomCompFunc = match.fun(randomCompFunc)
    environment(randomCompFunc) = environment()
    random_comps = aaply(1:iterations, 1, randomCompFunc, ..., .parallel = TRUE)
    significance <- sum(random_comps > obs_sim)/iterations
    return(c(SharedVariance = obs_sim, Prob = significance))
}
rand_x = RandomMatrix(39, 100, 1, 10)
rand_y = RandomMatrix(39, 100, 1, 10)
sig_dist_RandMat = Map(function(x, y) KrzProjSignificance(x, y, 'RandomMat'), rand_x, rand_y)
sig_dist_RandEV = Map(function(x, y) KrzProjSignificance(x, y, 'RandomEV'), rand_x, rand_y)
sig_dist_ShuffleEV = Map(function(x, y) KrzProjSignificance(x, y, 'ShuffleEV'), rand_x, rand_y)
sig_dist_ShufflePop = Map(function(x, y) KrzProjSignificance(x, y, 'ShufflePop'), rand_x, rand_x)
save.image("./power.Rdata")


## Null hipotesis is dissimilarity

## 1) Random Matrices

RandomMat <- function(index){
    rand_x <- RandomMatrix(n, variance = diag(x))
    rand_y <- RandomMatrix(n, variance = diag(y))
    KrzProjection(rand_x, rand_y)[[1]]
}

## 2) Random EigenVectors

RandomEV <- function(index, size = TRUE){
    if(size) random_eVecs = qr.Q(qr(matrix(c(rep(1, n), rnorm(n*(n-1))), n, n)))
    else     random_eVecs = qr.Q(qr(matrix(rnorm(n*n), n, n)))
    rand_x =  t(random_eVecs) %*% evals %*% random_eVecs
    KrzProjection(rand_x, y)[[1]]
}

## 3) Shuffled eigenvectors

ShuffleEV <- function(index){
    shuffle_x = sample(1:n, n)
    shuffle_y = sample(1:n, n)
    rand_x = t(eigen_x$vectors[,shuffle_x]) %*% diag(eigen_x$values) %*% eigen_x$vectors[,shuffle_x]
    rand_y = t(eigen_y$vectors[,shuffle_y]) %*% diag(eigen_y$values) %*% eigen_y$vectors[,shuffle_y]
    KrzProjection(rand_x, rand_y)[[1]]
}

## Null hipotesis is similarity

## 1) Random populations are shuffled

ShufflePop <- function(index, sample.x = 100, sample.y = 100){
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
