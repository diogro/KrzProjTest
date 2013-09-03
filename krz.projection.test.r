source("~/projects/lem-usp-R/matrix.func.r")
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
write.csv(sizematrices, "size-39t-eigen-vectors.csv")
write.csv(nosizematrices, "nosize-39t-eigen-vectors.csv")

