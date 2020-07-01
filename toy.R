
library(MASS)

source("foos.R")

n <- 500
p <- 300
alpha <- .05 # significance level
sigma <- 1 # noise level
maf <- 0.6 # minor allele frequency (MAF)


p1 <- 5 ### block size
rho <- .5
M <- matrix(rho, p1, p1)  
diag(M) <- 1
beta <- rep(0, p)


set.seed(1)
# threshold in power enhancement
# simulate snp matrix
G1 <- matrix(0, n, p)
for(i in 1:(p/p1)){
  G0 <- mvrnorm(n, rep(0, p1), M)
  G1[, (1 + (i - 1) * p1):(p1 * i)] <- G0
}
G <- matrix(1, n, p)
i1 <- which(G1 < qnorm(maf^2))
i3 <- which(G1 > qnorm(1 - (1 - maf)^2))
G[i1] <- 0
G[i3] <- 2
W <- rnorm(n, 0.1 * G[,1], 1)
X <- cbind(1, W)

delta <- caldelta(X, G, simple = TRUE, B = 1000, seed = 123)

e <- rnorm(n, 0, sigma)
Y <- 1 + W + G %*% beta + e
X <- cbind(1, W)
d <- ncol(X)

if(p >= n / 2){
  sigma2_est <- vrcv(Y, cbind(X[, -1], G) , Num = 10, 
                     seed = 1)
}else{
  fit <- lm(Y ~ X + G - 1)
  sigma2_est <- (summary(fit)$sigma)^2
}


ts <- compt(Y, X, G, sigma2 = sigma2_est, delta = delta)

# the p-value of the Th test
ts$ph

# the p-value of the T test
ts$phe




