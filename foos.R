
#**********************************************************
#
# calucate TH sign for the power enhancement component 
#
#**********************************************************

tsign <- function(x){
  ifelse(x >= 0, 1, -1)
}

#**********************************************************
#
# calucate delta_p for power enhancement 
#
#**********************************************************

caldelta <- function(X, G, B = 2000, seed = 1, simple = FALSE){
  
  set.seed(seed)
  n <- nrow(X)
  p <- ncol(G)
  
  if(!simple){
    n1 <- ceiling(n/2)
    
    bsf <- function(){
      idx <- sample(1:n, size = n1)
      X1 <- X[idx, ]
      G1 <- G[idx, ]
      e <- rnorm(n1)
      Px <- X1%*% solve(t(X1) %*% X1) %*% t(X1)
      IPx <- diag(1, n1) - Px
      halfA <- t(G1) %*% IPx
      vecD <- 1 / diag(halfA %*% t(halfA))
      halfA2 <- sweep(halfA, 1, sqrt(vecD), FUN = "*")
      betam <- halfA2 %*% e
      max(abs(betam))
    }
    tmp <- replicate(B, bsf())
  }else{
    tmp <- replicate(B, max(abs(rnorm(p))))
  }
  
  return(max(tmp))
  
}


#****************************************************************
#
# compute the TL, TH, T statistics in the paper
# the first column of X must be intercept 
# method = low, use TL; method = high, use TH
# ape: add power enhancement 
# delta: threshold in power enhancement
#
#****************************************************************

compt <- function(Y, X, G, sigma2 = NULL, method = "high", 
                  ape = TRUE, delta = NULL, 
                  boncor = TRUE, edgecor = FALSE, ...){
  
  n <- nrow(X)
  d <- ncol(X)
  p <- ncol(G)
  
  # use the estimator under the null if sigma2 not provided
  if(is.null(sigma2)){
    fit <- lm(Y ~ X - 1)
    sigma2 <- (summary(fit)$sigma)^2
  }
  
  Px <- X %*% solve(t(X) %*% X) %*% t(X)
  IPx <- diag(1, n) - Px
  halfA <- t(G) %*% IPx
  vecD <- 1 / diag(halfA %*% t(halfA))
  halfA2 <- sweep(halfA, 1, vecD, FUN = "*")
  A <- t(halfA) %*% halfA2
  
  diff <- NA # difference of p^-1||Sigma||_F^2 - p/n
  
  if(method == "low"){
    A1 <- A
    Q <- c(t(Y) %*% A1 %*% Y)
    sd1 <- sqrt(2 * sum(A1^2))
    TH <- (Q - p * sigma2)/(sd1 * sigma2)
  }else if(method == "high"){
    A1 <- A - p / (n - d) * IPx
    fnorm <- sum(A1^2)
    sd1 <- sqrt(2 * fnorm)
    Q1 <- c(t(Y) %*% A1 %*% Y)
    TH <- Q1 / (sd1 * sigma2)
    diff <- fnorm / p
  }
  ph <- 2 * pnorm(abs(TH), lower.tail = FALSE) # pvalue for TH
  
  
  Tpe <- NA # test statistic for power enhancement
  phe <- NA # p-value for Tpe
  pb <- NA # p-value for bonferroni correction
  if(ape){
    
    if(is.null(delta)){
      delta <- log(log(n)) * sqrt(log(p))
    }
    
    # marginal coefficients and sds
    betam <- halfA2 %*% Y
    sdm <- sqrt(sigma2 * vecD)
    
    hatS <- which(abs(betam) > sdm * delta)
    
    if(length(hatS) >= 1){
      T0 <- tsign(TH) * sqrt(p) * sum( betam[hatS]^2 / (sdm[hatS])^2 )
      Tpe <- TH + T0
    }else{
      Tpe <- TH
    }
    
    if(boncor){
      Tb <- max(abs(betam) / sdm)
      pb <- 2 * pnorm(Tb, lower.tail = FALSE) * p
    }
    
    phe <- 2 * pnorm(abs(Tpe), lower.tail = FALSE)
    
  }
  
  # edgeworth correction
  ph_edge <- rep(NA, 2)
  phe_edge <- rep(NA, 2)
  if(edgecor){
    
    THa <- abs(TH)
    M1 <- sum(diag(A1))
    M2 <- 2 * sum(A1^2)
    M3 <- 8 * sum(diag(A1 %*% A1 %*% A1))
    M4 <- 48 * sum(diag(A1 %*% A1 %*% A1 %*% A1))
    
    E1 <- pnorm(THa) + M3 * (1 - THa^2) * dnorm(THa) / 6 / (sqrt(M2))^3
    ph_edge[1] <- 2 * (1 - E1)
    
    E2 <- E1 - M4 * (THa^3 - 3 * THa) * dnorm(THa) / 24 / M2^2 - 
      M3^2 * (THa^5 - 10 * THa^3 + 15 * THa) * dnorm(THa) / 72 / M2^3
    
    ph_edge[2] <- 2 * (1 - E2)
    
    if(ape){
      if(length(hatS) >= 1){
        Tpea <- abs(Tpe)
        Ee1 <- pnorm(Tpea) + M3 * (1 - Tpea^2) * dnorm(Tpea) / 6 / (sqrt(M2))^3
        phe_edge[1] <- 2 * (1 - Ee1)
        
        Ee2 <- Ee1 - M4 * (Tpea^3 - 3 * Tpea) * dnorm(Tpea) / 24 / M2^2 - 
          M3^2 * (Tpea^5 - 10 * Tpea^3 + 15 * Tpea) * dnorm(Tpea) / 72 / M2^3
        
        phe_edge[2] <- 2 * (1 - Ee2)
        
        
      }else{
        phe_edge <- ph_edge
      }
    }
    
    
  }

  
  res <- list(TH = TH, ph = ph, sd1 = sd1, Tpe = Tpe, phe = phe, 
              pb = pb, ph_edge = ph_edge, phe_edge = phe_edge, 
              diff = diff)
  
  res
}


#*********************************************************************
#
# High dimensional variance estimation using refitted cross vadiation: 
# 1. Divide n samples into 2 sub-samples with sizes n1 = n / 2
# 2. Using SIS to select n1 * sis_prop variables if p > n
# No intercept is need in X
#
#**********************************************************************

# use replicate instead of for loop 
vrcv <- function(Y, X, Num = 10, sis_prop = 0.5, seed = NULL){
  
  n <- length(Y)
  p <- ncol(X)
  
  if(!is.null(seed) & is(seed, "numeric")){
    set.seed(seed)
  }
  
  n1 <- floor(n / 2)
  
  res <- replicate(Num, {
    cv1 <- sample(n, n1)
    X1 <- X[cv1, ]
    Y1 <- Y[cv1]
    X2 <- X[-cv1, ]
    Y2 <- Y[-cv1]

    
    # SIS
    cor1 <- cor(X1, Y1)
    cor2 <- cor(X2, Y2)
    
    n2 <- n1 * sis_prop
    
    ind1 <- order(abs(cor1), decreasing = TRUE)[1:n2]
    ind2 <- order(abs(cor2), decreasing = TRUE)[1:n2]
    
    X12 <- cbind(1, X1[, ind2])
    X22 <- cbind(1, X2[, ind1])
    
    # fit1 <- lm(Y1 ~ X12 - 1)
    # fit2 <- lm(Y2 ~ X22 - 1)
    # 
    # s1 <- (summary(fit1)$sigma)^2
    # s2 <- (summary(fit2)$sigma)^2
    
    IP1 <- diag(1, n1) - X12 %*% solve(t(X12) %*% X12) %*% t(X12)
    IP2 <- diag(1, n - n1) - X22 %*% solve(t(X22) %*% X22) %*% t(X22)
    
    s1 <- (t(Y1) %*% IP1 %*% Y1) / (n1 - n2 - 1)
    s2 <- (t(Y2) %*% IP2 %*% Y2) / (n - n1 - n2 - 1)
    
    (s1 + s2) / 2
    
  })
  
  mean(res)
  
}




