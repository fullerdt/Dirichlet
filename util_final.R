

sirt_digamma1 <- function(x, h=1e-3) #Hidden utility function. I couldn't access it through R 
  #but I found it in the package on Github. 
  #It's an implementation of the trigamma function (via numerical derivative of digamma)
{
  ( digamma(x+h) - digamma(x-h) ) / (2*h)
}


MLE <- function (x, weights = NULL, eps = 10^(-5), convcrit = 1e-05, 
                 maxit = 1000, oldfac = 0.3, progress = FALSE) #The mle estimations function
{
  N <- nrow(x)
  K <- ncol(x)
  x <- (x + eps)/(1 + 2 * eps)
  x <- x/rowSums(x)
  N <- nrow(x)
  if (is.null(weights)) {
    weights <- rep(1, N)
  }
  weights <- N * weights/sum(weights)
  log.pbar <- colMeans(weights * log(x))
  alphaprob <- colMeans(x * weights)
  p2 <- mean(x[, 1]^2 * weights)
  xsi <- (alphaprob[1] - p2)/(p2 - (alphaprob[1])^2)
  alpha <- xsi * alphaprob
  K1 <- matrix(1, K, K)
  conv <- 1
  iter <- 1
  while ((conv > convcrit) & (iter < maxit)) {    #This is the stepwise iteration that converges at the global max     
    alpha0 <- alpha
    g <- N * digamma(sum(alpha)) - N * digamma(alpha) + N * 
      log.pbar
    z <- N * sirt_digamma1(sum(alpha))
    H <- diag(-N * sirt_digamma1(alpha)) + z
    alpha <- alpha0 - solve(H, g)
    alpha[alpha < 0] <- 1e-10
    alpha <- alpha0 + oldfac * (alpha - alpha0)
    conv <- max(abs(alpha0 - alpha))
    if (progress) {
      print(paste(iter, sum(alpha), conv))
      utils::flush.console()
    }
    iter <- iter + 1
  }
  alpha0 <- sum(alpha)
  xsi <- alpha/alpha0
  res <- list(alpha = alpha, alpha0 = alpha0, xsi = xsi)
  return(res)
}


diri.nr <- function(x, type = 1, tol = 1e-07) {
  
  if ( type == 1 ) {
    runtime <- proc.time()
    n <- dim(x)[1]  ## the sample size
    p <- dim(x)[2]  ## dimensionality
    m <- Rfast::colmeans(x)
    zx <- t( log(x) )
    down <-  - sum( m * ( Rfast::rowmeans( zx ) - log(m) ) )
    sa <- 0.5 * (p - 1) / down  ## initial value for precision
    a1 <- sa * m  ## initial values
    gm <- Rfast::rowsums(zx)
    z <- n * digamma( sa )
    g <- z - n * digamma(a1) + gm
    qk <-  - n * trigamma(a1)
    b <- sum(g / qk) / ( 1/z - sum(1 / qk) )
    a2 <- a1 - (g - b)/qk
    i <- 2
    while ( sum( abs( a2 - a1 ) ) > tol ) {
      i <- i + 1
      a1 <- a2
      z <- n * digamma( sum(a1) )
      g <- z - n * digamma(a1) + gm
      qk <-  - n * trigamma(a1)
      b <- sum(g / qk) / ( 1/z - sum(1 / qk) )
      a2 <- a1 - (g - b) / qk
    }
    loglik <- n * lgamma( sum(a2) ) - n * sum( lgamma(a2) ) + sum( zx * (a2 - 1) )
    runtime <- proc.time() - runtime
    
    if ( is.null(colnames(x)) ) {
      names(a2) <- paste("X", 1:p, sep = "")
    } else  names(a2) <- colnames(x)
    res <- list(iter = i, loglik = loglik, param = a2, runtime = runtime)
    
  } else  res <- Rfast::diri.nr2(x, tol = tol)
  
  res
}

dirimean.test <- function(x, a) {
  ## x is the compositional data
  ## a is the hypothesized compositional mean vector
  n <- dim(x)[1]  ## sample size
  d <- dim(x)[2] #- 1
  
  if ( min(x) <= 0 || min(a) <= 0 ) {
    res <- paste("There are zeros in the data")
  } else {
    z <- t( log(x) )
    ## loglik is for the 'mle' type
    loglik <- function(phi) {
      phi <- exp(phi)
      ma <- phi * a
      n * lgamma( phi ) - n * sum( lgamma(ma) ) + sum( z * (ma - 1) )
    }
    ## phi under Ho
    phi <- sum(a)
    #    if ( phi == 1 ) {
    #      mod0 <- optimize(loglik, c(-20, 20), maximum = TRUE )
    #      ell0 <- mod0$objective
    #      phi0 <- exp( mod0$maximum )
    #      par0 <- phi0 * a
    #    } else if ( phi > 1 ) {
    ell0 <-  n * lgamma( phi ) - n * sum( lgamma(a) ) + sum( z * ( a - 1 )  )
    par0 <- a
    #    }
    ## parameters under H1
    mod1 <- diri.nr(x)
    ell1 <- mod1$loglik
    ## test statistic and p-value
    test <- 2 * (ell1 - ell0)
    pvalue <- pchisq(test, d, lower.tail = FALSE)
    param <- rbind(par0, mod1$param)
    rownames(param) <- c("Null", "Alternative")
    if ( is.null( colnames(x) ) ) {
      colnames(param) <- paste("X", 1:c( d ), sep = " " )
    } else  colnames(param) <- paste("X", 1:c( d ), sep = " " )
    
    lik <- c(ell0, ell1)
    names(lik) <- c("Null loglik", "Alternative loglik")
    info <- c(test, pvalue)
    names(info) <- c("Test", "p-value")
    res <- list(param = param, loglik = lik, info = info)
  }
  
  res
}

#function to find fisher information
DirFIMPO <- function(delta){
  FIMPO <- matrix(nrow=length(delta), ncol=length(delta), -trigamma(sum(delta)))
  diagonal <- c()
  for(i in 1:length(delta)){
    diagonal <- c(diagonal, trigamma(delta[i])-trigamma(sum(delta)))
  }
  diag(FIMPO) <- diagonal
  return(FIMPO)
}

#function to for Rao score test

RaoScore <- function(delta,X,Fisher) {
  #RaoScore = numeric(2) # when we want to create a function with multiple outcomes
  #names(RaoScore) = c("p-value","test_stat") #name outcomes
  n <- dim(X)[1]
  d <- dim(X) [2]
  G <- log(X)
  T <- colSums(G)
  gradient_list <- rep(0, d)
  for(i in 1:d){
    gradient_list[i] <- n %*% digamma(sum(delta)) - n %*% digamma(delta[i]) + T[i]}#Col:taxa,Row:sample
  gradient<- as.matrix(gradient_list)
  test_stat <- t(gradient) %*%  (solve(Fisher)/n) %*% gradient
  test_stat <- as.numeric(test_stat)
  p_val <- 1 - pchisq(test_stat, d)
  #RaoScore[1] = p_val
  #RaoScore[2] = test_stat
  #RaoScore
  return(p_val)
}

MLE_W <- function (x, weights = NULL, eps = 10^(-5), convcrit = 1e-05, 
                   maxit = 1000, oldfac = 0.3, progress = FALSE) #The mle estimations function
{
  N <- nrow(x)
  K <- ncol(x)
  
  x <- (x + eps)/(1 + 2 * eps)
  x <- x/rowSums(x)
  N <- nrow(x)
  if (is.null(weights)) {
    weights <- rep(1, N)
  }
  weights <- N * weights/sum(weights)
  log.pbar <- colMeans(weights * log(x))
  alphaprob <- colMeans(x * weights)
  p2 <- mean(x[, 1]^2 * weights)
  xsi <- (alphaprob[1] - p2)/(p2 - (alphaprob[1])^2)
  alpha <- xsi * alphaprob
  K1 <- matrix(1, K, K)
  conv <- 1
  iter <- 1
  while ((conv > convcrit) & (iter < maxit)) {    #This is the stepwise iteration that converges at the global max     
    alpha0 <- alpha
    g <- N * digamma(sum(alpha)) - N * digamma(alpha) + N * 
      log.pbar
    z <- N * sirt_digamma1(sum(alpha))
    H <- diag(-N * sirt_digamma1(alpha)) + z
    alpha <- alpha0 - solve(H, g)
    alpha[alpha < 0] <- 1e-10
    alpha <- alpha0 + oldfac * (alpha - alpha0)
    conv <- max(abs(alpha0 - alpha))
    if (progress) {
      print(paste(iter, sum(alpha), conv))
      utils::flush.console()
    }
    iter <- iter + 1
  }
  alpha0 <- sum(alpha)
  
  
  
  res <- list(alpha = alpha, alpha0 = alpha0)
  return(res)
}


#Function file to compute Fisher information
dirFIMPOW <- function(alpha){
  FIMPO <- matrix(nrow=length(alpha), ncol=length(alpha), -trigamma(sum(alpha)))
  diagonal <- c()
  for(i in 1:length(alpha)){
    diagonal <- c(diagonal, trigamma(alpha[i])-trigamma(sum(alpha)))
  }
  diag(FIMPO) <- diagonal
  return(FIMPO)
}

WaldTest = function(data,Delta)
  
{
  Delta_h<-MLE_W(data)$alpha
  Fisher<-dirFIMPOW(Delta_h)*nrow(data)
  
  WaldTest = numeric(1)
  names(WaldTest) = c("p-value")
  
  d = length(Delta) 
  
  #Asymvar=inv(Fisher)#Use library(matlib) to find inverse of matrix
  #Asymvar is the asymptotic covariance matrix is the inverse of Fisher information (Consistent estimator) divided by n. 
  
  W = ((t(Delta_h-Delta)) %*% (Fisher) %*%(Delta_h-Delta))
  W = as.numeric(W)
  pval = 1-pchisq(W,d)
  
  WaldTest[1] = pval
  WaldTest
  
} 


dirtestsim <- function(sim_param){
  del = sim_param[[1]]
  d = sim_param[[2]]
  alphas = sim_param[[3]]
  ns = sim_param[[4]]
  K = sim_param[[5]]
  
  sizes.A <- c()
  sizes.W <- c()
  sizes.R <- c()
  
  for(d0 in del){
    ress1.A <- list()
    ress2.A <- list()
    ress3.A <- list()
    ress1.W <- list()
    ress2.W <- list()
    ress3.W <- list()
    ress1.R <- list()
    ress2.R <- list()
    ress3.R <- list()
    delta <- rep(d0, d)
    for(n in ns){
      res.A <- c()
      res.W <- c()
      res.R <- c()
      for(k in 1:K){
        dirdata <- matrix(delta, nrow=n, ncol=length(delta), byrow=TRUE) %>%  #data matrix initialization
          sirt::dirichlet.simul()
        res.A <- c(res.A,dirimean.test(dirdata, delta)$info[2])
        res.W <- c(res.W,WaldTest(dirdata, Delta=delta))
        res.R <- c(res.R,RaoScore(X=dirdata, delta =delta, dirFIMPO(delta)))
      }
      
      ress1.A <- append(ress1.A, mean(res.A<=0.1))
      ress2.A <- append(ress2.A, mean(res.A<=0.05))
      ress3.A <- append(ress3.A, mean(res.A<=0.01))
      ress1.W <- append(ress1.W, mean(res.W<=0.1))
      ress2.W <- append(ress2.W, mean(res.W<=0.05))
      ress3.W <- append(ress3.W, mean(res.W<=0.01))
      ress1.R <- append(ress1.R, mean(res.R<=0.1))
      ress2.R <- append(ress2.R, mean(res.R<=0.05))
      ress3.R <- append(ress3.R, mean(res.R<=0.01))
    }
    
    sizes.A <- append(sizes.A,list(list(ress1.A, ress2.A, ress3.A)))
    sizes.W <- append(sizes.W,list(list(ress1.W, ress2.W, ress3.W)))
    sizes.R <- append(sizes.R,list(list(ress1.R, ress2.R, ress3.R)))
    
  }
  
  sizestest <- list(sizes.A, sizes.W, sizes.R)
}


dirplot <- function(file,
                    sizes, 
                    sim_param, 
                    type = c("method", "d"), 
                    delp = 0.5, 
                    method = "ALRT"){
  
  
  require(dplyr)
  require(magrittr)
  require(ggplot2)
  
  dels = sim_param[[1]]
  d = sim_param[[2]]
  alphas = sim_param[[3]]
  ns = sim_param[[4]]
  K = sim_param[[5]]
  
  if("method" %in% type){
    
    if(delp == 0.5){
      delc = 1
    }
    else if(delp == 1){
      delc = 2
    }
    else{
      delc = 3
    }
    
    pdf(paste0(file, "bymethod.pdf"), width=5, height = 3)
    
    for(i in 1:3){
      size <- unlist(c(sizes[[1]][[delc]][[i]], sizes[[2]][[delc]][[i]], sizes[[3]][[delc]][[i]]))
      delta <- c(rep("ALRT", length(ns)), rep("WT", length(ns)), rep("RT", length(ns)))
      nss <- c(ns, ns, ns)
      pd <- data.frame("ns"=nss, "delta" = as.factor(delta), "data" = size)
      er <- sqrt(alphas[i]*(1-alphas[i])/K)
      p <- ggplot(pd, aes(x=ns, y=data, color = delta)) +
        geom_point(position=position_dodge(5), size=0.8) +  
        geom_errorbar(aes(ymin=data-er, ymax=data+er, color=delta), width=.2,
                      position=position_dodge(5)) +
        scale_color_brewer(palette="Dark2") +
        geom_hline(yintercept = alphas[i]+0.2*alphas[i], color="red", linetype="dashed") +
        geom_hline(yintercept = alphas[i], color="black", linetype="dashed") +
        geom_hline(yintercept = alphas[i]-0.2*alphas[i], color="red", linetype="dashed") +
        xlab("n") +
        ylab("Size (SE)") +
        theme_bw()
      print(p)
    }
    dev.off()
  }
  
  if("d" %in% type){
    
    if(method == "ALRT"){
      delc = 1
    }
    else if(method == "WT"){
      delc = 2
    }
    else{
      delc = 3
    }
    
    
    pdf(paste0(file, "bydel.pdf"), width=5, height = 3)
    
    for(i in 1:3){
      size <- unlist(c(sizes[[delc]][[1]][[i]], sizes[[delc]][[2]][[i]], sizes[[delc]][[3]][[i]]))
      delta <- c(rep(0.5, length(ns)), rep(1, length(ns)), rep(5, length(ns)))
      nss <- c(ns, ns, ns)
      pd <- data.frame("ns"=nss, "delta" = as.factor(delta), "data" = size)
      er <- sqrt(alphas[i]*(1-alphas[i])/K)
      p <- ggplot(pd, aes(x=ns, y=data, color = delta)) +
        geom_point(position=position_dodge(5), size=0.8) +  
        geom_errorbar(aes(ymin=data-er, ymax=data+er, color=delta), width=.2,
                      position=position_dodge(5)) +
        scale_color_brewer(palette="Dark2") +
        geom_hline(yintercept = alphas[i]+0.2*alphas[i], color="red", linetype="dashed") +
        geom_hline(yintercept = alphas[i], color="black", linetype="dashed") +
        geom_hline(yintercept = alphas[i]-0.2*alphas[i], color="red", linetype="dashed") +
        xlab("n") +
        ylab("Size (SE)") +
        theme_bw()
      print(p)
    }
    #dev.off()
  }
  
  
}


DistanceFromSymmetry <- function(del0, del1){
  require(wavethresh)
  d <- length(del0)
  new0 <- c()
  new1 <- c()
  for(j in 1:d){
    new0 <- c(new0, del0[j]/(sum(del0)/d))
    new1 <- c(new1, del1[j]/(sum(del1)/d))
  }
  return(wavethresh::l2norm(new0, new1))
}
