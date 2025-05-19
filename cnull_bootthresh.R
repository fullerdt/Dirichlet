library(magrittr)
library(dplyr)


M=10000
K=500
d=2
del <- .1
alphas <- c(0.05)
ns <- c(10)
C=50
delta1 <- rep(4, d)
bounds <- c(1,3)


bootsall <- list()
for(d0 in del){
  delta <- rep(d0, d)
  boots <- c()
  for(n in ns){
    resQuant <- list()
    bootsperk <- list()
    for(k in 1:K){
      dirdata <- matrix(delta, nrow=n, ncol=length(delta), byrow=TRUE) %>%  
        sirt::dirichlet.simul()
      boot <- numeric(M)
      for (m in 1:M){
        z <- sum(log(dirdata))
        del_max<- optimize(fun_del, c(0,1000), tol = 0.0001,maximum = TRUE, n1=n , d1=d, z1=z)
        bootsample<- matrix(del_max$objective, nrow=n, ncol=length(delta), byrow=TRUE) %>%  
          sirt::dirichlet.simul()
        boot[m] <- newdirimean.test(bootsample)$info[1]
      }
      bootCutoff<-quantile(sort(boot),(1-alphas[1]))
      boots <- c(boots, bootCutoff)
    }
    bootsperk <- append(bootsperk, boots)
  }
  bootsall <- append(bootsall, list(list(bootsperk))) 
}   

save(bootsall,file="FINAL_Size_cutoffval_cnull_short_K500M10kns10del.1_d2.Rdata")



