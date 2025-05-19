library(dplyr)
library(magrittr)
source("newutil.R")

K=20000

storagenew <- list()

for(del in c(0.5, 1, 2, 3, 5, 10)){
  for(d in c(2)){
    
    delta=c(rep(del,d))
    
    N=300
    

    t300 <- dirtest(K=K,N=N,d=d, delta=delta)

    
    N=250
    

    t250 <- dirtest(K=K,N=N,d=d, delta=delta)

    
    N=200
    

    t200 <- dirtest(K=K,N=N,d=d, delta=delta)

    
    N=150
    

    t150 <- dirtest(K=K,N=N,d=d, delta=delta)

    
    N=100
    

    t100 <- dirtest(K=K,N=N,d=d, delta=delta)

    
    N=50
    
    

    t50 <- dirtest(K=K,N=N,d=d, delta=delta)

    
    N=45
    

    t45 <- dirtest(K=K,N=N,d=d, delta=delta)

    
    N=40
    

    t40 <- dirtest(K=K,N=N,d=d, delta=delta)

    
    N=35
    

    t35 <- dirtest(K=K,N=N,d=d, delta=delta)

    
    N=30
    

    t30 <- dirtest(K=K,N=N,d=d, delta=delta)

    
    N=25
    

    t25 <- dirtest(K=K,N=N,d=d, delta=delta)

    
    N=20
    

    t20 <- dirtest(K=K,N=N,d=d, delta=delta)

    
    N=15
    

    t15 <- dirtest(K=K,N=N,d=d, delta=delta)

    
    N=10
    

    t10 <- dirtest(K=K,N=N,d=d, delta=delta)

    
    N=5
    

    t5 <- dirtest(K=K,N=N,d=d, delta=delta)

    
    ts<-list(t5, t10, t15, t20, t25, t30, t35, t40, t45, t50, t100, t150, t200, t250, t300)
    
    storagenew <- append(storagenew, list(ts))
  }
}









