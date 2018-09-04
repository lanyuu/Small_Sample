# Monro Robbins for the second simulation
psi_ <- seq(-20,80,0.001)
alpha <- 0.01
x <- psi_
####Calculate the CI based on second method########
psi_j<-c()
diff<-c()
CI2 <- function(alpha){
  psi_j[1] <- x[which.min(abs(pnorm(r_root(psi_,y1,y2,y3))-(1-alpha)))]
  c <- 2/pnorm(r_root(psi_j[1],y1,y2,y3))
  for (i in 1:100){
    ooo<-optimize(LogLik_P, interval=c(-5, 50), maximum=TRUE,y1,y2,y3)
    psihat<-ooo$maximum
    A1=-(psihat+u)^2+(1+t)*psihat*(psihat+u)
    B1=(psihat+u)*(y1+2*y3+y2)-(1+t)*psihat*(y1+y3)
    C1=-y3*(y1+y3+y2)
    
    if ((B1^2-4*A1*C1) < 0 ){
      gammahat = (-B1)/(2*A1)
    } else {
      gammahat = (-B1+sqrt(B1^2-4*A1*C1))/(2*A1)
    }
    
    A2=(1+t)
    B2=(psihat*gammahat*(1+t)-(y1+y2))
    C2= -y2*psihat*gammahat
    
    if ((B2^2-4*A2*C2) < 0 ){
      betahat  =  (-B2)/(2*A2)
    } else {
      betahat  =  (-B2+(sqrt(B2^2-4*A2*C2)))/(2*A2)
    }
    gg<-gammahat*psi_j[i]+betahat
    
    if (gg<0){
      gg<-0
    }
    
    Y1 <- rpois(1,gg)
    Y2 <- rpois(1,betahat*t)
    Y3 <- rpois(1,gammahat*u)
    
    psi_hat<-optimize(LogLik_P,c(-5,50),Y1,Y2,Y3,maximum=TRUE,tol = 0.0000001)$maximum 
    log_psihat<-LogLik_P(psi_hat,Y1,Y2,Y3)
    logpsi<-LogLik_P(psi_j[i],Y1,Y2,Y3)
    
    if ((log_psihat-logpsi)<0){
      r1<- 0
    } else {
      r1<-sign(psi_hat-psi_j[i])*sqrt(2*(log_psihat-logpsi))
    }
    
    if ( r1 <= r_root(psi_j[i],y1,y2,y3)  ){
      psi_j[i+1] <- psi_j[i] + (c*alpha)/i
    } else {
      psi_j[i+1] <- psi_j[i] - (c*(1-alpha))/i
    }
    diff[i] = abs(psi_j[i+1] - psi_j[i])
  }
  
  return(list(psi_j,diff))
  #return(psi_j)
}

result2<-CI2(0.01)

###############Average the CI with 10 replications based on second method################################
rep=10
CI_mat2 <- matrix(data = NA, nrow = rep, ncol = 2)
Width_CI2<- function(alpha){
  for (i in 1:rep){
    result_L<-CI2(alpha)
    index_L<-which(abs(result_L[[2]] - 0.0001) == sort(abs(result_L[[2]] - 0.0001), decreasing = FALSE)[1])
    LL<-result_L[[1]][index_L]
    
    result_U<-CI2(1-alpha)
    index_U<-which(abs(result_U[[2]] - 0.0001) == sort(abs(result_U[[2]] - 0.0001), decreasing = FALSE)[1])
    UL<-result_U[[1]][index_U]
    CI_mat2[i,]<-c(LL,UL)
  }
  return(CI_mat2)
}


limits2<-Width_CI2(0.01)
CI_average2<-colMeans(limits2)


###########Check second simulation gives wider CIs###############################################
CI_plot2<-c()
CI_average2<-matrix(NA,50,2)
for (i in 1:50){
  CI_average2[i,]<-colMeans(Width_CI2(0.01))
  CI_plot2[i]<-CI_average2[i,2]-CI_average2[i,1]
}




