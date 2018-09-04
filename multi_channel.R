data_multi<-(matrix(data =c(1,7,5,15,50,1,5,12,17,55,2,4,2,19,60,2,7,9,21,65,1,9,6,
                           23,70,1,3,5,25,75,2,10,10,27,80,3,6,12,29,85,2,9,7,31,90,1,13,13,33,95),
                     byrow = TRUE, nrow = 10, ncol = 5))

# Multi-channel profile log-likelihood:
LogLik_PM <- function(psi,data_multi){
  gamma_hat<-c()
  beta_hat<-c()
  logLik<-c()
  for (i in 1:10) {
    y1<-data_multi[i,1]
    y2<-data_multi[i,2]
    y3<-data_multi[i,3]
    t<-data_multi[i,4]
    u<-data_multi[i,5]
    
    A1=-(psi+u)^2+(1+t)*psi*(psi+u)
    B1=(psi+u)*(y1+2*y3+y2)-(1+t)*psi*(y1+y3)
    C1=-y3*(y1+y3+y2)
    
    if ((B1^2-4*A1*C1) < 0 ){
      gamma_hat[i] = (-B1)/(2*A1)
    } else {
      gamma_hat[i] = (-B1+sqrt(B1^2-4*A1*C1))/(2*A1)
    }
    if (A1 == 0 ){
      gamma_hat[i] = (-C1)/(B1)
    } 

    A2=(1+t)
    B2=(psi*gamma_hat[i]*(1+t)-(y1+y2))
    C2= -y2*psi*gamma_hat[i]
    
    if ((B2^2-4*A2*C2) < 0 ){
      beta_hat[i]  =  (-B2)/(2*A2)
    } else {
      beta_hat[i]  =  (-B2+(sqrt(B2^2-4*A2*C2)))/(2*A2)
    }
  #}
    if (y1 == 0 & y2 != 0 & y3 != 0){
      logLik[i] = y2*log(beta_hat[i]*t)+y3*log(gamma_hat[i]*u)-(psi*gamma_hat[i]+beta_hat[i]+beta_hat[i]*t+gamma_hat[i]*u)
    } else if (y1 != 0 & y2 == 0 & y3 != 0){
      logLik[i] = y1*log(psi*gamma_hat[i]+beta_hat[i])+y3*log(gamma_hat[i]*u)-(psi*gamma_hat[i]+beta_hat[i]+beta_hat[i]*t+gamma_hat[i]*u)
    } else if (y1 != 0 & y2 != 0 & y3 == 0){
      logLik[i] = y1*log(psi*gamma_hat[i]+beta_hat[i])+y2*log(beta_hat[i]*t)-(psi*gamma_hat[i]+beta_hat[i]+beta_hat[i]*t+gamma_hat[i]*u)
    } else if (y1 == 0 & y2 == 0 & y3 != 0){
      logLik[i] = y3*log(gamma_hat[i]*u)-(psi*gamma_hat[i]+beta_hat[i]+beta_hat[i]*t+gamma_hat[i]*u)
    } else if (y1 != 0 & y2 == 0 & y3 == 0){
      logLik[i] = y1*log(psi*gamma_hat[i]+beta_hat[i])-(psi*gamma_hat[i]+beta_hat[i]+beta_hat[i]*t+gamma_hat[i]*u)
    } else if (y1 == 0 & y2 != 0 & y3 == 0){
      logLik[i] = y2*log(beta_hat[i]*t)-(psi*gamma_hat[i]+beta_hat[i]+beta_hat[i]*t+gamma_hat[i]*u)
    }
    else if (y1 == 0 & y2 != 0 & y3 == 0){
      logLik[i] = -(psi*gamma_hat[i]+beta_hat[i]+beta_hat[i]*t+gamma_hat[i]*u)
    }
    else{
      logLik[i] = y1*log(psi*gamma_hat[i]+beta_hat[i])+y2*log(beta_hat[i]*t)+y3*log(gamma_hat[i]*u)-(psi*gamma_hat[i]+beta_hat[i]+beta_hat[i]*t+gamma_hat[i]*u)
    }
  }
  logLik<-sum(logLik)
  return(logLik)
}

#psi_hat MLE：# 11.48716
optimize(LogLik_PM, interval=c(-5, 50), maximum=TRUE,tol = 0.0000001,data_multi)$maximum 

loglik_m<-c()
for (i in 1:length(seq(-3,50,0.01))){
  loglik_m[i]<-LogLik_PM(seq(-3,50,0.01)[i],data_multi)
}

#########Multi-channel profile log-likelihood###############################################################
plot(seq(-3,50,0.01),loglik_m,type="l",main='Multi-channel profile log-likelihood')
abline(v=optimize(LogLik_PM, data_multi,interval=c(-5, 50), maximum=TRUE)$maximum, col="blue")



########################################################################
r_rootM<-function(psi){
  r<-c()
  for (i in 1:length(psi)){
    psi_hat<-optimize(LogLik_PM, interval=c(-5, 50),data_multi, maximum=TRUE,tol = 0.0000001)$maximum 
    log_psihat<-LogLik_PM(psi_hat,data_multi)
    logpsi<-LogLik_PM(psi[i],data_multi)
    r[i]<-sign(psi_hat-psi[i])*sqrt(2*(log_psihat-logpsi))
}
  return(r)
}


#p-value of test for psi0=0: r(0)
1-pnorm(r_rootM(0))  #=1- 0.9999997 #3.123492e-07

#significance function for N(0,1) approximation of multi-channel:
psi_<-seq(-3,50,0.001)
plot(psi_,pnorm(r_rootM(psi_)),type='l',xlab='psi',ylab='significance',main='Multi-channel significance function')
abline(h=0.99, col="blue")
abline(h=0.5, col="blue")
abline(h=0.01, col="blue")
abline(v=0, col="blue")




























