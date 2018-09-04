y1 = 1
y2 = 8
y3 = 14
t = 27
u = 80
#psi=0

# Profile log-likelihood:
LogLik_P <- function(psi,y1,y2,y3){
  A1=-(psi+u)^2+(1+t)*psi*(psi+u)
  B1=(psi+u)*(y1+2*y3+y2)-(1+t)*psi*(y1+y3)
  C1=-y3*(y1+y3+y2)

  if ((B1^2-4*A1*C1) < 0 ){
    gamma_hat = (-B1)/(2*A1)
  } else {
    gamma_hat = (-B1+sqrt(B1^2-4*A1*C1))/(2*A1)
  }
  
  A2=(1+t)
  B2=(psi*gamma_hat*(1+t)-(y1+y2))
  C2= -y2*psi*gamma_hat
  
  if ((B2^2-4*A2*C2) < 0 ){
    beta_hat  =  (-B2)/(2*A2)
  } else {
    beta_hat  =  (-B2+(sqrt(B2^2-4*A2*C2)))/(2*A2)
  }
  
  if (y1 == 0 & y2 != 0 & y3 != 0){
    logLik = y2*log(beta_hat*t)+y3*log(gamma_hat*u)-(psi*gamma_hat+beta_hat+beta_hat*t+gamma_hat*u)
  } else if (y1 != 0 & y2 == 0 & y3 != 0){
    logLik = y1*log(psi*gamma_hat+beta_hat)+y3*log(gamma_hat*u)-(psi*gamma_hat+beta_hat+beta_hat*t+gamma_hat*u)
  } else if (y1 != 0 & y2 != 0 & y3 == 0){
    logLik = y1*log(psi*gamma_hat+beta_hat)+y2*log(beta_hat*t)-(psi*gamma_hat+beta_hat+beta_hat*t+gamma_hat*u)
  } else if (y1 == 0 & y2 == 0 & y3 != 0){
    logLik = y3*log(gamma_hat*u)-(psi*gamma_hat+beta_hat+beta_hat*t+gamma_hat*u)
  } else if (y1 != 0 & y2 == 0 & y3 == 0){
    logLik = y1*log(psi*gamma_hat+beta_hat)-(psi*gamma_hat+beta_hat+beta_hat*t+gamma_hat*u)
  } else if (y1 == 0 & y2 != 0 & y3 == 0){
    logLik = y2*log(beta_hat*t)-(psi*gamma_hat+beta_hat+beta_hat*t+gamma_hat*u)
  }
  else{
    logLik = y1*log(psi*gamma_hat+beta_hat)+y2*log(beta_hat*t)+y3*log(gamma_hat*u)-(psi*gamma_hat+beta_hat+beta_hat*t+gamma_hat*u)
  }
  return(logLik)
}

#Plot of profile log likelihood:
plot(seq(-3,50,0.01),LogLik_P(seq(-3,50,0.01),y1,y2,y3),type='l',xlab='psi',ylab='log-likelihood',main='profile log-likelihood')
abline(v=optimize(LogLik_P, interval=c(-5, 50), maximum=TRUE,y1,y2,y3)$maximum, col="blue")


# Signed root test statistic r:
r_root<-function(psi,y1,y2,y3){
  psi_hat<-optimize(LogLik_P, interval=c(-5, 50), maximum=TRUE,y1,y2,y3,tol = 0.0000001)$maximum 
  log_psihat<-LogLik_P(psi_hat,y1,y2,y3)
  logpsi<-LogLik_P(psi,y1,y2,y3)
  r<-sign(psi_hat-psi)*sqrt(2*(log_psihat-logpsi))
  return(r)
}

#p-value of test for psi0=0: r(0)
pnorm(r_root(0,y1,y2,y3))  #=1-0.8373245

#significance function for N(0,1) approximation:
psi_<-seq(-3,50,0.001)
plot(psi_,pnorm(r_root(psi_,y1,y2,y3)),type='l',xlab='psi',ylab='significance',main='significance function')
abline(h=0.99, col="blue")
abline(h=0.5, col="blue")
abline(h=0.01, col="blue")
abline(v=0, col="blue")

#Confidence limits for null hypothsis:
x<-psi_
x[which.min(abs(pnorm(r_root(psi_,y1,y2,y3))-0.99))] #confidence limits (psi=-2.644)
x[which.min(abs(pnorm(r_root(psi_,y1,y2,y3))-0.01))] # 33.835
x[which.min(abs(pnorm(r_root(psi_,y1,y2,y3))-0.5))]  # 4.021

33.835+2.644 = 36.479









