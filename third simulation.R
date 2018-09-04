######################################################################
#third simulation specifying (psi,beta,gamma)=(psihat,betahat,gammahat)
rep=10000
r3<-c()
distribution3<- function(y1,y2,y3){
  ooo<-optimize(LogLik_P, interval=c(-5, 50), maximum=TRUE,y1,y2,y3,tol = 0.0000001)
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
  gg<-gammahat*psihat+betahat
  
  if (gg<0){
    gg<-0
  }
  
  for (i in 1:rep){
    Y1 <- rpois(1,gg)
    Y2 <- rpois(1,betahat*t)
    Y3 <- rpois(1,gammahat*u)
    
    psi_hat<-optimize(LogLik_P,c(-5,50),Y1,Y2,Y3,maximum=TRUE,tol = 0.0000001)$maximum 
    log_psihat<-LogLik_P(psi_hat,Y1,Y2,Y3)
    logpsi<-LogLik_P(psihat,Y1,Y2,Y3)
    sign(psi_hat-psihat)*sqrt(2*(log_psihat-logpsi))
    if (log_psihat-logpsi< 0){
      r3[i]<-0
    }else{
      r3[i]<-sign(psi_hat-psihat)*sqrt(2*(log_psihat-logpsi))
    }
  }
  return(r3)
}
r3<-distribution3(y1,y2,y3) 
1-sum(r3 < r_root(0,y1,y2,y3))/rep



# Significance function:
s3<-c()
significance_3<-function(psi){
  for (i in 1:length(psi)){
    s3[i]<- sum(r3 < r_root(psi[i],y1,y2,y3))/10000
  }
  return(s3)
}
#s3<-significance_3(seq(-3,50,0.1)
################## third coverages##################################

r33<-c() # value of test statistic:
p33<-c()
rep=10000
#constrained MLE
psi = 0
beta = exp(1)^1.1
gamma = 1
t = 33
u = 100
#y=matrix(NA,rep,3)
for (i in 1:rep){
  # simulated distribution:
  y1 <- rpois(1,gamma*psi+beta)
  y2 <- rpois(1,beta*t)
  y3 <- rpois(1,gamma*u)
  #y[i,]<-c(y1,y2,y3)
  r33<-distribution3(y1,y2,y3)
  #sum(r33 < r_root(psi))/10000
  p33[i] <- sum(r33 < r_root(psi,y1,y2,y3))/rep
}
p33
