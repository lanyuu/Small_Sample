y1 = 1
y2 = 8
y3 = 14
t = 27
u = 80
##########################################################################
#Second simulation simulate specifying (psi,beta,gamma)=(psi_0,betahat,gammahat):2 overall MLE
# simulated distribution of r(psi_0)
r2<-c()
rep=10000
distribution2<- function(psi,y1,y2,y3){
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
  
  for (i in 1:rep){
    Y1 <- rpois(1,gammahat*psi+betahat)
    Y2 <- rpois(1,betahat*t)
    Y3 <- rpois(1,gammahat*u)
    psi_hat<-optimize(LogLik_P,c(-5,50),Y1,Y2,Y3,maximum=TRUE,tol = 0.0000001)$maximum 
    log_psihat<-LogLik_P(psi_hat,Y1,Y2,Y3)
    logpsi<-LogLik_P(psi,Y1,Y2,Y3)
    if (log_psihat-logpsi < 0){
      r2[i]<-0
    } else{
      r2[i]<-sign(psi_hat-psi)*sqrt(2*(log_psihat-logpsi))
    }
    
  }
  return(r2)
}



##########################################################################
#compare constrained_MLE of beta/gamma vs overall MLE of beta/gamma:
plot(density(distribution1(0,y1,y2,y3)), col = "red",main="(i) vs (ii)")
lines(density(distribution2(0,y1,y2,y3)), col = "green")

#plot the second significance function:
s2<-c()
significance_2<-function(psi){
  for (i in 1:length(psi)){
    s2[i]<- sum(distribution2(psi[i],y1=1,y2=8,y3=14) < r_root(psi[i],y1,y2,y3))/10000
  }
  return(s2)
}

########### second coverage###################

r22<-c() # value of test statistic:
p22<-c()
rep=10000
#constrained MLE
psi = 0
beta = exp(1)^1.1
gamma = 1
t = 33
u = 100

#y<-matrix(NA,rep,3)

for (i in 1:rep){
  # simulated distribution:
  y1 <- rpois(1,gamma*psi+beta)
  y2 <- rpois(1,beta*t)
  y3 <- rpois(1,gamma*u)
  #y[i,]<-c(y1,y2,y3)
  
  r22<-distribution2(psi,y1,y2,y3)
  p22[i] <- sum(r22 < r_root(psi,y1,y2,y3))/10000
}
which(is.na(p22)) 








