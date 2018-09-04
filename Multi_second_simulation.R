data2<-(matrix(data =c(1,7,5,15,50,1,5,12,17,55),byrow = TRUE, nrow = 2, ncol = 5))

#Initialize:
r2<-c()
#psi=0

gamma_hat<-c()
beta_hat<-c()
D<-matrix(data=NA,byrow = TRUE, nrow = 2, ncol = 5)
t = c(15,17)
u = c(50,55)
D[,4]<-t
D[,5]<-u

rep1=100 
check<-c()
distribution2_2<- function(psi,data2){
  for (ii in 1:rep1){
    for (i in 1:2) {
      y1<-data2[i,1]
      y2<-data2[i,2]
      y3<-data2[i,3]
      t<-data2[i,4]
      u<-data2[i,5]
      
      ooo<-optimize(LogLik_P2, interval=c(-5, 50), maximum=TRUE,tol = 0.0000001,data2)
      psihat<-ooo$maximum
      A1=-(psihat+u)^2+(1+t)*psihat*(psihat+u)
      B1=(psihat+u)*(y1+2*y3+y2)-(1+t)*psihat*(y1+y3)
      C1=-y3*(y1+y3+y2)
      
      if ((B1^2-4*A1*C1) < 0 ){
        gammahat[i] = (-B1)/(2*A1)
      } else {
        gammahat[i] = (-B1+sqrt(B1^2-4*A1*C1))/(2*A1)
      }
      if (A1 == 0 ){
        gammahat[i] = (-C1)/(B1)
      } 
    
      A2=(1+t)
      B2=(psihat*gammahat[i]*(1+t)-(y1+y2))
      C2= -y2*psihat*gammahat[i]
      
      if ((B2^2-4*A2*C2) < 0 ){
        betahat[i]  =  (-B2)/(2*A2)
      } else {
        betahat[i]  =  (-B2+(sqrt(B2^2-4*A2*C2)))/(2*A2)
      }

      y1 <- rpois(1,gammahat[i]*psi+betahat[i])
      y2 <- rpois(1,betahat[i]*t)
      y3 <- rpois(1,gammahat[i]*u)
      c(y1,y2,y3)
      if(y1==0&y2==0&y3==0){
        D[i,1:3]<-c(rpois(1,gammahat[i]*psi+betahat[i])+1,rpois(1,betahat[i]*t),rpois(1,gammahat[i]*u))
      } else{
        D[i,1:3]<-c(y1,y2,y3)
      }
    }
    
    psi_hat<-optimize(LogLik_P2, interval=c(-5, 50), D,maximum=TRUE,tol = 0.000000001)$maximum 
    log_psihat<-LogLik_P2(psi_hat,D)
    logpsi<-LogLik_P2(psi,D)
    
    if ((log_psihat-logpsi) < 0  ){
      r2[ii]<- 0
    } else {
      r2[ii]<-sign(psi_hat-psi)*sqrt(2*(log_psihat-logpsi))
    }
  }
  return(r2)
}
r2<-distribution2_2(0,data2)
which(is.na(r2)) 

plot(density(r2),main='simulated distribution of r(0) of 2-channel')
#pvalue of the simulation distribution of 2-channel:
1-sum(r2 < r_root2(0,data2))/rep1 
###################################
#plot the significance function:
s2<-c()
significance_22<-function(psi){
  for (i in 1:length(psi)){
    s2[i]<- sum(distribution2_2(psi[i],data2) < r_root2(psi[i],data2))/rep1
  }
  return(s2)
}

s2<-significance_22(seq(-min(betahat/gammahat),50,1))



############### second coverage #############################################

r22<-c() # value of test statistic:
p22<-c()

#constrained MLE
psi = 0
beta = c(0.2,1.1)
gamma = c(0.2,0.65)
t = c(15,33)
u = c(50,95)

data2<-matrix(data=NA,2,5)
data2[1:2,4]=t
data2[1:2,5]=u

R_<-c()
rep=100
#y<-matrix(data=0,nrow=rep,ncol=3)
for (i in 1:rep){
  for (ii in 1:2){
    y1[ii] <- rpois(1,gamma[ii]*psi+beta[ii])
    y2[ii] <- rpois(1,beta[ii]*t[ii])
    y3[ii] <- rpois(1,gamma[ii]*u[ii])
    data2[ii,1:3]<-c(y1[ii],y2[ii],y3[ii])
  }
  R_<-distribution2_2(psi,data2)
  p22[i] <- sum(R_< r_root2(psi,data2))/rep1
}







