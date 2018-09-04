data2<-(matrix(data =c(1,7,5,15,50,1,5,12,17,55),byrow = TRUE, nrow = 2, ncol = 5))

#Initialize:
r3<-c()
#psi=0

gamma_hat<-c()
beta_hat<-c()
D<-matrix(data=NA,byrow = TRUE, nrow = 2, ncol = 5)
t = c(15,17)
u = c(50,55)
D[,4]<-t
D[,5]<-u
gg<-c()

rep1=100
check<-c()
distribution3_2<- function(psi,data2){
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
      
      A2=(1+t)
      B2=(psihat*gammahat*(1+t)-(y1+y2))
      C2= -y2*psihat*gammahat
      
      if ((B2^2-4*A2*C2) < 0 ){
        betahat[i]  =  (-B2)/(2*A2)
      } else {
        betahat[i]  =  (-B2+(sqrt(B2^2-4*A2*C2)))/(2*A2)
      }
      gg[i]<-gammahat[i]*psihat+betahat[i]
      if (gg[i]<0){
        gg[i]<-0
      }
      
      Y1 <- rpois(1,gg[i])
      Y2 <- rpois(1,betahat[i]*t)
      Y3 <- rpois(1,gammahat[i]*u)
      c(y1,y2,y3)
      if(Y1==0&Y2==0&Y3==0){
        D[i,1:3]<-c(rpois(1,gg[i])+1,rpois(1,betahat[i]*t),rpois(1,gammahat[i]*u))
      } else{
        D[i,1:3]<-c(Y1,Y2,Y3)
      }
    }
    
    psi_hat<-optimize(LogLik_P2, interval=c(-5, 50), D,maximum=TRUE,tol = 0.000000001)$maximum 
    log_psihat<-LogLik_P2(psi_hat,D)
    logpsi<-LogLik_P2(psihat,D)
    
    if ((log_psihat-logpsi) < 0  ){
      r3[ii]<- 0
    } else {
      r3[ii]<-sign(psi_hat-psi)*sqrt(2*(log_psihat-logpsi))
    }
  }
  return(r3)
}
r3<-distribution3_2(0,data2)
which(is.na(r3)) 


###################################
#plot the significance function:
s3<-c()
significance_23<-function(psi){
  for (i in 1:length(psi)){
    s3[i]<- sum(distribution3_2(psi[i],data2) < r_root2(psi[i],data2))/rep1
  }
  return(s3)
}

s3<-significance_23(seq(-3,50,1))
############### second coverage#############################################

r33<-c() # value of test statistic:
p33<-c()

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
  R_<-distribution3_2(psi,data2)
  p33[i] <- sum(R_< r_root2(psi,data2))/rep1
}









