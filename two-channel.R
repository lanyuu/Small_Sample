#####two channel case #######
data2<-(matrix(data =c(1,7,5,15,50,1,5,12,17,55),byrow = TRUE, nrow = 2, ncol = 5))
# profile log-likelihood
LogLik_P2 <- function(psi,data2){
  gamma_hat<-c()
  beta_hat<-c()
  logLik<-c()
  for (i in 1:2) {
    y1<-data2[i,1]
    y2<-data2[i,2]
    y3<-data2[i,3]
    t<-data2[i,4]
    u<-data2[i,5]
    
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
    } else if (y1 == 0 & y2 == 0 & y3 == 0){
      logLik[i] = -(psi*gamma_hat[i]+beta_hat[i]+beta_hat[i]*t+gamma_hat[i]*u)
    } else{
      logLik[i] = y1*log(psi*gamma_hat[i]+beta_hat[i])+y2*log(beta_hat[i]*t)+y3*log(gamma_hat[i]*u)-(psi*gamma_hat[i]+beta_hat[i]+beta_hat[i]*t+gamma_hat[i]*u)
    }
  }
  logLik<-sum(logLik)
  return(logLik)
}

########## psi_hat MLE ############
optimize(LogLik_P2, interval=c(-5, 50), maximum=TRUE,tol = 0.0000001,data2)$maximum #3.688137

######### 2-channel profile log-likelihood #############
loglik_2<-c()
for (i in 1:length(seq(-3,50,0.01))){
  loglik_2[i]<-LogLik_P2(seq(-3,50,0.01)[i],data2)
}

plot(seq(-3,50,0.01),loglik_2,type="l",main='2-channel profile log-likelihood')
abline(v=optimize(LogLik_P2, data2,interval=c(-5, 50), maximum=TRUE)$maximum, col="blue")

######## R(psi) ###########
r_root2<-function(psi,data2){
  r<-c()
  for (i in 1:length(psi)){
    psi_hat<-optimize(LogLik_P2, interval=c(-5, 50),data2, maximum=TRUE,tol = 0.0000001)$maximum 
    log_psihat<-LogLik_P2(psi_hat,data2)
    logpsi<-LogLik_P2(psi[i],data2)
    if((log_psihat-logpsi)<0){
      r[i]<-0
    } else{
      r[i]<-sign(psi_hat-psi[i])*sqrt(2*(log_psihat-logpsi))
    }
  }
  return(r)
}


###################################################################################################################
####### R*of 2-channel####
gamma_hat<-c()
beta_hat<-c()
gammahat<-c()
betahat<-c()

# Modified likelihood root test statistic R*:
R_star2<-function(psi,data2){
  ### constrained MLEs ####
  for (i in 1:2) {
    y1<-data2[i,1]
    y2<-data2[i,2]
    y3<-data2[i,3]
    t<-data2[i,4]
    u<-data2[i,5]
    
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
  }
  #### overall MLEs ##
  ooo<-optimize(LogLik_P2, interval=c(-5, 50), maximum=TRUE,tol = 0.0000001,data2)
  psihat<-ooo$maximum
  for (i in 1:2) {
    y1<-data2[i,1]
    y2<-data2[i,2]
    y3<-data2[i,3]
    t<-data2[i,4]
    u<-data2[i,5]
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
    B2=(psihat*gammahat*(1+t)-(y1+y2))
    C2= -y2*psihat*gammahat
    
    if ((B2^2-4*A2*C2) < 0 ){
      betahat[i]  =  (-B2)/(2*A2)
    } else {
      betahat[i]  =  (-B2+(sqrt(B2^2-4*A2*C2)))/(2*A2)
    }
  }
  ############################################
  t<-data2[1:2,4]
  u<-data2[1:2,5]
  y1<-data2[1:2,1]
  y2<-data2[1:2,2]
  y3<-data2[1:2,3]
  ###########################################
  U11<- sum(gammahat*log(gammahat*psihat+betahat))-sum(gammahat*log(gamma_hat*psi+beta_hat))
  U21<- log(gammahat[1]*psihat+betahat[1])+t[1]*log(betahat[1]*t[1])-(log(gamma_hat[1]*psi+beta_hat[1])+t[1]*log(beta_hat[1]*t[1]))
  U31<- psihat*log(gammahat[1]*psihat+betahat[1])+u[1]*log(gammahat[1]*u[1])-( psihat*log(gamma_hat[1]*psi+beta_hat[1])+u[1]*log(gamma_hat[1]*u[1]))
  U41<- log(gammahat[2]*psihat+betahat[2])+t[2]*log(betahat[2]*t[2])-(log(gamma_hat[2]*psi+beta_hat[2])+t[2]*log(beta_hat[2]*t[2]))
  U51<- psihat*log(gammahat[2]*psihat+betahat[2])+u[2]*log(gammahat[2]*u[2])-(psihat*log(gamma_hat[2]*psi+beta_hat[2])+u[2]*log(gamma_hat[2]*u[2]))
  U12<- gammahat[1]/(gamma_hat[1]*psi+beta_hat[1])
  U22<- 1/(gamma_hat[1]*psi+beta_hat[1])+t[1]/beta_hat[1]
  U32<- psihat/(gamma_hat[1]*psi+beta_hat[1])
  U42=U52=0
  U13<- (psi*gammahat[1])/(gamma_hat[1]*psi+beta_hat[1])
  U23<-psi/(gamma_hat[1]*psi+beta_hat[1])
  U33<-(psi*psihat)/(gamma_hat[1]*psi+beta_hat[1])+u[1]/gamma_hat[1]
  U43=U53=0
  U14<- gammahat[2]/(gamma_hat[2]*psi+beta_hat[2])
  U24=U34=0
  U44<- 1/(gamma_hat[2]*psi+beta_hat[2])+t[2]/beta_hat[2]
  U54<- psihat/(gamma_hat[2]*psi+beta_hat[2])
  U15<- (psi*gammahat[2])/(gamma_hat[2]*psi+beta_hat[2])
  U25=U35=0
  U45<- psi/(gamma_hat[2]*psi+beta_hat[2])
  U55<-(psi*psihat)/(gamma_hat[2]*psi+beta_hat[2])+u[2]/gamma_hat[2]
  
  U<-matrix(c(U11,U12,U13,U14,U15, U21,U22,U23,U24,U25, 
              U31,U32,U33,U34,U35, U41,U42,U43,U44,U45, U51,U52,U53,U54,U55), nrow=5, byrow=TRUE)
  ###########################################
  L11<-sum((gammahat*gammahat)/(gammahat*psihat+betahat))
  L12<-gammahat[1]/(gammahat[1]*psihat+betahat[1])
  L13<-(psihat*gammahat[1])/(gammahat[1]*psihat+betahat[1])
  L14<-(gammahat[2])/(gammahat[2]*psihat+betahat[2])
  L15<-(psihat*gammahat[2])/(gammahat[2]*psihat+betahat[2])
  L21<-L12
  L22<-1/(gammahat[1]*psihat+betahat[1])+t[1]/betahat[1]
  L23<-psihat/(gammahat[1]*psihat+betahat[1])
  L24=L25=L34=L35=L42=L43=L52=L53=0
  L31=L13
  L32=L23
  L33<-(psihat*psihat)/(gammahat[1]*psihat+betahat[1])+u[1]/gammahat[1]
  L41=L14
  L44<-1/(gammahat[2]*psihat+betahat[2])+t[2]/betahat[2]
  L45<-psihat/(gammahat[2]*psihat+betahat[2])
  L51=L15
  L54=L45
  L55<-(psihat*psihat)/(gammahat[2]*psihat+betahat[2])+u[2]/gammahat[2]
  
  L<-matrix(c(L11,L12,L13,L14,L15, L21,L22,L23,L24,L25, 
              L31,L32,L33,L34,L35, L41,L42,L43,L44,L45, L51,L52,L53,L54,L55), nrow=5, byrow=TRUE)
  
  ###########################################
  J11<-sum((y1*gammahat*gammahat)/(gammahat*psihat+betahat))
  J12<-(y1[1]*gammahat[1])/(gammahat[1]*psihat+betahat[1])^2
  J13<-1-(y1[1]*betahat[1])/(gammahat[1]*psihat+betahat[1])^2
  J14<- (y1[2]*gammahat[2])/(gammahat[2]*psihat+betahat[2])^2
  J15<- 1-(y1[2]*betahat[2])/(gammahat[2]*psihat+betahat[2])^2
  J21<-J12
  J22<-(y1[1])/(gammahat[1]*psihat+betahat[1])^2+y2[1]/(betahat[1])^2
  J23<- (y1[1]*psihat)/(gammahat[1]*psihat+betahat[1])^2
  J24=J25=J34=J35=J42=J43=J52=J53=0
  J31=J13
  J32=J23
  J33<- (y1[1]*psihat*psihat)/(gammahat[1]*psihat+betahat[1])^2+y3[1]/(gammahat[1])^2
  J41=J14
  J44<- (y1[2])/(gammahat[2]*psihat+betahat[2])^2+y2[2]/(betahat[2])^2
  J45<- (y1[2]*psihat)/(gammahat[2]*psihat+betahat[2])^2
  J51=J15
  J54=J45
  J55<- (y1[2]*psihat*psihat)/(gammahat[2]*psihat+betahat[2])^2+y3[2]/(gammahat[2])^2
  
  J<-matrix(c(J11,J12,J13,J14,J15, J21,J22,J23,J24,J25, 
              J31,J32,J33,J34,J35, J41,J42,J43,J44,J45, J51,J52,J53,J54,J55), nrow=5, byrow=TRUE)
  
  ###########################################
  
  j11<-(y1[1])/(gamma_hat[1]*psi+beta_hat[1])^2+y2[1]/(beta_hat[1])^2
  j12<-(y1[1]*psi)/(gamma_hat[1]*psi+beta_hat[1])^2
  j13=j14=j23=j24=j31=j32=j41=j42=0
  j21<-j12
  j22<-(y1[1]*psi*psi)/(gamma_hat[1]*psi+beta_hat[1])^2+y3[1]/(gamma_hat[1])^2
  j33<-(y1[2])/(gamma_hat[2]*psi+beta_hat[2])^2+y2[2]/(beta_hat[2])^2
  j34<-(y1[2]*psi)/(gamma_hat[2]*psi+beta_hat[2])^2
  j43=j34
  j44<-(y1[2]*psi*psi)/(gamma_hat[2]*psi+beta_hat[2])^2+y3[2]/(gamma_hat[2])^2
  
  j<-matrix(c(j11,j12,j13,j14, j21,j22,j23,j24,
              j31,j32,j33,j34, j41,j42,j43,j44), nrow=4, byrow=TRUE)
  
  ############################################
  v<-(det(U)/det(L))*sqrt(det(J)/det(j))
  R<-r_root2(psi,data2)+(1/r_root2(psi,data2))*log(v/r_root2(psi,data2))
  return(R)
}
  
######### pvalue ############  
1-pnorm(R_star2(0,data2)) # 0.1067988

######## CI ########
psi_<-seq(-3,50,0.001)
x<-psi_
R_<-c()
for (i in 1:length(psi_)){
  R_[i]<-pnorm(R_star2(psi_[i],data2))
  U_R<-x[which.min(abs(R_-0.99))] 
  L_R<-x[which.min(abs(R_-0.01))]
}
#R_
U_R #-2.219
L_R #24.546
x[which.min(abs(R_-0.5))] #3.735
#L_R-U_R= 26.765

#significance functions:
plot(psi_,R_,type='l',lty=1,main='')
pnorm(R_star2(3.735,data2))

###############################################




















