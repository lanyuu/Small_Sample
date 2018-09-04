##################################################
########## R* ########################################
##################################################
y1 = 1
y2 = 8
y3 = 14
t = 27
u = 80

# Modified likelihood root test statistic R*:
R_star<-function(psi,y1,y2,y3){
  #constrained MLEs
  A1=-(psi+u)^2+(1+t)*psi*(psi+u)
  B1=(psi+u)*(y1+2*y3+y2)-(1+t)*psi*(y1+y3)
  C1=-y3*(y1+y3+y2)
  gamma_hat = (-B1+sqrt(B1^2-4*A1*C1))/(2*A1) 
  
  A2=(1+t)
  B2=(psi*gamma_hat*(1+t)-(y1+y2))
  C2= -y2*psi*gamma_hat
  beta_hat  =  (-B2+sqrt(B2^2-4*A2*C2))/(2*A2)
  
  #overall MLEs
  ooo<-optimize(LogLik_P, interval=c(-5, 50), maximum=TRUE,y1,y2,y3,tol=0.0000001)
  psihat<-ooo$maximum
  A11=-(psihat+u)^2+(1+t)*psihat*(psihat+u)
  B11=(psihat+u)*(y1+2*y3+y2)-(1+t)*psihat*(y1+y3)
  C11=-y3*(y1+y3+y2)
  

  gammahat = (-B11+sqrt(B11^2-4*A11*C11))/(2*A11)
  
  
  A22=(1+t)
  B22=(psihat*gammahat*(1+t)-(y1+y2))
  C22= -y2*psihat*gammahat
  
  if (y2==0){
    betahat = -B22/(A22)
  }else {
    betahat = (-B22+(sqrt(B22^2-4*A22*C22)))/(2*A22)
  }
  ############################################
  a1<-gammahat*log((gammahat*psihat+betahat)/(gamma_hat*psi+beta_hat))
  a2<-log((gammahat*psihat+betahat)/(gamma_hat*psi+beta_hat))+t*log(betahat/beta_hat)
  a3<-psihat*log((gammahat*psihat+betahat)/(gamma_hat*psi+beta_hat))+u*log(gammahat/gamma_hat)
  b11<-(gammahat)/(gamma_hat*psi+beta_hat)
  b12<-1/(gamma_hat*psi+beta_hat)+t/(beta_hat)
  b13<-(psihat)/(gamma_hat*psi+beta_hat)
  b21<-(gammahat*psi)/(gamma_hat*psi+beta_hat)
  b22<-(psi)/(gamma_hat*psi+beta_hat)
  b23<- (psihat*psi)/(gamma_hat*psi+beta_hat)+(u)/(gamma_hat) 
  L<- matrix(c(a1, a2, a3,
               b11, b12, b13,
               b21, b22, b23), nrow=3, byrow=TRUE)
  ############################################
  J11<- (gammahat*gammahat)/(gammahat*psihat+betahat)
  J12<- (gammahat)/(gammahat*psihat+betahat)
  J13<- (1-(betahat)/(gammahat*psihat+betahat))
  J21<- J12
  J22<- 1/(gammahat*psihat+betahat)+t/(betahat)
  J23<- psihat/(gammahat*psihat+betahat)
  J31<- J13
  J32<- J23
  J33<- (psihat*psihat)/(gammahat*psihat+betahat)+u/(gammahat)
  J <- matrix(c(J11,J12,J13,
                J21,J22,J23, 
                J31,J32,J33), nrow=3, byrow=TRUE)
  ############################################
  j22<-(gammahat*psihat+betahat)/((gamma_hat*psi+beta_hat)^2)+ (betahat*t)/((beta_hat)^2)
  j23<-((gammahat*psihat+betahat)*psi)/((gamma_hat*psi+beta_hat)^2)
  j32<- j23
  j33<-((gammahat*psihat+betahat)*psi*psi)/((gamma_hat*psi+beta_hat)^2)+(gammahat*u)/(gamma_hat*gamma_hat)
  jxx <- matrix(c(j22,j23, 
                  j32,j33), nrow=2, byrow=TRUE)
  
  ############################################
  v<-(det(L))/(sqrt(det(J)*det(jxx)))
  R<-r_root(psi,y1,y2,y3)+(1/r_root(psi,y1,y2,y3))*log(v/r_root(psi,y1,y2,y3))
  return(R)
}
1-pnorm(R_star(0,y1,y2,y3)) 

psi_<-seq(-3,50,0.01)
x<-psi_
R_<-c()
for (i in 1:length(psi_)){
  R_[i]<-pnorm(R_star(psi_[i],y1,y2,y3))
  U_R<-x[which.min(abs(R_-0.99))] 
  L_R<-x[which.min(abs(R_-0.01))]
}


#####first simulation################################
r_root<-function(psi,y1,y2,y3){
  psi_hat<-optimize(LogLik_P, interval=c(-5, 50), maximum=TRUE,y1,y2,y3,tol = 0.0000001)$maximum 
  log_psihat<-LogLik_P(psi_hat,y1,y2,y3)
  logpsi<-LogLik_P(psi,y1,y2,y3)
  if((log_psihat-logpsi)< 0 ){
    r<-0
  }else{
    r<-sign(psi_hat-psi)*sqrt(2*(log_psihat-logpsi))
  }
  return(r)
}
r_root(0,y1,y2,y3)
####################
R1<-c()
rep=100000
Y<-matrix(NA,rep,3)

distribution1_R <- function(psi,y1,y2,y3){
  for (i in 1:rep){
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
    
    
    Y1 <- rpois(1,gamma_hat*psi+beta_hat)
    Y2 <- rpois(1,beta_hat*t)
    Y3 <- rpois(1,gamma_hat*u)
    #c(Y1,Y2,Y3)
    
    Y[i,]<-c(Y1,Y2,Y3)

    if((Y1==0 & Y2==0 & Y3==0)){
      R1[i]<-0
    }else if(Y1==0 & Y2 !=0  & Y3 !=0){
      R1[i]<-r_root(psi,Y1,Y2,Y3)
    }else if(Y1==0 & Y2 ==0  & Y3 !=0){
      R1[i]<-r_root(psi,Y1,Y2,Y3)
    } else if (is.na(Y1)==TRUE & is.na(Y2)==TRUE & is.na(Y2)==TRUE){
      R1[i]<-0
    }else{
      R1[i]<-R_star(psi,Y1,Y2,Y3)
    }
  which(is.na(R1)) 
  }
  return(R1)
}

R1<-distribution1_R(psi=0,y1,y2,y3)

#plot the bootstrap significance function:
S1<-c()
significance_bR<-function(psi){
  for (i in 1:length(psi)){
    S1[i]<- sum(distribution1_R(psi[i],y1=1,y2=8,y3=14) < r_root(psi[i],y1,y2,y3))/rep
  }
  return(S1)
}

S1<-significance_bR(seq(-3,50,1))



