## bootstrap simulation steps of data D:
y1 = 1
y2 = 8
y3 = 14
t = 27
u = 80
#Initialize:
r1<-c()
#psi=0
rep=10000
distribution1<- function(psi,y1,y2,y3){
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
    
    psi_hat<-optimize(LogLik_P,c(-5,50),Y1,Y2,Y3,maximum=TRUE,tol = 0.0000001)$maximum 
    log_psihat<-LogLik_P(psi_hat,Y1,Y2,Y3)
    logpsi<-LogLik_P(psi,Y1,Y2,Y3)
    
    if ((log_psihat-logpsi)<0){
      r1[i]<- 0
    } else {
      r1[i]<-sign(psi_hat-psi)*sqrt(2*(log_psihat-logpsi))
    }
    
  }
return(r1)
}


r1<-distribution1(0,y1=1,y2=8,y3=14)   
which(is.na(r1)) 

plot(density(na.omit(r1)),main='Simulated distribution of r(0)',xlab='')

#pvalue of the simulation distribution:
1-sum(r1 < r_root(0,y1,y2,y3))/rep  

##########################################################################
#plot the bootstrap significance function:
s1<-c()
significance_b<-function(psi){
  for (i in 1:length(psi)){
    s1[i]<- sum(distribution1(psi[i],y1=1,y2=8,y3=14) < r_root(psi[i],y1,y2,y3))/10000
  }
  return(s1)
}

s1<-significance_b(seq(-3,50,0.5))

plot(seq(-3,50,0.5),s1,type='l',ylab='significance',xlab='psi',main='') #method1
lines(seq(-3,50,0.1),pnorm(r_root(seq(-3,50,0.1),y1,y2,y3)),type='l',lty=2,mian='') #unadjusted
lines(seq(-3,50,0.1),R_,type='l',lty=3,main='',lwd=1.3) #adjusted

###############confidence interval#############################################
psi1<-seq(-3,50,1)
L<-psi1[which.min(abs(s1 - 0.99))]
U<-psi1[which.min(abs(s1- 0.01))]

psi_l2<-seq(L-1,L+1,0.1)
L1<-psi_l2[which.min(abs(significance_b(psi_l2)- 0.99))]

psi_l3<-seq(L1-0.1,L1+0.1,0.01)
L2<-psi_l3[which.min(abs(significance_b(psi_l3)- 0.99))]

psi_l4<-seq(L2-0.01,L2+0.01,0.001)
L3<-psi_l4[which.min(abs(significance_b(psi_l4)- 0.99))]

U<-psi1[which.min(abs(s1- 0.01))]

psi_u2<-seq(U-1,U+1,0.1)
U1<-psi_u2[which.min(abs(significance_b(psi_u2)- 0.01))]

psi_u3<-seq(U1-0.1,U1+0.1,0.01)
U2<-psi_u3[which.min(abs(significance_b(psi_u3)- 0.01))]

psi_u4<-seq(U2-0.01,U2+0.01,0.001)
U3<-psi_u4[which.min(abs(significance_b(psi_u4)- 0.01))]


###############first coverage#############################################

r11<-c() # value of test statistic:
p11<-c()
rep=10000
#constrained MLE
psi = 0
beta = exp(1)^1.1
gamma = 1
t = 33
u = 100

y<-matrix(NA,rep,3)

for (i in 1:rep){
  # simulated distribution:
  y1 <- rpois(1,gamma*psi+beta)
  y2 <- rpois(1,beta*t)
  y3 <- rpois(1,gamma*u)
  y[i,]<-c(y1,y2,y3)
  
  r11<-distribution1(psi,y1,y2,y3)
  p11[i] <- sum(r11 < r_root(psi,y1,y2,y3))/10000 #significance functions
}
which(is.na(p11)) 




