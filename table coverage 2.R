########### r coverage ###########
psi = 0
beta = exp(1)^1.1
gamma = 1
t = 33
u = 100


# Tabel 1 for empirical converage for normal approximation:
x<-seq(-3,50,0.01)
j=0
#bound<-matrix(data=NA,10000,2)
for (i in 1:10000){
  y1 <- rpois(1,gamma*psi+beta)
  y2 <- rpois(1,beta*t)
  y3 <- rpois(1,gamma*u)
  psi_hat<-optimize(LogLik_P, interval=c(-5, 50), maximum=TRUE,y1,y2,y3,tol=0.000001)$maximum 
  log_psihat<-LogLik_P(psi_hat,y1,y2,y3)
  
  r<-sign(psi_hat-x)*sqrt(2*(log_psihat-LogLik_P(x,y1,y2,y3)))
  upper<- x[which.min(abs(pnorm(r)-0.99))] ##lower<- x[which.min(abs(pnorm(r)-0.99))]
  #bound[i,]<-c(lower,upper)
  if(psi<=upper){
    j=j+1
  }
}
j/10000 # for prob=0.01 check 0.99

####### R*coverage one channel ############
psi = 1
beta = exp(1)^1.1
gamma = 1
t = 33
u = 100

R_<-c()
rep=10000

for (i in 1:rep){
  y1 <- rpois(1,gamma*psi+beta)
  y2 <- rpois(1,beta*t)
  y3 <- rpois(1,gamma*u)
  if(y1==0){
      R_[i]<-r_root(psi,y1,y2,y3)
    }else {
    R_[i]<-(R_star(psi,y1,y2,y3))
    }
}

pstar<-c()
psi = 1
beta = exp(1)^1.1
gamma = 1
t = 33
u = 100


# Tabel 1 for empirical converage for normal approximation:
x<-seq(-3,50,0.01)
j=0
for (i in 1:10000){
  y1 <- rpois(1,gamma*psi+beta)
  y2 <- rpois(1,beta*t)
  y3 <- rpois(1,gamma*u)
  r<-r_root(psi)
  upper<- x[which.min(abs(pnorm(r)-0.01))]
  
  if( 1>upper){
    j=j+1
  }
}
########

psi = 0
beta = exp(1)^1.1
gamma = 1
t = 33
u = 100

R_<-c()
rep=10000

for (i in 1:rep){
  y1 <- rpois(1,gamma*psi+beta)
  y2 <- rpois(1,beta*t)
  y3 <- rpois(1,gamma*u)
  if(y1==0){
    R_[i]<-r_root(psi,y1,y2,y3)
  }else {
    R_[i]<-(R_star(psi,y1,y2,y3))
  }
}

