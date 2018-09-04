##########Monro_Robbins for first simulation####
psi_ <- seq(-20,80,0.001)
alpha <- 0.01
x <- psi_
####Calculate the CI based on first method#######
psi_j<-c()
diff<-c()
CI1 <- function(alpha){
  psi_j[1] <- x[which.min(abs(pnorm(r_root(psi_,y1,y2,y3))-(1-alpha)))]
  c <- 2/pnorm(r_root(psi_j[1],y1,y2,y3))
  for (i in 1:100){
    A1=-(psi_j[i]+u)^2+(1+t)*psi_j[i]*(psi_j[i]+u)
    B1=(psi_j[i]+u)*(y1+2*y3+y2)-(1+t)*psi_j[i]*(y1+y3)
    C1=-y3*(y1+y3+y2)
    
    if ((B1^2-4*A1*C1) < 0 ){
      gamma_hat = (-B1)/(2*A1)
    } else {
      gamma_hat = (-B1+sqrt(B1^2-4*A1*C1))/(2*A1)
    }
    
    A2=(1+t)
    B2=(psi_j[i]*gamma_hat*(1+t)-(y1+y2))
    C2= -y2*psi_j[i]*gamma_hat
    
    if ((B2^2-4*A2*C2) < 0 ){
      beta_hat  =  (-B2)/(2*A2)
    } else {
      beta_hat  =  (-B2+(sqrt(B2^2-4*A2*C2)))/(2*A2)
    }
    
    Y1 <- rpois(1,gamma_hat*psi_j[i]+beta_hat)
    Y2 <- rpois(1,beta_hat*t)
    Y3 <- rpois(1,gamma_hat*u)
    
    psi_hat<-optimize(LogLik_P,c(-5,50),Y1,Y2,Y3,maximum=TRUE,tol = 0.0000001)$maximum 
    log_psihat<-LogLik_P(psi_hat,Y1,Y2,Y3)
    logpsi<-LogLik_P(psi_j[i],Y1,Y2,Y3)
    
    if ((log_psihat-logpsi)<0){
      r1<- 0
    } else {
      r1<-sign(psi_hat-psi_j[i])*sqrt(2*(log_psihat-logpsi))
    }
    
    if ( r1 <= r_root(psi_j[i],y1,y2,y3)  ){
      psi_j[i+1] <- psi_j[i] + (c*alpha)/i
    } else {
      psi_j[i+1] <- psi_j[i] - (c*(1-alpha))/i
    }
    diff[i] = abs(psi_j[i+1] - psi_j[i])
  }
  
  return(list(psi_j,diff))
  #return(psi_j)
}

result1<-CI1(0.01)
###############Average the CI with 10 replications based on first method########
rep=10
CI_mat <- matrix(data = NA, nrow = rep, ncol = 2)
Width_CI1<- function(alpha){
  for (i in 1:rep){
    result_L<-CI1(alpha)
    index_L<-which(abs(result_L[[2]] - 0.0001) == sort(abs(result_L[[2]] - 0.0001), decreasing = FALSE)[1])
    LL<-result_L[[1]][index_L]
    
    result_U<-CI1(1-alpha)
    index_U<-which(abs(result_U[[2]] - 0.0001) == sort(abs(result_U[[2]] - 0.0001), decreasing = FALSE)[1])
    UL<-result_U[[1]][index_U]
    CI_mat[i,]<-c(LL,UL)
  }
  return(CI_mat)
}

limits1<-Width_CI1(0.01)
CI_average1<-colMeans(limits1)
CI_average1[2]-CI_average1[1]

###########Check fist simulation gives wider CIs######
CI_plot1<-c()
ite=50
CI_average1<-matrix(NA,ite,2)

for (i in 1:ite){
  CI_average1[i,]<-colMeans(Width_CI1(0.01))
  CI_plot1[i]<-CI_average1[i,2]-CI_average1[i,1]
}

