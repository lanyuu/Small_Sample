psi_ <- seq(-3,50,0.001)
alpha <- 0.01
x <- psi_

data2<-(matrix(data =c(1,7,5,15,50,1,5,12,17,55),byrow = TRUE, nrow = 2, ncol = 5))
y1<-data2[1:2,1]
y2<-data2[1:2,2]
y3<-data2[1:2,3]
t<-data2[1:2,4]
u<-data2[1:2,5]
gamma_hat<-c()
beta_hat<-c()
gammahat<-c()
betahat<-c()
D<-matrix(data=NA,byrow = TRUE, nrow = 2, ncol = 5)
D[,4]<-t
D[,5]<-u
####Calculate the CI based on first method#####
psi_j<-c()
diff<-c()
CI12 <- function(alpha){
  psi_j[1] <- x[which.min(abs(pnorm(r_root2(psi_,data2))-(1-alpha)))]
  c <- 2/pnorm(r_root2(psi_j[1],data2))
  for (i in 1:100){
    for (ii in 1:2) {
      y1<-data2[ii,1]
      y2<-data2[ii,2]
      y3<-data2[ii,3]
      t<-data2[ii,4]
      u<-data2[ii,5]
      
      A1=-(psi_j[i]+u)^2+(1+t)*psi_j[i]*(psi_j[i]+u)
      B1=(psi_j[i]+u)*(y1+2*y3+y2)-(1+t)*psi_j[i]*(y1+y3)
      C1=-y3*(y1+y3+y2)
      
      if ((B1^2-4*A1*C1) < 0 ){
        gamma_hat[ii] = (-B1)/(2*A1)
      } else {
        gamma_hat[ii] = (-B1+sqrt(B1^2-4*A1*C1))/(2*A1)
      }
      if (A1 == 0 ){
        gamma_hat[ii] = (-C1)/(B1)
      } 
      
      A2=(1+t)
      B2=(psi_j[i]*gamma_hat[ii]*(1+t)-(y1+y2))
      C2= -y2*psi_j[i]*gamma_hat[ii]
      
      if ((B2^2-4*A2*C2) < 0 ){
        beta_hat[ii]  =  (-B2)/(2*A2)
      } else {
        beta_hat[ii]  =  (-B2+(sqrt(B2^2-4*A2*C2)))/(2*A2)
      }
      
      Y1 <- rpois(1,gamma_hat[ii]*psi_j[i]+beta_hat[ii])
      Y2 <- rpois(1,beta_hat[ii]*t)
      Y3 <- rpois(1,gamma_hat[ii]*u)
      c(Y1,Y2,Y3)
      if(Y1==0&Y2==0&Y3==0){
        D[ii,1:3]<-c(rpois(1,gamma_hat[ii]*psi_j[i]+beta_hat[ii])+1,rpois(1,beta_hat[ii]*t),rpois(1,gamma_hat[ii]*u))
      } else{
        D[ii,1:3]<-c(Y1,Y2,Y3)
      }
    }
    
    psi_hat<-optimize(LogLik_P2, interval=c(-5, 50), D,maximum=TRUE,tol = 0.000000001)$maximum 
    log_psihat<-LogLik_P2(psi_hat,D)
    logpsi<-LogLik_P2(psi_j[i],D)
    
    if ((log_psihat-logpsi) < 0  ){
      r1<- 0
    } else {
      r1<-sign(psi_hat-psi_j[i])*sqrt(2*(log_psihat-logpsi))
    }
    
    if ( r1 <= r_root2(psi_j[i],data2)  ){
      psi_j[i+1] <- psi_j[i] + (c*alpha)/i
    } else {
      psi_j[i+1] <- psi_j[i] - (c*(1-alpha))/i
    }
    diff[i] = abs(psi_j[i+1] - psi_j[i])
  }
  
  return(list(psi_j,diff))
  #return(psi_j)
}

result1<-CI12(0.01)

############### Average the CI with 10 replications based on first method################################
rep=10
CI_mat <- matrix(data = NA, nrow = rep, ncol = 2)
Width_CI12<- function(alpha){
  for (i in 1:rep){
    result_L<-CI12(alpha)
    index_L<-which(abs(result_L[[2]] - 0.0001) == sort(abs(result_L[[2]] - 0.0001), decreasing = FALSE)[1])
    LL<-result_L[[1]][index_L]
    
    result_U<-CI12(1-alpha)
    index_U<-which(abs(result_U[[2]] - 0.0001) == sort(abs(result_U[[2]] - 0.0001), decreasing = FALSE)[1])
    UL<-result_U[[1]][index_U]
    CI_mat[i,]<-c(LL,UL)
  }
  return(CI_mat)
}

limits1<-Width_CI12(0.01)
CI_average1<-colMeans(limits1)
CI_average1
CI_average1[2]-CI_average1[1]
#one gives (-2.179041 35.890549)

###########Check fist simulation gives wider CIs###############################################
CI_plot1<-c()
ite=50
CI_average1<-matrix(NA,ite,2)

for (i in 1:ite){
  CI_average1[i,]<-colMeans(Width_CI12(0.01))
  CI_plot1[i]<-CI_average1[i,2]-CI_average1[i,1]
}

\


