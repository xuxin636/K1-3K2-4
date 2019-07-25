set.seed(1234)
library(mvtnorm)
cond <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
N=3000;J <- 20;K1 <- 3;K2 <- 4;K <- K1+K2;JJ_ind <- J/K2
Q <- cbind(matrix(rep(1,J*K1),nrow=J,ncol = K1),matrix(0,nrow=J,ncol=K2))
for(k in 1:K2){
  Q[(1+(k-1)*JJ_ind):(k*JJ_ind),k+K1] <- 1
}
Q[(J-1):J,K1] <- 0;Q[J,K1-1] <- 0#for identification
D_true <- cbind(rnorm(J,0,1),matrix(runif(J*K,0.5,2),nrow=J,ncol=K)*Q )
rho_true <- 0.3
sigma_true=diag(1,nrow=K);sigma_true[1:K1,1:K1][row(sigma_true[1:K1,1:K1])!=col(sigma_true[1:K1,1:K1])]=rho_true
theta_true<-rmvnorm(N,mean=rep(0,K),sigma=sigma_true)#specific
set.seed(cond)
response <- matrix(rbinom(N*J,1,1/(1+exp(-(cbind(rep(1,N),theta_true)%*%t(D_true))))),nrow=N,ncol=J)
S <- K1+1;KK <- 14;theta_min <- -4;theta_max <- 4;mm  <- seq(theta_min,theta_max,(theta_max-theta_min)/KK);

THETA_tuta <- matrix(0,nrow=length(mm)^S,ncol=S);
for(k in 1:S){THETA_tuta[,S-k+1] <- rep(c(rep(1,length(mm)^(k-1))%*%t(mm)),length(mm)^(S-k))}
THETA_tuta <- cbind(rep(1,nrow(THETA_tuta)),THETA_tuta)
D_initial <- cbind(sort(rnorm(J,0,1))[rank(colMeans(response))],matrix(runif(J*K,0.5,1.5),nrow=J,ncol=K)*Q)
D_initial[J,(K-1)] <- 0#for identification
response <- t(response)
A_0 <- t(D_initial)
THETA_tuta_12 <- t(matrix(t(THETA_tuta),nrow=ncol(THETA_tuta)*length(mm))[2:(K1+1),])
rho <- runif(1,0.25,0.75)
t_0 <- function(A_0,indicator){
  temp_0_2 <- THETA_tuta%*%rbind(A_0[1:(K1+1),indicator],matrix(A_0[(K1+2):(K+1),indicator][which(A_0[(K1+2):(K+1),indicator]!=0)],nrow=1))
  cc2 <- temp_0_2%*%response[indicator,]-rowSums(log(1+exp(temp_0_2)))-THETA_tuta[,ncol(THETA_tuta)]*THETA_tuta[,ncol(THETA_tuta)]/2
  return(list(temp_0_2,cc2))
};
lik <- function(A,post,res){
  ts <- THETA_tuta%*%c(A[1:(K1+1)],A[(K1+2):(K+1)][which(A[(K1+2):(K+1)]!=0)])
  likelihood <- sum((ts%*%t(res)-c(log(1+exp(ts))))*post)
  return(likelihood)
}
prox <- function(x,lammda){
  for(k in 1:length(x)){
    if(abs(x[k])<lammda){x[k] <- 0}
    else{
      if(x[k]>=lammda){x[k]=x[k]-lammda}
      else{x[k]=x[k]+lammda}
    }
  }
  
  return(x)
}

lammda <- 0;alpha <- 0.5

timestart <- Sys.time()
cc <- A_0
ind <- matrix(0,nrow=(K+2),ncol=J)
s <- length(mm);grad_temp <-matrix(0,nrow=(K+1),ncol=J)
stepsize <-  matrix(1/N,nrow=(K+1),ncol=J)
x3 <- rowSums(response);
temp <- list();cc1 <- list();u_temp <- array(dim=c(nrow(THETA_tuta),N,K2))
post <- x <- list()
for(m in 1:500){
  
  sigma=diag(1,nrow=K1);sigma[row(sigma)!=col(sigma)]=rho
  density <- (rowSums((THETA_tuta_12%*%solve(sigma))*THETA_tuta_12)/2)[c(rep(1,s)%*%t(1:nrow(THETA_tuta_12)))]
  
  for(k in 1:K2){
    indicator=which(A_0[k+K1+1,]!=0)
    s_temp <- t_0(A_0,indicator = indicator);
    cc1[[k]] <-s_temp[[2]]; 
    temp[[k]] <- t((1/(1+exp(-s_temp[[1]]))))
    u_temp[,,k] <- (log(matrix(colSums(matrix(exp(cc1[[k]]),nrow=s)),ncol=N)))[c(rep(1,s)%*%t(1:nrow(THETA_tuta_12))),]
  }
  
  for(k in 1:K2){
    indicator=which(Q[,k+K1]!=0);indicator_1 <- c(2:(K1+1),k+K1+1)
    post_temp <- exp(cc1[[k]]+apply(u_temp[,,-k],c(1,2),sum)-density)
    post[[k]] <- sweep(post_temp,2,colSums(post_temp),"/")
    x <- (t(THETA_tuta[,-1])%*%post[[k]])%*%t(response[indicator,])
    grad_temp[1,indicator] <- x3[indicator]-rowSums(temp[[k]]%*%post[[k]])
    for(mm in 1:(K1+1)){
      grad_temp[indicator_1[mm],indicator] <- x[mm,]-rowSums(temp[[k]]%*%(post[[k]]*THETA_tuta[,mm+1]))
    }
    
  }
  rho <- mean(colSums(post[[1]]*THETA_tuta[,2]*THETA_tuta[,3]))
  grad_temp[(2):(K+1),] <- grad_temp[(2):(K+1),]*t(Q)
  
  for(j in 1:J){
    indicator <- which(A_0[(K1+2):(K+1),j]!=0)
    lik_0 <- lik(A_0[,j],post[[indicator]],response[j,])
    for(k in 1:(K1+1)){
      tuta <- A_0[,j]
      tuta[k] <- prox(A_0[k,j]+grad_temp[k,j]*stepsize[k,j],lammda = lammda)
      while((lik(tuta,post[[indicator]],response[j,])-lik_0)<grad_temp[k,j]*(tuta[k]-A_0[k,j])-0.5*(tuta[k]-A_0[k,j])*(tuta[k]-A_0[k,j])/stepsize[k,j]){
        stepsize[k,j] <- stepsize[k,j]*alpha;
        tuta[k] <- prox(A_0[k,j]+grad_temp[k,j]*stepsize[k,j],lammda = lammda)
      }
      cc[k,j] <- tuta[k]
    }
    k <- indicator+(K1+1)
    tuta <- A_0[,j]
    tuta[k] <- prox(A_0[k,j]+grad_temp[k,j]*stepsize[k,j],lammda = lammda)
    while((lik(tuta,post[[indicator]],response[j,])-lik_0)<grad_temp[k,j]*(tuta[k]-A_0[k,j])-0.5*(tuta[k]-A_0[k,j])*(tuta[k]-A_0[k,j])/stepsize[k,j]){
      stepsize[k,j] <- stepsize[k,j]*alpha;
      tuta[k] <- prox(A_0[k,j]+grad_temp[k,j]*stepsize[k,j],lammda = lammda)
    }
    cc[k,j] <- tuta[k]
    
    
  }
  
  cc[(2):(K+1),] <- cc[(2):(K+1),]*t(Q)
  
  A_0 <- cc
}
timeend <- Sys.time() 
tt <- timeend-timestart
tt

ind[1:(K+1),] <- A_0;
ind[K+2,1] <- rho


write.csv(RESULT, file =paste0('K1=3K2=4',cond,'.csv'))

