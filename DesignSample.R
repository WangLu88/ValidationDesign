library(MASS)
t <- 1

roc_seq <- function(X_D1,X_D0,t){
  kk <- ceiling((1-t)*length(X_D0))
  S_inv_D0 <-  sort(X_D0)[kk]
  roc <- 1- mean(ifelse(X_D1<=S_inv_D0,1,0))
  D <- c(rep(1,length(X_D1)),rep(0,length(X_D0)))
  X <- c(X_D1,X_D0)
  fit <- glm(D ~ X, family = "binomial")
  ratio <- exp(fit$coefficients[1]+S_inv_D0*fit$coefficients[2])*(length(X_D0)/length(X_D1))
  sigma2 <- roc*(1-roc)/length(X_D1) + (ratio^2)*t*(1-t)/length(X_D0)
  result <- (roc-roc_0)/sqrt(sigma2)
  return(result)
}
## bound=1: P==1&stop==1 ##
## bound=2: P==0.5&stop==1 ##
## bound=3: P==1&stop==2 ##
## bound=4: P==0.5&stop==2 ##

n <- 240
mu_0 <- 0
sd_0 <- sd_1 <- 1
mu_1 <- 0.9
mu_null <- 0.6
roc_0 <- 1-pnorm(qnorm(1-t,mu_0,sd_0),mu_null,sd_1)
nrep <- 5000
alpha <- 0.05
max <- 50   ## available specimen volume
x_dim <- 3*max
result_r0.3 <- result_r0.5 <- matrix(NA,nrep,9)

set.seed(123)

for (k in 1:nrep){

## Generate Data ##  
n_D0 <- n_D1 <- n/2 
x_D0 <- mvrnorm(n_D0,rep(mu_0,x_dim),diag(rep(sd_0,x_dim)))
x_D1 <- mvrnorm(n_D1,rep(mu_1,x_dim),diag(rep(sd_1,x_dim)))

## Default Strategy ##
# K1 = 1 or 0, reject null or not
K1 <- rep(NA,max)
for (i in 1:max){
  z_K1 <- roc_seq(x_D1[,i],x_D0[,i],t)
  K1[i] <- ifelse(z_K1>=qnorm(1-alpha,0,1),1,0)
}

## Proposed Strategy ##

# K = 1/r, r is the proportion of samples used in the first stage #
for (K in 2:3) {
  
  if(K==3){
    ii1_case <- sample(1:n_D1,n_D1/K,replace=FALSE)
    ii1_ctrl <- sample(1:n_D0,n_D0/K,replace=FALSE)
    ii2_case <- sample((1:n_D1)[-ii1_case],n_D1/K,replace=FALSE)
    ii2_ctrl <- sample((1:n_D0)[-ii1_ctrl],n_D0/K,replace=FALSE)
    list_case <- list(sort(ii1_case),sort(ii2_case),(1:n_D1)[-c(ii1_case,ii2_case)])
    list_ctrl <- list(sort(ii1_ctrl),sort(ii2_ctrl),(1:n_D0)[-c(ii1_ctrl,ii2_ctrl)])
  }
  if(K==2){
    ii1_case <- sample(1:n_D1,n_D1/2,replace=FALSE)
    ii1_ctrl <- sample(1:n_D0,n_D0/2,replace=FALSE)
    list_case <- list(sort(ii1_case),(1:n_D1)[-ii1_case])
    list_ctrl <- list(sort(ii1_ctrl),(1:n_D0)[-ii1_ctrl])
  }

# Different types of boundary shapes and stopping rules #
  for (bound in 1:4){
    
    if (bound==1){
      boundary_r0.3 <- matrix(c(2.8596698, -0.9532233, 1.651031, 1.651031),2,2)
      boundary_r0.5 <- matrix(c(2.358497, 0, 1.667709, 1.667709),2,2)
    }
    if (bound==2){
      boundary_r0.3 <- matrix(c(1.8730601, 0.2897634, 1.87306, 1.87306),2,2)
      boundary_r0.5 <- matrix(c(1.8342166, 0.7597574, 1.834217, 1.834217),2,2)
    }
    if (bound==3){
      boundary_r0.3 <- matrix(c(Inf, -0.9610118, 1.642354, 1.642354),2,2)
      boundary_r0.5 <- matrix(c(Inf, -0.03144846, 1.633478, 1.633478),2,2)
    }
    if (bound==4){
      boundary_r0.3 <- matrix(c(Inf, 0.1131356, 1.589686, 1.589686),2,2)
      boundary_r0.5 <- matrix(c(Inf, 0.5658151, 1.577009, 1.577009),2,2)
    }
    if(K==3){boundary <- boundary_r0.3}
    if(K==2){boundary <- boundary_r0.5}
   
  # K3 = 1 or 0, reject null or not # 
  temp <- matrix(0,1,K)
  K3 <- rep(NA)
  for (i in 1:x_dim){
    aa <- apply(temp,2,sum)
    if (length(which(aa==min(aa)))==1) {ind <- which(aa==min(aa))}
    if (length(which(aa==min(aa)))!=1) {ind <- sample(which(aa==min(aa)),1,replace=FALSE)}
    if (min(aa)==max) {break}
    z_K3 <- roc_seq(x_D1[list_case[[ind]],i],x_D0[list_ctrl[[ind]],i],t)
    z_K3_final <- roc_seq(x_D1[,i],x_D0[,i],t)
    temp_new <- rep(1,K)
    K3_new <- rep(NA)
    if (z_K3>=boundary[1,1]) {
      temp_new[-ind] <- 0
      K3_new <- 1}
    if (z_K3<=boundary[2,1]) {
      temp_new[-ind] <- 0
      K3_new <- 0}
    if (z_K3>boundary[2,1]&z_K3<boundary[1,1]) {
      if (z_K3_final>=boundary[1,2]) {K3_new <- 1}
      if (z_K3_final<boundary[2,2]) {K3_new <- 0}
    } 
    temp <- rbind(temp,temp_new)
    K3 <- c(K3,K3_new)
    aa <- apply(temp,2,sum)
    if (min(aa)==max) {break}  
  }
  temp <- temp[-1,]
  K3 <- K3[-1]
  repeat{
    aa <- apply(temp,2,sum)
    if(max(aa)==max) {break}
    nn <- nrow(temp)
    if(max(aa)>max) {
      temp <- temp[-nn,]
      K3 <- K3[-nn]} 
  }
  
  if(K==3){
    result_r0.3[k,bound] <- length(K3)
    result_r0.3[k,bound+4] <- sum(K3)
    result_r0.3[k,9] <- sum(K1)}
  if(K==2){
    result_r0.5[k,bound] <- length(K3)
    result_r0.5[k,bound+4] <- sum(K3)
    result_r0.5[k,9] <- sum(K1)}
}
}

}
apply(result_r0.3,2,mean)
apply(result_r0.5,2,mean)
