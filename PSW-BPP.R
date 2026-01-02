library(mvtnorm)
library(invgamma)
library(MatchIt)
library(rstantools)
library(dplyr)
library(MASS)
library(distr)
library(stats)

start_time <- Sys.time()
rep=100
ATE.rep = c(); CP.ATE = c(); width.CI.rep = c()
var.rep = c(); CP.var = c(); var.width.CI.rep = c()
a1.rep = c(); a2.rep = c()
for(r in 1:rep){
  set.seed(1234+r)
  n=200 #current data sample size
  p=4 #without intercept
  n0=2000 #external data sample size
  b=c(2, 1, 1.5, -1.3)
  
  var_y=10; var_y0=10
  #current data
  mu=c(1,1.2,1.5,1.6)
  var_x=1
  cov_matrix <- matrix(0, p, p)
  diag(cov_matrix) <- var_x
  
  X=rmvnorm(n, mean=mu, sigma=cov_matrix) #design matrix without intercept
  colnames(X) =paste0("V", 1:p)
  X <- as.data.frame(X) %>%
    mutate(V1 = ifelse(V1 < median(V1), 0, 1)) #converting the continuous X1 covariate into binary
  
  e=rnorm(n,0, sd=sqrt(var_y))
  Y=as.matrix(X)%*%b+e
  
  #external data
  mu0=c(1,1,1) #no treatment variable is generated
  var_x0=1
  cov_matrix0 <- matrix(0, p-1, p-1)
  diag(cov_matrix0) <- var_x0
  
  X0=rmvnorm(n0, mean=mu0, sigma=cov_matrix0) #design matrix without intercept
  X0= cbind(rep(0, n0), X0) #treatment (only control arm) variable is included
  colnames(X0) =paste0("V", 1:p)
  
  e0 <- rnorm(n0, 0, sd=sqrt(var_y0))  # Error term
  Y0=as.matrix(X0)%*%b+e0
  
  study=c(rep(1,n), rep(0,n0)) #outcome for propensity score model; 1 for current data, 0 for external data
  fulldata=data.frame(cbind(rbind(X,X0),study)) #combined data (both current and external)
  colnames(fulldata)=c(paste0("X", 1:p),"study")
  p.score = fitted(glm(study ~ .-X1, data = fulldata, family = binomial)) 
  wght = pmin(1, (p.score/(1-p.score))) #Weight vector
  fulldata$wght <- wght
  fulldata$Y=c(Y,Y0)
  
  fulldata_c = subset(fulldata, study == 1) #current data (both treatment & control)
  fulldata_e = subset(fulldata, study == 0) #external control
  
  Xt = subset(fulldata_c, X1==1)[, 1:p] #current treatment data
  Yt = subset(fulldata_c, X1==1)[, "Y"]
  Xc = subset(fulldata_c, X1==0)[, 1:p] #current control data
  Yc = subset(fulldata_c, X1==0)[, "Y"]
  df_mean_c = length(Yc)-1 #t-df
  mu_mean_c = mean(Yc) #observed sample mean
  sigma_mean_c = sqrt(var(Yc)/length(Yc)) #standard error of the sample mean
  
  Xe = as.matrix(fulldata_e[,1:p]) #external control
  Ye = fulldata_e$Y
  df_mean_e = length(Ye)-1 #t-df
  mu_mean_e = mean(Ye) #observed sample mean
  sigma_mean_e = sqrt(var(Ye)/length(Ye)) #standard error of the sample mean
  
  t_tail <- function() {
    if(df_mean_e <= 0L) return(a1=0) # if there are less than or equal to one subject in the external data #then both are 0
    
    diff_mean <- (distr::Td(df_mean_c) * sigma_mean_c) - #distr::Td(df) is a Student-t distribution object (mean 0, symmetric).
      (distr::Td(df_mean_e) * sigma_mean_e)
    left_mean <- distr::p.l(diff_mean)(mu_mean_e - mu_mean_c) #distr::p.l(diff_mean) returns the lower-tail CDF function for diff_mean
    
    p_mean <- min(left_mean, 1-left_mean) # tail probability
    return(p_mean)
  }
  p_mean = t_tail()
  a1 = if_else(p_mean < 0.01, 0, p_mean)
  a1.rep = c(a1.rep, a1)
  
  #########################
  f_tail <- function(df1, df2, target_ratio) {
    if(df1 <= 2) return(0)
    mode <- ((df1-2)/df1) * (df2/(df2+2))
    
    df_val <- stats::df(target_ratio, df1, df2)
    f <- function(x) {
      (stats::df(x, df1, df2) - df_val)^2
    }
    if(target_ratio > mode) {
      q <- optim(par = mode/2, fn = f, method = "Brent",
                  lower = 0, upper = mode)
      if(!(near(q$value, 0))) stop("not convergent!")
      ans <- min(
        pf(q$par, df1, df2),
        1 - pf(target_ratio, df1, df2)
      )
    } else if (target_ratio < mode) {
      q <- optim(par = 2*mode - target_ratio, fn = f,
                  method = "Brent", lower = mode,
                  upper = qf(1-1e-6, df1, df2))
      if(!(near(q$value, 0))) stop("not convergent!")
      ans <- min(
        pf(target_ratio, df1, df2),
        1 - pf(q$par, df1, df2)
      )
    }
  }
  p_var <- f_tail(df1 = df_mean_e, df2 = df_mean_c, target_ratio = (var(Ye)/var(Yc)))
  a2 = if_else(p_var < 0.01, 0, p_var)
  a2.rep = c(a2.rep, a2)
  
  iter=5000
  beta=rep(1,p); sigsq=1.1
  
  #current data
  X=as.matrix(rbind(Xt, Xc))
  Y=c(Yt, Yc)
  
  #external data
  #Xe, Ye are already defined above
  We=diag(fulldata_e$wght) #PS-weight diagonal matrix
  We_sum=sum(fulldata_e$wght) #sum  of PS
  
  
  hatbeta_c=solve(t(X)%*%X)%*%t(X)%*%Y
  hatbeta_e = ginv(t(Xe)%*%We %*% Xe) %*% t(Xe)%*% We %*% Ye
  
  #Empty vectors and matrices
  post_sample=matrix(NA, iter, (p+1))
  
  #MCMC samples
  for(i in 1:iter){
    beta_var = sigsq*solve(t(X)%*%X+a1*t(Xe)%*%We%*%Xe)
    beta_mean = solve(t(X)%*%X+a1*t(Xe)%*%We%*%Xe) %*% (t(X)%*%X%*%hatbeta_c + a1*t(Xe)%*%We%*%Xe%*%hatbeta_e)
    beta = as.vector(rmvnorm(1, mean=beta_mean, sigma=beta_var))
    
    #sigsq_shape = (nrow(X)+(p*a1+(1-p)*a2)*We_sum)/2
    sigsq_shape = (nrow(X)+p*a1+(We_sum-p)*a2)/2
    sigsq_rate = (t(Y-X%*%hatbeta_c)%*%(Y-X%*%hatbeta_c)+a2*t(Ye-Xe%*%hatbeta_e)%*%We%*%(Ye-Xe%*%hatbeta_e))/2
    sigsq = rinvgamma(1,  shape=sigsq_shape, rate=sigsq_rate)
    post_sample[i,] = c(beta, sigsq)
  }
  trt.vec = post_sample[-(1:2000), 1]
  ATE = mean(trt.vec)
  ATE.rep = c(ATE.rep, ATE)
  ATE.quantiles = quantile(trt.vec, probs=c(0.025, 0.975)) #95% credible interval
  width.CI = ATE.quantiles[2]-ATE.quantiles[1]
  width.CI.rep = c(width.CI.rep, width.CI)
  
  true.trt=b[1] #true treatment effect
  
  if (true.trt >= ATE.quantiles[1] && true.trt <= ATE.quantiles[2]) {
    CP.ATE[r] <- 1
  } else {
    CP.ATE[r] <- 0
  }
  
  var.vec = post_sample[-(1:2000), (p+1)]
  var.mean = mean(var.vec)
  var.rep = c(var.rep, var.mean)
  var.quantiles = quantile(var.vec, probs=c(0.025, 0.975)) #95% credible interval
  var.width.CI = var.quantiles[2]-var.quantiles[1]
  var.width.CI.rep = c(var.width.CI.rep, var.width.CI)
  
  true.var.y=var_y #true variance
  
  if (true.var.y >= var.quantiles[1] && true.var.y <= var.quantiles[2]) {
    CP.var[r] <- 1
  } else {
    CP.var[r] <- 0
  }
}
time_taken <- Sys.time()-start_time

median(a1.rep); IQR(a1.rep) #median and IQR for a1
median(a2.rep); IQR(a2.rep) #median and IQR for a2

mean(ATE.rep); mean(true.trt-ATE.rep); mean(abs(true.trt-ATE.rep)) # mean, bias, and ABias of treatment effect
sqrt(mean((true.trt-ATE.rep)^2)); sd(ATE.rep); mean(width.CI.rep); sum(CP.ATE) #RMSE, SE, width, and CP of treatment effect
