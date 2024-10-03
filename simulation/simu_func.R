# simulation functions
## simulation function
simutwoX <- function(eff1 = 2,n = 500,p=500,q=10,r=5){
  U = matrix(rnorm(n*r),nrow = n,ncol = r)
  # x2 = matrix(rnorm(n,0,1),nrow = n,ncol = 1)
  sigma2 = invgamma::rinvgamma(p,3,2)
  if((q>=10)&(q>r)){
    alpha_ = rstiefel::rustiefel(q,r)%*%(sqrt(q)*diag(3-2*((1:r)-1)/(r-1)))
  }else{
    alpha_ = matrix(1/sqrt(r),nrow = q,ncol = r)
  }
  # alpha2 = matrix(2/sqrt(r),nrow = r,ncol = 1)
  # pvec = runif(p*q,0,0.3)
  # beta_ = matrix(eff1*sqrt(2)*rbinom(p*q,1,0.05),nrow = p,ncol = q)
  beta_ = matrix(0,nrow = p,ncol = q)
  for(i in 1:q) beta_[,i] = eff1*sqrt(2)*rbinom(p,1,0.15)
  # beta2 = matrix(eff2*sqrt(2)*rbinom(p,1,0.1),nrow = p,ncol = 1)
  gamma_ = rstiefel::rustiefel(p,r)%*%(sqrt(p)*diag(3-2*((1:r)-1)/(r-1))) # crucial step, random choice of such matrix results in bad behavior of all methods
  tauset_comp = apply(beta_,2,function(x) which(x!=0))
  X = U%*%t(alpha_)+matrix(rnorm(n*q),nrow = n,ncol = q) #;x2 = rbinom(n,1,plogis(U%*%alpha2+rnorm(n)))
  Y = X%*%t(beta_)+U%*%t(gamma_)+MASS::mvrnorm(n,rep(0,p),diag(sigma2))
  return(list(data = list(X = X,Y = Y, U = U),
              coef = list(beta = beta_,alpha = alpha_,
                          gamma = gamma_,tauset_comp = tauset_comp)))
}

# mixed effect simulation
simutwoX_mix <- function(eff1 = -1,eff2 = 5,n = 500,p=100,q=10,r=3){
  U = matrix(rnorm(n*r),nrow = n,ncol = r)
  # x2 = matrix(rnorm(n,0,1),nrow = n,ncol = 1)
  sigma2 = invgamma::rinvgamma(p,3,2)
  if((q>=10)&(q>r)){
    alpha_ = rstiefel::rustiefel(q,r)%*%(sqrt(q)*diag(3-2*((1:r)-1)/(r-1)))
  }else{
    alpha_ = matrix(1/sqrt(r),nrow = q,ncol = r)
  }
  # alpha2 = matrix(2/sqrt(r),nrow = r,ncol = 1)
  # pvec = runif(p*q,0,0.3)
  # beta_ = matrix(eff1*sqrt(2)*rbinom(p*q,1,0.05),nrow = p,ncol = q)
  beta_ = matrix(0,nrow = p,ncol = q)
  nonzero = 0.15*p
  eff1_count = floor(nonzero*0.6)
  for(i in 1:q) {
    indice = sample(1:p,nonzero)
    beta_[indice[1:eff1_count],i] = eff1*sqrt(2)
    beta_[indice[(eff1_count + 1):nonzero],i] = eff2*sqrt(2)
  }
  # beta2 = matrix(eff2*sqrt(2)*rbinom(p,1,0.1),nrow = p,ncol = 1)
  gamma_ = rstiefel::rustiefel(p,r)%*%(sqrt(p)*diag(3-2*((1:r)-1)/(r-1))) # crucial step, random choice of such matrix results in bad behavior of all methods
  tauset_comp = apply(beta_,2,function(x) which(x!=0))
  X = U%*%t(alpha_)+matrix(rnorm(n*q),nrow = n,ncol = q) #;x2 = rbinom(n,1,plogis(U%*%alpha2+rnorm(n)))
  Y = X%*%t(beta_)+U%*%t(gamma_)+MASS::mvrnorm(n,rep(0,p),diag(sigma2))
  return(list(data = list(X = X,Y = Y, U = U),
              coef = list(beta = beta_,alpha = alpha_,
                          gamma = gamma_,tauset_comp = tauset_comp)))
}



# for uniform residuals
simutwoX_unif <- function(eff1 = 2,n = 500,p=500,q=1,r=5){
  U = matrix(rnorm(n*r),nrow = n,ncol = r)
  # x2 = matrix(rnorm(n,0,1),nrow = n,ncol = 1)
  # sigma2 = invgamma::rinvgamma(p,3,2)
  alpha_ = matrix(1/sqrt(r),nrow = r,ncol = q)
  # alpha2 = matrix(2/sqrt(r),nrow = r,ncol = 1)
  # pvec = runif(p*q,0,0.3)
  # beta_ = matrix(eff1*sqrt(2)*rbinom(p*q,1,0.05),nrow = p,ncol = q)
  beta_ = matrix(0,nrow = p,ncol = q)
  for(i in 1:q) beta_[,i] = eff1*sqrt(2)*rbinom(p,1,0.15)
  # beta2 = matrix(eff2*sqrt(2)*rbinom(p,1,0.1),nrow = p,ncol = 1)
  gamma_ = rstiefel::rustiefel(p,r)%*%(sqrt(p)*diag(3-2*((1:r)-1)/(r-1))) # crucial step, random choice of such matrix results in bad behavior of all methods
  tauset_comp = apply(beta_,2,function(x) which(x!=0))
  X = U%*%(alpha_)+matrix(runif(n*q),nrow = n,ncol = q) #;x2 = rbinom(n,1,plogis(U%*%alpha2+rnorm(n)))
  Y = X%*%t(beta_)+U%*%t(gamma_)+matrix(runif(n*p),nrow = n,ncol = p)
  return(list(data = list(X = X,Y = Y, U = U),
              coef = list(beta = beta_,alpha = alpha_,
                          gamma = gamma_,tauset_comp = tauset_comp)))
}



######################################
# evaluating the estimation of S set
# parallel function
process_estS <- function(b) {
  # parameters
  # n = 500; p = 500
  n = n; p = p; q = q
  eff = 1; r = 3
  # function
  testHDX <- simutwoX(eff,n,p,q,r)
  tryCatch({
    if(length(unique(unlist(testHDX$coef$tauset_comp)))==ncol(testHDX$data$Y)){
      stop("No NCO")
    }else{
      greed = greedy_l0(testHDX$data$Y,q,testHDX$data$X)
      estS = greed$truest
      trueS = unique(unlist(testHDX$coef$tauset_comp))
      sens = length(intersect(setdiff(1:p,estS),setdiff(1:p,trueS)))/length(setdiff(1:p,trueS))
      FDR = length(intersect(setdiff(1:p,estS),trueS))/length(setdiff(1:p,estS))
      FDRp = length(intersect(setdiff(1:p,estS),trueS))/p
      return(list(greed,trueS,sens,FDR,FDRp))
    }
  },error = function(e){
    # cat("Error:", e$message, "\n")
    return(e$message)
  })
  # combined_diff$condition <- c(names(testcoef1), names(testcoef2))
  
}


# parallel for mixed effect
process_function_mix <- function(b) {
  # parameters
  n = 500; p = 100
  q = q; r = r
  # function
  testHDX <- simutwoX_mix(eff1=1,eff2=5,n,p,q,r)
  tryCatch({
    if(length(unique(unlist(testHDX$coef$tauset_comp)))==ncol(testHDX$data$Y)){
      stop("No NCO")
    }else{
      parallest <- debias_twostage(testHDX$data$X, testHDX$data$Y, testHDX$data$U, unique(unlist(testHDX$coef$tauset_comp)))
      deconfoundest <- debias_deconfound(testHDX$data$X, testHDX$data$Y, testHDX$data$U, unique(unlist(testHDX$coef$tauset_comp)))
      
      testcoef1 <- est_bias(parallest)
      testcoef2 <- est_bias(deconfoundest)
      
      combined_diff <- as.data.frame(rbind(do.call(rbind, lapply(testcoef1$MSE, mean,na.rm = T)),
                                           do.call(rbind, lapply(testcoef2$MSE, mean,na.rm = T))))
      infy_diff <- as.data.frame(rbind(do.call(rbind, lapply(testcoef1$infy_bias, mean,na.rm = T)),
                                       do.call(rbind, lapply(testcoef2$infy_bias, mean,na.rm = T))))
      return(list(parallest,deconfoundest,testcoef1,testcoef2,combined_diff,infy_diff,testHDX$coef$beta))
    }
  },error = function(e){
    # cat("Error:", e$message, "\n")
    return(e$message)
  })
  # combined_diff$condition <- c(names(testcoef1), names(testcoef2))
  
}


process_function_unif <- function(b) {
  # parameters
  n = 500; p = 100
  eff = eff; q = q; r = r
  # function
  testHDX <- simutwoX_unif(eff,n,p,q,r)
  tryCatch({
    if(length(unique(unlist(testHDX$coef$tauset_comp)))==ncol(testHDX$data$Y)){
      stop("No NCO")
    }else{
      parallest <- debias_twostage(testHDX$data$X, testHDX$data$Y, testHDX$data$U, unique(unlist(testHDX$coef$tauset_comp)))
      deconfoundest <- debias_deconfound(testHDX$data$X, testHDX$data$Y, testHDX$data$U, unique(unlist(testHDX$coef$tauset_comp)))
      
      testcoef1 <- est_bias(parallest)
      testcoef2 <- est_bias(deconfoundest)
      
      combined_diff <- as.data.frame(rbind(do.call(rbind, lapply(testcoef1$MSE, mean,na.rm = T)),
                                           do.call(rbind, lapply(testcoef2$MSE, mean,na.rm = T))))
      infy_diff <- as.data.frame(rbind(do.call(rbind, lapply(testcoef1$infy_bias, mean,na.rm = T)),
                                       do.call(rbind, lapply(testcoef2$infy_bias, mean,na.rm = T))))
      return(list(parallest,deconfoundest,testcoef1,testcoef2,combined_diff,infy_diff,testHDX$coef$beta))
    }
  },error = function(e){
    # cat("Error:", e$message, "\n")
    return(e$message)
  })
  # combined_diff$condition <- c(names(testcoef1), names(testcoef2))
  
}



