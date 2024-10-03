# real data debiasing
# debias function
debias_twostage_data <- function(X,Y,ptrue = 1:15){
  ## simple linear regression with X and U
  # beta_orc = lm(Y~X+U)$coef[2:(ncol(X)+1),]
  
  ## linear regression with only X
  beta_x = lm(Y~X)$coef[2:(ncol(X)+1),]
  naivep <- apply(Y,2,function(y_col) {
    summary(lm(y_col~X))$coef[2:(ncol(X)+1),4]
  })
  
  # factor analysis for X
  # hatr = nScree(cor(X))$Components$nkaiser
  
  pzero = setdiff(1:ncol(Y),ptrue)
  ## parallel Y 2023
  betal=pl=betarg=
    betargms=betalasso=epsilon=betalms=
    betascad = betax =  miaobeta = 
    matrix(0,ncol = ncol(Y),nrow = ncol(X))
  
  # for(i in 1:ncol(Y)){
  #   epsilon[,i] = lm(Y[,i]~X)$coef[2:(ncol(X)+1)]
  #   # multiple trt debias
  #   miaobeta[,i] = esti.null(Y[,i],X,v = NULL,hatr)
  # }
  
  for(i in ptrue){
    # step one
    # X = scale(X);Y = scale(Y)
    W = Y[,pzero]
    Zl = Y[,setdiff(ptrue,i)]
    modelone = lm(W~X+Zl)
    What = predict(modelone)
    # step two
    # modeltwo = lm(Y[,i]~X+What)
    # betal[,i] = modeltwo$coef[2:(ncol(X)+1)]
    # betax[,i] = esti.null(Y[,i],X,v = NULL,hatr)
    # pl[,i] = summary(modeltwo)$coef[2:(ncol(X)+1),4]
    param = cv.glmnet(as.matrix(cbind(X,What)),Y[,i],alpha = 0)
    betarg[,i] = coef(glmnet(as.matrix(cbind(X,What)),Y[,i],alpha = 0,lambda = param$lambda.min))[2:(ncol(X)+1)]
    # adaptive LASSO
    # betargms[,i] = (MASS::lm.ridge(Y[,i]~X+What,lambda = 0.01))$coef[1:ncol(X)]
    cv_model <- cv.glmnet(as.matrix(cbind(X,What)),Y[,i], alpha = 1)
    best_model <- as.vector(coef(glmnet(as.matrix(cbind(X,What)), Y[,i], alpha = 1, lambda = cv_model$lambda.min)))[-1]
    newx = sweep(as.matrix(cbind(X,What)), 2, best_model, '*')
    cv_modelstar <- cv.glmnet(newx,Y[,i], alpha = 1)
    betastar = coef(glmnet(newx, Y[,i], alpha = 1, lambda = cv_modelstar$lambda.min))[2:(ncol(X)+1)]
    betargms[,i] = c(best_model[1:(ncol(X))])*c(betastar)
    betalasso[,i] = best_model[1:(ncol(X))]
    # SCAD
    scad_model = cv.ncvreg(as.matrix(cbind(X,What)),Y[,i],seed = 233)
    scad_fit = (ncvreg(as.matrix(cbind(X,What)),Y[,i],lambda = scad_model$lambda))
    scad_lam <- scad_fit$lambda[which.min(BIC(scad_fit))]
    betascad[,i] = coef(scad_fit, lambda=scad_lam)[2:(ncol(X)+1)]
  }
  return(list(naive = beta_x,ridge = betarg,naivep = naivep,
              adlasso = betargms,lasso = betalasso,scad = betascad,direct = lm(Y~X)$coef[2:(ncol(X)+1),]))
}


# Wang's method
debias_deconfound_data <- function(X,Y,ptrue = 1:15){
  options(warn = -1)
  ## simple linear regression with X and U
  # beta_orc = lm(Y~X+U)$coef[2:(ncol(X)+1),]
  
  ## linear regression with only X
  # beta_x = lm(Y~X)$coef[2:(ncol(X)+1),]
  
  ## deconfounding test 2017
  ### negative control
  testx = qr(X)
  Yqr = qr.qty(testx,Y)
  testy = t(qr.Q(testx))%*%Y
  # fa
  fahat = psych::fa(cor(Yqr[(ncol(X)+1):nrow(Yqr),]),5)
  # use POET instead
  # fahat = POET::POET(Y = t(Yqr[2:nrow(Yqr),]),K = 5)
  gammahat = as.matrix(fahat$loadings)
  sigmahat = cor(fahat$residual) #fahat$SigmaU 
  sigmahatinv = diag(ncol(sigmahat))#solve(sigmahat)
  # alphahat_nc = solve(t(gammahat[pzero,])%*%sigmahatinv[pzero,pzero]%*%gammahat[pzero,])%*%t(gammahat[pzero,])%*%sigmahatinv[pzero,pzero]%*%(t(Yqr[1:ncol(X),pzero]))%*%t(solve(qr.R(testx)))
  # betahat_nc = (t(Y[1:ncol(X),ptrue]))%*%t(solve(qr.R(testx)))-gammahat[ptrue,]%*%alphahat_nc
  # betahat_nc_alt = lm(Y[1,ptrue]/testx$qr[1]~gammahat[ptrue,])$coef[1]
  
  ### sparsity
  robustmodel = apply(t(Yqr[1:ncol(X),])%*%t(solve(qr.R(testx))),2,function(x) rlm(x~gammahat)$coef)
  # scadmodel_selection = apply(t(Yqr[1:ncol(X),])%*%t(solve(qr.R(testx))),2,function(x) (ncvreg(unclass(gammahat),x,nlambda = 100)))
  # scad_lam <- lapply(scadmodel_selection, function(x) x$lambda[which.min(BIC(x))])
  # betascad = mapply(function(lam,model) coef(model, lambda=lam),scad_lam,scadmodel_selection)
  betahat_sps = t(Yqr[1:ncol(X),])%*%t(solve(qr.R(testx)))-gammahat%*%robustmodel[2:nrow(robustmodel),]
  # betahat_scad = t(Yqr[1:ncol(X),])%*%t(solve(qr.R(testx)))-gammahat%*%betascad[2:nrow(betascad),]
  # ## try LASSO for sparsity (LEAPP)
  
  return(list(robustrg = t(betahat_sps)))
}


esti.null <- function(y,x,v=NULL, nfact, dta=NULL, subset=NULL){
  if(is.null(dta)){
    if(is.null(subset)) subset <- 1:nrow(as.matrix(x))
    y <- as.matrix(y)[subset,]
    x <- as.matrix(x)[subset,]
  }else{
    if(is.null(subset)) subset <- 1:nrow(as.matrix(dta))
    y <- as.matrix(dta[subset,y])
    x <- as.matrix(dta[subset,x])
  }
  
  n <- nrow(as.matrix(x))
  dimx <- ncol(as.matrix(x))
  lambda <- log(n)/n
  
  # factor analysis
  if(is.null(v)){
    x.res <- x
    # ols Y~X
    lm.y <-lm(y~x) 
  }else{
    x.res <- lm(x~v)$res
    # ols Y~X+V
    lm.y <-lm(y~x+v) 
  }
  efa.fit <- as.matrix(factanal(x.res, factors=nfact)$loadings)
  halpha <- efa.fit*sqrt(diag(var(x.res)))
  hgamma <- solve(var(x.res))%*%halpha
  hxi <- lm.y$coef[2:(1+dimx)]
  # select confounded treatments
  confd.x <- apply(hgamma^2, 1, sum) > lambda
  if(sum(confd.x) < 2*nfact+1){
    hbeta <- hxi
    hdelta <- rep(0, nfact)
  }else{
    # estimate confounding bias using LMS
    delta.ini <-  lqs(y=hxi[confd.x], x=as.matrix(hgamma[confd.x,]), intercept=FALSE, method='lts')$coefficients
    # estimate effects
    beta.ini <- hxi - as.matrix(hgamma)%*%as.matrix(delta.ini)
    # selection of confounded null treatments
    confd.null.x <- confd.x
    confd.null.x[confd.x] <- order((beta.ini[confd.x])^2, decreasing=TRUE) > ((sum(confd.x)-nfact)/2)
    # update delta
    hgamma.confd.null <- as.matrix(hgamma[confd.null.x,])
    hxi.confd.null <-  hxi[confd.null.x]
    hdelta <- solve(t(hgamma.confd.null)%*%hgamma.confd.null) %*% t(hgamma.confd.null) %*% hxi.confd.null
    # final estimate
    hbeta <- hxi - as.matrix(hgamma)%*%hdelta
  }
  
  c(hbeta)
}




# estimation function
estbeta <- function(X,Y,X_demo = X_demo,B = 1000,rate = 0.9){
  # Y_FA <- get("Y_FA", envir = .GlobalEnv)
  # X_demo <- get(demo, envir = .GlobalEnv)
  Y_res = apply(Y,2,function(y){
    fit <- lm(y~X_demo)
    residuals(fit)
  })
  X_res = apply(X,2,function(y){
    fit <- lm(y~X_demo)
    residuals(fit)
  })
  greed1 = greedy_l0(Y_res,ncol(X_res),X_res,rate =rate)
  estS1 = greed1$truest
  parallest1 <- debias_twostage_data(X_res, Y_res, estS1)
  deconfoundest1 <- debias_deconfound_data(X_res, Y_res, estS1)
  
  # bootstrap interval
  process_function_boot_res <- function(b) {
    # parameters
    X = X_res; Y = Y_res; ptrue = estS1
    # bootstrap coverage probablity
    set.seed(233*b)
    boot_indices <- sample(1:nrow(X), replace = TRUE)
    X_boot <- X[boot_indices, , drop = FALSE]
    Y_boot <- Y[boot_indices, , drop = FALSE]
    
    # Apply the debias_twostage function to the bootstrap sample
    boot_estimates <- debias_twostage_data(X_boot, Y_boot, ptrue)
    boot_wang <- debias_deconfound_data(X_boot, Y_boot, ptrue)
    return(list(boot_estimates,boot_wang))
    
  }
  ncore = detectCores()
  cl <- makeCluster(ncore-2)
  clusterEvalQ(cl, {
    library(glmnet)
    library(MASS)
    library(ncvreg)
    library(nFactors)
  })
  clusterExport(cl, c("debias_twostage_data", "debias_deconfound_data", 
                      "X_res","Y_res","estS1","process_function_boot_res",
                      "esti.null"),envir = environment())
  boot_data_results_res <- parLapply(cl, 1:B, process_function_boot_res)
  stopCluster(cl)
  
  # evaluate the intervals
  bootmat <- do.call(rbind,lapply(boot_data_results_res,function(est) as.vector(est[[1]]$scad)))
  ci_lower <- apply(bootmat,2,quantile,probs=0.025)
  ci_upper <- apply(bootmat,2,quantile,probs=0.975)
  ci_lower2 <- apply(bootmat,2,quantile,probs=0.05)
  ci_upper2 <- apply(bootmat,2,quantile,probs=0.95)
  bootse <- apply(bootmat,2,sd)
  bootse_CI_scad_res = sapply(1:length(as.vector(parallest1$scad)),function(k){
    as.vector(parallest1$scad)[k]+c(-1,1)*1.96*bootse[k] 
  })
  bootse_CI_scad_90 = sapply(1:length(as.vector(parallest1$scad)),function(k){
    as.vector(parallest1$scad)[k]+c(-1,1)*1.64*bootse[k] 
  })
  boot_CI_scad_res = rbind(ci_lower,ci_upper)
  boot_CI_scad_90 = rbind(ci_lower2,ci_upper2)
  
  # wang's method
  bootmat <- do.call(rbind,lapply(boot_data_results_res,function(est) as.vector(est[[2]]$robustrg)))
  ci_lower <- apply(bootmat,2,quantile,probs=0.025)
  ci_upper <- apply(bootmat,2,quantile,probs=0.975)
  bootse <- apply(bootmat,2,sd)
  bootse_CI_wang_res = sapply(1:length(as.vector(deconfoundest1$robustrg)),function(k){
    as.vector(deconfoundest1$robustrg)[k]+c(-1,1)*1.96*bootse[k] 
  })
  boot_CI_wang_res = rbind(ci_lower,ci_upper)
  return(list(S = list(greed1,estS1),beta = parallest1, wang = deconfoundest1,
              bootstrap = list(boot = boot_data_results_res,
                               scad = list(bootse_CI_scad_res,
                                           boot_CI_scad_res,bootse_CI_scad_90,boot_CI_scad_90),
                               wang = list(bootse_CI_wang_res,
                                           boot_CI_wang_res))))
}

