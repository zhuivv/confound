# de-bias functions
require(MASS);require(ncvreg);require(glmnet);require(nFactors)

### parallel outcome two stage method
debias_twostage <- function(X,Y,U = testindp$data$U,ptrue = 1:15){
  ## simple linear regression with X and U
  beta_orc = lm(Y~X+U)$coef[2:(ncol(X)+1),]
  
  ## linear regression with only X
  beta_x = lm(Y~X)$coef[2:(ncol(X)+1),]
  
  # factor analysis for X
  hatr = nScree(cor(X))$Components$nkaiser
  
  pzero = setdiff(1:ncol(Y),ptrue)
  ## parallel Y 2023
  betal=pl=betarg=
    betargms=betalasso=epsilon=betalms=
    betascad = betax =  miaobeta = 
    matrix(0,ncol = ncol(Y),nrow = ncol(X))
  
  for(i in 1:ncol(Y)){
    epsilon[,i] = lm(Y[,i]~X)$coef[2:(ncol(X)+1)]
    # multiple trt debias
    miaobeta[,i] = esti.null(Y[,i],X,v = NULL,hatr)
  }
  
  for(i in ptrue){
    # step one
    # X = scale(X);Y = scale(Y)
    W = Y[,pzero]
    Zl = Y[,setdiff(ptrue,i)]
    modelone = lm(W~X+Zl)
    What = predict(modelone)
    # step two
    modeltwo = lm(Y[,i]~X+What)
    betal[,i] = modeltwo$coef[2:(ncol(X)+1)]
    betax[,i] = esti.null(Y[,i],X,v = NULL,hatr)
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
    # SE's
    # summary_fit <- summary(scad_fit)
    # lambda_index <- which.min(abs(summary_fit$lambda - scad_lam))
    # sdscad[,i] <- (summary_fit$se.asymptotic[lambda_index, ])[2:(ncol(X)+1)]
  }
  return(list(oracle = beta_orc,naive = beta_x,parall = betal,ridge = betarg,
              miaobeta = miaobeta,betax = betax,
              adlasso = betargms,lasso = betalasso,scad = betascad,direct = lm(Y~X)$coef[2:(ncol(X)+1),]))
}

# function for the null treatments estimation from Miao 2022
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

#### confounded testing methods (try other regression methods rather than robust regression)
debias_deconfound <- function(X,Y,U = testindp$data$U,ptrue = 1:15){
  options(warn = -1)
  ## simple linear regression with X and U
  beta_orc = lm(Y~X+U)$coef[2:(ncol(X)+1),]
  
  ## linear regression with only X
  # beta_x = lm(Y~X)$coef[2:(ncol(X)+1),]
  
  ## deconfounding test 2017
  ### negative control
  testx = qr(X)
  Yqr = qr.qty(testx,Y)
  testy = t(qr.Q(testx))%*%Y
  # fa
  hatr = nScree(cor(X))$Components$nkaiser
  fahat = psych::fa(cor(Yqr[(ncol(X)+1):nrow(Yqr),]),hatr)
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
  scadmodel_selection = apply(t(Yqr[1:ncol(X),])%*%t(solve(qr.R(testx))),2,function(x) (ncvreg(unclass(gammahat),x,nlambda = 100)))
  scad_lam <- lapply(scadmodel_selection, function(x) x$lambda[which.min(BIC(x))])
  betascad = mapply(function(lam,model) coef(model, lambda=lam),scad_lam,scadmodel_selection)
  betahat_sps = t(Yqr[1:ncol(X),])%*%t(solve(qr.R(testx)))-gammahat%*%robustmodel[2:nrow(robustmodel),]
  betahat_scad = t(Yqr[1:ncol(X),])%*%t(solve(qr.R(testx)))-gammahat%*%betascad[2:nrow(betascad),]
  # ## try LASSO for sparsity (LEAPP)
  
  return(list(oracle = beta_orc,robustrg = t(betahat_sps),scadrg = t(betahat_scad)))
}



