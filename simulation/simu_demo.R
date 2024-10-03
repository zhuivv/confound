# simulation demo
## evaluation functions
estS_eval <- function(input){
  return(list(
    sens_mean = mean(unlist(sapply(input[sapply(input, length) == 5],function(x) x[[3]])),na.rm = T),
    FDR_mean = mean(unlist(sapply(input[sapply(input, length) == 5],function(x) x[[4]])),na.rm = T),
    FDRp_mean = mean(unlist(sapply(input[sapply(input, length) == 5],function(x) x[[5]])),na.rm = T),
    sens_sd = sd(unlist(sapply(input[sapply(input, length) == 5],function(x) x[[3]])),na.rm = T),
    FDR_sd = sd(unlist(sapply(input[sapply(input, length) == 5],function(x) x[[4]])),na.rm = T),
    FDRp_sd = sd(unlist(sapply(input[sapply(input, length) == 5],function(x) x[[5]])),na.rm = T)
  ))
}

eval_beta <- function(input){
  len = 7
  MAE = do.call(cbind,lapply(input[sapply(input, length) == len],function(x) x[[5]]))
  infy = do.call(cbind,lapply(input[sapply(input, length) == len],function(x) x[[6]]))
  # coverage = do.call(cbind,lapply(input,function(x) x[[8]]))
  return(list( MAE_result = paste(rownames(input[[1]][[5]]),":",round(rowMeans(MAE)*100,2),"(",round(apply(MAE,1,sd)*100,3),")",sep = ""),
               infy_result = paste(rownames(input[[1]][[5]]),":",round(rowMeans(infy)*100,2),"(",round(apply(infy,1,sd)*100,2),")",sep = "")))
  # prob_result = paste(rownames(input[[1]][[5]]),":",round(rowMeans(coverage)*100,4),"(",round(apply(coverage,1,sd),5),")",sep = "")))
}



################################
## beta estimation with true S
library(parallel)
ncore = detectCores()
cl <- makeCluster(ncore-2)
clusterEvalQ(cl, {
  library(glmnet)
  library(MASS)
  library(ncvreg)
  library(nFactors)
})
# for genearl strong and weak effect:
eff = 5; q = 5; r = 3 
## adjust effect size eff, number of exposures q, number of unmeasured confounding r as needed
clusterExport(cl, c("simutwoX", "debias_twostage", "debias_deconfound", "est_bias","eff","q","r","process_function","esti.null"))  # add any other required objects
B <- 1:500  # Assuming B is defined somewhere in your code
set.seed(233);start_time = Sys.time()
bias_list_eff5_q5_r3 <- parLapply(cl, B, process_function)
Sys.time()-start_time
eval_beta(bias_list_eff5_q5_r3)

# for mixed effect:
q = 5; r = 3;p=100 # adjust as needed, p is the number of outcomes
clusterExport(cl, c("simutwoX_mix", "debias_twostage", "debias_deconfound", "est_bias","eff","q","r","process_function_mix","esti.null","bootstrap","calculate_coverage"))  # add any other required objects
# B <- 1:500  # Assuming B is defined somewhere in your code
set.seed(233);start_time = Sys.time()
bias_list_q5_r3_mix <- parLapply(cl, B, process_function_mix)
Sys.time()-start_time
eval_beta(bias_list_q5_r3_mix)

# for uniform residuals
# strong effect, 3 exposures, r = 3
eff = 5; q = 3; r = 3
clusterExport(cl, c("simutwoX_sparse", "debias_twostage", "debias_deconfound", "est_bias","eff","q","r","process_function_sparse","esti.null"))  # add any other required objects
# B <- 1:500  # Assuming B is defined somewhere in your code
set.seed(233);start_time = Sys.time()
bias_list_eff5_q3_r3_sparse <- parLapply(cl, B, process_function_sparse)
Sys.time()-start_time
eval_beta(bias_list_eff5_q3_r3_sparse)
stopCluster()


##############################
# beta estimation with estimated S
## also include S estimation evaluations
debias_twostage_fdr <- function(X,Y,U = testindp$data$U,ptrue = 1:15){
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
  
  for(i in 1:ncol(Y)){
    # step one
    # X = scale(X);Y = scale(Y)
    W = Y[,setdiff(pzero,i)]
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

process_functionS <- function(b) {
  # parameters
  set.seed(123+b)
  n = n; p = p;rate = rate
  eff1 = -1; eff2 = 5; q = q; r = 3
  # function
  testHDX <- simutwoX(eff2,n,p,q,r)
  tryCatch({
    if(length(unique(unlist(testHDX$coef$tauset_comp)))==ncol(testHDX$data$Y)){
      stop("No NCO")
    }else{
      greed = greedy_l0(testHDX$data$Y,ncol(testHDX$data$X),testHDX$data$X)
      estS = greed$truest
      # fnset = sample(trueS,round(length(setdiff(1:p,trueS))*rate))
      # parallest <- debias_twostage_fdr(testHDX$data$X, testHDX$data$Y, testHDX$data$U, trueS[!trueS%in%fnset])
      parallest <- debias_twostage_fdr(testHDX$data$X, testHDX$data$Y, testHDX$data$U, estS)
      deconfoundest <- debias_deconfound(testHDX$data$X, testHDX$data$Y, testHDX$data$U, estS)
      
      testcoef1 <- est_bias(parallest)
      testcoef2 <- est_bias(deconfoundest)
      
      combined_diff <- as.data.frame(rbind(do.call(rbind, lapply(testcoef1$MSE, mean,na.rm = T)),
                                           do.call(rbind, lapply(testcoef2$MSE, mean,na.rm = T))))
      infy_diff <- as.data.frame(rbind(do.call(rbind, lapply(testcoef1$infy_bias, mean,na.rm = T)),
                                       do.call(rbind, lapply(testcoef2$infy_bias, mean,na.rm = T))))
      
      # evaluation fo estS
      trueS = unique(unlist(testHDX$coef$tauset_comp))
      sens = length(intersect(setdiff(1:p,estS),setdiff(1:p,trueS)))/length(setdiff(1:p,trueS))
      FDR = length(intersect(setdiff(1:p,estS),trueS))/length(setdiff(1:p,estS))
      FDRp = length(intersect(setdiff(1:p,estS),trueS))/p
      return(list(estbeta = list(parallest,deconfoundest,testcoef1,
                                 testcoef2,combined_diff,infy_diff,testHDX$coef$beta),
                  estS = list(greed,trueS,sens,FDR,FDRp)))
      #, coverage_prob))
    }
  },error = function(e){
    # cat("Error:", e$message, "\n")
    return(e$message)
  })
  # combined_diff$condition <- c(names(testcoef1), names(testcoef2))
  
}

cl <- makeCluster(ncore-2)
clusterEvalQ(cl, {
  library(glmnet)
  library(MASS)
  library(ncvreg)
  library(nFactors)
})
q = 3; rate = 0.05; p = 100;n=100
clusterExport(cl, c("simutwoX", "debias_twostage_fdr","greedy_l0","model_list", "debias_deconfound", "est_bias","n","eff","q","r","rate","p","process_functionS","esti.null","bootstrap","calculate_coverage"))  # add any other required objects
B <- 1:100  # Assuming B is defined somewhere in your code
set.seed(233);start_time = Sys.time()
bias_n100_q3p100 <- parLapply(cl, B, process_functionS)
Sys.time()-start_time
# for sample size n = 100, 1000, the outcomes are bias_n100_q3p100 
## and bias_n1000_q3p100
## similar to the rest of the results
eval_beta(lapply(bias_n100_q3p100[which(lapply(bias_n100_q3p100,length)==2)], function(x) x$estbeta))
estS_eval(lapply(bias_n100_q3p100[which(lapply(bias_n100_q3p100,length)==2)],function(x) x$estS))





